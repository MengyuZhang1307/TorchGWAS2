#!/usr/bin/env python3

import sys
import os
import numpy as np
import torch
from tqdm import tqdm
import math


sys.path.append(os.path.dirname(os.path.abspath(__file__)))
# Add path for GEM module
sys.path.append("build")
sys.path.append("pymodules")

## function to read corrected residuals and C2
def read_correction_file(file_path: str):
    """
    @brief Read phenotype file with special format.

    @param file_path Path to the correction file.
    @return (header, c2, data_array)
        - header: list of phenotypes name (column names)
        - c2: numpy array of double (values from 2nd line after '#')
        - c_res: numpy 2D array of double (corrected residuals values)
    """
    c2 = np.empty(0, dtype=np.float64)
    c_res = []
    
    with open(file_path, "r") as f:
        header = f.readline().strip().split("\t")
        second_line = f.readline().strip().split("\t")
        if second_line[0].startswith("#"):
            c2 = np.array([float(val) for val in second_line if not val.startswith("#")],
              dtype=np.float64)
                          

        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                c_res.append([float(x) for x in parts[2:]])

    c_res = np.array(c_res, dtype=np.float64)

    return header, c2, c_res


def calc_t(corrected_res, geno, beta, gamma, sqrt_c2, ph_std):
    """
    Compute per-SNP stats using pre-normalized phenotypes and on-device sqrt(c2).

    corrected_res: (N, P)
    geno: (M, N)
    beta, gamma: (M, P) work buffers on same device as geno
    sqrt_c2_device: (P,) tensor on same device as inputs
    ph_std: (1, P) phenotype stds computed BEFORE standardizing corrected_res
    """
    N = corrected_res.shape[0]
    with torch.no_grad():
        # Compute stds (keep dims for broadcasting)
        # geno_std: (M, 1) across samples; ph_std: (1, P) provided from pre-standardized corrected_res
        geno_std = geno.std(1, keepdim=True, unbiased=False).clamp_min(1e-8)
        ph_std = ph_std.clamp_min(1e-8)

        # Row-wise center and scale genotypes using saved std (so we can reuse geno_std later)
        geno.sub_(geno.mean(1, keepdim=True)).div_(geno_std)

        # Score U = X^T Y written into beta (M,P); then convert to r = U/N
        torch.matmul(geno, corrected_res, out=beta)  # (M,N) @ (N,P) -> (M,P)
        beta.div_(N)  # now beta holds r (correlation) when inputs are standardized

        # Keep an unscaled copy of r for SE(r)
        r = beta.clone()  # (M,P)

        # As requested: additionally scale beta by phenotype and SNP stds
        # beta: (M,P) / (1,P) / (M,1) -> (M,P)
        beta.div_(ph_std)
        beta.div_(geno_std)

        # Compute SE from r, then scale SE by the same stds
        gamma.copy_(r)                # start from r
        gamma.pow_(2).sub_(1).div_(2 - N)  # (1 - r^2)/(N - 2)
        torch.sqrt(gamma, out=gamma)  # SE(r)
        gamma.div_(ph_std)            # adjust SE for phenotype scaling
        gamma.div_(geno_std)          # adjust SE for SNP scaling

        # Null-model calibration
        gamma.div_(sqrt_c2.unsqueeze(0))

        beta_coeffs = beta.cpu()
        se = gamma.cpu()
        t_stats = beta.div_(gamma).abs_().neg_().cpu()
        return t_stats, beta_coeffs, se


def run_gwas(runner, snps_per_chunk=1000, device='cuda',  compress=False):
    """
    Run GWAS using a pre-configured GEMRunner instance.
    """
    # Get C2 values from the runner (try fitting null model first)
    # c2_values = None
    # try:
    #     runner.run_fit_nullmodel()
    #     c2_values = runner.get_c2_values()
    #     print(f"Fitted null model and obtained C2 values for SE adjustment: {c2_values}")
    # except Exception as e:
    #     print(f"Could not fit null model or get C2 values: {e}")
    #     print("Proceeding without C2 adjustment")
    # if c2_values is None:
    #     return
    
    # Get corrected residuals from runner SHOULD GET CORRECTED SCALED RESIDUALS
    ph_headers, c2_values, corrected_res = read_correction_file("outAddlie.txt") # read corrected_res, c2 and ph_headers from file
    # corrected_res = torch.from_numpy(runner.get_phenotypes()).float()

    #intercept = torch.from_numpy(runner.get_covariates()).float()
    
    if device == 'cuda' and torch.cuda.is_available():
        device = torch.device('cuda')
    else:
        device = torch.device('cpu')

    # corrected_res = corrected_res.to(device)
    corrected_res = torch.from_numpy(corrected_res).float().to(device)
    #covariates = covariates.to(device)
    
    n_samples, n_corrected_res = corrected_res.shape

    # Center phenotypes with NaN-safe mean and replace NaNs
    #ph_mean = torch.nanmean(phenotypes, dim=0, keepdim=True)
    #phenotypes = phenotypes - ph_mean
    corrected_res = torch.nan_to_num(corrected_res, nan=0.0)

    #c = covariates.cpu().numpy()
    #c_mean = np.nanmean(c, axis=0, keepdims=True)
    #c_std = np.nanstd(c, axis=0, keepdims=True)
    #c_std[c_std < 1e-12] = 1.0
    #c = (c - c_mean) / c_std
    #c = np.nan_to_num(c, nan=0.0)
    #covarQ, _ = np.linalg.qr(c)
    #covarQ = torch.from_numpy(covarQ).float().to(device)
    
    #pheno_normalized = phenotypes - torch.matmul(covarQ, torch.matmul(covarQ.T, phenotypes))
    # Save phenotype std BEFORE normalization for scaling beta/gamma
    ph_std_pre = torch.std(corrected_res, dim=0, keepdim=True, unbiased=False).clamp_min(1e-8)
    corrected_res = corrected_res / ph_std_pre
    corrected_res = torch.nan_to_num(corrected_res, nan=0.0)
    
    # Precompute sqrt(c2) once on the target device (clamped for stability)
    sqrt_c2 = torch.from_numpy(np.asarray(c2_values)).to(device=device, dtype=torch.float32)
    sqrt_c2.sqrt_()

    # Start dosage streaming from runner
    queue = runner.start_dosage_stream(queue_capacity=20, snps_per_chunk=snps_per_chunk)
    
    all_t_stats = []
    all_beta = []
    all_se = []
    
    # Preallocate device buffers and reuse/slice for smaller final chunks
    geno_tensor = torch.empty(snps_per_chunk, n_samples, device=device)
    beta_tensor = torch.empty(snps_per_chunk, n_corrected_res, device=device)
    gamma_tensor = torch.empty(snps_per_chunk, n_corrected_res, device=device)
    

     # --- Prepare output files ---
    buffer_snps = 500_000
    out_prefix = "TGWAS"
    def _open(path, mode):
        if compress:
            return gzip.open(path + ".gz", mode + "t")
        return open(path, mode, buffering=1024*1024)

    handles = {}
    header_line = "BETA\tSE\tT_STAT\tP_VALUE\n"
    ph_headers = ph_headers[2:] 
    n_files = 600
    n_pheno = len(ph_headers)

    pheno_per_file = math.ceil(n_pheno / n_files)  # auto-calc phenos per file

    # --- Open grouped files ---
    handles = []
    for f in range(n_files):
        start = f * pheno_per_file
        end   = min((f + 1) * pheno_per_file, n_pheno)
        if start >= end:
            break
        group_headers = ph_headers[start:end]

        # Create file and write header
        path = f"{out_prefix}_block{f+1}.tsv"
        fh = _open(path, "w")
        header_cols = []
        for ph in group_headers:
            header_cols.append(f"BETA_{ph}")
            header_cols.append(f"SE_{ph}")
            header_cols.append(f"TSTAT_{ph}")
            header_cols.append(f"PVAL_{ph}")
        fh.write("\t".join(header_cols) + "\n")

        handles.append((start, end, fh))
    
    buffer = {fh: [] for _, _, fh in handles}
    snp_processed = 0
    for chunk_data in tqdm(queue, desc="Processing SNPs"):
            
        actual_snps = chunk_data.shape[0]
        
        if actual_snps == snps_per_chunk:
            geno = geno_tensor
            beta = beta_tensor
            gamma = gamma_tensor
        else:
            geno = geno_tensor[:actual_snps, :]
            beta = beta_tensor[:actual_snps, :]
            gamma = gamma_tensor[:actual_snps, :]

        # geno.copy_(chunk_data)
        # geno.copy_(torch.from_numpy(chunk_data).float().to(device))
        geno.copy_(torch.from_numpy(chunk_data).float())
        t_stats, beta_coeffs, se = calc_t(corrected_res, geno, beta, gamma, sqrt_c2, ph_std_pre)
        pvals = 2 * torch.special.ndtr(t_stats) #added by Sama
        # all_t_stats.append(t_stats)
        # all_beta.append(beta_coeffs)
        # all_se.append(se)
    

    # Convert to numpy
        t_np   = t_stats.cpu().numpy()
        b_np   = beta_coeffs.cpu().numpy()
        se_np  = se.cpu().numpy()
        p_np   = pvals.cpu().numpy()

        # Write per phenotypes
        for start, end, fh in handles:
            block_lines = []
            for r in range(actual_snps):
                row_vals = []
                for j in range(start, end):
                    row_vals.append(f"{b_np[r,j]:.6g}")
                    row_vals.append(f"{se_np[r,j]:.6g}")
                    row_vals.append(f"{t_np[r,j]:.6g}")
                    row_vals.append(f"{p_np[r,j]:.6g}")
                block_lines.append("\t".join(row_vals) + "\n")
            buffer[fh].extend(block_lines)

        snp_processed += actual_snps

        # Flush if buffer full
        if snp_processed >= buffer_snps:
            for fh, lines in buffer.items():
                if lines:
                    fh.writelines(lines)
                    buffer[fh] = []  # clear buffer
            snp_processed = 0

    # --- Final flush ---
    for fh, lines in buffer.items():
        if lines:
            fh.writelines(lines)
        fh.close()
    print(f"Wrote {len(ph_headers)} phenotype files with prefix {out_prefix}_*.tsv{'.gz' if compress else ''}")
    
    # t_statistics = torch.cat(all_t_stats, dim=0)
    # beta_coefficients = torch.cat(all_beta, dim=0)
    # standard_errors = torch.cat(all_se, dim=0)
    # p_values = 2 * torch.special.ndtr(t_statistics)
    # return {
    #     't_stats': t_statistics,
    #     'beta': beta_coefficients,
    #     'se': standard_errors,
    #     'p_values': p_values,
    #     'ph_headers': ph_headers,
    # }