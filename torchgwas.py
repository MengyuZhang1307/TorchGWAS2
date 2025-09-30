#!/usr/bin/env python3

import sys
import os
import numpy as np
import torch
from tqdm import tqdm
import math
import time
import pandas as pd


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
        # geno.sub_(geno.mean(1, keepdim=True)).div_(geno_std)
        geno_mean = geno.mean(1, keepdim=True)
        geno.sub_(geno_mean).div_(geno_std)
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
        return geno_mean, geno_std, t_stats, beta_coeffs, se


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
    
    ph_headers = ph_headers[2:]    
    # Preallocate device buffers and reuse/slice for smaller final chunks
    geno_tensor = torch.empty(snps_per_chunk, n_samples, device=device)
    beta_tensor = torch.empty(snps_per_chunk, n_corrected_res, device=device)
    gamma_tensor = torch.empty(snps_per_chunk, n_corrected_res, device=device)
    
    all_beta = []
    all_se = []
    all_t = []
    all_p = []
    all_mean = []
    all_std  = []
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
        mean, std, t_stats, beta_coeffs, se = calc_t(corrected_res, geno, beta, gamma, sqrt_c2, ph_std_pre)
        pvals = 2 * torch.special.ndtr(t_stats) #added by Sam

        # # Convert to numpy
        mean_np = mean.cpu().numpy()
        std_np = std.cpu().numpy()
        b_np   = beta_coeffs.cpu().numpy()
        se_np  = se.cpu().numpy()
        t_np   = t_stats.cpu().numpy()
        p_np   = pvals.cpu().numpy()

        all_mean.append(mean_np)
        all_std.append(std_np)
        all_beta.append(b_np)
        all_se.append(se_np)
        all_t.append(t_np)
        all_p.append(p_np)
        

    all_beta = np.vstack(all_beta)   # (num_snps, num_pheno)
    all_se   = np.vstack(all_se)
    all_t    = np.vstack(all_t)
    all_p    = np.vstack(all_p)
    all_mean = np.concatenate(all_mean) if all_mean else np.array([])
    all_std  = np.concatenate(all_std)  if all_std  else np.array([])
    all_stats = np.stack([all_beta, all_se, all_t, all_p], axis=2)
    num_snps, num_pheno = all_beta.shape
    print("num_snps", num_snps)
    print("num_pheno", num_pheno)
    # Flatten phenotypes → shape (num_snps, num_pheno*4)
    all_stats_2d = all_stats.reshape(num_snps, -1)
    # Add mean/std at the start → shape (num_snps, 2 + num_pheno*4)
    final_mat = np.concatenate(
        [all_mean.reshape(-1, 1), all_std.reshape(-1, 1), all_stats_2d],
        axis=1
    )
    headers = ["snp_mean", "snp_std"]
    for j in range(len(ph_headers)):
        headers.extend([f"beta_{j+1}", f"se_{j+1}", f"t_{j+1}", f"pval_{j+1}"])
    start_time_df = time.time()
    df = pd.DataFrame(final_mat, columns=headers)
    end_time_df = time.time()
    print("Wall time in seconds  for converting to df:", end_time_df - start_time_df)
    # Save once
    start_time_save = time.time()
    df.to_parquet("results_all_ADDLIE.parquet", engine="pyarrow", compression="snappy")
    end_time_save = time.time()
    print("Wall time in seconds for saving df:", end_time_save - start_time_save)
    df = pd.read_parquet("results_all_ADDLIE.parquet", engine="pyarrow")
    print(df.head())