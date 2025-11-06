#!/usr/bin/env python3

import sys
import os
import numpy as np
import torch
from tqdm import tqdm
import math
import time
import pandas as pd
import pyarrow as pa, pyarrow.parquet as pq
import gc

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
            if len(parts) >= 2:
                c_res.append([float(x) for x in parts[1:]])

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
        beta.mul_(ph_std)
        beta.div_(geno_std)

        # Compute SE from r, then scale SE by the same stds
        gamma.copy_(r)                # start from r
        gamma.pow_(2).sub_(1).div_(2 - N)  # (1 - r^2)/(N - 2)
        torch.sqrt(gamma, out=gamma)  # SE(r)
        gamma.mul_(ph_std)            # adjust SE for phenotype scaling
        gamma.div_(geno_std)          # adjust SE for SNP scaling

        # Null-model calibration
        gamma.div_(sqrt_c2.unsqueeze(0))

        # Extract final beta and SE BEFORE computing t-stats
        beta_coeffs = beta.clone().cpu()
        se = gamma.clone().cpu()
        t_stats = beta.div_(gamma).abs_().neg_().cpu()
        return geno_mean, geno_std, t_stats, beta_coeffs, se


def run_gwas(runner, out_file, snps_per_chunk=1000, device='cuda',  compress=False):
    """
    Run GWAS using a pre-configured GEMRunner instance.
    """
 
    ph_headers, c2_values, corrected_res = read_correction_file("intermediate_" + out_file) # read corrected_res, c2 and ph_headers from file
    
    if device == 'cuda' and torch.cuda.is_available():
        device = torch.device('cuda')
    else:
        device = torch.device('cpu')

    # corrected_res = corrected_res.to(device)
    corrected_res = torch.from_numpy(corrected_res).float().to(device)
    #covariates = covariates.to(device)
    
    n_samples, n_corrected_res = corrected_res.shape

    # Center phenotypes with NaN-safe mean and replace NaNs
    corrected_res = torch.nan_to_num(corrected_res, nan=0.0)

    # Save phenotype std BEFORE normalization for scaling beta/gamma
    ph_std_pre = torch.std(corrected_res, dim=0, keepdim=True, unbiased=False).clamp_min(1e-8)
    corrected_res = corrected_res / ph_std_pre
    corrected_res = torch.nan_to_num(corrected_res, nan=0.0)
    
    # Precompute sqrt(c2) once on the target device (clamped for stability)
    sqrt_c2 = torch.from_numpy(np.asarray(c2_values)).to(device=device, dtype=torch.float32)
    sqrt_c2.sqrt_()

    # Start dosage streaming from runner
    queue = runner.start_dosage_stream(queue_capacity=20, snps_per_chunk=snps_per_chunk)
    
    ph_headers = ph_headers[1:]    
    # Preallocate device buffers and reuse/slice for smaller final chunks
    geno_tensor = torch.empty(snps_per_chunk, n_samples, device=device)
    beta_tensor = torch.empty(snps_per_chunk, n_corrected_res, device=device)
    gamma_tensor = torch.empty(snps_per_chunk, n_corrected_res, device=device)
    
    buffer_size = 100_000

    headers = [
    "SNPID",
    "RSID",
    "CHR",
    "POS",
    "Non_Effect_Allele",
    "Effect_Allele",
    "N_Samples",
    "AF",
    "GV"
    ]

    for ph in ph_headers:
        headers.append(f"{ph}_BETA")
        headers.append(f"{ph}_SE")
        headers.append(f"{ph}_pvalue")

    rows_in_buffer = 0
    buffer = []
    writer = None
    out_path = "TGWAS_" + out_file
    if os.path.exists(out_path + ".parquet"):
        os.remove(out_path + ".parquet")

    start = time.time()

    for chunk_data, meta in tqdm(queue, desc="Processing SNPs"):
        actual_snps = chunk_data.shape[0]
        rows_in_buffer += actual_snps

        geno = geno_tensor[:actual_snps, :]
        beta = beta_tensor[:actual_snps, :]
        gamma = gamma_tensor[:actual_snps, :]

        geno.copy_(torch.from_numpy(chunk_data).float())
        mean, std, t_stats, beta_coeffs, se = calc_t(
            corrected_res, geno, beta, gamma, sqrt_c2, ph_std_pre
        )

        b_np = beta_coeffs.cpu().numpy().astype(np.float32)
        se_np = se.cpu().numpy().astype(np.float32)
        neg_log10_pval = (
            -(torch.log(torch.tensor(2.0)) + torch.special.log_ndtr(t_stats))
            / torch.log(torch.tensor(10.0))
        ).cpu().numpy().astype(np.float32)

        # keep order: BETA, SE, PVAL repeating per phenotype
        all_stats = np.stack([b_np, se_np, neg_log10_pval], axis=2)
        all_stats_2d = all_stats.reshape(b_np.shape[0], -1)

        # Convert metadata to Arrow arrays
        meta_arrays = []
        for k, v in meta.items():
            if isinstance(v[0], str):
                # v = [s.rstrip('\x00') for s in v]
                meta_arrays.append(pa.array(v, type=pa.string()))
            elif isinstance(v[0], (int, np.integer)):
                meta_arrays.append(pa.array(v, type=pa.int32()))
            else:
                meta_arrays.append(pa.array(v, type=pa.float32()))

        # Convert numeric results to Arrow arrays
 
        stat_arrays= [pa.array(col, type=pa.float32()) for col in all_stats_2d.T]
                        
                    
        table= pa.table(meta_arrays + stat_arrays, names=headers)
        buffer.append(table)

        # flush buffer
        if rows_in_buffer >= buffer_size:
            combined = pa.concat_tables(buffer)
            # write Feather (Arrow IPC v2)
            # feather.write_feather(combined, out_path + ".feather", compression="zstd")
            if writer is None:
                writer = pq.ParquetWriter(
                out_path + ".parquet", combined.schema, compression="snappy"
            )
            writer.write_table(combined)
            rows_in_buffer = 0
            buffer.clear()
            del combined, b_np, se_np, neg_log10_pval,
            all_stats, all_stats_2d, table, stat_arrays, meta_arrays
            gc.collect()
            torch.cuda.empty_cache()

    # --- Flush remaining ---
    if buffer:
        combined = pa.concat_tables(buffer)
        # feather.write_feather(combined, out_path + ".feather", compression="zstd")
        if writer is None:
                writer = pq.ParquetWriter(
                out_path + ".parquet", combined.schema, compression="snappy"
            )
        writer.write_table(combined)
        buffer.clear()
        del combined, b_np, se_np, neg_log10_pval,
        all_stats, all_stats_2d, table, stat_arrays, meta_arrays
        gc.collect()
        torch.cuda.empty_cache()
    if writer is not None:
        writer.close()
    end = time.time()
    print(f"time for chunk = {end - start:.2f}s")
