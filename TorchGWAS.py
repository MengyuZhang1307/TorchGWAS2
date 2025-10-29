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
import json

# sys.path.append(os.path.dirname(os.path.abspath(__file__)))
# Add path for GEM module
# sys.path.append("build_test")
# sys.path.append("pymodules")

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
    
  
    # # buffer_size = 500_000
    # buffer = []
    # beta_se_headers = []
    # for ph in ph_headers:
    #     beta_se_headers.append(f"{ph}_BETA")
    #     beta_se_headers.append(f"{ph}_SE")
    #     beta_se_headers.append(f"{ph}_pvalue")

    # start = time.time()

    # # Output and metadata setup
    # out_path = "TGWAS" + out_file
    # meta = {
    # "cols": len(ph_headers),
    # "total_cols": len(beta_se_headers),
    # "dtype": "float32",
    # "layout": "Beta_SE_rowmajor",
    # "headers": beta_se_headers,        
    # "buffer_flush_rows": 100_000,
    # "rows": 0
    # }

    # if os.path.exists(out_path):
    #     os.remove(out_path)
    # out = open(out_path, "ab")  # append binary

    # # Buffers
    # buffer = []
    # rows_in_buffer = 0
    # total_rows = 0

    # for chunk_data in tqdm(queue, desc="Processing SNPs"):
    #     actual_snps = chunk_data.shape[0]
    #     total_rows += actual_snps
    #     rows_in_buffer += actual_snps

    #     geno = geno_tensor[:actual_snps, :]
    #     beta = beta_tensor[:actual_snps, :]
    #     gamma = gamma_tensor[:actual_snps, :]

    #     geno.copy_(torch.from_numpy(chunk_data).float())
    #     mean, std, t_stats, beta_coeffs, se = calc_t(
    #         corrected_res, geno, beta, gamma, sqrt_c2, ph_std_pre
    #     )

    #     b_np = beta_coeffs.cpu().numpy().astype(np.float32)
    #     se_np = se.cpu().numpy().astype(np.float32)
    #     t_stats_bp = t_stats.numpy().astype(np.float32)
    #     # all_stats = np.stack([b_np, se_np], axis=2)
    #     # pvals = 2 * torch.special.ndtr(t_stats)
    #     neg_log10_pval = -(torch.log(torch.tensor(2.0)) + torch.special.log_ndtr(t_stats)) / torch.log(torch.tensor(10.0))
    #     all_stats = np.stack([b_np, se_np, neg_log10_pval], axis=2)
    #     all_stats_2d = all_stats.reshape(b_np.shape[0], -1)
    #     buffer.append(all_stats_2d)

    #     #  Flush every ~1M SNPs
    #     if rows_in_buffer >= meta["buffer_flush_rows"]:
    #         stack_buffer = np.vstack(buffer)
    #         stack_buffer.tofile(out)
    #         buffer.clear()
    #         rows_in_buffer = 0 
    #         out.flush()           


    # if rows_in_buffer > 0:
    #     stack_buffer = np.vstack(buffer)
    #     stack_buffer.tofile(out)
    #     buffer.clear()
    #     rows_in_buffer = 0  
    #     print("There are snps less than buffer size, add them to output")

    # out.close()

    # # Write metadata JSON
    # meta["rows"] = total_rows
    # with open(out_path + ".meta.json", "w") as f:
    #     json.dump(meta, f, indent=2)

    # end = time.time()
    # print(f" Done — {total_rows:,} SNPs written in {(end-start):.1f}s")

    # # Load metadata
    # with open(out_path + ".meta.json") as f:
    #     meta = json.load(f)

    # rows = meta["rows"]
    # cols = meta["cols"]

    # data = np.fromfile(out_path, dtype=np.float32)

    # # Reshape into (rows, cols*2)
    # data = data.reshape(rows, cols * 2)

    # # Print the first 5 rows
    # print("First 5 rows (Beta + SE + pvalue):")
    # print(data[:5, :20])   
    writer = None 
    buffer_size = 100_000
    rows_in_buffer = 0
    buffer = []
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

    start = time.time()

    # Output and metadata setup
    out_path = "TGWAS" + out_file
 
    if os.path.exists(out_path):
        os.remove(out_path)
    out = open(out_path, "ab")  # append binary

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
        t_stats_bp = t_stats.numpy().astype(np.float32)
        # all_stats = np.stack([b_np, se_np], axis=2)
        # pvals = 2 * torch.special.ndtr(t_stats)
        neg_log10_pval = -(torch.log(torch.tensor(2.0)) + torch.special.log_ndtr(t_stats)) / torch.log(torch.tensor(10.0))
        all_stats = np.stack([b_np, se_np, neg_log10_pval], axis=2)
        all_stats_2d = all_stats.reshape(b_np.shape[0], -1)
        # Convert metadata and stats to DataFrame
        df_meta = pd.DataFrame(meta)
        df_stats = pd.DataFrame(
            all_stats_2d,
            columns=[f"{ph}_BETA" for ph in ph_headers] +
                    [f"{ph}_SE" for ph in ph_headers] +
                    [f"{ph}_pvalue" for ph in ph_headers]
        )
        df = pd.concat([df_meta, df_stats], axis=1)

        buffer.append(df)

        
        if rows_in_buffer >= buffer_size:
            # Vertically stack all arrays from the buffer
            df = pd.concat(buffer, ignore_index=True)
            # Convert to Arrow Table and write to Parquet
            table = pa.Table.from_pandas(df, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(
                out_path + ".parquet", table.schema, compression="snappy"
            )
                
            writer.write_table(table)
            rows_in_buffer = 0
            # Clear buffer after writing
            buffer.clear()
    # flush remainder
    if buffer:
        df = pd.concat(buffer, ignore_index=True)
        # Create a DataFrame (column names already known)
        # df = pd.DataFrame(stacked, columns=headers)
        # Convert all numeric columns explicitly to float32 (for safety)
        # df = df.astype(np.float32, copy=False)
        table = pa.Table.from_pandas(df, preserve_index=False)
        if writer is None:
            writer = pq.ParquetWriter(
                out_path + ".parquet", table.schema, compression="snappy"
            )
        writer.write_table(table)

    # Close writer
    if writer:
        writer.close()
    
    end = time.time()
    print(f"time for chunck = {end - start}")

    #Read the parquet file head
    df_check = pd.read_parquet(out_path + ".parquet")

    # Show first 5 rows
    print(df_check.head())

    # Show the column names
    print(df_check.columns.tolist())