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
    sample_ids = []
    
    with open(file_path, "r") as f:
        header = f.readline().strip().split("\t")
        second_line = f.readline().strip().split("\t")
        if second_line[0].startswith("#"):
            c2 = np.array([float(val) for val in second_line if not val.startswith("#")],
              dtype=np.float64)
                          

        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                sample_ids.append(parts[0])
                c_res.append([float(x) for x in parts[1:]])

    c_res = np.array(c_res, dtype=np.float64)

    return header, c2, c_res, sample_ids


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


def run_gwas(runner, intermediate_file, TGWAS_file, snps_per_chunk=1000, device='cuda',  compress=False, regress_genotypes=True):
    """
    Run GWAS using a pre-configured GEMRunner instance.
    """
    # dir_name = os.path.dirname(out_file)
    # base_name = os.path.basename(out_file)
    # intermediate_file = os.path.join(dir_name, "intermediate_" + base_name)
    ph_headers, c2_values, corrected_res, resid_sample_ids = read_correction_file(intermediate_file) # read corrected_res, c2, ph_headers and sample ids from intermediate file
    
    if device == 'cuda' and torch.cuda.is_available():
        device = torch.device('cuda')
    else:
        device = torch.device('cpu')

    # corrected_res = corrected_res.to(device)
    corrected_res = torch.from_numpy(corrected_res).float().to(device)
    #covariates = covariates.to(device)
    
    n_samples, n_corrected_res = corrected_res.shape
    # Validate sample ids length matches residual rows
    if len(resid_sample_ids) != n_samples:
        print(f"Warning: intermediate file sample ID count ({len(resid_sample_ids)}) != residual rows ({n_samples}).")

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

    # If requested, prepare covariate projection to regress covariates out of genotypes (Check 1: before loop for each batches of SNPs)
    # cov_X = None
    start_readcov = time.time()
    cov_X = None
    proj_A = None
    if regress_genotypes:
        try:
            # applied_covs = None
            # # Primary: call runner-provided accessor if available
            # if hasattr(runner, "get_covariates"):
            #     cov_np = np.asarray(runner.get_covariates())  # shape (n_samples, n_covariates)
            #     applied_covs = None
            # else:
            # Fallback: read covariate file directly from runner options, or from python input.
            cov_file = runner.opt.cov_add
            cov_names = list(runner.opt.covariates) if hasattr(runner.opt, 'covariates') else []
            if not cov_names:
                # allow using all numeric columns if no explicit covariate names provided
                cov_names = []
            # choose delimiter from runner options if present, otherwise default to tab
            sep = getattr(runner.opt, 'cov_delim', '\t') or '\t'
            # Read with pandas and select requested covariate columns; assume sample order matches genotype order
            cov_df = pd.read_csv(cov_file, sep=sep)

            # Hard-coded covariate format: first column = family ID, second column = sample ID.
            # We will match residual/sample IDs using the sample ID (second) column and drop both
            # ID columns before selecting numeric covariates.
            cols = list(cov_df.columns)
            if len(cols) < 2:
                raise RuntimeError("Covariate file must have at least two columns: family ID and sample ID.")

            # Use the column (sample ID) for matching to intermediate residual sample IDs
            if hasattr(runner.opt, 'sampleid_header_name') and runner.opt.sampleid_header_name:
                sample_id_col = runner.opt.sampleid_header_name
                if sample_id_col not in cov_df.columns:
                    raise RuntimeError(f"Sample ID column '{sample_id_col}' not found in covariate file.")
                sample_ids_from_cov = cov_df[sample_id_col].astype(str).values
            else:
                raise RuntimeError(
                    "You must specify --sampleid-name for covariate file. "
                    "Default behavior of using the second column is disabled."
                )

            # Drop the two ID columns (FID and IID) to leave only covariate columns
            # cov_df2 = cov_df.drop(columns=[cols[0], cols[1]])

            # Select covariate columns
            if cov_names:
                try:
                    cov_sel = cov_df[cov_names]
                    applied_covs = list(cov_sel.columns)
                except Exception:
                    # fallback: take numeric columns
                    cov_sel = cov_df.select_dtypes(include=[np.number])
                    applied_covs = list(cov_sel.columns)
            else:
                cov_sel = cov_df.select_dtypes(include=[np.number])
                applied_covs = list(cov_sel.columns)

            # Reorder covariates to match residual/sample ID order from intermediate file (Check 2: match the Sample ID)
            try:
                # resid_sample_ids is read from the intermediate file earlier
                cov_sel = cov_sel.copy()
                cov_sel['_sample_id_for_match'] = sample_ids_from_cov
                cov_sel.set_index('_sample_id_for_match', inplace=True)
                # Reindex to the residual sample ID order; this will introduce NaN for missing rows
                cov_sel = cov_sel.reindex(resid_sample_ids)
                # If any missing after reindex, fail early
                if cov_sel.isnull().values.any():
                    missing = cov_sel.isnull().any(axis=1)
                    n_missing = int(missing.sum())
                    raise RuntimeError(f"Covariate file does not contain values for {n_missing} residual samples (after reindex).")
                # drop index and continue
                cov_sel.reset_index(drop=True, inplace=True)
            except Exception as e:
                raise

            cov_np = cov_sel.astype(float).values

            # Validate shape
            if cov_np.shape[0] != n_samples:
                # try transpose if user provided (n_covariates, n_samples)
                if cov_np.shape[1] == n_samples:
                    cov_np = cov_np.T
                else:
                    raise RuntimeError(f"Covariate shape mismatch: expected {n_samples} samples, got {cov_np.shape}")

            # Build design matrix with intercept
            intercept = np.ones((n_samples, 1), dtype=cov_np.dtype)
            cov_X_n = np.hstack([intercept, cov_np])  # shape (n_samples, p+1)
            # Compute projection coefficients matrix A = (X^T X)^{-1} X^T
            # XtX = cov_X.T @ cov_X
            # # use pseudo-inverse for numerical stability
            # inv_XtX = np.linalg.inv(XtX)
            # proj_A = inv_XtX @ cov_X.T  # shape (p+1, n_samples)
            ########################-Move to GPU-########################
            cov_X = torch.from_numpy(cov_X_n).to(device).float()
            print("cov_X device:", cov_X.device)

            # Compute XᵀX
            XtX = cov_X.T @ cov_X
            print("XtX device:", XtX.device)

            # Invert XᵀX
            inv_XtX = torch.linalg.inv(XtX)
            print("inv_XtX device:", inv_XtX.device)

            # Compute projection A
            proj_A = inv_XtX @ cov_X.T
            print("proj_A device:", proj_A.device)

            # Inform what covariates are applied (if known)
            if applied_covs is None:
                print("Info: genotype residualization: covariates obtained from runner.get_covariates()", flush=True)
            elif len(applied_covs) == 0:
                print("Info: genotype residualization: intercept only (no covariates)", flush=True)
            else:
                print(f"Info: genotype residualization will use covariates: {applied_covs}", flush=True)

        except Exception as e:
            print(f"Warning: failed to prepare covariate projection: {e}", flush=True)
            cov_X_n = None
            proj_A = None

    end_cov = time.time()
    print(f"time for reading covariate file and preparing projection = {end_cov - start_readcov:.2f}s")
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

    if os.path.exists(TGWAS_file):
        os.remove(TGWAS_file)

    start = time.time()

    for chunk_data, meta in tqdm(queue, desc="Processing SNPs"):
        actual_snps = chunk_data.shape[0]
        rows_in_buffer += actual_snps

        geno = geno_tensor[:actual_snps, :]
        beta = beta_tensor[:actual_snps, :]
        gamma = gamma_tensor[:actual_snps, :]

        # Regress covariates out of genotypes (Check 3: should operate on GPU)
        if proj_A is not None and chunk_data is not None:
            try:
                G = torch.from_numpy(chunk_data).float().to(device)  # shape (M, n_samples)
                print("G device:", G.device)
                coeffs = proj_A @ G.T
                print("coeffs device:", coeffs.device)
                fitted = cov_X @ coeffs
                G_resid = G - fitted.T
                geno.copy_(G_resid.float())
                print("geno device:", geno.device)

            except Exception as e:
                print(f"Warning: failed to regress covariates: {e}", flush=True)

        
        # if proj_A is not None and chunk_data is not None:
        #     try:
        #         G = np.asarray(chunk_data, dtype=np.float64)  # shape (M, n_samples)
        #         # coeffs: (p+1, M) = proj_A (p+1 x n_samples) @ G.T (n_samples x M)
        #         coeffs = proj_A @ G.T
        #         # fitted: (n_samples, M) = cov_X (n_samples x p+1) @ coeffs (p+1 x M)
        #         fitted = cov_X @ coeffs
        #         G_resid = G - fitted.T  # back to (M, n_samples)
        #         chunk_data = G_resid.astype(np.float32)
        #     except Exception as e:
        #         # fallback to original data on failure
        #         print(f"Warning: failed to regress covariates from genotypes for this chunk: {e}", flush=True)
        # geno.copy_(torch.from_numpy(chunk_data).float())
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
                TGWAS_file, combined.schema, compression="snappy"
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
                TGWAS_file, combined.schema, compression="snappy"
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
    print(f"time for calculating GWAS = {end - start:.2f}s")
