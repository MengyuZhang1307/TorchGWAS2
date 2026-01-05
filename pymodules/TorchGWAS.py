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


def remove_collinear_columns_np(X: np.ndarray, col_names=None):
    """
    Remove (near-)collinear columns using QR diag threshold.

    Parameters
    ----------
    X : (n, p) np.ndarray
        Design matrix.
    col_names : list[str] | None
        Optional column names length p.
    keep_first : bool
        If True, never drop column 0 (useful for intercept).

    Returns
    -------
    X_new : np.ndarray
    """
    print("Checking collinearity for the regressed covariates...\n")
    X = np.asarray(X)
    n, p = X.shape

    # QR decomposition (like HouseholderQR)
    # R is shape (min(n,p), p) in 'reduced' mode.
    _, R = np.linalg.qr(X, mode="reduced")

    # Make diagR length = p (pad zeros if p > n)
    diagR = np.zeros(p, dtype=float)
    d = min(n, p)
    if d > 0:
        diagR[:d] = np.abs(np.diag(R[:d, :d]))

    sqrtEps = np.sqrt(np.finfo(X.dtype).eps)
    maxdiag = diagR.max() if p > 0 else 0.0
    cutoff = maxdiag * sqrtEps

    dropped_idx = [j for j in range(p) if diagR[j] < cutoff]

    if dropped_idx:
        dropped_names = []
        for j in dropped_idx:
            if j == 0:
                # intercept got flagged 
                dropped_names.append("intercept")
            else:
                jj = j - 1
                if 0 <= jj < len(col_names):
                    dropped_names.append(col_names[jj])
        print(f"Dropped columns: {dropped_names}")
    else:
        print("Dropped columns: []")

    dropped_set = set(dropped_idx)
    keep_idx = [j for j in range(p) if j not in dropped_set]

    X_new = X[:, keep_idx]

    if col_names is not None:
        new_col_names = []
        for j in keep_idx:
            if j == 0:
                # skip intercept:
                continue
            else:
                jj = j - 1
                if 0 <= jj < len(col_names):
                    new_col_names.append(col_names[jj])

    return X_new, new_col_names


def has_duplicates(resid_sample_ids):
    return len(resid_sample_ids) != len(set(resid_sample_ids))


def fill_J(resid_sample_ids, cov_sample_ids, device="cpu", dtype=torch.float32):
    
    """
    resid_sample_ids: observation IDs in covariate order (after filtering)
    sample_id_for_G:  unique IDs in the SAME order as G columns
    """
    id2col = {sid: j for j, sid in enumerate(resid_sample_ids)}  # required for alignment to G
    n_obs = len(cov_sample_ids)
    n_unique = len(resid_sample_ids)

    col_idx = torch.tensor([id2col[sid] for sid in cov_sample_ids],
                           device=device, dtype=torch.long)
    row_idx = torch.arange(n_obs, device=device, dtype=torch.long)

    indices = torch.stack([row_idx, col_idx], dim=0)
    values  = torch.ones(n_obs, device=device, dtype=dtype)

    J = torch.sparse_coo_tensor(indices, values, (n_obs, n_unique)).coalesce()
    return J



def calc_cov_proj(runner_opt, resid_sample_ids, device):
    """
    geno_ids: unique IDs in the SAME order as G columns  (n_unique)
    obs_ids: observation IDs in covariate order (after filtering)
    Returns:
      proj_A: (p x n_unique) if no-dup case, else None
      cov_X:  (n_unique x p) if no-dup case, else (n_obs x p)
      has_dup: bool (dup IDs in cov file after filtering to genotype IDs)
      J: (n_obs x n_unique) if dup case, else None
      obs_ids: list[str]  (IDs in cov_X row order)
    """
    # Treat input resid_sample_ids as genotype sample IDs in G column order (unique)
    geno_id_set = set(resid_sample_ids)

    cov_file = runner_opt.cov_add
    cov_names = list(runner_opt.covariates) if hasattr(runner_opt, "covariates") else []
    sep = getattr(runner_opt, "cov_delim", "\t") or "\t"

    cov_df = pd.read_csv(cov_file, sep=sep)

    if not (hasattr(runner_opt, "sampleid_header_name") and runner_opt.sampleid_header_name):
        raise RuntimeError("You must specify --sampleid-name for covariate file.")
    sample_id_col = runner_opt.sampleid_header_name
    if sample_id_col not in cov_df.columns:
        raise RuntimeError(f"Sample ID column '{sample_id_col}' not found in covariate file.")

    sample_ids_from_cov = cov_df[sample_id_col].astype(str).to_numpy()

    if cov_names:
        cov_sel = cov_df[cov_names].copy()
    else:
        cov_sel = cov_df.select_dtypes(include=[np.number]).copy()

    # ---- filter cov rows to genotype IDs, KEEPING COV FILE ORDER ----
    keep_mask = pd.Series(sample_ids_from_cov).isin(geno_id_set).to_numpy()
    cov_sel = cov_sel.iloc[keep_mask].reset_index(drop=True)
    obs_ids = sample_ids_from_cov[keep_mask].tolist()

    has_dup = has_duplicates(obs_ids)

    # Build design matrix (in current obs order)
    cov_np = cov_sel.astype(float).to_numpy()
    intercept = np.ones((cov_np.shape[0], 1), dtype=cov_np.dtype)
    cov_X_n = np.hstack([intercept, cov_np])  # (n_rows, p)
    cov_X_n, new_cov_nam = remove_collinear_columns_np(cov_X_n, cov_names)
    # move to torch
    cov_X = torch.as_tensor(cov_X_n, device=device, dtype=torch.float32)

    if not has_dup:
        # reorder cov_X rows to match genotype order (resid_sample_ids)
        row_map = {sid: i for i, sid in enumerate(obs_ids)}  # unique now
        row_idx = [row_map[sid] for sid in resid_sample_ids]  # all should exist after filtering
        cov_X = cov_X[row_idx]  # (n_unique, p)

        # projection A = (X^T X)^-1 X^T  (more stable than inv: use solve)
        XtX = cov_X.T @ cov_X
        # # Invert XᵀX
        # inv_XtX = torch.linalg.inv(XtX)
        # # Compute projection A
        # proj_A = inv_XtX @ cov_X.T
        proj_A = torch.linalg.solve(XtX, cov_X.T)  # (p x n_unique)

        J = None
        return proj_A, cov_X, J, has_dup

    else:
        # duplicated obs IDs: keep cov order and build J to map obs->genotype columns
        J = fill_J(resid_sample_ids, obs_ids, device=device, dtype=torch.float32)
        proj_A = None
        return proj_A, cov_X, J, has_dup



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
    cov_X, proj_A, J, has_dup = None, None, None, None
    start_readcov = time.time()
    if regress_genotypes:
        proj_A, cov_X, J, has_dup = calc_cov_proj(
        runner.opt,
        resid_sample_ids,
        device)
    end_readcov = time.time()
    print(f"Time for reading covariate file and preparing projection = {end_readcov - start_readcov:.2f}s")
 
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

        if chunk_data is None:
            raise ValueError("chunk_data is None")

        G = torch.as_tensor(chunk_data, dtype=torch.float32, device=device)  # (M, n)

        if has_dup:
            J = J.coalesce()
            Jt = J.transpose(0, 1).coalesce()              # (n_uniq, n_obs)
            JT_X = torch.sparse.mm(Jt, cov_X)       # only supports(sp*dens) (n_uniq, p)
            S = JT_X.T                                     # (p, n_uniq)

            XTX = cov_X.T @ cov_X                          # (p, p)
            XTX_i_S = torch.linalg.solve(XTX, S)           # (p, n_uniq)
            counts = torch.sparse.sum(J, dim=0).to_dense() # (n_uniq,)
            GD = G * counts.unsqueeze(0)                       # (M, n_uniq)
            tmp = G @ JT_X                                 # (M, p)
            corr = tmp @ XTX_i_S                           # (M, n_uniq)

            geno.copy_(GD - corr) / counts.unsqueeze(0)                         # (M, n_uniq)

        else:
            coeffs = proj_A @ G.T
            fitted = cov_X @ coeffs
            geno.copy_(G - fitted.T)
        
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