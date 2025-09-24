#!/usr/bin/env python3

import sys
import os
import numpy as np
import torch
from tqdm import tqdm

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
        - header: list of strings (column names)
        - c2: numpy array of double (values from 2nd line after '#')
        - data_array: numpy 2D array of double (phenotype values)
    """
    c2 = np.empty(0, dtype=np.float64)
    c_res = []
    
    with open(file_path, "r") as f:
        header = f.readline().strip().split("\t")
        second_line = f.readline().strip().split("\t")
        if second_line[0].startswith("#"):
            c2 = [float(val) for val in second_line if not val.startswith("#")]
                          

        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                c_res.append([float(x) for x in parts[2:]])

    c_res = np.array(c_res, dtype=np.float64)

    return header, c2, c_res


def calc_t(pheno_normalized, geno, beta, gamma, sqrt_c2):
    """
    Compute per-SNP stats using pre-normalized phenotypes and on-device sqrt(c2).

    pheno_normalized: (N, P)
    geno: (M, N)
    beta, gamma: (M, P) work buffers on same device as geno
    sqrt_c2_device: (P,) tensor on same device as inputs
    """
    N = pheno_normalized.shape[0]
    with torch.no_grad():
        geno.sub_(geno.mean(1, keepdim=True)).div_(geno.std(1, keepdim=True))
        torch.matmul(geno, pheno_normalized, out=beta)
        beta.div_(N)
        gamma.copy_(beta)
        gamma.pow_(2).sub_(1).div_(2-N)
        torch.sqrt(gamma, out=gamma)
        gamma.div_(sqrt_c2.unsqueeze(0))
        beta_coeffs = beta.cpu()
        se = gamma.cpu()
        t_stats = beta.div_(gamma).abs_().neg_().cpu()
        return t_stats, beta_coeffs, se


def run_gwas(runner, snps_per_chunk=1000, device='cuda'):
    """
    Run GWAS using a pre-configured GEMRunner instance.
    """
    # Get C2 values from the runner (try fitting null model first)
    c2_values = None
    try:
        runner.run_fit_nullmodel()
        c2_values = runner.get_c2_values()
        print(f"Fitted null model and obtained C2 values for SE adjustment: {c2_values}")
    except Exception as e:
        print(f"Could not fit null model or get C2 values: {e}")
        print("Proceeding without C2 adjustment")
    if c2_values is None:
        return
    
    # Get phenotypes and covariates from runner
    phenotypes = torch.from_numpy(runner.get_phenotypes()).float()
    covariates = torch.from_numpy(runner.get_covariates()).float()
    
    if device == 'cuda' and torch.cuda.is_available():
        device = torch.device('cuda')
    else:
        device = torch.device('cpu')
    
    phenotypes = phenotypes.to(device)
    covariates = covariates.to(device)
    
    n_samples, n_phenotypes = phenotypes.shape
    
    # Center phenotypes with NaN-safe mean and replace NaNs
    ph_mean = torch.nanmean(phenotypes, dim=0, keepdim=True)
    phenotypes = phenotypes - ph_mean
    phenotypes = torch.nan_to_num(phenotypes, nan=0.0)
    
    c = covariates.cpu().numpy()
    c_mean = np.nanmean(c, axis=0, keepdims=True)
    c_std = np.nanstd(c, axis=0, keepdims=True)
    c_std[c_std < 1e-12] = 1.0
    c = (c - c_mean) / c_std
    c = np.nan_to_num(c, nan=0.0)
    covarQ, _ = np.linalg.qr(c)
    covarQ = torch.from_numpy(covarQ).float().to(device)
    
    pheno_normalized = phenotypes - torch.matmul(covarQ, torch.matmul(covarQ.T, phenotypes))
    pheno_normalized = pheno_normalized / torch.std(pheno_normalized, dim=0, keepdim=True)
    pheno_normalized = torch.nan_to_num(pheno_normalized, nan=0.0)
    
    # Precompute sqrt(c2) once on the target device (clamped for stability)
    sqrt_c2 = torch.from_numpy(np.asarray(c2_values)).to(device=device, dtype=torch.float32)
    sqrt_c2.sqrt_()

    # Start dosage streaming from runner
    queue = runner.start_dosage_stream(queue_capacity=10, snps_per_chunk=snps_per_chunk)
    
    all_t_stats = []
    all_beta = []
    all_se = []
    
    # Preallocate device buffers and reuse/slice for smaller final chunks
    geno_tensor = torch.empty(snps_per_chunk, n_samples, device=device)
    beta_tensor = torch.empty(snps_per_chunk, n_phenotypes, device=device)
    gamma_tensor = torch.empty(snps_per_chunk, n_phenotypes, device=device)
    
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
        
        geno.copy_(chunk_data)
        t_stats, beta_coeffs, se = calc_t(pheno_normalized, geno, beta, gamma, sqrt_c2)
        
        all_t_stats.append(t_stats)
        all_beta.append(beta_coeffs)
        all_se.append(se)
    
    t_statistics = torch.cat(all_t_stats, dim=0)
    beta_coefficients = torch.cat(all_beta, dim=0)
    standard_errors = torch.cat(all_se, dim=0)
    p_values = 2 * torch.special.ndtr(t_statistics)
    
    return {
        't_stats': t_statistics,
        'beta': beta_coefficients,
        'se': standard_errors,
        'p_values': p_values,
    }
