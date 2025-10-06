#to run on server $HOME/local/python3.11/bin/python3 Bgen.py or python3 Bgen.py 

import sys
import pandas as pd
#Set path
sys.path.append("build")  
sys.path.append("pymodules")  
import Mygen
from pymodules import ConfOpt
# from torchgwas import run_gwas
from torchgwas import read_correction_file, run_gwas
from Mygen import GEMRunner
import numpy as np
import time
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor

def report_array_info(name, arr):
    print(f"{name}:")
    print(f"  shape = {arr.shape}")
    print(f"  elements = {arr.size:,}")
    print(f"  memory ≈ {arr.nbytes / 1e6:.2f} MB\n")
#save results to files
def save_pheno(i, pheno_name, betas, ses, tstats, pvals):
    ph_df = pd.DataFrame({
        "beta": betas[:, i],
        "se": ses[:, i],
        "t_stat": tstats[:, i],
        "p_value": pvals[:, i],
    })
    ph_df.to_csv(f"rersults_{pheno_name}.csv", sep="\t", index=False)
    # df.to_parquet(f"results_{pheno_name}.parquet", engine="pyarrow", compression="zstd")

#if we have number of files as threads
def save_block(start, end, headers, betas, ses, tstats, pvals, block_id):
    # Build one DataFrame for this block of phenotypes
    cols = {}
    for i in range(start, end):
        ph = headers[i]
        cols[f"beta_{ph}"]  = betas[:, i]
        cols[f"se_{ph}"]    = ses[:, i]
        cols[f"t_stat_{ph}"] = tstats[:, i]
        cols[f"pval_{ph}"]  = pvals[:, i]

    ph_df = pd.DataFrame(cols)
    ph_df.to_csv(f"rersults_{block_id}.csv", sep="\t", index=False)
    # df.to_parquet(f"results_block_{block_id}.parquet", engine="pyarrow", compression="zstd")

start_time = time.time()

# opt = ConfOpt(pheno_add = "example/example.pheno-2id-repeated",
#             cov_add = "example/example.pheno",
#             pheno_delim = ',',
#             cov_delim = ',',
#             kin_add = "example/example.kinship",
#             kin_delim = ',',
#             kin_diag = 0.5,
#             geno_add = "example/example.bgen",
#             sample_add = "example/example.sample",
#             do_filters = False,
#             use_sample_file = True,
#             includeVariantFile = "",
#             stream_snps = 1,
#             sampleid_header_name = "sampleid",
#             random_slope_header_name = "",
#             covariates = ["cov3"],
#             exposures = ["cov1"],
#             interactions = [],
#             missing_key = "NA",
#             threads = 5, 
#             num_chunks = 5,
#             outfile = "outexample.txt")

# opt = ConfOpt(pheno_file = "example/example.pheno2-2id",
#             cov_file = "example/example.cov-2id",
#             delim_pheno = ',',
#             delim_cov = ',',
#             geno_file = "example/example.bgen",
#             sample_file = "example/example.sample",
#             do_filters = False,
#             use_sample_file = True,
#             includeVariantFile = "",
#             stream_snps = 1,
#             sampleid_header_name = "sampleid",
#             random_slope_header_name = "",
#             covariates = ["cov3"],
#             exposures = ["cov1"],
#             interactions = [],
#             missing_key = "NA",
#             threads = 5, 
#             num_chunks = 5,
#             outfile = "outexample.txt")

opt = ConfOpt(
    pheno_add = "data/T2_pheno_QT_repeated",
    cov_add = "data/T2_covar",
    pheno_delim = "\t",
    cov_delim = " ",
    geno_add = "data/all_filtered.bgen",
    sample_ = "data/MRI_samples_chr1.sample",
    do_filters = False,
    use_sample_file = True,
    includeVariantFile = "",
    stream_snps = 10000,
    sampleid_header_name = "IID",
    random_slope_header_name = "",
    covariates = ["PC1"],
    exposures = [],
    interactions = [],
    missing_key = "NA",
    kin_add = "data/kinship.txt",
    kin_delim = ' ',
    kin_diag = 0.5,
    threads = 90,
    num_chunks = 90,
    outfile = "outAddlie.txt"
)

runner = GEMRunner(opt.get())

# Run null model fitting
runner.run_fit_nullmodel()


print("Starting streaming dosage decode...")
run_gwas(runner, snps_per_chunk=1000, device='cuda')
end_time_gwas = time.time()
print("End of running Torch GWAS file\n")
print("Wall time in seconds :", end_time_gwas - start_time)
##calculate data type change
# start_dconv_time = time.time()
# t_stats = results['t_stats'].cpu().numpy()
# betas = results['beta'].cpu().numpy()
# ses = results['se'].cpu().numpy()
# pvals = results['p_values'].cpu().numpy()
# ph_headers = results['ph_headers'] [2:]  # <-- phenotype names
# end_dconv_time = time.time()
# print("End of data type conversion file\n")
# print("Wall time in seconds :", end_dconv_time - start_dconv_time)
# ## calculate size
# report_array_info("t_stats", t_stats)
# report_array_info("betas", betas)
# report_array_info("ses", ses)
# report_array_info("pvals", pvals)
# print("phenos name:", ph_headers)

# start_writing_time = time.time()

# # Example: 1280 phenotypes
# with ThreadPoolExecutor(max_workers=96) as executor:  # adjust threads
#     for i, pheno_name in enumerate(ph_headers):
#         executor.submit(save_pheno, i, pheno_name, betas, ses, t_stats, pvals)


# end_writing_time = time.time()
# print("End of writing 2 file\n")
# print("Wall time in seconds :", end_writing_time - start_writing_time)




# Divide columns across threads
# start_writing_time = time.time()

# n_pheno = len(ph_headers)
# n_threads = 96
# n_threads = min(n_threads, n_pheno)
# block_size = (n_pheno + n_threads - 1) // n_threads  # ceil division

# with ThreadPoolExecutor(max_workers=n_threads) as executor:
#     for block_id in range(n_threads):
#         start = block_id * block_size
#         end = min((block_id + 1) * block_size, n_pheno)
#         if start < end:  # only submit if block not empty
#             executor.submit(save_block, start, end, ph_headers, betas, ses, t_stats, pvals, block_id)

# with ProcessPoolExecutor(max_workers=48) as executor:  # match physical cores
#     for block_id in range(n_threads):
#         start = block_id * block_size
#         end = min((block_id + 1) * block_size, n_pheno)
#         if start < end:
#             executor.submit(save_block, start, end, ph_headers, betas, ses, t_stats, pvals, block_id)
# end_writing_time = time.time()

# print("End of writing 2 file\n")
# print("Wall time in seconds :", end_writing_time - start_writing_time)


