#to run on server $HOME/local/python3.11/bin/python3 Bgen.py or python3 Bgen.py 

import sys
import pandas as pd
#Set path
sys.path.append("build-test")  
sys.path.append("pymodules")  
# import Mygen
from pymodules import ConfOpt
# from torchgwas import run_gwas
from torchgwas import read_correction_file, run_gwas
from Mygen import GEMRunner
import numpy as np
import time
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor

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
    sample_add = "data/MRI_samples_chr1.sample",
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

# opt = ConfOpt(
#     pheno_add = "missing-simulationbyMengyu/phenotype_missing.txt",
#     cov_add = "missing-simulationbyMengyu/cov_missing.txt",
#     pheno_delim = "\t",
#     cov_delim = "\t",
#     geno_add = "missing-simulationbyMengyu/SA.bgen",
#     sample_add = "missing-simulationbyMengyu/SA.sample",
#     do_filters = False,
#     use_sample_file = True,
#     includeVariantFile = "",
#     stream_snps = 10000,
#     sampleid_header_name = "id",
#     random_slope_header_name = "",
#     covariates = ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"],
#     exposures = [],
#     interactions = [],
#     missing_key = "NA",
#     kin_add = "missing-simulationbyMengyu/kinfile_extended_fam.txt",
#     kin_delim = '\t',
#     kin_diag = 1,
#     threads = 90,
#     num_chunks = 90,
#     outfile = "missing-x1-x10.txt"
# )


runner = GEMRunner(opt.get())

# Run null model fitting
runner.run_fit_nullmodel()


print("Starting streaming dosage decode...")
run_gwas(runner, snps_per_chunk=1000, device='cuda')
end_time_gwas = time.time()
print("End of running Torch GWAS file\n")
print("Wall time in seconds :", end_time_gwas - start_time)
