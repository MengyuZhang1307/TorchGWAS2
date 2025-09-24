#to run on server $HOME/local/python3.11/bin/python3 Bgen.py or python3 Bgen.py 

import sys
#Set path
sys.path.append("build")  
sys.path.append("pymodules")  
import Mygen
from pymodules import ConfOpt
# from torchgwas import run_gwas
from torchgwas import read_correction_file
from Mygen import GEMRunner
import numpy as np
import time


# function to read corrected residuals and C2
import numpy as np


start_time = time.time()

opt = ConfOpt(pheno_file = "example/example.pheno2-2id",
            cov_file = "example/example.cov-2id",
            delim_pheno = ',',
            delim_cov = ',',
            geno_file = "example/example.bgen",
            sample_file = "example/example.sample",
            do_filters = False,
            use_sample_file = True,
            includeVariantFile = "",
            stream_snps = 1,
            sampleid_header_name = "sampleid",
            random_slope_header_name = "",
            covariates = ["cov3"],
            exposures = ["cov1"],
            interactions = [],
            missing_key = "NA",
            threads = 5, 
            num_chunks = 5,
            outfile = "outpy.txt")

# opt = ConfOpt(
#     pheno_file = "/HGCNT95FS/ADDLIE/Sama-GEM2/TORCH/T2_pheno_QT_repeated",
#     cov_file = "/HGCNT95FS/ADDLIE/Sama-GEM2/TORCH/T2_covar",
#     delim_pheno = "\t",
#     delim_cov = " ",
#     geno_file = "/HGCNT95FS/ADDLIE/Sama-GEM2/all_filtered.bgen",
#     sample_file = "/HGCNT95FS/ADDLIE/Sama-GEM2/MRI_samples_chr1.sample",
#     do_filters = False,
#     use_sample_file = True,
#     includeVariantFile = "",
#     stream_snps = 1000,
#     sampleid_header_name = "IID",
#     random_slope_header_name = "PC2",
#     covariates = ["PC1"],
#     exposures = ["SEX"],
#     interactions = [],
#     missing_key = "NA",
#     kin_path = "/HGCNT95FS/ADDLIE/Sama-GEM2/kinship.txt",
#     delim_k = ' ',
#     kin_diag = 0.5,
#     threads = 72,
#     num_chunks = 72,
#     outfile = "outAddlie.txt"
# )

runner = GEMRunner(opt.get())

# Run null model fitting
runner.run_fit_nullmodel()


print("Starting streaming dosage decode...")
# Capacity controls backpressure and is also used by C++ to choose thread count
queue_capacity = max(1, int(opt.threads))
snps_per_chunk = max(1, int(opt.stream_snps))
q = runner.start_dosage_stream(queue_capacity, snps_per_chunk)

total_rows = 0
first_chunk_shape = None
for chunk in q:
    if first_chunk_shape is None:
        first_chunk_shape = chunk.shape
        print("First chunk shape:", first_chunk_shape)
        print("First row sample:", chunk[0:1])
    total_rows += chunk.shape[0]
print("Total SNP rows streamed:", total_rows)
end_time = time.time()
print("End of reading Bgen file\n")
print("Wall time in seconds :", end_time - start_time)

header, c2, c_res = read_correction_file("outpy.txt")

print(c_res.shape)    
print(c_res.ndim)     
print(c_res.size)

