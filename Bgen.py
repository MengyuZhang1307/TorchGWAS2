#to run on server $HOME/local/python3.11/bin/python3 Bgen.py or python3 Bgen.py 

import sys
#Set path
sys.path.append("build-py")  
sys.path.append("pymodules")  
import Mygen
from pymodules import ConfOpt
from Mygen import GEMRunner
import numpy as np
import time
start_time = time.time()

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
#             threads = 2, 
#             num_chunks = 5,
#             outfile = "outpy.txt")

opt = ConfOpt(
    pheno_file = "/HGCNT95FS/ADDLIE/Sama-GEM2/TORCH/T2_pheno_QT_repeated",
    cov_file = "/HGCNT95FS/ADDLIE/Sama-GEM2/TORCH/T2_covar",
    delim_pheno = "\t",
    delim_cov = " ",
    geno_file = "/HGCNT95FS/ADDLIE/Sama-GEM2/all_filtered.bgen",
    sample_file = "/HGCNT95FS/ADDLIE/Sama-GEM2/MRI_samples_chr1.sample",
    do_filters = False,
    use_sample_file = True,
    includeVariantFile = "",
    stream_snps = 1000,
    sampleid_header_name = "IID",
    random_slope_header_name = "PC2",
    covariates = ["PC1"],
    exposures = ["SEX"],
    interactions = [],
    missing_key = "NA",
    kin_path = "/HGCNT95FS/ADDLIE/Sama-GEM2/kinship.txt",
    delim_k = ' ',
    kin_diag = 0.5,
    threads = 72,
    num_chunks = 32,
    outfile = "outAddlie.txt"
)

runner = GEMRunner(opt.get())

# NEW: Access raw phenotype and covariate data in BGEN sample order
try:
    print("Accessing phenotype and covariate data in BGEN sample order...")
    
    # Get data as NumPy arrays - all in BGEN sample order
    sample_ids = runner.get_sample_ids()           # BGEN sample order
    phenotypes = runner.get_phenotypes()           # BGEN sample order  
    covariates = runner.get_covariates()           # BGEN sample order
    
    # Verify alignment
    assert len(sample_ids) == phenotypes.shape[0] == covariates.shape[0], "Sample count mismatch!"
except RuntimeError as e:
    print(f"Could not access data: {e}")
except AssertionError as e:
    print(f"Sample alignment error: {e}")

# Run null model fitting
runner.run_fit_nullmodel()
# opt = GEMOptions()
# opt.pheno_file = "example/example.pheno2-2id"
# opt.cov_file = "example/example.cov-2id"
# opt.delim_pheno = ','  # 3rd argument
# opt.delim_cov = ','    # 4th argument
# opt.bgen_file = "example/example.bgen"
# opt.sample_file = "example/example.sample"
# opt.do_filters = False
# opt.use_sample_file = True
# opt.includeVariantFile = ""
# opt.stream_snps = 1
# opt.sampleid_header_name = "sampleid"
# opt.random_slope_header_name = ""
# opt.covariates = ["time", "cov3"]
# opt.exposures = ["cov1"]
# opt.interactions = ["time"]
# opt.missing_key = "NA"
# opt.threads = 2
# opt.outfile = "outpy.txt"


# opt.pheno_file = "/HGCNT95FS/ADDLIE/Sama-GEM2/TORCH/T2_pheno"
# opt.cov_file = "/HGCNT95FS/ADDLIE/Sama-GEM2/TORCH/T2_covar"
# opt.delim_pheno = " "  # 3rd argument
# opt.delim_cov = " "   # 4th argument
# opt.bgen_file = "/HGCNT95FS/ADDLIE/Sama-GEM2/all_filtered.bgen"
# opt.sample_file = "/HGCNT95FS/ADDLIE/Sama-GEM2/MRI_samples_chr1.sample"
# opt.do_filters = False
# opt.use_sample_file = True
# opt.includeVariantFile = ""
# opt.stream_snps = 1
# opt.sampleid_header_name = "IID"
# opt.random_slope_header_name = ""
# opt.covariates = ["PC2"]
# opt.exposures = ["PC1"]
# opt.interactions = []
# opt.missing_key = "NA"
# opt.threads = 72
# opt.outfile = "outAddlie.txt"



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




