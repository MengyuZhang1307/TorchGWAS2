#to run on server $HOME/local/python3.11/bin/python3 Bgen.py or python3 Bgen.py 

import sys
#Set path
sys.path.append("build")  
import csv
import Mygen
from Mygen import GEMOptions, run_bgen

opt = GEMOptions()
opt.pheno_file = "example/example.pheno2-2id"
opt.cov_file = "example/example.cov-2id"
opt.delim_pheno = ','  # 3rd argument
opt.delim_cov = ','    # 4th argument
opt.bgen_file = "example/example.bgen"
opt.sample_file = "example/example.sample"
opt.do_filters = False
opt.use_sample_file = True
opt.includeVariantFile = ""
opt.stream_snps = 1
opt.sampleid_header_name = "sampleid"
opt.random_slope_header_name = ""
opt.covariates = ["time", "cov3"]
opt.exposures = ["cov1"]
opt.interactions = ["time"]
opt.missing_key = "NA"
opt.threads = 2
opt.out_file = "outpy.txt"

results = run_bgen(opt)


results = run_bgen(opt)

with open("out_dosage.csv", "w", newline="") as f:
    writer = csv.writer(f)
    # Level 1: outermost (e.g., chromosomes or batches)
    for level1 in results:
        # Level 2: individual SNPs
        for level2 in level1:
            # Level 3: sample dosages
            for level3 in level2:
                writer.writerow(level3)  # Level 4: a list of floats
