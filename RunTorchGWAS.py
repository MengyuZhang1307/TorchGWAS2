#to run on server $HOME/local/python3.11/bin/python3 Bgen.py or python3 Bgen.py 

import pandas as pd
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), "pymodules"))
from pymodules import ConfOpt
# import Mygen
from pymodules import GEMRunner
from pymodules import run_gwas
from pymodules import parquet_to_text_duckdb
import numpy as np
import time
import argparse

def normalize_delim(s):
    # Convert to single-character delimiter that C++ expects
    if s in ["\\t", "t", "tab", "TAB", r"\t"]:
        return "\t"
    if s in ["\\0", "0", "space", " "]:
        return " "
    if s == ",":
        return ","
    return s

def parse_args():
    parser = argparse.ArgumentParser(description="Run TorchGWAS using GEM2 and Torch backend.")
    parser.add_argument("--pheno-file", type=str, required=True, help="Phenotype file path")
    parser.add_argument("--cov-file", type=str, required=True, help="Covariate file path")
    parser.add_argument("--bgen", type=str, required=True, help="genotype file path")
    parser.add_argument("--sample", type=str, required=True, help="BGEN Sample file path (optional)")
    parser.add_argument("--kin-file", type=str, default="", help="Kinship file path (optional)")
    parser.add_argument("--kin-diag", type=float, default=1.0, help="Diagonal value of " \
                        "kinship matrix that not accounting for inbreeding (Default: 1.0)")
    parser.add_argument("--pheno-delim", type=str, default=",", help="Phenotype file delimiter (default: comma)")
    parser.add_argument("--cov-delim", type=str, default=",", help="Covariate file delimiter (default: comma)")
    parser.add_argument("--kin-delim", type=str, default=",", help="Kinship file delimiter (default: comma)")
    parser.add_argument("--sampleid-name", type=str, help="sample ID herader name")
    parser.add_argument("--include-snp-file", type=str, default= "", help="Path to file containing a subset of variants in \
                        the specified genotype file to be used for analysis. The first line in this file is the header that specifies\
                        which variant identifier in the genotype file is used for ID matching. This must be 'snpid' (PLINK or BGEN)\
                        or 'rsid' (BGEN only). There should be one variantidentifier per line after the header.")
    parser.add_argument("--covar-names", nargs="+", help="Covariate names list")
    parser.add_argument("--random-slope-name", type=str, default = "", help="Column name in the covariate file that contains random slope (default: "").")
    parser.add_argument("--missing-value", type=str, default="NA", help="Indicates how missing values in the phenotype and covariate files are stored.")
    parser.add_argument("--threads", type=int, help="Number of threads")
    # parser.add_argument("--num-chunks", type=int, help="Number of chunks")
    parser.add_argument("--stream-snps", type=int, default=1000, help="Number of SNPs per chunk")
    parser.add_argument("--out", type=str, default="out.txt", help="Output file name")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cuda", help="Computation device (default: cuda)")
    parser.add_argument("--verbose", action="store_true", help="Print null model(default: False)")
    parser.add_argument("--convert", action="store_false", help="Convert binary to text file (default: True)")
    return parser.parse_args()

def main():
    args = parse_args()
    start_time = time.time()

    confopt = ConfOpt(
        pheno_add=args.pheno_file,
        pheno_delim=normalize_delim(args.pheno_delim),
        cov_add=args.cov_file,
        cov_delim=normalize_delim(args.cov_delim),
        geno_add=args.bgen,
        sample_add=args.sample,
        use_sample_file=bool(args.sample),
        do_filters=bool(args.include_snp_file),
        includeVariantFile=args.include_snp_file,
        stream_snps=args.stream_snps,
        sampleid_header_name= args.sampleid_name,
        covariates=args.covar_names,
        random_slope_header_name=args.random_slope_name,
        missing_key=args.missing_value,
        kin_add=args.kin_file,
        kin_delim=normalize_delim(args.kin_delim),
        kin_diag=args.kin_diag,
        threads=args.threads,
        # num_chunks=args.num_chunks,
        outfile=args.out
    )

    runner = GEMRunner(confopt.get())

    print("Running null model fitting ...")
    runner.run_fit_nullmodel()

    print("Starting dosage streaming and GWAS ...")
    run_gwas(runner, out_file=args.out, snps_per_chunk=args.stream_snps, device=args.device)

    end_time = time.time()
    print("\n TorchGWAS completed successfully.")
    print(f"Wall time: {(end_time - start_time):.2f} seconds")
    if args.convert:
        parquet_to_text_duckdb(input_file="TGWAS_" + args.out + ".parquet", output_file=args.out + "txt")

if __name__ == "__main__":
    main()



# start_time = time.time()

# # opt = ConfOpt(pheno_add = "example/example.pheno-2id-repeated",
# #             cov_add = "example/example.pheno",
# #             pheno_delim = ',',
# #             cov_delim = ',',
# #             kin_add = "example/example.kinship",
# #             kin_delim = ',',
# #             kin_diag = 0.5,
# #             geno_add = "example/example.bgen",
# #             sample_add = "example/example.sample",
# #             do_filters = False,
# #             use_sample_file = True,
# #             includeVariantFile = "",
# #             stream_snps = 1,
# #             sampleid_header_name = "sampleid",
# #             random_slope_header_name = "",
# #             covariates = ["cov3"],
# #             exposures = ["cov1"],
# #             interactions = [],
# #             missing_key = "NA",
# #             threads = 5, 
# #             num_chunks = 5,
# #             outfile = "outexample.txt")

# # opt = ConfOpt(pheno_file = "example/example.pheno2-2id",
# #             cov_file = "example/example.cov-2id",
# #             delim_pheno = ',',
# #             delim_cov = ',',
# #             geno_file = "example/example.bgen",
# #             sample_file = "example/example.sample",
# #             do_filters = False,
# #             use_sample_file = True,
# #             includeVariantFile = "",
# #             stream_snps = 1,
# #             sampleid_header_name = "sampleid",
# #             random_slope_header_name = "",
# #             covariates = ["cov3"],
# #             exposures = ["cov1"],
# #             interactions = [],
# #             missing_key = "NA",
# #             threads = 5, 
# #             num_chunks = 5,
# #             outfile = "outexample.txt")

# confopt = ConfOpt(
#     pheno_add = "data/T2_pheno_QT_repeated",
#     cov_add = "data/T2_covar",
#     pheno_delim = "\t",
#     cov_delim = " ",
#     geno_add = "data/all_filtered.bgen",
#     sample_add = "data/MRI_samples_chr1.sample",
#     do_filters = False,
#     use_sample_file = True,
#     includeVariantFile = "",
#     stream_snps = 1000,
#     sampleid_header_name = "IID",
#     random_slope_header_name = "",
#     covariates = ["PC1"],
#     exposures = [],
#     interactions = [],
#     missing_key = "NA",
#     kin_add = "data/kinship.txt",
#     kin_delim = ' ',
#     kin_diag = 0.5,
#     threads = 72,
#     num_chunks = 72,
#     outfile = "out.txt"
# )

# # opt = ConfOpt(
# #     pheno_add = "missing-simulationbyMengyu/phenotype_missing.txt",
# #     cov_add = "missing-simulationbyMengyu/cov_missing.txt",
# #     pheno_delim = "\t",
# #     cov_delim = "\t",
# #     geno_add = "missing-simulationbyMengyu/SA.bgen",
# #     sample_add = "missing-simulationbyMengyu/SA.sample",
# #     do_filters = False,
# #     use_sample_file = True,
# #     includeVariantFile = "",
# #     stream_snps = 10000,
# #     sampleid_header_name = "id",
# #     random_slope_header_name = "",
# #     covariates = ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"],
# #     exposures = [],
# #     interactions = [],
# #     missing_key = "NA",
# #     kin_add = "missing-simulationbyMengyu/kinfile_extended_fam.txt",
# #     kin_delim = '\t',
# #     kin_diag = 1,
# #     threads = 90,
# #     num_chunks = 90,
# #     outfile = "missing-x1-x10.txt"
# # )


# runner = GEMRunner(confopt.get()) #return opt obj from confopt obj

# # Run null model fitting
# runner.run_fit_nullmodel()


# print("Starting streaming dosage decode...")
# run_gwas(runner, runner.opt.outfile, snps_per_chunk=1000, device='cuda')
# end_time_gwas = time.time()
# print("End of running Torch GWAS file\n")
# print("Wall time in seconds :", end_time_gwas - start_time)
