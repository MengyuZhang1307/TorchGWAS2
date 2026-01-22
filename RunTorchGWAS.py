#to run on server $HOME/local/python3.11/bin/python3 Bgen.py or python3 Bgen.py 

import pandas as pd
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), "pymodules"))
from pymodules import ConfOpt
from pymodules import GEMRunner
from pymodules import run_gwas
from pymodules import parquet_to_text_duckdb
from pymodules import FDTee
import numpy as np
import time
import argparse
import logging
# from contextlib import redirect_stdout, redirect_stderr
import io
import tempfile
from pathlib import Path
import re
import traceback
# import threading
# import atexit
import faulthandler


def setup_step1_log(log_file, mode="a"):
    root = logging.getLogger()
    root.handlers.clear()
    root.setLevel(logging.INFO)
    root.propagate = False

    fh = logging.FileHandler(log_file, mode=mode, delay=False)
    fh.setFormatter(logging.Formatter("%(asctime)s [PY] %(levelname)s: %(message)s"))
    root.addHandler(fh)
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(logging.Formatter("[PY] %(levelname)s: %(message)s"))
    root.addHandler(sh)

_TEE = None

def setup_pipeline_log(log_path: str, mode: str = "a"):
    global _TEE

    if _TEE is None:
        # _TEE = FDTee(log_path, truncate=(mode == "w"))
        _TEE = FDTee(log_path, truncate=False, tee_to_terminal=True)

        root = logging.getLogger()
        root.setLevel(logging.INFO)
        root.handlers.clear()

        h = logging.StreamHandler(sys.stderr)  # goes through FDTee -> terminal + log
        fmt = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s",
                                "%Y-%m-%d %H:%M:%S")
        h.setFormatter(fmt)
        root.addHandler(h)

        def _excepthook(exc_type, exc, tb):
            logging.error("Uncaught Python exception:")
            logging.error("".join(traceback.format_exception(exc_type, exc, tb)).rstrip())
        sys.excepthook = _excepthook
    return log_path



def safe_stem(p: str) -> str:
    s = Path(p).stem
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s)

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
    parser.add_argument("--pheno-file", type=str, help="Phenotype file path (required for step1)")
    parser.add_argument("--cov-file", type=str, help="Covariate file path (required for step1 and step2)")
    parser.add_argument("--bgen", nargs="+", default=[], help="Step2: BGEN file(s).")
    parser.add_argument("--sample", nargs="+", default=[], help="Step2: SAMPLE file(s).")
    parser.add_argument("--kin-file", type=str, default="", help="Kinship file path (optional, required for step1 and step2 if using kinship)")
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
    parser.add_argument("--stream-snps", type=int, default=1000, help="Number of SNPs per chunk")
    parser.add_argument("--out", type=str, default="out.txt", help="Output file name")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cuda", help="Computation device (default: cuda)")
    parser.add_argument("--verbose", action="store_true", help="Print null model(default: False)")
    parser.add_argument("--convert", action="store_true", help="Convert binary to text file (default: True)")
    parser.add_argument("--step", choices=["all", "step1", "step2", "step3"], default="all",
                help=(
                    "Pipeline step to run:\n"
                    "all = perform all steps togethet\n"
                    "step1 = Fit the null model generate correction factors (write intermediate_*.txt)\n"
                    "step2 = TGWAS (write TGWAS_*.parquet)\n"
                    "step3 = Convert TGWAS_*.parquet to .txt\n"
                ))

    parser.add_argument("--parquet", nargs="+", default=[], help="Input TGWAS parquet file(s) for step3.")
    return parser.parse_args()

def build_logger_and_paths(args):
    """Only handles logger + dir/base names based on --out."""
    dir_name = os.path.dirname(args.out) or "."
    base_name = os.path.splitext(os.path.basename(args.out))[0]
    log_file = os.path.join(dir_name, base_name + ".log")
    return  dir_name, base_name

def build_conf_allsteps(args):
    """Config for all"""
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
        sampleid_header_name=args.sampleid_name,
        covariates=args.covar_names,
        random_slope_header_name=args.random_slope_name,
        missing_key=args.missing_value,
        kin_add=args.kin_file,
        kin_delim=normalize_delim(args.kin_delim),
        kin_diag=args.kin_diag,
        threads=args.threads,
        outfile=args.out,
        verbose=args.verbose,
    )
    return confopt

def build_conf_step1(args):
    """Config for step1: needs phenotype + covariates, genotype files, kin."""
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
        sampleid_header_name=args.sampleid_name,
        covariates=args.covar_names,
        random_slope_header_name=args.random_slope_name,
        missing_key=args.missing_value,
        kin_add=args.kin_file,
        kin_delim=normalize_delim(args.kin_delim),
        kin_diag=args.kin_diag,
        threads=args.threads,
        outfile=args.out,
        verbose=args.verbose,
    )
    return confopt


def build_conf_step2(args):
    """Config for step2: NO phenotype, but covariates, genotype files, kin."""
    confopt = ConfOpt(
        cov_add=args.cov_file,
        cov_delim=normalize_delim(args.cov_delim),
        geno_add=args.bgen,
        sample_add=args.sample,
        use_sample_file=bool(args.sample),
        do_filters=bool(args.include_snp_file),
        includeVariantFile=args.include_snp_file,
        stream_snps=args.stream_snps,
        sampleid_header_name=args.sampleid_name,
        covariates=args.covar_names,
        random_slope_header_name=args.random_slope_name,
        missing_key=args.missing_value,
        kin_add=args.kin_file,
        kin_delim=normalize_delim(args.kin_delim),
        kin_diag=args.kin_diag,
        threads=args.threads,
        outfile=args.out,
        verbose=args.verbose,
    )
    return confopt

def run_all(dir_name, base_name, args, log_file):
    # ------------------
    # STEP 1 (once)
    # ------------------
    logging.info("%s", "*" * 80)
    logging.info("STEP 1: Fitting null model")
    logging.info("%s", "*" * 80)
    sub_step1 = argparse.Namespace(**vars(args))
    sub_step1.bgen = args.bgen[0]
    sub_step1.sample = args.sample[0]
    conf_step1 = build_conf_step1(sub_step1)

    intermediate_file = os.path.join(dir_name, "intermediate_" + base_name + ".txt")
    runner = GEMRunner(conf_step1.get())        
    #Null model
    runner.run_fit_nullmodel()
    logging.info("Intermediate file (correction) path: %s", intermediate_file)
    setup_pipeline_log(log_file, mode="a")
    # ------------------
    # STEP 2 (loop)
    # ------------------
    print("*" * 80)
    print("STEP 2: Running Torch GWAS")
    print("*" * 80)
    for bgen_i, sample_i in zip(args.bgen, args.sample):
        sub_step2 = argparse.Namespace(**vars(args))
        sub_step2.bgen = bgen_i
        sub_step2.sample = sample_i

        base_i = safe_stem(bgen_i) + "_" + base_name

        conf_step2 = build_conf_step2(sub_step2)

        TGWAS_file = os.path.join(dir_name, base_i + ".parquet")
        runner = GEMRunner(conf_step2.get(), True) 
        run_gwas(
            runner,
            intermediate_file,               # correction file
            TGWAS_file,
            snps_per_chunk=args.stream_snps,
            device=args.device,
        )
        print(f"TGWAS parquet output: {TGWAS_file}")


    # ------------------
    # STEP 3
    # ------------------
    if args.convert:
        print("*" * 80)
        print("STEP 3: Converting binary to text")
        print("*" * 80)
        for bgen_i in args.bgen:
            base_i = safe_stem(bgen_i) + "_" + base_name
            parquet_file = os.path.join(dir_name, base_i + ".parquet")
            output_file = os.path.join(dir_name, base_i + ".txt")

            if not os.path.exists(parquet_file):
                logging.error(f"Parquet file not found, skipping: {parquet_file}")
                raise SystemExit(2)

            # print(f"Converting: {parquet_file} -> {output_file}")
            print(f"STEP 3: Converting {parquet_file} -> {output_file}")
            parquet_to_text_duckdb(parquet_file, output_file)


# Broken steps:
def run_step1(confopt, dir_name, base_name, args):
    """
    STEP 1:
      - GEMRunner init
      - run_fit_nullmodel
      - log intermediate + TGWAS filenames
      -run_gwas(runner, correction, TGWAS_file, ...)
    """
    intermediate_file = os.path.join(dir_name, "intermediate_" + base_name + ".txt")
    # 1) C++ init
    runner = GEMRunner(confopt.get())

    # 2) Null model
    runner.run_fit_nullmodel()

    # print(f"Intermediate file (correction) path: {intermediate_file}")    
    logging.info("Intermediate file (correction) path: %s", intermediate_file)

def run_step2(confopt, dir_name, base_i, base_name, args):
    """
    STEP 2:
      - GEMRunner init
      - run_gwas(runner, correction, parquet file, ...)
    """
    intermediate_file = os.path.join(dir_name, "intermediate_" + base_name + ".txt")
    TGWAS_file = os.path.join(dir_name, base_i + ".parquet")

    # print("STEP 2: Re-initializing GEMRunner and running GWAS/TGWAS...")
    print(f"TGWAS parquet output: {TGWAS_file}")

    runner = GEMRunner(confopt.get(), True) # True to match IDs for each genotype with intermediate file

    print("Starting GWAS/TGWAS with run_gwas...")

    run_gwas(
        runner,
        intermediate_file,               # correction file
        TGWAS_file,
        snps_per_chunk=args.stream_snps,
        device=args.device,
    )

def run_step3(args):
    """
    STEP 3:
        - Convert TGWAS_<base_name>.parquet -> <base_name>.txt
        - TGWAS_file (input parquet)
        - output_file (text)
        - logger
    """
    print("STEP 3: Converting binary to text")
    if not args.parquet:
        print("STEP 3 requires --parquet  (output of step2).")
        raise SystemExit(2)

    parquet_file = args.parquet
    # output_file = os.path.join(dir_name, base_name + ".txt")
    output_file = args.out 

    print(f"STEP 3: Converting {parquet_file} -> {output_file}")
    parquet_to_text_duckdb(parquet_file, output_file)


def main():
    global _TEE
    start_time = time.time()
    args = parse_args()
    crash_fp = open(args.out + ".crash.log", "w", buffering=1)
    faulthandler.enable(file=crash_fp, all_threads=True)
    crash_fp.write("\n==== crash log start ====\n")
    crash_fp.flush()
    try:
        dir_name, base_name = build_logger_and_paths(args)
        log_file = os.path.join(dir_name, base_name + ".log")
        if args.step == "all":
            setup_step1_log(log_file, mode="w")
            run_all(dir_name, base_name, args, log_file)

        # Step-specific requirements
        if args.step == "step1":
            setup_step1_log(log_file, mode="w")
            logging.info("%s", "*" * 80)
            logging.info("STEP 1: fitting null model...")
            if getattr(args, "convert", False):
                logging.warning("--convert is only used with --step all. Ignoring it for --step step1.")
            if not args.pheno_file:
                logging.error("STEP 1 requires --pheno-file.")
                raise SystemExit(2)
            if len(args.bgen) != 1:
                logging.error("STEP 1 requires exactly ONE --bgen file.")
                raise SystemExit(2)
            if len(args.sample) != 1:
                logging.error("STEP 1 requires exactly ONE --sample file.")
                raise SystemExit(2)

            args.bgen = args.bgen[0]
            args.sample = args.sample[0]
            confopt = build_conf_step1(args)
            run_step1(confopt, dir_name, base_name, args)

        elif args.step == "step2":
            setup_pipeline_log(log_file, mode="a")
            if getattr(args, "pheno_file", None):
                print("WARNING: --pheno-file is not used in step2; ignoring it for --step step2.", file=sys.stderr)
            if getattr(args, "convert", False):
                print("WARNING: --convert is only used with --step all. Ignoring it for --step step2.", file=sys.stderr)
            if len(args.bgen) != len(args.sample):
                print(f"--bgen count ({len(args.bgen)}) must match --sample count ({len(args.sample)}).")
                raise SystemExit(2)
            for bgen_i, sample_i in zip(args.bgen, args.sample):
                sub = argparse.Namespace(**vars(args))
                sub.bgen = bgen_i   
                sub.sample = sample_i  

                base_i = safe_stem(bgen_i) + "_" + base_name   # output: TGWAS_<base_i>.parquet
                print("*" * 80)
                print(f"STEP 2 batch item: bgen={bgen_i} -> sample={sample_i}")
                confopt = build_conf_step2(sub)
                run_step2(confopt, dir_name, base_i, base_name, sub)

        elif args.step == "step3":
            setup_pipeline_log(log_file, mode="a")
            for pq in args.parquet:
                sub = argparse.Namespace(**vars(args))
                sub.parquet = pq   # use string per run 
                sub.out = str(Path(pq).with_suffix(".txt"))
                print("*" * 80)
                run_step3(sub)

        end_time = time.time()
        print("\nTorchGWAS pipeline step completed successfully.")
        print(f"Wall time: {(end_time - start_time):.2f} seconds")
    except SystemExit:
        # keep SystemExit behavior (argparse / your raise SystemExit)
        raise
    except Exception:
        # THIS will now go into the log because tee is still active
        logging.exception("Uncaught Python exception:")
        raise SystemExit(1)
    
    finally:
        try:
            logging.shutdown() 
        finally:
            if _TEE is not None:
                _TEE.close()
                _TEE = None

if __name__ == "__main__":
    main()
