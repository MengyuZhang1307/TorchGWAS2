#to run on server $HOME/local/python3.11/bin/python3 Bgen.py or python3 Bgen.py 

import pandas as pd
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), "pymodules"))
from pymodules import ConfOpt
from pymodules import GEMRunner
from pymodules import run_gwas
from pymodules import parquet_to_text_duckdb
import numpy as np
import time
import argparse
import logging
from contextlib import redirect_stdout, redirect_stderr
import io
import tempfile
from pathlib import Path
import re
import traceback

def _ensure_parent_dir(path: str):
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)


def setup_pipeline_log(out_prefix: str, step_name: str, mode: str = "a") -> str:
    """
    mode:
      - "w" => create/truncate (Step 1)
      - "a" => append (Step 2/3)
    """
    if mode not in ("w", "a"):
        raise ValueError("mode must be 'w' or 'a'")

    log_path = out_prefix + ".log"
    log_dir = os.path.dirname(log_path)
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)

    # If Step 1: truncate first
    if mode == "w":
        with open(log_path, "w"):
            pass

    # OS-level redirect (append from now on)
    fd = os.open(log_path, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o644)
    os.dup2(fd, 2)  # stderr
    os.dup2(fd, 1)  # stdout
    os.close(fd)

    # Re-wrap Python streams (line-buffered)
    sys.stdout = os.fdopen(1, "w", buffering=1)
    sys.stderr = os.fdopen(2, "w", buffering=1)

    # Python logging -> stderr (now the log file)
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    root.handlers.clear()

    h = logging.StreamHandler(sys.stderr)
    fmt = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s", "%Y-%m-%d %H:%M:%S")
    h.setFormatter(fmt)
    root.addHandler(h)

    logging.info("=" * 72)
    logging.info(step_name)
    logging.info(f"LOG FILE: {log_path}  (mode={mode})")
    logging.info("=" * 72)

    # Uncaught exceptions -> log
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

def run_all(dir_name, base_name, args):
    if args.step == "all":
        # ------------------
        # STEP 1 (once)
        # ------------------
        
        print("STEP 1: Fitting null model")

        sub_step1 = argparse.Namespace(**vars(args))
        sub_step1.bgen = args.bgen[0]
        sub_step1.sample = args.sample[0]

        conf_step1 = build_conf_step1(sub_step1)
        # run_step1(conf_step1, logger, dir_name, base_name, sub_step1)
        intermediate_file = os.path.join(dir_name, "intermediate_" + base_name + ".txt")
        runner = GEMRunner(conf_step1.get())        
        #Null model
        runner.run_fit_nullmodel()
        # ------------------
        # STEP 2 (loop)
        # ------------------
        print("STEP 2: Running Torch GWAS")
        for bgen_i, sample_i in zip(args.bgen, args.sample):
            sub_step2 = argparse.Namespace(**vars(args))
            sub_step2.bgen = bgen_i
            sub_step2.sample = sample_i

            base_i = safe_stem(bgen_i) + "_" + base_name

            conf_step2 = build_conf_step2(sub_step2)
            # run_step2(confopt, logger, dir_name, base_i, sub)
            TGWAS_file = os.path.join(dir_name, base_i + ".parquet")
            runner = GEMRunner(conf_step2.get(), True) 
            run_gwas(
                runner,
                intermediate_file,               # correction file
                TGWAS_file,
                snps_per_chunk=args.stream_snps,
                device=args.device,
            )


        # ------------------
        # STEP 3
        # ------------------
        if args.convert:
            print("STEP 3: Converting binary to text")
            for bgen_i in args.bgen:
                base_i = safe_stem(bgen_i) + "_" + base_name
                parquet_file = os.path.join(dir_name, base_i + ".parquet")
                txt_file = os.path.join(dir_name, base_i + ".txt")

                if not os.path.exists(parquet_file):
                    raise SystemExit(f"Parquet file not found, skipping: {parquet_file}")

                print(f"Converting: {parquet_file} -> {txt_file}")
                parquet_to_text_duckdb(parquet_file, txt_file)


# Broken steps:
def run_step1(confopt, logger, dir_name, base_name, args):
    """
    STEP 1:
      - GEMRunner init
      - run_fit_nullmodel
      - log intermediate + TGWAS filenames
      -run_gwas(runner, correction, TGWAS_file, ...)
    """
    intermediate_file = os.path.join(dir_name, "intermediate_" + base_name + ".txt")
    # TGWAS_file = os.path.join(dir_name, "TGWAS_" + base_name + ".parquet")

    logger.info("STEP 1: Initializing GEMRunner and fitting null model...")

    # 1) C++ init
    runner = GEMRunner(confopt.get())

    # 2) Null model
    logger.info("Running null model fitting ...")
    runner.run_fit_nullmodel()

    logger.info(f"Intermediate file (correction) path: {intermediate_file}")
    # logger.info(f"TGWAS parquet file (for step2/step3): {TGWAS_file}")
    logger.info(f"Run step2 with: --correction-add {intermediate_file}")

 
def run_step2(confopt, logger, dir_name, base_name, args):
    """
    STEP 2:
      - GEMRunner init
      - run_gwas(runner, correction, parquet file, ...)
    """
    intermediate_file = os.path.join(dir_name, "intermediate_" + base_name + ".txt")
    TGWAS_file = os.path.join(dir_name, base_name + ".parquet")

    logger.info("STEP 2: Re-initializing GEMRunner and running GWAS/TGWAS...")
    logger.info(f"Using correction (intermediate) file: {intermediate_file}")
    logger.info(f"TGWAS parquet output: {TGWAS_file}")

    runner = GEMRunner(confopt.get(), True) # True to match IDs for each genotype with intermediate file

    logger.info("Starting GWAS/TGWAS with run_gwas...")

    run_gwas(
        runner,
        intermediate_file,               # correction file
        TGWAS_file,
        snps_per_chunk=args.stream_snps,
        device=args.device,
    )

def run_step3(logger, args):
    """
    STEP 3:
        - Convert TGWAS_<base_name>.parquet -> <base_name>.txt
        - TGWAS_file (input parquet)
        - output_file (text)
        - logger
    """
    if not args.parquet:
        raise SystemExit("STEP 2 requires --parquet  (output of step2).")

    parquet_file = args.parquet
    # output_file = os.path.join(dir_name, base_name + ".txt")
    output_file = args.out 

    logger.info(f"STEP 3: Converting {parquet_file} -> {output_file} ...")

    parquet_to_text_duckdb(parquet_file, output_file)


def main():
    start_time = time.time()
    args = parse_args()
    dir_name, base_name = build_logger_and_paths(args)
    if args.step == "all":
        setup_pipeline_log(args.out, "STEP 1: FIT NULL MODEL", mode="w")
        run_all(dir_name, base_name, args)

    # Step-specific requirements


    if args.step == "step1":
        setup_pipeline_log(args.out, "STEP 1: FIT NULL MODEL", mode="w")
        if not args.pheno_file:
            raise SystemExit("STEP 1 requires --pheno-file.")
        if len(args.bgen) != 1:
            raise SystemExit("STEP 1 requires exactly ONE --bgen file.")
        if len(args.sample) != 1:
            raise SystemExit("STEP 1 requires exactly ONE --sample file.")

    # convert list -> string for pybind GEMOptions
        args.bgen = args.bgen[0]
        args.sample = args.sample[0]
        confopt = build_conf_step1(args)
        run_step1(confopt, logger, dir_name, base_name, args)

    elif args.step == "step2":
        setup_pipeline_log(args.out, "STEP 1: FIT NULL MODEL", mode="a")
        if len(args.bgen) != len(args.sample):
            raise SystemExit(f"--bgen count ({len(args.bgen)}) must match --sample count ({len(args.sample)}).")
        for bgen_i, sample_i in zip(args.bgen, args.sample):
            sub = argparse.Namespace(**vars(args))
            sub.bgen = bgen_i   
            sub.sample = sample_i  

            base_i = safe_stem(bgen_i) + "_" + base_name   # output: TGWAS_<base_i>.parquet
            logger.info("*" * 80)
            logger.info(f"STEP 2 batch item: bgen={bgen_i} -> sample={sample_i}")
            confopt = build_conf_step2(sub)
            run_step2(confopt, logger, dir_name, base_i, sub)

    elif args.step == "step3":
        setup_pipeline_log(args.out, "STEP 1: FIT NULL MODEL", mode="a")
        for pq in args.parquet:
            sub = argparse.Namespace(**vars(args))

            # keep parquet as a single file for this run
            sub.parquet = pq   # use string per run 

            # output name: same stem, .txt
            sub.out = str(Path(pq).with_suffix(".txt"))

            logger.info("*" * 80)
            run_step3(logger, sub)

    end_time = time.time()
    print("\nTorchGWAS pipeline step completed successfully.")
    print(f"Wall time: {(end_time - start_time):.2f} seconds")

if __name__ == "__main__":
    main()
