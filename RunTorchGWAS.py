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

def safe_stem(p: str) -> str:
    s = Path(p).stem
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s)
class CaptureCStdout:
    def __enter__(self):
        self._orig_stdout_fd = sys.stdout.fileno()

        # Save a duplicate of the original file descriptor
        self._saved_stdout_fd = os.dup(self._orig_stdout_fd)

        # Create a temporary file to capture output
        self._tmpfile = tempfile.TemporaryFile(mode="w+b")

        # Redirect stdout to the temporary file
        os.dup2(self._tmpfile.fileno(), self._orig_stdout_fd)

        return self

    def __exit__(self, exc_type, exc_value, traceback):
        # Restore the original stdout
        os.dup2(self._saved_stdout_fd, self._orig_stdout_fd)

        # Read captured output
        self._tmpfile.seek(0)
        self.output = self._tmpfile.read().decode()

        # Cleanup
        self._tmpfile.close()
        os.close(self._saved_stdout_fd)
class CaptureCStderr:
    def __enter__(self):
        self.fd = sys.stderr.fileno()
        self.saved_fd = os.dup(self.fd)
        self.tmp = tempfile.TemporaryFile(mode="w+b")
        os.dup2(self.tmp.fileno(), self.fd)
        return self

    def __exit__(self, *args):
        os.dup2(self.saved_fd, self.fd)
        os.close(self.saved_fd)
        self.tmp.seek(0)
        self.output = self.tmp.read().decode()
        self.tmp.close()

def setup_logger(out_path, step, truncate=False):
    """
    Create a logger that prints to both file and console.
    If truncate=True, overwrite the log file.
    If truncate=False, append to the existing file.
    """
    logger = logging.getLogger("TGWAS")
    logger.setLevel(logging.INFO)

    # IMPORTANT: avoid duplicate handlers if called multiple times
    if logger.handlers:
        logger.handlers.clear()

    fmt = logging.Formatter(
        "[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    mode = "w" if truncate or step == "all" else "a"
    fh = logging.FileHandler(out_path, mode=mode)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    ch = logging.StreamHandler()
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    return logger

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
    parser.add_argument("--correction-add", dest="correction_add", type=str, default="", help="Path to intermediate file produced by step1 (used in step2).")
    parser.add_argument("--parquet", nargs="+", default=[], help="Input TGWAS parquet file(s) for step3.")
    return parser.parse_args()

def build_logger_and_paths(args):
    """Only handles logger + dir/base names based on --out."""
    dir_name = os.path.dirname(args.out) or "."
    base_name = os.path.splitext(os.path.basename(args.out))[0]
    log_file = os.path.join(dir_name, base_name + ".log")

    # For step1 we truncate; for step2/3 we append to the same log
    truncate = (args.step == "step1")

    logger = setup_logger(log_file, args.step, truncate=truncate)
    return logger, dir_name, base_name

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
def run_all(confopt, logger, dir_name, base_name, args):
    """
    All STEP:
      - GEMRunner init
      - run_fit_nullmodel
      - log intermediate + TGWAS filenames
      - Convert TGWAS_<base_name>.parquet -> <base_name>.txt
    """
    intermediate_file = os.path.join(dir_name, "intermediate_" + base_name + ".txt")
    TGWAS_file = os.path.join(dir_name, "TGWAS_" + base_name + ".parquet")

    logger.info("STEP 1: Initializing GEMRunner and fitting null model...")

    # 1) C++ init
    # with CaptureCStdout() as cap_init, CaptureCStderr() as cap_init_err:
    #     runner = GEMRunner(confopt.get())
    runner = GEMRunner(confopt.get(), False)
    # init_output = (cap_init.output + "\n" + cap_init_err.output).strip()
    # if init_output:
    #     logger.info("\n********** C++ Initialization Output **********\n" + init_output)

    # 2) Null model
    logger.info("Running null model fitting ...")
    # with CaptureCStdout() as cap_out, CaptureCStderr() as cap_err:
    #     runner.run_fit_nullmodel()
    runner.run_fit_nullmodel()
    # merged = (cap_out.output + "\n" + cap_err.output).strip()
    # if merged:
    #     logger.info(
    #         "\n****************************** C++ Null Model Output ******************************\n"
    #         + merged
    #     )
    # else:
    #     logger.info("No C++ output captured from null model.")

    # logger.info(f"Intermediate file (correction) path: {intermediate_file}")
    # # logger.info(f"TGWAS parquet file (for step2/step3): {TGWAS_file}")
    # logger.info(f"Run step2 with: --correction-add {intermediate_file}")
    # logger.info("STEP 2: Re-initializing GEMRunner and running GWAS/TGWAS...")
    # logger.info(f"Using correction (intermediate) file: {intermediate_file}")
    # logger.info(f"TGWAS parquet output: {TGWAS_file}")
    # logger.info("Starting GWAS/TGWAS with run_gwas...")
    # cxx_buffer = io.StringIO()
    # with redirect_stdout(cxx_buffer), redirect_stderr(cxx_buffer):
    #     run_gwas(
    #         runner,
    #         intermediate_file,               # correction file
    #         TGWAS_file,
    #         snps_per_chunk=args.stream_snps,
    #         device=args.device,
    #     )
    run_gwas(
            runner,
            intermediate_file,               # correction file
            TGWAS_file,
            snps_per_chunk=args.stream_snps,
            device=args.device,
        )
    # captured_output = cxx_buffer.getvalue().strip()
    # if captured_output:
    #     logger.info(
    #         "\n****************************** TGWAS Output ******************************\n"
    #         + captured_output)
    output_file = os.path.join(dir_name, base_name + ".txt")
    # logger.info(f"STEP 3: Converting {TGWAS_file} -> {output_file} ...")
    # cxx_buffer.seek(0)
    # cxx_buffer.truncate(0)
    # if args.convert:
    #     with redirect_stdout(cxx_buffer), redirect_stderr(cxx_buffer):
    #         parquet_to_text_duckdb(TGWAS_file, output_file)
    if args.convert:
        parquet_to_text_duckdb(TGWAS_file, output_file)
    # captured_output_conversion = cxx_buffer.getvalue().strip()
    # if captured_output_conversion:
    #     logger.info(
    #         "\n****************************** Conversion of Output Binary to Text ******************************\n"
    #         + captured_output_conversion)


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
    with CaptureCStdout() as cap_init, CaptureCStderr() as cap_init_err:
        runner = GEMRunner(confopt.get())

    init_output = (cap_init.output + "\n" + cap_init_err.output).strip()
    if init_output:
        logger.info("\n********** C++ Initialization Output **********\n" + init_output)

    # 2) Null model
    logger.info("Running null model fitting ...")
    with CaptureCStdout() as cap_out, CaptureCStderr() as cap_err:
        runner.run_fit_nullmodel()

    merged = (cap_out.output + "\n" + cap_err.output).strip()
    if merged:
        logger.info(
            "\n****************************** C++ Null Model Output ******************************\n"
            + merged
        )
    else:
        logger.info("No C++ output captured from null model.")

    logger.info(f"Intermediate file (correction) path: {intermediate_file}")
    # logger.info(f"TGWAS parquet file (for step2/step3): {TGWAS_file}")
    logger.info(f"Run step2 with: --correction-add {intermediate_file}")

 
def run_step2(confopt, logger, dir_name, base_name, args):
    """
    STEP 2:
      - GEMRunner init
      - run_gwas(runner, correction, TGWAS_file, ...)
    """
    if not args.correction_add:
        raise SystemExit("STEP 2 requires --correction <intermediate_file> (output of step1).")

    intermediate_file = args.correction_add                     # this *is* the correction file
    TGWAS_file = os.path.join(dir_name, base_name + ".parquet")

    logger.info("STEP 2: Re-initializing GEMRunner and running GWAS/TGWAS...")
    logger.info(f"Using correction (intermediate) file: {intermediate_file}")
    logger.info(f"TGWAS parquet output: {TGWAS_file}")

    with CaptureCStdout() as cap_init, CaptureCStderr() as cap_init_err:
        runner = GEMRunner(confopt.get(), True) # True to match IDs for eacg genotype with intermediate file

    init_output = (cap_init.output + "\n" + cap_init_err.output).strip()
    if init_output:
        logger.info(
            "\n********** C++ Initialization Output (Step 2) **********\n"
            + init_output
        )

    cxx_buffer = io.StringIO()
    logger.info("Starting GWAS/TGWAS with run_gwas...")

    with redirect_stdout(cxx_buffer), redirect_stderr(cxx_buffer):
        run_gwas(
            runner,
            intermediate_file,               # correction file
            TGWAS_file,
            snps_per_chunk=args.stream_snps,
            device=args.device,
        )

    captured_output = cxx_buffer.getvalue().strip()
    if captured_output:
        logger.info(
            "\n****************************** TGWAS Output ******************************\n"
            + captured_output)

def run_step3(logger, args):
    """
    STEP 3:
        - Convert TGWAS_<base_name>.parquet -> <base_name>.txt
        - TGWAS_file (input parquet)
        - output_file (text)
        - logger
    """
    if not args.parquet:
        raise SystemExit("STEP 2 requires --parquet <TGWAS_file> (output of step2).")

    TGWAS_file = args.parquet
    # output_file = os.path.join(dir_name, base_name + ".txt")
    output_file = args.out 
    cxx_buffer = io.StringIO()
    logger.info(f"STEP 3: Converting {TGWAS_file} -> {output_file} ...")

    with redirect_stdout(cxx_buffer), redirect_stderr(cxx_buffer):
        parquet_to_text_duckdb(TGWAS_file, output_file)

    captured_output_conversion = cxx_buffer.getvalue().strip()
    if captured_output_conversion:
        logger.info(
            "\n****************************** Conversion of Output Binary to Text ******************************\n"
            + captured_output_conversion)

def main():
    start_time = time.time()
    args = parse_args()
    logger, dir_name, base_name = build_logger_and_paths(args)
    if args.step == "all":
        confopt = build_conf_allsteps(args)
        run_all(confopt, logger, dir_name, base_name, args)

    # Step-specific requirements
    if args.step == "step1" and not args.pheno_file:
        raise SystemExit("STEP 1 requires --pheno-file.")


    if args.step == "step1":
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
        for pq in args.parquet:
            sub = argparse.Namespace(**vars(args))

            # keep parquet as a single file for this run
            sub.parquet = pq   # use string per run 

            # output name: same stem, .txt
            sub.out = str(Path(pq).with_suffix(".txt"))

            logger.info("*" * 80)
            run_step3(logger, sub)

    end_time = time.time()
    logger.info("\nTorchGWAS pipeline step completed successfully.")
    logger.info(f"Wall time: {(end_time - start_time):.2f} seconds")

if __name__ == "__main__":
    main()
