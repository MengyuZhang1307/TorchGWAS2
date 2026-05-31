#!/usr/bin/env python3
"""
Smoke test for PLINK BED input: Step 1 (null model) + Step 2 (dosage stream).
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "pymodules"))

from Mygen import GEMOptions, GEMRunner
from TorchGWAS import run_gwas

EXAMPLE = os.path.join(os.path.dirname(__file__), "example")
CORR    = "/tmp/test_plink_corr.txt"
OUT     = "/tmp/test_plink_out.parquet"

def make_opt(geno_file, corr_file=CORR):
    opt = GEMOptions()
    opt.pheno_add            = os.path.join(EXAMPLE, "example.pheno2-2id")
    opt.pheno_delim          = ','
    opt.cov_add              = os.path.join(EXAMPLE, "example.cov-2id")
    opt.cov_delim            = ','
    opt.geno_add             = geno_file
    opt.sampleid_header_name = "ID"
    opt.covariates           = ["cov1", "cov3"]
    opt.missing_key          = "NA"
    opt.threads              = 4
    opt.stream_snps          = 50
    opt.out_file             = "/tmp/test_plink_out.txt"
    opt.log_file             = "/tmp/test_plink.log"
    opt.null_log_file        = "/tmp/test_plink_null.log"
    opt.corr_file            = corr_file
    return opt

def test_bed_full():
    print("=" * 60)
    print("TEST: BED - Step 1 (null model)")
    print("=" * 60)
    bed = os.path.join(EXAMPLE, "example.bed")
    opt = make_opt(bed)

    # Step 1
    runner = GEMRunner(opt, False)
    runner.run_fit_nullmodel()
    print("  Step 1 OK")

    print()
    print("=" * 60)
    print("TEST: BED - Step 2 (dosage stream -> GWAS)")
    print("=" * 60)

    # Step 2: re-init with match_ids=True to align samples with correction file
    runner2 = GEMRunner(opt, True)
    run_gwas(runner2, CORR, OUT, snps_per_chunk=50, device="cpu")
    print(f"  Step 2 OK  ->  {OUT}")

if __name__ == "__main__":
    try:
        test_bed_full()
    except Exception as e:
        import traceback
        print(f"FAILED: {e}")
        traceback.print_exc()
