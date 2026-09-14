#!/usr/bin/env python3
"""
scripts/run_pipeline.py
Master End-to-End Execution Script for Targeted Nanopore Analysis (nCATS).
Executes QC, Alignment Diagnostics, Variant Calling, Phasing, Methylation, STR/SV, and Dashboard.
"""

import os
import sys
import subprocess
import time

def run_step(step_name, command):
    print("\n" + "="*70)
    print(f"▶ EXECUTING: {step_name}")
    print(f"  Command: {command}")
    print("="*70)
    start_time = time.time()
    
    res = subprocess.run(command, shell=True)
    elapsed = time.time() - start_time
    
    if res.returncode != 0:
        print(f"❌ ERROR in step: {step_name} (Exit code: {res.returncode})")
        sys.exit(res.returncode)
    else:
        print(f"✔ COMPLETED: {step_name} in {elapsed:.2f} seconds.")

def main():
    print("**********************************************************************")
    print("  TARGETED NANOPORE GENOMIC & EPIGENETIC ANALYSIS (nCATS) PIPELINE   ")
    print("**********************************************************************")
    total_start = time.time()

    # Step 1: Quality Control
    run_step("1. Raw Read Quality Control", "python3 qc.py")

    # Step 2: Target Coverage
    run_step("2. Target Coverage and Depth Evaluation", "python3 plot_coverage.py")

    # Step 3: Variant Calling with Dual-Strand Filtering
    run_step("3. Long-Read Variant Calling & Dual-Strand Filter", "python3 scripts/variant_caller.py")

    # Step 4: Long-Read Haplotype Phasing & Haplotagging
    run_step("4. Read-Backed Haplotype Phasing & BAM Haplotagging", "python3 scripts/phasing.py")

    # Step 5: Native CpG Methylation & WGBS Validation
    run_step("5. Native CpG 5mC Epigenetic Profiling & WGBS Benchmark", "python3 scripts/methylation.py")

    # Step 6: Forensic STR Profiling & Structural Variant Calling
    run_step("6. Forensic STR & Structural Variant Analysis", "python3 scripts/str_sv_analyzer.py")

    # Step 7: Publication Visualizations & Multi-Omic Dashboard
    run_step("7. Multi-Omic Publication Visualizations & Dashboard", "python3 scripts/multiomic_dashboard.py")

    total_elapsed = time.time() - total_start
    print("\n" + "*"*70)
    print(f"🎉 ENTIRE PIPELINE COMPLETED SUCCESSFULLY IN {total_elapsed:.2f} SECONDS!")
    print("   All results, VCFs, TSV tables, and 300 DPI visualizations are saved in results/")
    print("*"*70)

if __name__ == "__main__":
    main()
