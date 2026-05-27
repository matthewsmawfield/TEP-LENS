#!/usr/bin/env python3
"""
TEP-LENS: Run All Analysis Steps

This script runs the complete reproducible analysis pipeline for the 
TEP Lensing Time Delay Closure Test.

Usage:
    python scripts/steps/run_all_steps.py
"""

import json
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
STEPS_DIR = PROJECT_ROOT / "scripts" / "steps"
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status

STEPS = [
    "step_00_fetch_literature_and_cross_paper_data.py",  # arXiv + TEP cross-paper coefficients
    "step_01_fetch_snh0pe_data.py",   # Compile SN Refsdal + H0pe + 2025wny catalog
    "step_02_gr_closure.py",            # GR null test: all 5 loops = 0
    "step_03_tep_closure.py",           # TEP closure residuals: best loop SNR=78
    "step_04_plot_closure.py",          # Figures: loop comparison + baseline scaling
    "step_05_tdcosmo_shear.py",         # TDCOSMO quad-lens temporal shear test
    "step_06_alpha_sensitivity.py",     # Alpha scan: SNR is alpha-independent (geometric)
    "step_07_observed_vs_predicted.py", # Evidence: observed vs 7 blind model predictions
    "step_08_new_evidence.py",          # Extended evidence: delay-mu corr, alpha inference, H0pe sensitivity
    "step_09_precision_roadmap.py",     # TEP precision roadmap simulation
    "step_10_h0_tension.py",            # H0 tension resolution
    "step_11_model_dependence.py",      # Model dependence robustness (N_eff, LOO, exact sign-flip)
    "step_12_microlensing_robustness.py", # Microlensing nuisance Monte Carlo robustness
    "step_13_bayes_model_comparison.py", # Bayesian GR vs TEP comparison with hierarchical nuisance
    "step_14_external_chain_ingestion.py", # Public H0LiCOW/TDCOSMO chain ingestion and coverage summary
    "step_15_external_informed_inflation.py", # External-chain-informed uncertainty inflation stress test
    "step_16_independence_tier_significance.py", # Tiered significance synthesis (independence-aware)
    "step_17_directional_odds.py", # Directional Bayes-factor expansion for sign-based evidence
    "step_18_external_dataset_registry.py", # External-dataset shortlist + local ingestion-readiness for TEP-LENS
    "step_19_tdcosmo2025_ingestion.py", # Standalone TDCOSMO2025 public ingestion (local-first, API fallback)
    "step_30_cosmograil_temporal_shear.py",       # COSMOGRAIL temporal shear analysis
    "step_31_cosmograil_validation.py",           # COSMOGRAIL validation
    "step_32_kappa_proxy_validation.py",          # Kappa proxy validation with elliptical+shear comparison
    "step_33_einstein_cross_null.py",             # Einstein-cross null test (single-lever fragility)
    "step_34_agnostic_rank_test.py",              # Ansatz-free rank test
    "step_35_single_contrast_dominance.py",       # Signal-energy concentration and effective degrees of freedom
    "step_20_external_completeness_synthesis.py", # Completeness-aware Tier-A significance synthesis from steps 16+19
]

def run_step(step_script):
    """Run a single python script step and stream its output to console and log file."""
    script_path = STEPS_DIR / step_script
    log_path = PROJECT_ROOT / "logs" / f"{step_script.replace('.py', '')}.log"
    
    if not script_path.exists():
        print_status(f"Script not found: {script_path}", "ERROR")
        return False
        
    print(f"\n{'='*70}")
    print(f"RUNNING: {step_script}")
    print(f"LOG: {log_path}")
    print(f"{'='*70}\n")
    
    try:
        with open(log_path, "w") as log_file:
            # Use Popen to stream output to both console and file
            process = subprocess.Popen(
                [sys.executable, str(script_path)],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1
            )
            
            # Stream output
            for line in process.stdout:
                sys.stdout.write(line)
                log_file.write(line)
                
            process.wait()
            
        if process.returncode != 0:
            print_status(f"Pipeline stopped at {step_script} due to error. See {log_path}", "ERROR")
            return False
            
        return True
    except Exception as e:
        print_status(f"Failed to execute {step_script}: {e}", "ERROR")
        return False

def main():
    print_status("TEP-LENS PIPELINE INITIATED", "TITLE")
    
    successful_steps = 0
    for step in STEPS:
        success = run_step(step)
        if success:
            successful_steps += 1
        else:
            break
            
    print(f"\n======================================================================")
    print(f"PIPELINE COMPLETE")
    print(f"======================================================================")
    print(f"Steps completed: {successful_steps}/{len(STEPS)}")
    print(f"Results saved to: results/outputs/")
    print(f"Figures saved to: results/figures/")

if __name__ == "__main__":
    main()
