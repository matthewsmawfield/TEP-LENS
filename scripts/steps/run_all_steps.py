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
    "step_01_fetch_snh0pe_data.py",
    "step_02_gr_closure.py",
    "step_03_tep_closure.py",
    "step_04_plot_closure.py"
]

def run_step(step_script):
    """Run a single python script step and stream its output."""
    script_path = STEPS_DIR / step_script
    if not script_path.exists():
        print_status(f"Script not found: {script_path}", "ERROR")
        return False
        
    print(f"\n======================================================================")
    print(f"RUNNING: {step_script}")
    print(f"======================================================================\n")
    
    try:
        result = subprocess.run(
            [sys.executable, str(script_path)],
            check=True,
            text=True
        )
        return True
    except subprocess.CalledProcessError as e:
        print_status(f"Pipeline stopped at {step_script} due to error.", "ERROR")
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
