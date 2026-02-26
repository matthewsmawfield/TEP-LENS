#!/usr/bin/env python3
"""
TEP-LENS: Step 01 - Fetch SN H0pe Lensing Time Delay Data

This script initializes the fundamental observables for the SN H0pe 
triply imaged supernova system (Pierel et al. 2024, Chen et al. 2024).
"""

import json
import os
import sys
from pathlib import Path

# Add project root to path for imports
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status

STEP_NUM = "01"

def safe_json_default(obj):
    """Handle numpy types in JSON serialization."""
    if hasattr(obj, 'item'):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

def main():
    print_status(f"STEP {STEP_NUM}: Initializing SN H0pe Data Structures", "TITLE")
    
    # Create required directories
    data_dir = PROJECT_ROOT / "data" / "raw" / "snh0pe"
    data_dir.mkdir(parents=True, exist_ok=True)
    out_dir = PROJECT_ROOT / "results" / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)

    print_status("Compiling literature data for SN H0pe (PLCK G165.7+67.0)")
    
    # -------------------------------------------------------------
    # LITERATURE DATA FOR SN H0PE
    # -------------------------------------------------------------
    snh0pe_data = {
        "metadata": {
            "name": "SN H0pe",
            "host_cluster": "PLCK G165.7+67.0",
            "z_lens": 0.351,
            "z_src": 1.783,
            "images": ["A", "B", "C"],
            "notes": "Data sourced from Pierel et al. 2024 (ApJ 967, 50) and Frye et al. 2024."
        },
        "time_delays": {
            "pierel_2024_photometric": {
                "dt_AB": {"value": -116.6, "err_plus": 10.8, "err_minus": 9.3},
                "dt_CB": {"value": -48.6, "err_plus": 3.6, "err_minus": 4.0},
                "dt_AC_derived": {"value": -68.0, "err": 11.5}
            }
        },
        "lens_model": {
            "mu_A": {"value": 5.4, "err": 0.5},
            "mu_B": {"value": 2.5, "err": 0.3},
            "mu_C": {"value": 2.0, "err": 0.2}
        }
    }

    # Save to raw
    raw_path = data_dir / "snh0pe_literature_data.json"
    with open(raw_path, 'w') as f:
        json.dump(snh0pe_data, f, indent=2, default=safe_json_default)
        
    print_status(f"Saved primary catalog to {raw_path}")
    
    # Save standard pipeline output
    output_path = out_dir / f"step_{STEP_NUM}_snh0pe_data.json"
    with open(output_path, 'w') as f:
        json.dump({"step": STEP_NUM, "status": "success", "data": snh0pe_data}, f, indent=2, default=safe_json_default)
        
    print_status(f"Results saved to {output_path}")
    print_status(f"Step {STEP_NUM} complete.")

if __name__ == "__main__":
    main()
