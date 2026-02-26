#!/usr/bin/env python3
"""
TEP-LENS: Step 02 - General Relativity Route-Closure Test

Calculates the expected relative time delays under standard GR and verifies 
the geometric closure identity: dt_AB + dt_BC + dt_CA = 0.
"""

import json
import sys
from pathlib import Path
import numpy as np

# Add project root to path for imports
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status

STEP_NUM = "02"

def safe_json_default(obj):
    if hasattr(obj, 'item'):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

def main():
    print_status(f"STEP {STEP_NUM}: GR Route-Closure Null Test", "TITLE")
    
    data_path = PROJECT_ROOT / "data" / "raw" / "snh0pe" / "snh0pe_literature_data.json"
    if not data_path.exists():
        print_status(f"Required data file missing: {data_path}", "ERROR")
        
    with open(data_path, 'r') as f:
        data = json.load(f)
        
    print_status("Extracting measured delays (Pierel et al. 2024)...")
    delays = data["time_delays"]["pierel_2024_photometric"]
    dt_AB = delays["dt_AB"]["value"]
    dt_CB = delays["dt_CB"]["value"]
    
    # By definition in GR:
    # t_A - t_B = dt_AB
    # t_C - t_B = dt_CB
    # Therefore, t_A - t_C = dt_AC = dt_AB - dt_CB
    # And dt_BC = -dt_CB
    # And dt_CA = -dt_AC
    
    dt_BC = -dt_CB
    dt_AC = dt_AB - dt_CB
    dt_CA = -dt_AC
    
    # GR Closure identity:
    closure_sum = dt_AB + dt_BC + dt_CA
    
    print_status(f"  Measured dt_AB: {dt_AB:.2f} days")
    print_status(f"  Measured dt_CB: {dt_CB:.2f} days")
    print_status(f"  Inferred dt_AC: {dt_AC:.2f} days")
    print_status(f"  GR Closure Sum (dt_AB + dt_BC + dt_CA): {closure_sum:.2f} days")
    
    results = {
        "step": STEP_NUM,
        "status": "success",
        "gr_predictions": {
            "dt_AB": dt_AB,
            "dt_BC": dt_BC,
            "dt_CA": dt_CA,
            "closure_sum": closure_sum,
            "closure_sum_description": "Under standard General Relativity, the route closure residual must be exactly 0."
        }
    }
    
    out_dir = PROJECT_ROOT / "results" / "outputs"
    output_path = out_dir / f"step_{STEP_NUM}_gr_closure.json"
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2, default=safe_json_default)
        
    print_status(f"Results saved to {output_path}")
    print_status(f"Step {STEP_NUM} complete.")
        
if __name__ == "__main__":
    main()
