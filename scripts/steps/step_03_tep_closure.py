#!/usr/bin/env python3
"""
TEP-LENS: Step 03 - TEP Route-Closure Test

Calculates the predicted TEP closure residual for SN H0pe.

Under TEP, the effective speed of light scales as c' = c / Gamma_t.
This means the geometric travel time dt_geom is stretched by Gamma_t.

For a closure test, we use the fact that the light paths traverse different 
gravitational potential depths through the lens cluster.

Assuming we have macro-model magnifications mu_A, mu_B, mu_C, these 
correlate with the local potential depth at the image positions.
For a first-order prediction, Gamma_t scales with the projected mass density, 
which correlates with magnification.
Gamma_i = 1 + alpha * log10(mu_i)
"""

import json
import sys
from pathlib import Path
import numpy as np

# Add project root to path for imports
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status

STEP_NUM = "03"

def safe_json_default(obj):
    if hasattr(obj, 'item'):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

def main():
    print_status(f"STEP {STEP_NUM}: TEP Route-Closure Prediction", "TITLE")
    
    data_path = PROJECT_ROOT / "data" / "raw" / "snh0pe" / "snh0pe_literature_data.json"
    if not data_path.exists():
        print_status(f"Required data file missing: {data_path}", "ERROR")
        
    with open(data_path, 'r') as f:
        data = json.load(f)
        
    print_status("Loading observed delays and lens model magnifications...")
    delays = data["time_delays"]["pierel_2024_photometric"]
    dt_AB_obs = delays["dt_AB"]["value"]
    dt_CB_obs = delays["dt_CB"]["value"]
    
    # Deriving dt_BC and dt_CA for a loop: A -> B -> C -> A
    dt_AB = dt_AB_obs
    dt_BC = -dt_CB_obs
    dt_CA = -(dt_AB_obs - dt_CB_obs)
    
    alpha_tep = 0.05  # Representative coupling magnitude for demonstration
    
    mu_A = data["lens_model"]["mu_A"]["value"]
    mu_B = data["lens_model"]["mu_B"]["value"]
    mu_C = data["lens_model"]["mu_C"]["value"]
    
    print_status("Computing temporal shear factors...")
    Gamma_A = 1 + alpha_tep * np.log10(mu_A)
    Gamma_B = 1 + alpha_tep * np.log10(mu_B)
    Gamma_C = 1 + alpha_tep * np.log10(mu_C)
    
    print_status(f"  Path A (mu={mu_A}): Gamma = {Gamma_A:.4f}")
    print_status(f"  Path B (mu={mu_B}): Gamma = {Gamma_B:.4f}")
    print_status(f"  Path C (mu={mu_C}): Gamma = {Gamma_C:.4f}")
    
    R_closure = (Gamma_A - 1) * dt_AB + (Gamma_B - 1) * dt_BC + (Gamma_C - 1) * dt_CA
    
    print_status(f"  Predicted TEP Closure Residual: {R_closure:.2f} days")
    
    results = {
        "step": STEP_NUM,
        "status": "success",
        "tep_predictions": {
            "Gamma_A": Gamma_A,
            "Gamma_B": Gamma_B,
            "Gamma_C": Gamma_C,
            "residual_days": R_closure,
            "alpha_assumed": alpha_tep
        }
    }
    
    out_dir = PROJECT_ROOT / "results" / "outputs"
    output_path = out_dir / f"step_{STEP_NUM}_tep_closure.json"
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2, default=safe_json_default)
        
    print_status(f"Results saved to {output_path}")
    print_status(f"Step {STEP_NUM} complete.")
        
if __name__ == "__main__":
    main()
