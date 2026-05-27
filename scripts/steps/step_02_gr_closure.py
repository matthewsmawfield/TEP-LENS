#!/usr/bin/env python3
"""
TEP-LENS: Step 02 - General Relativity Algebraic Loop Sum (SN Refsdal)

With 5 images (S1-S4, SX) and 4 measured pairwise delays relative to S1,
we can construct multiple independent 3-image algebraic loops.

Under any theory with globally assignable arrival times per image,
any closed loop sum is identically zero:
  L_loop(i,j,k) = dt_ij + dt_jk + dt_ki = 0

This is a purely algebraic identity, not a testable observable. The genuine
TEP test is the blind-prediction residual (Step 07), not a closure violation.

Loops constructed here from the 4 measured Kelly+2023 delays:
  - Loop A: S1 -> S2 -> S3 -> S1  (inner Einstein cross triplet)
  - Loop B: S1 -> S2 -> S4 -> S1  (inner cross triplet)
  - Loop C: S1 -> S3 -> S4 -> S1  (inner cross triplet)
  - Loop D: S1 -> S2 -> SX -> S1  (cross-to-arc loop)
  - Loop E: S1 -> S4 -> SX -> S1  (cross-to-arc loop)
"""

import json
import sys
from pathlib import Path
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status

STEP_NUM = "02"

def safe_json_default(obj):
    if hasattr(obj, 'item'):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

def loop_sum(dt_ij, dt_jk, dt_ki):
    """Algebraic loop sum. Identically zero under any theory with globally assignable arrival times."""
    return dt_ij + dt_jk + dt_ki

def propagate_error(*errs):
    """Quadrature error propagation."""
    return float(np.sqrt(sum(e**2 for e in errs)))

def main():
    print_status(f"STEP {STEP_NUM}: GR Algebraic Loop Sum — SN Refsdal", "TITLE")

    catalog_path = PROJECT_ROOT / "data" / "raw" / "sn_lensing" / "lensed_sn_catalog.json"
    if not catalog_path.exists():
        print_status(f"Required catalog missing: {catalog_path}. Run step_01 first.", "ERROR")

    with open(catalog_path, 'r') as f:
        catalog = json.load(f)

    print_status("Extracting SN Refsdal time delays (Kelly et al. 2023, ApJ 948, 93)...")
    delays = catalog["sn_refsdal"]["time_delays_days"]

    # Absolute arrival times relative to S1 (positive = arrives later than S1)
    dt = {
        "S1": 0.0,
        "S2": delays["dt_S2_S1"]["value"],
        "S3": delays["dt_S3_S1"]["value"],
        "S4": delays["dt_S4_S1"]["value"],
        "SX": delays["dt_SX_S1"]["value"],
    }
    err = {
        "S1": 0.0,
        "S2": delays["dt_S2_S1"]["err"],
        "S3": delays["dt_S3_S1"]["err"],
        "S4": delays["dt_S4_S1"]["err"],
        "SX": delays["dt_SX_S1"]["err"],
    }

    for img, val in dt.items():
        if img != "S1":
            print_status(f"  dt_{img}_S1 = {val:+.1f} ± {err[img]:.1f} days")

    # ------------------------------------------------------------------
    # Build closure loops: dt_ij = dt_j - dt_i (time of j minus time of i)
    # Loop closure:  dt_ij + dt_jk + dt_ki = (dt_j-dt_i)+(dt_k-dt_j)+(dt_i-dt_k) = 0 exactly
    # ------------------------------------------------------------------
    def dt_ij(i, j):
        return dt[j] - dt[i], propagate_error(err[i], err[j])

    loops = {
        "S1_S2_S3": ["S1", "S2", "S3"],
        "S1_S2_S4": ["S1", "S2", "S4"],
        "S1_S3_S4": ["S1", "S3", "S4"],
        "S1_S2_SX": ["S1", "S2", "SX"],
        "S1_S4_SX": ["S1", "S4", "SX"],
    }

    loop_results = {}
    print_status("Computing GR algebraic loop sums for all loops...")

    for name, (i, j, k) in loops.items():
        d_ij, e_ij = dt_ij(i, j)
        d_jk, e_jk = dt_ij(j, k)
        d_ki, e_ki = dt_ij(k, i)
        R = loop_sum(d_ij, d_jk, d_ki)
        R_err = propagate_error(e_ij, e_jk, e_ki)
        loop_results[name] = {
            "images": [i, j, k],
            "dt_ij": d_ij, "dt_jk": d_jk, "dt_ki": d_ki,
            "loop_sum_gr": R,
            "loop_sum_err": R_err,
            "gr_null_description": "Identically 0 under GR by construction."
        }
        print_status(f"  Loop {i}-{j}-{k}: L_GR = {R:+.4f} days (identically zero by arithmetic)")

    print_status("GR null test: all loops confirm L_GR = 0.000 days (exact).")

    results = {
        "step": STEP_NUM,
        "status": "success",
        "system": "SN Refsdal (MACS J1149.6+2223)",
        "reference": "Kelly et al. 2023, ApJ 948, 93",
        "gr_algebraic_loop_loops": loop_results,
        "summary": (
            "Under any theory with globally assignable arrival times, the algebraic "
            "loop sum is identically zero for any combination of measured delays. "
            "Five independent loops are computed from the four Kelly+2023 measured delays. "
            "The measured delay errors propagate to give the observational sensitivity for each loop."
        )
    }

    out_dir = PROJECT_ROOT / "results" / "outputs"
    output_path = out_dir / f"step_{STEP_NUM}_gr_closure.json"
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2, default=safe_json_default)

    print_status(f"Results saved to {output_path}")
    print_status(f"Step {STEP_NUM} complete.")

if __name__ == "__main__":
    main()
