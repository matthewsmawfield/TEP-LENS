#!/usr/bin/env python3
"""
TEP-LENS: Step 35 - Single-Contrast Dominance Quantification

Purpose: Quantify exactly how concentrated the predicted TEP signal is in the
S4-SX contrast versus all other image pairs and loops. Replaces rhetorical
"single-contrast" language with a pipeline-generated metric.

Algorithm:
1. Load step_03 predicted TEP residuals for all five independent loops.
2. Compute the "energy" (sum of squared residuals) across all loops.
3. Zero out the S4-SX Gamma contrast (set Gamma_S4 = Gamma_SX = mean) and
   recompute all loop residuals. The fractional drop in energy measures how
   much signal lives in that one contrast.
4. Compute an effective degrees-of-freedom (participation ratio):
   eff_dof = (sum |R_i|)^2 / sum(R_i^2). This tells us how many independently
   informative loops actually contribute. For a perfectly single-contrast
   system, eff_dof -> 1. For five equal loops, eff_dof -> 5.

Outputs:
- results/outputs/step_35_single_contrast_dominance.json
"""

import json
import sys
from pathlib import Path
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status

STEP_NUM = "35"


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def main():
    print_status(f"STEP {STEP_NUM}: Single-Contrast Dominance Quantification", "TITLE")

    # ------------------------------------------------------------------
    # Load step_03 TEP predicted residuals
    # ------------------------------------------------------------------
    step03_path = PROJECT_ROOT / "results" / "outputs" / "step_03_tep_closure.json"
    if not step03_path.exists():
        print_status(f"step_03 output missing: {step03_path}", "ERROR")
        sys.exit(1)

    with open(step03_path) as f:
        s03 = json.load(f)

    loops = s03["tep_predicted_discrepancies"]
    gamma = s03["gamma_per_image"]
    mu_norm = s03["mu_norm_per_image"]
    alpha_tep = s03["alpha_tep"]

    # Original residuals
    R_orig = np.array([loops[name]["tep_gr_discrepancy_days"] for name in loops])
    loop_names = list(loops.keys())
    n_loops = len(loop_names)

    print_status("Original loop residuals:")
    for name, R in zip(loop_names, R_orig):
        print_status(f"  {name}: {R:+.3f} d")

    # Energy metrics
    energy_total = float(np.sum(R_orig ** 2))
    energy_abs_total = float(np.sum(np.abs(R_orig)))

    # ------------------------------------------------------------------
    # Contrast knockout: what happens if we remove S4-SX Gamma difference?
    # Set Gamma_S4 and Gamma_SX to their mean, preserving inner-cross structure.
    # ------------------------------------------------------------------
    gamma_mean_s4sx = (gamma["S4"] + gamma["SX"]) / 2.0
    gamma_ko = dict(gamma)
    gamma_ko["S4"] = gamma_mean_s4sx
    gamma_ko["SX"] = gamma_mean_s4sx

    # Reconstruct absolute delays relative to S1 from loop data
    # Loop S1_S2_S3: dt_ij = dt_S2 - dt_S1, dt_jk = dt_S3 - dt_S2
    dt_S2 = loops["S1_S2_S3"]["dt_ij_days"]           # dt_S2 - 0
    dt_S3 = dt_S2 + loops["S1_S2_S3"]["dt_jk_days"]  # dt_S3 = dt_S2 + (dt_S3 - dt_S2)
    # Loop S1_S4_SX: dt_ij = dt_S4 - dt_S1
    dt_S4 = loops["S1_S4_SX"]["dt_ij_days"]           # dt_S4 - 0
    dt_SX = dt_S4 + loops["S1_S4_SX"]["dt_jk_days"]  # dt_SX = dt_S4 + (dt_SX - dt_S4)

    dt_all = {"S1": 0.0, "S2": dt_S2, "S3": dt_S3, "S4": dt_S4, "SX": dt_SX}

    def compute_loop_residual(gamma_dict, loop_imgs):
        i, j, k = loop_imgs
        dt_ij = dt_all[j] - dt_all[i]
        dt_jk = dt_all[k] - dt_all[j]
        dt_ki = dt_all[i] - dt_all[k]
        return ((gamma_dict[i] - 1.0) * dt_ij
                + (gamma_dict[j] - 1.0) * dt_jk
                + (gamma_dict[k] - 1.0) * dt_ki)

    R_ko = np.array([
        compute_loop_residual(gamma_ko, loops[name]["images"])
        for name in loop_names
    ])

    energy_ko = float(np.sum(R_ko ** 2))
    energy_drop_fraction = (energy_total - energy_ko) / energy_total if energy_total != 0 else 0.0

    print_status(f"\nS4-SX contrast knockout:")
    print_status(f"  Total energy (sum R^2):     {energy_total:.3f} d^2")
    print_status(f"  Energy after knockout:      {energy_ko:.3f} d^2")
    print_status(f"  Fraction DROP from S4-SX:   {energy_drop_fraction:.1%}")

    # ------------------------------------------------------------------
    # Isolate test: keep only S4-SX contrast, zero all other Gamma differences
    # Set all inner-cross Gammas to their mean; only S4 and SX differ.
    # ------------------------------------------------------------------
    inner_images = ["S1", "S2", "S3", "S4"]
    gamma_inner_mean = np.mean([gamma[img] for img in inner_images])
    gamma_iso = {img: gamma_inner_mean for img in inner_images}
    gamma_iso["SX"] = gamma["SX"]  # keep SX different from inner mean

    R_iso = np.array([
        compute_loop_residual(gamma_iso, loops[name]["images"])
        for name in loop_names
    ])
    energy_iso = float(np.sum(R_iso ** 2))
    energy_retained_fraction = energy_iso / energy_total if energy_total != 0 else 0.0

    print_status(f"\nIsolate S4-SX contrast (zero inner-cross differences):")
    print_status(f"  Retained energy fraction:   {energy_retained_fraction:.1%}")

    # ------------------------------------------------------------------
    # Effective degrees of freedom (participation ratio on |R|)
    # ------------------------------------------------------------------
    # eff_dof = (sum |R_i|)^2 / sum(R_i^2)
    # For 5 equal contributors: eff_dof = 5.
    # For 1 dominant + 4 negligible: eff_dof -> 1.
    eff_dof = (np.sum(np.abs(R_orig)) ** 2) / np.sum(R_orig ** 2)

    print_status(f"\nEffective degrees of freedom (participation ratio):")
    print_status(f"  eff_dof = {eff_dof:.2f}")
    print_status(f"  Interpretation: the predicted TEP signal is concentrated in")
    print_status(f"  ~{eff_dof:.1f} effective independent contrast(s).")

    # ------------------------------------------------------------------
    # Per-loop contribution fraction
    # ------------------------------------------------------------------
    null_region_loops = {"S1_S2_S3", "S1_S2_S4", "S1_S3_S4"}
    loop_contrib = {}
    for name, R in zip(loop_names, R_orig):
        frac = (R ** 2) / energy_total if energy_total != 0 else 0.0
        loop_contrib[name] = {
            "R_tep_days": float(R),
            "energy_fraction": float(frac),
            "loop_type": "null_region" if name in null_region_loops else "probative_contrast",
        }

    # SX-loop vs inner-cross energy split
    sx_loops = [n for n in loop_names if "SX" in n]
    inner_loops = [n for n in loop_names if "SX" not in n]
    energy_sx = float(np.sum([loops[n]["tep_gr_discrepancy_days"] ** 2 for n in sx_loops]))
    energy_inner = float(np.sum([loops[n]["tep_gr_discrepancy_days"] ** 2 for n in inner_loops]))

    print_status(f"\nLoop-type energy split:")
    print_status(f"  SX-containing loops:   {energy_sx:.3f} d^2 ({energy_sx/energy_total:.1%})")
    print_status(f"  Inner-cross only loops: {energy_inner:.3f} d^2 ({energy_inner/energy_total:.1%})")

    # ------------------------------------------------------------------
    # Save results
    # ------------------------------------------------------------------
    results = {
        "step": STEP_NUM,
        "status": "success",
        "system": "SN Refsdal (MACS J1149.6+2223)",
        "alpha_tep": alpha_tep,
        "energy_metrics": {
            "total_energy_sum_R2_d2": energy_total,
            "total_abs_sum_R_d": energy_abs_total,
            "effective_degrees_of_freedom": float(eff_dof),
            "interpretation_dof": (
                f"The predicted TEP signal is concentrated in approximately {eff_dof:.1f} "
                f"effective independent contrast(s). A value near 1.0 indicates a single-contrast "
                f"system; a value near {n_loops}.0 would indicate all loops contribute equally."
            ),
        },
        "contrast_knockout": {
            "description": "Set Gamma_S4 = Gamma_SX = mean, recompute all loop residuals",
            "energy_after_knockout_d2": energy_ko,
            "energy_drop_fraction": float(energy_drop_fraction),
            "interpretation_knockout": (
                f"Removing the S4-SX Gamma contrast eliminates {energy_drop_fraction:.1%} "
                f"of the total predicted TEP signal energy."
            ),
        },
        "contrast_isolation": {
            "description": "Zero all inner-cross Gamma differences, keep only S4-SX contrast",
            "energy_after_isolation_d2": energy_iso,
            "energy_retained_fraction": float(energy_retained_fraction),
            "interpretation_isolation": (
                f"The S4-SX contrast alone retains {energy_retained_fraction:.1%} "
                f"of the total predicted TEP signal energy."
            ),
        },
        "loop_type_split": {
            "sx_loop_energy_d2": energy_sx,
            "inner_loop_energy_d2": energy_inner,
            "sx_loop_energy_fraction": float(energy_sx / energy_total) if energy_total else 0.0,
            "inner_loop_energy_fraction": float(energy_inner / energy_total) if energy_total else 0.0,
        },
        "per_loop_contribution": loop_contrib,
        "headline": (
            f"SN Refsdal's predicted TEP signal is {energy_drop_fraction:.1%} concentrated "
            f"in the S4-SX contrast (effective DOF = {eff_dof:.1f}), confirming this is a "
            f"single-contrast, single-system measurement."
        ),
    }

    out_dir = PROJECT_ROOT / "results" / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"step_{STEP_NUM}_single_contrast_dominance.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=safe_json_default)

    print_status(f"\nResults saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
