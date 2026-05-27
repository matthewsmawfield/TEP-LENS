#!/usr/bin/env python3
"""
TEP-LENS: Step 36 - Forward-Prediction Falsification Thresholds for SN 2025wny

Purpose: Compute pre-registered falsification thresholds for SN 2025wny,
the next resolved multiply-imaged supernova. The proxy model predicts a
non-zero blind-prediction residual for any system where images sample
different potential depths. This script quantifies:

1. Predicted residual range given plausible magnification and geometry priors.
2. Required delay measurement precision for 3-sigma and 5-sigma tests.
3. Explicit falsification condition: what observed residual would reject
   the linear log-magnification ansatz at 95% CL?

These thresholds are pre-registered in the manuscript (§4.9) and are
computed deterministically by this script; they will not be adjusted
post-measurement.

Assumptions (from Johansson et al. 2025, ApJ 995, L17):
- Four images (A-D) in an Einstein-cross geometry.
- Brightest image magnification: mu ~ 15-80 (discovery-paper estimate 20-50,
  but we broaden to 15-80 to be conservative).
- Faintest image magnification: mu ~ 2-15 (consistent with quad-lens flux ratios).
- Longest pairwise delay baseline: 10-150 days (galaxy-scale lens at z_l=0.375).
- Proxy coupling: alpha_lens = -0.055 +/- 0.044 (empirically determined from
  SN Refsdal; propagated as a Gaussian prior).

Algorithm:
1. Monte Carlo over magnification priors, delay baselines, and alpha uncertainty.
2. For each draw, compute the proxy-model predicted residual for the
   brightest-faintest contrast loop.
3. Report percentile envelopes and required precisions.
4. Output falsification threshold: the residual magnitude below which the
   proxy model is excluded at 95% CL for a given measurement precision.
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status, safe_json_default

STEP_NUM = "36"




def main():
    print_status(
        f"STEP {STEP_NUM}: Forward-Prediction Falsification Thresholds — SN 2025wny",
        "TITLE",
    )

    rng = np.random.default_rng(42)
    n_mc = 100000

    # ------------------------------------------------------------------
    # Priors (physically motivated, conservative)
    # ------------------------------------------------------------------
    # Brightest image magnification: discovery paper estimates 20-50;
    # we broaden to 15-80 to include systematic uncertainty.
    mu_bright = rng.uniform(15.0, 80.0, size=n_mc)
    # Faintest image magnification: typical quad-lens flux ratios give ~2-15.
    mu_faint = rng.uniform(2.0, 15.0, size=n_mc)
    # Intermediate images: draw between faint and bright.
    mu_mid1 = mu_faint + rng.uniform(0.0, 1.0, size=n_mc) * (mu_bright - mu_faint)
    mu_mid2 = mu_faint + rng.uniform(0.0, 1.0, size=n_mc) * (mu_bright - mu_faint)

    # Mean magnification across the four-image system
    mu_mean = (mu_bright + mu_mid1 + mu_mid2 + mu_faint) / 4.0

    # Longest delay baseline: galaxy-scale quad lenses typically 10-150 days.
    # SN 2025wny is at z_l = 0.375, so delays are shorter than cluster-scale Refsdal.
    dt_baseline = rng.uniform(10.0, 150.0, size=n_mc)

    # Empirical coupling with uncertainty from SN Refsdal
    alpha = rng.normal(loc=-0.055, scale=0.044, size=n_mc)

    # ------------------------------------------------------------------
    # Compute proxy-model predicted residual for brightest-faintest loop
    # ------------------------------------------------------------------
    # Gamma_t = 1 + alpha * log10(mu / mu_mean)
    gamma_bright = 1.0 + alpha * np.log10(mu_bright / mu_mean)
    gamma_faint = 1.0 + alpha * np.log10(mu_faint / mu_mean)

    # The fractional TEP shift between bright and faint images
    delta_gamma = gamma_bright - gamma_faint

    # Predicted residual: R = delta_gamma * dt_baseline (approximate loop form)
    # For a 4-image loop with two intermediate images, the exact loop residual
    # is smaller than the simple product because the intermediate images
    # partially cancel. We apply a geometric correction factor drawn from
    # the distribution of typical quad-lens loop structures.
    loop_geometry_factor = rng.uniform(0.6, 1.0, size=n_mc)
    R_pred = delta_gamma * dt_baseline * loop_geometry_factor

    # The residual is predicted to be positive (faint image late, bright image early)
    # under the measured negative alpha. We report the positive envelope.
    # For alpha draws that cross zero, the sign flips; we clip to physical range.
    R_pred = np.where(alpha < 0, np.abs(R_pred), -np.abs(R_pred))
    # We care about the magnitude for falsification, so take absolute values
    # for the threshold computation, but keep sign information for prediction.
    R_mag = np.abs(R_pred)

    # ------------------------------------------------------------------
    # Summary statistics
    # ------------------------------------------------------------------
    R_p16 = float(np.percentile(R_mag, 16))
    R_p50 = float(np.percentile(R_mag, 50))
    R_p84 = float(np.percentile(R_mag, 84))
    R_p95 = float(np.percentile(R_mag, 95))
    R_p99 = float(np.percentile(R_mag, 99))

    # Fraction of draws with predicted residual > 1 day (detectable threshold)
    frac_gt_1d = float(np.mean(R_mag > 1.0))
    frac_gt_5d = float(np.mean(R_mag > 5.0))

    print_status(f"Predicted residual magnitude distribution (N={n_mc}):")
    print_status(f"  16th percentile: {R_p16:.2f} d")
    print_status(f"  50th percentile: {R_p50:.2f} d")
    print_status(f"  84th percentile: {R_p84:.2f} d")
    print_status(f"  95th percentile: {R_p95:.2f} d")
    print_status(f"  99th percentile: {R_p99:.2f} d")
    print_status(f"  P(R > 1 d): {frac_gt_1d:.2%}")
    print_status(f"  P(R > 5 d): {frac_gt_5d:.2%}")

    # ------------------------------------------------------------------
    # Required measurement precision
    # ------------------------------------------------------------------
    # For a z-test: z = R_pred / sigma_dt
    # Required sigma_dt = R_pred / z_target
    z_3sigma = 3.0
    z_5sigma = 5.0

    sigma_req_3sigma = R_p50 / z_3sigma
    sigma_req_5sigma = R_p50 / z_5sigma

    # Conservative: use 16th percentile residual (smaller signal)
    sigma_req_3sigma_conservative = R_p16 / z_3sigma
    sigma_req_5sigma_conservative = R_p16 / z_5sigma

    print_status(f"\nRequired delay precision (longest baseline):")
    print_status(f"  For 3-sigma test (median signal): sigma_dt < {sigma_req_3sigma:.2f} d")
    print_status(f"  For 5-sigma test (median signal): sigma_dt < {sigma_req_5sigma:.2f} d")
    print_status(f"  For 3-sigma test (conservative, 16th pct): sigma_dt < {sigma_req_3sigma_conservative:.2f} d")
    print_status(f"  For 5-sigma test (conservative, 16th pct): sigma_dt < {sigma_req_5sigma_conservative:.2f} d")

    # ------------------------------------------------------------------
    # Falsification threshold
    # ------------------------------------------------------------------
    # If observed |R_obs| < threshold at measurement precision sigma_dt,
    # exclude proxy model at 95% CL.
    # For a 2-sigma exclusion: threshold = 2 * sigma_dt
    # We report the threshold as a function of sigma_dt.
    sigma_dt_grid = np.array([0.5, 1.0, 1.5, 2.0, 3.0, 5.0])
    falsification_thresholds = []
    for sdt in sigma_dt_grid:
        # 2-sigma exclusion: |R_obs| < 2*sigma_dt excludes the model
        # if the predicted residual is larger than this.
        # Fraction of MC draws excluded at this precision:
        excluded_frac = float(np.mean(R_mag < 2.0 * sdt))
        falsification_thresholds.append({
            "sigma_dt_days": float(sdt),
            "exclusion_threshold_2sigma_days": float(2.0 * sdt),
            "fraction_of_prior_predictions_excluded": excluded_frac,
            "interpretation": (
                f"At sigma_dt={sdt:.1f}d, a residual |R_obs|<{2.0*sdt:.1f}d "
                f"excludes {excluded_frac:.1%} of the prior prediction envelope."
            ),
        })

    # Pre-registered falsification condition for manuscript §4.9
    # Use sigma_dt = 2.0 d as the reference precision mentioned in the text.
    sigma_ref = 2.0
    threshold_ref = 2.0 * sigma_ref  # 2-sigma
    excluded_frac_ref = float(np.mean(R_mag < threshold_ref))

    print_status(f"\nPre-registered falsification condition (sigma_dt = {sigma_ref} d):")
    print_status(f"  If |R_obs| < {threshold_ref:.1f} d, the proxy model is excluded at 95% CL.")
    print_status(f"  This threshold excludes {excluded_frac_ref:.1%} of the prior prediction envelope.")

    # ------------------------------------------------------------------
    # Save results
    # ------------------------------------------------------------------
    results = {
        "step": STEP_NUM,
        "status": "success",
        "system": "SN 2025wny (z_s=2.011, z_l=0.375)",
        "priors": {
            "mu_bright_range": [15.0, 80.0],
            "mu_faint_range": [2.0, 15.0],
            "dt_baseline_range_days": [10.0, 150.0],
            "alpha_prior": {"mean": -0.055, "sigma": 0.044},
            "loop_geometry_factor_range": [0.6, 1.0],
            "n_monte_carlo": n_mc,
            "seed": 42,
        },
        "predicted_residual_distribution_days": {
            "p16": R_p16,
            "p50": R_p50,
            "p84": R_p84,
            "p95": R_p95,
            "p99": R_p99,
            "fraction_gt_1d": frac_gt_1d,
            "fraction_gt_5d": frac_gt_5d,
        },
        "required_precision_days": {
            "for_3sigma_median_signal": sigma_req_3sigma,
            "for_5sigma_median_signal": sigma_req_5sigma,
            "for_3sigma_conservative_16th_pct": sigma_req_3sigma_conservative,
            "for_5sigma_conservative_16th_pct": sigma_req_5sigma_conservative,
        },
        "falsification_thresholds_vs_precision": falsification_thresholds,
        "pre_registered_falsification_condition": {
            "sigma_dt_reference_days": sigma_ref,
            "exclusion_threshold_2sigma_days": threshold_ref,
            "fraction_of_prior_excluded": excluded_frac_ref,
            "statement": (
                f"If an independent blind-prediction residual for the longest-baseline "
                f"loop in SN 2025wny is consistent with zero at the 2-sigma level "
                f"(|R_obs| < {threshold_ref:.1f} d) and the delay precision satisfies "
                f"sigma_dt < {sigma_ref:.1f} d, the linear log-magnification ansatz "
                f"is excluded at 95% confidence for that system geometry."
            ),
        },
        "headline": (
            f"SN 2025wny predicted proxy-model residual: {R_p50:.1f} "
            f"[{R_p16:.1f}, {R_p84:.1f}] days (16th-84th percentile). "
            f"Required precision for 3-sigma test: sigma_dt < {sigma_req_3sigma:.1f} days. "
            f"Falsification threshold (sigma_dt={sigma_ref}d): |R_obs| < {threshold_ref:.1f}d."
        ),
    }

    out_dir = PROJECT_ROOT / "results" / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"step_{STEP_NUM}_falsification_forward_prediction.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=safe_json_default)

    print_status(f"\nResults saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
