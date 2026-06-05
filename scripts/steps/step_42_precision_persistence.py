#!/usr/bin/env python3
"""
TEP-LENS: Step 42 - Precision-Persistence Test

Purpose: Test the leading conventional explanation for the SN Refsdal SX
under-prediction — that the blind GR models simply had immature/low-precision
mass models, and a sufficiently precise GR model would converge to the observed
376 d with zero residual.

Method: regress the per-model blind-prediction residual (obs - pred) on the
per-model uncertainty sigma_model, and extrapolate to the high-precision limit
(sigma_model -> 0). Two competing hypotheses:

  (immaturity / GR)   : residual -> 0 as sigma_model -> 0  (intercept consistent with 0)
  (irreducible / TEP) : residual -> r0 > 0                 (intercept significantly > 0)

If the highest-precision models retain a positive residual and the sigma->0
intercept is inconsistent with zero, the model-immaturity explanation is
disfavoured: the under-prediction is not an artifact of low blind-model
precision.

HONEST-REPORTING CONTRACT: reports the fitted intercept and its significance
regardless of sign. Caveats (small N, shared lens inputs, intercept is an
extrapolation, post-blind models see SX) are recorded in the output.

Inputs : results/outputs/step_07_observed_vs_predicted.json
Outputs: results/outputs/step_42_precision_persistence.json
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status

STEP_NUM = "42"


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def main():
    print_status(f"STEP {STEP_NUM}: Precision-Persistence Test — SN Refsdal SX", "TITLE")

    s07 = json.load(open(PROJECT_ROOT / "results" / "outputs" / "step_07_observed_vs_predicted.json"))
    models = s07["per_model_results"]

    names = [m["name"] for m in models]
    resid = np.array([m["delta_obs_minus_pred_days"] for m in models], dtype=float)
    sig_model = np.array([m["sigma_model_days"] for m in models], dtype=float)
    blind = np.array([m["blind"] for m in models], dtype=bool)

    print_status("Per-model (sigma_model, residual):")
    for n, s, r, b in zip(names, sig_model, resid, blind):
        print_status(f"  {n:12s} sigma={s:5.1f} d  resid={r:+6.1f} d  blind={b}")

    # Weighted linear regression resid = r0 + slope * sigma_model.
    # Weights = 1/sigma_model^2 (more precise models constrain the intercept more).
    # The intercept r0 is the residual extrapolated to a perfectly precise GR model.
    def weighted_linfit(x, y, w):
        W = np.sum(w)
        mx = np.sum(w * x) / W
        my = np.sum(w * y) / W
        Sxx = np.sum(w * (x - mx) ** 2)
        Sxy = np.sum(w * (x - mx) * (y - my))
        slope = Sxy / Sxx
        intercept = my - slope * mx
        # Parameter covariance (assuming weights are inverse-variance up to a scale).
        n = len(x)
        resid_fit = y - (intercept + slope * x)
        dof = max(n - 2, 1)
        chi2 = np.sum(w * resid_fit ** 2)
        scale = chi2 / dof  # reduced chi2; inflate errors if model misfit
        var_int = scale * (1.0 / W + mx ** 2 / Sxx)
        var_slope = scale * (1.0 / Sxx)
        return (float(intercept), float(np.sqrt(var_int)),
                float(slope), float(np.sqrt(var_slope)),
                float(chi2), int(dof))

    for subset_name, mask in (("all_8", np.ones_like(blind)), ("blind_7", blind)):
        x = sig_model[mask]
        y = resid[mask]
        w = 1.0 / x ** 2  # precision weighting
        intercept, sig_int, slope, sig_slope, chi2, dof = weighted_linfit(x, y, w)
        z_int = intercept / sig_int if sig_int > 0 else float("nan")
        p_int = float(stats.norm.sf(z_int)) if np.isfinite(z_int) else float("nan")

        # Most-precise model residual (direct, model-independent check).
        i_best = int(np.argmin(x))
        best_name = np.array(names)[mask][i_best]
        best_sigma = float(x[i_best]); best_resid = float(y[i_best])

        print_status(f"\n[{subset_name}] precision-weighted fit  resid = r0 + slope*sigma:")
        print_status(f"  intercept r0 (sigma->0) = {intercept:+.2f} +/- {sig_int:.2f} d  "
                     f"(z={z_int:+.2f}, one-sided p={p_int:.3f})")
        print_status(f"  slope                   = {slope:+.3f} +/- {sig_slope:.3f} d per d")
        print_status(f"  most-precise model      = {best_name} (sigma={best_sigma:.0f} d, resid={best_resid:+.1f} d)")

        out_key = subset_name
        globals().setdefault("_fit", {})[out_key] = {
            "intercept_r0_days": intercept,
            "intercept_sigma_days": sig_int,
            "intercept_z": z_int,
            "intercept_p_one_sided": p_int,
            "slope_days_per_day": slope,
            "slope_sigma": sig_slope,
            "chi2": chi2,
            "dof": dof,
            "most_precise_model": {"name": str(best_name),
                                   "sigma_model_days": best_sigma,
                                   "residual_days": best_resid},
        }

    fits = globals()["_fit"]
    r0_blind = fits["blind_7"]["intercept_r0_days"]
    p_blind = fits["blind_7"]["intercept_p_one_sided"]

    if r0_blind > 0 and p_blind < 0.05:
        verdict = ("IMMATURITY DISFAVOURED: the sigma->0 intercept is positive and "
                   "inconsistent with zero, so the SX under-prediction does not vanish in "
                   "the high-precision limit. A perfectly precise GR model of this type "
                   "would still under-predict — consistent with a genuine (TEP-like) "
                   "residual rather than blind-model immaturity.")
    elif r0_blind > 0:
        verdict = ("SUGGESTIVE BUT NOT DECISIVE: the sigma->0 intercept is positive "
                   f"(r0={r0_blind:+.1f} d) but only at p={p_blind:.2f}; the immaturity "
                   "explanation is weakened, not excluded, at current N and precision.")
    else:
        verdict = ("IMMATURITY NOT EXCLUDED: the intercept is consistent with (or below) "
                   "zero, so a high-precision GR model could plausibly remove the residual.")

    print_status("\n" + verdict)

    caveats = [
        "Small N (7 blind / 8 total) and shared lens inputs limit the regression.",
        "The intercept is an extrapolation beyond the data's precision range.",
        "Residual and sigma_model may be mildly correlated across method families.",
        "Grillo+2024 is a post-blind precision update; it still retains a ~+14 d residual, "
        "which is the single strongest data point against the immaturity explanation, but "
        "post-blind models are not independent confirmations.",
        "Schuldt et al. 2026 reach sub-percent statistical precision on dt_SX:S1 but are "
        "post-observation (use SX), so they cannot serve as a blind high-precision anchor; "
        "they bound what precision is achievable, not the blind residual.",
    ]

    headline = (
        f"Precision-persistence (blind 7): sigma->0 residual intercept r0 = "
        f"{r0_blind:+.1f} +/- {fits['blind_7']['intercept_sigma_days']:.1f} d "
        f"(one-sided p={p_blind:.3f}). Most-precise blind model "
        f"({fits['blind_7']['most_precise_model']['name']}, "
        f"sigma={fits['blind_7']['most_precise_model']['sigma_model_days']:.0f} d) retains "
        f"a {fits['blind_7']['most_precise_model']['residual_days']:+.0f} d residual."
    )

    results = {
        "step": STEP_NUM,
        "status": "success",
        "description": (
            "Tests whether the SN Refsdal SX blind under-prediction vanishes in the "
            "high-precision (sigma_model -> 0) limit. A positive intercept disfavours the "
            "model-immaturity (GR) explanation."),
        "fits": fits,
        "verdict": verdict,
        "caveats": caveats,
        "headline": headline,
    }

    out_path = PROJECT_ROOT / "results" / "outputs" / f"step_{STEP_NUM}_precision_persistence.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=safe_json_default)
    print_status(f"\nResults saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
