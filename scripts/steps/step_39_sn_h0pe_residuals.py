#!/usr/bin/env python3
"""
TEP-LENS: Step 39 - SN H0pe Blind-Prediction Residual Test

SN H0pe is a triply-imaged Type Ia supernova in the galaxy cluster
PLCK G165.7+67.0 (z_lens=0.351, z_src=1.783). Pierel et al. (2024, ApJ 967, 50)
measured photometric time delays between images A, B, and C from JWST light curves.
Pascale et al. (2025, ApJ, arXiv:2403.18902) published blind lens-model predictions
from seven independent modeling teams, constrained by cluster lensing features
without using measured time delays. Strict blinding protocols were followed,
with some models receiving post-unblinding corrections to time-delay sampling
errors (the underlying mass models remained blind).

This step computes observed-minus-predicted residuals for the two measured
delay pairs (AB, CB), tests the GR null for each pair, and reports an
overall system summary. Because SN H0pe has three images, a loop-closure
TEP test is theoretically possible, but the required geometric delays
from lens models are not available in the published tables. We therefore
use the same single-pair proxy-model residual as Step 38:
    R_TEP ~ alpha * log10(mu_i/mu_j) * dt_GR
which for H0pe's moderate magnification contrast predicts residuals of
order ~1-2 d — small compared to per-model scatter (~5-50 d).
"""

import json
import sys
from pathlib import Path
import numpy as np
from scipy import stats as scipy_stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status
from scripts.utils.tep_config import ALPHA_PROXY, SIGMA_ALPHA_PROXY

STEP_NUM = "39"


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def analyze_delay_pair(models, obs_value, obs_err_plus, obs_err_minus, pair_name):
    """Compute residuals and tests for a single delay pair."""
    obs_err_sym = (obs_err_plus + obs_err_minus) / 2.0
    n_total = len(models)

    results = []
    for m in models:
        delta = obs_value - m["dt_pred"]
        sigma_m = (m["err_plus"] + m["err_minus"]) / 2.0
        sigma_tot = float(np.sqrt(obs_err_sym**2 + sigma_m**2))
        z = delta / sigma_tot

        results.append(
            {
                "name": m["name"],
                "dt_pred_days": m["dt_pred"],
                "err_plus": m["err_plus"],
                "err_minus": m["err_minus"],
                "sigma_model_days": sigma_m,
                "delta_obs_minus_pred_days": float(delta),
                "sigma_total_days": sigma_tot,
                "z_score": float(z),
            }
        )

    deltas = np.array([r["delta_obs_minus_pred_days"] for r in results])
    sigma_tots = np.array([r["sigma_total_days"] for r in results])
    weights = 1.0 / sigma_tots**2
    w_sum = weights.sum()

    R_weighted = float((weights * deltas).sum() / w_sum)
    sigma_R = float(1.0 / np.sqrt(w_sum))
    z_weighted = R_weighted / sigma_R

    R_mean = float(np.mean(deltas))
    R_std = float(np.std(deltas))
    n_positive = int(np.sum(deltas > 0))

    p_binom = float(
        scipy_stats.binomtest(n_positive, n_total, 0.5, alternative="greater").pvalue
    )
    z_binom = float((n_positive - n_total * 0.5) / np.sqrt(n_total * 0.25))

    wilcoxon = scipy_stats.wilcoxon(deltas, alternative="two-sided")
    p_wilcoxon = float(wilcoxon.pvalue)

    chi2_gr = float(np.sum((deltas / sigma_tots) ** 2))

    return {
        "pair": pair_name,
        "observed": {
            "value_days": obs_value,
            "err_plus": obs_err_plus,
            "err_minus": obs_err_minus,
        },
        "n_models": n_total,
        "per_model": results,
        "weighted_mean_residual": {
            "R_obs_days": R_weighted,
            "sigma_R_obs_days": sigma_R,
            "z_score": z_weighted,
        },
        "unweighted": {
            "mean_days": R_mean,
            "std_days": R_std,
            "n_positive": n_positive,
        },
        "binomial": {
            "n_positive": n_positive,
            "n_total": n_total,
            "p_value": p_binom,
            "z_approx": z_binom,
        },
        "wilcoxon": {
            "p_value_two_sided": p_wilcoxon,
            "n_nonzero": int(np.sum(deltas != 0)),
        },
        "chi2_gr_null": chi2_gr,
    }


def sign_label(value, tol=1e-12):
    """Return a stable sign label for directional checks."""
    if value > tol:
        return "positive"
    if value < -tol:
        return "negative"
    return "zero"


def directional_sign_test(result, predicted_residual):
    """Binomial sign test for residuals matching the predicted TEP sign."""
    predicted_sign = sign_label(predicted_residual)
    if predicted_sign == "zero":
        return {
            "predicted_sign": predicted_sign,
            "n_matching": 0,
            "n_total": 0,
            "p_value_one_sided": 1.0,
            "z_approx": 0.0,
            "note": "Predicted residual is numerically zero; sign test is undefined.",
        }

    residuals = np.array(
        [r["delta_obs_minus_pred_days"] for r in result["per_model"]],
        dtype=float,
    )
    nonzero = residuals[residuals != 0]
    if predicted_sign == "positive":
        n_matching = int(np.sum(nonzero > 0))
    else:
        n_matching = int(np.sum(nonzero < 0))
    n_total = int(len(nonzero))

    if n_total == 0:
        p_value = 1.0
        z_approx = 0.0
    else:
        p_value = float(
            scipy_stats.binomtest(
                n_matching, n_total, 0.5, alternative="greater"
            ).pvalue
        )
        z_approx = float((n_matching - n_total * 0.5) / np.sqrt(n_total * 0.25))

    return {
        "predicted_sign": predicted_sign,
        "n_matching": n_matching,
        "n_total": n_total,
        "p_value_one_sided": p_value,
        "z_approx": z_approx,
    }


def central_prediction(models, excluded_model_names=None):
    """Central model delay for TEP prediction, with explicit exclusions."""
    excluded = set(excluded_model_names or [])
    selected = [m for m in models if m["name"] not in excluded]
    if not selected:
        raise ValueError("No models available for central prediction.")
    return {
        "mean_days": float(np.mean([m["dt_pred"] for m in selected])),
        "median_days": float(np.median([m["dt_pred"] for m in selected])),
        "n_models": len(selected),
        "excluded_model_names": sorted(excluded),
        "method": (
            "mean excluding explicitly zero-weighted extreme outliers"
            if excluded
            else "mean over all models"
        ),
    }


def main():
    print_status(
        f"STEP {STEP_NUM}: SN H0pe Blind-Prediction Residual Test", "TITLE"
    )

    # Load alpha from Step 07
    s07_path = (
        PROJECT_ROOT / "results" / "outputs" / "step_07_observed_vs_predicted.json"
    )
    if s07_path.exists():
        with open(s07_path) as f:
            s07 = json.load(f)
        bai = s07.get("bootstrap_alpha_inference", {})
        alpha_mean = float(bai.get("alpha_mean", -0.055))
        alpha_std = float(bai.get("alpha_std", 0.011))
        print_status(
            f"Loaded alpha from Step 07: {alpha_mean:+.4f} +/- {alpha_std:.4f}"
        )
    else:
        alpha_mean = ALPHA_PROXY
        alpha_std = SIGMA_ALPHA_PROXY
        print_status(
            f"Step 07 output not found; falling back to alpha = {ALPHA_PROXY} +/- {SIGMA_ALPHA_PROXY}",
            "WARN",
        )

    # ------------------------------------------------------------------
    # Load observed delays and magnifications from canonical source
    # ------------------------------------------------------------------
    lit_path = PROJECT_ROOT / "data" / "raw" / "snh0pe" / "snh0pe_literature_data.json"
    with open(lit_path) as f:
        lit = json.load(f)

    # Use Pierel+2024 photometric delays as canonical (first measured delays)
    td = lit["time_delays"]["pierel_2024_photometric"]
    obs_ab = float(td["dt_AB"]["value"])
    obs_ab_err_plus = float(td["dt_AB"]["err_plus"])
    obs_ab_err_minus = float(td["dt_AB"]["err_minus"])
    obs_cb = float(td["dt_CB"]["value"])
    obs_cb_err_plus = float(td["dt_CB"]["err_plus"])
    obs_cb_err_minus = float(td["dt_CB"]["err_minus"])

    # Use photometric magnifications for TEP proxy prediction
    phot = lit["photometry"]
    mu_A = float(phot["mu_A"]["value"])
    mu_B = float(phot["mu_B"]["value"])
    mu_C = float(phot["mu_C"]["value"])

    print_status(
        f"Observed Delta_t(A,B) = {obs_ab:+.1f} "
        f"(+{obs_ab_err_plus}/-{obs_ab_err_minus}) d (Pierel+2024)"
    )
    print_status(
        f"Observed Delta_t(C,B) = {obs_cb:+.1f} "
        f"(+{obs_cb_err_plus}/-{obs_cb_err_minus}) d (Pierel+2024)"
    )

    # ------------------------------------------------------------------
    # Blind model predictions (Pascale+2025, Table 1)
    # All seven models were built blind to time-delay measurements.
    # Some models had post-unblinding corrections to sampling (see tablenotes
    # in Pascale+2025); the underlying mass models remained blind.
    # WSLAP+ is an extreme outlier for both delays and was given zero weight
    # in the H0 inference; we include it for completeness.
    # ------------------------------------------------------------------

    models_ab = [
        {"name": "GLAFIC", "dt_pred": -105.18, "err_plus": 5.16, "err_minus": 7.89},
        {"name": "Zitrin-analytic", "dt_pred": -105.52, "err_plus": 5.09, "err_minus": 6.32},
        {"name": "LENSTOOL", "dt_pred": -102.68, "err_plus": 6.47, "err_minus": 6.80},
        {"name": "MARS", "dt_pred": -136.27, "err_plus": 20.23, "err_minus": 20.32},
        {"name": "Chen-2020", "dt_pred": -112.25, "err_plus": 6.43, "err_minus": 6.57},
        {"name": "WSLAP+", "dt_pred": -273.35, "err_plus": 95.30, "err_minus": 95.30},
        {"name": "Zitrin-LTM", "dt_pred": -96.45, "err_plus": 5.41, "err_minus": 5.85},
    ]

    models_cb = [
        {"name": "GLAFIC", "dt_pred": -50.71, "err_plus": 3.37, "err_minus": 5.09},
        {"name": "Zitrin-analytic", "dt_pred": -41.06, "err_plus": 14.27, "err_minus": 13.19},
        {"name": "LENSTOOL", "dt_pred": -54.09, "err_plus": 3.18, "err_minus": 3.76},
        {"name": "MARS", "dt_pred": -63.72, "err_plus": 29.24, "err_minus": 28.46},
        {"name": "Chen-2020", "dt_pred": -53.35, "err_plus": 2.71, "err_minus": 2.99},
        {"name": "WSLAP+", "dt_pred": +342.75, "err_plus": 92.54, "err_minus": 92.54},
        {"name": "Zitrin-LTM", "dt_pred": -27.64, "err_plus": 6.52, "err_minus": 4.80},
    ]

    print_status(
        f"Observed magnifications: mu_A={mu_A}, mu_B={mu_B}, mu_C={mu_C} "
        f"(Pierel+2024 photometry)"
    )

    # ------------------------------------------------------------------
    # Analyse each delay pair
    # ------------------------------------------------------------------
    result_ab = analyze_delay_pair(
        models_ab, obs_ab, obs_ab_err_plus, obs_ab_err_minus, "AB"
    )
    result_cb = analyze_delay_pair(
        models_cb, obs_cb, obs_cb_err_plus, obs_cb_err_minus, "CB"
    )

    # ------------------------------------------------------------------
    # WSLAP+ leave-out sensitivity: WSLAP+ is an extreme outlier
    # (given zero weight in the official H0pe H0 inference). The leave-out
    # subset answers the referee question: does the directional consistency
    # claim survive when the largest outlier is dropped?
    # ------------------------------------------------------------------
    models_ab_no_wslap = [m for m in models_ab if m["name"] != "WSLAP+"]
    models_cb_no_wslap = [m for m in models_cb if m["name"] != "WSLAP+"]
    result_ab_no_wslap = analyze_delay_pair(
        models_ab_no_wslap, obs_ab, obs_ab_err_plus, obs_ab_err_minus, "AB_no_WSLAP"
    )
    result_cb_no_wslap = analyze_delay_pair(
        models_cb_no_wslap, obs_cb, obs_cb_err_plus, obs_cb_err_minus, "CB_no_WSLAP"
    )
    print_status(
        f"\n[WSLAP+ leave-out] AB: n_neg={sum(1 for r in result_ab_no_wslap['per_model'] if r['delta_obs_minus_pred_days'] < 0)}/6, "
        f"R_obs = {result_ab_no_wslap['weighted_mean_residual']['R_obs_days']:+.2f} +/- "
        f"{result_ab_no_wslap['weighted_mean_residual']['sigma_R_obs_days']:.2f} d"
    )
    print_status(
        f"[WSLAP+ leave-out] CB: n_neg={sum(1 for r in result_cb_no_wslap['per_model'] if r['delta_obs_minus_pred_days'] < 0)}/6, "
        f"R_obs = {result_cb_no_wslap['weighted_mean_residual']['R_obs_days']:+.2f} +/- "
        f"{result_cb_no_wslap['weighted_mean_residual']['sigma_R_obs_days']:.2f} d"
    )

    print_status(f"\n--- Delay pair AB ---")
    print_status(
        f"Weighted mean residual = {result_ab['weighted_mean_residual']['R_obs_days']:+.2f} "
        f"+/- {result_ab['weighted_mean_residual']['sigma_R_obs_days']:.2f} d"
    )
    print_status(
        f"Binomial: {result_ab['binomial']['n_positive']}/7 positive, "
        f"p = {result_ab['binomial']['p_value']:.4f}"
    )
    print_status(
        f"Wilcoxon (two-sided): p = {result_ab['wilcoxon']['p_value_two_sided']:.4f}"
    )

    print_status(f"\n--- Delay pair CB ---")
    print_status(
        f"Weighted mean residual = {result_cb['weighted_mean_residual']['R_obs_days']:+.2f} "
        f"+/- {result_cb['weighted_mean_residual']['sigma_R_obs_days']:.2f} d"
    )
    print_status(
        f"Binomial: {result_cb['binomial']['n_positive']}/7 positive, "
        f"p = {result_cb['binomial']['p_value']:.4f}"
    )
    print_status(
        f"Wilcoxon (two-sided): p = {result_cb['wilcoxon']['p_value_two_sided']:.4f}"
    )

    # ------------------------------------------------------------------
    # TEP proxy-model predicted residuals
    # ------------------------------------------------------------------
    # WSLAP+ is documented as an extreme outlier and was given zero weight in
    # the official H0pe H0 inference. Use the WSLAP+-excluded central delay for
    # the proxy prediction, while retaining all-model residuals above.
    central_ab = central_prediction(models_ab, excluded_model_names=["WSLAP+"])
    central_cb = central_prediction(models_cb, excluded_model_names=["WSLAP+"])
    all_model_central_ab = central_prediction(models_ab)
    all_model_central_cb = central_prediction(models_cb)
    mean_pred_ab = central_ab["mean_days"]
    mean_pred_cb = central_cb["mean_days"]

    log_mu_ratio_ab = float(np.log10(mu_A / mu_B))
    log_mu_ratio_cb = float(np.log10(mu_C / mu_B))

    R_tep_ab = alpha_mean * log_mu_ratio_ab * mean_pred_ab
    R_tep_cb = alpha_mean * log_mu_ratio_cb * mean_pred_cb
    R_tep_ab_all_model_mean = (
        alpha_mean * log_mu_ratio_ab * all_model_central_ab["mean_days"]
    )
    R_tep_cb_all_model_mean = (
        alpha_mean * log_mu_ratio_cb * all_model_central_cb["mean_days"]
    )

    direction_ab = directional_sign_test(result_ab, R_tep_ab)
    direction_cb = directional_sign_test(result_cb, R_tep_cb)
    direction_ab_no_wslap = directional_sign_test(result_ab_no_wslap, R_tep_ab)
    direction_cb_no_wslap = directional_sign_test(result_cb_no_wslap, R_tep_cb)

    print_status(f"\nTEP proxy predicted residuals:")
    print_status(
        f"  AB: {R_tep_ab:+.3f} d (log10(mu_A/mu_B)={log_mu_ratio_ab:+.3f}, "
        f"dt_mean={mean_pred_ab:+.1f} d; WSLAP+ excluded)"
    )
    print_status(
        f"  CB: {R_tep_cb:+.3f} d (log10(mu_C/mu_B)={log_mu_ratio_cb:+.3f}, "
        f"dt_mean={mean_pred_cb:+.1f} d; WSLAP+ excluded)"
    )
    print_status(
        f"  All-model mean diagnostic: AB={R_tep_ab_all_model_mean:+.3f} d, "
        f"CB={R_tep_cb_all_model_mean:+.3f} d"
    )
    print_status(
        f"  Directional matches: AB={direction_ab['n_matching']}/{direction_ab['n_total']} "
        f"{direction_ab['predicted_sign']}, CB={direction_cb['n_matching']}/"
        f"{direction_cb['n_total']} {direction_cb['predicted_sign']}"
    )

    # ------------------------------------------------------------------
    # Overall system summary
    # ------------------------------------------------------------------
    # Average the two pairs' weighted mean residuals (they are correlated
    # through shared image B, so this is a summary, not an independent measure).
    R_sys_mean = float(
        np.mean(
            [
                result_ab["weighted_mean_residual"]["R_obs_days"],
                result_cb["weighted_mean_residual"]["R_obs_days"],
            ]
        )
    )

    print_status(f"\nSystem summary:")
    print_status(f"  Mean residual across pairs: {R_sys_mean:+.2f} d")

    # ------------------------------------------------------------------
    # Save outputs
    # ------------------------------------------------------------------
    out = {
        "step": STEP_NUM,
        "system": "SN H0pe",
        "status": "success",
        "description": (
            "Blind-prediction residual test for SN H0pe (PLCK G165.7+67.0). "
            "Seven independent lens models from Pascale+2025 compared against "
            "Pierel+2024 observed delays for two image pairs (AB, CB)."
        ),
        "alpha_prior": {"mean": alpha_mean, "sigma": alpha_std},
        "observed": {
            "dt_AB_days": obs_ab,
            "dt_AB_err_plus": obs_ab_err_plus,
            "dt_AB_err_minus": obs_ab_err_minus,
            "dt_CB_days": obs_cb,
            "dt_CB_err_plus": obs_cb_err_plus,
            "dt_CB_err_minus": obs_cb_err_minus,
            "mu_A": mu_A,
            "mu_B": mu_B,
            "mu_C": mu_C,
            "reference": "Pierel et al. 2024, ApJ 967, 50",
        },
        "delay_pair_AB": result_ab,
        "delay_pair_CB": result_cb,
        "wslap_plus_leave_out": {
            "description": (
                "WSLAP+ is an extreme outlier for both H0pe delays (predicted "
                "-273 d and +343 d versus observed -116 d and -49 d) and was "
                "given zero weight in the official H0pe H0 inference. This "
                "subset analysis recomputes all directional and weighted-mean "
                "statistics with WSLAP+ removed, to test whether the published "
                "directional consistency claim survives removal of the dominant outlier."
            ),
            "delay_pair_AB_no_WSLAP": result_ab_no_wslap,
            "delay_pair_CB_no_WSLAP": result_cb_no_wslap,
            "directional_sign_tests": {
                "AB_no_WSLAP": direction_ab_no_wslap,
                "CB_no_WSLAP": direction_cb_no_wslap,
            },
        },
        "tep_prediction": {
            "R_tep_AB_days": R_tep_ab,
            "R_tep_CB_days": R_tep_cb,
            "alpha_used": alpha_mean,
            "log10_mu_ratio_AB": log_mu_ratio_ab,
            "log10_mu_ratio_CB": log_mu_ratio_cb,
            "dt_mean_pred_AB": mean_pred_ab,
            "dt_mean_pred_CB": mean_pred_cb,
            "central_prediction_AB": central_ab,
            "central_prediction_CB": central_cb,
            "all_model_mean_diagnostic": {
                "R_tep_AB_days": R_tep_ab_all_model_mean,
                "R_tep_CB_days": R_tep_cb_all_model_mean,
                "central_prediction_AB": all_model_central_ab,
                "central_prediction_CB": all_model_central_cb,
                "note": (
                    "All-model mean is retained as a diagnostic only because WSLAP+ "
                    "is an explicitly documented extreme outlier for H0pe and was "
                    "zero-weighted in the official H0 inference."
                ),
            },
            "directional_sign_tests": {
                "AB": direction_ab,
                "CB": direction_cb,
            },
            "note": (
                "Single-pair TEP residuals are sub-day to ~2 d, smaller than "
                "per-model scatter (~5-50 d). Primary central predictions exclude "
                "WSLAP+, which is an explicitly documented extreme outlier and was "
                "zero-weighted in the official H0 inference. No loop-closure test "
                "is performed because geometric delays from lens models are not "
                "published."
            ),
        },
        "system_summary": {
            "mean_residual_across_pairs_days": R_sys_mean,
            "n_delay_pairs": 2,
            "n_models_per_pair": 7,
        },
        "limitations": [
            "Two delay pairs share image B; they are not statistically independent.",
            "TEP predicted residuals are ~1-2 d, far below per-model scatter.",
            "WSLAP+ model is an extreme outlier for both delays (given zero weight in H0 inference).",
            "Some models had post-unblinding corrections to sampling; values used are final corrected.",
            "No loop-closure geometric delays published; single-pair proxy test only.",
        ],
    }

    out_path = (
        PROJECT_ROOT
        / "results"
        / "outputs"
        / f"step_{STEP_NUM}_sn_h0pe_residuals.json"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=safe_json_default)
    print_status(f"\nResults saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
