#!/usr/bin/env python3
"""
TEP-LENS: Step 38 - SN Encore Blind-Prediction Residual Test

SN Encore is a multiply-imaged Type Ia supernova in the galaxy
cluster MACS J0138.2-2155. Pierel et al. (2026, ApJ, arXiv:2509.12301) measured
the time delay between images 1b and 1a from JWST light-curve photometry.
Suyu et al. (2025/2026, A&A, arXiv:2509.12319) published blind lens-model
predictions from seven independent modeling teams (eight total models, with
one team contributing two variants). The models were constrained by cluster
lensing features and SN image positions but NOT by the measured time delays
or magnifications (Pierel et al. 2026, Section 6.2, confirms blinding).

This step computes observed-minus-predicted residuals for the single measured
delay pair (1b,1a) and tests the GR null (residuals symmetric about zero).
Because SN Encore has only two resolved images (1a, 1b), there is no closed
three-image loop for a direct loop-closure TEP test. The TEP proxy-model
predicted residual for a single delay pair scales as
    R_TEP ~ alpha * log10(mu_1/mu_2) * dt_GR
which for Encore's modest magnification contrast (mu_1b/mu_1a ~ 1.5) and
delay baseline (~37 d) yields a predicted shift of less than 1 day — well below
the per-model scatter. SN Encore therefore serves primarily as an independent
consistency check rather than a high-SNR evidence strand.
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

STEP_NUM = "38"


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def main():
    print_status(
        f"STEP {STEP_NUM}: SN Encore Blind-Prediction Residual Test", "TITLE"
    )

    # ------------------------------------------------------------------
    # Load observed delay from canonical catalog
    # ------------------------------------------------------------------
    catalog_path = PROJECT_ROOT / "data" / "raw" / "sn_lensing" / "lensed_sn_catalog.json"
    with open(catalog_path) as f:
        catalog = json.load(f)

    encore = catalog["sn_encore"]
    td = encore["time_delays_days"]["dt_1b_1a"]
    obs_value = float(td["value"])
    obs_err_plus = float(td.get("err_plus", td.get("err", 0.0)))
    obs_err_minus = float(td.get("err_minus", td.get("err", 0.0)))
    obs_err_sym = (obs_err_plus + obs_err_minus) / 2.0

    print_status(
        f"Observed Delta_t(1b,1a) = {obs_value:+.1f} "
        f"(+{obs_err_plus}/-{obs_err_minus}) d ("
        f"{encore['metadata']['notes']})"
    )

    # ------------------------------------------------------------------
    # Blind model predictions (Suyu+2025/2026, Table E.1)
    # All models are pre-blind (constrained before time-delay measurement).
    # Values extracted from arXiv:2509.12319, Appendix E, Table E.1.
    # The Lenstool I and WSLAP+ models show large scatter; we retain them
    # for statistical completeness but flag them as outliers.
    # ------------------------------------------------------------------
    models = [
        {
            "name": "glafic",
            "team": "Oguri",
            "ref": "Suyu+2025, Table E.1 (glafic)",
            "dt_pred": -32.4,
            "err_plus": 2.2,
            "err_minus": 2.5,
            "method": "GLAFIC parametric",
            "blind": True,
        },
        {
            "name": "GLEE",
            "team": "Acebron/Grillo",
            "ref": "Suyu+2025, Table E.1 (GLEE)",
            "dt_pred": -37.1,
            "err_plus": 2.7,
            "err_minus": 2.6,
            "method": "GLEE parametric",
            "blind": True,
        },
        {
            "name": "GLEE-baseline",
            "team": "Acebron/Grillo",
            "ref": "Suyu+2025, Table E.1 (GLEE-baseline)",
            "dt_pred": -36.9,
            "err_plus": 2.6,
            "err_minus": 3.0,
            "method": "GLEE parametric (baseline)",
            "blind": True,
        },
        {
            "name": "Lenstool I",
            "team": "Acebron",
            "ref": "Suyu+2025, Table E.1 (Lenstool I)",
            "dt_pred": -75.0,
            "err_plus": 55.0,
            "err_minus": 54.0,
            "method": "Lenstool parametric",
            "blind": True,
        },
        {
            "name": "Lenstool II",
            "team": "Acebron",
            "ref": "Suyu+2025, Table E.1 (Lenstool II)",
            "dt_pred": -35.6,
            "err_plus": 3.1,
            "err_minus": 4.3,
            "method": "Lenstool parametric (alt. scaling)",
            "blind": True,
        },
        {
            "name": "MrMARTIAN",
            "team": "Furtak/Zitrin",
            "ref": "Suyu+2025, Table E.1 (MrMARTIAN)",
            "dt_pred": -40.6,
            "err_plus": 9.8,
            "err_minus": 2.9,
            "method": "MrMARTIAN free-form",
            "blind": True,
        },
        {
            "name": "WSLAP+",
            "team": "Diego/Sendra",
            "ref": "Suyu+2025, Table E.1 (WSLAP+)",
            "dt_pred": -112.0,
            "err_plus": 32.0,
            "err_minus": 32.0,
            "method": "WSLAP+ hybrid",
            "blind": True,
        },
        {
            "name": "Zitrin-analytic",
            "team": "Zitrin",
            "ref": "Suyu+2025, Table E.1 (Zitrin-analytic)",
            "dt_pred": -40.2,
            "err_plus": 7.5,
            "err_minus": 11.1,
            "method": "Zitrin-LTM analytic",
            "blind": True,
        },
    ]

    n_blind = sum(1 for m in models if m["blind"])
    n_total = len(models)
    print_status(
        f"Loaded {n_blind} blind models ({n_total} total) from Suyu+2025 Table E.1"
    )

    # ------------------------------------------------------------------
    # TEP proxy-model predicted residual for this single delay pair
    # For a two-image delay, the proxy correction is:
    #   R_TEP = alpha * [log10(mu_1a) - log10(mu_1b)] * dt_GR
    # Using observed magnifications from catalog and mean predicted delay.
    # ------------------------------------------------------------------
    mu_1a_obs = float(encore["magnification_proxies"]["mu_absolute"]["1a"]["value"])
    mu_1b_obs = float(encore["magnification_proxies"]["mu_absolute"]["1b"]["value"])
    dt_mean_pred = float(np.mean([m["dt_pred"] for m in models]))

    log_mu_ratio = np.log10(mu_1a_obs / mu_1b_obs)

    # Load alpha from Step 07 bootstrap inference
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

    R_tep_prediction = alpha_mean * log_mu_ratio * dt_mean_pred
    R_tep_unit = R_tep_prediction / alpha_mean if alpha_mean != 0 else 0.0

    print_status(
        f"TEP proxy predicted residual: {R_tep_prediction:+.3f} d "
        f"(alpha={alpha_mean}, log10(mu_ratio)={log_mu_ratio:+.3f}, dt_mean={dt_mean_pred:+.1f} d)"
    )
    print_status(
        f"Sensitivity: {abs(R_tep_unit):.3f} d per unit alpha (small because "
        f"magnification contrast is modest)"
    )

    # ------------------------------------------------------------------
    # Per-model residuals
    # ------------------------------------------------------------------
    print_status("\nPer-model observed - predicted residuals:")
    print_status(
        f"{'Model':<14} {'Pred':>7} {'Delta':>8} {'sigma_m':>8} {'sigma_tot':>10} {'z':>7}"
    )

    results_per_model = []
    for m in models:
        delta = obs_value - m["dt_pred"]
        sigma_m = (m["err_plus"] + m["err_minus"]) / 2.0
        sigma_tot = float(np.sqrt(obs_err_sym**2 + sigma_m**2))
        z = delta / sigma_tot
        alpha_inferred = delta / R_tep_unit if R_tep_unit != 0 else None

        results_per_model.append(
            {
                "name": m["name"],
                "team": m["team"],
                "ref": m["ref"],
                "method": m["method"],
                "blind": m["blind"],
                "dt_pred_days": m["dt_pred"],
                "err_plus": m["err_plus"],
                "err_minus": m["err_minus"],
                "sigma_model_days": sigma_m,
                "delta_obs_minus_pred_days": float(delta),
                "sigma_total_days": sigma_tot,
                "z_score": float(z),
                "alpha_inferred": (
                    float(alpha_inferred) if alpha_inferred is not None else None
                ),
            }
        )

        print_status(
            f"{m['name']:<14} {m['dt_pred']:>+7.1f}  {delta:>+8.1f}  "
            f"{sigma_m:>8.1f}  {sigma_tot:>10.1f}  {z:>+7.2f}"
        )

    # ------------------------------------------------------------------
    # Weighted and unweighted statistics
    # ------------------------------------------------------------------
    deltas = np.array([r["delta_obs_minus_pred_days"] for r in results_per_model])
    sigma_tots = np.array([r["sigma_total_days"] for r in results_per_model])
    weights = 1.0 / sigma_tots**2
    w_sum = weights.sum()

    R_obs_weighted = float((weights * deltas).sum() / w_sum)
    sigma_R_obs = float(1.0 / np.sqrt(w_sum))
    z_weighted = R_obs_weighted / sigma_R_obs

    R_obs_mean = float(np.mean(deltas))
    R_obs_std = float(np.std(deltas))
    n_positive = int(np.sum(deltas > 0))
    n_negative = int(np.sum(deltas < 0))

    print_status(
        f"\nWeighted mean residual R_obs = {R_obs_weighted:+.2f} "
        f"+/- {sigma_R_obs:.2f} d  (z = {z_weighted:+.2f})"
    )
    print_status(
        f"Unweighted mean: {R_obs_mean:+.2f} d, std = {R_obs_std:.2f} d"
    )
    print_status(f"Positive residuals: {n_positive}/{n_total}, negative: {n_negative}/{n_total}")

    # ------------------------------------------------------------------
    # TEP consistency test
    # ------------------------------------------------------------------
    tep_residual = R_obs_weighted - R_tep_prediction
    z_tep = tep_residual / sigma_R_obs
    alpha_inferred_wmean = R_obs_weighted / R_tep_unit if R_tep_unit != 0 else None

    print_status(f"\nTEP consistency test:")
    print_status(
        f"  R_obs (weighted)  = {R_obs_weighted:+.2f} +/- {sigma_R_obs:.2f} d"
    )
    print_status(f"  R_TEP_pred        = {R_tep_prediction:+.3f} d")
    print_status(f"  Tension           = {z_tep:+.2f} sigma")
    if alpha_inferred_wmean is not None:
        print_status(
            f"  Inferred alpha    = {alpha_inferred_wmean:.4f} "
            f"(if R_obs attributed entirely to TEP; highly uncertain)"
        )

    # GR null test
    z_gr = R_obs_weighted / sigma_R_obs
    print_status(f"\nGR null test:")
    print_status(
        f"  R_obs = {R_obs_weighted:+.2f} +/- {sigma_R_obs:.2f} d"
    )
    print_status(f"  Tension with GR (R=0): {z_gr:+.2f} sigma")

    # ------------------------------------------------------------------
    # Binomial sign test
    # ------------------------------------------------------------------
    p_binom = float(
        scipy_stats.binomtest(n_positive, n_total, 0.5, alternative="greater").pvalue
    )
    z_binom_approx = float((n_positive - n_total * 0.5) / np.sqrt(n_total * 0.25))
    print_status(f"\nBinomial sign test (all {n_total} models):")
    print_status(f"  Positive residuals: {n_positive}/{n_total}")
    print_status(f"  One-sided p (H0: p=0.5): p = {p_binom:.4f}")
    print_status(f"  Equivalent z: z = {z_binom_approx:+.2f}")

    # Wilcoxon signed-rank test
    wilcoxon_result = scipy_stats.wilcoxon(deltas, alternative="two-sided")
    p_wilcoxon = float(wilcoxon_result.pvalue)
    n_nonzero = int(np.sum(deltas != 0))
    print_status(
        f"\nWilcoxon signed-rank (all {n_total}, two-sided): p = {p_wilcoxon:.6f}"
    )

    # ------------------------------------------------------------------
    # Chi-squared model comparison: GR (R=0) vs TEP (R=R_tep_prediction)
    # ------------------------------------------------------------------
    chi2_gr = float(np.sum(((deltas - 0) / sigma_tots) ** 2))
    chi2_tep = float(
        np.sum(((deltas - R_tep_prediction) / sigma_tots) ** 2)
    )
    delta_chi2 = chi2_gr - chi2_tep

    # Both are fixed predictions (0 free params), not nested models.
    # Report individual goodness-of-fit p-values; delta_chi2 has no formal p-value.
    p_chi2_gr = float(scipy_stats.chi2.sf(chi2_gr, df=n_total))
    p_chi2_tep = float(scipy_stats.chi2.sf(chi2_tep, df=n_total))

    print_status(f"\nChi-squared model comparison ({n_total} models):")
    print_status(f"  chi^2 under GR  (R=0):            {chi2_gr:.3f}  (p={p_chi2_gr:.3f})")
    print_status(
        f"  chi^2 under TEP (R={R_tep_prediction:.2f} d):  {chi2_tep:.3f}  (p={p_chi2_tep:.3f})"
    )
    print_status(f"  Delta chi^2 = {delta_chi2:+.3f}")
    print_status(f"  NOTE: Both are fixed predictions (0 free params); delta_chi2 has no formal p-value.")

    # ------------------------------------------------------------------
    # Save outputs
    # ------------------------------------------------------------------
    out = {
        "step": STEP_NUM,
        "system": "SN Encore",
        "status": "success",
        "description": (
            "Blind-prediction residual test for SN Encore (MACS J0138-2155). "
            "Eight independent lens models from Suyu+2025 compared against "
            "Pierel+2026 observed delay."
        ),
        "observed": {
            "dt_1b_1a_days": obs_value,
            "err_plus": obs_err_plus,
            "err_minus": obs_err_minus,
            "reference": "Pierel et al. 2026, ApJ, Table 5",
        },
        "magnifications_observed": {
            "mu_1a": mu_1a_obs,
            "mu_1b": mu_1b_obs,
            "reference": "Pierel et al. 2026",
        },
        "models": results_per_model,
        "weighted_mean_residual": {
            "R_obs_days": R_obs_weighted,
            "sigma_R_obs_days": sigma_R_obs,
            "z_score": z_weighted,
            "alpha_inferred": (
                float(alpha_inferred_wmean) if alpha_inferred_wmean is not None else None
            ),
        },
        "unweighted_statistics": {
            "mean_days": R_obs_mean,
            "std_days": R_obs_std,
            "n_positive": n_positive,
            "n_negative": n_negative,
        },
        "tep_prediction": {
            "R_tep_prediction_days": R_tep_prediction,
            "alpha_used": alpha_mean,
            "log10_mu_ratio": float(log_mu_ratio),
            "dt_mean_pred_days": float(dt_mean_pred),
            "R_tep_unit_days_per_alpha": float(R_tep_unit),
            "note": (
                "Single-pair TEP residual is ~<1 d because magnification contrast "
                "is modest. No loop-closure test possible with only two images."
            ),
        },
        "gr_null_test": {
            "z_score": z_gr,
            "p_value_two_sided": float(2 * scipy_stats.norm.sf(abs(z_gr))),
        },
        "binomial_sign_test": {
            "n_positive": n_positive,
            "n_total": n_total,
            "p_value_one_sided": p_binom,
            "z_approx": z_binom_approx,
        },
        "wilcoxon_test": {
            "p_value_two_sided": p_wilcoxon,
            "n_nonzero": n_nonzero,
        },
        "chi2_comparison": {
            "chi2_gr": chi2_gr,
            "chi2_tep": chi2_tep,
            "delta_chi2": delta_chi2,
            "p_value_chi2_gr": p_chi2_gr,
            "p_value_chi2_tep": p_chi2_tep,
            "note": "Both GR and TEP are fixed predictions (0 free params), not nested models. Delta chi^2 has no formal p-value.",
        },
        "limitations": [
            "Only two resolved images (1a, 1b); no three-image loop for direct TEP closure test.",
            "TEP predicted residual is <1 d, far smaller than per-model scatter (~3-50 d).",
            "System serves as independent consistency check, not high-SNR evidence strand.",
            "Models with very large uncertainties (Lenstool I, WSLAP+) dominate unweighted scatter but are downweighted in the precision-weighted mean.",
        ],
    }

    out_path = (
        PROJECT_ROOT
        / "results"
        / "outputs"
        / f"step_{STEP_NUM}_sn_encore_residuals.json"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=safe_json_default)
    print_status(f"\nResults saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
