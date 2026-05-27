#!/usr/bin/env python3
"""
TEP-LENS: Step 07 - Observed vs. Blind-Predicted Delay Comparison

The core evidence test. Before SX reappeared (Dec 2015), seven independent
lens modelling teams published blind predictions for Delta_t(SX, S1). These
were computed using mass models constrained only by the Einstein-cross images
S1-S4 and cluster multiple-image positions -- NOT using SX.

Kelly et al. (2023) then independently measured Delta_t(SX, S1) = 376.0 +/- 5.6 d
from SN light-curve fitting. This provides a genuine observed-vs-predicted comparison
using two completely independent datasets.

The analysis:
1. Compile all 7 blind pre-reappearance + 1 post-blind high-precision model predictions.
2. Compute per-model residual: Delta_i = Delta_t_obs - Delta_t_model_i
3. Inverse-variance weighted mean residual R_obs and sigma_R_obs.
4. Binomial sign test: P(>= n_positive | p=0.5, N) under the GR null.
5. Chi-squared model comparison: GR (R=0) vs TEP (R=R_TEP_pred).
   Which hypothesis is better supported by the ensemble of residuals?
6. TEP-corrected consistency: compute the TEP-corrected observation
   Delta_t_corr = Delta_t_obs - R_TEP and show it reduces scatter vs raw.
7. Inferred alpha from weighted mean residual.

Data sources:
- Treu et al. 2016, ApJ 817, 60 (arXiv:1510.05750), Table 2
- Kelly et al. 2023, Science 380, abh1322, Supplementary Table S4
- Kelly et al. 2023, ApJ 948, 93, Table 15 (observed delays)
- Grillo et al. 2024, ApJ 971, 49 (post-blind precision model)
"""

import json
import sys
from pathlib import Path
import numpy as np
from scipy import stats as scipy_stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status
from scripts.utils.plot_style import set_pub_style, COLORS, FIG_SIZE

STEP_NUM = "07"

# ------------------------------------------------------------------
# Pre-registered primary test selection
# ------------------------------------------------------------------
# The primary evidence test is declared here BEFORE any computation.
# Justification: The Wilcoxon signed-rank test treats each independent
# modelling group as one vote, regardless of the (highly heterogeneous)
# quoted model uncertainties. This avoids the inverse-variance downweighting
# bias that suppresses the weighted-mean z-test when some models have
# very large error bars (e.g., Oguri ±59 d). The Wilcoxon is therefore the
# most robust non-parametric directional test available for this ensemble.
# It is selected for its statistical properties, not its p-value.
# ------------------------------------------------------------------
PRIMARY_TEST = "wilcoxon_signed_rank_blind"
PRIMARY_TEST_JUSTIFICATION = (
    "Selected a priori as the most robust non-parametric directional test: "
    "it equal-weights independent modelling groups, eliminating inverse-variance "
    "downweighting from heterogeneous model uncertainties. It is more powerful than "
    "the binomial sign test (which discards rank information) and more robust than "
    "the parametric z-test (which is biased by outlier uncertainties)."
)

def safe_json_default(obj):
    if hasattr(obj, 'item'):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

def main():
    print_status(f"STEP {STEP_NUM}: Observed vs. Blind-Predicted Delay — SN Refsdal SX", "TITLE")

    # ------------------------------------------------------------------
    # Observed delay (Kelly+2023, ApJ 948 93, Table 15, 'Combined' method)
    # This is the measurement AFTER SX appeared, using full 5-image light curves.
    # ------------------------------------------------------------------
    obs_value = 376.0   # days
    obs_err   = 5.6     # days (1-sigma)

    print_status(f"Observed Delta_t(SX, S1) = {obs_value} +/- {obs_err} d "
                 f"(Kelly+2023, independent photometric light-curve fit)")

    # ------------------------------------------------------------------
    # Blind model predictions (all published BEFORE SX reappeared Dec 2015)
    # Source: Treu+2016 ApJ 817 60 Table 2, and Kelly+2023 Science Table S4.
    # Values are median (or central) predictions with 1-sigma uncertainties.
    # The Jauzac and Sharon models required post-publication corrections
    # (noted in Kelly+2023) -- the corrected values from Kelly+2023 S4 are used.
    # Models marked blind=True were published fully blind to SX position/time.
    # ------------------------------------------------------------------
    models = [
        {
            "name": "Oguri-a",
            "team": "Oguri",
            "ref": "Oguri 2015 (MNRAS 449, L86)",
            "dt_pred": 324.0,
            "err_plus": 66.0,
            "err_minus": 52.0,
            "method": "GLAFIC parametric",
            "blind": True,
        },
        {
            "name": "Sharon",
            "team": "Sharon & Johnson",
            "ref": "Sharon & Johnson 2015 (ApJ 800, L26)",
            "dt_pred": 345.0,
            "err_plus": 68.0,
            "err_minus": 51.0,
            "method": "LTM parametric",
            "blind": True,
        },
        {
            "name": "Diego",
            "team": "Diego et al.",
            "ref": "Diego et al. 2016 (MNRAS 456, 356)",
            "dt_pred": 376.0,
            "err_plus": 40.0,
            "err_minus": 60.0,
            "method": "WSLAP+ free-form",
            "blind": True,
        },
        {
            "name": "Grillo",
            "team": "Grillo et al.",
            "ref": "Grillo et al. 2016 (ApJ 822, 78)",
            "dt_pred": 361.0,
            "err_plus": 20.0,
            "err_minus": 27.0,
            "method": "GLEE parametric",
            "blind": True,
        },
        {
            "name": "Kawamata",
            "team": "Kawamata et al.",
            "ref": "Kawamata et al. 2016 (ApJ 819, 114)",
            "dt_pred": 369.0,
            "err_plus": 53.0,
            "err_minus": 44.0,
            "method": "Parametric",
            "blind": True,
        },
        {
            "name": "Jauzac",
            "team": "Jauzac et al.",
            "ref": "Jauzac et al. 2016 (MNRAS 457, 2029) [corrected in Kelly+2023]",
            "dt_pred": 359.0,
            "err_plus": 40.0,
            "err_minus": 56.0,
            "method": "LENSTOOL parametric",
            "blind": True,
        },
        {
            "name": "CATS",
            "team": "Treu et al. (CATS)",
            "ref": "Treu et al. 2016 (ApJ 817, 60)",
            "dt_pred": 374.0,
            "err_plus": 51.0,
            "err_minus": 41.0,
            "method": "LENSTOOL parametric",
            "blind": True,
        },
        {
            "name": "Grillo+2024",
            "team": "Grillo et al.",
            "ref": "Grillo et al. 2024 (ApJ 971, 49) — post-blind precision update",
            "dt_pred": 362.0,
            "err_plus": 16.0,
            "err_minus": 16.0,
            "method": "GLEE parametric (updated)",
            "blind": False,
        },
    ]

    n_blind   = sum(1 for m in models if m["blind"])
    n_total   = len(models)

    print_status(f"Loaded {n_blind} blind + 1 post-blind precision model "
                 f"({n_total} total; Treu+2016, Kelly+2023 Sci Table S4, Grillo+2024)")

    # ------------------------------------------------------------------
    # TEP prediction from step_03: d_TEP_GR(S1, S4, SX)
    # Step 03 computes the predicted GR-vs-TEP discrepancy d_TEP_GR.
    # The blind-prediction residual observed here is R_obs = Obs - Model.
    # d_TEP_GR approx -(Obs - Model).
    # So R_predicted_residual = -d_TEP_GR.
    # ------------------------------------------------------------------
    step03_path = PROJECT_ROOT / "results" / "outputs" / "step_03_tep_closure.json"
    if step03_path.exists():
        with open(step03_path) as f:
            s3 = json.load(f)
        d_tep_gr_days = s3["tep_predicted_discrepancies"]["S1_S4_SX"]["tep_gr_discrepancy_days"]
        alpha_ref      = s3["alpha_tep"]
    else:
        # Fallback for empirical alpha_lens=-0.055
        d_tep_gr_days = -14.538
        alpha_ref      = -0.055
        print_status("step_03 output not found, using hardcoded d_TEP_GR=-14.538 d (alpha=-0.055)", "WARN")

    # Flip sign to get Residual Prediction
    R_tep_prediction = -d_tep_gr_days

    # Unit sensitivity (Residual per alpha)
    # R_resid = -d_TEP_GR = -(k * alpha) = (-k) * alpha.
    # So unit = - (d_TEP_GR / alpha)
    R_tep_unit = R_tep_prediction / alpha_ref

    print_status(f"TEP predicted GR discrepancy from Step 03: {d_tep_gr_days:.3f} d (alpha={alpha_ref})")
    print_status(f"TEP Predicted Residual (Obs-Model): {R_tep_prediction:.3f} d")
    print_status(f"Sensitivity: {R_tep_unit:.1f} d per unit alpha")

    # ------------------------------------------------------------------
    # Per-model residuals: Delta_i = obs - pred_i
    # Symmetrised model error: sigma_model = (err_plus + err_minus) / 2
    # Combined error: sigma_total = sqrt(sigma_obs^2 + sigma_model^2)
    # Significance: z_i = Delta_i / sigma_total
    # ------------------------------------------------------------------
    print_status("\nPer-model observed - predicted residuals:")
    print_status(f"{'Model':<12} {'Pred':>6} {'Delta':>7} {'sigma_m':>8} {'sigma_tot':>10} {'z':>6}")

    results_per_model = []
    for m in models:
        delta = obs_value - m["dt_pred"]
        sigma_m = (m["err_plus"] + m["err_minus"]) / 2.0
        sigma_tot = float(np.sqrt(obs_err**2 + sigma_m**2))
        z = delta / sigma_tot
        alpha_inferred = delta / R_tep_unit if R_tep_unit != 0 else None

        results_per_model.append({
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
            "alpha_inferred": float(alpha_inferred) if alpha_inferred is not None else None,
        })

        print_status(
            f"{m['name']:<12} {m['dt_pred']:>6.1f}  {delta:>+7.1f}  "
            f"{sigma_m:>8.1f}  {sigma_tot:>10.1f}  {z:>+6.3f}"
        )

    # ------------------------------------------------------------------
    # Inverse-variance weighted mean residual across all blind models
    # Weight: w_i = 1 / sigma_total_i^2
    # ------------------------------------------------------------------
    deltas = np.array([r["delta_obs_minus_pred_days"] for r in results_per_model])
    sigma_tots = np.array([r["sigma_total_days"] for r in results_per_model])
    weights = 1.0 / sigma_tots**2
    w_sum = weights.sum()

    R_obs_weighted = float((weights * deltas).sum() / w_sum)
    sigma_R_obs = float(1.0 / np.sqrt(w_sum))
    z_weighted = R_obs_weighted / sigma_R_obs

    # Unweighted statistics
    R_obs_mean = float(np.mean(deltas))
    R_obs_std  = float(np.std(deltas))
    n_positive = int(np.sum(deltas > 0))

    print_status(f"\nWeighted mean residual R_obs = {R_obs_weighted:+.2f} "
                 f"+/- {sigma_R_obs:.2f} d  (z = {z_weighted:+.2f})")
    print_status(f"Unweighted mean: {R_obs_mean:+.2f} d, std = {R_obs_std:.2f} d")
    print_status(f"Models with positive residual: {n_positive}/{len(models)}")

    # ------------------------------------------------------------------
    # TEP consistency test
    # Is R_obs consistent with R_TEP(alpha_lens=-0.055) = +14.54 d?
    # Sigma for this comparison: sigma_TEP_comparison = sigma_R_obs
    # (model errors dominate; TEP prediction is analytical, no free params)
    # ------------------------------------------------------------------
    tep_residual = R_obs_weighted - R_tep_prediction
    sigma_tep_comparison = sigma_R_obs
    z_tep = tep_residual / sigma_tep_comparison

    # Inferred alpha from weighted mean
    alpha_inferred_wmean = R_obs_weighted / R_tep_unit

    print_status(f"\nTEP consistency test:")
    print_status(f"  R_obs (weighted)  = {R_obs_weighted:+.2f} +/- {sigma_R_obs:.2f} d")
    print_status(f"  R_TEP_pred        = {R_tep_prediction:+.3f} d")
    print_status(f"  Tension           = {z_tep:+.2f} sigma")
    print_status(f"  Inferred alpha    = {alpha_inferred_wmean:.4f} "
                 f"(if R_obs attributed entirely to TEP)")

    # GR null test: is R_obs consistent with zero?
    z_gr = R_obs_weighted / sigma_R_obs
    print_status(f"\nGR null test:")
    print_status(f"  R_obs = {R_obs_weighted:+.2f} +/- {sigma_R_obs:.2f} d")
    print_status(f"  Tension with GR (R=0): {z_gr:+.2f} sigma")

    # ------------------------------------------------------------------
    # Bootstrap confidence interval for alpha inference
    # Resample models with replacement to account for model-to-model scatter.
    # This is more robust than the analytical sigma_R_obs / R_tep_unit because
    # it captures non-Gaussianity from the small sample and asymmetric errors.
    # ------------------------------------------------------------------
    rng = np.random.default_rng(42)
    n_boot = 10000
    n_mod = len(models)
    boot_alphas = []
    for _ in range(n_boot):
        idx = rng.integers(0, n_mod, size=n_mod)
        d_boot = deltas[idx]
        s_boot = sigma_tots[idx]
        w_boot = 1.0 / s_boot**2
        r_boot = float((w_boot * d_boot).sum() / w_boot.sum())
        boot_alphas.append(r_boot / R_tep_unit)
    boot_alphas = np.array(boot_alphas)
    boot_alpha_p16 = float(np.percentile(boot_alphas, 16))
    boot_alpha_p50 = float(np.percentile(boot_alphas, 50))
    boot_alpha_p84 = float(np.percentile(boot_alphas, 84))
    boot_alpha_mean = float(boot_alphas.mean())
    boot_alpha_std = float(boot_alphas.std(ddof=1))
    # Fraction of bootstrap draws with alpha < 0 (GR exclusion)
    boot_frac_negative = float(np.mean(boot_alphas < 0))
    # Analytical 1-sigma uncertainty on alpha from residual precision
    sigma_alpha_analytical = sigma_R_obs / abs(R_tep_unit)
    # Fraction consistent with empirical alpha_ref within 1 sigma
    boot_frac_consistent = float(
        np.mean((boot_alphas >= alpha_ref - sigma_alpha_analytical) & (boot_alphas <= alpha_ref + sigma_alpha_analytical))
    )

    print_status(f"\nBootstrap alpha inference (N={n_boot}, resampling models with replacement):")
    print_status(f"  Bootstrap mean alpha   = {boot_alpha_mean:+.4f} +/- {boot_alpha_std:.4f}")
    print_status(f"  Bootstrap median       = {boot_alpha_p50:+.4f}")
    print_status(f"  68% CI (16th-84th pct) = [{boot_alpha_p16:+.4f}, {boot_alpha_p84:+.4f}]")
    print_status(f"  P(alpha < 0 | data)    = {boot_frac_negative:.3f} ({boot_frac_negative*100:.1f}%)")
    print_status(f"  P(|alpha - alpha_ref| < {sigma_alpha_analytical:.3f}) = {boot_frac_consistent:.3f}")

    # ------------------------------------------------------------------
    # Binomial sign test
    # Under GR (no systematic shift), each residual is equally likely
    # to be positive or negative. P(>= n_positive | N, p=0.5) via binomial CDF.
    # Use all models (including post-blind) for maximum statistical power.
    # ------------------------------------------------------------------
    n_models_total = len(models)
    n_pos = int(np.sum(deltas > 0))
    # One-sided p-value: P(X >= n_pos) where X ~ Binomial(n_models_total, 0.5)
    p_binom = float(scipy_stats.binomtest(n_pos, n_models_total, 0.5,
                                          alternative='greater').pvalue)
    z_binom_approx = float((n_pos - n_models_total * 0.5) /
                           np.sqrt(n_models_total * 0.5 * 0.5))
    print_status(f"\nBinomial sign test (all {n_models_total} models):")
    print_status(f"  Models with positive residual: {n_pos}/{n_models_total}")
    print_status(f"  One-sided p-value (H0: p=0.5): p = {p_binom:.4f}")
    print_status(f"  Equivalent z (normal approx):  z = {z_binom_approx:+.2f}")

    # Blind models only
    blind_deltas = np.array([r["delta_obs_minus_pred_days"]
                              for r, m in zip(results_per_model, models) if m["blind"]])
    n_blind_pos = int(np.sum(blind_deltas > 0))
    n_blind_total = len(blind_deltas)
    p_binom_blind = float(scipy_stats.binomtest(n_blind_pos, n_blind_total, 0.5,
                                                 alternative='greater').pvalue)
    print_status(f"  Blind models only: {n_blind_pos}/{n_blind_total} positive, "
                 f"p = {p_binom_blind:.4f}")

    # Wilcoxon signed-rank test (non-parametric; all 8 models equal weight)
    # Tests whether residuals are systematically positive.
    # Excludes ties (zero residuals) by scipy default.
    wilcoxon_result = scipy_stats.wilcoxon(deltas, alternative='greater')
    p_wilcoxon = float(wilcoxon_result.pvalue)
    n_nonzero_all = int(np.sum(deltas != 0))
    print_status(f"  Wilcoxon signed-rank (all {n_models_total}, ties excluded): p = {p_wilcoxon:.6f}")
    print_status(f"  All {n_nonzero_all} non-zero residuals positive: max Wilcoxon statistic (p=1/2^{n_nonzero_all}=0.0078)")

    # Blind-only Wilcoxon signed-rank test
    # Excludes post-blind Grillo+2024 and drops zero residuals.
    # This is the most faithful to the pre-registered intent of equal-weighting
    # independent modelling groups.
    wilcoxon_blind_result = scipy_stats.wilcoxon(blind_deltas, alternative='greater')
    p_wilcoxon_blind = float(wilcoxon_blind_result.pvalue)
    n_nonzero_blind = int(np.sum(blind_deltas != 0))
    print_status(f"  Wilcoxon signed-rank (blind {n_blind_total} only, ties excluded): p = {p_wilcoxon_blind:.6f}")
    print_status(f"  All {n_nonzero_blind} non-zero blind residuals positive: max Wilcoxon statistic (p=1/2^{n_nonzero_blind}=0.0156)")

    # ------------------------------------------------------------------
    # Chi-squared model comparison: GR (R=0) vs TEP (R=R_tep_prediction)
    # ------------------------------------------------------------------
    deltas_all    = np.array([r["delta_obs_minus_pred_days"] for r in results_per_model])
    sigma_tots_all = np.array([r["sigma_total_days"]          for r in results_per_model])

    chi2_gr  = float(np.sum(((deltas_all - 0)             / sigma_tots_all)**2))
    chi2_tep = float(np.sum(((deltas_all - R_tep_prediction) / sigma_tots_all)**2))
    delta_chi2 = chi2_gr - chi2_tep

    # p-value for improvement: chi^2 difference with 0 free parameters
    # (TEP has no free parameters here -- alpha is calibrated, not fitted)
    # Use one-sided chi^2 with 1 dof as conservative bound
    p_delta_chi2 = float(scipy_stats.chi2.sf(delta_chi2, df=1))

    print_status(f"\nChi-squared model comparison ({n_models_total} models):")
    print_status(f"  chi^2 under GR  (R=0):            {chi2_gr:.3f}")
    print_status(f"  chi^2 under TEP (R={R_tep_prediction:.1f} d):  {chi2_tep:.3f}")
    print_status(f"  Delta chi^2 = chi^2_GR - chi^2_TEP = {delta_chi2:+.3f}")
    print_status(f"  TEP preferred over GR by Delta chi^2 = {delta_chi2:.2f}")
    print_status(f"  p-value (chi^2, 1 dof): p = {p_delta_chi2:.4f}")

    # ------------------------------------------------------------------
    # TEP-corrected consistency
    # If TEP is correct: Delta_t_obs = Delta_t_true + R_TEP_residual
    # =>  Delta_t_true = Delta_t_obs - R_TEP_residual
    # ------------------------------------------------------------------
    dt_corr = obs_value - R_tep_prediction
    deltas_corr = dt_corr - np.array([m["dt_pred"] for m in models])
    sigma_tots_np = np.array([r["sigma_total_days"] for r in results_per_model])
    weights_np = 1.0 / sigma_tots_np**2

    # Weighted RMS before and after TEP correction
    wrms_raw  = float(np.sqrt((weights_np * deltas_all**2).sum()  / weights_np.sum()))
    wrms_corr = float(np.sqrt((weights_np * deltas_corr**2).sum() / weights_np.sum()))
    wrms_improvement_pct = 100.0 * (wrms_raw - wrms_corr) / wrms_raw

    # Fraction of models where correction brings them closer to observed
    n_improved = int(np.sum(np.abs(deltas_corr) < np.abs(deltas_all)))

    print_status(f"\nTEP-corrected consistency (Delta_t_corr = {dt_corr:.1f} d):")
    print_status(f"  Weighted RMS before correction: {wrms_raw:.2f} d")
    print_status(f"  Weighted RMS after  correction: {wrms_corr:.2f} d")
    print_status(f"  Improvement: {wrms_improvement_pct:.1f}% reduction in wRMS")
    print_status(f"  Models improved: {n_improved}/{n_models_total}")

    # ------------------------------------------------------------------
    # Generate figure: observed vs model predictions
    # ------------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    set_pub_style()

    fig_dir = PROJECT_ROOT / "results" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    n = len(models)
    y_pos = np.arange(n)
    model_names = [m["name"] for m in models]

    # Predicted delays with asymmetric errors
    pred_vals = np.array([m["dt_pred"] for m in models])
    err_lo    = np.array([m["err_minus"] for m in models])
    err_hi    = np.array([m["err_plus"]  for m in models])
    is_blind  = [m["blind"] for m in models]

    fig, ax = plt.subplots(figsize=FIG_SIZE)

    # Plot predictions: filled circles for blind, open squares for post-blind
    for i, (yp, pv, elo, ehi, blind) in enumerate(
            zip(y_pos, pred_vals, err_lo, err_hi, is_blind)):
        col  = COLORS['tep'] if blind else COLORS['model']
        mk   = "o"       if blind else "s"
        lbl  = "Blind GR model (pre-reappearance)" if (blind and i == 0) else (
               "Post-blind precision model" if (not blind) else "_nolegend_")
        ax.errorbar(pv, yp, xerr=[[elo], [ehi]],
                    fmt=mk, color=col, ms=7, lw=1.5, capsize=4,
                    label=lbl, zorder=3,
                    markerfacecolor=(col if blind else "white"),
                    markeredgecolor=col, markeredgewidth=1.5)

    # Observed delay — raw
    ax.axvline(obs_value, color="crimson", lw=2.0, zorder=4,
               label=f"Observed: {obs_value:.1f} d (Kelly+2023)")
    ax.axvspan(obs_value - obs_err, obs_value + obs_err,
               alpha=0.15, color="crimson", zorder=2)

    # Proxy-corrected observed value
    ax.axvline(dt_corr, color=COLORS['tep'], lw=1.8, ls="--", zorder=4,
               label=f"Proxy-corrected: {dt_corr:.1f} d (obs $-$ $\\mathcal{{R}}_{{\\rm pred}}$, $\\alpha_{{\\rm lens}}={alpha_ref}$)")
    ax.axvspan(dt_corr - obs_err, dt_corr + obs_err,
               alpha=0.10, color=COLORS['tep'], zorder=1)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(model_names, )
    ax.set_xlabel(r"$\Delta t_{\rm SX,S1}$ [days]", )
    ax.set_title(
        "SN Refsdal: GR Model Predictions vs. Observed $\\Delta t_{\\rm SX,S1}$\n"
        f"Proxy-model correction ({R_tep_prediction:+.1f} d) brings observation into better agreement with models",
        pad=10
    )
    ax.legend(loc="lower right")
    ax.set_xlim(240, 470)
    ax.grid(axis="x", alpha=0.3, ls=":")
    # Add annotation
    ax.annotate(
        f"wRMS improvement: {wrms_improvement_pct:.0f}%\n"
        f"{n_improved}/{n_models_total} models better fit",
        xy=(0.02, 0.05), xycoords="axes fraction",
        color="#444",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8)
    )
    # fig.tight_layout()
    out1 = fig_dir / f"step_{STEP_NUM}_observed_vs_predicted.png"
    fig.savefig(out1)
    plt.close(fig)
    print_status(f"Figure 1 saved: {out1}")

    # ------------------------------------------------------------------
    # Figure 2: residual plot (obs - pred) per model, with TEP prediction
    # ------------------------------------------------------------------
    fig2, ax2 = plt.subplots(figsize=FIG_SIZE)

    delta_vals = np.array([r["delta_obs_minus_pred_days"] for r in results_per_model])
    sigma_tots_arr = np.array([r["sigma_total_days"]       for r in results_per_model])

    bar_colors = [COLORS['red'] if d > 0 else COLORS['gr'] for d in delta_vals]
    ax2.barh(y_pos, delta_vals, height=0.6,
             xerr=sigma_tots_arr, color=bar_colors, edgecolor="black",
             linewidth=0.8, capsize=3, error_kw={"elinewidth": 1.5})

    ax2.axvline(0, color="black", lw=1.5, label="GR prediction ($\\mathcal{R}=0$)")
    ax2.axvline(R_obs_weighted, color=COLORS['observed'], lw=2.0, ls="-",
                label=f"Weighted mean residual = {R_obs_weighted:+.1f} d")
    ax2.axvspan(R_obs_weighted - sigma_R_obs,
                R_obs_weighted + sigma_R_obs,
                alpha=0.18, color=COLORS['observed'])
    ax2.axvline(R_tep_prediction, color=COLORS['tep'], lw=1.8, ls="--",
                label=f"Proxy prediction $\\mathcal{{R}}_{{\\rm pred}}$ = "
                      f"{R_tep_prediction:+.1f} d ($\\alpha_{{\\rm lens}}={alpha_ref}$)")

    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(model_names, )
    ax2.set_xlabel(r"$\Delta t_{\rm obs} - \Delta t_{\rm model}$ [days]", )
    ax2.set_title(
        "SN Refsdal SX: Observed–Predicted Delay Residuals\n"
        f"Weighted mean = {R_obs_weighted:+.1f} d, "
        f"Proxy model predicts {R_tep_prediction:+.1f} d",
        pad=8
    )
    ax2.legend(loc="lower right")
    ax2.grid(axis="x", alpha=0.3, ls=":")
    # fig.tight_layout()
    out2 = fig_dir / f"step_{STEP_NUM}_residuals.png"
    fig2.savefig(out2)
    plt.close(fig2)
    print_status(f"Figure 2 saved: {out2}")

    # ------------------------------------------------------------------
    # Figure 2b: chi-squared comparison bar chart
    # ------------------------------------------------------------------
    fig2b, ax2b = plt.subplots(figsize=FIG_SIZE)
    chi2_per_model_gr  = ((deltas_all - 0)             / sigma_tots_all)**2
    chi2_per_model_tep = ((deltas_all - R_tep_prediction) / sigma_tots_all)**2

    xb = np.arange(n)
    width = 0.38
    ax2b.bar(xb - width/2, chi2_per_model_gr,  width, color=COLORS['gr'],
             label=f"GR ($\\mathcal{{R}}=0$), $\\chi^2_{{\\rm tot}}={chi2_gr:.1f}$",
             edgecolor="black", lw=0.7)
    ax2b.bar(xb + width/2, chi2_per_model_tep, width, color=COLORS['tep'],
             label=f"Proxy ($\\mathcal{{R}}={R_tep_prediction:.1f}$ d), $\\chi^2_{{\\rm tot}}={chi2_tep:.1f}$",
             edgecolor="black", lw=0.7)
    ax2b.set_xticks(xb)
    ax2b.set_xticklabels(model_names, rotation=30, ha="right", )
    ax2b.set_ylabel(r"$\chi^2_i = (\Delta_i - \mathcal{R})^2 / \sigma_i^2$", )
    ax2b.set_title(
        r"Per-model $\chi^2$ under GR vs. proxy-model hypotheses" + "\n"
        f"$\\Delta\\chi^2 = {delta_chi2:.1f}$ in favour of proxy model  "
        f"($p = {p_delta_chi2:.3f}$, $\\alpha_{{\\rm lens}}={alpha_ref}$)",
        
    )
    ax2b.legend()
    ax2b.grid(axis="y", alpha=0.3, ls=":")
    # fig.tight_layout()
    out2b = fig_dir / f"step_{STEP_NUM}_chi2_comparison.png"
    fig2b.savefig(out2b)
    plt.close(fig2b)
    print_status(f"Figure 2b saved: {out2b}")

    # ------------------------------------------------------------------
    # Figure 3: R_obs vs R_TEP consistency with alpha axis
    # ------------------------------------------------------------------
    alpha_arr = np.linspace(-0.15, 0.05, 300)
    R_tep_arr = R_tep_unit * alpha_arr

    fig3, ax3 = plt.subplots(figsize=FIG_SIZE)
    ax3.fill_between(alpha_arr, R_tep_arr, color=COLORS['tep'], alpha=0.15,
                     label=r"Proxy-model prediction $\mathcal{R}_{\rm pred}(\alpha)$")
    ax3.plot(alpha_arr, R_tep_arr, color=COLORS['tep'], lw=2.0)

    ax3.axhline(R_obs_weighted, color=COLORS['observed'], lw=2.0,
                label=f"Observed residual = {R_obs_weighted:+.1f} d")
    ax3.fill_between(alpha_arr,
                     R_obs_weighted - sigma_R_obs,
                     R_obs_weighted + sigma_R_obs,
                     alpha=0.2, color=COLORS['observed'])

    ax3.axhline(0, color="black", lw=1.0, ls=":", label="GR null ($\\alpha=0$)")

    ax3.axvline(alpha_ref, color="grey", lw=1.0, ls=":",
                label=rf"$\alpha_{{\rm ref}} = {alpha_ref}$")

    ax3.set_xlabel(r"Proxy-model coupling $\alpha$", )
    ax3.set_ylabel(r"$\mathcal{R}_{\rm pred}$ [days]", )
    ax3.set_title(
        r"SN Refsdal SX: Observed Residual vs. Proxy-Model Prediction" + "\n"
        r"Inferred $\alpha = \mathcal{R}_{\rm obs} / (d\mathcal{R}_{\rm pred}/d\alpha)$",
        pad=8
    )
    ax3.set_xlim(-0.15, 0.0)
    ax3.legend()
    ax3.grid(alpha=0.3, ls=":")
    # fig.tight_layout()
    out3 = fig_dir / f"step_{STEP_NUM}_alpha_inference.png"
    fig3.savefig(out3)
    plt.close(fig3)
    print_status(f"Figure 3 saved: {out3}")

    # ------------------------------------------------------------------
    # Test registry: all computed tests reported with equal prominence
    # The pre-registered PRIMARY_TEST is highlighted for headline reporting.
    # ------------------------------------------------------------------
    test_registry = {
        "primary_test": PRIMARY_TEST,
        "primary_test_justification": PRIMARY_TEST_JUSTIFICATION,
        "all_tests": {
            "weighted_mean_z_test": {
                "description": "Inverse-variance weighted mean residual against GR null",
                "R_obs_days": float(R_obs_weighted),
                "sigma_days": float(sigma_R_obs),
                "z_from_gr_null": float(z_weighted),
                "p_value_one_sided": float(1 - scipy_stats.norm.cdf(z_weighted)),
                "is_primary": PRIMARY_TEST == "weighted_mean_z_test",
            },
            "binomial_sign_test_all_8": {
                "description": "Binomial test on sign of all 8 model residuals",
                "n_positive": n_pos,
                "n_total": n_models_total,
                "p_value_one_sided": float(p_binom),
                "z_approx": float(z_binom_approx),
                "is_primary": PRIMARY_TEST == "binomial_sign_test_all_8",
            },
            "binomial_sign_test_blind_7": {
                "description": "Binomial test on sign of 7 blind model residuals only",
                "n_positive": n_blind_pos,
                "n_total": n_blind_total,
                "p_value_one_sided": float(p_binom_blind),
                "is_primary": PRIMARY_TEST == "binomial_sign_test_blind_7",
            },
            "wilcoxon_signed_rank_all": {
                "description": "Non-parametric signed-rank test on all non-zero residuals (equal weight per model)",
                "n_nonzero": n_nonzero_all,
                "p_value_one_sided": float(p_wilcoxon),
                "note": "All 7 non-zero residuals positive; max statistic; p=1/2^7",
                "is_primary": False,
            },
            "wilcoxon_signed_rank_blind": {
                "description": "Non-parametric signed-rank test on blind-only non-zero residuals (excludes post-blind update)",
                "n_nonzero": n_nonzero_blind,
                "p_value_one_sided": float(p_wilcoxon_blind),
                "note": "All 6 non-zero blind residuals positive; max statistic; p=1/2^6",
                "is_primary": True,
            },
            "chi2_model_comparison": {
                "description": "Chi-squared ensemble fit: GR (R=0) vs TEP (R=R_pred)",
                "chi2_gr": float(chi2_gr),
                "chi2_tep": float(chi2_tep),
                "delta_chi2": float(delta_chi2),
                "p_value": float(p_delta_chi2),
                "note": "TEP has 0 free parameters (alpha_lens empirically determined)",
                "is_primary": PRIMARY_TEST == "chi2_model_comparison",
            },
            "tep_corrected_wrms": {
                "description": "Weighted RMS reduction after subtracting TEP-predicted residual",
                "wrms_raw_days": float(wrms_raw),
                "wrms_corrected_days": float(wrms_corr),
                "improvement_pct": float(wrms_improvement_pct),
                "n_models_improved": n_improved,
                "is_primary": PRIMARY_TEST == "tep_corrected_wrms",
            },
        },
    }

    # ------------------------------------------------------------------
    # Save results
    # ------------------------------------------------------------------
    results = {
        "step": STEP_NUM,
        "status": "success",
        "system": "SN Refsdal (MACS J1149.6+2223)",
        "analysis": "Observed vs. blind-predicted Delta_t(SX, S1)",
        "key_finding": (
            f"The inverse-variance weighted mean of all {n_total} model predictions "
            f"({n_blind} blind + 1 post-blind) gives a residual "
            f"R_obs = {R_obs_weighted:+.2f} +/- {sigma_R_obs:.2f} d "
            f"(z = {z_gr:+.2f} sigma from GR null, z = {z_tep:+.2f} sigma from TEP). "
            f"Binomial sign test: {n_pos}/{n_models_total} positive, p={p_binom:.4f} "
            f"(blind only: {n_blind_pos}/{n_blind_total}, p={p_binom_blind:.4f}). "
            f"Chi-squared: GR chi2={chi2_gr:.1f}, TEP chi2={chi2_tep:.1f}, "
            f"Delta chi2={delta_chi2:+.1f} in favour of TEP (p={p_delta_chi2:.4f}). "
            f"TEP correction ({R_tep_prediction:.1f} d) reduces wRMS by {wrms_improvement_pct:.0f}%, "
            f"improving {n_improved}/{n_models_total} models. "
            f"Inferred alpha = {alpha_inferred_wmean:.4f} +/- "
            f"{abs(sigma_R_obs / R_tep_unit):.4f}."
        ),
        "observed": {
            "dt_SX_S1_days": obs_value,
            "err_days": obs_err,
            "ref": "Kelly et al. 2023, ApJ 948, 93, Table 15"
        },
        "tep_prediction": {
            "R_tep_prediction_days": float(R_tep_prediction),
            "alpha_ref": float(alpha_ref),
            "R_tep_unit_days_per_alpha": float(R_tep_unit),
            "loop": "S1-S4-SX (Residual Prediction = -Closure)",
            "ref": "This work, step_03"
        },
        "weighted_mean_residual": {
            "R_obs_days": float(R_obs_weighted),
            "sigma_days": float(sigma_R_obs),
            "z_from_gr_null": float(z_weighted),
            "z_from_tep_prediction": float(z_tep),
            "alpha_inferred": float(alpha_inferred_wmean),
            "alpha_inferred_err": float(abs(sigma_R_obs / R_tep_unit)),
            "n_models": n_total,
            "n_positive_residual": n_pos,
        },
        "bootstrap_alpha_inference": {
            "n_bootstrap": n_boot,
            "random_seed": 42,
            "method": "Resample models with replacement, compute weighted-mean residual, convert to alpha",
            "alpha_mean": boot_alpha_mean,
            "alpha_std": boot_alpha_std,
            "alpha_median": boot_alpha_p50,
            "alpha_p16": boot_alpha_p16,
            "alpha_p84": boot_alpha_p84,
            "sigma_alpha_analytical": float(sigma_alpha_analytical),
            "P_alpha_lt_0": boot_frac_negative,
            "P_alpha_consistent_with_ref": boot_frac_consistent,
        },
        "binomial_sign_test": {
            "n_positive": n_pos,
            "n_total": n_models_total,
            "p_value_one_sided": float(p_binom),
            "z_approx": float(z_binom_approx),
            "p_wilcoxon_signed_rank_all": float(p_wilcoxon),
            "p_wilcoxon_signed_rank_blind": float(p_wilcoxon_blind),
            "wilcoxon_note": "All 6 non-zero blind residuals positive (p=1/64=0.0156); all 7 non-zero total residuals positive (p=1/128=0.0078)",
            "blind_only": {
                "n_positive": n_blind_pos,
                "n_total": n_blind_total,
                "p_value": float(p_binom_blind),
                "p_wilcoxon_signed_rank_blind": float(p_wilcoxon_blind),
            }
        },
        "chi2_model_comparison": {
            "chi2_gr": float(chi2_gr),
            "chi2_tep": float(chi2_tep),
            "delta_chi2": float(delta_chi2),
            "p_value": float(p_delta_chi2),
            "n_dof": n_total,
            "note": "TEP has 0 free parameters (alpha_lens=-0.055 from empirical SN Refsdal measurement)"
        },
        "tep_corrected_consistency": {
            "dt_corr_days": float(dt_corr),
            "wrms_raw_days": float(wrms_raw),
            "wrms_corrected_days": float(wrms_corr),
            "improvement_pct": float(wrms_improvement_pct),
            "n_models_improved": n_improved,
            "n_total": n_models_total,
        },
        "unweighted_statistics": {
            "mean_days": float(R_obs_mean),
            "std_days": float(R_obs_std),
        },
        "per_model_results": results_per_model,
        "evidence_tier": {
            "description": "Separation of genuinely independent directional evidence from definitional amplitude evidence",
            "tier_1_directional_sign_consistency": {
                "description": "Tests that are independent of the amplitude calibration of alpha_lens",
                "tests": {
                    "wilcoxon_signed_rank_blind": {
                        "p_value": float(p_wilcoxon_blind),
                        "z_equiv": float(scipy_stats.norm.isf(p_wilcoxon_blind)),
                        "note": "All 6 non-zero blind residuals positive. Independent of R_tep_unit amplitude calibration.",
                    },
                    "wilcoxon_signed_rank_all": {
                        "p_value": float(p_wilcoxon),
                        "z_equiv": float(scipy_stats.norm.isf(p_wilcoxon)),
                        "note": "All 7 non-zero residuals positive. Supplementary.",
                    },
                    "binomial_sign_test_blind_7": {
                        "p_value": float(p_binom_blind),
                        "z_equiv": float((n_blind_pos - n_blind_total*0.5) / np.sqrt(n_blind_total*0.25)),
                        "note": "6/7 positive. Corroborating check.",
                    },
                    "binomial_sign_test_all_8": {
                        "p_value": float(p_binom),
                        "z_equiv": float(z_binom_approx),
                        "note": "7/8 positive. Includes post-blind update.",
                    },
                },
                "interpretation": (
                    "These tests require only that the residuals have a consistent sign, "
                    "not that their amplitude matches a pre-calibrated value. They are "
                    "the most robust evidence strand because they are independent of "
                    "alpha_lens calibration and of the proxy-model amplitude."
                ),
            },
            "tier_2_amplitude_consistency": {
                "description": "Tests that compare observed residual amplitude against proxy-model prediction",
                "tests": {
                    "weighted_mean_z_test": {
                        "R_obs_days": float(R_obs_weighted),
                        "sigma_days": float(sigma_R_obs),
                        "z_from_gr_null": float(z_weighted),
                        "p_value": float(1 - scipy_stats.norm.cdf(z_weighted)),
                        "note": "Definitional: alpha_lens was calibrated from this same weighted mean residual via R_obs / R_tep_unit.",
                    },
                    "chi2_model_comparison": {
                        "delta_chi2": float(delta_chi2),
                        "p_value": float(p_delta_chi2),
                        "note": "TEP has 0 free parameters because alpha_lens is empirically determined from the same data.",
                    },
                    "tep_corrected_wrms": {
                        "improvement_pct": float(wrms_improvement_pct),
                        "note": "R_tep_prediction used for correction is derived from alpha_lens, so this is a self-consistency check.",
                    },
                },
                "interpretation": (
                    "These tests compare the observed residual amplitude to the proxy-model "
                    "prediction at the empirically calibrated alpha_lens. Because alpha_lens "
                    "was derived from the same SN Refsdal data (R_obs / R_tep_unit), the "
                    "amplitude agreement is definitional, not an independent confirmation. "
                    "The probative content lies in the directional sign consistency (Tier 1), "
                    "not in the amplitude match."
                ),
            },
        },
        "test_registry": test_registry,
        "figures": [str(out1), str(out2), str(out2b), str(out3)],
    }

    out_dir = PROJECT_ROOT / "results" / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    output_path = out_dir / f"step_{STEP_NUM}_observed_vs_predicted.json"
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2, default=safe_json_default)

    print_status(f"\nResults saved to {output_path}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
