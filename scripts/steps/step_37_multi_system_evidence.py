#!/usr/bin/env python3
"""
TEP-LENS: Step 37 - Multi-System Evidence Accumulation Projection

Purpose: Quantify how the TEP evidence strengthens as independent multiply-imaged
supernova systems are added. The current SN Refsdal measurement is a single-system,
single-contrast test (D_eff ~ 2.0). No single system can yield decisive evidence
because the proxy model is phenomenological and alpha_proxy was calibrated on the
data it tests. This step projects the number of independent systems required for
decisive significance under valid statistical combination rules.

Valid combination methods (independent systems only):
1. Stouffer's z-method for independent directional tests.
2. Inverse-variance weighted alpha inference.
3. Binomial across systems: P(k of N positive | GR).

Invalid methods (not used):
- Fisher's method: requires independent p-values from uncorrelated tests.
  Within each system, the tests are correlated; across systems they are not,
  but Fisher assumes p-values test the same null hypothesis in the same way.
- Naive p-value multiplication: ignores test correlation and differing power.

Assumptions:
- Each new system is statistically independent of SN Refsdal.
- Per-system precision is characterised by blind-model uncertainty sigma_model.
- The proxy-model coupling alpha_proxy is common across all lens systems.
- Long-baseline contrast dominates the signal (proven by Step 35 for Refsdal).
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status

STEP_NUM = "37"


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def stouffer_combine(z_scores):
    """Stouffer's method for combining independent z-scores.
    z_total = sum(z_i) / sqrt(N). Valid when each z_i tests the same
    directional hypothesis on independent data."""
    z_scores = np.array(z_scores)
    n = len(z_scores)
    if n == 0:
        return {"z_total": 0.0, "p_total": 1.0, "n_systems": 0}
    z_total = float(z_scores.sum() / np.sqrt(n))
    p_total = float(stats.norm.sf(z_total))
    return {"z_total": z_total, "p_total": p_total, "n_systems": n}


def inverse_variance_alpha(alphas, sigma_alphas):
    """Combine independent alpha inferences from multiple systems."""
    alphas = np.array(alphas)
    sigma_alphas = np.array(sigma_alphas)
    weights = 1.0 / sigma_alphas**2
    w_sum = weights.sum()
    alpha_combined = float((weights * alphas).sum() / w_sum)
    sigma_combined = float(1.0 / np.sqrt(w_sum))
    z_vs_zero = alpha_combined / sigma_combined
    # One-sided p-value for alpha < 0 (predicted direction).
    # norm.cdf(z_obs) = P(Z < z_obs) gives the correct left-tail probability.
    p_vs_zero = float(stats.norm.cdf(z_vs_zero))
    return {
        "alpha_combined": alpha_combined,
        "sigma_combined": sigma_combined,
        "z_vs_zero": z_vs_zero,
        "p_vs_zero": p_vs_zero,
        "n_systems": len(alphas),
    }


def binomial_cross_system(n_positive, n_total):
    """Binomial test: P(>= n_positive | p=0.5, N=n_total).
    Under GR, each independent system has 50% chance of positive residual."""
    p = float(stats.binomtest(n_positive, n_total, 0.5, alternative="greater").pvalue)
    z = float((n_positive - n_total * 0.5) / np.sqrt(n_total * 0.25))
    return {"n_positive": n_positive, "n_total": n_total, "p_value": p, "z_approx": z}


def load_headline_significance(step16):
    """Return the current correlation-aware headline test with legacy fallback."""
    if "headline_significance" in step16:
        headline = step16["headline_significance"]
        return {
            "test": headline.get("test", "headline significance"),
            "z_score": float(headline["z_score"]),
            "p_value": float(headline["p_value"]),
        }

    hierarchy = step16.get("test_hierarchy", {})
    if "tier_1b_primary_correlation_aware" in hierarchy:
        headline = hierarchy["tier_1b_primary_correlation_aware"]
        return {
            "test": headline.get("test", "Exact family-sign-flip test blind"),
            "z_score": float(headline["z_score"]),
            "p_value": float(headline["p_value"]),
        }

    tests = step16.get("individual_tests", {})
    if "exact_permutation" in tests:
        headline = tests["exact_permutation"]
        return {
            "test": "Exact family-sign-flip test blind",
            "z_score": float(headline["z"]),
            "p_value": float(headline["p"]),
        }

    raise KeyError(
        "Step 16 output does not contain headline_significance, "
        "test_hierarchy.tier_1b_primary_correlation_aware, or individual_tests.exact_permutation"
    )


def per_system_z(sigma_model_days, R_tep_true=14.538, n_blind_models=7):
    """Expected z-score per system given model uncertainty.
    z = R_tep / (sigma_model / sqrt(n_models))
    """
    sigma_mean = sigma_model_days / np.sqrt(n_blind_models)
    z = R_tep_true / sigma_mean
    return z


def main():
    print_status(
        f"STEP {STEP_NUM}: Multi-System Evidence Accumulation Projection",
        "TITLE",
    )

    rng = np.random.default_rng(42)

    # ------------------------------------------------------------------
    # Load current Refsdal template
    # ------------------------------------------------------------------
    s07_path = PROJECT_ROOT / "results" / "outputs" / "step_07_observed_vs_predicted.json"
    s16_path = PROJECT_ROOT / "results" / "outputs" / "step_16_tier_significance.json"

    with open(s07_path) as f:
        s07 = json.load(f)
    with open(s16_path) as f:
        s16 = json.load(f)

    # Template alpha from Refsdal
    alpha_ref = float(s07["weighted_mean_residual"]["alpha_inferred"])
    sigma_alpha_ref = float(s07["weighted_mean_residual"]["alpha_inferred_err"])
    r_tep_true = float(s07["tep_prediction"]["R_tep_prediction_days"])

    # Template headline z from the current correlation-aware primary test.
    headline = load_headline_significance(s16)
    z_template = float(headline["z_score"])
    p_template = float(headline["p_value"])
    headline_test = headline["test"]

    print_status(f"Template system (SN Refsdal SX):")
    print_status(f"  alpha_inferred = {alpha_ref:+.4f} +/- {sigma_alpha_ref:.4f}")
    print_status(f"  R_tep = {r_tep_true:.3f} d")
    print_status(f"  {headline_test}: z = {z_template:.2f}, p = {p_template:.4f}")

    # ------------------------------------------------------------------
    # 1. Stouffer projection: identical systems
    # ------------------------------------------------------------------
    max_n = 20
    stouffer_results = []
    for n in range(1, max_n + 1):
        # Assume each new system matches Refsdal precision
        z_per = z_template
        z_total = z_per * np.sqrt(n)
        p_total = float(stats.norm.sf(z_total))
        stouffer_results.append({
            "n_systems": n,
            "z_per_system": float(z_per),
            "z_combined": float(z_total),
            "p_combined": p_total,
            "reaches_3sigma": z_total >= 3.0,
            "reaches_5sigma": z_total >= 5.0,
        })

    # Find thresholds
    n_for_3sigma = next((r["n_systems"] for r in stouffer_results if r["reaches_3sigma"]), None)
    n_for_5sigma = next((r["n_systems"] for r in stouffer_results if r["reaches_5sigma"]), None)

    print_status(f"\nStouffer projection (identical Refsdal-precision systems):")
    print_status(f"  3-sigma evidence:  N = {n_for_3sigma} systems")
    print_status(f"  5-sigma discovery: N = {n_for_5sigma} systems")

    # ------------------------------------------------------------------
    # 2. Inverse-variance alpha projection
    # ------------------------------------------------------------------
    alpha_proj = []
    for n in range(1, max_n + 1):
        alphas = [alpha_ref] * n
        sigmas = [sigma_alpha_ref] * n
        comb = inverse_variance_alpha(alphas, sigmas)
        alpha_proj.append({
            "n_systems": n,
            "alpha_combined": comb["alpha_combined"],
            "sigma_combined": comb["sigma_combined"],
            "z_vs_zero": comb["z_vs_zero"],
            "p_vs_zero": comb["p_vs_zero"],
            "reaches_3sigma": abs(comb["z_vs_zero"]) >= 3.0,
            "reaches_5sigma": abs(comb["z_vs_zero"]) >= 5.0,
        })

    n_alpha_3sig = next((r["n_systems"] for r in alpha_proj if r["reaches_3sigma"]), None)
    n_alpha_5sig = next((r["n_systems"] for r in alpha_proj if r["reaches_5sigma"]), None)

    print_status(f"\nAlpha-inference projection (identical systems):")
    print_status(f"  Combined sigma_alpha = {sigma_alpha_ref:.4f} / sqrt(N)")
    print_status(f"  3-sigma exclusion of alpha=0: N = {n_alpha_3sig}")
    print_status(f"  5-sigma exclusion of alpha=0: N = {n_alpha_5sig}")

    # ------------------------------------------------------------------
    # 3. Precision-vs-N trade-off
    # ------------------------------------------------------------------
    sigma_models = np.linspace(5.0, 60.0, 30)
    n_systems_grid = np.arange(1, max_n + 1)

    precision_grid = np.zeros((len(sigma_models), len(n_systems_grid)))
    for i, sigma_m in enumerate(sigma_models):
        for j, n in enumerate(n_systems_grid):
            z_per = per_system_z(sigma_m, r_tep_true)
            z_total = z_per * np.sqrt(n)
            precision_grid[i, j] = z_total

    # Find contour lines
    def find_contour(grid, target, xs, ys):
        points = []
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1]):
                if grid[i, j] >= target:
                    points.append((float(xs[j]), float(ys[i])))
        return points

    contour_3sig = find_contour(precision_grid, 3.0, n_systems_grid, sigma_models)
    contour_5sig = find_contour(precision_grid, 5.0, n_systems_grid, sigma_models)

    # ------------------------------------------------------------------
    # 4. Binomial cross-system: how many systems need positive residuals?
    # ------------------------------------------------------------------
    binomial_proj = []
    for n in range(1, max_n + 1):
        # Under TEP: P(positive) = 1.0 (from Refsdal's 6/6 blind positive)
        # Under GR: P(positive) = 0.50
        # For k positive out of n, binomial p = P(X >= k | p=0.5, n)
        for k in range(n, 0, -1):
            b = binomial_cross_system(k, n)
            if b["p_value"] <= 0.05:
                binomial_proj.append({
                    "n_systems": n,
                    "min_positive_for_05": k,
                    "p_value_at_k": b["p_value"],
                    "z_approx": b["z_approx"],
                })
                break
        else:
            binomial_proj.append({
                "n_systems": n,
                "min_positive_for_05": n + 1,  # impossible
                "p_value_at_k": 1.0,
                "z_approx": 0.0,
            })

    # ------------------------------------------------------------------
    # 5. Designated multi-system thresholds
    # ------------------------------------------------------------------
    prospective_thresholds = {
        "stouffer_3sigma_n_systems": n_for_3sigma,
        "stouffer_5sigma_n_systems": n_for_5sigma,
        "alpha_3sigma_n_systems": n_alpha_3sig,
        "alpha_5sigma_n_systems": n_alpha_5sig,
        "note": (
            "These thresholds assume each new system provides an independent blind-prediction "
            "residual test with precision equivalent to SN Refsdal (blind model uncertainty "
            "sigma_model ~ 16-60 d, 6 blind models, R_tep ~ 14.5 d). Degraded precision in "
            "new systems increases the required N proportionally. These thresholds are "
            "designated prospectively for future systems and should not be adjusted post-measurement."
        ),
    }

    print_status(f"\nProspective thresholds:")
    print_status(f"  {prospective_thresholds['note']}")

    # ------------------------------------------------------------------
    # Save outputs
    # ------------------------------------------------------------------
    out = {
        "step": STEP_NUM,
        "status": "success",
        "description": "Multi-system evidence accumulation projection",
        "template_system": {
            "system": "SN Refsdal SX",
            "alpha_inferred": alpha_ref,
            "sigma_alpha": sigma_alpha_ref,
            "R_tep_days": r_tep_true,
            "headline_test": headline_test,
            "z_headline": z_template,
            "p_headline": p_template,
        },
        "stouffer_projection": {
            "method": "Stouffer's z-method for independent systems",
            "assumption": f"Each system provides identical-precision {headline_test}",
            "n_for_3sigma": n_for_3sigma,
            "n_for_5sigma": n_for_5sigma,
            "per_system_results": stouffer_results,
        },
        "alpha_inference_projection": {
            "method": "Inverse-variance weighted alpha across independent systems",
            "assumption": "Each system provides independent alpha inference with same precision as Refsdal",
            "n_for_3sigma": n_alpha_3sig,
            "n_for_5sigma": n_alpha_5sig,
            "per_system_results": alpha_proj,
        },
        "binomial_cross_system": {
            "method": "Binomial sign test across independent systems",
            "assumption": "Each system yields a positive or negative residual under GR null",
            "results": binomial_proj,
        },
        "precision_vs_n_grid": {
            "sigma_model_days": [float(s) for s in sigma_models],
            "n_systems": [int(n) for n in n_systems_grid],
            "z_combined_grid": precision_grid.tolist(),
            "contour_3sigma_points": contour_3sig,
            "contour_5sigma_points": contour_5sig,
        },
        "prospective_thresholds": prospective_thresholds,
    }

    out_path = PROJECT_ROOT / "results" / "outputs" / f"step_{STEP_NUM}_multi_system_evidence.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=safe_json_default)
    print_status(f"\nResults saved to {out_path}")

    # ------------------------------------------------------------------
    # Figure: precision vs N contour plot
    # ------------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 6))

        # Contour plot
        X, Y = np.meshgrid(n_systems_grid, sigma_models)
        cs = ax.contour(X, Y, precision_grid, levels=[3, 5], colors=["C0", "C1"], linewidths=2)
        ax.clabel(cs, inline=True, fontsize=10, fmt="z = %1.0f$\\sigma$")

        ax.set_xlabel("Number of independent systems ($N$)")
        ax.set_ylabel("Per-system model uncertainty ($\\sigma_{\\rm model}$, days)")
        ax.set_title("Multi-System Evidence Accumulation (Stouffer Combination)")
        ax.set_xlim(1, max_n)
        ax.set_ylim(60, 5)
        ax.grid(True, alpha=0.3)

        # Add Refsdal point
        sigma_refsdal = 20.0  # approximate median
        ax.scatter([1], [sigma_refsdal], color="red", s=100, zorder=5, label="SN Refsdal (template)")
        ax.legend(loc="upper right")

        fig.tight_layout()
        fig_path = (
            PROJECT_ROOT / "results" / "figures" / f"step_{STEP_NUM}_multi_system_evidence.png"
        )
        fig.savefig(fig_path, dpi=300, bbox_inches="tight")
        print_status(f"Figure saved to {fig_path}")
        plt.close(fig)
    except Exception as e:
        print_status(f"Figure generation skipped: {e}", "WARNING")

    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
