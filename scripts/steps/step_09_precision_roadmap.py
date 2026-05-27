#!/usr/bin/env python3
"""
TEP-LENS: Step 09 - Precision Roadmap

Simulates the statistical significance of the TEP vs GR test as lens modeling
uncertainties shrink. Current uncertainties are large (16-60 days). Near-future
cluster models (e.g., from JWST data + deep MUSE spectroscopy) aim for < 5 days.

This script computes:
1. Expected z-score (from inverse-variance weighted mean) as a function of sigma_model.
2. Expected fraction of positive residuals (binomial test) as a function of sigma_model.
3. Decisive falsification threshold.
"""

import json
import sys
from pathlib import Path
import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status

STEP_NUM = "09"

def safe_json_default(obj):
    if hasattr(obj, 'item'):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

def main():
    print_status(f"STEP {STEP_NUM}: Precision Roadmap Simulation", "TITLE")

    # Current true signal (from SN Refsdal SX, alpha_lens=-0.055)
    R_tep_true = 14.538  # days
    
    # Range of future model uncertainties
    sigma_models = np.logspace(0.3, 1.8, 50) # From ~2 days to ~63 days
    
    # Number of independent models
    n_models = 8
    
    results_list = []
    
    for sigma in sigma_models:
        # 1. z-score of weighted mean
        # Weighted mean error for N models each with error sigma
        # sigma_mean = sigma / sqrt(N)
        sigma_mean = sigma / np.sqrt(n_models)
        expected_z = R_tep_true / sigma_mean
        
        # 2. Binomial test expected p-value
        # Probability that a single model residual is > 0
        # If true mean is 14.538 and scatter is sigma
        p_positive = stats.norm.sf(0, loc=R_tep_true, scale=sigma)
        
        # Expected number of positive residuals
        expected_positive = n_models * p_positive
        
        # Probability of getting >= 7/8 positive
        prob_7_or_8 = stats.binom.sf(6, n_models, p_positive)
        # Probability of getting exactly 8/8
        prob_8 = stats.binom.pmf(8, n_models, p_positive)
        
        results_list.append({
            "sigma_model_days": float(sigma),
            "expected_z_score": float(expected_z),
            "p_positive_per_model": float(p_positive),
            "prob_7_or_8_positive": float(prob_7_or_8),
            "prob_8_positive": float(prob_8)
        })

    # Find the threshold where expected z-score crosses 3 and 5
    sigma_3sigma = None
    sigma_5sigma = None
    for r in results_list:
        if r["expected_z_score"] >= 3.0 and sigma_3sigma is None:
            sigma_3sigma = r["sigma_model_days"]
        if r["expected_z_score"] >= 5.0 and sigma_5sigma is None:
            sigma_5sigma = r["sigma_model_days"]
            
    # For exactly 3 sigma:
    # 3 = 14.538 / (sigma / sqrt(8))
    # sigma = 14.538 * sqrt(8) / 3 = 13.70 days
    exact_sigma_3sigma = R_tep_true * np.sqrt(n_models) / 3.0
    exact_sigma_5sigma = R_tep_true * np.sqrt(n_models) / 5.0

    # Falsification thresholds: if GR is true (R=0), what precision is needed
    # to exclude TEP (R=R_tep_true) at a given confidence level?
    # sigma_mean = sigma / sqrt(n_models)
    # z_excl = R_tep_true / sigma_mean = R_tep_true * sqrt(n_models) / sigma
    # For 3-sigma exclusion of TEP: sigma = R_tep_true * sqrt(n_models) / 3
    # This is the same formula as detection, but now asking: if GR is true,
    # at what precision does the observation become incompatible with TEP?
    # (Because the two hypotheses are separated by R_tep_true days.)
    # The detection and exclusion thresholds are symmetric for a single test.
    sigma_excl_3sigma = exact_sigma_3sigma
    sigma_excl_5sigma = exact_sigma_5sigma

    # More useful: given current R_obs = 14.6 +/- 11.6 d, what residual would
    # falsify TEP at N-sigma?  (Threshold for |R_obs - R_tep_true| > N * sigma_R_obs)
    current_sigma_R = 11.6  # from step_07 weighted mean
    falsif_2sigma = 2.0 * current_sigma_R
    falsif_3sigma = 3.0 * current_sigma_R
    falsif_5sigma = 5.0 * current_sigma_R

    print_status(f"Thresholds for N={n_models} independent models:")
    print_status(f"  3-sigma evidence requires per-model err < {exact_sigma_3sigma:.1f} d")
    print_status(f"  5-sigma discovery requires per-model err < {exact_sigma_5sigma:.1f} d")
    print_status(f"Falsification thresholds (if GR is true, exclude TEP):")
    print_status(f"  Same sigma thresholds apply (symmetric for single test)")
    print_status(f"  At current precision (sigma_R={current_sigma_R:.1f} d):")
    print_status(f"    TEP excluded at 2-sigma if |R_obs - {R_tep_true:.1f}| > {falsif_2sigma:.1f} d")
    print_status(f"    TEP excluded at 3-sigma if |R_obs - {R_tep_true:.1f}| > {falsif_3sigma:.1f} d")

    # Plot the precision roadmap
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from scripts.utils.plot_style import set_pub_style, COLORS, FIG_SIZE
        
        set_pub_style()
        fig, ax1 = plt.subplots(figsize=FIG_SIZE)
        
        ax1.plot([r["sigma_model_days"] for r in results_list], 
                 [r["expected_z_score"] for r in results_list], 
                 color=COLORS["accent"], lw=2, label="Expected $z$-score")
        
        ax1.axhline(3, color="black", ls="--", alpha=0.5)
        ax1.axhline(5, color="black", ls=":", alpha=0.5)
        ax1.axvline(exact_sigma_3sigma, color="gray", ls="--", alpha=0.5)
        ax1.axvline(exact_sigma_5sigma, color="gray", ls=":", alpha=0.5)
        
        ax1.text(60, 3.2, r"3$\sigma$ Evidence", va="bottom", ha="right")
        ax1.text(60, 5.2, r"5$\sigma$ Discovery", va="bottom", ha="right")
        
        ax1.set_xlabel(r"Per-Model Uncertainty $\sigma_{\rm model}$ [days]")
        ax1.set_ylabel(r"Expected Ensemble Significance ($z$-score)")
        ax1.set_xlim(65, 2) # Reversed axis to show decreasing uncertainty
        ax1.set_ylim(0, 10)
        
        ax1.set_title(f"TEP Precision Roadmap (SN Refsdal, 8 Models)")
        # fig.tight_layout()
        
        out_dir = PROJECT_ROOT / "results" / "figures"
        out_dir.mkdir(exist_ok=True, parents=True)
        out_fig = out_dir / f"step_{STEP_NUM}_precision_roadmap.png"
        fig.savefig(out_fig)
        plt.close(fig)
        print_status(f"Figure saved to {out_fig}")
        
    except Exception as e:
        print_status(f"Failed to plot: {e}")

    out_json = PROJECT_ROOT / "results" / "outputs" / f"step_{STEP_NUM}_precision_roadmap.json"
    with open(out_json, "w") as f:
        json.dump({
            "step": STEP_NUM,
            "R_tep_true_days": R_tep_true,
            "n_models": n_models,
            "thresholds": {
                "sigma_3sigma_days": exact_sigma_3sigma,
                "sigma_5sigma_days": exact_sigma_5sigma
            },
            "falsification": {
                "note": "Symmetric thresholds: if GR is true, same sigma_model excludes TEP at same confidence.",
                "current_precision": {
                    "sigma_R_obs_days": current_sigma_R,
                    "tep_excluded_2sigma_if_abs_residual_deviation_gt_days": falsif_2sigma,
                    "tep_excluded_3sigma_if_abs_residual_deviation_gt_days": falsif_3sigma,
                    "tep_excluded_5sigma_if_abs_residual_deviation_gt_days": falsif_5sigma,
                },
                "sigma_model_for_tep_exclusion": {
                    "3sigma_days": sigma_excl_3sigma,
                    "5sigma_days": sigma_excl_5sigma,
                }
            },
            "data": results_list
        }, f, indent=2, default=safe_json_default)
    print_status(f"Results saved to {out_json}")
    print_status(f"Step {STEP_NUM} complete.")

if __name__ == "__main__":
    main()
