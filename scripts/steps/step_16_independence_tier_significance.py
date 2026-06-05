#!/usr/bin/env python3
"""
TEP-LENS: Step 16 - Independence-Tier Significance

Evaluates the overall significance of the evidence.
Because the various tests for SN Refsdal (Wilcoxon sign test, weighted mean,
Pearson correlation, alpha inference) all rely fundamentally on the single
anomalous SX arrival time, they are highly correlated.

Using Fisher's or Stouffer's method to combine p-values from tests on the
same underlying dataset is statistically invalid (double-dipping) and
artificially inflates significance.

The most defensible approach is to select the single most robust non-parametric
test as the headline significance for the system, supported by the consistency
of the other metrics.
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status
from scripts.utils.plot_style import set_pub_style, COLORS, FIG_SIZE

STEP_NUM = "16"


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def main():
    print_status(f"STEP {STEP_NUM}: Independence-Tier Significance", "TITLE")

    s07_path = PROJECT_ROOT / "results" / "outputs" / "step_07_observed_vs_predicted.json"
    s08_path = PROJECT_ROOT / "results" / "outputs" / "step_08_new_evidence.json"
    s11_path = PROJECT_ROOT / "results" / "outputs" / "step_11_model_dependence.json"
    s15_path = PROJECT_ROOT / "results" / "outputs" / "step_15_external_inflation.json"

    with open(s07_path) as f:
        s07 = json.load(f)
    with open(s08_path) as f:
        s08 = json.load(f)
    with open(s11_path) as f:
        s11 = json.load(f)
    with open(s15_path) as f:
        s15 = json.load(f)

    # Primary independence test: Wilcoxon signed-rank on blind models only
    p_sign = float(s07["binomial_sign_test"]["p_wilcoxon_signed_rank_blind"])
    # Primary correlation-aware bound: exact family-sign-flip test.
    p_perm = float(s11["permutation_test_wilcoxon"]["blind_7"]["p_value_one_sided"])
    # Sensitivity exploration: block-bootstrap (demoted from operational primary)
    p_boot = float(s11["block_bootstrap_wilcoxon"]["blind_7"]["p_median"])
    p_ext = float(s15["scenarios"]["kappa50"]["p_z_one_sided"])
    p_pearson = float(s08["test_A_delay_mu_correlation"]["pearson_p_onesided"])
    p_alpha = float(s08["test_C_alpha_inference"]["p_vs_gr_null_onesided"])

    # Calculate z-scores for individual tests
    z_sign = float(stats.norm.isf(p_sign))
    z_perm = float(stats.norm.isf(p_perm))
    z_boot = float(stats.norm.isf(p_boot))
    z_ext = float(stats.norm.isf(p_ext))
    z_pearson = float(stats.norm.isf(p_pearson))
    z_alpha = float(stats.norm.isf(p_alpha))

    print_status("Individual Correlated Tests (SN Refsdal SX):")
    print_status(f"  Wilcoxon signed-rank (blind): z={z_sign:.2f}, p={p_sign:.4f} (Independence primary)")
    boot = s11["block_bootstrap_wilcoxon"]["blind_7"]

    print_status(f"  Exact family-sign-flip (blind): z={z_perm:.2f}, p={p_perm:.4f} (Correlation-aware primary)")
    print_status(
        f"  Block-bootstrap (blind):       z={z_boot:.2f}, p={p_boot:.4f} "
        f"[{boot['p_16']:.4f}, {boot['p_84']:.4f}] (Sensitivity)"
    )
    print_status(f"  Pearson delay-mu:             z={z_pearson:.2f}, p={p_pearson:.4f}")
    print_status(f"  Ext-informed w-mean:          z={z_ext:.2f}, p={p_ext:.4f}")
    print_status(f"  Alpha vs zero:              z={z_alpha:.2f}, p={p_alpha:.4f}")
    print_status("\nNOTE: These tests are fundamentally correlated. Combining them")
    print_status("via Fisher/Stouffer is invalid. The conservative significance")
    print_status(f"is the most robust single test: z={z_perm:.2f} (Exact family-sign-flip).")

    out_fig = None
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        set_pub_style()

        labels = [
            "Exact Perm\n(Primary)",
            "Wilcoxon\n(Indep.)",
            "Block-Boot\n(Sens.)",
            "Pearson\nCorr.",
            "Ext-Inf\nMean",
            "Alpha vs\nZero",
        ]
        zvals = [z_perm, z_sign, z_boot, z_pearson, z_ext, z_alpha]
        bar_colors = [COLORS['red'], COLORS['accent'], COLORS['model'], COLORS['tep'], COLORS['tep'], COLORS['tep']]

        fig, ax = plt.subplots(figsize=FIG_SIZE)
        ax.bar(labels, zvals, color=bar_colors)
        ax.axhline(1.96, ls="--", color="black", lw=1.0, label="95% CI (1.96σ)")
        ax.axhline(3.0, ls=":", color="black", lw=1.0, label="99.7% CI (3.0σ)")
        ax.set_ylabel("Equivalent z (one-sided)")
        ax.set_title("Step 16: Correlated SN Refsdal Evidence Strands")
        ax.legend(loc="upper right")
        ax.grid(axis="y", alpha=0.3, ls=":")
        fig.tight_layout()

        out_fig = PROJECT_ROOT / "results" / "figures" / f"step_{STEP_NUM}_tier_significance.png"
        fig.savefig(out_fig)
        plt.close(fig)
        print_status(f"Figure saved to {out_fig}")
    except Exception as e:
        print_status(f"Plotting failed: {e}", "ERROR")

    out = {
        "step": STEP_NUM,
        "status": "success",
        "methodology_note": "Fisher/Stouffer combination removed due to invalid double-dipping on correlated data.",
        "headline_significance": {
            "test": "Exact family-sign-flip test blind",
            "z_score": z_perm,
            "p_value": p_perm,
            "note": "Correlation-aware primary result; kept for downstream compatibility.",
        },
        "test_hierarchy": {
            "tier_1a_primary_independence": {
                "test": "Wilcoxon signed-rank blind",
                "z_score": z_sign,
                "p_value": p_sign,
            },
            "tier_1b_primary_correlation_aware": {
                "test": "Exact family-sign-flip test blind",
                "z_score": z_perm,
                "p_value": p_perm,
            },
            "tier_2_sensitivity_exploration": {
                "test": "Block-bootstrap Wilcoxon blind",
                "z_score": z_boot,
                "p_value": p_boot,
                "caveat": "Demoted from operational primary; retained for sensitivity exploration.",
            },
        },
        "individual_tests": {
            "exact_permutation": {"z": z_perm, "p": p_perm},
            "wilcoxon": {"z": z_sign, "p": p_sign},
            "block_bootstrap": {"z": z_boot, "p": p_boot},
            "pearson": {"z": z_pearson, "p": p_pearson},
            "external_mean": {"z": z_ext, "p": p_ext},
            "alpha_inference": {"z": z_alpha, "p": p_alpha}
        },
        "figure": str(out_fig) if out_fig else None,
    }

    out_path = PROJECT_ROOT / "results" / "outputs" / f"step_{STEP_NUM}_tier_significance.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=safe_json_default)

    print_status(f"Results saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")

if __name__ == "__main__":
    main()
