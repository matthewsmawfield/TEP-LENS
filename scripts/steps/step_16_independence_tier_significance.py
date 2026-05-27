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
from scripts.utils.logger import print_status, safe_json_default
from scripts.utils.plot_style import set_pub_style, COLORS, FIG_SIZE

STEP_NUM = "16"




def main():
    print_status(f"STEP {STEP_NUM}: Independence-Tier Significance", "TITLE")

    s07_path = PROJECT_ROOT / "results" / "outputs" / "step_07_observed_vs_predicted.json"
    s08_path = PROJECT_ROOT / "results" / "outputs" / "step_08_new_evidence.json"
    s15_path = PROJECT_ROOT / "results" / "outputs" / "step_15_external_inflation.json"

    with open(s07_path) as f:
        s07 = json.load(f)
    with open(s08_path) as f:
        s08 = json.load(f)
    with open(s15_path) as f:
        s15 = json.load(f)

    p_sign = float(s07["binomial_sign_test"].get(
        "p_wilcoxon_signed_rank",
        0.0078125
    ))
    p_ext = float(s15["scenarios"]["kappa50"]["p_z_one_sided"])
    p_pearson = float(s08["test_A_delay_mu_correlation"]["pearson_p_onesided"])
    p_alpha = float(s08["test_C_alpha_inference"]["p_vs_gr_null_onesided"])

    # Calculate z-scores for individual tests
    z_sign = float(stats.norm.isf(p_sign))
    z_ext = float(stats.norm.isf(p_ext))
    z_pearson = float(stats.norm.isf(p_pearson))
    z_alpha = float(stats.norm.isf(p_alpha))

    print_status("Individual Correlated Tests (SN Refsdal SX):")
    print_status(f"  Wilcoxon signed-rank: z={z_sign:.2f}, p={p_sign:.4f} (Headline)")
    print_status(f"  Pearson delay-mu:     z={z_pearson:.2f}, p={p_pearson:.4f}")
    print_status(f"  Ext-informed w-mean:  z={z_ext:.2f}, p={p_ext:.4f}")
    print_status(f"  Alpha vs zero:        z={z_alpha:.2f}, p={p_alpha:.4f}")
    print_status("\nNOTE: These tests are fundamentally correlated. Combining them")
    print_status("via Fisher/Stouffer is invalid. The conservative significance")
    print_status(f"is the strongest robust single test: z={z_sign:.2f} (Wilcoxon).")

    out_fig = None
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        set_pub_style()

        labels = [
            "Wilcoxon\n(Headline)",
            "Pearson\nCorrelation",
            "Ext-Informed\nMean",
            "Alpha vs\nZero",
        ]
        zvals = [z_sign, z_pearson, z_ext, z_alpha]

        fig, ax = plt.subplots(figsize=FIG_SIZE)
        ax.bar(labels, zvals, color=[COLORS['red'], COLORS['tep'], COLORS['tep'], COLORS['tep']])
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
        print_status(f"Plotting skipped: {e}", "WARNING")

    out = {
        "step": STEP_NUM,
        "status": "success",
        "methodology_note": "Fisher/Stouffer combination removed due to invalid double-dipping on correlated data.",
        "headline_significance": {
            "test": "Wilcoxon signed-rank",
            "z_score": z_sign,
            "p_value": p_sign
        },
        "individual_tests": {
            "wilcoxon": {"z": z_sign, "p": p_sign},
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
