#!/usr/bin/env python3
"""
TEP-LENS: Step 17 - Directional Odds Expansion

Adds a Bayes-factor directional-odds layer for sign-based evidence.
This complements p-value reporting with an odds interpretation while
remaining transparent about data dependence and sample size.
"""

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.special import beta as beta_fn
from scipy.special import betainc

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status

STEP_NUM = "17"


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def method_family(method_name: str) -> str:
    m = method_name.lower()
    if "glafic" in m:
        return "GLAFIC"
    if "ltm" in m:
        return "LTM"
    if "wslap" in m:
        return "WSLAP+"
    if "glee" in m:
        return "GLEE"
    if "lenstool" in m:
        return "LENSTOOL"
    return "OTHER"


def directional_bf_from_k_n(k_pos: int, n_total: int):
    """
    Compare:
      H0: sign-positive probability p = 0.5
      H1: p ~ Uniform(0.5, 1.0)  [directional alternative]

    Data likelihood for k positives out of n:
      L0 = 0.5^n
      L1 = 2 * Integral_{0.5}^1 p^k (1-p)^(n-k) dp
         = 2 * Beta(k+1, n-k+1) * (1 - I_0.5(k+1, n-k+1))
    """
    if n_total <= 0:
        return {
            "k_positive": int(k_pos),
            "n_total": int(n_total),
            "binomial_p_one_sided": None,
            "likelihood_h0": None,
            "likelihood_h1": None,
            "bf10_directional": None,
            "posterior_h1_equal_prior": None,
        }

    a = k_pos + 1
    b = n_total - k_pos + 1

    l0 = 0.5**n_total
    tail_from_half = 1.0 - betainc(a, b, 0.5)
    l1 = 2.0 * beta_fn(a, b) * tail_from_half

    bf10 = l1 / l0 if l0 > 0 else np.inf
    post_h1 = bf10 / (1.0 + bf10)

    p_binom = float(stats.binomtest(k_pos, n_total, 0.5, alternative="greater").pvalue)

    return {
        "k_positive": int(k_pos),
        "n_total": int(n_total),
        "binomial_p_one_sided": p_binom,
        "likelihood_h0": float(l0),
        "likelihood_h1": float(l1),
        "bf10_directional": float(bf10),
        "posterior_h1_equal_prior": float(post_h1),
    }


def scenario_from_signs(sign_values: np.ndarray):
    nz = sign_values[np.abs(sign_values) > 1e-12]
    k = int(np.sum(nz > 0))
    n = int(nz.size)
    out = directional_bf_from_k_n(k, n)
    out["n_zero_excluded"] = int(sign_values.size - n)
    return out


def main():
    print_status(f"STEP {STEP_NUM}: Directional Odds Expansion", "TITLE")

    s07_path = PROJECT_ROOT / "results" / "outputs" / "step_07_observed_vs_predicted.json"
    with open(s07_path) as f:
        s07 = json.load(f)

    per_model = s07["per_model_results"]

    deltas_all = np.array([float(r["delta_obs_minus_pred_days"]) for r in per_model], dtype=float)
    blind_mask = np.array([bool(r.get("blind", False)) for r in per_model], dtype=bool)

    # Scenario 1: all models, non-zero residuals only (Wilcoxon-compatible sign set)
    all_nonzero = scenario_from_signs(deltas_all)

    # Scenario 2: blind models only, non-zero residuals
    blind_nonzero = scenario_from_signs(deltas_all[blind_mask])

    # Scenario 3: method-family collapsed residual signs (inverse-variance weighted mean per family)
    fam_records = {}
    for rec in per_model:
        fam = method_family(str(rec.get("method", "")))
        fam_records.setdefault(fam, []).append(rec)

    fam_names = sorted(fam_records.keys())
    fam_means = []
    fam_debug = []
    for fam in fam_names:
        items = fam_records[fam]
        vals = np.array([float(x["delta_obs_minus_pred_days"]) for x in items], dtype=float)
        sigs = np.array([float(x["sigma_total_days"]) for x in items], dtype=float)
        w = 1.0 / np.maximum(sigs, 1e-9) ** 2
        mean_fam = float(np.sum(w * vals) / np.sum(w))
        fam_means.append(mean_fam)
        fam_debug.append(
            {
                "family": fam,
                "n_models": int(len(items)),
                "weighted_mean_residual_days": mean_fam,
            }
        )

    family_nonzero = scenario_from_signs(np.array(fam_means, dtype=float))

    print_status(
        "All-model nonzero sign set: "
        f"k={all_nonzero['k_positive']}/{all_nonzero['n_total']}, "
        f"p={all_nonzero['binomial_p_one_sided']:.4f}, "
        f"BF10={all_nonzero['bf10_directional']:.2f}"
    )
    print_status(
        "Blind-only nonzero sign set: "
        f"k={blind_nonzero['k_positive']}/{blind_nonzero['n_total']}, "
        f"p={blind_nonzero['binomial_p_one_sided']:.4f}, "
        f"BF10={blind_nonzero['bf10_directional']:.2f}"
    )
    print_status(
        "Method-family collapsed sign set: "
        f"k={family_nonzero['k_positive']}/{family_nonzero['n_total']}, "
        f"p={family_nonzero['binomial_p_one_sided']:.4f}, "
        f"BF10={family_nonzero['bf10_directional']:.2f}"
    )

    out_fig = None
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from scripts.utils.plot_style import set_pub_style, COLORS, FIG_SIZES, save_fig

        labels = ["All nonzero", "Blind nonzero", "Family collapsed"]
        bfs = [
            all_nonzero["bf10_directional"],
            blind_nonzero["bf10_directional"],
            family_nonzero["bf10_directional"],
        ]
        log10_bf = [math.log10(max(b, 1e-12)) for b in bfs]

        set_pub_style()
        fig, ax = plt.subplots(figsize=FIG_SIZES["web_standard"])
        bar_colors = [COLORS["tep"], COLORS["observed"], COLORS["model"]]
        ax.bar(labels, log10_bf, color=bar_colors, edgecolor=COLORS["text"], linewidth=0.6)
        ax.axhline(math.log10(3.0), ls="--", color=COLORS["gr"], lw=1.0)
        ax.axhline(math.log10(10.0), ls=":", color=COLORS["red"], lw=1.0)
        ax.set_ylabel(r"$\log_{10}(\mathrm{BF}_{10})$ for directional sign model")
        ax.set_title("Step 17: Directional Bayes-factor odds by evidence subset")

        for i, (bf, y) in enumerate(zip(bfs, log10_bf)):
            ax.text(i, y + 0.04, f"BF={bf:.1f}", ha="center", va="bottom")

        out_fig = PROJECT_ROOT / "results" / "figures" / f"step_{STEP_NUM}_directional_odds.png"
        save_fig(fig, out_fig)
        print_status(f"Figure saved to {out_fig}")
    except Exception as e:
        print_status(f"Plotting failed: {e}", "ERROR")

    out = {
        "step": STEP_NUM,
        "status": "success",
        "description": (
            "Directional Bayes-factor expansion for sign-based evidence. "
            "H0: p(sign+)=0.5; H1: p(sign+)~Uniform(0.5,1)."
        ),
        "scenarios": {
            "all_models_nonzero": all_nonzero,
            "blind_models_nonzero": blind_nonzero,
            "method_family_collapsed": family_nonzero,
        },
        "method_family_details": fam_debug,
        "figure": str(out_fig) if out_fig else None,
    }

    out_path = PROJECT_ROOT / "results" / "outputs" / f"step_{STEP_NUM}_directional_odds.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=safe_json_default)

    print_status(f"Results saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
