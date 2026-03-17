#!/usr/bin/env python3
"""
TEP-LENS: Step 11 - Model Dependence Robustness

Purpose:
1) Quantify sensitivity of headline evidence to inter-model dependence.
2) Run leave-one-out (LOO) stress tests.
3) Add exact/permutation-style tests for sign and weighted-mean residual.

This step does not introduce new astrophysical assumptions; it stress-tests
inference built from step_07 outputs.
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status

STEP_NUM = "11"


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def method_family(method_str: str) -> str:
    s = method_str.lower()
    if "lenstool" in s:
        return "LENSTOOL"
    if "glee" in s:
        return "GLEE"
    if "glafic" in s:
        return "GLAFIC"
    if "ltm" in s:
        return "LTM"
    if "wslap" in s:
        return "WSLAP+"
    return "OTHER"


def weighted_mean_and_sigma(values: np.ndarray, sigmas: np.ndarray):
    w = 1.0 / (sigmas**2)
    mu = float(np.sum(w * values) / np.sum(w))
    sigma_mu = float(1.0 / np.sqrt(np.sum(w)))
    return mu, sigma_mu


def exact_sign_flip_pvalue(values: np.ndarray, sigmas: np.ndarray):
    """
    Exact two-sided style tail under symmetric sign-flip null for weighted mean.
    Uses all 2^N sign assignments on |values|.
    """
    abs_vals = np.abs(values)
    obs_mu, _ = weighted_mean_and_sigma(values, sigmas)

    n = len(values)
    tails = 0
    total = 1 << n

    for mask in range(total):
        signs = np.ones(n)
        for i in range(n):
            if (mask >> i) & 1:
                signs[i] = -1.0
        trial = signs * abs_vals
        mu_trial, _ = weighted_mean_and_sigma(trial, sigmas)
        if mu_trial >= obs_mu:
            tails += 1

    p_one_sided = tails / total
    return float(p_one_sided)


def beta_binom_tail_p(n_pos: int, n_total: int, rho: float):
    """
    Correlated sign-test null via beta-binomial with p=0.5 and ICC=rho.
    rho=0 recovers ordinary binomial.
    """
    if rho <= 1e-12:
        return float(stats.binomtest(n_pos, n_total, 0.5, alternative="greater").pvalue)

    # For beta-binomial: rho = 1 / (alpha + beta + 1) when alpha=beta for p=0.5.
    # Solve alpha+beta = 1/rho - 1, and alpha=beta.
    a_plus_b = (1.0 / rho) - 1.0
    alpha = 0.5 * a_plus_b
    beta = 0.5 * a_plus_b

    xs = np.arange(n_pos, n_total + 1)
    return float(np.sum(stats.betabinom.pmf(xs, n_total, alpha, beta)))


def main():
    print_status(f"STEP {STEP_NUM}: Model Dependence Robustness", "TITLE")

    step07_path = PROJECT_ROOT / "results" / "outputs" / "step_07_observed_vs_predicted.json"
    with open(step07_path) as f:
        s07 = json.load(f)

    models = s07["per_model_results"]
    deltas = np.array([m["delta_obs_minus_pred_days"] for m in models], dtype=float)
    sigmas = np.array([m["sigma_total_days"] for m in models], dtype=float)
    names = [m["name"] for m in models]

    n_total = len(models)
    n_pos = int(np.sum(deltas > 0))

    # Baseline metrics
    mu_all, sigma_all = weighted_mean_and_sigma(deltas, sigmas)
    z_all = mu_all / sigma_all
    p_binom_all = float(stats.binomtest(n_pos, n_total, 0.5, alternative="greater").pvalue)
    p_flip_exact = exact_sign_flip_pvalue(deltas, sigmas)
    wsr = stats.wilcoxon(deltas, alternative="greater", zero_method="pratt", mode="exact")
    p_wilcoxon = float(wsr.pvalue)

    print_status(f"Baseline weighted mean residual: {mu_all:+.2f} ± {sigma_all:.2f} d (z={z_all:.2f})")
    print_status(f"Baseline sign test: {n_pos}/{n_total} positive, p={p_binom_all:.4f}")
    print_status(f"Exact sign-flip weighted-mean p(one-sided): {p_flip_exact:.4f}")
    print_status(f"Wilcoxon signed-rank p(one-sided): {p_wilcoxon:.4f}")

    # Method-family effective sample size (Kish-style weighting)
    fams = [method_family(m.get("method", "")) for m in models]
    fam_counts = {}
    for fam in fams:
        fam_counts[fam] = fam_counts.get(fam, 0) + 1

    # Each model weighted by inverse family size
    fam_weights = np.array([1.0 / fam_counts[f] for f in fams], dtype=float)
    neff_family = float((np.sum(fam_weights) ** 2) / np.sum(fam_weights**2))

    print_status(f"Method-family effective sample size (Kish proxy): N_eff={neff_family:.2f} from N={n_total}")

    # Correlation sensitivity curve for sign-test p-values
    rho_grid = np.linspace(0.0, 0.8, 17)
    corr_curve = []
    for rho in rho_grid:
        p_corr = beta_binom_tail_p(n_pos, n_total, float(rho))
        neff_rho = n_total / (1.0 + (n_total - 1.0) * rho)
        corr_curve.append(
            {
                "rho": float(rho),
                "n_eff": float(neff_rho),
                "p_beta_binom_tail": float(p_corr),
                "z_equiv": float(stats.norm.isf(min(max(p_corr, 1e-12), 1 - 1e-12))),
            }
        )

    # Leave-one-out stress tests
    loo = []
    R_tep = s07["tep_prediction"]["R_tep_prediction_days"]
    for i in range(n_total):
        mask = np.ones(n_total, dtype=bool)
        mask[i] = False

        d_sub = deltas[mask]
        s_sub = sigmas[mask]
        n_sub = len(d_sub)
        n_pos_sub = int(np.sum(d_sub > 0))

        mu_sub, sigma_sub = weighted_mean_and_sigma(d_sub, s_sub)
        z_sub = mu_sub / sigma_sub
        p_sign_sub = float(stats.binomtest(n_pos_sub, n_sub, 0.5, alternative="greater").pvalue)

        chi2_gr = float(np.sum((d_sub / s_sub) ** 2))
        chi2_tep = float(np.sum(((d_sub - R_tep) / s_sub) ** 2))
        dchi2 = chi2_gr - chi2_tep

        loo.append(
            {
                "left_out": names[i],
                "weighted_mean_days": float(mu_sub),
                "weighted_sigma_days": float(sigma_sub),
                "z_from_gr": float(z_sub),
                "n_positive": int(n_pos_sub),
                "n_total": int(n_sub),
                "p_sign": float(p_sign_sub),
                "delta_chi2_gr_minus_tep": float(dchi2),
            }
        )

    p_vals = np.array([x["p_sign"] for x in loo])
    z_vals = np.array([x["z_from_gr"] for x in loo])
    dchi2_vals = np.array([x["delta_chi2_gr_minus_tep"] for x in loo])

    loo_summary = {
        "p_sign_min": float(np.min(p_vals)),
        "p_sign_max": float(np.max(p_vals)),
        "z_min": float(np.min(z_vals)),
        "z_max": float(np.max(z_vals)),
        "delta_chi2_min": float(np.min(dchi2_vals)),
        "delta_chi2_max": float(np.max(dchi2_vals)),
        "n_cases_tep_better": int(np.sum(dchi2_vals > 0)),
        "n_cases_total": int(len(dchi2_vals)),
    }

    print_status("LOO robustness summary:")
    print_status(f"  Sign-test p range: {loo_summary['p_sign_min']:.4f} to {loo_summary['p_sign_max']:.4f}")
    print_status(f"  Weighted-mean z range: {loo_summary['z_min']:.2f} to {loo_summary['z_max']:.2f}")
    print_status(
        f"  Δχ² (GR-TEP) range: {loo_summary['delta_chi2_min']:+.2f} to {loo_summary['delta_chi2_max']:+.2f}"
    )

    # Plot
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from scripts.utils.plot_style import set_pub_style, COLORS, FIG_SIZE

        set_pub_style()

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=FIG_SIZE)

        # Correlation sensitivity
        ax1.plot(
            [c["rho"] for c in corr_curve],
            [c["p_beta_binom_tail"] for c in corr_curve],
            marker="o",
            lw=1.8,
        )
        ax1.axhline(p_binom_all, color="gray", ls="--", lw=1.2, label="Independent binomial p")
        ax1.set_xlabel("Assumed inter-model ICC (rho)")
        ax1.set_ylabel("One-sided sign-test p-value")
        ax1.set_title("Dependence sensitivity (beta-binomial)")
        ax1.set_ylim(0, 1)
        ax1.grid(alpha=0.3, ls=":")
        ax1.legend()

        # LOO sign p-values
        order = np.argsort(p_vals)
        ax2.barh(np.arange(n_total), p_vals[order], color="#4c78a8")
        ax2.axvline(0.05, color="black", ls="--", lw=1.0)
        ax2.set_yticks(np.arange(n_total))
        ax2.set_yticklabels([loo[i]["left_out"] for i in order], )
        ax2.set_xlabel("LOO one-sided sign-test p-value")
        ax2.set_title("Leave-one-out stress test")
        ax2.grid(alpha=0.3, ls=":", axis="x")

        fig.tight_layout()
        out_fig = PROJECT_ROOT / "results" / "figures" / f"step_{STEP_NUM}_model_dependence.png"
        fig.savefig(out_fig)
        plt.close(fig)
        print_status(f"Figure saved to {out_fig}")
    except Exception as e:
        print_status(f"Plotting failed: {e}", "ERROR")
        out_fig = None

    # Save
    out = {
        "step": STEP_NUM,
        "status": "success",
        "baseline": {
            "weighted_mean_days": float(mu_all),
            "weighted_sigma_days": float(sigma_all),
            "z_from_gr": float(z_all),
            "n_positive": int(n_pos),
            "n_total": int(n_total),
            "p_binomial_independent": float(p_binom_all),
            "p_weighted_mean_signflip_exact": float(p_flip_exact),
            "p_wilcoxon_signed_rank": float(p_wilcoxon),
        },
        "effective_sample_size": {
            "n_total": int(n_total),
            "method_family_counts": fam_counts,
            "n_eff_family_kish": float(neff_family),
        },
        "dependence_sensitivity_curve": corr_curve,
        "leave_one_out": loo,
        "leave_one_out_summary": loo_summary,
        "figure": str(out_fig) if out_fig else None,
    }

    out_path = PROJECT_ROOT / "results" / "outputs" / f"step_{STEP_NUM}_model_dependence.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=safe_json_default)

    print_status(f"Results saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
