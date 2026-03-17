#!/usr/bin/env python3
"""
TEP-LENS: Step 12 - Microlensing Nuisance Robustness

Monte Carlo propagation of flux-ratio perturbations to test how strongly
microlensing-scale systematics can change key TEP inferences.

Outputs per perturbation level:
- Distribution of R_TEP(alpha=-0.05) for S1-S4-SX loop
- Probability TEP still improves chi^2 vs GR
- Probability inferred alpha remains < 0
- Robustness interval for R_TEP and alpha_inferred
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status

STEP_NUM = "12"


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def compute_r_tep_alpha05(mu_norm):
    """
    Compute S1-S4-SX closure residual under alpha=-0.05 using the same loop
    decomposition used by prior steps.
    """
    alpha = -0.05
    gamma = {k: 1.0 + alpha * np.log10(v) for k, v in mu_norm.items()}

    # Loop (S1, S4, SX) using arrival times relative to S1:
    # t(S1)=0, t(S4)=20.3, t(SX)=376.0
    dt_14 = 20.3
    dt_4x = 355.7
    dt_x1 = -376.0

    r = (gamma["S1"] - 1.0) * dt_14 + (gamma["S4"] - 1.0) * dt_4x + (gamma["SX"] - 1.0) * dt_x1
    return float(r)


def main():
    print_status(f"STEP {STEP_NUM}: Microlensing Robustness Monte Carlo", "TITLE")

    step07_path = PROJECT_ROOT / "results" / "outputs" / "step_07_observed_vs_predicted.json"
    with open(step07_path) as f:
        s07 = json.load(f)

    catalog_path = PROJECT_ROOT / "data" / "raw" / "sn_lensing" / "lensed_sn_catalog.json"
    with open(catalog_path) as f:
        catalog = json.load(f)

    fluxes = catalog["sn_refsdal"]["magnification_proxies"]["flux_total"]

    mu_base = {k: float(fluxes[k]["value"]) for k in ["S1", "S2", "S3", "S4", "SX"]}
    mu_err = {k: float(fluxes[k]["err"]) for k in ["S1", "S2", "S3", "S4", "SX"]}

    deltas = np.array([m["delta_obs_minus_pred_days"] for m in s07["per_model_results"]], dtype=float)
    sigmas = np.array([m["sigma_total_days"] for m in s07["per_model_results"]], dtype=float)

    r_obs = float(s07["weighted_mean_residual"]["R_obs_days"])
    r_tep_nominal = float(s07["tep_prediction"]["R_tep_alpha05_days"])

    rng = np.random.default_rng(42)
    n_draws = 20000
    frac_levels = [0.10, 0.20, 0.30]

    summaries = []

    for frac in frac_levels:
        r_draws = np.zeros(n_draws)
        alpha_inf_draws = np.zeros(n_draws)
        dchi2_draws = np.zeros(n_draws)

        for i in range(n_draws):
            mu_pert = {}
            for k in mu_base:
                # Combine quoted photometric error and nuisance microlensing fraction.
                sigma_abs = np.sqrt(mu_err[k] ** 2 + (frac * mu_base[k]) ** 2)
                draw = rng.normal(mu_base[k], sigma_abs)
                mu_pert[k] = max(draw, 1e-4)

            mu_mean = np.mean(list(mu_pert.values()))
            mu_norm = {k: mu_pert[k] / mu_mean for k in mu_pert}

            r_tep = compute_r_tep_alpha05(mu_norm)
            r_draws[i] = r_tep

            r_unit = r_tep / 0.05 if abs(r_tep) > 1e-12 else np.nan
            alpha_inf_draws[i] = r_obs / r_unit if np.isfinite(r_unit) else np.nan

            chi2_gr = np.sum((deltas / sigmas) ** 2)
            chi2_tep = np.sum(((deltas - r_tep) / sigmas) ** 2)
            dchi2_draws[i] = chi2_gr - chi2_tep

        finite_alpha = np.isfinite(alpha_inf_draws)
        alpha_clean = alpha_inf_draws[finite_alpha]

        summary = {
            "microlensing_fraction": float(frac),
            "n_draws": int(n_draws),
            "R_tep_median_days": float(np.median(r_draws)),
            "R_tep_p16_days": float(np.percentile(r_draws, 16)),
            "R_tep_p84_days": float(np.percentile(r_draws, 84)),
            "R_tep_p025_days": float(np.percentile(r_draws, 2.5)),
            "R_tep_p975_days": float(np.percentile(r_draws, 97.5)),
            "P_R_tep_positive": float(np.mean(r_draws > 0)),
            "alpha_inferred_median": float(np.median(alpha_clean)),
            "alpha_inferred_p16": float(np.percentile(alpha_clean, 16)),
            "alpha_inferred_p84": float(np.percentile(alpha_clean, 84)),
            "P_alpha_inferred_positive": float(np.mean(alpha_clean > 0)),
            "delta_chi2_median": float(np.median(dchi2_draws)),
            "delta_chi2_p16": float(np.percentile(dchi2_draws, 16)),
            "delta_chi2_p84": float(np.percentile(dchi2_draws, 84)),
            "P_delta_chi2_gt_0": float(np.mean(dchi2_draws > 0)),
        }
        summaries.append(summary)

        print_status(
            f"frac={frac:.2f}: R_TEP={summary['R_tep_median_days']:.2f} "
            f"[{summary['R_tep_p16_days']:.2f}, {summary['R_tep_p84_days']:.2f}] d; "
            f"P(Δχ²>0)={summary['P_delta_chi2_gt_0']:.3f}"
        )

    # Plot
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))

        x = [s["microlensing_fraction"] for s in summaries]
        r_med = [s["R_tep_median_days"] for s in summaries]
        r_lo = [s["R_tep_median_days"] - s["R_tep_p16_days"] for s in summaries]
        r_hi = [s["R_tep_p84_days"] - s["R_tep_median_days"] for s in summaries]

        axs[0].errorbar(x, r_med, yerr=[r_lo, r_hi], marker="o", capsize=4, lw=1.8)
        axs[0].axhline(r_tep_nominal, ls="--", color="gray", lw=1.2, label="Nominal R_TEP")
        axs[0].set_xlabel("Microlensing nuisance fraction")
        axs[0].set_ylabel("R_TEP (days)")
        axs[0].set_title("R_TEP robustness to flux perturbations")
        axs[0].grid(alpha=0.3, ls=":")
        axs[0].legend()

        p_chi = [s["P_delta_chi2_gt_0"] for s in summaries]
        p_alpha = [s["P_alpha_inferred_positive"] for s in summaries]
        axs[1].plot(x, p_chi, marker="o", lw=1.8, label="P(TEP better fit: Δχ²>0)")
        axs[1].plot(x, p_alpha, marker="s", lw=1.8, label="P(alpha_inferred>0)")
        axs[1].set_ylim(0, 1)
        axs[1].set_xlabel("Microlensing nuisance fraction")
        axs[1].set_ylabel("Robustness probability")
        axs[1].set_title("Inference robustness")
        axs[1].grid(alpha=0.3, ls=":")
        axs[1].legend()

        # fig.tight_layout()
        out_fig = PROJECT_ROOT / "results" / "figures" / f"step_{STEP_NUM}_microlensing_robustness.png"
        fig.savefig(out_fig)
        plt.close(fig)
        print_status(f"Figure saved to {out_fig}")
    except Exception as e:
        print_status(f"Plotting failed: {e}", "ERROR")
        out_fig = None

    out = {
        "step": STEP_NUM,
        "status": "success",
        "nominal": {
            "R_obs_weighted_days": r_obs,
            "R_tep_alpha05_nominal_days": r_tep_nominal,
        },
        "design": {
            "draws_per_level": n_draws,
            "microlensing_fraction_levels": frac_levels,
            "note": "Flux perturbations combine quoted photometric errors and fractional nuisance terms.",
        },
        "summaries": summaries,
        "figure": str(out_fig) if out_fig else None,
    }

    out_path = PROJECT_ROOT / "results" / "outputs" / f"step_{STEP_NUM}_microlensing_robustness.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=safe_json_default)

    print_status(f"Results saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
