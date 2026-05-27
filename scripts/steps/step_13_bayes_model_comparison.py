#!/usr/bin/env python3
"""
TEP-LENS: Step 13 - Bayesian GR vs Proxy-Model Comparison

Implements a conservative empirical-Bayesian comparison using the per-model
residuals from step_07:

Data model:
  delta_i ~ Normal(mu_shift, sqrt(sigma_i^2 + tau^2))

Where:
  - sigma_i are reported model+measurement errors from step_07
  - tau is extra between-model overdispersion (hierarchical nuisance)
  - mu_shift is hypothesis-dependent:
      H_GR         : mu_shift = mu_bias
      H_TEP_fixed  : mu_shift = mu_bias + R_TEP(alpha=-0.055)
      H_TEP_free   : mu_shift = mu_bias + alpha * R_unit

Priors (transparent and weakly informative):
  mu_bias ~ Normal(0, 40 d)
  tau ~ HalfNormal(20 d)
  alpha ~ Normal(-0.055, 0.05)  [proxy-model free-alpha model only; centered on empirical SN Refsdal measurement]

Outputs:
  - Marginal log-evidence per model (numerical quadrature / grid integration)
  - Bayes factors BF(proxy_fixed/GR), BF(proxy_free/GR)
  - Posterior summaries for tau and alpha (free model)
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.special import logsumexp

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status

STEP_NUM = "13"


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def loglike_residuals(deltas, sigmas, mu_shift, tau):
    var = sigmas**2 + tau**2
    return -0.5 * np.sum(np.log(2.0 * np.pi * var) + ((deltas - mu_shift) ** 2) / var)


def loglike_mu_grid(deltas, sigmas, mu_grid, tau, offset=0.0):
    """Vectorized log-likelihood on a 1D mu grid for fixed tau and additive offset."""
    var = sigmas**2 + tau**2
    log_norm = np.sum(np.log(2.0 * np.pi * var))
    mu_eff = mu_grid + offset
    resid = deltas[None, :] - mu_eff[:, None]
    quad = np.sum((resid**2) / var[None, :], axis=1)
    return -0.5 * (log_norm + quad)


def evaluate_scenario(
    deltas,
    sigmas,
    r_tep,
    r_unit,
    mu_grid,
    tau_grid,
    alpha_grid,
    mu_prior_loc,
    mu_prior_scale,
    tau_prior_scale,
    alpha_prior_loc,
    alpha_prior_scale,
):
    d_mu = float(mu_grid[1] - mu_grid[0])
    d_tau = float(tau_grid[1] - tau_grid[0])
    d_alpha = float(alpha_grid[1] - alpha_grid[0])

    mu_prior = stats.norm.logpdf(mu_grid, loc=mu_prior_loc, scale=mu_prior_scale)
    tau_prior = stats.halfnorm.logpdf(tau_grid, loc=0.0, scale=tau_prior_scale)
    alpha_prior = stats.norm.logpdf(alpha_grid, loc=alpha_prior_loc, scale=alpha_prior_scale)

    # H_GR
    log_post_gr = np.empty((len(tau_grid), len(mu_grid)))
    for it, tau in enumerate(tau_grid):
        ll = loglike_mu_grid(deltas, sigmas, mu_grid, tau, offset=0.0)
        log_post_gr[it, :] = ll + mu_prior + tau_prior[it]
    logZ_gr = logsumexp(log_post_gr) + np.log(d_mu) + np.log(d_tau)

    # H_TEP_fixed
    log_post_tf = np.empty((len(tau_grid), len(mu_grid)))
    for it, tau in enumerate(tau_grid):
        ll = loglike_mu_grid(deltas, sigmas, mu_grid, tau, offset=r_tep)
        log_post_tf[it, :] = ll + mu_prior + tau_prior[it]
    logZ_tf = logsumexp(log_post_tf) + np.log(d_mu) + np.log(d_tau)

    # H_TEP_free
    log_terms = []
    log_marg_tau = np.full(len(tau_grid), -np.inf)
    log_marg_alpha = np.full(len(alpha_grid), -np.inf)

    for it, tau in enumerate(tau_grid):
        shifts = alpha_grid * r_unit
        var = sigmas**2 + tau**2
        log_norm = np.sum(np.log(2.0 * np.pi * var))
        mu_eff = shifts[:, None] + mu_grid[None, :]
        resid = deltas[None, None, :] - mu_eff[:, :, None]
        quad = np.sum((resid**2) / var[None, None, :], axis=2)
        ll_mat = -0.5 * (log_norm + quad)
        mat = ll_mat + mu_prior[None, :] + alpha_prior[:, None] + tau_prior[it]

        lt_tau = logsumexp(mat) + np.log(d_alpha) + np.log(d_mu)
        log_terms.append(lt_tau)
        log_marg_tau[it] = lt_tau

        log_alpha_given_tau = logsumexp(mat, axis=1) + np.log(d_mu)
        log_marg_alpha = np.logaddexp(log_marg_alpha, log_alpha_given_tau)

    logZ_tfree = logsumexp(np.array(log_terms)) + np.log(d_tau)

    # Posterior summaries
    log_marg_tau -= logsumexp(log_marg_tau)
    p_tau = np.exp(log_marg_tau)

    log_marg_alpha -= logsumexp(log_marg_alpha)
    p_alpha = np.exp(log_marg_alpha)

    tau_mean = float(np.sum(tau_grid * p_tau))
    tau_p16 = float(np.interp(0.16, np.cumsum(p_tau), tau_grid))
    tau_p84 = float(np.interp(0.84, np.cumsum(p_tau), tau_grid))

    alpha_mean = float(np.sum(alpha_grid * p_alpha))
    alpha_p16 = float(np.interp(0.16, np.cumsum(p_alpha), alpha_grid))
    alpha_p84 = float(np.interp(0.84, np.cumsum(p_alpha), alpha_grid))

    log_bf_tf_gr = float(logZ_tf - logZ_gr)
    log_bf_tfree_gr = float(logZ_tfree - logZ_gr)

    bf_tf_gr = float(np.exp(min(log_bf_tf_gr, 700.0)))
    bf_tfree_gr = float(np.exp(min(log_bf_tfree_gr, 700.0)))

    return {
        "logZ_gr": float(logZ_gr),
        "logZ_tep_fixed": float(logZ_tf),
        "logZ_tep_free": float(logZ_tfree),
        "log_bf_tep_fixed_over_gr": float(log_bf_tf_gr),
        "log_bf_tep_free_over_gr": float(log_bf_tfree_gr),
        "bf_tep_fixed_over_gr": float(bf_tf_gr),
        "bf_tep_free_over_gr": float(bf_tfree_gr),
        "tau_mean_days": tau_mean,
        "tau_p16_days": tau_p16,
        "tau_p84_days": tau_p84,
        "alpha_mean": alpha_mean,
        "alpha_p16": alpha_p16,
        "alpha_p84": alpha_p84,
        "posterior_tau_grid": p_tau,
        "posterior_alpha_grid": p_alpha,
    }


def main():
    print_status(f"STEP {STEP_NUM}: Bayesian GR vs TEP Comparison", "TITLE")

    step07_path = PROJECT_ROOT / "results" / "outputs" / "step_07_observed_vs_predicted.json"
    with open(step07_path) as f:
        s07 = json.load(f)

    deltas = np.array([m["delta_obs_minus_pred_days"] for m in s07["per_model_results"]], dtype=float)
    sigmas = np.array([m["sigma_total_days"] for m in s07["per_model_results"]], dtype=float)

    # TEP prediction (fixed alpha model)
    # The step_07 JSON output key was updated from R_tep_alpha05_days to R_tep_prediction_days
    if "R_tep_prediction_days" in s07["tep_prediction"]:
        r_tep = float(s07["tep_prediction"]["R_tep_prediction_days"])
    else:
        r_tep = float(s07["tep_prediction"]["R_tep_alpha05_days"]) # Fallback if using old file
    
    # Calculate R_tep_unit (days per unit alpha)
    # R_tep = R_tep_unit * alpha_ref
    alpha_ref = float(s07["tep_prediction"]["alpha_ref"])
    r_unit = float(s07["tep_prediction"]["R_tep_unit_days_per_alpha"])

    # Priors / grids
    mu_grid = np.linspace(-80.0, 80.0, 401)   # days
    tau_grid = np.linspace(0.0, 60.0, 301)    # days
    # Alpha grid for negative coupling (empirical lensing-sector value)
    alpha_grid = np.linspace(-0.20, 0.10, 401)

    scenarios = {
        "baseline": {
            "mu_prior_loc": 0.0,
            "mu_prior_scale": 40.0,
            "tau_prior_scale": 20.0,
            "alpha_prior_loc": 0.0,  # GR-centered null prior for free-alpha model
            "alpha_prior_scale": 0.15,
            "note": "Proper GR-null prior for free-alpha model (alpha ~ N(0, 0.15)); fixed-alpha model tests alpha=-0.055 directly.",
        },
        "h0pe2025_informed": {
            "mu_prior_loc": 8.0,
            "mu_prior_scale": 50.0,
            "tau_prior_scale": 25.0,
            "alpha_prior_loc": 0.0,  # GR-centered null prior for free-alpha model
            "alpha_prior_scale": 0.20,
            "note": (
                "Sensitivity prior set informed by reported H0pe lens-model bias direction "
                "(2510.07637): allows broader/positive residual-bias nuisance without forcing it."
                " Free-alpha prior remains GR-centered to avoid circularity."
            ),
        },
    }

    scenario_results = {}
    for name, cfg in scenarios.items():
        res = evaluate_scenario(
            deltas,
            sigmas,
            r_tep,
            r_unit,
            mu_grid,
            tau_grid,
            alpha_grid,
            cfg["mu_prior_loc"],
            cfg["mu_prior_scale"],
            cfg["tau_prior_scale"],
            cfg["alpha_prior_loc"],
            cfg["alpha_prior_scale"],
        )
        scenario_results[name] = res

        print_status(f"[{name}] logZ(GR)        = {res['logZ_gr']:.2f}")
        print_status(f"[{name}] logZ(TEP fixed) = {res['logZ_tep_fixed']:.2f}")
        print_status(f"[{name}] logZ(TEP free)  = {res['logZ_tep_free']:.2f}")
        print_status(
            f"[{name}] log BF fixed/GR = {res['log_bf_tep_fixed_over_gr']:+.2f} "
            f"(BF={res['bf_tep_fixed_over_gr']:.3g})"
        )
        print_status(
            f"[{name}] log BF free/GR  = {res['log_bf_tep_free_over_gr']:+.2f} "
            f"(BF={res['bf_tep_free_over_gr']:.3g})"
        )
        print_status(
            f"[{name}] Posterior tau   = {res['tau_mean_days']:.2f} "
            f"[{res['tau_p16_days']:.2f}, {res['tau_p84_days']:.2f}] d"
        )
        print_status(
            f"[{name}] Posterior alpha = {res['alpha_mean']:.4f} "
            f"[{res['alpha_p16']:.4f}, {res['alpha_p84']:.4f}]"
        )

    baseline = scenario_results["baseline"]

    # Plot
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from scripts.utils.plot_style import set_pub_style, COLORS, FIG_SIZES, save_fig

        set_pub_style()
        fig, axs = plt.subplots(1, 2, figsize=FIG_SIZES["web_two_panel"])

        axs[0].plot(tau_grid, baseline["posterior_tau_grid"], lw=1.8, color=COLORS["gr"])
        axs[0].axvline(baseline["tau_mean_days"], ls="--", lw=1.0, color=COLORS["observed"])
        axs[0].set_xlabel("tau (extra dispersion, days)")
        axs[0].set_ylabel("Posterior density (arb)")
        axs[0].set_title("Posterior for hierarchical overdispersion")

        axs[1].plot(alpha_grid, baseline["posterior_alpha_grid"], lw=1.8, color=COLORS["tep"])
        axs[1].axvline(-0.055, ls="--", lw=1.0, color=COLORS["observed"], label="Empirical alpha = -0.055")
        axs[1].axvline(0.0, ls=":", lw=1.2, color=COLORS["red"], label="GR null (alpha = 0)")
        axs[1].set_xlabel("alpha")
        axs[1].set_ylabel("Posterior density (arb)")
        axs[1].set_title("Posterior alpha (TEP free model)")
        axs[1].legend()

        out_fig = PROJECT_ROOT / "results" / "figures" / f"step_{STEP_NUM}_bayes_comparison.png"
        save_fig(fig, out_fig)
        print_status(f"Figure saved to {out_fig}")
    except Exception as e:
        print_status(f"Plotting failed: {e}", "ERROR")
        out_fig = None

    out = {
        "step": STEP_NUM,
        "status": "success",
        "priors": {
            "baseline": {
                "mu_bias": "Normal(0, 40 d)",
                "tau": "HalfNormal(20 d)",
                "alpha": "Normal(0, 0.15) [TEP-free model; GR-centered null prior]",
            },
            "h0pe2025_informed": {
                "mu_bias": "Normal(+8, 50 d)",
                "tau": "HalfNormal(25 d)",
                "alpha": "Normal(0, 0.20) [TEP-free model; GR-centered null prior]",
            },
        },
        "scenarios": {
            k: {
                "note": scenarios[k]["note"],
                "evidence": {
                    "logZ_gr": scenario_results[k]["logZ_gr"],
                    "logZ_tep_fixed": scenario_results[k]["logZ_tep_fixed"],
                    "logZ_tep_free": scenario_results[k]["logZ_tep_free"],
                    "log_bf_tep_fixed_over_gr": scenario_results[k]["log_bf_tep_fixed_over_gr"],
                    "log_bf_tep_free_over_gr": scenario_results[k]["log_bf_tep_free_over_gr"],
                    "bf_tep_fixed_over_gr": scenario_results[k]["bf_tep_fixed_over_gr"],
                    "bf_tep_free_over_gr": scenario_results[k]["bf_tep_free_over_gr"],
                },
                "posterior": {
                    "tau_mean_days": scenario_results[k]["tau_mean_days"],
                    "tau_p16_days": scenario_results[k]["tau_p16_days"],
                    "tau_p84_days": scenario_results[k]["tau_p84_days"],
                    "alpha_mean": scenario_results[k]["alpha_mean"],
                    "alpha_p16": scenario_results[k]["alpha_p16"],
                    "alpha_p84": scenario_results[k]["alpha_p84"],
                },
            }
            for k in scenario_results
        },
        "figure": str(out_fig) if out_fig else None,
    }

    out_path = PROJECT_ROOT / "results" / "outputs" / f"step_{STEP_NUM}_bayes_model_comparison.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=safe_json_default)

    print_status(f"Results saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
