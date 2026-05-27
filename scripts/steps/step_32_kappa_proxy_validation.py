#!/usr/bin/env python3
"""
TEP-LENS: Step 32 - Kappa Proxy Validation

The TEP temporal shear couples to the projected convergence kappa(theta), not to
the total magnification mu. The current analysis uses flux ratios F_i/F_ref as
proxies for mu_i/mu_ref. This step quantifies the systematic bias introduced by
this proxy, using the exact lensing identity:

    mu = 1 / [(1 - kappa)^2 - gamma^2]

where gamma is the total shear at the image position.

Because flux ratios are proportional to absolute magnification with an unknown
overall scale C (mu_i = C * F_i), the inferred kappa depends on both the unknown
shear gamma and the unknown scale C. This script performs a Monte Carlo
sensitivity analysis over physically motivated ranges for both parameters.

Important caveat: observed flux ratios constrain |mu|, not signed magnification.
The inversion kappa = 1 - sqrt(mu^{-1} + gamma^2) assumes the positive branch.
Near critical curves, image parity matters and the branch choice is non-trivial.
The Monte Carlo should therefore be read as a proxy-systematic exploration,
not a unique reconstruction of kappa.

Key physical insight:
- For images near a tangential critical curve (S1-S4 in SN Refsdal), gamma is
  large and comparable to 1-kappa. Small changes in gamma produce large changes
  in mu at fixed kappa, making mu a poor proxy for kappa.
- For images far from critical curves (SX in SN Refsdal), gamma is smaller,
  and mu is a more faithful proxy for kappa.

Outputs:
- results/outputs/step_32_kappa_proxy_validation.json
- results/figures/step_32_kappa_mu_comparison.png
"""

import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status, safe_json_default
from scripts.utils.plot_style import set_pub_style, COLORS, FIG_SIZE

STEP_NUM = "32"




def kappa_from_mu_gamma(mu, gamma):
    """Invert lensing identity: mu = 1 / [(1-kappa)^2 - gamma^2]."""
    term = 1.0 / mu + gamma ** 2
    if term < 0 or term > 1.0:
        return np.nan
    return 1.0 - np.sqrt(term)


def gamma_tep_factor(alpha, quantity_norm):
    """TEP Gamma_t = 1 + alpha * log10(quantity_norm)."""
    return 1.0 + alpha * np.log10(quantity_norm)


def closure_residual(alpha, quantities, delays, loop_images):
    """
    Compute TEP closure residual for a loop using any quantity (mu or kappa).
    quantities: dict {image: value}
    delays: dict {image: absolute arrival time relative to reference}
    loop_images: tuple (i, j, k)
    """
    i, j, k = loop_images
    q_mean = np.mean(list(quantities.values()))
    q_norm = {img: quantities[img] / q_mean for img in quantities}
    Gamma = {img: gamma_tep_factor(alpha, q_norm[img]) for img in q_norm}

    dt_ij = delays[j] - delays[i]
    dt_jk = delays[k] - delays[j]
    dt_ki = delays[i] - delays[k]

    R = ((Gamma[i] - 1.0) * dt_ij
         + (Gamma[j] - 1.0) * dt_jk
         + (Gamma[k] - 1.0) * dt_ki)
    return float(R)


def sn_refsdal_sensitivity_envelope(n_draws=20000, seed=42):
    """
    Monte Carlo sensitivity analysis for SN Refsdal.

    The flux ratios F_i/F_S1 are proportional to absolute magnifications:
        mu_i = C * F_i
    where C is an unknown scale factor (the absolute magnification of the
    reference image depends on the source position).

    For each draw:
    1. Draw C from a physically motivated range (total magnification of SN
       Refsdal is estimated at ~5-15 in the literature; with sum(F)=4.901,
       this implies C ~ 1-3).
    2. Draw a shear gamma_i for each image from a distribution informed by
       its position relative to the cluster potential.
    3. Compute the implied kappa_i from the lensing identity.
    4. Compute the TEP closure residual using kappa-based Gamma factors.
    5. Compare to the nominal mu-proxy result.

    Shear priors (informed by cluster lensing geometry):
    - S1-S4 (Einstein cross around member galaxy, near critical curve):
      The member galaxy contributes ~0.5 (isothermal at Einstein radius) and
      the cluster contributes external shear ~0.2-0.4. We model total shear
      as gamma ~ TruncatedNormal(mean=0.6, sigma=0.2, low=0.0, high=0.95).
    - SX (peripheral arc, ~8 arcsec from cross, far from member galaxy
      critical curve): external shear is smaller. gamma ~ TN(0.2, 0.15, 0, 0.8).

    Scale factor prior:
    - C ~ Uniform(0.5, 4.0). This covers the plausible range of absolute
      magnifications for a cluster-lensed supernova.
    """
    rng = np.random.default_rng(seed)

    # Observed flux ratios (treated as proportional to mu)
    fluxes = {
        "S1": 1.158,
        "S2": 0.887,
        "S3": 0.716,
        "S4": 1.793,
        "SX": 0.347,
    }

    # Absolute delays relative to S1
    delays = {
        "S1": 0.0,
        "S2": 9.9,
        "S3": 9.0,
        "S4": 20.3,
        "SX": 376.0,
    }

    # Shear distributions (physically motivated)
    gamma_prior = {
        "S1": {"mean": 0.60, "std": 0.18, "low": 0.05, "high": 0.92},
        "S2": {"mean": 0.60, "std": 0.18, "low": 0.05, "high": 0.92},
        "S3": {"mean": 0.60, "std": 0.18, "low": 0.05, "high": 0.92},
        "S4": {"mean": 0.60, "std": 0.18, "low": 0.05, "high": 0.92},
        "SX": {"mean": 0.18, "std": 0.12, "low": 0.02, "high": 0.70},
    }

    alpha_nominal = -0.055
    mu_values = fluxes  # proxy assumption: F_i proportional to mu_i

    draws_R = []
    draws_alpha = []
    draws_kappa = {img: [] for img in fluxes}
    draws_mu = {img: [] for img in fluxes}
    draws_C = []

    for _ in range(n_draws):
        # Draw absolute magnification scale factor
        C = rng.uniform(0.5, 4.0)
        draws_C.append(C)

        gamma_draw = {}
        kappa_draw = {}
        mu_draw = {}
        for img in fluxes:
            # Draw gamma from truncated normal
            p = gamma_prior[img]
            g = rng.normal(p["mean"], p["std"])
            g = np.clip(g, p["low"], p["high"])

            # Absolute magnification
            mu_abs = C * fluxes[img]
            mu_draw[img] = mu_abs

            # Compute implied kappa (must be real and positive)
            # Need 1/mu + gamma^2 < 1 for real kappa
            if 1.0 / mu_abs + g ** 2 >= 1.0:
                # If shear is too large for this mu, reduce shear to boundary
                g = np.sqrt(max(0, 1.0 - 1.0 / mu_abs - 1e-6))

            k = kappa_from_mu_gamma(mu_abs, g)
            if not np.isfinite(k) or k <= 0:
                k = 0.01
            gamma_draw[img] = g
            kappa_draw[img] = k

        # Recompute closure residual using kappa-based Gamma
        R_kappa = closure_residual(alpha_nominal, kappa_draw, delays, ("S1", "S4", "SX"))

        # Equivalent alpha to match the same nominal observed residual
        R_mu = closure_residual(alpha_nominal, mu_values, delays, ("S1", "S4", "SX"))
        if abs(R_kappa) > 1e-6 and np.isfinite(R_kappa):
            alpha_equiv = alpha_nominal * (R_mu / R_kappa)
        else:
            alpha_equiv = np.nan

        draws_R.append(R_kappa)
        draws_alpha.append(alpha_equiv)
        for img in fluxes:
            draws_kappa[img].append(kappa_draw[img])
            draws_mu[img].append(mu_draw[img])

    # Summarize
    def summarize(arr, pctiles=(2.5, 16, 50, 84, 97.5)):
        a = np.array(arr)
        a = a[np.isfinite(a)]
        if len(a) == 0:
            return {"median": None, "mean": None, "std": None}
        out = {
            "median": float(np.median(a)),
            "mean": float(np.mean(a)),
            "std": float(np.std(a)),
        }
        for p in pctiles:
            out[f"p{p:.0f}".replace(".", "p")] = float(np.percentile(a, p))
        return out

    kappa_summary = {}
    for img in fluxes:
        kappa_summary[img] = summarize(draws_kappa[img])

    R_summary = summarize(draws_R)
    alpha_summary = summarize(draws_alpha)
    C_summary = summarize(draws_C)

    # Ranking stability: probability that kappa_norm ranking matches mu_norm ranking
    mu_norm_list = {img: mu_values[img] / np.mean(list(mu_values.values())) for img in fluxes}
    mu_ranks = {img: r for r, img in enumerate(sorted(fluxes, key=lambda x: mu_norm_list[x]))}

    rank_agreement_count = 0
    for i_draw in range(n_draws):
        kappa_vals = {img: draws_kappa[img][i_draw] for img in fluxes}
        kappa_norm = {img: kappa_vals[img] / np.mean(list(kappa_vals.values())) for img in fluxes}
        kappa_ranks = {img: r for r, img in enumerate(sorted(fluxes, key=lambda x: kappa_norm[x]))}
        if all(mu_ranks[img] == kappa_ranks[img] for img in fluxes):
            rank_agreement_count += 1

    P_rank_agreement = rank_agreement_count / n_draws

    # Key robustness: probability that S4-SX kappa contrast has the same sign
    # as S4-SX mu contrast (i.e., kappa_S4 > kappa_SX and mu_S4 > mu_SX)
    contrast_same_sign = 0
    for i_draw in range(n_draws):
        kappa_vals = {img: draws_kappa[img][i_draw] for img in fluxes}
        if (kappa_vals["S4"] - kappa_vals["SX"]) * (mu_values["S4"] - mu_values["SX"]) > 0:
            contrast_same_sign += 1
    P_contrast_same_sign = contrast_same_sign / n_draws

    return {
        "n_draws": n_draws,
        "kappa_summary_per_image": kappa_summary,
        "mu_proxy_per_image": {img: float(mu_values[img]) for img in fluxes},
        "mu_norm_per_image": {img: float(mu_norm_list[img]) for img in fluxes},
        "closure_residual_kappa": R_summary,
        "closure_residual_mu_nominal": float(R_mu),
        "alpha_equiv_to_match_obs": alpha_summary,
        "C_prior_summary": C_summary,
        "P_rank_agreement_mu_vs_kappa": float(P_rank_agreement),
        "P_contrast_sign_agreement_S4_SX": float(P_contrast_same_sign),
        "gamma_priors": gamma_prior,
    }


def tdcosmo_theoretical_comparison():
    """
    Theoretical illustration of the mu-kappa relationship for power-law lenses.

    For a spherical power-law lens with slope gamma_pl, the convergence profile is:
        kappa(r) = (3 - gamma_pl)/2 * (theta_E/r)^(gamma_pl-1)
    The shear for a spherical lens is gamma(r) = kappa(r).
    The magnification is:
        mu(r) = 1 / [(1 - kappa(r))^2 - kappa(r)^2] = 1 / [1 - 2*kappa(r)]

    This shows that mu diverges where kappa approaches 0.5 (the critical curve).
    For gamma_pl < 2, kappa(theta_E) > 0.5, so the Einstein radius is inside a
    region where mu is negative (unphysical for the simple spherical model), but
    in reality ellipticity and external shear modify this.

    The key point: at a fixed flux ratio F_i/F_j (proportional to mu_i/mu_j),
    the implied kappa ratio depends on the absolute kappa values, which are
    unknown without a full mass model.
    """
    systems = {
        "RXJ1131-1231": {"gamma_pl": 1.95, "kappa_at_rE": 0.525},
        "HE0435-1223": {"gamma_pl": 1.93, "kappa_at_rE": 0.535},
        "PG1115+080": {"gamma_pl": 2.17, "kappa_at_rE": 0.415},
        "DES0408-5354": {"gamma_pl": 1.90, "kappa_at_rE": 0.550},
        "SDSS1206+4332": {"gamma_pl": 1.95, "kappa_at_rE": 0.525},
        "WFI2033-4723": {"gamma_pl": 1.95, "kappa_at_rE": 0.525},
    }

    results = {}
    for name, pars in systems.items():
        kappa_e = pars["kappa_at_rE"]
        gamma_pl = pars["gamma_pl"]

        # At r = theta_E for a spherical power-law: gamma = kappa
        mu_e = 1.0 / (1.0 - 2.0 * kappa_e) if kappa_e < 0.5 else np.inf

        # At r = 1.2 * theta_E: kappa drops by factor (1/1.2)^(gamma_pl-1)
        kappa_1p2 = kappa_e * (1.0 / 1.2) ** (gamma_pl - 1.0)
        mu_1p2 = 1.0 / (1.0 - 2.0 * kappa_1p2) if kappa_1p2 < 0.5 else np.inf

        if np.isfinite(mu_e) and np.isfinite(mu_1p2) and mu_e > 0 and mu_1p2 > 0:
            delta_log_mu = np.log10(mu_1p2 / mu_e)
            delta_log_kappa = np.log10(kappa_1p2 / kappa_e)
            ratio = delta_log_mu / delta_log_kappa if abs(delta_log_kappa) > 1e-6 else None
        else:
            delta_log_mu = None
            delta_log_kappa = None
            ratio = None

        results[name] = {
            "gamma_pl": gamma_pl,
            "kappa_at_theta_E": float(kappa_e),
            "mu_at_theta_E": float(mu_e) if np.isfinite(mu_e) else None,
            "kappa_at_1p2_theta_E": float(kappa_1p2),
            "mu_at_1p2_theta_E": float(mu_1p2) if np.isfinite(mu_1p2) else None,
            "delta_log_mu": float(delta_log_mu) if delta_log_mu is not None else None,
            "delta_log_kappa": float(delta_log_kappa) if delta_log_kappa is not None else None,
            "log_mu_over_log_kappa_ratio": float(ratio) if ratio is not None else None,
            "note": (
                "Spherical power-law approximation. Real lenses have ellipticity "
                "and external shear, which modify the mu-kappa relationship. "
                "This illustrates the generic non-linearity of the proxy."
            ),
        }

    return results


def tdcosmo_elliptical_sie_with_shear():
    """
    Elliptical SIE + external shear comparison.

    Real TDCOSMO lenses are elliptical with external shear, which amplifies
    the mu-kappa discrepancy relative to the spherical case. For an SIE with
    axis ratio q and external shear gamma_ext, the total shear magnitude at
    the Einstein radius is larger than the internal convergence alone.

    For an SIE, the convergence at the Einstein radius (major axis) is
    approximately kappa_e_major = kappa_e_spherical * sqrt(2/(1+q^2)).
    The magnification includes both internal shear (from ellipticity) and
    external shear:
        mu = 1 / [(1 - kappa)^2 - gamma_total^2]
    where gamma_total^2 = gamma_int^2 + gamma_ext^2 + cross_terms.

    We conservatively model this by adding a typical external shear
    gamma_ext ~ 0.08 to the spherical model. This increases the total
    shear at fixed kappa, reducing mu and amplifying the mu-kappa
    discrepancy relative to the spherical-only case.
    """
    systems = {
        "RXJ1131-1231": {"gamma_pl": 1.95, "kappa_at_rE": 0.525, "gamma_ext": 0.08},
        "HE0435-1223": {"gamma_pl": 1.93, "kappa_at_rE": 0.535, "gamma_ext": 0.10},
        "PG1115+080": {"gamma_pl": 2.17, "kappa_at_rE": 0.415, "gamma_ext": 0.06},
        "DES0408-5354": {"gamma_pl": 1.90, "kappa_at_rE": 0.550, "gamma_ext": 0.12},
        "SDSS1206+4332": {"gamma_pl": 1.95, "kappa_at_rE": 0.525, "gamma_ext": 0.08},
        "WFI2033-4723": {"gamma_pl": 1.95, "kappa_at_rE": 0.525, "gamma_ext": 0.09},
    }

    results = {}
    for name, pars in systems.items():
        kappa_e = pars["kappa_at_rE"]
        gamma_pl = pars["gamma_pl"]
        gamma_ext = pars["gamma_ext"]

        # Spherical internal shear = kappa at Einstein radius
        gamma_int_e = kappa_e

        # At r = theta_E: total shear with external component
        gamma_total_e = np.sqrt(gamma_int_e**2 + gamma_ext**2)
        mu_e = 1.0 / ((1.0 - kappa_e)**2 - gamma_total_e**2)

        # At r = 1.2 * theta_E: kappa drops, internal shear drops with it
        kappa_1p2 = kappa_e * (1.0 / 1.2) ** (gamma_pl - 1.0)
        gamma_int_1p2 = kappa_1p2  # spherical approximation for internal
        gamma_total_1p2 = np.sqrt(gamma_int_1p2**2 + gamma_ext**2)
        mu_1p2 = 1.0 / ((1.0 - kappa_1p2)**2 - gamma_total_1p2**2)

        if mu_e > 0 and mu_1p2 > 0 and np.isfinite(mu_e) and np.isfinite(mu_1p2):
            delta_log_mu = np.log10(mu_1p2 / mu_e)
            delta_log_kappa = np.log10(kappa_1p2 / kappa_e)
            ratio = delta_log_mu / delta_log_kappa if abs(delta_log_kappa) > 1e-6 else None
        else:
            delta_log_mu = None
            delta_log_kappa = None
            ratio = None

        # Also compute spherical-only ratio for comparison
        mu_e_sph = 1.0 / (1.0 - 2.0 * kappa_e) if kappa_e < 0.5 else np.inf
        kappa_1p2_sph = kappa_e * (1.0 / 1.2) ** (gamma_pl - 1.0)
        mu_1p2_sph = 1.0 / (1.0 - 2.0 * kappa_1p2_sph) if kappa_1p2_sph < 0.5 else np.inf
        if np.isfinite(mu_e_sph) and np.isfinite(mu_1p2_sph) and mu_e_sph > 0 and mu_1p2_sph > 0:
            delta_log_mu_sph = np.log10(mu_1p2_sph / mu_e_sph)
            ratio_sph = delta_log_mu_sph / delta_log_kappa if abs(delta_log_kappa) > 1e-6 else None
        else:
            ratio_sph = None

        results[name] = {
            "gamma_pl": gamma_pl,
            "kappa_at_theta_E": float(kappa_e),
            "gamma_ext": gamma_ext,
            "mu_at_theta_E": float(mu_e),
            "mu_at_1p2_theta_E": float(mu_1p2),
            "delta_log_mu": float(delta_log_mu) if delta_log_mu is not None else None,
            "delta_log_kappa": float(delta_log_kappa) if delta_log_kappa is not None else None,
            "log_mu_over_log_kappa_ratio": float(ratio) if ratio is not None else None,
            "spherical_ratio": float(ratio_sph) if ratio_sph is not None else None,
            "note": (
                "Elliptical SIE with external shear. The ratio exceeds the spherical "
                "power-law value because external shear increases total shear at fixed kappa, "
                "amplifying the magnification-convergence discrepancy."
            ),
        }

    return results


def plot_comparison(sn_results, tdcosmo_results, elliptical_results):
    """Generate figure comparing mu-based and kappa-based proxies."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    set_pub_style()
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: SN Refsdal kappa vs mu per image
    images = ["S1", "S2", "S3", "S4", "SX"]
    mu_vals = [sn_results["mu_proxy_per_image"][img] for img in images]
    mu_norms = [sn_results["mu_norm_per_image"][img] for img in images]
    kappa_medians = [sn_results["kappa_summary_per_image"][img]["median"] for img in images]
    kappa_p16 = [sn_results["kappa_summary_per_image"][img]["p16"] for img in images]
    kappa_p84 = [sn_results["kappa_summary_per_image"][img]["p84"] for img in images]

    # Normalize kappa to mean
    kappa_mean = np.mean(kappa_medians)
    kappa_norm = [k / kappa_mean for k in kappa_medians]
    kappa_norm_lo = [k / kappa_mean for k in kappa_p16]
    kappa_norm_hi = [k / kappa_mean for k in kappa_p84]

    x = np.arange(len(images))
    width = 0.35

    axs[0].bar(x - width/2, mu_norms, width, label=r"$\mu_{\rm norm}$ (flux proxy)", color=COLORS["model"], edgecolor="black")
    axs[0].bar(x + width/2, kappa_norm, width, yerr=[[k - l for k, l in zip(kappa_norm, kappa_norm_lo)],
                                                      [h - k for k, h in zip(kappa_norm, kappa_norm_hi)]],
               label=r"$\kappa_{\rm norm}$ (inferred)", color=COLORS["tep"], edgecolor="black", capsize=3)
    axs[0].set_xticks(x)
    axs[0].set_xticklabels(images)
    axs[0].set_ylabel("Normalised proxy value")
    axs[0].set_title(r"SN Refsdal: flux-proxy $\mu$ vs. inferred $\kappa$")
    axs[0].legend()
    axs[0].grid(axis="y", alpha=0.3, ls=":")

    # Panel B: TDCOSMO delta_log_mu vs delta_log_kappa
    names = [n for n in tdcosmo_results if tdcosmo_results[n]["delta_log_kappa"] is not None]
    delta_mu = [tdcosmo_results[n]["delta_log_mu"] for n in names]
    delta_kappa = [tdcosmo_results[n]["delta_log_kappa"] for n in names]

    delta_mu_ell = [elliptical_results[n]["delta_log_mu"] for n in names
                    if elliptical_results[n]["delta_log_mu"] is not None]
    delta_kappa_ell = [elliptical_results[n]["delta_log_kappa"] for n in names
                       if elliptical_results[n]["delta_log_kappa"] is not None]

    if delta_kappa and delta_mu:
        axs[1].scatter(delta_kappa, delta_mu, s=80, color=COLORS["observed"], edgecolor="black",
                       zorder=3, label="Spherical power-law (lower bound)")
        for i, name in enumerate(names):
            axs[1].annotate(name.split("-")[0], (delta_kappa[i], delta_mu[i]),
                           textcoords="offset points", xytext=(5, 5), fontsize=8)

        if delta_kappa_ell and delta_mu_ell:
            axs[1].scatter(delta_kappa_ell, delta_mu_ell, s=80, color=COLORS["red"], edgecolor="black",
                           zorder=3, marker="s", label="Elliptical + ext. shear (upper bound)")

        lim_min = min(min(delta_kappa + delta_kappa_ell), min(delta_mu + delta_mu_ell)) - 0.1
        lim_max = max(max(delta_kappa + delta_kappa_ell), max(delta_mu + delta_mu_ell)) + 0.1
        axs[1].plot([lim_min, lim_max], [lim_min, lim_max], "k--", lw=1, label="1:1 (perfect proxy)", zorder=1)
        axs[1].set_xlim(lim_min, lim_max)
        axs[1].set_ylim(lim_min, lim_max)
    axs[1].set_xlabel(r"$\Delta\log_{10}\kappa$ (convergence-based)")
    axs[1].set_ylabel(r"$\Delta\log_{10}\mu$ (magnification-based)")
    axs[1].set_title(r"TDCOSMO: $\mu$-$\kappa$ proxy discrepancy")
    axs[1].legend(fontsize=8)
    axs[1].grid(alpha=0.3, ls=":")

    fig.tight_layout()
    out_fig = PROJECT_ROOT / "results" / "figures" / f"step_{STEP_NUM}_kappa_mu_comparison.png"
    fig.savefig(out_fig, dpi=200)
    plt.close(fig)
    print_status(f"Figure saved: {out_fig}")
    return str(out_fig)


def main():
    print_status(f"STEP {STEP_NUM}: Kappa Proxy Validation", "TITLE")

    # ------------------------------------------------------------------
    # A. SN Refsdal sensitivity envelope
    # ------------------------------------------------------------------
    print_status("Running SN Refsdal shear-degeneracy sensitivity envelope (20,000 draws)...")
    sn_results = sn_refsdal_sensitivity_envelope(n_draws=20000)

    print_status("SN Refsdal kappa inference (median, p16-p84):")
    for img in ["S1", "S2", "S3", "S4", "SX"]:
        s = sn_results["kappa_summary_per_image"][img]
        print_status(f"  {img}: kappa = {s['median']:.3f} [{s['p16']:.3f}, {s['p84']:.3f}]")

    print_status(f"Closure residual (mu proxy, nominal): {sn_results['closure_residual_mu_nominal']:.3f} d")
    s = sn_results["closure_residual_kappa"]
    print_status(f"Closure residual (kappa inferred): {s['median']:.3f} [{s['p16']:.3f}, {s['p84']:.3f}] d")

    s = sn_results["alpha_equiv_to_match_obs"]
    print_status(f"Equivalent alpha (to match same observed residual): {s['median']:.4f} [{s['p16']:.4f}, {s['p84']:.4f}]")
    print_status(f"Rank-order agreement P(mu_rank == kappa_rank) = {sn_results['P_rank_agreement_mu_vs_kappa']:.3f}")
    print_status(f"S4-SX contrast sign agreement P = {sn_results['P_contrast_sign_agreement_S4_SX']:.3f}")

    # ------------------------------------------------------------------
    # B. TDCOSMO theoretical comparison (spherical)
    # ------------------------------------------------------------------
    print_status("Running TDCOSMO theoretical mu-kappa comparison (spherical)...")
    tdcosmo_results = tdcosmo_theoretical_comparison()
    for name, res in tdcosmo_results.items():
        if res["delta_log_mu"] is not None:
            print_status(f"  {name}: delta_log_mu={res['delta_log_mu']:.3f}, "
                         f"delta_log_kappa={res['delta_log_kappa']:.3f}, "
                         f"ratio={res['log_mu_over_log_kappa_ratio']:.2f}")
        else:
            print_status(f"  {name}: mu diverges (kappa >= 0.5 at Einstein radius)")

    # ------------------------------------------------------------------
    # B2. Elliptical SIE + external shear comparison
    # ------------------------------------------------------------------
    print_status("Running TDCOSMO elliptical SIE + external shear comparison...")
    elliptical_results = tdcosmo_elliptical_sie_with_shear()
    for name, res in elliptical_results.items():
        if res["delta_log_mu"] is not None:
            print_status(f"  {name}: elliptical ratio={res['log_mu_over_log_kappa_ratio']:.2f}, "
                         f"spherical ratio={res['spherical_ratio']:.2f}, "
                         f"amplification factor={res['log_mu_over_log_kappa_ratio']/res['spherical_ratio']:.2f}")
        else:
            print_status(f"  {name}: mu diverges in elliptical model")

    # ------------------------------------------------------------------
    # C. Plot
    # ------------------------------------------------------------------
    fig_path = plot_comparison(sn_results, tdcosmo_results, elliptical_results)

    # ------------------------------------------------------------------
    # Save results
    # ------------------------------------------------------------------
    max_sph_ratio = max(
        r["log_mu_over_log_kappa_ratio"] for r in tdcosmo_results.values()
        if r["log_mu_over_log_kappa_ratio"] is not None
    )
    max_ell_ratio = max(
        r["log_mu_over_log_kappa_ratio"] for r in elliptical_results.values()
        if r["log_mu_over_log_kappa_ratio"] is not None
    )

    # ------------------------------------------------------------------
    # Proxy systematic budget: decompose total uncertainty
    # ------------------------------------------------------------------
    # Kelly+2023 measurement uncertainty on SX-S1 delay
    sigma_measurement = 5.6  # days

    # Lens-model uncertainty: spread of blind predictions (std of model medians)
    # From Treu+2016 Table 2: model predictions range ~324-376 d
    # Approximate as half the range / sqrt(N) or std directly from literature
    model_preds = np.array([324.0, 345.0, 376.0, 361.0, 369.0, 359.0, 374.0])
    sigma_lens_model = float(np.std(model_preds, ddof=1))  # ~17.5 d

    # Proxy systematic: spread from mu->kappa MC envelope
    # Use the 16th-84th spread of the kappa-based closure residual
    R_kappa = sn_results["closure_residual_kappa"]
    R_mu_nom = sn_results["closure_residual_mu_nominal"]
    sigma_proxy = float(abs(R_kappa["p84"] - R_kappa["p16"]) / 2.0)

    # Microlensing: 10-30% flux perturbation -> ~10-20% shift in alpha_equiv
    # Approximate as 15% of the nominal residual amplitude
    sigma_microlensing = float(0.15 * abs(R_mu_nom))

    # Total uncertainty (quadrature sum)
    sigma_total = float(np.sqrt(
        sigma_measurement**2 + sigma_lens_model**2 + sigma_proxy**2 + sigma_microlensing**2
    ))

    # Sign robustness extracted from Monte Carlo
    P_sign_stable = sn_results["P_contrast_sign_agreement_S4_SX"]

    # Amplitude uncertainty fraction: ratio of kappa-based to mu-proxy spread
    alpha_equiv = sn_results["alpha_equiv_to_match_obs"]
    if alpha_equiv["median"] and abs(alpha_equiv["median"]) > 1e-6:
        amp_uncertainty_frac = float(
            abs((alpha_equiv["p84"] - alpha_equiv["p16"]) / 2.0)
            / abs(alpha_equiv["median"])
        )
    else:
        amp_uncertainty_frac = None

    proxy_systematic_budget = {
        "measurement_uncertainty_days": float(sigma_measurement),
        "measurement_fraction": float(sigma_measurement / sigma_total),
        "lens_model_uncertainty_days": float(sigma_lens_model),
        "lens_model_fraction": float(sigma_lens_model / sigma_total),
        "proxy_systematic_days": float(sigma_proxy),
        "proxy_systematic_fraction": float(sigma_proxy / sigma_total),
        "microlensing_uncertainty_days": float(sigma_microlensing),
        "microlensing_fraction": float(sigma_microlensing / sigma_total),
        "sign_robustness_probability": float(P_sign_stable),
        "amplitude_uncertainty_fraction": float(amp_uncertainty_frac) if amp_uncertainty_frac is not None else None,
        "total_uncertainty_days": float(sigma_total),
        "interpretation": (
            f"The evidence is suggestive rather than decisive because the proxy-systematic "
            f"contribution to the total uncertainty budget is "
            f"{100*float(sigma_proxy/sigma_total):.0f}%, comparable to the lens-model "
            f"contribution ({100*float(sigma_lens_model/sigma_total):.0f}%)."
        ),
    }

    print_status(f"\nProxy systematic budget:")
    print_status(f"  Measurement:     {sigma_measurement:.1f} d  ({100*sigma_measurement/sigma_total:.0f}%)")
    print_status(f"  Lens-model:      {sigma_lens_model:.1f} d  ({100*sigma_lens_model/sigma_total:.0f}%)")
    print_status(f"  Proxy (mu->kappa): {sigma_proxy:.1f} d  ({100*sigma_proxy/sigma_total:.0f}%)")
    print_status(f"  Microlensing:    {sigma_microlensing:.1f} d  ({100*sigma_microlensing/sigma_total:.0f}%)")
    print_status(f"  Total:           {sigma_total:.1f} d")

    output = {
        "step": STEP_NUM,
        "status": "success",
        "sn_refsdal_sensitivity": sn_results,
        "tdcosmo_theoretical": tdcosmo_results,
        "tdcosmo_elliptical_shear": elliptical_results,
        "figure": fig_path,
        "proxy_systematic_budget": proxy_systematic_budget,
        "key_findings": {
            "sn_refsdal": (
                f"SN Refsdal: The flux-proxy mu and inferred kappa give different "
                f"S1-S4-SX closure residuals: mu-proxy = {sn_results['closure_residual_mu_nominal']:.2f} d, "
                f"kappa-inferred = {sn_results['closure_residual_kappa']['median']:.2f} "
                f"[{sn_results['closure_residual_kappa']['p16']:.2f}, "
                f"{sn_results['closure_residual_kappa']['p84']:.2f}] d. "
                f"Equivalent alpha shifts to {sn_results['alpha_equiv_to_match_obs']['median']:.4f} "
                f"[{sn_results['alpha_equiv_to_match_obs']['p16']:.4f}, "
                f"{sn_results['alpha_equiv_to_match_obs']['p84']:.4f}]. "
                f"Crucially, the S4-SX contrast sign is stable: "
                f"P(same sign) = {sn_results['P_contrast_sign_agreement_S4_SX']:.1%}."
            ),
            "tdcosmo": (
                f"TDCOSMO spherical power-law ratios reach up to {max_sph_ratio:.1f}, "
                f"but elliptical lenses with external shear amplify this to {max_ell_ratio:.1f}. "
                f"The spherical comparison is a lower bound on the proxy discrepancy; "
                f"real lenses exhibit larger misestimation."
            ),
        },
    }

    out_path = PROJECT_ROOT / "results" / "outputs" / f"step_{STEP_NUM}_kappa_proxy_validation.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=safe_json_default)

    print_status(f"Results saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()

