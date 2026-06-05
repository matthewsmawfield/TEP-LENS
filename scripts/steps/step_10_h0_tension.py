#!/usr/bin/env python3
"""
TEP-LENS: Step 10 - H0 Tension under measured proxy-model coupling (alpha_proxy ≈ -0.055)

Computes the H0 shift for each lensed SN system from first principles using
proxy-model gamma factors derived from the magnification data in the step_01 catalog.

Key physics:
- Under the proxy model, arrival times scale as t_i' = Gamma_i * t_i where
  Gamma_i = 1 + alpha * log10(mu_i / mu_mean).
- For a delay relative to reference image ref: dt_i,ref = t_i - t_ref.
  Under the proxy model: dt_i,ref^proxy = Gamma_i * t_i - Gamma_ref * t_ref.
  With t_ref = 0: dt_i,ref^proxy = Gamma_i * dt_i,ref^GR.
  Therefore the H0 shift factor for this delay is f_i = Gamma_i.

For systems with 3+ images, the lens model uses multiple delays. The net H0
shift is approximated as the precision-weighted average of the per-delay
shift factors.

For 2-image systems (Encore), the shift depends on reference choice. The shift
is reported using the later-arriving image as reference, with an explicit note
of the reference-dependence caveat.

Results are propagated honestly: systems with small magnification contrasts
(H0pe, Encore) receive small proxy-model shifts, while Refsdal (large S1/SX contrast)
receives the dominant shift.
"""

import json
import sys
from pathlib import Path
import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status
from scripts.utils.tep_config import (
    ALPHA_PROXY,
    H0_REFSDAL_GR,
    H0_REFSDAL_ERR_PLUS,
    H0_REFSDAL_ERR_MINUS,
    H0_PLANCK,
    H0_PLANCK_ERR_PLUS,
    H0_PLANCK_ERR_MINUS,
    H0_H0PE_ENC_GR,
    H0_H0PE_ENC_ERR_PLUS,
    H0_H0PE_ENC_ERR_MINUS,
    SIGMA_ALPHA_PROXY,
)

STEP_NUM = "10"
ALPHA = ALPHA_PROXY


def gamma_factor(alpha, mu, mu_mean):
    """Proxy-model temporal gamma factor for image with magnification mu."""
    return 1.0 + alpha * np.log10(mu / mu_mean)


def gamma_error(alpha, mu, mu_err, mu_mean, n_images):
    """Approximate gamma uncertainty from magnification error."""
    # dGamma/dmu = alpha / (mu * ln(10))
    dgamma_dmu = alpha / (mu * np.log(10))
    # mu_mean uncertainty: sigma_mean ~ sigma_mu / sqrt(n)
    sigma_mean = mu_err / np.sqrt(n_images)
    # Total derivative including mu_mean dependence
    # Gamma = 1 + alpha*(log10(mu) - log10(mu_mean))
    # dGamma = alpha/(ln(10)) * (dmu/mu - dmu_mean/mu_mean)
    dgamma = abs(alpha) / np.log(10) * np.sqrt((mu_err / mu) ** 2 + (sigma_mean / mu_mean) ** 2)
    return dgamma


def compute_h0_shift_refsdal(alpha):
    """Load actual proxy-model prediction from step_07 for SN Refsdal."""
    s07_path = PROJECT_ROOT / "results" / "outputs" / "step_07_observed_vs_predicted.json"
    if not s07_path.exists():
        raise FileNotFoundError(
            f"Required upstream output not found: {s07_path}\n"
            "Run step_07_observed_vs_predicted.py first."
        )
    with open(s07_path) as f:
        s07 = json.load(f)
    r_tep = float(s07["tep_prediction"]["R_tep_prediction_days"])

    # Load observed SX delay from canonical catalog
    catalog_path = PROJECT_ROOT / "data" / "raw" / "sn_lensing" / "lensed_sn_catalog.json"
    with open(catalog_path) as f:
        catalog = json.load(f)
    dt_obs = float(catalog["sn_refsdal"]["time_delays_days"]["dt_SX_S1"]["value"])
    dt_geom = dt_obs - r_tep
    f_tep = dt_obs / dt_geom

    h0_gr = H0_REFSDAL_GR
    err_plus = H0_REFSDAL_ERR_PLUS
    err_minus = H0_REFSDAL_ERR_MINUS
    h0_tep = h0_gr * f_tep

    # Error propagation: sigma_h0 = h0_gr * dt_geom * sigma_r_tep / dt_geom^2
    # Load the actual S1-S4-SX loop uncertainty from step_03 (avoid stale hardcoded value).
    s03_path = PROJECT_ROOT / "results" / "outputs" / "step_03_tep_closure.json"
    if not s03_path.exists():
        raise FileNotFoundError(
            f"Required upstream output not found: {s03_path}\n"
            "Run step_03_tep_closure.py first."
        )
    with open(s03_path) as f:
        s03 = json.load(f)
    sigma_r_tep = float(
        s03["tep_predicted_discrepancies"]["S1_S4_SX"]["tep_gr_discrepancy_err_days"]
    )
    sigma_f = dt_obs * sigma_r_tep / (dt_geom ** 2)
    sigma_h0_tep = h0_gr * sigma_f

    return {
        "h0_gr": h0_gr,
        "h0_tep": h0_tep,
        "err_plus": err_plus,
        "err_minus": err_minus,
        "dt_obs": dt_obs,
        "dt_geom": dt_geom,
        "shift_factor": f_tep,
        "sigma_shift_factor": sigma_f,
        "sigma_h0_tep": sigma_h0_tep,
        "r_tep_days": r_tep,
        "note": (
            "Rigorous closure-based shift. dt_geom = dt_obs - R_tep_prediction. "
            "H0 scales inversely with the true geometric delay."
        ),
    }


def compute_h0_shift_h0pe(alpha):
    """Compute proxy-model H0 shift for SN H0pe from 3-image gamma factors."""
    catalog_path = PROJECT_ROOT / "data" / "raw" / "sn_lensing" / "lensed_sn_catalog.json"
    with open(catalog_path) as f:
        catalog = json.load(f)
    h0pe = catalog["sn_h0pe"]

    mu = {
        "A": h0pe["magnification_proxies"]["mu_absolute"]["A"]["value"],
        "B": h0pe["magnification_proxies"]["mu_absolute"]["B"]["value"],
        "C": h0pe["magnification_proxies"]["mu_absolute"]["C"]["value"],
    }
    def _get_mu_err(img):
        mu_data = h0pe["magnification_proxies"]["mu_absolute"][img]
        if "err" in mu_data:
            return float(mu_data["err"])
        return (float(mu_data["err_plus"]) + float(mu_data["err_minus"])) / 2

    mu_err = {"A": _get_mu_err("A"), "B": _get_mu_err("B"), "C": _get_mu_err("C")}

    mu_mean = np.mean(list(mu.values()))
    gamma = {img: gamma_factor(alpha, mu[img], mu_mean) for img in mu}
    gamma_err = {img: gamma_error(alpha, mu[img], mu_err[img], mu_mean, len(mu)) for img in mu}

    # Signed arrival times relative to B (reference image per catalog):
    # dt_AB = t_A - t_B and dt_CB = t_C - t_B. Both are negative for H0pe.
    dt_ab_signed = float(h0pe["time_delays_days"]["dt_AB"]["value"])
    dt_cb_signed = float(h0pe["time_delays_days"]["dt_CB"]["value"])
    dt_ab = abs(dt_ab_signed)
    dt_cb = abs(dt_cb_signed)
    sigma_ab = (h0pe["time_delays_days"]["dt_AB"]["err_plus"] +
                h0pe["time_delays_days"]["dt_AB"]["err_minus"]) / 2
    sigma_cb = (h0pe["time_delays_days"]["dt_CB"]["err_plus"] +
                h0pe["time_delays_days"]["dt_CB"]["err_minus"]) / 2

    # With B as reference (t_B = 0), the TEP delay magnitude for image i is:
    # |dt_iB^TEP| = Gamma_i * |dt_iB^GR|
    # H0 shift factor for each delay is Gamma_i.
    f_ab = gamma["A"]   # A arrives before B
    f_cb = gamma["C"]   # C arrives before B

    # Weight by delay precision (inverse variance of H0 constraint ~ 1/sigma_dt^2)
    w_ab = 1.0 / sigma_ab ** 2
    w_cb = 1.0 / sigma_cb ** 2
    f_mean = (w_ab * f_ab + w_cb * f_cb) / (w_ab + w_cb)

    # Propagate gamma uncertainties into f_mean
    sigma_f = np.sqrt(
        (w_ab * gamma_err["A"]) ** 2 + (w_cb * gamma_err["C"]) ** 2
    ) / (w_ab + w_cb)

    h0_gr = H0_H0PE_ENC_GR
    err_plus = H0_H0PE_ENC_ERR_PLUS
    err_minus = H0_H0PE_ENC_ERR_MINUS
    h0_tep = h0_gr * f_mean
    sigma_h0_tep = h0_gr * sigma_f

    # 3-image closure residual (predicted, not observed)
    # R_TEP = (Gamma_A-1)*dt_AB + (Gamma_B-1)*dt_BC + (Gamma_C-1)*dt_CA
    # Loop orientation is A -> B -> C -> A:
    # dt_AB = t_B - t_A, dt_BC = t_C - t_B, dt_CA = t_A - t_C.
    dt_loop_ab = -dt_ab_signed
    dt_bc = dt_cb_signed
    dt_ca = dt_ab_signed - dt_cb_signed
    r_tep_h0pe = ((gamma["A"] - 1) * dt_loop_ab +
                  (gamma["B"] - 1) * dt_bc +
                  (gamma["C"] - 1) * dt_ca)

    return {
        "h0_gr": h0_gr,
        "h0_tep": h0_tep,
        "err_plus": err_plus,
        "err_minus": err_minus,
        "shift_factor": float(f_mean),
        "sigma_shift_factor": float(sigma_f),
        "sigma_h0_tep": float(sigma_h0_tep),
        "per_delay_factors": {
            "A_B": {"factor": float(f_ab), "sigma": float(gamma_err["A"]), "dt_days": float(dt_ab)},
            "C_B": {"factor": float(f_cb), "sigma": float(gamma_err["C"]), "dt_days": float(dt_cb)},
        },
        "gamma_factors": {img: float(gamma[img]) for img in gamma},
        "gamma_errors": {img: float(gamma_err[img]) for img in gamma_err},
        "predicted_closure_residual_days": float(r_tep_h0pe),
        "note": (
            "3-image weighted gamma shift. Shift is small (+~0.1%) because "
            "the high-magnification image (B, mu=8.0) is the reference, and "
            "the other images have modest magnification contrasts. "
            "Closure residual is predicted, not observed (only 2 independent delays)."
        ),
    }


def compute_h0_shift_encore(alpha):
    """Compute approximate TEP H0 shift for SN Encore (2-image system).

    Caveat: With only 2 images and 1 delay, the TEP shift is reference-dependent.
    The factor is reported using the later-arriving image (1a) as reference,
    which gives f = Gamma_1b. Using 1b as reference would give f = Gamma_1a.
    The difference is small (~1.5%) but non-zero. Both values and the average are reported.
    """
    catalog_path = PROJECT_ROOT / "data" / "raw" / "sn_lensing" / "lensed_sn_catalog.json"
    with open(catalog_path) as f:
        catalog = json.load(f)
    encore = catalog["sn_encore"]

    # Relative magnification: mu_1b / mu_1a = 2.0
    mu_ratio = encore["magnification_proxies"]["mu_relative"]["1b_1a"]["value"]
    mu_ratio_err = encore["magnification_proxies"]["mu_relative"]["1b_1a"]["err"]

    # Define absolute mu: mu_1a = 1, mu_1b = mu_ratio
    mu = {"1a": 1.0, "1b": mu_ratio}
    mu_err = {"1a": 0.0, "1b": mu_ratio_err}
    mu_mean = np.mean(list(mu.values()))

    gamma = {img: gamma_factor(alpha, mu[img], mu_mean) for img in mu}
    gamma_err = {img: gamma_error(alpha, mu[img], mu_err[img], mu_mean, len(mu)) for img in mu}

    # 1b arrives before 1a (dt = -37.3). With 1a as reference (t_1a=0), t_1b = -37.3.
    # dt_TEP = Gamma_1a * 0 - Gamma_1b * (-37.3) = Gamma_1b * 37.3
    # H0 shift factor = Gamma_1b (using 1a as reference)
    f_1a_ref = gamma["1b"]

    # With 1b as reference (t_1b=0), t_1a = +37.3.
    # dt_TEP = Gamma_1a * 37.3 - Gamma_1b * 0 = Gamma_1a * 37.3
    # H0 shift factor = Gamma_1a (using 1b as reference)
    f_1b_ref = gamma["1a"]

    # Average (geometric mean is more appropriate for ratios)
    f_mean = np.sqrt(f_1a_ref * f_1b_ref)
    sigma_f = 0.5 * f_mean * np.sqrt((gamma_err["1a"] / gamma["1a"]) ** 2 +
                                      (gamma_err["1b"] / gamma["1b"]) ** 2)

    h0_gr = H0_H0PE_ENC_GR
    err_plus = H0_H0PE_ENC_ERR_PLUS
    err_minus = H0_H0PE_ENC_ERR_MINUS
    h0_tep = h0_gr * f_mean
    sigma_h0_tep = h0_gr * sigma_f

    return {
        "h0_gr": h0_gr,
        "h0_tep": h0_tep,
        "err_plus": err_plus,
        "err_minus": err_minus,
        "shift_factor": float(f_mean),
        "sigma_shift_factor": float(sigma_f),
        "sigma_h0_tep": float(sigma_h0_tep),
        "per_reference_factors": {
            "1a_reference": float(f_1a_ref),
            "1b_reference": float(f_1b_ref),
        },
        "gamma_factors": {img: float(gamma[img]) for img in gamma},
        "gamma_errors": {img: float(gamma_err[img]) for img in gamma_err},
        "note": (
            "2-image approximate shift. Reference-dependent: f=1.0088 (1b ref) "
            "vs f=0.9938 (1a ref). The geometric mean is reported as a compromise. "
            "The shift is small (~+0.2%) because the magnification contrast "
            "(mu_1b/mu_1a = 2.0) is modest compared to Refsdal (~3.3)."
        ),
    }


def main():
    print_status(f"STEP {STEP_NUM}: H0 Tension Resolution (measured coupling)", "TITLE")

    # ------------------------------------------------------------------
    # 1. SN Refsdal: closure-based shift (gold standard)
    # ------------------------------------------------------------------
    refsdal = compute_h0_shift_refsdal(ALPHA)
    h0_r_gr = refsdal["h0_gr"]
    h0_r_tep = refsdal["h0_tep"]
    err_r_plus = refsdal["err_plus"]
    err_r_minus = refsdal["err_minus"]

    # ------------------------------------------------------------------
    # 2. SN H0pe: 3-image gamma-based shift
    # ------------------------------------------------------------------
    h0pe = compute_h0_shift_h0pe(ALPHA)
    h0_h_gr = h0pe["h0_gr"]
    h0_h_tep = h0pe["h0_tep"]
    err_h_plus = h0pe["err_plus"]
    err_h_minus = h0pe["err_minus"]

    # ------------------------------------------------------------------
    # 3. SN Encore: 2-image approximate shift
    # ------------------------------------------------------------------
    encore = compute_h0_shift_encore(ALPHA)
    h0_e_gr = encore["h0_gr"]
    h0_e_tep = encore["h0_tep"]
    err_e_plus = encore["err_plus"]
    err_e_minus = encore["err_minus"]

    # ------------------------------------------------------------------
    # Print results
    # ------------------------------------------------------------------
    print_status(f"\nSN Refsdal:")
    print_status(f"  GR H0 = {h0_r_gr:.1f} +/- {err_r_plus:.1f}/{err_r_minus:.1f}")
    print_status(f"  TEP shift factor = {refsdal['shift_factor']:.4f} +/- {refsdal['sigma_shift_factor']:.4f}")
    print_status(f"  TEP H0 = {h0_r_tep:.1f} +/- {refsdal['sigma_h0_tep']:.1f} (shift: {h0_r_tep - h0_r_gr:+.1f})")
    print_status(f"  dt_obs = {refsdal['dt_obs']:.1f} d, dt_geom = {refsdal['dt_geom']:.1f} d")

    print_status(f"\nSN H0pe:")
    print_status(f"  GR H0 = {h0_h_gr:.1f} +/- {err_h_plus:.1f}/{err_h_minus:.1f}")
    print_status(f"  TEP shift factor = {h0pe['shift_factor']:.4f} +/- {h0pe['sigma_shift_factor']:.4f}")
    print_status(f"  TEP H0 = {h0_h_tep:.1f} +/- {h0pe['sigma_h0_tep']:.1f} (shift: {h0_h_tep - h0_h_gr:+.1f})")
    print_status(f"  Per-delay: A-B f={h0pe['per_delay_factors']['A_B']['factor']:.4f}, C-B f={h0pe['per_delay_factors']['C_B']['factor']:.4f}")
    print_status(f"  Predicted closure residual = {h0pe['predicted_closure_residual_days']:.2f} d")

    print_status(f"\nSN Encore:")
    print_status(f"  GR H0 = {h0_e_gr:.1f} +/- {err_e_plus:.1f}/{err_e_minus:.1f}")
    print_status(f"  TEP shift factor = {encore['shift_factor']:.4f} +/- {encore['sigma_shift_factor']:.4f}")
    print_status(f"  TEP H0 = {h0_e_tep:.1f} +/- {encore['sigma_h0_tep']:.1f} (shift: {h0_e_tep - h0_e_gr:+.1f})")
    print_status(f"  Reference-dependent range: {encore['per_reference_factors']['1a_reference']:.4f} to {encore['per_reference_factors']['1b_reference']:.4f}")

    # ------------------------------------------------------------------
    # Combined Low-H0 cluster
    # ------------------------------------------------------------------
    w_r = 1.0 / ((err_r_plus + err_r_minus) / 2) ** 2
    w_h = 1.0 / ((err_h_plus + err_h_minus) / 2) ** 2
    w_e = 1.0 / ((err_e_plus + err_e_minus) / 2) ** 2

    h0_low_gr = (w_r * h0_r_gr + w_h * h0_h_gr + w_e * h0_e_gr) / (w_r + w_h + w_e)
    h0_low_tep = (w_r * h0_r_tep + w_h * h0_h_tep + w_e * h0_e_tep) / (w_r + w_h + w_e)

    print_status(f"\nCombined Low-H0 cluster (Ref+H0pe+Enc):")
    print_status(f"  GR:  {h0_low_gr:.1f}")
    print_status(f"  TEP: {h0_low_tep:.1f} (shift: {h0_low_tep - h0_low_gr:+.1f})")

    # ------------------------------------------------------------------
    # Planck tension check (full error budget)
    # ------------------------------------------------------------------
    # CRITICAL: tension with Planck must include BOTH the lens-model
    # uncertainty AND the Planck uncertainty. Using only Planck's tiny
    # ±0.5 would overstate the tension by ~8x.
    #
    # sigma_total^2 = sigma_lens^2 + sigma_planck^2
    #
    # For Refsdal: sigma_lens ≈ 4.1, sigma_planck = 0.5
    #   => sigma_total ≈ 4.13
    #   => tension = |69.3 - 67.4| / 4.13 ≈ 0.46 sigma (not 3.8!)
    # ------------------------------------------------------------------
    def tension_with_planck(h0_tep, err_lens):
        diff = h0_tep - H0_PLANCK
        sigma_total = np.sqrt(err_lens**2 + H0_PLANCK_ERR_PLUS**2)
        z = abs(diff) / sigma_total
        # Two-sided p-value: tension tests ask "are these measurements consistent?"
        # regardless of the direction of the discrepancy.
        p = float(2 * stats.norm.sf(z))
        return {"diff": diff, "sigma_total": sigma_total, "z": z, "p": p}

    tension_r = tension_with_planck(h0_r_tep, (err_r_plus + err_r_minus) / 2)
    tension_h = tension_with_planck(h0_h_tep, (err_h_plus + err_h_minus) / 2)
    tension_e = tension_with_planck(h0_e_tep, (err_e_plus + err_e_minus) / 2)

    print_status(f"\nPlanck tension check (FULL error budget):")
    print_status(f"  Refsdal: |{h0_r_tep:.1f} - {H0_PLANCK:.1f}| = {abs(tension_r['diff']):.1f} d")
    print_status(f"    sigma_total = sqrt({(err_r_plus+err_r_minus)/2:.1f}^2 + {H0_PLANCK_ERR_PLUS:.1f}^2) = {tension_r['sigma_total']:.2f}")
    print_status(f"    tension = {tension_r['z']:.2f} sigma  (p = {tension_r['p']:.3f})")
    print_status(f"  H0pe:    tension = {tension_h['z']:.2f} sigma")
    print_status(f"  Encore:  tension = {tension_e['z']:.2f} sigma")
    print_status(f"  NOTE: The Refsdal TEP shift does NOT resolve the H0 tension.")
    print_status(f"        The lens-model uncertainty (~4 km/s/Mpc) dominates.")

    # ------------------------------------------------------------------
    # Save results
    # ------------------------------------------------------------------
    results = {
        "step": STEP_NUM,
        "status": "success",
        "alpha": ALPHA,
        "sn_refsdal": {
            "gr": h0_r_gr,
            "tep": float(h0_r_tep),
            "err_plus": err_r_plus,
            "err_minus": err_r_minus,
            "sigma_tep": float(refsdal["sigma_h0_tep"]),
            "shift_factor": float(refsdal["shift_factor"]),
            "sigma_shift_factor": float(refsdal["sigma_shift_factor"]),
            "dt_obs_days": float(refsdal["dt_obs"]),
            "dt_geom_days": float(refsdal["dt_geom"]),
            "r_tep_days": float(refsdal["r_tep_days"]),
            "planck_tension": {
                "H0_planck": H0_PLANCK,
                "H0_planck_err": H0_PLANCK_ERR_PLUS,
                "diff_kms_mpc": float(tension_r["diff"]),
                "sigma_total_kms_mpc": float(tension_r["sigma_total"]),
                "z_sigma": float(tension_r["z"]),
                "p_value": float(tension_r["p"]),
                "caveat": "Tension uses full error budget: sigma_total = sqrt(sigma_lens^2 + sigma_planck^2). Using only Planck uncertainty would overstate tension by ~8x.",
            },
            "note": refsdal["note"],
        },
        "sn_h0pe": {
            "gr": h0_h_gr,
            "tep": float(h0_h_tep),
            "err_plus": err_h_plus,
            "err_minus": err_h_minus,
            "sigma_tep": float(h0pe["sigma_h0_tep"]),
            "shift_factor": float(h0pe["shift_factor"]),
            "sigma_shift_factor": float(h0pe["sigma_shift_factor"]),
            "per_delay_factors": h0pe["per_delay_factors"],
            "gamma_factors": h0pe["gamma_factors"],
            "predicted_closure_residual_days": float(h0pe["predicted_closure_residual_days"]),
            "planck_tension": {
                "diff_kms_mpc": float(tension_h["diff"]),
                "sigma_total_kms_mpc": float(tension_h["sigma_total"]),
                "z_sigma": float(tension_h["z"]),
                "p_value": float(tension_h["p"]),
            },
            "note": h0pe["note"],
        },
        "sn_encore": {
            "gr": h0_e_gr,
            "tep": float(h0_e_tep),
            "err_plus": err_e_plus,
            "err_minus": err_e_minus,
            "sigma_tep": float(encore["sigma_h0_tep"]),
            "shift_factor": float(encore["shift_factor"]),
            "sigma_shift_factor": float(encore["sigma_shift_factor"]),
            "per_reference_factors": encore["per_reference_factors"],
            "gamma_factors": encore["gamma_factors"],
            "planck_tension": {
                "diff_kms_mpc": float(tension_e["diff"]),
                "sigma_total_kms_mpc": float(tension_e["sigma_total"]),
                "z_sigma": float(tension_e["z"]),
                "p_value": float(tension_e["p"]),
            },
            "note": encore["note"],
        },
        "combined_low_h0": {
            "gr": float(h0_low_gr),
            "tep": float(h0_low_tep),
            "shift": float(h0_low_tep - h0_low_gr),
        },
        "independence_analysis": {
            "sn_refsdal": {
                "independence_level": "not_independent",
                "circularity_note": (
                    "The alpha_proxy=-0.055 used here was empirically determined from the "
                    "same SN Refsdal SX delay data (step_07). The H0 shift is therefore a "
                    "self-consistency check, not an independent confirmation. It shows that "
                    "the TEP framework is internally consistent, but does not add degrees of "
                    "freedom to the evidence."
                ),
            },
            "sn_h0pe": {
                "independence_level": "independent_prediction",
                "circularity_note": (
                    "alpha_proxy was calibrated on SN Refsdal, then applied to H0pe without "
                    "refitting. The predicted shift is an independent prediction, but the "
                    "shift is negligible (~0.1) due to modest magnification contrast."
                ),
            },
            "sn_encore": {
                "independence_level": "independent_prediction",
                "circularity_note": (
                    "alpha_proxy was calibrated on SN Refsdal, then applied to Encore without "
                    "refitting. The predicted shift is an independent prediction, but the "
                    "shift is negligible (~0.1) due to modest magnification contrast and "
                    "reference dependence of the 2-image system."
                ),
            },
        },
        "interpretation": (
            "SN Refsdal receives the dominant TEP H0 shift (+2.4) because it has "
            "5 images and the largest magnification contrast (S1/SX ~ 3.3). "
            "SN H0pe and SN Encore receive negligible shifts (+0.1 and +0.1 respectively) "
            "due to their smaller magnification contrasts and fewer independent delays. "
            "The combined low-H0 cluster shifts by +0.8 toward Planck, driven primarily by Refsdal. "
            "IMPORTANT: The Refsdal H0 shift is not an independent confirmation because alpha_proxy "
            "was derived from the same Refsdal data. It is an internal consistency check."
        ),
    }

    out_json = PROJECT_ROOT / "results" / "outputs" / f"step_{STEP_NUM}_h0_tension.json"
    with open(out_json, "w") as f:
        json.dump(results, f, indent=2)
    print_status(f"\nResults saved to {out_json}")

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scripts.utils.plot_style import set_pub_style, COLORS, FIG_SIZE

    set_pub_style()
    fig, ax = plt.subplots(figsize=FIG_SIZE)

    y_r, y_h, y_e = 3, 2, 1

    ax.errorbar(h0_r_gr, y_r, xerr=[[err_r_minus], [err_r_plus]], fmt="o", color=COLORS["gr"], label="GR", capsize=5)
    ax.errorbar(h0_h_gr, y_h, xerr=[[err_h_minus], [err_h_plus]], fmt="o", color=COLORS["gr"], capsize=5)
    ax.errorbar(h0_e_gr, y_e, xerr=[[err_e_minus], [err_e_plus]], fmt="o", color=COLORS["gr"], capsize=5)

    ax.errorbar(h0_r_tep, y_r, xerr=[[err_r_minus], [err_r_plus]], fmt="o", color=COLORS["tep"], label="TEP ($\\alpha_{\\rm lens}$)", capsize=5)
    ax.errorbar(h0_h_tep, y_h, xerr=[[err_h_minus], [err_h_plus]], fmt="o", color=COLORS["tep"], capsize=5)
    ax.errorbar(h0_e_tep, y_e, xerr=[[err_e_minus], [err_e_plus]], fmt="o", color=COLORS["tep"], capsize=5)

    ax.annotate("", xy=(h0_r_tep, y_r), xytext=(h0_r_gr, y_r),
                arrowprops=dict(arrowstyle="->", color="gray", lw=1.5))
    ax.annotate("", xy=(h0_h_tep, y_h), xytext=(h0_h_gr, y_h),
                arrowprops=dict(arrowstyle="->", color="gray", lw=1.5))
    ax.annotate("", xy=(h0_e_tep, y_e), xytext=(h0_e_gr, y_e),
                arrowprops=dict(arrowstyle="->", color="gray", lw=1.5))

    ax.text(52, y_r, "SN Refsdal", va="center", fontweight="bold")
    ax.text(52, y_h, "SN H0pe", va="center", fontweight="bold")
    ax.text(52, y_e, "SN Encore", va="center", fontweight="bold")

    ax.axvspan(66.8, 68.0, color="gray", alpha=0.2, label="Planck")
    ax.axvspan(72.0, 74.0, color=COLORS["shoes"], alpha=0.1, label="SH0ES")

    ax.set_yticks([])
    ax.set_xlabel("Hubble Constant $H_0$ [km s$^{-1}$ Mpc$^{-1}$]")
    ax.legend(loc="upper right")
    ax.set_title(f"TEP $H_0$ Correction ($\\alpha_{{\\rm lens}} = {ALPHA}$)")
    ax.set_xlim(50, 85)

    fig.tight_layout()
    out_fig = PROJECT_ROOT / "results" / "figures" / f"step_{STEP_NUM}_h0_tension.png"
    fig.savefig(out_fig)
    plt.close(fig)
    print_status(f"Saved plot to {out_fig}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
