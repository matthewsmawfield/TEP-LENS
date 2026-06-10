#!/usr/bin/env python3
"""
TEP-LENS: Step 03 - Proxy-Model GR Discrepancy / Sensitivity (SN Refsdal)

Under the log-magnification phenomenological proxy model, the effective transit
time of light scales with local magnification Gamma_t:

    t_obs_i = t_geom_i * Gamma_t(i)

where Gamma_t(i) = 1 + alpha * log10(mu_i), and mu_i is the lensing
magnification at image i (used as a proxy for the projected potential depth).
In the absence of a solved TEP lensing transfer function for kappa(theta)
or psi(theta), this is a first-order phenomenological proxy, not a derived
fundamental coupling.

For a closed loop i->j->k->i, the proxy-model predicted GR discrepancy is:

    d_proxy_GR = (Gamma_i - 1)*dt_ij + (Gamma_j - 1)*dt_jk + (Gamma_k - 1)*dt_ki

This is non-zero if and only if the Gamma_t values differ between images.
The MSD cannot mimic this, because any global mass sheet scales all delays
by the same factor, leaving the algebraic loop sum unchanged at zero under both
GR and the proxy model.

Magnification proxies: Kelly+2023 total relative flux ratios (F_i/F_ref),
which are direct proportional measures of the absolute magnification at
each image position.
"""

import json
import sys
from pathlib import Path
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status
from scripts.utils.tep_config import ALPHA_PROXY, SIGMA_ALPHA_PROXY

STEP_NUM = "03"

def safe_json_default(obj):
    if hasattr(obj, 'item'):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

def propagate_error(*errs):
    return float(np.sqrt(sum(e**2 for e in errs)))

def main():
    print_status(f"STEP {STEP_NUM}: Proxy-Model Predicted GR Discrepancy — SN Refsdal", "TITLE")

    catalog_path = PROJECT_ROOT / "data" / "raw" / "sn_lensing" / "lensed_sn_catalog.json"
    if not catalog_path.exists():
        print_status(f"Required catalog missing: {catalog_path}. Run step_01 first.", "ERROR")

    with open(catalog_path, 'r') as f:
        catalog = json.load(f)

    refsdal = catalog["sn_refsdal"]
    delays = refsdal["time_delays_days"]
    fluxes = refsdal["magnification_proxies"]["flux_total"]

    # Absolute arrival times relative to S1
    dt = {
        "S1": 0.0,
        "S2": delays["dt_S2_S1"]["value"],
        "S3": delays["dt_S3_S1"]["value"],
        "S4": delays["dt_S4_S1"]["value"],
        "SX": delays["dt_SX_S1"]["value"],
    }
    dt_err = {
        "S1": 0.0,
        "S2": delays["dt_S2_S1"]["err"],
        "S3": delays["dt_S3_S1"]["err"],
        "S4": delays["dt_S4_S1"]["err"],
        "SX": delays["dt_SX_S1"]["err"],
    }

    # Relative magnification proxies (flux ratio relative to mean)
    mu_rel = {img: fluxes[img]["value"] for img in fluxes}
    mu_ref = np.mean(list(mu_rel.values()))  # normalise to mean flux
    mu_norm = {img: mu_rel[img] / mu_ref for img in mu_rel}

    alpha_tep = ALPHA_PROXY  # Nominal illustrative coupling; not a pre-observation forecast

    # Gamma_t(i) = 1 + alpha * log10(mu_norm_i)
    Gamma = {img: 1.0 + alpha_tep * np.log10(mu_norm[img]) for img in mu_norm}

    print_status(f"Proxy-model coupling alpha = {alpha_tep} (Expansion/Screening)")
    print_status("Temporal shear factors Gamma_t per image:")
    for img in ["S1", "S2", "S3", "S4", "SX"]:
        print_status(f"  {img}: mu_rel={mu_rel[img]:.3f}, mu_norm={mu_norm[img]:.3f}, "
                     f"Gamma={Gamma[img]:.5f}")

    # ------------------------------------------------------------------
    # Closure loops
    # R_TEP = sum over edges (i->j) of: (Gamma_i - 1) * dt_ij
    # where dt_ij = dt[j] - dt[i]
    # ------------------------------------------------------------------
    loops = {
        "S1_S2_S3": ["S1", "S2", "S3"],
        "S1_S2_S4": ["S1", "S2", "S4"],
        "S1_S3_S4": ["S1", "S3", "S4"],
        "S1_S2_SX": ["S1", "S2", "SX"],
        "S1_S4_SX": ["S1", "S4", "SX"],
    }

    loop_results = {}
    print_status("Computing proxy-model predicted GR discrepancies...")

    for name, (i, j, k) in loops.items():
        dt_ij = dt[j] - dt[i]
        dt_jk = dt[k] - dt[j]
        dt_ki = dt[i] - dt[k]

        d_tep_gr = ((Gamma[i] - 1) * dt_ij
                 + (Gamma[j] - 1) * dt_jk
                 + (Gamma[k] - 1) * dt_ki)

        # Error propagation via partial derivatives w.r.t. FREE delays.
        # All delays are measured relative to image i (dt_i = 0, sigma_i = 0).
        # d = (G_i-1)*(dt_j-dt_i) + (G_j-1)*(dt_k-dt_j) + (G_k-1)*(dt_i-dt_k)
        # dd/d(dt_j) = (G_i-1) - (G_j-1) = G_i - G_j
        # dd/d(dt_k) = (G_j-1) - (G_k-1) = G_j - G_k
        # This avoids double-counting shared delay uncertainties in the loop edges.
        # CAVEAT: If delays are pairwise differences from a common reference (S1),
        # uncertainties may be correlated. Here we treat them as independent,
        # which is conservative when Gamma differences are small.
        dd_ddt_j = Gamma[i] - Gamma[j]
        dd_ddt_k = Gamma[j] - Gamma[k]
        d_err = propagate_error(
            dd_ddt_j * dt_err[j],
            dd_ddt_k * dt_err[k],
        )

        # Signal-to-noise of predicted discrepancy vs measurement sensitivity
        snr = abs(d_tep_gr) / d_err if d_err > 0 else 0.0

        loop_results[name] = {
            "images": [i, j, k],
            "Gamma": {i: float(Gamma[i]), j: float(Gamma[j]), k: float(Gamma[k])},
            "dt_ij_days": float(dt_ij),
            "dt_jk_days": float(dt_jk),
            "dt_ki_days": float(dt_ki),
            "tep_gr_discrepancy_days": float(d_tep_gr),
            "tep_gr_discrepancy_err_days": float(d_err),
            "snr": float(snr),
        }

        print_status(f"  Loop {i}-{j}-{k}: d_proxy_GR = {d_tep_gr:+.3f} ± {d_err:.3f} days "
                     f"(SNR = {snr:.2f})")

    # Summary statistics
    all_d = [v["tep_gr_discrepancy_days"] for v in loop_results.values()]
    all_snr = [v["snr"] for v in loop_results.values()]
    best_loop = max(loop_results, key=lambda k: loop_results[k]["snr"])

    print_status(f"Best loop by SNR: {best_loop} "
                 f"(SNR={loop_results[best_loop]['snr']:.2f})")

    # ------------------------------------------------------------------
    # Proxy validation: load step_32 kappa-proxy sensitivity if available
    # ------------------------------------------------------------------
    step32_path = PROJECT_ROOT / "results" / "outputs" / "step_32_kappa_proxy_validation.json"
    proxy_validation = None
    if step32_path.exists():
        try:
            with open(step32_path) as f:
                proxy_validation = json.load(f)["sn_refsdal_sensitivity"]
        except Exception:
            pass

    results = {
        "step": STEP_NUM,
        "status": "success",
        "system": "SN Refsdal (MACS J1149.6+2223)",
        "reference": "Kelly et al. 2023, ApJ 948, 93",
        "alpha_tep": alpha_tep,
        "gamma_per_image": {img: float(Gamma[img]) for img in Gamma},
        "mu_norm_per_image": {img: float(mu_norm[img]) for img in mu_norm},
        "tep_predicted_discrepancies": loop_results,
        "best_loop_by_snr": best_loop,
        "mean_tep_gr_discrepancy_days": float(np.mean(all_d)),
        "max_snr": float(max(all_snr)),
        "proxy_validation": proxy_validation,
    }

    out_dir = PROJECT_ROOT / "results" / "outputs"
    output_path = out_dir / f"step_{STEP_NUM}_tep_closure.json"
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2, default=safe_json_default)

    print_status(f"Results saved to {output_path}")
    print_status(f"Step {STEP_NUM} complete.")

if __name__ == "__main__":
    main()
