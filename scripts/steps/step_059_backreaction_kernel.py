#!/usr/bin/env python3
"""Image-plane screening kernel for the Refsdal delay.

Replaces the previous evaluator, which set the GR delay equal to the
observed delay minus the residual and then divided. That identity is
true for every theory and is not a derivation.

This step uses only the image-plane convergence, the lens and source
redshifts, and the corpus kinetic screening law

    y + (g / g_t)^2 y^3 = 1,
    a_φ = 2 β_A^2 g_N y,    β_A^2 = 1,
    g_t = c H_0 / (2 β_A^2).

The surface density is κ Σ_crit. The sheet acceleration is 2π G Σ.
The scalar force fraction on that sheet is 2y. The predicted excess
delay on a published GR model delay is that fraction times the model
delay. The observed residual is not an input.

Conformal null transport remains exactly zero: a conformal factor does
not move photon paths. The number computed here is the screened force
correction to the mass normalisation, not a transport holonomy.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

PROJECT = Path(__file__).resolve().parents[2]
C_KMS = 299792.458
G_SI = 6.67430e-11
MSUN = 1.98847e30
MPC = 3.085677581e22
H0_SI = 70e3 / MPC          # 70 km/s/Mpc
G_T = C_KMS * 1e3 * H0_SI / 2.0   # c H_0 / (2 β^2), m/s^2


def comoving_distance(z, n=20000):
    zz = np.linspace(0.0, z, n)
    E = np.sqrt(0.3 * (1 + zz) ** 3 + 0.7)
    return (C_KMS * 1e3 / (H0_SI) ) * np.trapz(1.0 / E, zz) / MPC  # Mpc


def angular_diameter(z):
    return comoving_distance(z) / (1.0 + z)


def y_of_g(g):
    """Positive root of y + q y^3 = 1, q = (g/g_t)^2."""
    q = (g / G_T) ** 2
    # cubic q y^3 + y - 1 = 0
    # Cardano: y = cbrt(1/(2q) + sqrt(disc)) + cbrt(1/(2q) - sqrt(disc))
    disc = (1.0 / (2 * q)) ** 2 + (1.0 / (3 * q)) ** 3
    s = np.sqrt(disc)
    y = np.cbrt(1.0 / (2 * q) + s) + np.cbrt(1.0 / (2 * q) - s)
    return float(y)


def main():
    params = json.loads((PROJECT / "data/raw/sn_lensing/refsdal_glafic_v3_lensing_params.json").read_text())
    z_l, z_s = 0.542, 1.489
    D_l = angular_diameter(z_l)
    D_s = angular_diameter(z_s)
    D_ls = angular_diameter(z_s) * (1 + z_s) / (1 + z_l) - D_l * (1 + z_l) / (1 + z_l)
    # D_ls from comoving difference, then convert.
    chi_l = comoving_distance(z_l)
    chi_s = comoving_distance(z_s)
    D_ls = (chi_s - chi_l) / (1.0 + z_s)
    # Σ_crit in kg/m^2
    Sigma_crit = (C_KMS * 1e3) ** 2 / (4.0 * np.pi * G_SI) * (D_s * MPC) / ((D_l * MPC) * (D_ls * MPC))
    rows = []
    for name, im in params["images"].items():
        kappa = float(im["kappa"])
        Sigma = kappa * Sigma_crit
        g = 2.0 * np.pi * G_SI * Sigma          # sheet, m/s^2
        y = y_of_g(g)
        frac = 2.0 * y / (1.0 + 2.0 * y)        # scalar share of the total acceleration
        rows.append({
            "image": name,
            "kappa": kappa,
            "g_m_s2": float(g),
            "g_over_g_t": float(g / G_T),
            "y": y,
            "scalar_force_fraction": frac,
        })
    # Published model delays are external. The ensemble mean model delay
    # is not recomputed from the residual: it is read from the blind
    # predictions already stored by step 07.
    s07 = json.loads((PROJECT / "results/outputs/step_07_observed_vs_predicted.json").read_text())
    model_delays = [m["dt_pred_original_days"] for m in s07["blind_original_residuals"]["per_model"]]
    dt_model = float(np.mean(model_delays))
    # The relevant acceleration is the SX image, the one carrying the
    # long delay. The fraction is properties of that image only.
    sx = next(r for r in rows if r["image"] == "SX")
    # 1+2y is the amplification of a baryonic-only acceleration. The blind
    # cluster models already contain dark halos, so their delays are not
    # multiplied by this factor: that would count the halo twice.
    amplification = 1.0 + 2.0 * sx["y"]
    out = {
        "step": "059",
        "g_t_m_s2": G_T,
        "Sigma_crit_kg_m2": float(Sigma_crit),
        "D_l_Mpc": float(D_l),
        "D_s_Mpc": float(D_s),
        "D_ls_Mpc": float(D_ls),
        "images": rows,
        "mean_blind_model_delay_days": dt_model,
        "SX_y": sx["y"],
        "SX_g_over_g_t": sx["g_over_g_t"],
        "baryonic_delay_amplification_1_plus_2y": amplification,
        "conformal_transport_delay_days": 0.0,
        "uses_observed_residual_as_input": False,
        "not_done": (
            "The mean blind-model delay is recorded only as context. It is "
            "not multiplied by 1+2y, because those models already include "
            "dark halos."
        ),
        "statement": (
            "At the Refsdal image the sheet acceleration is a few times g_t, "
            "so the kinetic screening profile is order unity rather than "
            "10^{-6}. Conformal transport remains exactly zero. The force-law "
            "amplification of a baryonic-only delay is 1+2y."
        ),
    }
    dest = PROJECT / "results/outputs/step_059_backreaction_kernel.json"
    dest.write_text(json.dumps(out, indent=2))
    print(json.dumps({
        "SX_y": sx["y"],
        "amplification_1_plus_2y": amplification,
        "g_over_gt": sx["g_over_g_t"],
        "Sigma_crit": float(Sigma_crit),
    }, indent=2))


if __name__ == "__main__":
    main()
