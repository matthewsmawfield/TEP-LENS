#!/usr/bin/env python3
"""
TEP-LENS: Step 44 - Direct-Convergence Loop Residual (model kappa / inv-kappa, not flux proxy)

Purpose: replace the flux-magnification proxy with the actual macro-model
convergence kappa (and parity-signed magnification mu) at the SN Refsdal image
positions, from the v3 GLAFIC model tabulated by Kelly et al. 2023. This is the
de-systematized test the referee asked for: it uses the quantity TEP is argued
to couple to (potential depth) rather than the microlensing-vulnerable flux
magnification.

Because beta_A = -1 locks the TEP scalar field to track the Newtonian potential
(Phi < 0 -> phi > 0 -> A = exp(-|phi|) < 1), deep potential means slower clocks.
For an isothermal-like profile the potential scales as ~1/kappa, so the
physically-motivated tracer is 1/kappa, not kappa itself.

It computes the proxy-model loop residual R_TEP(S1,S4,SX) four ways:
  (A) flux-ratio proxy        : Gamma = 1 + alpha*log10(|F|_norm)   [the paper's nominal]
  (B) parity-signed model mu  : Gamma = 1 + alpha*log10(|mu|_norm)  [model magnitudes, no microlensing]
  (C) model convergence kappa : Gamma = 1 + alpha*log10(kappa_norm) [density, NOT the potential]
  (D) model inv-kappa         : Gamma = 1 + alpha*log10((1/kappa)_norm) [potential proxy for beta_A=-1]

and compares the predicted sign/magnitude to the observed blind residual
(+30.1 d, i.e. R_TEP/GR ~ -14.5 d in closure convention).

HONEST-REPORTING CONTRACT: this step reports the residual under each tracer
regardless of whether it supports TEP. The GLAFIC convergence places SX at the
HIGHEST kappa (deepest), inverting the flux-proxy ordering, so the convergence-
based residual can have the opposite sign to the flux-based one. That outcome,
if it stands, means the headline sign agreement is specific to the flux proxy.

Inputs : data/raw/sn_lensing/refsdal_glafic_v3_lensing_params.json
         results/outputs/step_07_observed_vs_predicted.json (observed residual)
Outputs: results/outputs/step_44_direct_kappa_residual.json
"""

import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status
from scripts.utils.tep_config import ALPHA_PROXY

STEP_NUM = "44"
LOOP = ("S1", "S4", "SX")


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def loop_residual(alpha, quantity, delays, loop=LOOP):
    """R = sum (Gamma_i - 1) * dt around the loop; Gamma = 1 + alpha*log10(q/qbar)."""
    qbar = np.mean(list(quantity.values()))
    qn = {im: quantity[im] / qbar for im in quantity}
    G = {im: 1.0 + alpha * np.log10(qn[im]) for im in qn}
    i, j, k = loop
    return float((G[i] - 1.0) * (delays[j] - delays[i])
                 + (G[j] - 1.0) * (delays[k] - delays[j])
                 + (G[k] - 1.0) * (delays[i] - delays[k])), qn, G


def main():
    print_status(f"STEP {STEP_NUM}: Direct-Convergence Loop Residual (model kappa)", "TITLE")

    gl = json.load(open(PROJECT_ROOT / "data" / "raw" / "sn_lensing" / "refsdal_glafic_v3_lensing_params.json"))
    imgs = gl["images"]
    delays = gl["delays_days_rel_S1"]

    kappa = {im: imgs[im]["kappa"] for im in imgs}
    inv_kappa = {im: 1.0 / imgs[im]["kappa"] for im in imgs}
    mu_abs = {im: abs(imgs[im]["mu_signed"]) for im in imgs}
    # Flux-ratio proxy values (the paper's nominal F_i/F_ref).
    flux = {"S1": 1.158, "S2": 0.887, "S3": 0.716, "S4": 1.793, "SX": 0.347}

    s07 = json.load(open(PROJECT_ROOT / "results" / "outputs" / "step_07_observed_vs_predicted.json"))
    R_obs = float(s07["weighted_mean_residual"]["R_obs_days"])              # +30.1 (obs - model)
    R_tep_pred_obs_sign = -1.0  # closure convention: predicted (obs-model) ~ -R_closure

    print_status(f"GLAFIC v3 kappa: " + ", ".join(f"{im}={kappa[im]:.3f}" for im in imgs))
    print_status(f"  -> highest kappa: {max(kappa, key=kappa.get)} ({max(kappa.values()):.3f}); "
                 f"lowest: {min(kappa, key=kappa.get)} ({min(kappa.values()):.3f})")

    results_by_tracer = {}
    for name, q in (("flux_ratio_proxy", flux), ("model_mu_abs", mu_abs),
                     ("model_kappa", kappa), ("model_inv_kappa", inv_kappa)):
        R, qn, G = loop_residual(ALPHA_PROXY, q, delays)
        # predicted observed residual (obs - model) ~ -R_closure
        R_pred_obs = -R
        dGamma_S4_SX = G["S4"] - G["SX"]
        sign_matches_obs = bool(np.sign(R_pred_obs) == np.sign(R_obs))
        results_by_tracer[name] = {
            "R_closure_days": R,
            "R_predicted_obs_minus_model_days": R_pred_obs,
            "dGamma_S4_minus_SX": float(dGamma_S4_SX),
            "q_norm": {im: float(qn[im]) for im in qn},
            "predicted_sign_matches_observed": sign_matches_obs,
        }
        print_status(f"\n[{name}]")
        print_status(f"  q_norm S4={qn['S4']:.3f}  SX={qn['SX']:.3f}   dGamma(S4-SX)={dGamma_S4_SX:+.4f}")
        print_status(f"  R_closure={R:+.2f} d  -> predicted (obs-model)={R_pred_obs:+.2f} d  "
                     f"(observed={R_obs:+.2f} d)  sign-match: {sign_matches_obs}")

    flux_R = results_by_tracer["flux_ratio_proxy"]["R_predicted_obs_minus_model_days"]
    inv_kappa_R = results_by_tracer["model_inv_kappa"]["R_predicted_obs_minus_model_days"]
    kappa_R = results_by_tracer["model_kappa"]["R_predicted_obs_minus_model_days"]
    # Compare the physically-motivated tracer (1/kappa for beta_A=-1) against flux proxy
    sign_flip_phys = bool(np.sign(flux_R) != np.sign(inv_kappa_R))

    if sign_flip_phys:
        verdict = (
            "WARNING: under the physically-motivated INV-KAPPA tracer (1/kappa, appropriate for "
            f"beta_A=-1 potential coupling), the predicted residual is {inv_kappa_R:+.1f} d — "
            f"OPPOSITE in sign to the flux-proxy prediction ({flux_R:+.1f} d) and to the "
            f"observed blind residual ({R_obs:+.1f} d). The GLAFIC v3 maps place SX at the "
            "HIGHEST kappa (deepest potential in 1/kappa terms), which inverts the flux-proxy "
            "ordering. The headline sign agreement therefore depends on using flux magnification "
            "as the tracer; it does not survive substitution of the potential-proportional 1/kappa. "
            "The blind-residual FACT (models under-predict SX) is unchanged, but TEP's claim to "
            "PREDICT its sign via a potential coupling is not supported by these GLAFIC values."
        )
    else:
        verdict = (
            f"Inv-kappa residual ({inv_kappa_R:+.1f} d) shares the sign of the flux-proxy "
            f"prediction ({flux_R:+.1f} d) and the observation ({R_obs:+.1f} d); the sign "
            "evidence survives substitution of the 1/kappa potential proxy (amplitude differs)."
        )
    print_status("\n" + verdict)

    results = {
        "step": STEP_NUM,
        "status": "success",
        "description": ("Loop residual under flux-ratio proxy vs parity-signed model "
                        "magnification vs model convergence, using GLAFIC v3 values "
                        "(Kelly+2023). Tests whether the headline sign survives the mu->kappa "
                        "substitution the referee requires."),
        "alpha_proxy_ref": ALPHA_PROXY,
        "observed_blind_residual_days": R_obs,
        "glafic_v3_kappa": kappa,
        "kappa_ordering": {"highest": max(kappa, key=kappa.get), "lowest": min(kappa, key=kappa.get)},
        "results_by_tracer": results_by_tracer,
        "sign_flip_flux_vs_kappa": bool(np.sign(flux_R) != np.sign(kappa_R)),
        "sign_flip_flux_vs_inv_kappa": sign_flip_phys,
        "verdict": verdict,
        "provenance": gl["provenance"],
        "caveats": [
            "1/kappa is a proxy for potential depth under the isothermal assumption (beta_A=-1); "
            "the true TEP coupling variable is |Phi|, i.e. the lensing potential psi, which "
            "would be the cleanest tracer but is not tabulated in Kelly+2023.",
            "Total kappa includes the cluster-member-galaxy contribution at S1-S4; a cluster-"
            "only potential decomposition could shift the comparison and should be checked.",
            "GLAFIC v3 values are web-transcribed (see provenance); confirm against archived "
            "maps before publication.",
        ],
    }
    out_path = PROJECT_ROOT / "results" / "outputs" / f"step_{STEP_NUM}_direct_kappa_residual.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=safe_json_default)
    print_status(f"\nResults saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
