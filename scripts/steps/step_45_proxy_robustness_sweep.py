#!/usr/bin/env python3
"""
TEP-LENS: Step 45 - Proxy-Robustness Sweep (Model-File + Theory-Space MC)

Purpose: Quantify the robustness of the headline directional evidence
(all blind residuals positive) to the choice of proxy tracer.

Two complementary analyses:

A. MODEL-FILE DISCOVERY
   Evaluate every published lens-model parameter file under flux, |mu|,
   and kappa tracers. Limited by data availability (only GLAFIC v3
   currently publishes image-position kappa values).

B. THEORY-SPACE MONTE CARLO (20,000 draws)
   Replicate step_32's shear-degeneracy envelope but compute the
   sign-match fraction: what fraction of physically allowed kappa
   configurations predict the same residual sign as observed?
   This samples the full theory space rather than counting models.

Tracers evaluated:
  (A) flux_ratio_proxy  : q = F_i / F_ref
  (B) model_mu_abs      : q = |mu_i|
  (C) model_kappa       : q = kappa_i

HONEST-REPORTING CONTRACT: this step does not rescue TEP if a tracer fails.
It measures how proxy-dependent the evidence is.

Inputs:
  - data/raw/sn_lensing/*_lensing_params.json
  - results/outputs/step_07_observed_vs_predicted.json
  - data/raw/sn_lensing/lensed_sn_catalog.json

Output:
  - results/outputs/step_45_proxy_robustness_sweep.json
"""

import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status
from scripts.utils.tep_config import ALPHA_PROXY

STEP_NUM = "45"
LOOP = ("S1", "S4", "SX")


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def loop_residual(alpha, quantity, delays, loop=LOOP):
    """
    Proxy-model loop residual.
    R = sum_{edges} (Gamma_i - 1) * dt_ij
    with Gamma_i = 1 + alpha * log10(q_i / qbar)
    """
    qbar = np.mean(list(quantity.values()))
    qn = {im: quantity[im] / qbar for im in quantity}
    G = {im: 1.0 + alpha * np.log10(qn[im]) for im in qn}
    i, j, k = loop
    R = (
        (G[i] - 1.0) * (delays[j] - delays[i])
        + (G[j] - 1.0) * (delays[k] - delays[j])
        + (G[k] - 1.0) * (delays[i] - delays[k])
    )
    return float(R), qn, G


def kappa_from_mu_gamma(mu, gamma):
    """Invert lensing identity: mu = 1 / [(1-kappa)^2 - gamma^2]."""
    term = 1.0 / mu + gamma ** 2
    if term < 0 or term > 1.0:
        return np.nan
    return 1.0 - np.sqrt(term)


def load_observed_residual():
    """Load the observed blind-prediction residual from step_07."""
    s07_path = PROJECT_ROOT / "results" / "outputs" / "step_07_observed_vs_predicted.json"
    with open(s07_path) as f:
        s07 = json.load(f)
    return float(s07["weighted_mean_residual"]["R_obs_days"])


def discover_model_files():
    """Find all lens-model parameter JSONs in data/raw/sn_lensing/."""
    data_dir = PROJECT_ROOT / "data" / "raw" / "sn_lensing"
    files = sorted(data_dir.glob("*_lensing_params.json"))
    return files


def process_model(model_path, observed_residual, flux_proxy):
    """Evaluate a single lens model under all three tracers."""
    with open(model_path) as f:
        model = json.load(f)

    model_name = model_path.stem.replace("_lensing_params", "")
    imgs = model.get("images", {})

    # Use model delays if present, otherwise canonical catalog
    delays = model.get("delays_days_rel_S1", {})
    if not delays or not all(im in delays for im in LOOP):
        # Fall back to catalog
        catalog_path = PROJECT_ROOT / "data" / "raw" / "sn_lensing" / "lensed_sn_catalog.json"
        with open(catalog_path) as f:
            catalog = json.load(f)
        refsdal = catalog["sn_refsdal"]
        delays = {
            "S1": 0.0,
            "S2": refsdal["time_delays_days"]["dt_S2_S1"]["value"],
            "S3": refsdal["time_delays_days"]["dt_S3_S1"]["value"],
            "S4": refsdal["time_delays_days"]["dt_S4_S1"]["value"],
            "SX": refsdal["time_delays_days"]["dt_SX_S1"]["value"],
        }

    # Build tracer dictionaries for images present in the model
    present = [im for im in LOOP if im in imgs]
    if len(present) < 3:
        print_status(
            f"Model {model_name} missing images for loop {LOOP}; skipping.",
            "WARNING",
        )
        return None

    # Tracer (A): flux proxy (always available from catalog)
    flux_sub = {im: flux_proxy[im] for im in present}

    # Tracer (B): parity-signed |mu|
    mu_abs = {}
    for im in present:
        mu_val = imgs[im].get("mu_signed")
        if mu_val is None:
            mu_abs = None
            break
        mu_abs[im] = abs(mu_val)

    # Tracer (C): kappa
    kappa = {}
    for im in present:
        k_val = imgs[im].get("kappa")
        if k_val is None:
            kappa = None
            break
        kappa[im] = k_val

    results = {"model_name": model_name}

    for tracer_name, q in (
        ("flux_ratio_proxy", flux_sub),
        ("model_mu_abs", mu_abs),
        ("model_kappa", kappa),
    ):
        if q is None:
            results[tracer_name] = {
                "available": False,
                "reason": "Required quantity not present in model file",
            }
            continue

        R, qn, G = loop_residual(ALPHA_PROXY, q, delays, loop=tuple(present))
        pred_obs = -R  # closure convention -> observed residual
        sign_match = bool(np.sign(pred_obs) == np.sign(observed_residual))

        results[tracer_name] = {
            "available": True,
            "R_closure_days": R,
            "R_predicted_obs_minus_model_days": pred_obs,
            "observed_residual_days": observed_residual,
            "sign_match": sign_match,
            "q_norm": {im: float(qn[im]) for im in qn},
            "Gamma": {im: float(G[im]) for im in G},
            "S4_SX_ordering": {
                "q_S4": float(q.get("S4", np.nan)),
                "q_SX": float(q.get("SX", np.nan)),
                "S4_gt_SX": bool(q.get("S4", 0) > q.get("SX", 0)),
            },
        }

    return results


def run_theory_space_mc(observed_residual, n_draws=20000, seed=20260607):
    """
    Replicate step_32's shear-degeneracy MC but collect sign-match stats.
    For each draw, infer kappa from flux + drawn shear, compute S1-S4-SX
    residual sign, and report the fraction matching the observed residual.
    """
    rng = np.random.default_rng(seed)

    fluxes = {"S1": 1.158, "S2": 0.887, "S3": 0.716, "S4": 1.793, "SX": 0.347}
    delays = {"S1": 0.0, "S2": 9.9, "S3": 9.0, "S4": 20.3, "SX": 376.0}
    gamma_prior = {
        "S1": {"mean": 0.60, "std": 0.18, "low": 0.05, "high": 0.92},
        "S2": {"mean": 0.60, "std": 0.18, "low": 0.05, "high": 0.92},
        "S3": {"mean": 0.60, "std": 0.18, "low": 0.05, "high": 0.92},
        "S4": {"mean": 0.60, "std": 0.18, "low": 0.05, "high": 0.92},
        "SX": {"mean": 0.18, "std": 0.12, "low": 0.02, "high": 0.70},
    }

    sign_match_flags = []
    pred_obs_vals = []

    for _ in range(n_draws):
        C = rng.uniform(0.5, 4.0)
        kappa_draw = {}
        for img in fluxes:
            p = gamma_prior[img]
            g = rng.normal(p["mean"], p["std"])
            g = np.clip(g, p["low"], p["high"])
            mu_abs = C * fluxes[img]
            if 1.0 / mu_abs + g ** 2 >= 1.0:
                g = np.sqrt(max(0, 1.0 - 1.0 / mu_abs - 1e-6))
            k = kappa_from_mu_gamma(mu_abs, g)
            if not np.isfinite(k) or k <= 0:
                k = 0.01
            kappa_draw[img] = k

        R, _, _ = loop_residual(ALPHA_PROXY, kappa_draw, delays, LOOP)
        pred_obs = -R
        sign_match_flags.append(np.sign(pred_obs) == np.sign(observed_residual))
        pred_obs_vals.append(pred_obs)

    pred_obs_vals = np.array(pred_obs_vals)
    return {
        "n_draws": n_draws,
        "seed": seed,
        "fraction_sign_match": float(np.mean(sign_match_flags)),
        "n_sign_match": int(np.sum(sign_match_flags)),
        "predicted_obs_residual": {
            "median": float(np.median(pred_obs_vals)),
            "mean": float(np.mean(pred_obs_vals)),
            "std": float(np.std(pred_obs_vals, ddof=1)),
            "p16": float(np.percentile(pred_obs_vals, 16)),
            "p84": float(np.percentile(pred_obs_vals, 84)),
            "p2.5": float(np.percentile(pred_obs_vals, 2.5)),
            "p97.5": float(np.percentile(pred_obs_vals, 97.5)),
        },
    }


def main():
    print_status(
        f"STEP {STEP_NUM}: Proxy-Robustness Sweep — SN Refsdal", "TITLE"
    )

    observed_residual = load_observed_residual()
    print_status(f"Observed residual (step_07): R_obs = {observed_residual:+.2f} d")

    # Canonical flux proxy from catalog
    catalog_path = PROJECT_ROOT / "data" / "raw" / "sn_lensing" / "lensed_sn_catalog.json"
    with open(catalog_path) as f:
        catalog = json.load(f)
    refsdal = catalog["sn_refsdal"]
    flux_proxy = {
        im: refsdal["magnification_proxies"]["flux_total"][im]["value"]
        for im in ["S1", "S2", "S3", "S4", "SX"]
    }

    model_files = discover_model_files()
    print_status(f"Discovered {len(model_files)} lens-model parameter file(s)")

    model_results = []
    for mf in model_files:
        res = process_model(mf, observed_residual, flux_proxy)
        if res is not None:
            model_results.append(res)
            print_status(f"  {res['model_name']}: processed")

    # Aggregate robustness across models
    tracer_names = ["flux_ratio_proxy", "model_mu_abs", "model_kappa"]
    robustness = {}
    for tname in tracer_names:
        n_available = sum(1 for r in model_results if r.get(tname, {}).get("available", False))
        n_match = sum(
            1
            for r in model_results
            if r.get(tname, {}).get("available", False) and r[tname]["sign_match"]
        )
        if n_available > 0:
            fraction = n_match / n_available
        else:
            fraction = None

        robustness[tname] = {
            "n_models_available": n_available,
            "n_models_sign_match": n_match,
            "fraction_sign_match": float(fraction) if fraction is not None else None,
            "interpretation": (
                f"{n_match}/{n_available} models predict the same residual sign "
                f"as the observed blind residual under the {tname} tracer."
            ),
        }

    # Per-tracer amplitude summary
    amplitude_summary = {}
    for tname in tracer_names:
        preds = [
            r[tname]["R_predicted_obs_minus_model_days"]
            for r in model_results
            if r.get(tname, {}).get("available", False)
        ]
        if preds:
            std_val = float(np.std(preds, ddof=1)) if len(preds) > 1 else None
            amplitude_summary[tname] = {
                "median_days": float(np.median(preds)),
                "mean_days": float(np.mean(preds)),
                "std_days": std_val,
                "min_days": float(np.min(preds)),
                "max_days": float(np.max(preds)),
            }
        else:
            amplitude_summary[tname] = {"available": False}

    # S4-SX ordering stability (fraction of models where S4 > SX)
    ordering_stability = {}
    for tname in tracer_names:
        flags = [
            r[tname]["S4_SX_ordering"]["S4_gt_SX"]
            for r in model_results
            if r.get(tname, {}).get("available", False)
        ]
        if flags:
            ordering_stability[tname] = {
                "fraction_S4_gt_SX": float(np.mean(flags)),
                "n_models": len(flags),
            }
        else:
            ordering_stability[tname] = {"available": False}

    # ------------------------------------------------------------------
    # B. Theory-space Monte Carlo (the scientifically meaningful part)
    # ------------------------------------------------------------------
    print_status("\nRunning theory-space Monte Carlo (20,000 kappa draws)...")
    mc_results = run_theory_space_mc(observed_residual, n_draws=20000, seed=20260607)
    p = mc_results["predicted_obs_residual"]
    print_status(
        f"  Theory-space sign-match: {mc_results['n_sign_match']}/{mc_results['n_draws']} "
        f"({mc_results['fraction_sign_match']:.1%})"
    )
    print_status(
        f"  Predicted residual (kappa MC): median={p['median']:+.2f} d, "
        f"mean={p['mean']:+.2f} d, std={p['std']:.2f} d"
    )
    print_status(
        f"  95% CI: [{p['p2.5']:+.2f}, {p['p97.5']:+.2f}] d"
    )

    results = {
        "step": STEP_NUM,
        "status": "success",
        "description": (
            "Systematic evaluation of the proxy-model predicted residual sign "
            "under multiple lens-model tracers. Measures how proxy-dependent "
            "the headline directional evidence is."
        ),
        "alpha_proxy": ALPHA_PROXY,
        "observed_residual_days": observed_residual,
        "loop_images": list(LOOP),
        "n_models_evaluated": len(model_results),
        "per_model_results": model_results,
        "robustness": robustness,
        "amplitude_summary": amplitude_summary,
        "S4_SX_ordering_stability": ordering_stability,
        "theory_space_mc": mc_results,
        "interpretation": (
            "The headline sign evidence (observed residual positive) is robust "
            "only if the majority of physically-motivated tracers predict the "
            "same sign. If convergence-based tracers systematically flip the sign, "
            "the evidence is specific to the flux proxy and cannot be claimed as "
            "a first-principles test of TEP without a solved lensing transfer function."
        ),
    }

    out_dir = PROJECT_ROOT / "results" / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"step_{STEP_NUM}_proxy_robustness_sweep.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=safe_json_default)

    print_status(f"\nRobustness summary:")
    for tname, rob in robustness.items():
        frac = rob["fraction_sign_match"]
        if frac is not None:
            print_status(f"  {tname}: sign-match = {rob['n_models_sign_match']}/{rob['n_models_available']} ({frac:.0%})")
        else:
            print_status(f"  {tname}: no models available")
    print_status(
        f"  theory_space_mc (kappa, 20k draws): sign-match = "
        f"{mc_results['n_sign_match']}/{mc_results['n_draws']} "
        f"({mc_results['fraction_sign_match']:.1%})"
    )

    print_status(f"\nResults saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
