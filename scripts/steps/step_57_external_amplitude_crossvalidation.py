#!/usr/bin/env python3
"""TEP-LENS Step 57 — external amplitude cross-validation.

This step asks the non-circular amplitude question for issue 19-1: does a
log-magnification response coefficient inferred without SN Refsdal predict the
Refsdal residual?  SN Encore and the primary SN H0pe AB contrast form the
external training set; Refsdal is held out.

The calculation also corrects a covariance error in the legacy cross-system
summaries.  A measured delay is common to every lens-model residual for a
given system, so its uncertainty must be added once after combining the model
predictions.  It must not be divided by sqrt(N_models).

No synthetic data are used.  Inputs are the archived outputs of Steps 07, 38,
and 39, which trace to the published Refsdal, Encore, and H0pe measurements and
blind lens-model predictions.
"""

import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status

STEP_NUM = "57"


def fixed_effect_prediction(models, excluded=()):
    """Combine model predictions only; observational error is handled later."""
    excluded = set(excluded)
    selected = [m for m in models if m["name"] not in excluded]
    values = np.array([m["dt_pred_days"] for m in selected], dtype=float)
    errors = np.array([m["sigma_model_days"] for m in selected], dtype=float)
    weights = 1.0 / errors**2
    mean = float(np.sum(weights * values) / np.sum(weights))
    sigma = float(1.0 / np.sqrt(np.sum(weights)))
    return mean, sigma, [m["name"] for m in selected]


def system_residual(observed, obs_sigma, models, excluded=()):
    pred, pred_sigma, names = fixed_effect_prediction(models, excluded=excluded)
    return {
        "observed_delay_days": float(observed),
        "observed_sigma_days": float(obs_sigma),
        "combined_model_prediction_days": pred,
        "combined_model_sigma_days": pred_sigma,
        "residual_days": float(observed - pred),
        "residual_sigma_days": float(np.hypot(obs_sigma, pred_sigma)),
        "models_used": names,
        "excluded_models": list(excluded),
        "covariance_rule": (
            "The common observed-delay error is added once after combining model "
            "predictions; it is not treated as independent across models."
        ),
    }


def fit_through_origin(residuals, sigmas, sensitivities):
    y = np.asarray(residuals, dtype=float)
    s = np.asarray(sigmas, dtype=float)
    x = np.asarray(sensitivities, dtype=float)
    w = 1.0 / s**2
    coefficient = float(np.sum(w * x * y) / np.sum(w * x**2))
    sigma = float(1.0 / np.sqrt(np.sum(w * x**2)))
    chi2 = float(np.sum(w * (y - coefficient * x) ** 2))
    return coefficient, sigma, chi2


def main():
    print_status(f"STEP {STEP_NUM}: External amplitude cross-validation", "TITLE")

    outdir = PROJECT_ROOT / "results" / "outputs"
    s07 = json.load(open(outdir / "step_07_observed_vs_predicted.json"))
    s38 = json.load(open(outdir / "step_38_sn_encore_residuals.json"))
    s39 = json.load(open(outdir / "step_39_sn_h0pe_residuals.json"))

    # Recompute one system-level residual per independent supernova.  The full
    # corrected Refsdal compilation is retained to address the +30 d issue.
    ref = system_residual(
        s07["observed"]["dt_SX_S1_days"],
        s07["observed"]["err_days"],
        s07["per_model_results"],
    )
    enc_obs = s38["observed"]
    encore = system_residual(
        enc_obs["dt_1b_1a_days"],
        0.5 * (enc_obs["err_plus"] + enc_obs["err_minus"]),
        s38["models"],
    )
    hab_obs = s39["delay_pair_AB"]["observed"]
    h0pe_ab = system_residual(
        hab_obs["value_days"],
        0.5 * (hab_obs["err_plus"] + hab_obs["err_minus"]),
        s39["delay_pair_AB"]["per_model"],
        excluded=("WSLAP+",),
    )

    sensitivities = {
        "Refsdal": float(s07["proxy_sensitivity"]["R_tep_unit_days_per_kappa"]),
        "Encore": float(s38["tep_prediction"]["R_tep_unit_days_per_kappa"]),
        "H0pe-AB": float(
            s39["tep_prediction"]["R_tep_AB_days"] / s39["tep_prediction"]["kappa_used"]
        ),
    }
    systems = {"Refsdal": ref, "Encore": encore, "H0pe-AB": h0pe_ab}
    for name, data in systems.items():
        data["unit_sensitivity_days_per_kappa"] = sensitivities[name]

    # Hold Refsdal out: train on one contrast from each independent external SN.
    train_names = ("Encore", "H0pe-AB")
    k_ext, k_ext_sigma, chi2_ext = fit_through_origin(
        [systems[n]["residual_days"] for n in train_names],
        [systems[n]["residual_sigma_days"] for n in train_names],
        [sensitivities[n] for n in train_names],
    )
    ref_pred = k_ext * sensitivities["Refsdal"]
    ref_pred_sigma = abs(sensitivities["Refsdal"]) * k_ext_sigma
    ref_pull = (ref["residual_days"] - ref_pred) / np.hypot(
        ref["residual_sigma_days"], ref_pred_sigma
    )

    # Reverse check: calibrate on Refsdal and predict the two external systems.
    k_ref = ref["residual_days"] / sensitivities["Refsdal"]
    k_ref_sigma = ref["residual_sigma_days"] / abs(sensitivities["Refsdal"])
    reverse = {}
    for name in train_names:
        prediction = k_ref * sensitivities[name]
        prediction_sigma = abs(sensitivities[name]) * k_ref_sigma
        pull = (systems[name]["residual_days"] - prediction) / np.hypot(
            systems[name]["residual_sigma_days"], prediction_sigma
        )
        reverse[name] = {
            "prediction_days": float(prediction),
            "prediction_sigma_days": float(prediction_sigma),
            "observed_residual_days": systems[name]["residual_days"],
            "pull_sigma": float(pull),
        }

    # Compare the proxy scaling against a generic delay-fraction scaling using
    # one statistically independent contrast per system.
    names = ("Refsdal", "Encore", "H0pe-AB")
    aligned = np.array([abs(systems[n]["residual_days"]) for n in names])
    sigma = np.array([systems[n]["residual_sigma_days"] for n in names])
    proxy_basis = np.array([abs(sensitivities[n]) for n in names])
    delay_basis = np.array([376.0, 39.8, 116.6])
    _, _, chi2_proxy = fit_through_origin(aligned, sigma, proxy_basis)
    _, _, chi2_fractional = fit_through_origin(aligned, sigma, delay_basis)
    aic_proxy = chi2_proxy + 2.0
    aic_fractional = chi2_fractional + 2.0

    verdict = (
        f"After holding Refsdal out and treating each observed-delay uncertainty as a "
        f"shared system-level error, Encore plus H0pe-AB give kappa_lens = "
        f"{k_ext:+.3f} +/- {k_ext_sigma:.3f}.  This predicts a Refsdal residual of "
        f"{ref_pred:+.1f} +/- {ref_pred_sigma:.1f} d, compared with the covariance-corrected "
        f"{ref['residual_days']:+.1f} +/- {ref['residual_sigma_days']:.1f} d; the difference "
        f"is {abs(ref_pull):.2f} sigma.  The external prediction is therefore statistically "
        f"consistent but too imprecise to close the amplitude test, and its central value is "
        f"larger than Refsdal by a factor {abs(ref_pred/ref['residual_days']):.1f}.  A generic "
        f"fractional-delay model remains modestly preferred over the log-magnification proxy "
        f"(Delta AIC = {aic_proxy-aic_fractional:.2f})."
    )
    print_status(verdict)

    result = {
        "step": STEP_NUM,
        "status": "success",
        "description": "Leave-one-system-out amplitude cross-validation with shared-observation covariance.",
        "systems": systems,
        "external_training": {
            "systems": list(train_names),
            "kappa_lens": k_ext,
            "sigma_kappa_lens": k_ext_sigma,
            "chi2": chi2_ext,
        },
        "held_out_refsdal_prediction": {
            "predicted_residual_days": float(ref_pred),
            "predicted_sigma_days": float(ref_pred_sigma),
            "observed_residual_days": ref["residual_days"],
            "observed_sigma_days": ref["residual_sigma_days"],
            "pull_sigma": float(ref_pull),
        },
        "reverse_refsdal_calibration": {
            "kappa_lens": float(k_ref),
            "sigma_kappa_lens": float(k_ref_sigma),
            "external_predictions": reverse,
        },
        "model_comparison": {
            "one_contrast_per_independent_system": list(names),
            "tep_proxy": {"chi2": chi2_proxy, "aic": aic_proxy},
            "uniform_fractional_delay": {"chi2": chi2_fractional, "aic": aic_fractional},
            "delta_aic_tep_minus_fractional": float(aic_proxy - aic_fractional),
        },
        "verdict": verdict,
    }
    path = outdir / "step_57_external_amplitude_crossvalidation.json"
    with open(path, "w") as f:
        json.dump(result, f, indent=2)
    print_status(f"Results saved to {path}")


if __name__ == "__main__":
    main()
