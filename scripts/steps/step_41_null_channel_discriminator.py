#!/usr/bin/env python3
"""
TEP-LENS: Step 41 - Null-Channel Discriminator (TEP vs. generic modelling bias)

Purpose: Directly confront the strongest mundane alternative to the TEP
interpretation — that the same-sign blind-prediction residuals across SN
Refsdal, Encore, and H0pe reflect a *generic* lens-model bias (e.g. an
H0/mass-sheet-like fractional offset, or a constant additive offset) rather
than potential-dependent temporal shear.

Method: in the TEP-predicted-positive sign frame, fit the cross-system
observed residuals to four competing one-parameter (or zero-parameter)
models and compare via chi-squared and AIC:

  (N) Null            : r_i = 0                       (GR; no systematic)
  (K) Constant offset : r_i = k                       (additive model bias)
  (F) Uniform fraction: r_i = c * |dt|_i              (H0 / MSD-like rescaling)
  (T) TEP proxy       : r_i = alpha * s_i             (s_i = R_tep_pred_i/|alpha_proxy|)

The discriminating question is whether the residual scales with the
TEP magnification-contrast sensitivity (model T) or merely with the delay
baseline (model F) / a constant (model K). This is the cross-system analogue
of a placebo test: model F/K carry NO potential-depth information, so if they
fit as well as T, the amplitude is not TEP-probative and only the *sign*
pattern (handled separately) carries evidence.

HONEST-REPORTING CONTRACT: this script reports whichever model the data
prefer, including outcomes that disfavour the TEP amplitude interpretation.
With only 3-4 contrasts the comparison is weak; the value is (i) quantifying
the TEP-vs-bias degeneracy now, and (ii) providing the discriminator that
sharpens as more long-baseline systems are added.

Inputs : results/outputs/step_40_cross_system_trio.json
Outputs: results/outputs/step_41_null_channel_discriminator.json
"""

import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status
from scripts.utils.tep_config import ALPHA_PROXY

STEP_NUM = "41"


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def fit_one_param(r, sigma, basis):
    """Weighted least-squares for r_i = theta * basis_i. Returns theta, chi2."""
    w = 1.0 / sigma**2
    theta = float(np.sum(w * basis * r) / np.sum(w * basis**2))
    chi2 = float(np.sum(w * (r - theta * basis) ** 2))
    return theta, chi2


def main():
    print_status(f"STEP {STEP_NUM}: Null-Channel Discriminator (TEP vs generic bias)", "TITLE")

    s40_path = PROJECT_ROOT / "results" / "outputs" / "step_40_cross_system_trio.json"
    if not s40_path.exists():
        raise FileNotFoundError(f"Required upstream output not found: {s40_path}")
    s40 = json.load(open(s40_path))

    # Per-contrast inputs (R_obs, sigma, R_tep_pred, baseline |dt|).
    baselines = {"SN Refsdal": 376.0, "SN Encore": 39.8}
    contrasts = []
    for s in s40["systems"]:
        name = s["name"]
        if name == "SN H0pe":
            for pair, bl in (("AB", 116.6), ("CB", 48.6)):
                pd = s["pair_details"][pair]
                contrasts.append({
                    "label": f"H0pe-{pair}",
                    "R_obs": float(pd["R_obs"]),
                    "sigma": float(s40["systems"][[x["name"] for x in s40["systems"]].index("SN H0pe")]
                                   ["pair_details"][pair].get("sigma_R", pd.get("sigma_R", np.nan))),
                    "R_tep_pred": float(pd["R_tep_pred"]),
                    "baseline": bl,
                })
        else:
            contrasts.append({
                "label": name.replace("SN ", ""),
                "R_obs": float(s["R_obs_days"]),
                "sigma": float(s["sigma_R_days"]),
                "R_tep_pred": float(s["R_tep_pred_days"]),
                "baseline": baselines[name],
            })

    # H0pe pair sigmas are not always stored in step_40; backfill from step_39.
    s39_path = PROJECT_ROOT / "results" / "outputs" / "step_39_sn_h0pe_residuals.json"
    if s39_path.exists():
        s39 = json.load(open(s39_path))
        sig_map = {}
        for pair_key, lab in (("delay_pair_AB", "H0pe-AB"), ("delay_pair_CB", "H0pe-CB")):
            wm = s39.get(pair_key, {}).get("weighted_mean", {}) or s39.get(pair_key, {})
            # try common locations
            for k in ("sigma_R_obs_days", "sigma_days", "sigma_R_days"):
                if k in wm:
                    sig_map[lab] = float(wm[k]); break
        for c in contrasts:
            if c["label"] in sig_map and (np.isnan(c["sigma"]) or c["sigma"] == 0):
                c["sigma"] = sig_map[c["label"]]

    # Hard fallback (documented): published weighted-mean sigmas.
    fallback_sigma = {"H0pe-AB": 5.13, "H0pe-CB": 2.67}
    for c in contrasts:
        if np.isnan(c["sigma"]) or c["sigma"] == 0:
            c["sigma"] = fallback_sigma.get(c["label"], c["sigma"])

    # Align to TEP-predicted-positive frame: r_i = R_obs_i * sign(R_tep_pred_i)
    # so the TEP-predicted direction is +. pred_aligned = |R_tep_pred|.
    for c in contrasts:
        s_dir = np.sign(c["R_tep_pred"]) or 1.0
        c["r_aligned"] = c["R_obs"] * s_dir
        c["pred_aligned"] = abs(c["R_tep_pred"])
        c["sign_match"] = bool(c["r_aligned"] > 0)

    labels = [c["label"] for c in contrasts]
    r = np.array([c["r_aligned"] for c in contrasts])
    sigma = np.array([c["sigma"] for c in contrasts])
    dts = np.array([c["baseline"] for c in contrasts])
    s_tep = np.array([c["pred_aligned"] for c in contrasts]) / abs(ALPHA_PROXY)  # TEP unit sensitivity

    print_status(f"Contrasts: {labels}")
    print_status(f"r_aligned (TEP+ frame): {np.round(r,2).tolist()}")
    print_status(f"sigma:                  {np.round(sigma,2).tolist()}")

    n = len(r)
    w = 1.0 / sigma**2

    # Model N: null
    chi2_N = float(np.sum(w * r**2)); k_N = 0
    # Model K: constant offset
    k_const = float(np.sum(w * r) / np.sum(w)); chi2_K = float(np.sum(w * (r - k_const) ** 2)); k_K = 1
    # Model F: uniform fractional (H0/MSD-like)
    c_frac, chi2_F = fit_one_param(r, sigma, dts); k_F = 1
    # Model T: TEP proxy
    alpha_fit, chi2_T = fit_one_param(r, sigma, s_tep); k_T = 1

    def aic(chi2, k):
        return chi2 + 2 * k

    models = {
        "null_GR":            {"params": k_N, "chi2": chi2_N, "aic": aic(chi2_N, k_N),
                               "best_fit": {}},
        "constant_offset":    {"params": k_K, "chi2": chi2_K, "aic": aic(chi2_K, k_K),
                               "best_fit": {"k_days": k_const}},
        "uniform_fractional": {"params": k_F, "chi2": chi2_F, "aic": aic(chi2_F, k_F),
                               "best_fit": {"c_fractional": c_frac,
                                            "note": "r_i = c*|dt|_i; H0/mass-sheet-like rescaling, no potential-depth info"}},
        "tep_proxy":          {"params": k_T, "chi2": chi2_T, "aic": aic(chi2_T, k_T),
                               "best_fit": {"alpha_fit": alpha_fit, "alpha_proxy_ref": ALPHA_PROXY}},
    }
    aic_min = min(m["aic"] for m in models.values())
    for name, m in models.items():
        m["delta_aic"] = float(m["aic"] - aic_min)

    ranked = sorted(models.items(), key=lambda kv: kv[1]["aic"])
    best_name = ranked[0][0]

    print_status("\nModel comparison (lower AIC = preferred):")
    for name, m in ranked:
        print_status(f"  {name:20s} chi2={m['chi2']:6.2f} (dof={n-m['params']})  "
                     f"AIC={m['aic']:6.2f}  dAIC={m['delta_aic']:5.2f}")

    # Amplitude-excess diagnostic: observed / TEP-predicted per contrast.
    excess = []
    for c in contrasts:
        ratio = c["r_aligned"] / c["pred_aligned"] if c["pred_aligned"] != 0 else float("nan")
        excess.append({"label": c["label"], "r_obs_aligned": c["r_aligned"],
                       "r_tep_pred_aligned": c["pred_aligned"], "obs_over_pred": ratio})
    print_status("\nAmplitude-excess (obs / TEP-pred) per contrast:")
    for e in excess:
        print_status(f"  {e['label']:10s} obs/pred = {e['obs_over_pred']:.1f}")

    # Sign test (the residual evidence channel, independent of amplitude).
    n_match = int(np.sum([c["sign_match"] for c in contrasts]))
    from scipy import stats as st
    p_sign = float(st.binomtest(n_match, n, 0.5, alternative="greater").pvalue)

    interpretation = []
    if models["uniform_fractional"]["aic"] <= models["tep_proxy"]["aic"]:
        interpretation.append(
            "AMPLITUDE NOT TEP-PROBATIVE: a uniform fractional (H0/mass-sheet-like) "
            "offset fits the cross-system residual magnitudes at least as well as the "
            "TEP magnification proxy. The TEP proxy under-predicts the Encore/H0pe "
            "residual magnitudes (they are several-day, while TEP predicts sub-2 d), "
            "so amplitude cannot distinguish TEP from a generic baseline-scaling bias.")
    else:
        interpretation.append(
            "TEP proxy is the AIC-preferred amplitude model across contrasts.")
    interpretation.append(
        f"Sign channel: {n_match}/{n} contrasts match the TEP-predicted direction "
        f"(one-sided binomial p = {p_sign:.3f}). Under a uniform fractional or constant "
        f"offset, the per-contrast SIGN is not predicted by potential-depth ordering, so "
        f"the directional match — not the amplitude — remains the carrier of TEP evidence.")
    interpretation.append(
        "With only %d contrasts this discriminator is weak; it sharpens as additional "
        "independent long-baseline systems are added (the same N that drives the "
        "Stouffer projection in Step 37)." % n)

    headline = (
        f"Cross-system discriminator (n={n}): AIC-preferred amplitude model = '{best_name}'. "
        f"TEP proxy dAIC = {models['tep_proxy']['delta_aic']:.2f}, "
        f"uniform-fractional dAIC = {models['uniform_fractional']['delta_aic']:.2f}. "
        f"Sign match {n_match}/{n} (p={p_sign:.3f}). "
        f"Amplitude does not favour TEP over a generic baseline-scaling bias; the sign "
        f"pattern is the residual evidence."
    )
    print_status("\n" + headline)

    results = {
        "step": STEP_NUM,
        "status": "success",
        "description": (
            "Cross-system discriminator separating the TEP magnification proxy from "
            "generic lens-model bias models (null, constant offset, uniform fractional "
            "H0/MSD-like). Amplitude vs sign decomposition."),
        "alpha_proxy_ref": ALPHA_PROXY,
        "contrasts": contrasts,
        "model_comparison": models,
        "aic_preferred_model": best_name,
        "amplitude_excess": excess,
        "sign_channel": {"n_match": n_match, "n_total": n, "p_one_sided": p_sign},
        "interpretation": interpretation,
        "headline": headline,
    }

    out_dir = PROJECT_ROOT / "results" / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"step_{STEP_NUM}_null_channel_discriminator.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=safe_json_default)
    print_status(f"\nResults saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
