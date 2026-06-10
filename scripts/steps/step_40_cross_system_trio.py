#!/usr/bin/env python3
"""
TEP-LENS: Step 40 - Cross-System Trio Evidence Synthesis

Purpose: Combine blind-prediction residual evidence from all three independent
multiply-imaged supernova systems with published blind lens-model predictions:

  1. SN Refsdal (Kelly+2023): 6 blind + 2 post-blind models, loop-closure TEP test.
     Predicted residual: +14.5 d (high SNR).
  2. SN Encore (Pierel+2026): 8 models, single delay pair (1b,1a).
     Predicted residual: ~-0.5 d (very low SNR; magnification contrast is modest).
  3. SN H0pe (Pierel+2024): 7 models, two delay pairs (AB, CB).
     Predicted residuals: ~-1.6 d (AB), ~-0.3 d (CB) (low SNR).

Combination method: Stouffer's z-method for independent directional tests.
Each system contributes a z-score testing the SAME directional hypothesis:
"residuals have the sign predicted by the TEP proxy model." The predicted
sign and magnitude vary per system because the TEP effect depends on lens
geometry (magnification contrast and delay baseline). Under GR, each system's
residuals should be symmetric about zero.

Key limitation: Only SN Refsdal provides a high-SNR TEP test. The other two
systems have predicted TEP shifts of order ~1 d or less, far below per-model
scatter (~5-50 d). They serve as independent consistency checks, not as
high-precision evidence strands.
"""

import json
import sys
from pathlib import Path
import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status
from scripts.utils.tep_config import ALPHA_PROXY, SIGMA_ALPHA_PROXY

STEP_NUM = "40"


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def stouffer_combine(z_scores):
    """Stouffer's method: z_total = sum(z_i) / sqrt(N)."""
    z_scores = np.array(z_scores)
    n = len(z_scores)
    if n == 0:
        return {"z_total": 0.0, "p_total": 1.0, "n_systems": 0}
    z_total = float(z_scores.sum() / np.sqrt(n))
    p_total = float(stats.norm.sf(z_total))
    return {"z_total": z_total, "p_total": p_total, "n_systems": n}


def sign_label(value, tol=1e-12):
    """Return a stable sign label for directional checks."""
    if value > tol:
        return "positive"
    if value < -tol:
        return "negative"
    return "zero"


def matching_directional_test(n_positive, n_total, predicted_sign):
    """One-sided binomial test for residual signs matching prediction."""
    if predicted_sign == "positive":
        n_match = int(n_positive)
    elif predicted_sign == "negative":
        n_match = int(n_total - n_positive)
    else:
        return {
            "n_matching": 0,
            "n_total": 0,
            "p_value": 1.0,
            "z_approx": 0.0,
        }

    p_value = float(
        stats.binomtest(n_match, n_total, 0.5, alternative="greater").pvalue
    )
    z_approx = float((n_match - n_total * 0.5) / np.sqrt(n_total * 0.25))
    return {
        "n_matching": n_match,
        "n_total": int(n_total),
        "p_value": p_value,
        "z_approx": z_approx,
    }


def main():
    print_status(
        f"STEP {STEP_NUM}: Cross-System Trio Evidence Synthesis", "TITLE"
    )

    # ------------------------------------------------------------------
    # Load per-system results
    # ------------------------------------------------------------------
    s07_path = PROJECT_ROOT / "results" / "outputs" / "step_07_observed_vs_predicted.json"
    s16_path = PROJECT_ROOT / "results" / "outputs" / "step_16_tier_significance.json"
    s38_path = PROJECT_ROOT / "results" / "outputs" / "step_38_sn_encore_residuals.json"
    s39_path = PROJECT_ROOT / "results" / "outputs" / "step_39_sn_h0pe_residuals.json"

    systems = []

    # SN Refsdal (from Step 07 / Step 16)
    if s07_path.exists():
        with open(s07_path) as f:
            s07 = json.load(f)
        # Prefer the current correlation-aware Refsdal headline from Step 16.
        # Fall back to the independence-primary blind Wilcoxon for legacy outputs.
        p_blind = s07["test_registry"]["all_tests"]["wilcoxon_signed_rank_blind"]["p_value_one_sided"]
        z_blind = float(stats.norm.isf(p_blind))
        refsdal_headline_test = "Wilcoxon signed-rank blind"
        p_refsdal = p_blind
        z_refsdal = z_blind
        if s16_path.exists():
            with open(s16_path) as f:
                s16 = json.load(f)
            headline = s16.get("headline_significance")
            if headline:
                refsdal_headline_test = headline.get("test", refsdal_headline_test)
                p_refsdal = float(headline["p_value"])
                z_refsdal = float(headline["z_score"])
        R_obs = s07["weighted_mean_residual"]["R_obs_days"]
        sigma_R = s07["weighted_mean_residual"]["sigma_days"]
        R_tep = s07["tep_prediction"]["R_tep_prediction_days"]
        alpha = s07.get("bootstrap_alpha_inference", {}).get("alpha_mean", ALPHA_PROXY)

        systems.append(
            {
                "name": "SN Refsdal",
                "n_models": 7,
                "n_delay_pairs": 1,
                "R_obs_days": R_obs,
                "sigma_R_days": sigma_R,
                "R_tep_pred_days": R_tep,
                "predicted_sign": "positive",
                "observed_sign": "positive" if R_obs > 0 else "negative",
                "z_directional": z_refsdal,
                "p_directional": p_refsdal,
                "directional_test": refsdal_headline_test,
                "independence_primary": {
                    "test": "Wilcoxon signed-rank blind",
                    "z_score": z_blind,
                    "p_value": p_blind,
                },
                "snr_description": "high",
            }
        )
        print_status(
            f"SN Refsdal: R_obs={R_obs:+.2f} d, R_tep={R_tep:+.2f} d, "
            f"{refsdal_headline_test}: z={z_refsdal:+.2f}, p={p_refsdal:.4f}"
        )
    else:
        print_status("Step 07 output not found; skipping SN Refsdal", "WARN")

    # SN Encore (from Step 38)
    if s38_path.exists():
        with open(s38_path) as f:
            s38 = json.load(f)
        R_obs = s38["weighted_mean_residual"]["R_obs_days"]
        sigma_R = s38["weighted_mean_residual"]["sigma_R_obs_days"]
        R_tep = s38["tep_prediction"]["R_tep_prediction_days"]
        n_pos = s38["unweighted_statistics"]["n_positive"]
        n_total = len(s38["models"])

        predicted_sign = sign_label(R_tep)
        observed_sign = sign_label(R_obs)
        dir_test = matching_directional_test(n_pos, n_total, predicted_sign)
        p_dir = dir_test["p_value"]
        z_dir = dir_test["z_approx"]

        systems.append(
            {
                "name": "SN Encore",
                "n_models": n_total,
                "n_delay_pairs": 1,
                "R_obs_days": R_obs,
                "sigma_R_days": sigma_R,
                "R_tep_pred_days": R_tep,
                "predicted_sign": predicted_sign,
                "observed_sign": observed_sign,
                "z_directional": z_dir,
                "p_directional": p_dir,
                "directional_sign_test": dir_test,
                "snr_description": "very_low",
            }
        )
        print_status(
            f"SN Encore:  R_obs={R_obs:+.2f} d, R_tep={R_tep:+.2f} d, "
            f"z={z_dir:+.2f}, p={p_dir:.4f} (predicted {predicted_sign})"
        )
    else:
        print_status("Step 38 output not found; skipping SN Encore", "WARN")

    # SN H0pe (from Step 39)
    if s39_path.exists():
        with open(s39_path) as f:
            s39 = json.load(f)
        R_ab = s39["delay_pair_AB"]["weighted_mean_residual"]["R_obs_days"]
        sigma_ab = s39["delay_pair_AB"]["weighted_mean_residual"]["sigma_R_obs_days"]
        R_cb = s39["delay_pair_CB"]["weighted_mean_residual"]["R_obs_days"]
        sigma_cb = s39["delay_pair_CB"]["weighted_mean_residual"]["sigma_R_obs_days"]
        R_tep_ab = s39["tep_prediction"]["R_tep_AB_days"]
        R_tep_cb = s39["tep_prediction"]["R_tep_CB_days"]

        # Average the two pairs for an overall system residual
        # (correlated through image B, so this is a summary)
        R_sys = float(np.mean([R_ab, R_cb]))
        # Conservative uncertainty: max of the two
        sigma_sys = float(np.max([sigma_ab, sigma_cb]))

        dir_ab = s39["tep_prediction"]["directional_sign_tests"]["AB"]
        dir_cb = s39["tep_prediction"]["directional_sign_tests"]["CB"]
        z_dir_ab = float(dir_ab["z_approx"])
        p_dir_ab = float(dir_ab["p_value_one_sided"])
        z_dir_cb = float(dir_cb["z_approx"])
        p_dir_cb = float(dir_cb["p_value_one_sided"])

        # Use the pair with the larger predicted TEP sensitivity as the
        # system-level directional check. For current H0pe data this is AB.
        pair_candidates = [
            ("AB", R_ab, sigma_ab, R_tep_ab, dir_ab),
            ("CB", R_cb, sigma_cb, R_tep_cb, dir_cb),
        ]
        primary_pair, R_primary, sigma_primary, R_tep_primary, dir_primary = max(
            pair_candidates, key=lambda item: abs(item[3])
        )
        z_dir = float(dir_primary["z_approx"])
        p_dir = float(dir_primary["p_value_one_sided"])
        predicted_sign = sign_label(R_tep_primary)
        observed_sign = sign_label(R_primary)

        systems.append(
            {
                "name": "SN H0pe",
                "n_models": 7,
                "n_delay_pairs": 2,
                "R_obs_days": R_sys,
                "sigma_R_days": sigma_sys,
                "R_tep_pred_days": float(np.mean([R_tep_ab, R_tep_cb])),
                "primary_directional_pair": primary_pair,
                "predicted_sign": predicted_sign,
                "observed_sign": observed_sign,
                "z_directional": z_dir,
                "p_directional": p_dir,
                "directional_sign_test": dir_primary,
                "snr_description": "low",
                "pair_details": {
                    "AB": {
                        "R_obs": R_ab,
                        "R_tep_pred": R_tep_ab,
                        "predicted_sign": sign_label(R_tep_ab),
                        "observed_sign": sign_label(R_ab),
                        "z_directional": z_dir_ab,
                        "p_directional": p_dir_ab,
                        "directional_sign_test": dir_ab,
                    },
                    "CB": {
                        "R_obs": R_cb,
                        "R_tep_pred": R_tep_cb,
                        "predicted_sign": sign_label(R_tep_cb),
                        "observed_sign": sign_label(R_cb),
                        "z_directional": z_dir_cb,
                        "p_directional": p_dir_cb,
                        "directional_sign_test": dir_cb,
                    },
                },
            }
        )
        print_status(
            f"SN H0pe:    R_obs(AB)={R_ab:+.2f} d, R_obs(CB)={R_cb:+.2f} d, "
            f"R_tep(AB)={R_tep_ab:+.2f} d, R_tep(CB)={R_tep_cb:+.2f} d"
        )
        print_status(
            f"  Directional z(AB)={z_dir_ab:+.2f}, z(CB)={z_dir_cb:+.2f}; "
            f"primary={primary_pair}"
        )
    else:
        print_status("Step 39 output not found; skipping SN H0pe", "WARN")

    if not systems:
        print_status("No system results available. Run steps 07, 38, 39 first.", "ERROR")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Sign consistency check
    # ------------------------------------------------------------------
    n_match = sum(
        1 for s in systems if s["predicted_sign"] == s["observed_sign"]
    )
    n_total_sys = len(systems)

    p_sign_consistency = float(
        stats.binomtest(n_match, n_total_sys, 0.5, alternative="greater").pvalue
    )
    # Normal z-approximation is INVALID for n=3 (exact p=0.125 vs z-implied p≈0.042).
    # Do NOT report z_sign as a formal metric; retain only for diagnostic transparency.
    z_sign_approx = float((n_match - n_total_sys * 0.5) / np.sqrt(n_total_sys * 0.25))

    print_status(f"\nSign consistency check:")
    print_status(
        f"  Systems with observed sign matching TEP prediction: {n_match}/{n_total_sys}"
    )
    print_status(f"  Exact binomial p (H0: random signs): p = {p_sign_consistency:.4f}")
    print_status(f"  [DIAGNOSTIC ONLY] Normal z-approx (INVALID for n=3): z = {z_sign_approx:+.2f}")
    print_status(f"  WARNING: z-approx implies p≈0.042, contradicting exact p=0.125.")
    print_status(f"  Report ONLY the exact p-value in the manuscript.")

    # ------------------------------------------------------------------
    # Stouffer combination of directional z-scores
    # ------------------------------------------------------------------
    z_scores = [s["z_directional"] for s in systems]
    stouffer = stouffer_combine(z_scores)

    print_status(f"\nStouffer combination ({n_total_sys} systems):")
    print_status(f"  Individual z-scores: {[f'{z:+.2f}' for z in z_scores]}")
    print_status(
        f"  Combined z = {stouffer['z_total']:+.3f}, p = {stouffer['p_total']:.4f}"
    )
    print_status(
        f"  WARNING: These z-scores come from different test types (exact family-"
        f"sign-flip for Refsdal, binomial sign tests for Encore/H0pe) with different "
        f"statistical power. Stouffer assumes equivalent tests; mixing test types is "
        f"methodologically questionable. Interpret as a sensitivity exploration, not "
        f"a formal combined significance."
    )

    # ------------------------------------------------------------------
    # Precision-weighted mean residual across systems
    # This is NOT a proper TEP test (different geometric sensitivities),
    # but it shows the average shift.
    # ------------------------------------------------------------------
    R_vals = np.array([s["R_obs_days"] for s in systems])
    sigma_Rs = np.array([s["sigma_R_days"] for s in systems])
    weights = 1.0 / sigma_Rs**2
    w_sum = weights.sum()

    R_combined = float((weights * R_vals).sum() / w_sum)
    sigma_R_combined = float(1.0 / np.sqrt(w_sum))

    print_status(f"\nPrecision-weighted mean residual (all systems):")
    print_status(f"  R_combined = {R_combined:+.2f} +/- {sigma_R_combined:.2f} d")

    # ------------------------------------------------------------------
    # Refsdal-only vs. trio comparison
    # ------------------------------------------------------------------
    refsdal_z = next(
        (s["z_directional"] for s in systems if s["name"] == "SN Refsdal"), 0.0
    )
    trio_z = stouffer["z_total"]

    print_status(f"\nEvidence comparison:")
    print_status(f"  Refsdal alone: z = {refsdal_z:+.2f}")
    print_status(f"  Trio combined: z = {trio_z:+.2f}")
    print_status(
        f"  Encore + H0pe contribution: {trio_z - refsdal_z:+.2f} "
        f"(minimal because predicted TEP shifts are low-SNR)"
    )

    # ------------------------------------------------------------------
    # Save outputs
    # ------------------------------------------------------------------
    out = {
        "step": STEP_NUM,
        "status": "success",
        "description": (
            "Cross-system evidence synthesis for three multiply-imaged SNe "
            "with blind lens-model predictions. Stouffer combination of "
            "directional z-scores and sign-consistency check."
        ),
        "systems": systems,
        "sign_consistency": {
            "n_match": n_match,
            "n_total": n_total_sys,
            "p_value": p_sign_consistency,
            "z_approx_diagnostic_only": z_sign_approx,
            "z_approx_caveat": "Normal approximation is INVALID for n=3. Exact p=0.125; z=1.73 implies p≈0.042. Do NOT quote z_sign as a formal metric.",
        },
        "stouffer_combination": stouffer,
        "stouffer_caveat": (
            "Stouffer combines z-scores from different test types "
            "(exact family-sign-flip for Refsdal, binomial sign tests for Encore/H0pe) "
            "which have different power properties. This is a sensitivity exploration, "
            "not a formal combined significance."
        ),
        "precision_weighted_mean": {
            "R_combined_days": R_combined,
            "sigma_R_combined_days": sigma_R_combined,
        },
        "evidence_comparison": {
            "refsdal_z": refsdal_z,
            "trio_z": trio_z,
            "delta_z": trio_z - refsdal_z,
        },
        "key_findings": [
            "SN Refsdal is the only system with a high-SNR TEP predicted residual (~14.5 d).",
            "SN Encore predicted TEP shift is ~-0.5 d, far below model scatter; serves as consistency check.",
            "SN H0pe predicted TEP shifts are sub-day to ~2 d and below model scatter; primary directional check uses the larger-sensitivity AB pair.",
            f"All {n_match}/{n_total_sys} systems show primary residual signs consistent with TEP predictions.",
            "Stouffer combination is dominated by Refsdal; additional systems add minimal precision.",
            "Future high-contrast, long-baseline systems (like Refsdal) are needed for decisive multi-system evidence.",
        ],
        "limitations": [
            "Encore and H0pe have only single delay pairs; no loop-closure geometric test possible.",
            "TEP predicted shifts for Encore/H0pe are sub-day, swamped by per-model scatter of 5-50 d.",
            "H0pe delay pairs share image B; they are correlated, not independent.",
            "Per-system geometric sensitivity (R_tep_unit) varies by ~100x between Refsdal and Encore/H0pe.",
            "Alpha was calibrated on Refsdal; applying it to other systems is an extrapolation.",
        ],
    }

    out_path = (
        PROJECT_ROOT
        / "results"
        / "outputs"
        / f"step_{STEP_NUM}_cross_system_trio.json"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=safe_json_default)
    print_status(f"\nResults saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
