#!/usr/bin/env python3
"""
TEP-LENS: Step 15 - External-Informed Uncertainty Inflation Test

Uses externally ingested H0LiCOW/TDCOSMO chain summaries (step_14) to define
an empirical uncertainty-inflation prior, then stress-tests the SN Refsdal
model-comparison metrics under this externally informed inflation.
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status

STEP_NUM = "15"


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def weighted_mean_sigma(vals: np.ndarray, sigmas: np.ndarray):
    w = 1.0 / sigmas**2
    m = float(np.sum(w * vals) / np.sum(w))
    s = float(1.0 / np.sqrt(np.sum(w)))
    return m, s


def evaluate_with_inflation(deltas, sigmas, r_tep, kappa):
    sig_eff = sigmas * np.sqrt(1.0 + kappa**2)

    m, s = weighted_mean_sigma(deltas, sig_eff)
    z = m / s
    p_z_one = float(stats.norm.sf(z))

    n_pos = int(np.sum(deltas > 0))
    n_tot = len(deltas)
    p_sign = float(stats.binomtest(n_pos, n_tot, 0.5, alternative="greater").pvalue)

    chi2_gr = float(np.sum((deltas / sig_eff) ** 2))
    chi2_tep = float(np.sum(((deltas - r_tep) / sig_eff) ** 2))
    dchi2 = chi2_gr - chi2_tep

    return {
        "kappa": float(kappa),
        "weighted_mean_days": float(m),
        "weighted_sigma_days": float(s),
        "z_from_gr": float(z),
        "p_z_one_sided": float(p_z_one),
        "p_sign_one_sided": float(p_sign),
        "delta_chi2_gr_minus_tep": float(dchi2),
    }


def main():
    print_status(f"STEP {STEP_NUM}: External-Informed Inflation", "TITLE")

    s07_path = PROJECT_ROOT / "results" / "outputs" / "step_07_observed_vs_predicted.json"
    s14_path = PROJECT_ROOT / "results" / "outputs" / "step_14_external_chain_ingestion.json"

    with open(s07_path) as f:
        s07 = json.load(f)
    with open(s14_path) as f:
        s14 = json.load(f)

    deltas = np.array([m["delta_obs_minus_pred_days"] for m in s07["per_model_results"]], dtype=float)
    sigmas = np.array([m["sigma_total_days"] for m in s07["per_model_results"]], dtype=float)
    r_tep = float(s07["tep_prediction"]["R_tep_prediction_days"])

    cvs = []
    for rec in s14.get("manifest", []):
        summ = rec.get("summary")
        if not rec.get("downloaded") or not summ:
            continue
        med = summ.get("col0_median")
        std = summ.get("col0_std")
        if med is None or std is None:
            continue
        if not np.isfinite(med) or not np.isfinite(std):
            continue
        if abs(med) < 1e-12:
            continue
        cvs.append(abs(std / med))

    if len(cvs) == 0:
        raise FileNotFoundError(
            f"No valid external-chain CV data found in {step14_path}.\n"
            "Run step_14_external_chain_ingestion.py first."
        )

    cvs = np.array(cvs, dtype=float)
    k16, k50, k84 = np.percentile(cvs, [16, 50, 84])

    print_status(
        f"External chain CV prior (from step_14): kappa16={k16:.3f}, "
        f"kappa50={k50:.3f}, kappa84={k84:.3f}"
    )

    grid = np.linspace(max(0.0, k16 * 0.5), max(k84 * 1.5, k50 + 0.05), 25)
    curve = [evaluate_with_inflation(deltas, sigmas, r_tep, float(k)) for k in grid]

    at16 = evaluate_with_inflation(deltas, sigmas, r_tep, float(k16))
    at50 = evaluate_with_inflation(deltas, sigmas, r_tep, float(k50))
    at84 = evaluate_with_inflation(deltas, sigmas, r_tep, float(k84))

    print_status(
        f"kappa50 scenario: z={at50['z_from_gr']:.2f}, p_z={at50['p_z_one_sided']:.3f}, "
        f"Δχ²={at50['delta_chi2_gr_minus_tep']:+.2f}"
    )

    out_fig = None
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from scripts.utils.plot_style import set_pub_style, COLORS, FIG_SIZES, save_fig

        set_pub_style()
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=FIG_SIZES["web_two_panel"])

        ax1.plot([c["kappa"] for c in curve], [c["z_from_gr"] for c in curve], lw=1.8, color=COLORS["gr"])
        ax1.axvline(k50, ls="--", color=COLORS["observed"], lw=1.0, label="External median kappa")
        ax1.set_xlabel("External-informed inflation kappa")
        ax1.set_ylabel("Weighted-mean z (vs GR)")
        ax1.legend()

        ax2.plot(
            [c["kappa"] for c in curve],
            [c["delta_chi2_gr_minus_tep"] for c in curve],
            lw=1.8,
            color=COLORS["tep"],
        )
        ax2.axhline(0, ls="--", color=COLORS["observed"], lw=1.0)
        ax2.axvline(k50, ls="--", color=COLORS["observed"], lw=1.0)
        ax2.set_xlabel("External-informed inflation kappa")
        ax2.set_ylabel(r"$\Delta\chi^2$ (GR - TEP)")

        out_fig = PROJECT_ROOT / "results" / "figures" / f"step_{STEP_NUM}_external_inflation.png"
        save_fig(fig, out_fig)
        print_status(f"Figure saved to {out_fig}")
    except Exception as e:
        print_status(f"Plotting failed: {e}", "ERROR")

    out = {
        "step": STEP_NUM,
        "status": "success",
        "external_cv_values": cvs.tolist(),
        "kappa_percentiles": {"p16": float(k16), "p50": float(k50), "p84": float(k84)},
        "scenarios": {"kappa16": at16, "kappa50": at50, "kappa84": at84},
        "curve": curve,
        "figure": str(out_fig) if out_fig else None,
    }

    out_path = PROJECT_ROOT / "results" / "outputs" / f"step_{STEP_NUM}_external_inflation.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=safe_json_default)

    print_status(f"Results saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
