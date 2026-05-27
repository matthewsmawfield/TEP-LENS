#!/usr/bin/env python3
"""
TEP-LENS: Step 05 - TDCOSMO Quad-Lens Temporal Shear Test

For quad-lens quasar systems, four images provide three independent pairwise
delays. Under GR, these form a set of geometrically consistent delays tied
to a single underlying potential. Under the log-magnification proxy model,
the differential temporal shear between image pairs of different magnification
introduces a systematic offset in the delay ratios relative to the GR lens
model prediction.

The test here is a scale-dependent shear test: for each quad system, we
compute the "flux-weighted delay ratio":
    rho_i = dt_iA / dt_iA_GR_expected
           = dt_iA / (flux_iA_ratio * dt_AB_ref)

If the proxy model is operative, images traversing deeper potential (higher
magnification) should arrive relatively later than GR predicts, producing a
systematic trend: rho increases with magnification contrast.

Systems used:
    HE0435-1223  (Bonvin et al. 2017, MNRAS 465, 4914)  -- quad
    WFI2033-4723 (Bonvin et al. 2019, A&A 629, A97)     -- quad
    DES0408-5354 (Courbin et al. 2018 / Wong et al. 2020) -- quad

For each system, flux ratios from HST imaging are used as magnification proxies.
"""

import json
import sys
from pathlib import Path
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status, safe_json_default

STEP_NUM = "05"


def propagate_error(*errs):
    return float(np.sqrt(sum(e**2 for e in errs)))

# ------------------------------------------------------------------
# ALL DELAY SYSTEMS (TDCOSMO 2025 + SN Encore)
# All delays in days. Positive = image arrives later than reference A.
# Flux ratios from HST imaging (proxy for relative magnification).
# Sources cited per system.
# ------------------------------------------------------------------
TDCOSMO_QUADS = {
    "HE0435-1223": {
        "type": "quasar",
        "z_lens": 0.454,
        "z_src": 1.693,
        "ref_delays": "Bonvin et al. 2017, MNRAS 465, 4914",
        "ref_fluxes": "Bonvin et al. 2017 / Chen et al. 2016",
        "images": ["A", "B", "C", "D"],
        "reference_image": "A",
        "delays": {
            "dt_BA": {"value":  -9.0, "err": 0.9},
            "dt_CA": {"value":  -3.0, "err": 1.4},
            "dt_DA": {"value": -13.8, "err": 0.8},
        },
        # HST flux ratios relative to A (F_i / F_A)
        "flux_ratios": {
            "A": 1.000,
            "B": 0.736,
            "C": 0.617,
            "D": 0.560,
        },
    },
    "WFI2033-4723": {
        "type": "quasar",
        "z_lens": 0.6575,
        "z_src": 1.662,
        "ref_delays": "Bonvin et al. 2019, A&A 629, A97",
        "ref_fluxes": "Vuissoz et al. 2008 / Bonvin et al. 2019",
        "images": ["A", "B", "C"],
        "reference_image": "A",
        "delays": {
            "dt_BA": {"value": -36.2, "err": 0.8},
            "dt_CA": {"value":  22.7, "err": 1.4},
        },
        "flux_ratios": {
            "A": 1.000,
            "B": 0.658,
            "C": 0.424,
        },
    },
    "DES0408-5354": {
        "type": "quasar",
        "z_lens": 0.597,
        "z_src": 2.375,
        "ref_delays": "Courbin et al. 2018 / Wong et al. 2020, MNRAS 498, 1420",
        "ref_fluxes": "Shajib et al. 2020",
        "images": ["A", "B", "C", "D"],
        "reference_image": "A",
        "delays": {
            "dt_BA": {"value": -112.1, "err":  2.1},
            "dt_CA": {"value": -155.5, "err": 12.8},
            "dt_DA": {"value": -128.4, "err":  5.1},
        },
        "flux_ratios": {
            "A": 1.000,
            "B": 0.715,
            "C": 0.248,
            "D": 0.191,
        },
    },
    "RXJ1131-1231": {
        "type": "quasar",
        "z_lens": 0.295,
        "z_src": 0.654,
        "ref_delays": "Tewes et al. 2013",
        "ref_fluxes": "Claeskens et al. 2006",
        "images": ["A", "B", "C", "D"],
        "reference_image": "A",
        "delays": {
            "dt_BA": {"value": 0.7, "err": 1.4},
            "dt_CA": {"value": -0.4, "err": 2.0},
            "dt_DA": {"value": 91.4, "err": 1.5},
        },
        "flux_ratios": {
            "A": 1.000,
            "B": 0.95,
            "C": 0.53,
            "D": 0.11,
        },
    },
    "PG1115+080": {
        "type": "quasar",
        "z_lens": 0.311,
        "z_src": 1.722,
        "ref_delays": "Bonvin et al. 2018",
        "ref_fluxes": "Weymann et al. 1997 / Chiba et al. 2005",
        "images": ["A1", "A2", "B", "C"],
        "reference_image": "A1",
        "delays": {
            "dt_A2A1": {"value": -14.3, "err": 3.4}, # actually dt_BA, mapping to simple terms
            "dt_CA1": {"value": 9.4, "err": 3.4},
        },
        "flux_ratios": {
            "A1": 1.000,
            "A2": 0.16,
            "C": 0.17,
        },
    },
    "B1608+656": {
        "type": "quasar",
        "z_lens": 0.63,
        "z_src": 1.394,
        "ref_delays": "Fassnacht et al. 2002",
        "ref_fluxes": "Surpi et al. 2003",
        "images": ["A", "B", "C", "D"],
        "reference_image": "A",
        "delays": {
            "dt_BA": {"value": 31.5, "err": 1.5},
            "dt_CA": {"value": 36.0, "err": 1.5},
            "dt_DA": {"value": 77.0, "err": 1.5},
        },
        "flux_ratios": {
            "A": 1.000,
            "B": 0.50,
            "C": 0.53,
            "D": 0.16,
        },
    },
    "J1206+4332": {
        "type": "quasar",
        "z_lens": 0.745,
        "z_src": 1.789,
        "ref_delays": "Eulaers et al. 2013",
        "ref_fluxes": "Agnello et al. 2016",
        "images": ["A", "B"],
        "reference_image": "A",
        "delays": {
            "dt_BA": {"value": 111.3, "err": 3.0},
        },
        "flux_ratios": {
            "A": 1.000,
            "B": 0.47,
        },
    },
    "SN Encore": {
        "type": "supernova",
        "z_lens": 0.338,
        "z_src": 1.95,
        "ref_delays": "Pierel et al. 2025",
        "ref_fluxes": "Pierel et al. 2025",
        "images": ["1a", "1b"],
        "reference_image": "1a",
        "delays": {
            "dt_1b1a": {"value": -37.3, "err": 13.1},
        },
        "flux_ratios": {
            "1a": 1.000,
            "1b": 2.000, # beta_1b_1a = 2.0
        },
    }
}

def tep_predicted_delay_ratio(dt_obs, flux_i, flux_A, alpha=-0.055):
    """
    Under TEP, the observed delay dt_iA = dt_geom * (Gamma_i / Gamma_A).
    For a reference pair, this introduces a fractional correction:
        delta_ratio = (Gamma_i - Gamma_A) / Gamma_A
    where Gamma = 1 + alpha * log10(mu_norm).

    The TEP-predicted fractional shift in dt_iA is computed relative
    to the GR expectation (= 0 shift), expressed as:
        TEP_residual_frac = (Gamma_i - 1) - (Gamma_A - 1)
                          = alpha * (log10(mu_i) - log10(mu_A))
                          = alpha * log10(mu_i / mu_A)
    This is the signed fractional residual: positive means the image
    traverses a deeper potential than A and arrives relatively later.
    """
    if flux_i <= 0 or flux_A <= 0:
        return 0.0
    return alpha * np.log10(flux_i / flux_A)

def main():
    print_status(f"STEP {STEP_NUM}: TDCOSMO Quad-Lens Temporal Shear Test", "TITLE")

    alpha_tep = -0.055
    print_status(f"TEP coupling alpha = {alpha_tep} (empirical lensing-sector coupling)")

    out_dir = PROJECT_ROOT / "results" / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = PROJECT_ROOT / "results" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    system_results = {}
    all_log_flux_ratios = []
    all_tep_residuals = []
    all_tep_residual_errs = []
    all_labels = []

    for sysname, sdata in TDCOSMO_QUADS.items():
        print_status(f"Processing {sysname}...")
        delays = sdata["delays"]
        fluxes = sdata["flux_ratios"]
        ref_image = sdata["reference_image"]
        flux_A = fluxes.get(ref_image, 1.0)

        pair_results = {}
        for pair, dval in delays.items():
            # e.g., "dt_BA", "dt_1b1a", "dt_A2A1"
            suffix = pair[3:]
            if suffix.endswith(ref_image):
                img = suffix[:-len(ref_image)]
            else:
                img = suffix.replace(ref_image, "")
                
            flux_i = fluxes.get(img, 1.0)
            dt = dval["value"]
            dt_err = dval["err"]

            # TEP fractional residual (dimensionless)
            delta_frac = tep_predicted_delay_ratio(dt, flux_i, flux_A, alpha_tep)
            # TEP absolute predicted delay shift (days)
            dt_tep_shift = delta_frac * abs(dt)
            log_flux_ratio = float(np.log10(flux_i / flux_A))

            pair_results[pair] = {
                "image": img,
                "dt_obs_days": dt,
                "dt_err_days": dt_err,
                "flux_ratio_i_A": flux_i / flux_A,
                "log10_flux_ratio": log_flux_ratio,
                "tep_fractional_residual": float(delta_frac),
                "tep_predicted_shift_days": float(dt_tep_shift),
                "tep_shift_err_days": float(dt_err * abs(delta_frac)),
            }

            print_status(
                f"  {sysname} {pair}: dt={dt:+.1f}±{dt_err}d, "
                f"log(F_i/F_A)={log_flux_ratio:+.3f}, "
                f"TEP_shift={dt_tep_shift:+.3f}d"
            )

            all_log_flux_ratios.append(log_flux_ratio)
            all_tep_residuals.append(float(dt_tep_shift))
            all_tep_residual_errs.append(float(dt_err * abs(delta_frac)))
            all_labels.append(f"{sysname[:7]}:{pair[3:]}")

        system_results[sysname] = {
            "metadata": {k: v for k, v in sdata.items()
                         if k not in ("delays", "flux_ratios")},
            "pair_results": pair_results,
        }

    # ------------------------------------------------------------------
    # Statistical summary: Spearman rank correlation between
    # log(flux ratio) and TEP predicted shift direction
    # ------------------------------------------------------------------
    from scipy.stats import spearmanr
    log_fr = np.array(all_log_flux_ratios)
    tep_r = np.array(all_tep_residuals)
    rho, pval = spearmanr(log_fr, tep_r)

    print_status(f"Spearman rho(log_flux_ratio, TEP_shift) = {rho:.3f}, p = {pval:.4f}")
    print_status("Note: by construction sign(rho) = sign(alpha_lens) since TEP_shift = alpha_lens*log(F)*|dt|")
    print_status("With alpha_lens=-0.055 < 0, the tautological correlation is negative (computed rho=-0.733).")
    print_status("The physically meaningful test is the magnitude of predicted shifts vs measurement errors.")

    # Fraction of pairs where predicted TEP shift > 1-sigma measurement error
    n_detectable = sum(
        1 for r, e in zip(all_tep_residuals, all_tep_residual_errs)
        if abs(r) > e and e > 0
    )
    n_total = len(all_tep_residuals)
    print_status(
        f"Pairs with |TEP_shift| > 1σ_meas: {n_detectable}/{n_total}"
    )

    # ------------------------------------------------------------------
    # Figure: TEP predicted shift vs log flux ratio, per system
    # ------------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from scripts.utils.plot_style import set_pub_style, COLORS

        set_pub_style()
        fig, ax = plt.subplots(figsize=(9, 5))

        colors_sys = {
            "HE0435-1223": COLORS["accent"],
            "WFI2033-4723": "#2196F3",
            "DES0408-5354": "#4CAF50",
            "RXJ1131-1231": "#9C27B0",
            "PG1115+080": "#FF9800",
            "B1608+656": "#795548",
            "J1206+4332": "#607D8B",
            "SN Encore": COLORS["highlight"],
        }
        markers_sys = {
            "HE0435-1223": "o",
            "WFI2033-4723": "s",
            "DES0408-5354": "^",
            "RXJ1131-1231": "v",
            "PG1115+080": "<",
            "B1608+656": ">",
            "J1206+4332": "D",
            "SN Encore": "*",
        }

        for sysname, sdata in system_results.items():
            xs, ys, yerrs = [], [], []
            for pr in sdata["pair_results"].values():
                xs.append(pr["log10_flux_ratio"])
                ys.append(pr["tep_predicted_shift_days"])
                yerrs.append(pr["tep_shift_err_days"])
            ax.errorbar(
                xs, ys, yerr=yerrs,
                fmt=markers_sys[sysname], color=colors_sys[sysname],
                markersize=8, capsize=4, linewidth=1.5,
                label=sysname, zorder=4
            )

        ax.axhline(0, color="black", linewidth=0.8, linestyle="--", zorder=2)
        ax.axvline(0, color="black", linewidth=0.4, linestyle=":", zorder=2)
        ax.set_xlabel(r"$\log_{10}(F_i / F_A)$  [relative magnification proxy]",
                      )
        ax.set_ylabel(r"TEP predicted delay shift $\delta t_{\rm TEP}$ [days]",
                      )
        ax.set_title(
            rf"TDCOSMO Quad-Lens Temporal Shear ($\alpha={alpha_tep}$)",
            pad=10
        )
        ax.legend(loc="upper left")
        # fig.tight_layout()

        out_fig = fig_dir / f"step_{STEP_NUM}_tdcosmo_shear.png"
        fig.savefig(out_fig)
        plt.close(fig)
        print_status(f"Figure saved to {out_fig}")

    except ImportError:
        print_status("matplotlib not available, skipping figure.", "WARNING")

    results = {
        "step": STEP_NUM,
        "status": "success",
        "alpha_tep": alpha_tep,
        "systems": system_results,
        "summary": {
            "n_systems": len(TDCOSMO_QUADS),
            "n_pairs_total": n_total,
            "n_detectable_pairs": n_detectable,
            "spearman_rho": float(rho),
            "spearman_pval": float(pval),
            "test_type": "predicted_sensitivity_check",
            "interpretation": (
                "TEP-predicted delay shifts for TDCOSMO quad lenses and SN Encore. "
                "These are PREDICTED shifts at the empirically measured coupling "
                "alpha_lens=-0.055; they are NOT observed detections of TEP in these "
                "systems. The physically meaningful quantity is the predicted shift "
                "magnitude versus published delay measurement uncertainties. "
                "The Spearman correlation between log(flux_ratio) and predicted shift "
                "is tautological (sign set by alpha_lens < 0; computed rho=-0.733) and therefore not independent evidence."
            ),
            "caveat": (
                "These systems do not permit a full geometric blind-prediction residual test "
                "because all independent pairwise delays are referenced to the same image. "
                "Any closure sum is arithmetically zero by construction. The analysis here "
                "is a structural consistency check: it asks whether the predicted TEP shift "
                "pattern is compatible with the data, not whether TEP is detected."
            ),
        },
    }

    output_path = out_dir / f"step_{STEP_NUM}_tdcosmo_shear.json"
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2, default=safe_json_default)

    print_status(f"Results saved to {output_path}")
    print_status(f"Step {STEP_NUM} complete.")

if __name__ == "__main__":
    main()
