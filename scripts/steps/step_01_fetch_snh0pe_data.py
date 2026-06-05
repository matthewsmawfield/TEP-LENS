#!/usr/bin/env python3
"""
TEP-LENS: Step 01 - Compile Lensed Supernova Literature Data

Primary system: SN Refsdal (MACS J1149.6+2223), the only lensed SN with
five resolved images and precision-measured time delays (Kelly et al. 2023,
ApJ 948, 93). Four independent pairwise delays relative to S1 are available,
enabling multiple independent closure loops.

Secondary system: SN H0pe (PLCK G165.7+67.0), three images (Pierel et al. 2024).

Prediction target: SN 2025wny (z=2.011), four images, no time delays yet
(Johansson et al. 2025, ApJ 995, L17).
"""

import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status

STEP_NUM = "01"

def safe_json_default(obj):
    if hasattr(obj, 'item'):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

def main():
    print_status(f"STEP {STEP_NUM}: Compiling Lensed Supernova Literature Data", "TITLE")

    data_dir = PROJECT_ROOT / "data" / "raw" / "sn_lensing"
    data_dir.mkdir(parents=True, exist_ok=True)
    out_dir = PROJECT_ROOT / "results" / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # PRIMARY: SN REFSDAL (Kelly et al. 2023, ApJ 948, 93)
    # ------------------------------------------------------------------
    # Time delays are relative to S1 (S1 is the reference image, Dt > 0
    # means image arrives later than S1).
    # Source: Table 15, "Combined" method median values.
    # Symmetric 1-sigma errors from (84th-16th percentile)/2.
    # Magnification ratios (relative total flux) from same paper.
    # S1 flux (arbitrary units) taken as reference = 1.158 (total F).
    # Full-flux column: S1=1.158, S2=0.887, S3=0.716, S4=1.793, SX=~0.348
    # (SX magnification ratio mu_SX/mu_S1 ~ 0.30 from abstract).
    # ------------------------------------------------------------------
    print_status("Loading SN Refsdal data (Kelly et al. 2023)...")

    sn_refsdal = {
        "metadata": {
            "name": "SN Refsdal",
            "host_cluster": "MACS J1149.6+2223",
            "z_lens": 0.542,
            "z_src": 1.489,
            "images": ["S1", "S2", "S3", "S4", "SX"],
            "reference_image": "S1",
            "notes": (
                "Kelly et al. 2023, ApJ 948, 93. Time delays from Table 15 "
                "'Combined' method. Flux ratios from photometric analysis. "
                "S1-S4 form an Einstein cross; SX is a separate host arc image "
                "~8 arcsec away, appearing ~376 days after S1."
            )
        },
        "time_delays_days": {
            "dt_S2_S1": {"value":   9.9, "err": 4.0, "ref": "Kelly+2023 Table 15"},
            "dt_S3_S1": {"value":   9.0, "err": 4.2, "ref": "Kelly+2023 Table 15"},
            "dt_S4_S1": {"value":  20.3, "err": 6.4, "ref": "Kelly+2023 Table 15"},
            "dt_SX_S1": {"value": 376.0, "err": 5.6, "ref": "Kelly+2023 Table 15"}
        },
        "magnification_proxies": {
            "note": (
                "Total flux ratios (F_i / F_ref) used as relative magnification "
                "proxies. Values from Kelly+2023 photometric analysis. "
                "Absolute magnifications uncertain; ratios are robust."
            ),
            "flux_total": {
                "S1": {"value": 1.158, "err": 0.009},
                "S2": {"value": 0.887, "err": 0.007},
                "S3": {"value": 0.716, "err": 0.006},
                "S4": {"value": 1.793, "err": 0.017},
                "SX": {"value": 0.347, "err": 0.017}
            }
        }
    }

    # ------------------------------------------------------------------
    # SECONDARY: SN H0pe (Pierel et al. 2024, ApJ 967, 50)
    # Canonical photometric delays and magnifications used in the
    # blind-prediction residual test (Step 39) and manuscript.
    # Grayling et al. 2025 (BayeSN-TD) values are preserved in
    # data/raw/snh0pe/snh0pe_literature_data.json for reference.
    # ------------------------------------------------------------------
    print_status("Loading SN H0pe data (Pierel et al. 2024)...")

    sn_h0pe = {
        "metadata": {
            "name": "SN H0pe",
            "host_cluster": "PLCK G165.7+67.0",
            "z_lens": 0.351,
            "z_src": 1.783,
            "images": ["A", "B", "C"],
            "reference_image": "B",
            "notes": (
                "Pierel et al. 2024 (ApJ 967, 50). Photometric time delays and "
                "magnifications from JWST light-curve analysis."
            )
        },
        "time_delays_days": {
            "dt_AB": {"value": -116.6, "err_plus": 10.8, "err_minus": 9.3,
                      "ref": "Pierel+2024"},
            "dt_CB": {"value":  -48.6, "err_plus": 3.6, "err_minus": 4.0,
                      "ref": "Pierel+2024"}
        },
        "magnification_proxies": {
            "note": "Absolute magnifications from Pierel+2024 photometry.",
            "mu_absolute": {
                "A": {"value": 4.43, "err": 0.5},
                "B": {"value": 8.00, "err": 0.8},
                "C": {"value": 6.43, "err": 0.6}
            }
        }
    }

    # ------------------------------------------------------------------
    # SECONDARY: SN Encore (Pierel et al. 2026, ApJ, arXiv:2509.12301)
    # ------------------------------------------------------------------
    print_status("Loading SN Encore data (Pierel et al. 2026)...")

    sn_encore = {
        "metadata": {
            "name": "SN Encore",
            "host_cluster": "MACS J0138-2155",
            "z_lens": 0.338,
            "z_src": 1.95,
            "images": ["1a", "1b", "1c"],
            "reference_image": "1a",
            "notes": "Pierel et al. 2026 (ApJ, arXiv:2509.12301). Table 5 final systematic-corrected delay."
        },
        "time_delays_days": {
            "dt_1b_1a": {"value": -39.8, "err_plus": 3.9, "err_minus": 3.3,
                      "ref": "Pierel+2026 Table 5"}
        },
        "magnification_proxies": {
            "note": "Absolute magnifications from Pierel+2026 photometry (Table 5).",
            "mu_absolute": {
                "1a": {"value": 21.8, "err": 2.0},
                "1b": {"value": 32.4, "err": 3.0}
            },
            "mu_relative": {
                "1b_1a": {"value": 1.49, "err": 0.15}
            }
        }
    }

    # ------------------------------------------------------------------
    # PREDICTION TARGET: SN 2025wny (Johansson et al. 2025, ApJ 995, L17)
    # ------------------------------------------------------------------
    print_status("Loading SN 2025wny discovery data (Johansson et al. 2025)...")

    sn_2025wny = {
        "metadata": {
            "name": "SN 2025wny",
            "host_galaxy_lens": "PS1J0716+3821 (z_lens=0.3754)",
            "z_lens": 0.3754,
            "z_src": 2.011,
            "images": ["A", "B", "C", "D"],
            "image_separation_arcsec": 1.7,
            "notes": (
                "Johansson et al. 2025, ApJ 995, L17; Taubenberger et al. 2025 (HOLISMOKES XIX). "
                "First lensed SLSN-I. Four images in Einstein cross pattern. "
                "Host galaxy PS1J0716+3821 known quadruple since Cañameras et al. 2020. "
                "Time delays not yet measured. HST (PID 17611) and JWST (PID 5564) "
                "follow-up ongoing (PI: Goobar). High magnification contrast makes this "
                "the most promising future TEP test after Refsdal."
            )
        },
        "time_delays_days": {
            "note": "Not yet measured. Ongoing photometric monitoring by Maidanak, Lulin, COLIBRI, Wendelstein. HST/JWST follow-up scheduled.",
            "post_hoc_predictions": {
                "reference": "Witt-Wynne model, Wynne & Schechter 2018, as applied in arXiv:2605.11090. NOTE: post-hoc, not blind.",
                "dt_CA_days": {"value": -20, "note": "A trails C by ~20 d"},
                "dt_BA_days": {"value": -30, "note": "A trails B by ~30 d"},
                "dt_DA_days": {"value": -175, "note": "A trails D by ~175 d"}
            }
        },
        "magnification_proxies": {
            "mu_A_estimate": {"value_min": 20, "value_max": 50,
                              "ref": "Johansson+2025 from light curve comparison"},
            "note": "Very high magnification contrast expected. If mu_A ~ 20-50 and other images are mu ~ 1-5, contrast ratio could exceed 10x, comparable to Refsdal S4-SX."
        },
        "monitoring_status": {
            "time_delays_measured": False,
            "blind_lens_models_available": False,
            "hst_program": "PID 17611 (PI: Goobar)",
            "jwst_program": "PID 5564 (PI: Goobar)",
            "priority_for_tep": "HIGH — 4 images, high magnification contrast, long baseline (~175 d predicted)"
        }
    }

    # ------------------------------------------------------------------
    # VENUS DISCOVERY: SN Eos (Coulter et al. 2026)
    # ------------------------------------------------------------------
    print_status("Loading SN Eos discovery data (Coulter et al. 2026)...")

    sn_eos = {
        "metadata": {
            "name": "SN Eos",
            "host_cluster": "MACS J1931.8-2635",
            "z_lens": 0.35,
            "z_src": 5.133,
            "images": ["1", "2"],
            "notes": (
                "Coulter et al. 2026 (arXiv). Multiply-imaged Type II SN at z=5.13, "
                "discovered by VENUS in JWST/NIRCam imaging. Two confirmed images. "
                "Farthest spectroscopically confirmed SN ever discovered. "
                "Not viable for loop-closure TEP test (only 2 images). "
                "VENUS lens model (arXiv:2602.14074) available but no blind time-delay predictions."
            )
        },
        "time_delays_days": {
            "note": "Not yet measured. Only two images confirmed; no delay measurement possible yet."
        },
        "magnification_proxies": {
            "note": "Magnification not yet measured from light curves."
        },
        "monitoring_status": {
            "time_delays_measured": False,
            "blind_lens_models_available": False,
            "priority_for_tep": "LOW — only 2 images, no loop closure possible"
        }
    }

    catalog = {
        "sn_refsdal": sn_refsdal,
        "sn_h0pe": sn_h0pe,
        "sn_encore": sn_encore,
        "sn_2025wny": sn_2025wny,
        "sn_eos": sn_eos
    }

    raw_path = data_dir / "lensed_sn_catalog.json"
    with open(raw_path, 'w') as f:
        json.dump(catalog, f, indent=2, default=safe_json_default)
    print_status(f"Saved raw catalog to {raw_path}")

    output_path = out_dir / f"step_{STEP_NUM}_lensed_sn_catalog.json"
    with open(output_path, 'w') as f:
        json.dump({"step": STEP_NUM, "status": "success",
                   "primary_system": "SN Refsdal",
                   "n_images_primary": 5,
                   "n_measured_delays": 4,
                   "catalog": catalog},
                  f, indent=2, default=safe_json_default)

    print_status(f"Results saved to {output_path}")
    print_status(f"Step {STEP_NUM} complete.")

if __name__ == "__main__":
    main()
