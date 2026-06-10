#!/usr/bin/env python3
"""
TEP-LENS: Step 38 - CosmoGRAIL Cross-System Directional Consistency Check

Purpose: Test whether observed time delays in 18 CosmoGRAIL-monitored quadruply-
and triply-imaged quasar lens systems are directionally consistent with the TEP
proxy-model prediction that more magnified images arrive earlier (alpha_proxy < 0).

This is NOT an independent evidence strand. The test is confounded by lens
geometry: in strong lensing, more magnified images can naturally have either
shorter or longer geometric delays depending on the source position relative to
the caustic. The blind-prediction residual test (Step 07) remains the only
pre-specified primary evidence test because it compares observed delays to
lens-model-predicted geometric delays, removing the geometric confound.

What this test provides:
- A cross-system consistency check: if the TEP effect is real, it should
  manifest across many independent lens systems, not just SN Refsdal.
- A directional bias test: under GR, the delay-magnification correlation
  should average to zero across many systems (geometry is random). Under TEP,
  there should be a systematic bias.
- An upper limit: if no directional bias is seen across 18 systems, this
  constrains the amplitude of any universal magnification-dependent delay effect.

Method:
1. For each CosmoGRAIL system, load the light curve RDB file.
2. Compute mean flux per image from the photometric data (flux ~ 10^(-0.4*mag)).
   Time-averaged flux ratios approximate relative lensing magnifications.
3. Extract observed time delays between image pairs from the Step 30 output.
4. For each pair (i,j), compute the TEP-predicted sign:
   TEP(alpha<0) predicts: brighter image arrives EARLIER.
   So if mu_i > mu_j, TEP predicts dt_ji > 0 (j arrives later than i).
   Agreement condition: sign(dt_ij) == sign(mu_j - mu_i).
5. Count agreements across all pairs in all systems.
6. Binomial test under GR null: P(agree) = 0.5 by symmetry.
7. Per-system Kendall/Spearman rank correlation between delay and magnification.
8. Stouffer combination of per-system z-scores.

Caveats (must be stated clearly in output):
- Quasar fluxes include microlensing and intrinsic variability. Mean fluxes are
  imperfect magnification proxies.
- Lens geometry can accidentally align with the TEP prediction in some systems
  and against it in others. The geometric effect does not average to exactly
  zero; there may be weak systematic trends (e.g., cusp configurations).
- This test is a consistency check, not independent evidence. It uses the same
  alpha_proxy value calibrated on SN Refsdal.
- The systems are quasars, not supernovae. Different astrophysical systematics.
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

STEP_NUM = "38b"


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def parse_rdb_mean_fluxes(rdb_path, expected_images):
    """Parse a CosmoGRAIL RDB file and return mean flux per image.
    
    Handles two formats:
    - Files with header line (column names like mag_A, magerr_A)
    - Files without header (assumes standard column order)
    """
    with open(rdb_path, "r") as f:
        lines = f.readlines()

    # Detect header
    header_line = None
    type_line_idx = None
    for i, line in enumerate(lines):
        if "mag_" in line and "magerr_" in line:
            header_line = line.strip().split()
            type_line_idx = i + 1
            break
        elif "====" in line and "mag" in line:
            type_line_idx = i
            break

    if header_line is None:
        # No explicit header; infer from data structure
        # Try to find first non-comment, non-blank data line
        for i, line in enumerate(lines):
            if line.strip() and not line.strip().startswith("#"):
                parts = line.strip().split()
                n_cols = len(parts)
                # Standard format: time, (mag, err) per image, telescope
                # n_data_cols = n_cols - 1 (telescope) or n_cols
                n_data_cols = n_cols - 1 if parts[-1].isalpha() or len(parts[-1]) > 4 else n_cols
                n_images = (n_data_cols - 1) // 2
                header_line = ["time"] + [f"mag_{chr(65+j)}" for j in range(n_images)] + \
                              [f"magerr_{chr(65+j)}" for j in range(n_images)]
                if n_data_cols % 2 == 0:
                    header_line.append("telescope")
                break
        if header_line is None:
            return None

    # Map image labels to column indices
    mag_cols = {}
    for img in expected_images:
        mag_key = f"mag_{img}"
        if mag_key in header_line:
            mag_cols[img] = header_line.index(mag_key)

    if not mag_cols:
        return None

    # Read data
    fluxes = {img: [] for img in mag_cols}
    start_idx = type_line_idx + 1 if type_line_idx is not None else 0
    for line in lines[start_idx:]:
        parts = line.strip().split()
        if not parts:
            continue
        for img, col_idx in mag_cols.items():
            if col_idx < len(parts):
                try:
                    mag = float(parts[col_idx])
                    flux = 10.0 ** (-0.4 * mag)
                    fluxes[img].append(flux)
                except ValueError:
                    continue

    if not any(fluxes.values()):
        return None

    # Return mean flux per image
    mean_fluxes = {}
    for img, flist in fluxes.items():
        if flist:
            mean_fluxes[img] = float(np.mean(flist))

    return mean_fluxes


def tep_predicted_sign(dt_ij, mu_i, mu_j, alpha=-0.055):
    """Return +1 if the pair agrees with TEP prediction, -1 if disagrees, 0 if tie.
    
    TEP(alpha<0): brighter image (higher mu) arrives earlier.
    If mu_i > mu_j: TEP predicts i arrives earlier, so dt_ij = t_i - t_j < 0.
    Agreement: sign(dt_ij) == sign(mu_j - mu_i)
    """
    if dt_ij == 0 or mu_i == mu_j:
        return 0
    obs_sign = np.sign(dt_ij)
    tep_sign = np.sign(mu_j - mu_i)
    if obs_sign == tep_sign:
        return 1
    else:
        return -1


def main():
    print_status(
        f"STEP {STEP_NUM}: CosmoGRAIL Cross-System Directional Consistency Check",
        "TITLE",
    )

    # ------------------------------------------------------------------
    # Load Step 30 outputs (observed delays)
    # ------------------------------------------------------------------
    s30_path = PROJECT_ROOT / "results" / "outputs" / "step_30_cosmograil_temporal_shear.json"
    if not s30_path.exists():
        print_status(f"Step 30 output not found: {s30_path}", "ERROR")
        return

    with open(s30_path) as f:
        s30 = json.load(f)

    # ------------------------------------------------------------------
    # Process each CosmoGRAIL system
    # ------------------------------------------------------------------
    data_dir = PROJECT_ROOT / "data" / "cosmograil"
    rdb_files = sorted(data_dir.glob("*.rdb"))

    system_results = []
    all_pair_results = []
    skipped_systems = []

    for sys_name, sys_data in s30["systems"].items():
        # Extract base system name (remove band suffix if present)
        # e.g., "HE0435_R" -> "HE0435", "DESJ0408" -> "DESJ0408"
        base_name = sys_name.split("_")[0]

        # Find matching RDB file
        rdb_match = None
        for rdb in rdb_files:
            rdb_stem = rdb.stem.upper()
            if base_name.upper() in rdb_stem:
                rdb_match = rdb
                break

        if rdb_match is None:
            skipped_systems.append(f"{sys_name}: no RDB file")
            continue

        images = sys_data.get("image_labels", [])
        pairs = sys_data.get("pairs", {})

        if len(images) < 2:
            skipped_systems.append(f"{sys_name}: only {len(images)} images")
            continue

        # Parse mean fluxes
        mean_fluxes = parse_rdb_mean_fluxes(rdb_match, images)
        if mean_fluxes is None or len(mean_fluxes) < 2:
            skipped_systems.append(f"{sys_name}: could not parse fluxes")
            continue

        # Normalise to mean flux = 1.0 (same convention as Step 03)
        mu_mean = np.mean(list(mean_fluxes.values()))
        mu_norm = {img: mean_fluxes[img] / mu_mean for img in mean_fluxes}

        # Compute Gamma_t per image
        alpha = ALPHA_PROXY
        Gamma = {img: 1.0 + alpha * np.log10(mu_norm[img]) for img in mu_norm}

        # Process each pair
        pair_results = []
        for pair_key, pair_data in pairs.items():
            imgs = pair_key.split("-")
            if len(imgs) != 2:
                continue
            img_i, img_j = imgs[0], imgs[1]

            if img_i not in mu_norm or img_j not in mu_norm:
                continue

            delay_info = pair_data.get("broadband", {})
            if "delay_days" not in delay_info:
                continue

            dt_ij = float(delay_info["delay_days"])
            dt_err = float(delay_info.get("uncertainty_days", 0.0))

            mu_i = mu_norm[img_i]
            mu_j = mu_norm[img_j]

            agree = tep_predicted_sign(dt_ij, mu_i, mu_j, alpha)

            pair_results.append({
                "pair": pair_key,
                "dt_ij_days": dt_ij,
                "dt_err_days": dt_err,
                "mu_i": float(mu_i),
                "mu_j": float(mu_j),
                "Gamma_i": float(Gamma[img_i]),
                "Gamma_j": float(Gamma[img_j]),
                "agrees_with_tep": int(agree),
            })

            all_pair_results.append({
                "system": sys_name,
                "pair": pair_key,
                "dt_ij_days": dt_ij,
                "mu_i": float(mu_i),
                "mu_j": float(mu_j),
                "agrees_with_tep": int(agree),
            })

        if not pair_results:
            skipped_systems.append(f"{sys_name}: no valid pairs")
            continue

        # Per-system statistics
        n_agree = sum(1 for p in pair_results if p["agrees_with_tep"] == 1)
        n_disagree = sum(1 for p in pair_results if p["agrees_with_tep"] == -1)
        n_total = n_agree + n_disagree

        # Directional agreement fraction per system
        # NOTE: P(agree | GR) is NOT necessarily 0.5. Lens geometry can create
        # systematic correlations between magnification and delay (e.g., images near
        # the Einstein radius are both more magnified and have shorter delays).
        # The binomial p=0.5 null is a simplifying assumption; the true null is
        # lens-model-dependent and generally unknown without blind predictions.
        if n_total > 0:
            p_binom = float(stats.binomtest(n_agree, n_total, 0.5, alternative="greater").pvalue)
            z_binom = float((n_agree - n_total * 0.5) / np.sqrt(n_total * 0.25))
        else:
            p_binom = 1.0
            z_binom = 0.0

        # Rank correlation: delay rank vs (reverse) magnification rank
        # For images, compute mean delay relative to reference (first image)
        ref_img = images[0]
        delays_rel = {ref_img: 0.0}
        for pair_key, pair_data in pairs.items():
            imgs = pair_key.split("-")
            if len(imgs) != 2:
                continue
            img_i, img_j = imgs[0], imgs[1]
            if img_i not in mu_norm or img_j not in mu_norm:
                continue
            delay_info = pair_data.get("broadband", {})
            if "delay_days" not in delay_info:
                continue
            dt_ij = float(delay_info["delay_days"])
            # Build relative delays (approximate)
            if img_i == ref_img and img_j not in delays_rel:
                delays_rel[img_j] = dt_ij
            elif img_j == ref_img and img_i not in delays_rel:
                delays_rel[img_i] = -dt_ij

        # If we have delays for all images, compute rank correlation
        tau = None
        tau_p = None
        common_imgs = [img for img in images if img in delays_rel and img in mu_norm]
        if len(common_imgs) >= 3:
            delay_vals = [delays_rel[img] for img in common_imgs]
            mu_vals = [mu_norm[img] for img in common_imgs]
            # TEP predicts negative correlation (higher mu -> shorter delay)
            try:
                tau, tau_p = stats.kendalltau(delay_vals, mu_vals)
                tau = float(tau)
                tau_p = float(tau_p)
            except Exception:
                tau = None
                tau_p = None

        system_results.append({
            "system": sys_name,
            "n_images": len(images),
            "n_pairs": n_total,
            "n_agree": n_agree,
            "n_disagree": n_disagree,
            "p_binomial": p_binom,
            "z_binomial": z_binom,
            "kendall_tau": tau,
            "kendall_p": tau_p,
            "pairs": pair_results,
            "mean_fluxes": {img: float(mean_fluxes[img]) for img in mean_fluxes},
            "mu_norm": {img: float(mu_norm[img]) for img in mu_norm},
        })

    # ------------------------------------------------------------------
    # Cross-system combination
    # ------------------------------------------------------------------
    n_systems = len(system_results)
    total_pairs = sum(r["n_pairs"] for r in system_results)
    total_agree = sum(r["n_agree"] for r in system_results)
    total_disagree = sum(r["n_disagree"] for r in system_results)

    print_status(f"\nProcessed {n_systems} CosmoGRAIL systems")
    print_status(f"  Total valid pairs: {total_pairs}")
    print_status(f"  Agree with TEP: {total_agree}")
    print_status(f"  Disagree with TEP: {total_disagree}")
    print_status(f"  Skipped: {len(skipped_systems)}")

    # Binomial across all pairs (p=0.5 null is approximate; see caveat above)
    if total_pairs > 0:
        p_all = float(stats.binomtest(total_agree, total_pairs, 0.5, alternative="greater").pvalue)
        z_all = float((total_agree - total_pairs * 0.5) / np.sqrt(total_pairs * 0.25))
    else:
        p_all = 1.0
        z_all = 0.0

    print_status(f"\nCross-system directional agreement (all pairs):")
    print_status(f"  {total_agree}/{total_pairs} pairs show brighter image arriving earlier")
    print_status(f"  Agreement fraction: {total_agree/total_pairs:.1%}")
    print_status(f"  Approximate binomial (p=0.5 null): p = {p_all:.4f}, z = {z_all:+.2f}")
    print_status(f"  [CAVEAT] p=0.5 is approximate; lens geometry can bias the true null.")

    # Stouffer combination of per-system z-scores
    z_scores = [r["z_binomial"] for r in system_results if r["n_pairs"] > 0]
    if z_scores:
        z_stouffer = float(np.sum(z_scores) / np.sqrt(len(z_scores)))
        p_stouffer = float(stats.norm.sf(z_stouffer))
    else:
        z_stouffer = 0.0
        p_stouffer = 1.0

    print_status(f"\nStouffer combination ({len(z_scores)} systems):")
    print_status(f"  z = {z_stouffer:+.2f}, p = {p_stouffer:.4f}")

    # Kendall tau summary (reference-dependent; not a robust evidence metric)
    tau_vals = [r["kendall_tau"] for r in system_results if r["kendall_tau"] is not None]
    if tau_vals:
        mean_tau = float(np.mean(tau_vals))
        # Do not compute p-value: tau is strongly reference-dependent and
        # the null distribution is unknown due to geometric confounding.
    else:
        mean_tau = None

    if mean_tau is not None:
        print_status(f"\nMean Kendall tau (delay vs magnification, reference-dependent):")
        print_status(f"  tau = {mean_tau:+.3f}")
        print_status(f"  [NOTE] Tau sign depends on choice of reference image.")
        print_status(f"  [NOTE] Lens geometry dominates; no valid null hypothesis without lens models.")

    # ------------------------------------------------------------------
    # Save outputs
    # ------------------------------------------------------------------
    out = {
        "step": STEP_NUM,
        "status": "success",
        "description": (
            "Cross-system directional consistency check using CosmoGRAIL quasar lens "
            "systems. Tests whether brighter images systematically arrive earlier, as "
            "predicted by the TEP proxy model with alpha_proxy < 0."
        ),
        "caveats": [
            "This is a consistency check, not an independent evidence strand.",
            "Lens geometry can confound the delay-magnification relationship; the blind-prediction residual test (Step 07) is the only pre-specified primary test because it compares observed delays to lens-model-predicted geometric delays.",
            "The binomial p=0.5 null is approximate. Lens geometry can create systematic correlations between magnification and delay (e.g., images near the Einstein radius are both more magnified and have shorter delays). The true null is lens-model-dependent and generally unknown without blind predictions.",
            "Quasar fluxes include microlensing and intrinsic variability; mean fluxes are imperfect magnification proxies.",
            "Systems are quasars, not supernovae; different astrophysical systematics apply.",
            f"alpha_proxy = {ALPHA_PROXY} was calibrated on SN Refsdal and applied here without refitting.",
            "Kendall tau is reference-dependent and its sign depends on the choice of reference image; it is not reported as an evidence metric.",
        ],
        "n_systems_processed": n_systems,
        "n_systems_skipped": len(skipped_systems),
        "skipped_systems": skipped_systems,
        "total_pairs": total_pairs,
        "total_agree": total_agree,
        "total_disagree": total_disagree,
        "cross_system_binomial": {
            "p_value": p_all,
            "z_score": z_all,
            "n_agree": total_agree,
            "n_total": total_pairs,
        },
        "stouffer_combination": {
            "z_score": z_stouffer,
            "p_value": p_stouffer,
            "n_systems": len(z_scores),
        },
        "kendall_summary": {
            "mean_tau": mean_tau,
            "n_systems_with_tau": len(tau_vals) if tau_vals else 0,
            "note": "Tau is reference-dependent; no valid null hypothesis without lens models.",
        },
        "per_system_results": system_results,
        "all_pair_results": all_pair_results,
    }

    out_path = PROJECT_ROOT / "results" / "outputs" / f"step_{STEP_NUM}_cosmograil_cross_system.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=safe_json_default)
    print_status(f"\nResults saved to {out_path}")

    # ------------------------------------------------------------------
    # Figure: per-system agreement fraction
    # ------------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from scripts.utils.plot_style import COLORS

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # Left panel: agreement fraction per system
        ax1 = axes[0]
        names = [r["system"] for r in system_results]
        agree_frac = [r["n_agree"] / r["n_pairs"] if r["n_pairs"] > 0 else 0 for r in system_results]
        colors = [COLORS['tep'] if f > 0.5 else COLORS['null'] for f in agree_frac]
        ax1.barh(range(len(names)), agree_frac, color=colors, alpha=0.7)
        ax1.axvline(0.5, color="black", linestyle="--", label="GR null (50%)")
        ax1.set_yticks(range(len(names)))
        ax1.set_yticklabels(names, fontsize=7)
        ax1.set_xlabel("Fraction of pairs agreeing with TEP prediction")
        ax1.set_xlim(0, 1)
        ax1.set_title("Per-system directional agreement")
        ax1.legend(loc="lower right")

        # Right panel: Kendall tau distribution
        ax2 = axes[1]
        if tau_vals:
            ax2.hist(tau_vals, bins=np.arange(-1, 1.1, 0.2), edgecolor="black", alpha=0.7)
            ax2.axvline(0.0, color="black", linestyle="--", label="GR null (tau=0)")
            ax2.axvline(mean_tau, color="C0", linestyle="-", linewidth=2, label=f"Mean tau = {mean_tau:+.3f}")
            ax2.set_xlabel("Kendall tau (delay vs magnification)")
            ax2.set_ylabel("Number of systems")
            ax2.set_title("Cross-system rank correlation")
            ax2.legend()

        fig.suptitle(
            f"CosmoGRAIL Cross-System Check: {total_agree}/{total_pairs} pairs agree "
            f"(p={p_all:.3f}, z={z_all:+.2f})\n"
            "CAVEAT: Confounded by lens geometry; not an independent evidence strand",
            fontsize=10,
        )
        fig.tight_layout()

        fig_path = (
            PROJECT_ROOT / "results" / "figures" / f"step_{STEP_NUM}_cross_system.png"
        )
        fig.savefig(fig_path, dpi=300, bbox_inches="tight")
        print_status(f"Figure saved to {fig_path}")
        plt.close(fig)
    except Exception as e:
        print_status(f"Figure generation skipped: {e}", "WARNING")

    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
