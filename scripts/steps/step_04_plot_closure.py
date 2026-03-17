#!/usr/bin/env python3
"""
TEP-LENS: Step 04 - Plot Route Closure Residuals (SN Refsdal)

Generates two publication-quality figures:
  Figure 1: Per-loop comparison of GR null (0) vs TEP predicted residual,
            with 1-sigma measurement sensitivity shown as error shading.
  Figure 2: TEP residual magnitude vs loop time-baseline, showing the
            SX loops dominate due to the 376-day baseline to image SX.
"""

import json
import sys
from pathlib import Path
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status
from scripts.utils.plot_style import set_pub_style, COLORS, FIG_SIZE

STEP_NUM = "04"

# Readable loop labels for figures
LOOP_LABELS = {
    "S1_S2_S3": "S1–S2–S3",
    "S1_S2_S4": "S1–S2–S4",
    "S1_S3_S4": "S1–S3–S4",
    "S1_S2_SX": "S1–S2–SX",
    "S1_S4_SX": "S1–S4–SX",
}

def main():
    print_status(f"STEP {STEP_NUM}: Generating Route Closure Figures — SN Refsdal", "TITLE")

    tep_path = PROJECT_ROOT / "results" / "outputs" / "step_03_tep_closure.json"
    if not tep_path.exists():
        print_status(f"Missing {tep_path}. Run step_03 first.", "ERROR")
    with open(tep_path, 'r') as f:
        tep_data = json.load(f)

    alpha = tep_data["alpha_tep"]
    loops = tep_data["tep_closure_loops"]

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print_status("matplotlib not installed, skipping plot generation.", "WARNING")
        return

    set_pub_style()

    out_dir = PROJECT_ROOT / "results" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Figure 1: Per-loop GR vs TEP bar chart with measurement sensitivity
    # ------------------------------------------------------------------
    loop_names = list(loops.keys())
    labels = [LOOP_LABELS[n] for n in loop_names]
    R_tep = [loops[n]["R_tep_days"] for n in loop_names]
    R_err = [loops[n]["R_tep_err_days"] for n in loop_names]
    snr = [loops[n]["snr"] for n in loop_names]

    x = np.arange(len(loop_names))
    width = 0.35

    fig, ax = plt.subplots(figsize=FIG_SIZE)

    # GR null bars (always 0, shown as baseline)
    bars_gr = ax.bar(x - width / 2, [0.0] * len(x), width,
                     color=COLORS['gray'], edgecolor='black', linewidth=1.2,
                     label=r'GR ($\Lambda$CDM): $\mathcal{R}=0$', zorder=3)

    # TEP predicted residual bars
    bars_tep = ax.bar(x + width / 2, R_tep, width,
                      color=COLORS['accent'], edgecolor='black', linewidth=1.2,
                      label=rf'TEP ($\alpha={alpha}$): $\mathcal{{R}}_{{TEP}}$', zorder=3)

    # Measurement 1-sigma sensitivity as horizontal error spans on GR side
    for xi, (err, label) in enumerate(zip(R_err, labels)):
        ax.errorbar(xi - width / 2, 0.0, yerr=err,
                    fmt='none', color='black', capsize=5, capthick=1.5,
                    linewidth=1.5, zorder=4, label='1$\\sigma$ sensitivity' if xi == 0 else '')

    # Annotate SNR on TEP bars
    for xi, (r, s) in enumerate(zip(R_tep, snr)):
        va = 'bottom' if r >= 0 else 'top'
        offset = 0.005 if r >= 0 else -0.005
        ax.text(xi + width / 2, r + offset,
                f'SNR={s:.2f}', ha='center', va=va, ,
                color='black', fontweight='bold')

    ax.axhline(0, color='black', linewidth=0.8, linestyle='-', zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, )
    ax.set_ylabel(r'Closure Residual $\mathcal{R}_{\rm closure}$ [days]', )
    ax.set_title(
        'SN Refsdal: Route-Closure Test — GR Null vs TEP Prediction\n'
        r'(Kelly et al. 2023, $\alpha_{\rm TEP}=0.05$)',
        , pad=12
    )
    ax.legend(loc='upper left')
    ax.set_xlabel('Closure Loop (image triplet)', )

    # fig.tight_layout()
    out_path1 = out_dir / f"step_{STEP_NUM}_closure_residual.png"
    fig.savefig(out_path1)
    plt.close(fig)
    print_status(f"Figure 1 saved to {out_path1}")

    # ------------------------------------------------------------------
    # Figure 2: TEP residual magnitude vs loop time-baseline (|dt_SX_S1|)
    # ------------------------------------------------------------------
    baselines = [
        max(abs(loops[n]["dt_ij_days"]),
            abs(loops[n]["dt_jk_days"]),
            abs(loops[n]["dt_ki_days"]))
        for n in loop_names
    ]

    fig2, ax2 = plt.subplots(figsize=FIG_SIZE)
    sc = ax2.scatter(baselines, [abs(r) for r in R_tep], s=120,
                     c=snr, cmap='plasma', edgecolors='black', linewidth=1.2,
                     zorder=4)
    cbar = plt.colorbar(sc, ax=ax2)
    cbar.set_label('Predicted SNR', )

    for xi, (b, r, lab) in enumerate(zip(baselines, R_tep, labels)):
        ax2.annotate(lab, (b, abs(r)), textcoords='offset points',
                     xytext=(6, 4), )

    ax2.set_xlabel('Maximum pairwise delay in loop [days]', )
    ax2.set_ylabel(r'$|\mathcal{R}_{\rm TEP}|$ [days]', )
    ax2.set_title(
        r'TEP Residual Magnitude vs Loop Time-Baseline (SN Refsdal, $\alpha=-0.05$)',
        , pad=10
    )
    # fig.tight_layout()
    out_path2 = out_dir / f"step_{STEP_NUM}_baseline_vs_residual.png"
    fig2.savefig(out_path2)
    plt.close(fig2)
    print_status(f"Figure 2 saved to {out_path2}")

    print_status(f"Step {STEP_NUM} complete.")

if __name__ == "__main__":
    main()
