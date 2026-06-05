#!/usr/bin/env python3
"""
Generate publication-quality figures for COSMOGRAIL temporal shear analysis.
"""

import json
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.plot_style import set_pub_style, COLORS, FIG_SIZES, save_fig

set_pub_style()
STEP_NUM = "32b"


def load_results():
    """Load temporal shear results."""
    results_path = Path(__file__).parent.parent.parent / 'results' / 'outputs' / 'step_30_cosmograil_temporal_shear.json'
    with open(results_path) as f:
        return json.load(f)


def collect_pair_rows(data):
    """Collect finite temporal-shear pair summaries from step 30 output."""
    rows = []
    for sys_id, sys_data in data['systems'].items():
        for pair_id, pair_data in sys_data['pairs'].items():
            gamma = pair_data.get('gamma', {})
            value = gamma.get('value')
            if value is None or not np.isfinite(value):
                continue

            uncertainty = gamma.get('uncertainty')
            sigma = abs(value / uncertainty) if uncertainty and uncertainty > 0 else 0.0
            broadband_delay = pair_data.get('broadband', {}).get('delay_days')
            intercept = gamma.get('intercept')
            rows.append({
                'system': sys_id,
                'pair': pair_id,
                'gamma': float(value),
                'uncertainty': float(uncertainty) if uncertainty and np.isfinite(uncertainty) else np.nan,
                'sigma': float(sigma),
                'broadband_delay': float(broadband_delay) if broadband_delay is not None and np.isfinite(broadband_delay) else np.nan,
                'intercept': float(intercept) if intercept is not None and np.isfinite(intercept) else np.nan,
                'multiscale': pair_data.get('multiscale', {}),
            })
    return rows


def plot_gamma_distribution(data, output_dir):
    """Plot distribution of gamma values across all pairs."""
    fig, ax = plt.subplots(figsize=FIG_SIZES["web_standard"])

    rows = collect_pair_rows(data)
    if not rows:
        ax.text(0.5, 0.5, 'No finite temporal-shear measurements', ha='center', va='center')
        ax.axis('off')
        png_path = output_dir / f'step_{STEP_NUM}_cosmograil_gamma_significance.png'
        pdf_path = output_dir / f'step_{STEP_NUM}_cosmograil_gamma_significance.pdf'
        save_fig(fig, png_path, close=False)
        fig.savefig(pdf_path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print(f"Saved: {png_path.name} / {pdf_path.name}")
        return {"png": str(png_path), "pdf": str(pdf_path)}
    
    # Sort by sigma
    rows = sorted(rows, key=lambda row: row['sigma'], reverse=True)
    sigmas = [row['sigma'] for row in rows]
    labels = [f"{row['system']}\n{row['pair']}" for row in rows]
    colors = [
        COLORS['significant'] if row['sigma'] >= 3 else
        COLORS['marginal'] if row['sigma'] >= 2 else
        COLORS['null']
        for row in rows
    ]
    
    x = np.arange(len(rows))
    bars = ax.bar(x, sigmas, color=colors, edgecolor='black', linewidth=0.5)
    
    # Add significance thresholds
    ax.axhline(3, color=COLORS['red'], linestyle='--', linewidth=1, label='3σ threshold')
    ax.axhline(2, color=COLORS['marginal'], linestyle='--', linewidth=1, label='2σ threshold')
    
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Significance (|Γ|/σ)')
    ax.set_xlabel('Lens System / Image Pair')
    ax.set_title('COSMOGRAIL Temporal Shear Detection Significance')
    ax.legend(loc='upper right')
    ax.set_ylim(0, max(3.2, max(sigmas) * 1.1))
    
    # Annotate top 3
    for i in range(min(3, len(sigmas))):
        if sigmas[i] >= 3:
            ax.annotate(f'{sigmas[i]:.1f}σ', (x[i], sigmas[i] + 0.2), 
                       ha='center', fontsize=8, fontweight='bold')
    
    png_path = output_dir / f'step_{STEP_NUM}_cosmograil_gamma_significance.png'
    pdf_path = output_dir / f'step_{STEP_NUM}_cosmograil_gamma_significance.pdf'
    save_fig(fig, png_path, close=False)
    fig.savefig(pdf_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Saved: {png_path.name} / {pdf_path.name}")
    return {"png": str(png_path), "pdf": str(pdf_path)}


def plot_multiscale_delays(data, output_dir):
    """Plot delay vs timescale for the strongest finite measurements."""
    fig, axes = plt.subplots(2, 3, figsize=FIG_SIZES["web_quad"], sharey=False)
    axes = axes.flatten()

    tau_values = [5, 10, 20, 40, 80, 160]
    rows = [
        row for row in collect_pair_rows(data)
        if sum(
            1 for tau in tau_values
            if row['multiscale'].get(str(tau), {}).get('delay_days') is not None
            and np.isfinite(row['multiscale'].get(str(tau), {}).get('delay_days'))
        ) >= 2
    ]
    rows = sorted(rows, key=lambda row: row['sigma'], reverse=True)[:5]

    if not rows:
        axes[0].text(0.5, 0.5, 'No multiscale delay fits passed quality cuts', ha='center', va='center')
        for ax in axes:
            ax.axis('off')

    for idx, row in enumerate(rows):
        ax = axes[idx]

        delays = []
        valid_log_tau = []

        for tau in tau_values:
            ms = row['multiscale'].get(str(tau), {})
            d = ms.get('delay_days')
            if d is not None and np.isfinite(d):
                delays.append(d)
                valid_log_tau.append(np.log10(tau))
        
        if len(delays) < 2:
            continue
        
        # Plot data points
        color = COLORS['significant'] if row['sigma'] >= 3 else COLORS['marginal'] if row['sigma'] >= 2 else COLORS['null']
        ax.scatter(valid_log_tau, delays, s=60, c=color,
                  edgecolor='black', linewidth=0.5, zorder=5)
        
        # Fit line
        gamma = row['gamma']
        intercept = row['intercept']
        if np.isfinite(intercept):
            x_fit = np.linspace(min(valid_log_tau) - 0.1, max(valid_log_tau) + 0.1, 100)
            y_fit = gamma * x_fit + intercept
            ax.plot(x_fit, y_fit, color=COLORS['fit'], linewidth=2,
                    label=f'Γ = {gamma:.0f} days/decade')
        
        # Broadband delay reference
        bb_delay = row['broadband_delay']
        if bb_delay is not None and np.isfinite(bb_delay):
            ax.axhline(bb_delay, color=COLORS['gr'], linestyle=':', linewidth=1, alpha=0.7)
        
        ax.set_xlabel('log₁₀(τ) [days]')
        if idx % 3 == 0:
            ax.set_ylabel('Time Delay [days]')
        ax.set_title(f"{row['system']} {row['pair']}")
        if ax.get_legend_handles_labels()[0]:
            ax.legend(loc='best', fontsize=8)
    
    # Hide unused subplot
    for ax in axes[len(rows):]:
        ax.axis('off')
    
    png_path = output_dir / f'step_{STEP_NUM}_cosmograil_multiscale_delays.png'
    pdf_path = output_dir / f'step_{STEP_NUM}_cosmograil_multiscale_delays.pdf'
    save_fig(fig, png_path, close=False)
    fig.savefig(pdf_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Saved: {png_path.name} / {pdf_path.name}")
    return {"png": str(png_path), "pdf": str(pdf_path)}


def plot_summary_panel(data, output_dir):
    """Create a summary panel figure."""
    fig = plt.figure(figsize=FIG_SIZES["web_tall"])
    
    # Panel A: Gamma vs broadband delay
    ax1 = fig.add_subplot(2, 2, 1)
    
    rows = collect_pair_rows(data)
    gammas = [row['gamma'] for row in rows]
    bb_delays = [row['broadband_delay'] for row in rows]
    sigmas = [row['sigma'] for row in rows]
    
    colors = [COLORS['significant'] if s >= 3 else COLORS['marginal'] if s >= 2 else COLORS['null'] 
              for s in sigmas]
    
    finite_bb = np.isfinite(bb_delays) & np.isfinite(gammas)
    ax1.scatter(np.array(bb_delays)[finite_bb], np.array(gammas)[finite_bb],
                c=np.array(colors, dtype=object)[finite_bb],
                s=60, edgecolor='black', linewidth=0.5)
    ax1.axhline(0, color=COLORS['gr'], linestyle='--', linewidth=1)
    ax1.set_xlabel('Broadband Delay [days]')
    ax1.set_ylabel('Γ [days/decade]')
    ax1.set_title('A. Temporal Shear vs Broadband Delay')
    
    # Panel B: Histogram of gamma
    ax2 = fig.add_subplot(2, 2, 2)
    valid_gammas = [g for g in gammas if np.isfinite(g)]
    if valid_gammas:
        ax2.hist(valid_gammas, bins=15, color=COLORS['null'], edgecolor='black', alpha=0.7)
    ax2.axvline(0, color=COLORS['red'], linestyle='--', linewidth=2, label='GR prediction (Γ=0)')
    ax2.set_xlabel('Γ [days/decade]')
    ax2.set_ylabel('Count')
    ax2.set_title('B. Distribution of Temporal Shear')
    ax2.legend()
    
    # Panel C: Significance histogram
    ax3 = fig.add_subplot(2, 2, 3)
    if sigmas:
        ax3.hist(sigmas, bins=np.arange(0, max(8, max(sigmas) + 1), 0.5),
                color=COLORS['significant'], edgecolor='black', alpha=0.7)
    ax3.axvline(3, color=COLORS['red'], linestyle='--', linewidth=2, label='3σ threshold')
    ax3.axvline(2, color=COLORS['marginal'], linestyle='--', linewidth=2, label='2σ threshold')
    ax3.set_xlabel('Significance (|Γ|/σ)')
    ax3.set_ylabel('Count')
    ax3.set_title('C. Detection Significance Distribution')
    ax3.legend()
    
    # Panel D: Summary text
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.axis('off')
    
    n_sig_3 = sum(1 for s in sigmas if s >= 3)
    n_sig_2 = sum(1 for s in sigmas if s >= 2)
    n_total = len(sigmas)
    top_rows = sorted(rows, key=lambda row: row['sigma'], reverse=True)[:3]
    if top_rows:
        top_lines = []
        for row in top_rows:
            uncertainty = row['uncertainty']
            unc_text = f" +/- {uncertainty:.1f}" if np.isfinite(uncertainty) else ""
            top_lines.append(
                f"  - {row['system']} {row['pair']}: "
                f"Gamma = {row['gamma']:.1f}{unc_text} ({row['sigma']:.2f} sigma)"
            )
        top_text = "\n".join(top_lines)
    else:
        top_text = "  - No finite Gamma measurements"

    summary = data.get('summary', {})
    status = summary.get('interpretation', 'No interpretation available.')
    
    summary_text = f"""
COSMOGRAIL Temporal Shear Analysis
----------------------------------

Systems analyzed: {summary.get('n_systems', len(data['systems']))}
Image pairs: {n_total}

Significant detections:
  - >3 sigma: {n_sig_3} pairs
  - >2 sigma: {n_sig_2} pairs

Strongest finite measurements:
{top_text}

Interpretation:
Under GR, Gamma = 0 (constant delay).
Under TEP-GL, Gamma != 0 (scale-dependent).

Status: {status}
"""
    
    ax4.text(0.05, 0.95, summary_text, transform=ax4.transAxes, 
            fontsize=10, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor=COLORS['light_gray'], alpha=0.6, edgecolor=COLORS['gr']))
    ax4.set_title('D. Summary')
    
    png_path = output_dir / f'step_{STEP_NUM}_cosmograil_summary_panel.png'
    pdf_path = output_dir / f'step_{STEP_NUM}_cosmograil_summary_panel.pdf'
    save_fig(fig, png_path, close=False)
    fig.savefig(pdf_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Saved: {png_path.name} / {pdf_path.name}")
    return {"png": str(png_path), "pdf": str(pdf_path)}


def main():
    project_root = Path(__file__).parent.parent.parent
    output_dir = project_root / 'results' / 'figures'
    manifest_dir = project_root / 'results' / 'outputs'
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_dir.mkdir(parents=True, exist_ok=True)
    
    print("Loading COSMOGRAIL temporal shear results...")
    data = load_results()
    
    print("Generating figures...")
    figures = {
        "gamma_significance": plot_gamma_distribution(data, output_dir),
        "multiscale_delays": plot_multiscale_delays(data, output_dir),
        "summary_panel": plot_summary_panel(data, output_dir),
    }

    rows = collect_pair_rows(data)
    manifest = {
        "step": STEP_NUM,
        "status": "success",
        "description": "Publication figures for COSMOGRAIL temporal-shear diagnostics.",
        "input": str(project_root / 'results' / 'outputs' / 'step_30_cosmograil_temporal_shear.json'),
        "figures": figures,
        "summary": {
            "n_systems": data.get("summary", {}).get("n_systems"),
            "n_pairs_valid": len(rows),
            "n_significant_2sigma": sum(1 for row in rows if row["sigma"] >= 2),
            "n_significant_3sigma": sum(1 for row in rows if row["sigma"] >= 3),
            "interpretation": data.get("summary", {}).get("interpretation"),
        },
    }
    manifest_path = manifest_dir / f"step_{STEP_NUM}_temporal_shear_figures.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"Saved manifest: {manifest_path}")
    
    print("\nDone!")


if __name__ == '__main__':
    main()
