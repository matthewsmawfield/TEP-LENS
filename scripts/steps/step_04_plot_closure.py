#!/usr/bin/env python3
"""
TEP-LENS: Step 04 - Plot Route Closure Residuals

Generates a visualization comparing the General Relativity null closure 
expectation vs the predicted TEP geometric residual.
"""

import json
import sys
from pathlib import Path
import numpy as np

# Add project root to path for imports
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status
from scripts.utils.plot_style import set_pub_style, COLORS, FIG_SIZE

STEP_NUM = "04"

def main():
    print_status(f"STEP {STEP_NUM}: Generating Route Closure Figures", "TITLE")
    
    # Load GR results
    gr_path = PROJECT_ROOT / "results" / "outputs" / "step_02_gr_closure.json"
    if not gr_path.exists():
        print_status(f"Missing {gr_path}", "ERROR")
    with open(gr_path, 'r') as f:
        gr_data = json.load(f)["gr_predictions"]
        
    # Load TEP results
    tep_path = PROJECT_ROOT / "results" / "outputs" / "step_03_tep_closure.json"
    if not tep_path.exists():
        print_status(f"Missing {tep_path}", "ERROR")
    with open(tep_path, 'r') as f:
        tep_data = json.load(f)["tep_predictions"]
        
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print_status("matplotlib not installed, skipping plot generation.", "WARNING")
        return
        
    set_pub_style()
    
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    
    # Simple bar chart comparing the two residuals
    models = ['GR ($\\Lambda$CDM)', f'TEP ($\\alpha={tep_data["alpha_assumed"]}$)']
    residuals = [gr_data["closure_sum"], tep_data["residual_days"]]
    
    colors = [COLORS['gray'], COLORS['accent']]
    
    ax.bar(models, residuals, color=colors, width=0.5, edgecolor='black', linewidth=1.5)
    
    ax.axhline(0, color='black', linewidth=1, linestyle='-')
    ax.set_ylabel(r"Closure Residual $\mathcal{R}_{\rm closure}$ [days]", fontsize=14)
    ax.set_title("SN H0pe Triply-Imaged Route Closure Test", fontsize=15, pad=15)
    
    # Add value annotations
    for i, v in enumerate(residuals):
        v_str = f"{v:.2f}" if v != 0 else "0.00"
        offset = -0.5 if v < 0 else 0.5
        va = 'top' if v < 0 else 'bottom'
        ax.text(i, v + offset, f"{v_str} days", ha='center', va=va, fontsize=12, fontweight='bold',
                color=colors[i])
                
    # Add a note about geometric consistency
    ax.text(0.05, 0.95, r"$\mathcal{R} = \Delta t_{AB} + \Delta t_{BC} + \Delta t_{CA}$", 
            transform=ax.transAxes, fontsize=14, verticalalignment='top',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor=COLORS['gray']))
    
    fig.tight_layout()
    
    out_dir = PROJECT_ROOT / "results" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"step_{STEP_NUM}_closure_residual.png"
    
    fig.savefig(out_path)
    plt.close(fig)
    
    print_status(f"Figure saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")

if __name__ == "__main__":
    main()
