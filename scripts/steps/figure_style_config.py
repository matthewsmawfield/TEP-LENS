"""
Unified figure styling configuration for TEP-LENS manuscript.
Matches the manuscript's warm/cosmic color palette.
"""

import matplotlib.pyplot as plt
import matplotlib as mpl

# Color palette matching manuscript CSS
COLORS = {
    # Primary manuscript colors
    'purple_deep': '#301E30',    # Deep Purple/Black
    'purple_dark': '#551F34',    # Dark Red/Purple
    'purple_mid': '#4B3A55',     # Muted Purple
    'purple_accent': '#A14856',  # Red/Pink Accent
    'gray_purple': '#6A5A70',    # Muted Purple/Gray
    
    # Data visualization colors
    'gr_null': '#6A5A70',        # Gray-purple for GR/null hypothesis
    'tep_prediction': '#A14856', # Red/pink accent for TEP predictions
    'observed': '#551F34',       # Dark red/purple for observations
    'model_ensemble': '#4B3A55', # Muted purple for model ensemble
    
    # Supporting colors
    'positive': '#A14856',       # Red/pink for positive evidence
    'neutral': '#6A5A70',        # Gray-purple for neutral/reference
    'background': '#F9F7F8',     # Warm white background
    'text': '#201520',           # Warm black text
}

# Figure styling parameters
FIGURE_CONFIG = {
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'axes.edgecolor': COLORS['gray_purple'],
    'axes.labelcolor': COLORS['text'],
    'axes.linewidth': 1.2,
    'axes.grid': True,
    'grid.alpha': 0.15,
    'grid.color': COLORS['gray_purple'],
    'grid.linestyle': '-',
    'grid.linewidth': 0.5,
    'xtick.color': COLORS['text'],
    'ytick.color': COLORS['text'],
    'text.color': COLORS['text'],
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'DejaVu Serif'],
    'font.size': 11,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 13,
    'lines.linewidth': 1.8,
    'lines.markersize': 8,
}

def apply_style():
    """Apply the unified TEP-LENS figure style."""
    mpl.rcParams.update(FIGURE_CONFIG)

def get_color_palette():
    """Return the standard color palette for multi-series plots."""
    return [
        COLORS['tep_prediction'],  # Primary: TEP/positive
        COLORS['gr_null'],         # Secondary: GR/null
        COLORS['observed'],        # Tertiary: Observations
        COLORS['model_ensemble'],  # Quaternary: Models
        COLORS['purple_accent'],   # Accent
    ]

def save_figure(fig, filename, dpi=300, bbox_inches='tight'):
    """
    Save figure with consistent settings.
    
    Parameters:
    -----------
    fig : matplotlib.figure.Figure
        Figure to save
    filename : str
        Output filename (should include path)
    dpi : int
        Resolution (default 300 for publication quality)
    bbox_inches : str
        Bounding box setting (default 'tight')
    """
    fig.savefig(filename, dpi=dpi, bbox_inches=bbox_inches, 
                facecolor='white', edgecolor='none')
    print(f"✓ Saved: {filename}")
