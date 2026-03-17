# Re-usable matplotlib styling configurations
# Standardized to use blues/bluey-grey color scheme and consistent dimensions
import matplotlib.pyplot as plt
import matplotlib as mpl

def set_pub_style():
    plt.style.use('default')
    plt.rcParams.update({
        'figure.figsize': (10, 6),
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'STIXGeneral', 'DejaVu Serif', 'serif'],
        'mathtext.fontset': 'dejavuserif',
        'font.size': 11,
        'axes.labelsize': 12,
        'axes.titlesize': 13,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'axes.grid': True,
        'grid.alpha': 0.2,
        'grid.linestyle': '-',
        'grid.linewidth': 0.5,
        'grid.color': '#475569',
        'axes.linewidth': 1.2,
        'axes.edgecolor': '#475569',
        'axes.labelcolor': '#0f172a',
        'xtick.color': '#0f172a',
        'ytick.color': '#0f172a',
        'xtick.major.width': 1.2,
        'ytick.major.width': 1.2,
        'lines.linewidth': 2.0,
        'lines.markersize': 7,
        'figure.dpi': 300,
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'savefig.dpi': 300,
        'savefig.format': 'png',
        'savefig.facecolor': 'white',
        'text.color': '#0f172a'
    })
    
    # Set color cycle to blues/slate
    mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=[
        '#2563eb', '#64748b', '#60a5fa', '#1e3a8a', '#b43b4e'
    ])

# Unified color palette matching blues/bluey-grey scheme
COLORS = {
    # Legacy names (for backward compatibility)
    "primary": "#1e3a8a",      # Dark Blue
    "accent": "#2563eb",       # Royal Blue
    "highlight": "#2563eb",    # Royal Blue
    "gray": "#64748b",         # Slate Gray
    "light_gray": "#cbd5e1",   # Light Slate
    "success": "#60a5fa",      # Light Blue
    "info": "#475569",         # Dark Slate
    
    # Semantic names for clarity
    "tep": "#2563eb",          # TEP predictions (Royal Blue)
    "gr": "#64748b",           # GR/null hypothesis (Slate Gray)
    "observed": "#1e3a8a",     # Observations (Dark Blue)
    "model": "#475569",        # Model ensemble (Dark Slate)
    "text": "#0f172a",         # Text color (Very dark blue-gray)
    "background": "#ffffff",   # Background
    "red": "#b43b4e"           # Contrast accent if needed
}

FIG_SIZE = (10, 6)
