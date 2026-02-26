# Re-usable matplotlib styling configurations
import matplotlib.pyplot as plt

def set_pub_style():
    plt.style.use('default')
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'STIXGeneral', 'DejaVu Serif', 'serif'],
        'font.size': 12,
        'axes.labelsize': 14,
        'axes.titlesize': 14,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 12,
        'axes.grid': True,
        'grid.alpha': 0.3,
        'grid.linestyle': '--',
        'axes.linewidth': 1.2,
        'xtick.major.width': 1.2,
        'ytick.major.width': 1.2,
        'lines.linewidth': 2.0,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.format': 'png'
    })

COLORS = {
    "primary": "#2c3e50",
    "accent": "#e74c3c",
    "highlight": "#f39c12",
    "gray": "#7f8c8d",
    "light_gray": "#bdc3c7",
    "success": "#27ae60",
    "info": "#2980b9"
}

FIG_SIZE = (8, 6)
