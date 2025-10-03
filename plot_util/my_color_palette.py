import matplotlib as mpl
from matplotlib import font_manager

# Define main colors
# Define improved main colors with better contrast
MAIN_COLORS = {
    1: "#5D3A9B",  # Deeper purple for more contrast
    2: "#8B4FAD",  # More saturated lavender
    3: "#E63976",  # Stronger pink
    4: "#FFA500",  # More vibrant orange
}

# Define additional color palettes
SECONDARY_COLORS = {
    1: "#E41A1C",  # strong red
    2: "#377EB8",  # vivid blue
    3: "#4DAF4A",  # green
    4: "#984EA3",  # purple
    5: "#FF7F00",  # orange
    6: "#A65628",  # brown
    7: "#F781BF",  # pink
    8: "#999999",  # gray (neutral, good contrast line)
}
GOOD_COLOR = "#00A79D"
BAD_COLOR = "#EC008C"

# Define font settings
DEFAULT_FONT = "Arial"  # Change to preferred font

LABEL_SIZE = 15
TICK_SIZE = 15
LEGEND_SIZE = 15

def apply_style():
    params = {
        'axes.labelsize': LABEL_SIZE,
        'axes.titlesize': LABEL_SIZE,
        'xtick.labelsize': TICK_SIZE,
        'ytick.labelsize': TICK_SIZE,
        'legend.fontsize': LEGEND_SIZE,
        'legend.title_fontsize': LEGEND_SIZE,
        'figure.titlesize': LABEL_SIZE,
        'figure.titleweight': 'bold',
        'text.usetex': True,  # Use LaTeX rendering
        'font.family': 'sans-serif',  # Ensures sans-serif is used
        'font.sans-serif': ['Ariel'],  # Correct place for Helvetica
        'savefig.bbox': 'tight',
        'savefig.format': 'pdf',
        'savefig.dpi': 300,
    }

    mpl.rcParams.update(params)

FONT = font_manager.FontProperties(family='Ariel', size=LABEL_SIZE)