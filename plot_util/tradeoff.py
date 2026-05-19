import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

from tools.utils import load_data_all
from my_color_palette import apply_style, MAIN_COLORS, FONT

from constants import PLOT_OUT_PATH

# ── Experiment numbers — update to match your data ───────────────────────────
PRIVACY_EXPERIMENT        = 7   # fixed-privacy optimisation
UTILITY_EXPERIMENT        = 10   # fixed-utility optimisation
RECOMB_PRIVACY_EXPERIMENT = 7   # RECOMB fixed-privacy optimisation
RECOMB_UTILITY_EXPERIMENT = 10   # RECOMB fixed-utility optimisation


def _plot_tradeoff_panel(ax, aggre_df, population):
    """Plot a single privacy-utility tradeoff curve for one population."""
    subset = aggre_df if population == "All" else aggre_df[aggre_df["population"] == population]
    df = subset.select_dtypes(include=[np.number]).copy()
    df["utility_loss_norm"] = df["utility_loss"] / df["true_max_utility_loss"]
    df["pmi_gain_norm"]     = df["pmi_gain"]     / df["max_pmi_gain"]

    g_mean = df.groupby("capacity").mean()
    g_max  = df.groupby("capacity").max()
    g_min  = df.groupby("capacity").min()

    x      = g_mean["utility_loss_norm"]
    y_mean = 1 - g_mean["pmi_gain_norm"]
    y_max  = 1 - g_max["pmi_gain_norm"]
    y_min  = 1 - g_min["pmi_gain_norm"]

    max_idx = g_mean["utility_loss_norm"].idxmax()
    ref_x   = g_mean.loc[max_idx, "utility_loss_norm"]
    ref_y   = 1 - g_mean.loc[max_idx, "pmi_gain_norm"]

    ax.plot(x, y_mean, color=MAIN_COLORS[1], marker=".")
    ax.fill_between(x, y_max, y_min, color=MAIN_COLORS[1], alpha=0.2)
    ax.hlines(y=ref_y, xmin=ref_x, xmax=1, linestyles="dashed", colors=MAIN_COLORS[1])
    ax.set_title(population)
    ax.set_xlim([0, 1])
    ax.set_xlabel("Utility Loss")
    ax.set_ylabel("Privacy Risk")


def plot_figure_A(
    privacy_exp=PRIVACY_EXPERIMENT,
    utility_exp=UTILITY_EXPERIMENT,
    figsize=(15, 10),
    output_prefix="",
):
    """2x5 grid showing privacy-utility tradeoff across populations for two
    optimisation strategies (fixed utility top row, fixed privacy bottom row)."""
    print(f"Using experiment #{privacy_exp} (fixed privacy) and #{utility_exp} (fixed utility)")
    aggre_df_utility = load_data_all(utility_exp)
    aggre_df_privacy = load_data_all(privacy_exp)

    populations = [
        "All",
        "East Asian Ancestry",
        "American Ancestry",
        "South Asian Ancestry",
        "African Ancestry",
    ]

    fig, axes = plt.subplots(
        2, len(populations),
        sharex=True, sharey=True,
        figsize=figsize,
        gridspec_kw={"hspace": 0.3},
    )
    axes = np.atleast_2d(axes)

    fig.text(0.5, 0.915, "Fixed Utility ($\\eta$)",    ha="center", va="bottom", fontproperties=FONT)
    fig.text(0.5, 0.48,  "Fixed Privacy ($\\epsilon$)", ha="center", va="bottom", fontproperties=FONT)

    for col, pop in enumerate(populations):
        _plot_tradeoff_panel(axes[0, col], aggre_df_utility, pop)
        _plot_tradeoff_panel(axes[1, col], aggre_df_privacy, pop)

    plt.tight_layout()
    plt.savefig(f"{PLOT_OUT_PATH}/{output_prefix}M2A_privacy_vs_utility_all.pdf")
    plt.close()


def plot_figure_A_recomb(
    privacy_exp=RECOMB_PRIVACY_EXPERIMENT,
    utility_exp=RECOMB_UTILITY_EXPERIMENT,
):
    """Compact 1x2 tradeoff figure with zoom insets (RECOMB submission)."""
    print(f"Using experiment #{privacy_exp} (fixed privacy) and #{utility_exp} (fixed utility)")
    aggre_df_utility = load_data_all(utility_exp)
    aggre_df_privacy = load_data_all(privacy_exp)

    fig, axes = plt.subplots(
        1, 2,
        sharex=True, sharey=True,
        figsize=(5, 2),
        gridspec_kw={"hspace": 0.3},
    )

    def _plot_with_bounds(ax, x, y_mean, y_max, y_min, ref_x, ref_y):
        ax.plot(x, y_mean, color=MAIN_COLORS[1], marker=".", markersize=4, zorder=2)
        ax.fill_between(x, y_min, y_max, color=MAIN_COLORS[1], alpha=0.2, zorder=1)
        ax.plot(x, y_max, color=MAIN_COLORS[2], linewidth=0.5, alpha=0.5, zorder=2)
        ax.plot(x, y_min, color=MAIN_COLORS[2], linewidth=0.5, alpha=0.5, zorder=2)
        ax.hlines(y=ref_y, xmin=ref_x, xmax=1, linestyles="dashed", colors=MAIN_COLORS[1], zorder=3)

    datasets   = [aggre_df_utility, aggre_df_privacy]
    titles     = ["Fixed Utility ($\\eta$)", "Fixed Privacy ($\\epsilon$)"]
    zoom_spans = [0.001, 0.0001]

    for i, (df, title, span) in enumerate(zip(datasets, titles, zoom_spans)):
        ax = axes[i]
        df_num = df.select_dtypes(include=[np.number]).copy()
        df_num["utility_loss_norm"] = df_num["utility_loss"] / df_num["true_max_utility_loss"]
        df_num["pmi_gain_norm"]     = df_num["pmi_gain"]     / df_num["max_pmi_gain"]

        stats = df_num.groupby("capacity").agg(["mean", "max", "min"])
        x_raw = stats["utility_loss_norm"]["mean"]
        idx   = np.argsort(x_raw)
        x       = x_raw.iloc[idx].to_numpy()
        y_mean  = (1 - stats["pmi_gain_norm"]["mean"]).iloc[idx].to_numpy()
        y_max   = (1 - stats["pmi_gain_norm"]["min"]).iloc[idx].to_numpy()
        y_min   = (1 - stats["pmi_gain_norm"]["max"]).iloc[idx].to_numpy()

        max_util_idx = x_raw.idxmax()
        ref_x = x_raw.loc[max_util_idx]
        ref_y = (1 - stats["pmi_gain_norm"]["mean"]).loc[max_util_idx]

        _plot_with_bounds(ax, x, y_mean, y_max, y_min, ref_x, ref_y)
        ax.set_title(title)
        ax.set_xlim([0, 1])
        ax.set_xlabel("Utility Loss")
        ax.set_ylabel("Privacy Risk")

        # Zoom inset centred near x = 0.2
        ax_ins = inset_axes(ax, width="50%", height="50%", loc=1, borderpad=1)
        _plot_with_bounds(ax_ins, x, y_mean, y_max, y_min, ref_x, ref_y)
        closest = np.argmin(np.abs(x - 0.2))
        x_c, y_c = x[closest], y_mean[closest]
        zoom_factor = int(1.0 / (span * 2))
        if zoom_factor % 10 == 9:
            zoom_factor += 1
        ax_ins.set_title(f"{zoom_factor}x Zoom", fontsize=15, pad=2)
        ax_ins.set_xlim(x_c - span, x_c + span)
        ax_ins.set_ylim(y_c - span, y_c + span)
        ax_ins.tick_params(labelleft=False, labelbottom=False)
        mark_inset(ax, ax_ins, loc1=2, loc2=4, fc="none", ec="0.5", linestyle="--")

    plt.tight_layout()
    plt.savefig(f"{PLOT_OUT_PATH}/RECOMB_M2A_privacy_vs_utility_all_recomb.png")
    plt.savefig(f"{PLOT_OUT_PATH}/RECOMB_M2A_privacy_vs_utility_all_recomb.pdf")
    plt.close()


if __name__ == "__main__":
    apply_style()
    # Original paper figures
    plot_figure_A()
    # RECOMB paper figures
    plot_figure_A(RECOMB_PRIVACY_EXPERIMENT, RECOMB_UTILITY_EXPERIMENT, figsize=(15, 5), output_prefix="RECOMB_")
    plot_figure_A_recomb()
