import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
from tools.utils import load_data_all
from my_color_palette import apply_style, SECONDARY_COLORS, FONT
from plot_helpers import add_marker, find_when_gapscore_is_neg

from constants import PLOT_OUT_PATH

# ── Experiment numbers — update to match your data ───────────────────────────
ONE_PRIVACY_BASELINE = 5    # experiment where the individual is removed from the pangenome (baseline)
MAIN_EXPERIMENT      = 15   # main PanMixer obfuscation experiment
HAVE_MARKER          = True

def plot_M4():
      # Downstream analysis on utility of individuals in the pangenome
    data_df = load_data_all(MAIN_EXPERIMENT)

    data_df["ld"] = data_df["ld_sums"] / data_df["ld_counts"]

    data_df["privacy_loss"] = 1 - data_df["pmi_gain_normalized"]

    lowest_gapscore_neg = find_when_gapscore_is_neg(data_df, "genotypes_score")
    data_df["utility_loss"] = data_df["utility_loss"] / data_df["true_max_utility_loss"]

    privacy_baseline = load_data_all(ONE_PRIVACY_BASELINE)
    privacy_baseline["ld"] = privacy_baseline["ld_sums"] / privacy_baseline["ld_counts"]

    subplots = ["All", "SNPs Only",  "SNPs (MAF $<$ 0.05)", "LD"]
    subplot_columns = ["wd_af", "wd_af_snp_only", "wd_af_snps_0_05", "ld"]
    af_indices = {0, 1, 2}  # AF plots that get a baseline subplot above

    #for subplot_column in subplot_columns:
    #    if subplot_column not in data_df.columns:
    #        print(f"Column {subplot_column} not found in data")
    #        return

    series = ["All", "East Asian Ancestry", "American Ancestry", "South Asian Ancestry", "African Ancestry"]

    fig = plt.figure(figsize=(20, 8))
    gs = gridspec.GridSpec(2, len(subplots), figure=fig, height_ratios=[1, 6], hspace=0.08)
    top_axes = [None] * len(subplots)
    axes = []
    for i in range(len(subplots)):
        if i in af_indices:
            ax_top = fig.add_subplot(gs[0, i])
            top_axes[i] = ax_top
            ax_main = fig.add_subplot(gs[1, i], sharex=ax_top)
        else:
            ax_main = fig.add_subplot(gs[1, i])
        axes.append(ax_main)

    line_markers = [
        ("o", 0),   # circle every 13 points
        ("s", 4),   # square every 17 points
        ("*", 8),   # triangle up every 19 points
        ("D", 12),   # triangle down every 23 points
        ("D", 16),   # diamond every 29 points
        ("P", 31),   # plus-filled every 31 points
        ("X", 37),   # X-filled every 37 points
        ("*", 41),   # star every 41 points
        ("+", 43),   # plus sign every 43 points
        ("x", 47),   # cross every 47 points
        ("<", 53),   # triangle left every 53 points
        (">", 59),   # triangle right every 59 points
        ("1", 61),   # tri-down tick
        ("2", 67),   # tri-up tick
        ("3", 71),   # tri-left tick
        ("4", 73),   # tri-right tick
    ]

    legend_handles = []

    for i, (subplot, y_col) in enumerate(zip(subplots, subplot_columns)):
        ax = axes[i]
        ax_top = top_axes[i]

        min_x = float("inf")
        max_x = 0

        for j, serie in enumerate(series):
            if serie != "All":
                data_df_subset = data_df[data_df["population"] == serie].select_dtypes(include=[np.number])
            else:
                data_df_subset = data_df.select_dtypes(include=[np.number])

            utility_mean = data_df_subset.groupby("capacity").mean()

            x_pmi = utility_mean["privacy_loss"]

            min_x = min(min_x, x_pmi.min())
            max_x = max(max_x, x_pmi.max())

            if serie == "All":
                line, = ax.plot(x_pmi, (utility_mean[y_col]), label=serie, color="black", zorder = 30)
            else:
                marker_type, spacing_every = line_markers[j-1]
                line, = ax.plot(x_pmi, (utility_mean[y_col]), label=serie, color=SECONDARY_COLORS[j + 1], zorder = 3, marker=marker_type, markevery = (spacing_every, 16), markersize=15)
            if i == 0:
                legend_handles.append(line)

            if serie == "All" and HAVE_MARKER:
                add_marker(ax, x_pmi, (utility_mean[y_col]), target_x = lowest_gapscore_neg)

        ax.set_xlabel("Privacy Loss")
        if ax_top is not None:
            ax_top.set_title(subplot)
        else:
            ax.set_title(subplot)

        # Set y-axis label only for first subplot or for LD individually
        if i == 0 or y_col == "ld":
            ax.set_ylabel("Wasserstein divergence")

        # Add baseline
        privacy_baseline = privacy_baseline.select_dtypes(include=[np.number])
        privacy_baseline_mean = privacy_baseline.groupby("capacity").mean()
        loss = privacy_baseline_mean[y_col].iloc[0]
        print(f"For {y_col} the loss is {loss}")

        if ax_top is not None:
            # Plot the removed-individual baseline in the small top subplot (broken y-axis)
            line = ax_top.hlines((loss), min_x, max_x, color="black", linestyle="--", zorder=5, label="Removed individual from pangenome")
            margin = max(loss * 0.4, 1e-6)
            ax_top.set_ylim(loss - margin, loss + margin)
            ax_top.set_yticks([loss])
            ax_top.set_yticklabels([f"{loss:.4f}"])
            ax_top.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
            ax_top.spines['bottom'].set_visible(False)
            ax.spines['top'].set_visible(False)
            # Diagonal break marks
            d = 0.012
            kwargs = dict(transform=ax_top.transAxes, color='k', clip_on=False, linewidth=1)
            ax_top.plot((-d, +d), (-3*d, +3*d), **kwargs)
            ax_top.plot((1 - d, 1 + d), (-3*d, +3*d), **kwargs)
            kwargs = dict(transform=ax.transAxes, color='k', clip_on=False, linewidth=1)
            ax.plot((-d, +d), (1 - d/2, 1 + d/2), **kwargs)
            ax.plot((1 - d, 1 + d), (1 - d/2, 1 + d/2), **kwargs)
        else:
            line = ax.hlines((loss), min_x, max_x, color="black", linestyle="--", zorder=5, label="Removed individual from pangenome")
        if i == 0:
            legend_handles.append(line)

    # Legend
    series.append("Removed individual from pangenome")
    if HAVE_MARKER:
        legend_handles.append(plt.Line2D([0], [0], color="blue", marker='o', linestyle='None', label="Loss when linking fails"))
        series.append("Loss when Linking Fails")

    fig.legend(handles=legend_handles, labels=series, loc='lower center', ncol=len(series), bbox_to_anchor=(0.5, 0))

    # Flip x-axis
    for ax in axes[1:3]:
        ax.sharey   (axes[0])
    
    for ax in axes[0:3]:
        ax.set_ylim(0, 0.01)
    axes[3].set_ylim(0,0.03)

    axes[3].set_ylabel("Euclidean Distance Between LD Decay Vectors")

    for ax in axes:
        ax.set_xlim(1, 0.0001)
        ax.set_xscale("log")
    for ax_top in top_axes:
        if ax_top is not None:
            ax_top.set_xscale("log")
            ax_top.set_xlim(1, 0.0001)

    plt.tight_layout(rect=[0, 0.1, 1, 0.95])
    plt.savefig(f"{PLOT_OUT_PATH}/M4_v2_downstream_analysis.pdf")
    plt.savefig(f"{PLOT_OUT_PATH}/M4_v2_downstream_analysis.png", dpi=300)
    plt.close()

    
if __name__ == "__main__":
    apply_style()
    plot_M4()
