import numpy as np
import pickle
import os
import matplotlib.pyplot as plt
import pandas as pd
import json
from tools.utils import load_data_all
from matplotlib.ticker import ScalarFormatter

from my_color_palette import apply_style, SECONDARY_COLORS, FONT

ONE_PRIVACY_BASELINE = 5 
MAIN_EXPERIMENT = 2

HAVE_MARKER = True

from constants import (
    PLOT_OUT_PATH,
)

def add_marker(ax, x_values, y_values, target_x=0.01):
    x_values = x_values.to_numpy()
    y_values = y_values.to_numpy()
    # Find indices of the two closest points to target_x
    closest_points = np.argsort(np.abs(x_values - target_x))[:2]
    
    # Extract the x and y values of the two closest points
    x1, x2 = x_values[closest_points]
    y1, y2 = y_values[closest_points]
    
    # Calculate the y-value at target_x using linear interpolation
    target_y = y1 + (y2 - y1) * (target_x - x1) / (x2 - x1)
    
    # Plot the marker at the target point
    ax.scatter(target_x, target_y, color="blue", zorder=5)
    ax.vlines(target_x, 0, target_y, color="blue", linestyle="--", zorder=5)
    ax.annotate(f"{target_y:.3f}", 
                (target_x, target_y), 
                textcoords="offset points", 
                xytext=(-35, 0), 
                ha='center', fontsize=15,
                zorder=5,
                fontproperties=FONT)

def find_when_gapscore_is_neg(df_all):
    unique_individuals = df_all["subject"].unique()

    lowest_privacy_loss = []

    for subject in unique_individuals:
        df_subject = df_all[df_all["subject"] == subject]
        gap_scores_negative = df_subject[df_subject["genotypes_score"] < 0]
        #print(subject, 1 - gap_scores_negative["pmi_gain_normalized"].min())

        print(len(gap_scores_negative))

        lowest_privacy_loss.append(gap_scores_negative["privacy_loss"].max())

    print(lowest_privacy_loss)
    return min(lowest_privacy_loss)

def _line_style_cycle():
    # Enough distinct styles for multiple series
    return ["-", "--", "-.", ":", (0, (5, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1))]

def plot_M4():
      # Downstream analysis on utility of individuals in the pangenome
    data_df = load_data_all(MAIN_EXPERIMENT)

    data_df["privacy_loss"] = 1 - (data_df["pmi_gain"] / data_df["max_pmi_gain"])

    lowest_gapscore_neg = find_when_gapscore_is_neg(data_df)
    print(lowest_gapscore_neg)
    

    data_df["privacy_loss"] = 1 - data_df["pmi_gain_normalized"]
    data_df["utility_loss"] = data_df["utility_loss"] / data_df["true_max_utility_loss"]
    privacy_baseline = load_data_all(ONE_PRIVACY_BASELINE)

    subplots = ["All", "SNPs Only",  "SNPs (MAF $<$ 0.05)", "LD"]
    subplot_columns = ["wd_af", "wd_af_snp_only", "wd_af_snps_0_05", "ld_euclidean"]

    #for subplot_column in subplot_columns:
    #    if subplot_column not in data_df.columns:
    #        print(f"Column {subplot_column} not found in data")
    #        return

    series = ["All", "East Asian Ancestry", "American Ancestry", "South Asian Ancestry", "African Ancestry"]

    fig, axes = plt.subplots(1, len(subplots), figsize=(20, 8), sharex=True)  # remove sharey=True

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

        min_x = float("inf")
        max_x = 0

        for j, serie in enumerate(series):
            if serie != "All":
                data_df_subset = data_df[data_df["population"] == serie].select_dtypes(include=[np.number])
            else:
                data_df_subset = data_df.select_dtypes(include=[np.number])

            utility_mean = data_df_subset.groupby("capacity").mean()
            utility_max = data_df_subset.groupby("capacity").max()
            utility_min = data_df_subset.groupby("capacity").min()

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
        ax.set_title(subplot)

        # Set y-axis label only for first subplot or for LD individually
        if i == 0 or y_col == "ld":
            ax.set_ylabel("Wasserstein divergence")

        # Add baseline
        privacy_baseline = privacy_baseline.select_dtypes(include=[np.number])
        privacy_baseline_mean = privacy_baseline.groupby("capacity").mean()
        loss = privacy_baseline_mean[y_col].iloc[0]
        print(f"For {y_col} the loss is {loss}")
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
        ax.sharey(axes[0])
    
    for ax in axes[0:3]:
        ax.set_ylim(0, 0.01)
    axes[3].set_ylim(0,0.7)

    axes[3].set_ylabel("Euclidean Distance Between LD Decay Vectors")

    for ax in axes:
        ax.set_xlim(1, 0.0001)
        ax.set_xscale("log")

    plt.tight_layout(rect=[0, 0.1, 1, 0.95])
    plt.savefig(f"{PLOT_OUT_PATH}/M4_v2_downstream_analysis.pdf")
    plt.savefig(f"{PLOT_OUT_PATH}/M4_v2_downstream_analysis.png", dpi=300)
    plt.close()

    
if __name__ == "__main__":
    apply_style()
    plot_M4()
