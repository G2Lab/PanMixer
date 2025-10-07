import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd

from tools.utils import load_data_all
from my_color_palette import apply_style, MAIN_COLORS, GOOD_COLOR, BAD_COLOR
from constants import PLOT_OUT_PATH

GAPSCORE_RESULTS = 2

def plot_privacy_loss_two_boxes_per_metric():
    # Load once
    df = load_data_all(GAPSCORE_RESULTS).copy()

    # --- Define privacy loss on y-axis ---
    df["privacy_loss"] = 1 - (df["pmi_gain"] / df["max_pmi_gain"])

    # remove subject in 1000 genomes dataset
    df = df[df["subject"] != "HG02145"]

    # We'll split by linkability using gap_score sign
    if "gap_score" not in df.columns:
        raise ValueError("Expected column 'gap_score' not found.")

    # Metrics to display (two boxes per metric)
    metrics = [
        "genotypes_score",
        "haplotype_score",
        "haplotype_both_score",
        "no_weight_score",
    ]
    # Keep only metrics present
    metrics = [m for m in metrics if m in df.columns]
    if not metrics:
        raise ValueError("None of the requested metric columns were found in the dataframe.")

    # Prepare box data
    box_data = []      # [cannot_link_metric1, can_link_metric1, cannot_link_metric2, can_link_metric2, ...]
    positions = []     # x positions
    labels = ["Genotypes", "Haplotypes", "Haplotypes (both)", "No weight"]
    labels = labels[:len(metrics)]

    group_spacing = 1.5   # space between metric groups
    within_offset = 0.25  # half-distance between the two boxes in a group
    box_width = 0.35

    # Build data per metric
    for i, metric in enumerate(metrics):
        # rows valid for this metric + have privacy_loss and gap_score
        dfi = df.dropna(subset=[metric, "privacy_loss", metric]).copy()
        # define linkability split by gap_score sign
        cannot_link = dfi[dfi[metric] < 0]["privacy_loss"].values
        can_link    = dfi[dfi[metric] > 0]["privacy_loss"].values

        # Append in fixed order so colors/legend are consistent
        base = i * group_spacing + 1  # base x for group i (1-indexed for nicer tick positions)
        positions += [base - within_offset, base + within_offset]
        box_data += [cannot_link, can_link]

    # Plot
    fig, ax = plt.subplots(figsize=(11, 6))
    bp = ax.boxplot(
        box_data,
        positions=positions,
        patch_artist=True,
        widths=box_width,
        showfliers=False,
    )

    # Color: cannot_link (GOOD_COLOR), can_link (BAD_COLOR) alternating
    colors = []
    for _ in range(len(metrics)):
        colors += [GOOD_COLOR, BAD_COLOR]
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)

    # X ticks at group centers
    centers = [1 + i * group_spacing for i in range(len(metrics))]
    ax.set_xticks(centers, labels, rotation=0, ha="center")

    # Y axis: privacy loss in [0, 1]
    ax.set_ylabel("Privacy Risk")
    ax.set_yscale("log")   # <-- make y-axis logarithmic

       # Legend
    legend_handles = [
        plt.Line2D([0], [0], color=GOOD_COLOR, lw=6, label="gap score $<$ 0 (Cannot Link)"),
        plt.Line2D([0], [0], color=BAD_COLOR, lw=6, label="gap score $>$ 0 (Can Link)"),
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),   # push legend below the plot
        frameon=False,
        ncol=2,                        # two columns
    )

    plt.tight_layout()
    os.makedirs(PLOT_OUT_PATH, exist_ok=True)
    outpath = f"{PLOT_OUT_PATH}/Supplementary_PrivacyLoss_TwoBoxesPerMetric.pdf"
    plt.savefig(outpath)
    plt.close()
    print(f"Saved: {outpath}")

if __name__ == "__main__":
    apply_style()
    plot_privacy_loss_two_boxes_per_metric()
