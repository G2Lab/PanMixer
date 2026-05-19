import numpy as np
import pickle
import os
import matplotlib.pyplot as plt
import pandas as pd
import json
from tools.utils import load_data_all
from matplotlib.ticker import ScalarFormatter
from my_color_palette import apply_style, SECONDARY_COLORS, FONT

# ── Experiment numbers — update to match your data ───────────────────────────
SCALING_EXPERIMENT = 32   # multi-target scaling experiment (all 22 chromosomes)

HAVE_MARKER = True

from constants import (
    PLOT_OUT_PATH,
    EXPERIMENT_PATH
)

read_subjects = {
    "HG00138": "EUR",
    "HG00635": "EAS",
    "HG01112": "AMR",
    "HG01600": "EAS",
    "HG02698": "SAS",
    "NA12778": "EUR",
    "NA18853": "AFR",
}

def _line_style_cycle():
    # Enough distinct styles for multiple series
    return ["-", "--", "-.", ":", (0, (5, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1))]

def plot_scaling():
      # Downstream analysis on utility of individuals in the pangenome
    dfs = [
        pd.read_csv(f"{EXPERIMENT_PATH}/exp_{SCALING_EXPERIMENT}/data_chr{i}.csv") for i in range(1, 23)
    ]
    ld_sums = np.sum([dfs[i]["ld_sums"] for i in range(22)], axis=0)
    ld_counts = np.sum([dfs[i]["ld_counts"] for i in range(22)], axis=0)
    ld_loss = ld_sums / ld_counts

    wd_af = np.mean([dfs[i]["wd_af"] for i in range(22)], axis=0)
    wd_af_snp_only = np.mean([dfs[i]["wd_af_snp_only"] for i in range(22)], axis=0)
    wd_af_complex = np.mean([dfs[i]["wd_af_complex"] for i in range(22)], axis=0)
    wd_af_snps_0_01 = np.mean([dfs[i]["wd_af_snps_0_01"] for i in range(22)], axis=0)
    wd_af_snps_0_05 = np.mean([dfs[i]["wd_af_snps_0_05"] for i in range(22)], axis=0)
    wd_af_snps_0_1 = np.mean([dfs[i]["wd_af_snps_0_1"] for i in range(22)], axis=0)
    wd_af_snps_0_5 = np.mean([dfs[i]["wd_af_snps_0_5"] for i in range(22)], axis=0)

    data_df = pd.DataFrame({
        "num_targets": dfs[0]["num_targets"],
        "wd_af": wd_af,
        "wd_af_snp_only": wd_af_snp_only,
        "wd_af_complex": wd_af_complex,
        "wd_af_snps_0_01": wd_af_snps_0_01,
        "wd_af_snps_0_05": wd_af_snps_0_05,
        "wd_af_snps_0_1": wd_af_snps_0_1,
        "wd_af_snps_0_5": wd_af_snps_0_5,
        "ld_loss": ld_loss,
    })

    panmixer_chr21_df = pd.read_csv(f"{EXPERIMENT_PATH}/exp_{SCALING_EXPERIMENT}/data_chr21.csv")

    subplots_top_row = ["All", "SNPs Only", "LD"]
    subplot_columns_top_row = ["wd_af", "wd_af_snp_only", "ld_loss"]
    y_labels_top_row = ["Wasserstein divergence", "Wasserstein divergence", "LD Loss"]

    subplots_bottom_row = ["Perfect", "Gapless", "MAPQ60"]
    subplot_columns_bottom_row = ["perfect", "gapless", "mapq60"]

    subplots = subplots_top_row + subplots_bottom_row
    subplot_columns = subplot_columns_top_row + subplot_columns_bottom_row

    METRICS = {
        "perfect": ("total_perfect", "total_aligned"),
        "gapless": ("total_gapless_softclips_allowed", "total_aligned"),
        "mapq60":  ("mapping_quality_max_60_reads",        "total_aligned"),
    }

    panmixer_chr21_df["perfect"] = panmixer_chr21_df["HG00138_total_perfect"] / panmixer_chr21_df["HG00138_total_aligned"]
    panmixer_chr21_df["gapless"] = panmixer_chr21_df["HG00138_total_gapless_softclips_allowed"] / panmixer_chr21_df["HG00138_total_aligned"]
    panmixer_chr21_df["mapq60"] = panmixer_chr21_df["HG00138_mapping_quality_max_60_reads"] / panmixer_chr21_df["HG00138_total_aligned"]

    panmixer = {
        "perfect": 100 * panmixer_chr21_df["perfect"],
        "gapless": 100 * panmixer_chr21_df["gapless"],
        "mapq60": 100 * panmixer_chr21_df["mapq60"],
    }



    #for subplot_column in subplot_columns:
    #    if subplot_column not in data_df.columns:
    #        print(f"Column {subplot_column} not found in data")
    #        return

    series = ["All", "EAS", "AMR", "SAS", "AFR"]
    series_full = ["All", "East Asian Ancestry", "American Ancestry", "South Asian Ancestry", "African Ancestry"]

    fig, axes = plt.subplots(2, 3, figsize=(10, 8))

    line_markers = [
        ("*", 0),   # circle every 13 points
        ("s", 8),   # square every 17 points
        ("o", 16),   # triangle up every 19 points
        ("D", 24),   # triangle down every 23 points
        ("D", 32),   # diamond every 29 points
        ("P", 50),   # plus-filled every 31 points
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

    x = data_df["num_targets"]

    for i, (subplot, y_col) in enumerate(zip(subplots, subplot_columns)):
        ax = axes[i // 3][i % 3]
        
        # 1. Determine Y data based on row index
        if i < 3:
            y = data_df[y_col]
        else:
            y = panmixer[y_col]

        # 2. Create a temporary DataFrame to handle the grouping
        #    This aligns the x (with duplicates) with the corresponding y
        temp_df = pd.DataFrame({'x': x, 'y': y})
        
        # 3. Group by 'x' and calculate mean, min, and max
        #    sort_index() ensures the line is drawn from left to right correctly
        grouped = temp_df.groupby('x')['y'].agg(['mean', 'min', 'max']).sort_index()

        # 4. Plot the Mean Line
        ax.plot(grouped.index, grouped['mean'], label='Mean', linewidth=2)

        # 5. Plot the Shaded Region (Min to Max)
        ax.fill_between(
            grouped.index, 
            grouped['min'], 
            grouped['max'], 
            alpha=0.3,       # Transparency of shading
            label='Range'    # Optional label
        )

        ax.set_title(subplot)
        ax.set_xlabel("Number of Target Individuals")
        
        if i < 3:
            ax.set_ylabel(y_labels_top_row[i])
        if i == 3:
            ax.set_ylabel("Percent of Reads")

    plt.tight_layout(rect=[0, 0.1, 1, 0.95])
    plt.savefig(f"{PLOT_OUT_PATH}/RECOMB_supplementary_scaling.pdf")
    plt.savefig(f"{PLOT_OUT_PATH}/RECOMB_supplementary_scaling.png", dpi=300)
    plt.close()

    
if __name__ == "__main__":
    apply_style()
    plot_scaling()
