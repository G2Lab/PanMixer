import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score

from tools.utils import load_data_all
from my_color_palette import apply_style, MAIN_COLORS, GOOD_COLOR, BAD_COLOR


from constants import (
    PLOT_OUT_PATH,
    SUBPOPS,
    SUBPOPS_SMALL,
)

# ── Experiment numbers — update to match your data ───────────────────────────
GAPSCORE_RESULTS = 7   # experiment containing gap-score / linking results

def plot_combo_privacy_utility():
    data_df = load_data_all(GAPSCORE_RESULTS)
    data_df["privacy_gain"] = 1 - data_df["pmi_gain_normalized"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5, 3), sharey=False)

    conditions = {
        "PanMixer": {
            "link_column": "genotypes_score",
            "pg_column": "privacy_gain",
            "ul_column": "utility_loss_normalized",
        }
    }

    def prepare_boxplot_data(metric_col):
        box_data = []
        positions = []
        xtick_labels = []
        pos_counter = 0

        for subpop in SUBPOPS:
            if subpop == "All":
                data_df_subset = data_df
            else:
                data_df_subset = data_df[data_df["population"] == subpop]

            for cond_name, cond_info in conditions.items():
                link_col = cond_info["link_column"]
                metric = cond_info[metric_col]

                df_cond = data_df_subset.dropna(subset=[metric, link_col]).copy()
                df_cond["link_outcome"] = (df_cond[link_col] > 0).astype(int)

                can_link = df_cond[df_cond["link_outcome"] == 1][metric]
                cannot_link = df_cond[df_cond["link_outcome"] == 0][metric]

                positions.append(pos_counter)
                positions.append(pos_counter + 0.5)
                box_data.append(can_link)
                box_data.append(cannot_link)

                pos_counter += 1.0

            pos_counter += 0.8
            xtick_labels.append(subpop)

        return box_data, positions, xtick_labels

    # Plot Privacy Loss
    privacy_data, privacy_pos, _ = prepare_boxplot_data("pg_column")
    bp1 = ax1.boxplot(privacy_data, positions=privacy_pos, patch_artist=True, widths=0.3, showfliers=False)
    for patch, color in zip(bp1['boxes'], [BAD_COLOR, GOOD_COLOR] * len(privacy_data)):
        patch.set_facecolor(color)

    xticks = SUBPOPS_SMALL

    ax1.set_xticks([np.mean(privacy_pos[i*2:i*2+2]) for i in range(len(SUBPOPS))], labels=SUBPOPS_SMALL)
    ax1.set_xticklabels(xticks, rotation=45, ha="right")
    ax1.set_ylabel("Privacy Risk")
    ax1.set_yscale("log")
    #ax1.set_title("Privacy Loss")

    # Plot Utility Loss
    utility_data, utility_pos, _ = prepare_boxplot_data("ul_column")
    bp2 = ax2.boxplot(utility_data, positions=utility_pos, patch_artist=True, widths=0.3, showfliers=False)
    for patch, color in zip(bp2['boxes'], [BAD_COLOR, GOOD_COLOR] * len(utility_data)):
        patch.set_facecolor(color)

    #ax2.axhline(y=baseline_utility_mean, color='gray', linestyle='--', label='Unique Baseline')

    ax2.set_xticks([np.mean(utility_pos[i*2:i*2+2]) for i in range(len(SUBPOPS))], labels=SUBPOPS_SMALL)
    ax2.set_xticklabels(xticks, rotation=45, ha="right")
    ax2.set_ylabel("Utility Loss")
    #ax2.set_title("Utility Loss")

    # Shared Legend
    legend_handles = [
        plt.Line2D([0], [0], color=BAD_COLOR, lw=4, label="Can Link"),
        plt.Line2D([0], [0], color=GOOD_COLOR, lw=4, label="Cannot Link"),
        #plt.Line2D([0], [0], color='gray', linestyle='--', lw=2, label="Unique Baseline")
    ]
    fig.legend(handles=legend_handles,
           loc="lower center",
           bbox_to_anchor=(0.5, -0.1),
           ncol=3,
           frameon=False)

    plt.tight_layout()
    plt.savefig(f"{PLOT_OUT_PATH}/RECOMB_M3_Privacy_Utility_Boxplots_SideBySide.png")
    plt.savefig(f"{PLOT_OUT_PATH}/RECOMB_M3_Privacy_Utility_Boxplots_SideBySide.pdf")
    plt.close()


if __name__ == "__main__":
    apply_style()
    plot_combo_privacy_utility()