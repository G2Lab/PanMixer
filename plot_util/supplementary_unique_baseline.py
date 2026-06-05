import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score

from tools.common.utils import load_data_all
from plot_util.my_color_palette import apply_style, MAIN_COLORS, GOOD_COLOR, BAD_COLOR
import matplotlib.ticker as mtick
import matplotlib.ticker as ticker

from constants import (
    PLOT_OUT_PATH,
    SUBPOPS,
    SUBPOPS_SMALL,
)

# ── Experiment numbers — update to match your data ───────────────────────────
GAPSCORE_RESULTS  = 27   # experiment containing gap-score / linking results
BASELINE_UNIQUE   = 19   # experiment where only unique-sequence variants are obfuscated
UNEDITED_BASELINE = 6    # experiment with the original, unedited pangenome

# Define the formatter function
def human_format(x, pos):
    if abs(x) >= 1e6:
        return f'{x*1e-6:.1f}M'.replace('.0M', 'M')
    if abs(x) >= 1e3:
        return f'{x*1e-3:.1f}k'.replace('.0k', 'k')
    return f'{x:g}'

def plot_supplementary_unique_baseline():
    panmixer_df = load_data_all(GAPSCORE_RESULTS)
    unique_df = load_data_all(BASELINE_UNIQUE)
    unedited_df = load_data_all(UNEDITED_BASELINE)

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))

    #x-axis pmi
    #y-axis gapscore genotypes

    panmixer_df["pmi_gain_norm"] = panmixer_df["pmi_gain"] / panmixer_df["max_pmi_gain"]
    unedited_df["pmi_gain_norm"] = unedited_df["pmi_gain"] / unedited_df["max_pmi_gain"]

    panmixer_df["privacy_risk"] = 1 - panmixer_df["pmi_gain_norm"]
    unedited_df["privacy_risk"] = 1 - unedited_df["pmi_gain_norm"]

    assert "genotypes_score" in panmixer_df.columns
    assert "genotypes_score" in unique_df.columns
    #assert "genotypes_score" in unedited_df.columns

    #ax.scatter(unedited_df["pmi_gain_norm"], unedited_df["genotypes_score"], label="Unedited")
    # scatter for poitns where privacy risk > 0

    #strange bug with subjects column


    above_0_scores = panmixer_df[panmixer_df["genotypes_score"] > 0]
    below_0_scores = panmixer_df[panmixer_df["genotypes_score"] <= 0]

    #remove capacity = 1
    above_0_scores = above_0_scores[above_0_scores["capacity"] > 0.000000000001]

    #average by capacity
    above_0_scores_grouped = above_0_scores.groupby("capacity").mean(numeric_only=True)
    below_0_scores_grouped = below_0_scores.groupby("capacity").mean(numeric_only=True)

    ax.scatter(above_0_scores_grouped["privacy_risk"], above_0_scores_grouped["genotypes_score"], label="PanMixer (Can Link; Gap Score $\geq$ 0)", color=BAD_COLOR)
    ax.scatter(below_0_scores_grouped["privacy_risk"], below_0_scores_grouped["genotypes_score"], label="PanMixer (Cannot Link; Gap Score $<$ 0)", color=GOOD_COLOR)

    #draw horizontal lines for unique baseline since it does not have a PMI
    ax.axhline(y=unique_df["genotypes_score"].mean(), color='gray', linestyle='--', label='Removed Unique Nodes')
    
    #add an x mark for unedited

    privacy_risk_unedited = 1
    gap_score_unedited = unedited_df["genotypes_score"].mean()
    ax.scatter(privacy_risk_unedited, gap_score_unedited, color='black', marker='x', label='Original')

    #draw a sharp line at 0
    ax.axhline(y=0, color='black', linestyle='-', label = "Gap Score = 0")
    
    ax.set_xlabel("Privacy Risk")
    ax.set_ylabel("Gap Score")

    ax.xaxis.set_major_formatter(ticker.FuncFormatter(human_format))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(human_format))

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3), 
              ncol=1, frameon=True, borderaxespad=0.)
    plt.tight_layout()
    
    os.makedirs(PLOT_OUT_PATH, exist_ok=True)
    plt.savefig(f"{PLOT_OUT_PATH}/RECOMB_Supplementary_Unique_Baseline.pdf")
    plt.savefig(f"{PLOT_OUT_PATH}/RECOMB_Supplementary_Unique_Baseline.png")
    plt.close()





if __name__ == "__main__":
    apply_style()
    plot_supplementary_unique_baseline()