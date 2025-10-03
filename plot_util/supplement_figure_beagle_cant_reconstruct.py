import numpy as np
import pandas as pd
from tools.utils import load_data
import json
import os

import matplotlib.pyplot as plt

from my_color_palette import apply_style, MAIN_COLORS, FONT

from constants import (
    PLOT_OUT_PATH,
    EXPERIMENT_PATH,
)

MAIN_EXPERIMENT = 2
CHROMOSOME = 21

def add_accuracy_stats(df, experiment_number, chromosome, overwrite=False):
    if not overwrite and "obfuscated_accuracy" in df.columns:
        return

    if not os.path.exists(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr21/0/accuracy_results.json"):
        print("No accuracy stats found")
        return

    subsets = ["all","complex", "snp_only", "snp_only_maf_0.05", "snp_only_maf_0.5"]
    methods = ["obfuscated_accuracy", "beagle_accuracy", "major_allele_accuracy", "reconstruction_accuracy"]

    accuracy_results = {method: {subset: [] for subset in subsets} for method in methods}

    nan_count = 0
    for i in range(len(df)):
        file_path = f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/chr{chromosome}/{i}/accuracy_results.json"

        if os.path.exists(file_path):
            with open(file_path, "r") as f:
                stats = json.load(f)
            for method in methods:
                for subset in subsets:
                    accuracy_results[method][subset].append(stats[subset][method])
        else:
            for method in methods:
                for subset in subsets:
                    accuracy_results[method][subset].append(np.nan)
                    nan_count += 1

    for method in methods:
        for subset in subsets:
            df[f"{method}_{subset}"] = accuracy_results[method][subset]

    print(f"Accuracy stats total NaNs inserted: {nan_count}")

def read_data():
    data = load_data(MAIN_EXPERIMENT, CHROMOSOME)
    add_accuracy_stats(data, MAIN_EXPERIMENT, CHROMOSOME)

    return data

def plot_beagle_accuracy(experiment):
    data_df = read_data()

    data_df["num_obfuscated_variants"] = data_df["number_of_changed_alleles"]

    # Series to plot
    series = ["beagle_accuracy_all", "major_allele_accuracy_all"]
    nice_names = ["Beagle", "Major Allele Baseline"]

    # Filter necessary columns
    data_df = data_df[["num_obfuscated_variants", "capacity"] + series]

    # Average by capacity
    data_df_averaged = data_df.groupby("capacity").mean().reset_index()

    # Create a single plot
    fig, ax = plt.subplots(figsize=(9, 6))

    # Plot each series on the same axes
    for i, serie in enumerate(series):
        if serie == "major_allele_accuracy_all":
            ax.plot(
                data_df_averaged["num_obfuscated_variants"],
                data_df_averaged[serie],
                label=nice_names[i],
                color="black",
                linestyle="--",
                linewidth=2.5,
                zorder=5
            )
        else:
            ax.plot(
                data_df_averaged["num_obfuscated_variants"],
                data_df_averaged[serie],
                label=nice_names[i],
                color=MAIN_COLORS[i + 1],
                alpha=0.7,
                linewidth=2
            )

    # Highlight regions where Beagle is below the major allele accuracy
    major_allele = data_df_averaged["major_allele_accuracy_all"]
    for serie in ["beagle_accuracy_all"]:
        below = data_df_averaged[serie] < major_allele
        ax.fill_between(
            data_df_averaged["num_obfuscated_variants"],
            data_df_averaged[serie],
            major_allele,
            where=below,
            color="gray",
            alpha=0.2,
            interpolate=True
        )

    ax.set_xlabel("\# of Obfuscated Variants")
    ax.set_ylabel("Imputation Accuracy")
    ax.legend()
    ax.grid(True)

    plt.tight_layout()
    print("Saved to", f"{PLOT_OUT_PATH}/beagle_accuracy_experiment_{experiment}.pdf")
    plt.savefig(f"{PLOT_OUT_PATH}/beagle_accuracy_experiment_{experiment}.pdf", dpi=300)
    plt.close()

def plot_beagle_accuracy(experiment):
    data_df = read_data()

    data_df["num_obfuscated_variants"] = data_df["number_of_changed_alleles"]

    # Series to plot
    series = ["reconstruction_accuracy_all"]
    nice_names = ["Beagle"]

    # Filter necessary columns
    data_df = data_df[["num_obfuscated_variants", "capacity"] + series]

    # Average by capacity
    data_df_averaged = data_df.groupby("capacity").mean().reset_index()

    # Create a single plot
    fig, ax = plt.subplots(figsize=(9, 6))

    # Plot each series on the same axes
    for i, serie in enumerate(series):
        ax.plot(
            data_df_averaged["num_obfuscated_variants"],
            data_df_averaged[serie],
            label=nice_names[i],
            color=MAIN_COLORS[i + 1],
            alpha=0.7,
            linewidth=2
        )

    ax.set_xlabel("\# of Obfuscated Variants")
    ax.set_ylabel("Reconstruction Accuracy")
    ax.legend()
    ax.grid(True)

    plt.tight_layout()
    plt.savefig(f"{PLOT_OUT_PATH}/beagle_reconstruction_accuracy_{experiment}.pdf", dpi=300)
    plt.close()

def plot_beagle_delta_vs_major(experiment):
    df = read_data()
    df["num_obfuscated_variants"] = df["number_of_changed_alleles"]
    df = df[["num_obfuscated_variants","capacity","beagle_accuracy_all","major_allele_accuracy_all"]].dropna()

    df["delta"] = df["beagle_accuracy_all"] - df["major_allele_accuracy_all"]
    grouped = df.groupby("capacity")

    # Bootstrap CIs per capacity
    rng = np.random.default_rng(0)
    rows = []
    B = 2000
    for cap, g in grouped:
        vals = g["delta"].to_numpy()
        if len(vals) == 0: 
            continue
        boot = []
        for _ in range(B):
            samp = rng.choice(vals, size=len(vals), replace=True)
            boot.append(np.mean(samp))
        boot = np.array(boot)
        rows.append({
            "capacity": cap,
            "mean": np.mean(vals),
            "lo": np.percentile(boot, 2.5),
            "hi": np.percentile(boot, 97.5),
            "n": len(vals),
            "num_obfuscated_variants": g["num_obfuscated_variants"].mean()
        })
    ci = pd.DataFrame(rows).sort_values("num_obfuscated_variants")

    fig, ax = plt.subplots(figsize=(9,6))
        # Beagle reconstruction curve + shading
    line_beagle, = ax.plot(ci["num_obfuscated_variants"], ci["mean"],
                           color=MAIN_COLORS[1], linewidth=2)
    ax.fill_between(ci["num_obfuscated_variants"], ci["lo"], ci["hi"],
                    alpha=0.2, color=MAIN_COLORS[1])

    # Baseline line
    line_baseline = ax.axhline(0, linestyle="--", linewidth=1.5, color="black")

    # Second legend (baseline)
    ax.legend([line_beagle,line_baseline], ["Beagle Reconstruction", "Major Allele baseline"],
              loc="lower left")

    ax.set_xlabel("\# of Obfuscated Variants")
    ax.set_ylabel("Accuracy difference")
    #ax.set_title(f"Beagle underperforms the Major Allele baseline in {below}/{len(ci)} capacity bins")
    ax.grid(True)
    plt.tight_layout()
    print(f"{PLOT_OUT_PATH}/beagle_delta_major_{experiment}.pdf")
    plt.savefig(f"{PLOT_OUT_PATH}/beagle_delta_major_{experiment}.pdf", dpi=300)
    plt.close()

if __name__ == "__main__":
    apply_style()
    plot_beagle_accuracy(MAIN_EXPERIMENT)
    plot_beagle_delta_vs_major(MAIN_EXPERIMENT)
