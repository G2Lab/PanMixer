import numpy as np
import pickle
import os
import matplotlib.pyplot as plt
import pandas as pd
import json
from matplotlib.lines import Line2D
from tools.utils import load_data
from my_color_palette import apply_style, MAIN_COLORS, SECONDARY_COLORS, FONT

from constants import (
    PLOT_OUT_PATH,
)

READ_SUBJECTS = ["HG00346", "HG02688", "HG00473", "HG01250", "HG02703"]
READ_POPULATIONS = ["EUR", "SAS", "EAS", "AMR", "AFR"]

REMOVED_POPULATION = ["EAS", "AMR", "AFR", "SAS"]

read_subjects = {
    "HG00138": "EUR",
    "HG00635": "EAS",
    "HG01112": "AMR",
    "HG01600": "EAS",
    "HG02698": "SAS",
    "NA12778": "EUR",
    "NA18853": "AFR",
}

def plot_A():
    empty_chr21_df = pd.read_csv("./experiments/exp_5/data_chr21.csv")
    master_chr21_df = pd.read_csv("./experiments/exp_6/data_chr21.csv")
    panmixer_chr21_df = pd.read_csv("./experiments/exp_2/data_chr21.csv")

    keep_rows = [21, 72, 126, 172]
    panmixer_chr21_df = panmixer_chr21_df.iloc[keep_rows]


    # Which subjects?
    # read_subjects = ["ind1", "ind2", ...]  # whatever yours are

    # Define which columns form each metric: numerator vs denominator suffixes
    # (Adjust these suffixes if your actual column names differ.)
    METRICS = {
        "perfect": ("total_perfect", "total_aligned"),
        "gapless": ("total_gapless_softclips_allowed", "total_aligned"),
        "mapq60":  ("mapping_quality_max_60_reads",        "total_aligned"),
    }

    def collect_metric_points(df: pd.DataFrame, read_subjects, metrics=METRICS):
        out = {m: [] for m in metrics}
        for subj in read_subjects:
            for m, (num_sfx, den_sfx) in metrics.items():
                num_col = f"{subj}_{num_sfx}"
                den_col = f"{subj}_{den_sfx}"

                # Coerce to numeric in case CSV has strings; non-numeric -> NaN
                num = pd.to_numeric(df.get(num_col), errors="coerce")
                den = pd.to_numeric(df.get(den_col), errors="coerce")

                # Compute ratio safely and drop NaNs / inf
                ratio = num / den
                ratio = ratio[ratio.notna() & pd.notna(den) & (den != 0)]


                # Extend the list with all valid points from this subject
                out[m].extend(ratio.tolist())
        return out

    empty_data_pts   = collect_metric_points(empty_chr21_df, read_subjects)
    master_data_pts  = collect_metric_points(master_chr21_df, read_subjects)
    panmixer_data_pts = collect_metric_points(panmixer_chr21_df, read_subjects)


    removed = {
        "perfect": 100 * np.array(empty_data_pts["perfect"]),
        "gapless": 100 * np.array(empty_data_pts["gapless"]),
        "mapq60": 100 * np.array(empty_data_pts["mapq60"]),
    }
    panmixer = {
        "perfect": 100 * np.array(panmixer_data_pts["perfect"]),
        "gapless": 100 * np.array(panmixer_data_pts["gapless"]),
        "mapq60": 100 * np.array(panmixer_data_pts["mapq60"]),
    }
    unedited = {
        "perfect": 100 * np.array(master_data_pts["perfect"]),
        "gapless": 100 * np.array(master_data_pts["gapless"]),
        "mapq60": 100 * np.array(master_data_pts["mapq60"]),
    }

    # Organize into rows
    conditions = ["Panmixer ($\epsilon_{\mathrm{private}}$)", "Original"]
    metrics = ["perfect", "gapless", "mapq60"]
    metrics_clean_name = ["Perfect", "Gapless", "MAPQ 60"]

    data = {
        "Removed": removed,
        "Panmixer ($\epsilon_{\mathrm{private}}$)": panmixer,
        "Original": unedited,
    }

    fig, axes = plt.subplots(
        nrows=1, ncols=3, sharey=True, figsize=(10, 3)
    )

    for ax, metric, metric_clean in zip(axes, metrics, metrics_clean_name):
        # Collect data for all three conditions
        metric_data = [data[c][metric] for c in conditions]

        """
        parts = ax.violinplot(
            metric_data,
            vert=False,
            showmeans=True,
            showmedians=False,
            showextrema=True,
            widths=0.9,
        )

        # Make violins a bit transparent; keep default color
        for pc in parts['bodies']:
            pc.set_alpha(0.6)

        # y positions are 1..N in violinplot
        y_pos = np.arange(1, len(conditions) + 1)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(conditions)

        ax.set_title(f"{metric_clean}")
        ax.set_xlabel("percent of reads")
        ax.grid(axis="x", linestyle="--", alpha=0.5)
        """

        # Horizontal boxplot
        bp = ax.boxplot(
            metric_data,
            vert=False,   # horizontal
            patch_artist=True,
            labels=conditions,
        )

        print(f"Mean {metric_clean}", np.mean(metric_data[0]), np.mean(metric_data[1]))
        
        ax.set_title(metric_clean)
        ax.set_xlabel(f"Percent of Reads")
        #ax.set_xlim(0, 100)  # adjust based on your data
        ax.grid(axis="x", linestyle="--", alpha=0.5)

    # Shared y-label
    plt.tight_layout()
    plt.savefig(f"{PLOT_OUT_PATH}/M5A.pdf")
    plt.close()



def plot_B():
    # correlation utility metric vs mappe dperfectly, mapped gapless, map
    panmixer_chr21_df = pd.read_csv("./experiments/exp_2/data_chr21.csv")
    panmixer_chr21_df = panmixer_chr21_df[panmixer_chr21_df["subject"] == "HG00438"]

    # Metric definitions: (numerator_col_suffix, denominator_col_suffix)
    METRICS = {
        "perfect": ("total_perfect", "total_aligned"),
        "gapless": ("total_gapless_softclips_allowed", "total_aligned"),
        "mapq60":  ("mapping_quality_max_60_reads", "total_aligned"),
    }

    corr_perfect = []
    corr_gapless = []
    corr_mapq60 = []

    # Utility-loss column name

    for index,row in panmixer_chr21_df.iterrows():
        util_loss = row["utility_loss_normalized"]
        
        for read_subject in [list(read_subjects.keys())[0]]:
            if pd.isna(row[f"{read_subject}_total_aligned"]):
                continue

            ta = pd.to_numeric(row.get(f"{read_subject}_total_aligned"))

            tp = pd.to_numeric(row.get(f"{read_subject}_total_perfect"), errors="coerce")
            tg = pd.to_numeric(row.get(f"{read_subject}_total_gapless_softclips_allowed"), errors="coerce")
            tm = pd.to_numeric(row.get(f"{read_subject}_mapping_quality_max_60_reads"), errors="coerce")

            if not pd.isna(tp):
                corr_perfect.append((util_loss, 100.0 * tp / ta))
            if not pd.isna(tg):
                corr_gapless.append((util_loss, 100.0 * tg / ta))
            if not pd.isna(tm):
                corr_mapq60.append((util_loss, 100.0 * tm / ta))

    # Helper to compute r and r^2 safely
    def safe_r2(pairs):
        if not pairs:
            return np.nan, np.nan, np.array([]), np.array([])
        arr = np.asarray(pairs, dtype=float)
        x = arr[:, 0]
        y = arr[:, 1]
        m = ~np.isnan(x) & ~np.isnan(y)
        x = x[m]; y = y[m]
        if x.size < 3 or np.isclose(np.std(x, ddof=1), 0) or np.isclose(np.std(y, ddof=1), 0):
            return np.nan, np.nan, x, y
        r = np.corrcoef(x, y)[0, 1]
        return r, r*r, x, y

    r_perfect, r2_perfect, x_p, y_p   = safe_r2(corr_perfect)
    r_gapless, r2_gapless, x_g, y_g   = safe_r2(corr_gapless)
    r_mapq60, r2_mapq60, x_m, y_m     = safe_r2(corr_mapq60)

    # Plot: 3 subplots, shared x = utility_loss
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharex=True)
    panels = [
        ("Perfect", x_p, y_p, r2_perfect),
        ("Gapless", x_g, y_g, r2_gapless),
        ("MAPQ 60",  x_m, y_m, r2_mapq60),
    ]

    plt.ticklabel_format(style='plain', axis='y')

    for ax, (name, x, y, r2) in zip(axes, panels):
        if x.size:
            ax.scatter(x, y, alpha=0.6)
            # Least-squares line if x has spread
            if np.std(x) > 0 and x.size >= 2:
                slope, intercept = np.polyfit(x, y, 1)
                xline = np.linspace(x.min(), x.max(), 100)
                yline = slope * xline + intercept
                ax.plot(xline, yline)
            if not np.isnan(r2):
                ax.text(0.70, 0.95, fr"$R^2$ = {r2:.3f}",
                        transform=ax.transAxes, va="top", ha="left")
        ax.set_title(name)
        ax.set_xlabel("Utility Loss")
        ax.grid(True, linestyle="--", alpha=0.4)

    axes[2].yaxis.get_major_formatter().set_scientific(False)

    plt.ticklabel_format(useOffset=False, style='plain', axis='y')


    axes[0].set_ylabel("Percent of Reads")
    fig.tight_layout()
    print("finished plots to", f"{PLOT_OUT_PATH}/M5B.pdf")
    plt.savefig(f"{PLOT_OUT_PATH}/M5B.pdf")
    plt.close()

if __name__ == "__main__":
    apply_style()

    plot_A()
    plot_B()