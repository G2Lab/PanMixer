import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

from tools.common.utils import load_data
from plot_util.my_color_palette import apply_style, MAIN_COLORS, SECONDARY_COLORS, FONT
from plot_util.plot_helpers import collect_metric_points, READ_SUBJECTS, PANGENOME_SUBJECTS

from constants import PLOT_OUT_PATH

# ── Experiment numbers — update to match your data ───────────────────────────
GIRAFFE_EXPERIMENT  = 7   # PanMixer read-mapping experiment
CORR_EXPERIMENT     = 7   # experiment used for utility-vs-mapping correlation
#PRIVATE_EXPERIMENT  = 27   # RECOMB PanMixer experiment
#EMPTY_EXPERIMENT    = 5    # subject-removed baseline
NO_EDIT_EXPERIMENT = 3    # unedited baseline (full pangenome)
NO_EDIT_FULL_EXPERIMENT = 6


def _select_panmixer_rows(panmixer_exp):
    """Select one row per subject: highest utility_loss among rows where both
    genotype and haplotype gap scores are negative (i.e. private)."""
    all_df = pd.read_csv(f"./experiments/exp_{panmixer_exp}/data_all.csv")
    mask = (all_df["genotypes_score"] < 0) & (all_df["haplotype_score"] < 0)
    filtered = all_df[mask]
    return filtered.groupby("subject")["utility_loss"].idxmax().values


def plot_A(panmixer_exp=GIRAFFE_EXPERIMENT, figsize=(12, 5), output_prefix=""):
    """Horizontal paired boxplots: Panmixer vs Original for each pangenome type.

    1 row × 3 columns (Perfect / Gapless / MAPQ 60).  Each subplot has 3 paired
    rows (Full, Filtered, Personalized) with Original and Panmixer side-by-side,
    distinguished by color.
    """
    master_full_df = pd.read_csv(f"./experiments/exp_{NO_EDIT_FULL_EXPERIMENT}/data_chr21.csv")
    master_df      = pd.read_csv(f"./experiments/exp_{NO_EDIT_EXPERIMENT}/data_chr21.csv")
    panmixer_df    = pd.read_csv(f"./experiments/exp_{panmixer_exp}/data_chr21.csv")
    keep_rows      = _select_panmixer_rows(panmixer_exp)
    panmixer_df    = panmixer_df.iloc[keep_rows]

    metrics        = ["perfect", "gapless", "mapq60"]
    metrics_labels = ["Perfect", "Gapless", "MAPQ 60"]

    # (label, col_prefix, original_df)
    pangenome_types = [
        ("Full",         "",              master_full_df),
        ("Filtered",     "filtered_",     master_df),
        ("Personalized", "personalized_", master_df),
    ]

    COLOR_ORIGINAL = MAIN_COLORS[1]
    COLOR_PANMIXER = SECONDARY_COLORS[1]

    fig, axes = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=figsize)

    # Positions: pairs at (1,2), (4,5), (7,8) — gap between groups
    positions_orig    = [2, 5, 8]
    positions_panmix  = [1, 4, 7]
    group_centers     = [1.5, 4.5, 7.5]

    for ax, metric, mlabel in zip(axes, metrics, metrics_labels):
        orig_data = []
        panmix_data = []

        for pg_label, col_prefix, orig_df in pangenome_types:
            matched_reads = {
                s: r for s, r in READ_SUBJECTS.items()
                if f"{col_prefix}{s}_total_aligned" in orig_df.columns
                and f"{col_prefix}{s}_total_aligned" in panmixer_df.columns
            }
            orig_pts = collect_metric_points(orig_df, matched_reads, col_prefix=col_prefix)
            panmix_pts = collect_metric_points(panmixer_df, matched_reads, col_prefix=col_prefix)
            orig_arr = 100 * np.array(orig_pts[metric])
            panmix_arr = 100 * np.array(panmix_pts[metric])
            orig_data.append(orig_arr)
            panmix_data.append(panmix_arr)
            print(f"  {mlabel:8s} | {pg_label:13s} | Original mean = {orig_arr.mean():.3f}%"
                  f" (n={orig_arr.size}) | Panmixer mean = {panmix_arr.mean():.3f}%"
                  f" (n={panmix_arr.size}) | Δ = {panmix_arr.mean() - orig_arr.mean():+.3f}%")

        box_kw = dict(vert=False, patch_artist=True, widths=0.8,
                      showfliers=True, flierprops=dict(markersize=4))

        bp_orig = ax.boxplot(orig_data, positions=positions_orig, **box_kw)
        bp_panmix = ax.boxplot(panmix_data, positions=positions_panmix, **box_kw)

        for patch in bp_orig["boxes"]:
            patch.set_facecolor(COLOR_ORIGINAL)
            patch.set_alpha(0.7)
        for patch in bp_panmix["boxes"]:
            patch.set_facecolor(COLOR_PANMIXER)
            patch.set_alpha(0.7)
        for bp in [bp_orig, bp_panmix]:
            for element in ["whiskers", "caps", "medians"]:
                for item in bp[element]:
                    item.set_color("black")

        # Separators between pangenome type groups
        ax.axhline(y=3, color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
        ax.axhline(y=6, color="gray", linestyle="--", linewidth=0.5, alpha=0.5)

        ax.set_title(mlabel)
        ax.grid(axis="x", linestyle="--", alpha=0.4)

    # Y-axis labels: pangenome type names at group centers
    pg_labels = [pg[0] for pg in pangenome_types]
    for ax in axes:
        ax.set_yticks(group_centers)
        ax.set_yticklabels(pg_labels)
        ax.invert_yaxis()

    # Legend
    legend_handles = [
        Line2D([0], [0], color=COLOR_ORIGINAL, marker="s", linestyle="None",
               markersize=10, label="Original", alpha=0.7),
        Line2D([0], [0], color=COLOR_PANMIXER, marker="s", linestyle="None",
               markersize=10, label="PanMixer", alpha=0.7),
    ]
    axes[-1].legend(handles=legend_handles, loc="lower right", framealpha=0.9)

    axes[1].set_xlabel("Percent of Reads")
    plt.tight_layout()
    plt.savefig(f"{PLOT_OUT_PATH}/{output_prefix}M5A.pdf")
    plt.savefig(f"{PLOT_OUT_PATH}/{output_prefix}M5A.png")
    plt.close()


def plot_B(
    panmixer_exp=CORR_EXPERIMENT,
    pangenome_subject=PANGENOME_SUBJECTS[0],
    read_subject=list(READ_SUBJECTS.keys())[0],
    figsize=(12, 4),
    output_prefix="",
):
    """Scatter plots of utility loss vs read-mapping quality metrics."""
    df = pd.read_csv(f"./experiments/exp_{panmixer_exp}/data_chr21.csv")
    df = df[df["subject"] == pangenome_subject]

    print(f"Generating plots: pangenome={pangenome_subject}, read={read_subject}")

    corr_perfect, corr_gapless, corr_mapq60 = [], [], []

    for _, row in df.iterrows():
        util_loss = row["utility_loss_normalized"]
        if pd.isna(row.get(f"{read_subject}_total_aligned")):
            continue
        ta = pd.to_numeric(row.get(f"{read_subject}_total_aligned"))
        tp = pd.to_numeric(row.get(f"{read_subject}_total_perfect"), errors="coerce")
        tg = pd.to_numeric(row.get(f"{read_subject}_total_gapless_softclips_allowed"), errors="coerce")
        tm = pd.to_numeric(row.get(f"{read_subject}_mapping_quality_max_60_reads"), errors="coerce")
        if not pd.isna(tp): corr_perfect.append((util_loss, 100.0 * tp / ta))
        if not pd.isna(tg): corr_gapless.append((util_loss, 100.0 * tg / ta))
        if not pd.isna(tm): corr_mapq60.append((util_loss, 100.0 * tm / ta))

    def safe_r2(pairs):
        if not pairs:
            return np.nan, np.nan, np.array([]), np.array([])
        arr = np.asarray(pairs, dtype=float)
        x, y = arr[:, 0], arr[:, 1]
        m = ~np.isnan(x) & ~np.isnan(y)
        x, y = x[m], y[m]
        if x.size < 3 or np.isclose(np.std(x, ddof=1), 0) or np.isclose(np.std(y, ddof=1), 0):
            return np.nan, np.nan, x, y
        r = np.corrcoef(x, y)[0, 1]
        return r, r * r, x, y

    r_p, r2_p, x_p, y_p = safe_r2(corr_perfect)
    r_g, r2_g, x_g, y_g = safe_r2(corr_gapless)
    r_m, r2_m, x_m, y_m = safe_r2(corr_mapq60)

    fig, axes = plt.subplots(1, 3, figsize=figsize, sharex=True)
    panels = [("Perfect", x_p, y_p, r2_p), ("Gapless", x_g, y_g, r2_g), ("MAPQ 60", x_m, y_m, r2_m)]

    for ax, (name, x, y, r2) in zip(axes, panels):
        if x.size:
            ax.scatter(x, y, alpha=0.6, s=50)
            if np.std(x) > 0 and x.size >= 2:
                slope, intercept = np.polyfit(x, y, 1)
                xline = np.linspace(x.min(), x.max(), 100)
                ax.plot(xline, slope * xline + intercept)
            if not np.isnan(r2):
                if name == "MAPQ 60":
                    ax.text(0.70, 0.95, fr"$R^2$ = {r2:.3f}", transform=ax.transAxes, va="top", ha="right")
                else:
                    ax.text(0.70, 0.95, fr"$R^2$ = {r2:.3f}", transform=ax.transAxes, va="top", ha="left")

        ax.set_title(name)
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.ticklabel_format(style="plain", axis="y")
        ax.yaxis.get_major_formatter().set_useOffset(False)
        ax.yaxis.get_major_formatter().set_scientific(False)

    axes[0].set_ylabel("Percent of Reads")
    axes[1].set_xlabel("Utility Loss")
    fig.tight_layout()
    print(f"Saved to {PLOT_OUT_PATH}/{output_prefix}M5B_{pangenome_subject}_{read_subject}.png")
    plt.savefig(f"{PLOT_OUT_PATH}/{output_prefix}M5B_{pangenome_subject}_{read_subject}.png")
    plt.close()


if __name__ == "__main__":
    apply_style()
    # Original paper figures
    plot_A(GIRAFFE_EXPERIMENT)
    for pg_subj in PANGENOME_SUBJECTS:
        for rd_subj in READ_SUBJECTS:
            plot_B(CORR_EXPERIMENT, pangenome_subject=pg_subj, read_subject=rd_subj)
    # RECOMB figures (disabled — PRIVATE_EXPERIMENT not set)
    # plot_A(PRIVATE_EXPERIMENT, figsize=(5, 1), output_prefix="RECOMB_")
    # for pg_subj in PANGENOME_SUBJECTS:
    #     for rd_subj in READ_SUBJECTS:
    #         plot_B(CORR_EXPERIMENT, pangenome_subject=pg_subj, read_subject=rd_subj,
    #                figsize=(8, 2), output_prefix="RECOMB_")
