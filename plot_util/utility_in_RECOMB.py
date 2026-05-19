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
MAIN_EXPERIMENT      = 7    # main PanMixer obfuscation experiment
HAVE_MARKER          = True

CORRELATIONS_OUT = PLOT_OUT_PATH + "/correlations_"


def _prepare_main_df():
    """Load the main experiment and add the derived columns shared by both plots."""
    data_df = load_data_all(MAIN_EXPERIMENT)
    data_df["ld_loss"] = data_df["ld_sums"] / data_df["ld_counts"]
    data_df["privacy_loss"] = 1 - data_df["pmi_gain_normalized"]
    data_df["utility_loss"] = data_df["utility_loss"] / data_df["true_max_utility_loss"]
    return data_df


def _prepare_baseline_df():
    """Load the removed-individual baseline and add the derived columns."""
    privacy_baseline = load_data_all(ONE_PRIVACY_BASELINE)
    privacy_baseline["ld_loss"] = privacy_baseline["ld_sums"] / privacy_baseline["ld_counts"]
    return privacy_baseline


# ─────────────────────────────────────────────────────────────────────────────
# M4: privacy-loss-vs-utility-loss line plots with broken-axis baselines
# ─────────────────────────────────────────────────────────────────────────────
def plot_M4(data_df=None, privacy_baseline=None):
    if data_df is None:
        data_df = _prepare_main_df()
    if privacy_baseline is None:
        privacy_baseline = _prepare_baseline_df()

    lowest_gapscore_neg = find_when_gapscore_is_neg(data_df, "genotypes_score")
    print(lowest_gapscore_neg)

    subplots = ["All", "SNPs Only",  "SNPs (MAF $<$ 0.05)", "LD"]
    subplot_columns = ["wd_af", "wd_af_snp_only", "wd_af_snps_0_05", "ld_loss"]
    af_indices = {0, 1, 2}  # AF plots that get a baseline subplot above

    series = ["All", "EAS", "AMR", "SAS", "AFR"]
    series_full = ["All", "East Asian Ancestry", "American Ancestry", "South Asian Ancestry", "African Ancestry"]

    fig = plt.figure(figsize=(18, 6))
    gs = gridspec.GridSpec(2, len(subplots), figure=fig, height_ratios=[1, 6], hspace=0.18, wspace=0.7)
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

    for i, (subplot, y_col) in enumerate(zip(subplots, subplot_columns)):
        ax = axes[i]
        ax_top = top_axes[i]

        min_x = float("inf")
        max_x = 0

        for j, serie in enumerate(series):
            if serie != "All":
                data_df_subset = data_df[data_df["population"] == series_full[j]].select_dtypes(include=[np.number])
            else:
                data_df_subset = data_df.select_dtypes(include=[np.number])

            utility_mean = data_df_subset.groupby("capacity").mean()

            x_pmi = utility_mean["privacy_loss"]

            min_x = min(min_x, x_pmi.min())
            max_x = max(max_x, x_pmi.max())

            if serie == "All":
                line, = ax.plot(x_pmi, (utility_mean[y_col]), label=serie, color="black", zorder=30)
            else:
                marker_type, spacing_every = line_markers[j-1]
                if j == 1:
                    color = SECONDARY_COLORS[j + 2]
                elif j == 2:
                    color = SECONDARY_COLORS[j]
                else:
                    color = SECONDARY_COLORS[j + 1]

                line, = ax.plot(x_pmi, (utility_mean[y_col]), label=serie, color=color, zorder=3, marker=marker_type, markevery=(spacing_every, 32), markersize=15)
            if i == 0:
                legend_handles.append(line)

            if serie == "All" and HAVE_MARKER:
                add_marker(ax, x_pmi, (utility_mean[y_col]), target_x=lowest_gapscore_neg)

        ax.set_xlabel("Privacy Risk")
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
    series.append("Removed")
    if HAVE_MARKER:
        legend_handles.append(plt.Line2D([0], [0], color="blue", marker='o', linestyle='None', label="Loss when linking fails"))
        series.append("$\\epsilon_{\\mathrm{private}}$")

    fig.legend(handles=legend_handles, labels=series, loc='lower center', ncol=len(series), bbox_to_anchor=(0.5, -0.08))

    # Flip x-axis
    for ax in axes[1:3]:
        ax.sharey(axes[0])

    for ax in axes[0:3]:
        ax.set_ylim(0, 0.01)
    axes[3].set_ylim(0, 0.02)

    axes[3].set_ylabel("LD Loss")

    for ax in axes:
        ax.set_xlim(1, 0.0001)
        ax.set_xscale("log")
    for ax_top in top_axes:
        if ax_top is not None:
            ax_top.set_xscale("log")
            ax_top.set_xlim(1, 0.0001)

    plt.tight_layout(rect=[0, 0.12, 1, 0.95])
    plt.savefig(f"{PLOT_OUT_PATH}/RECOMB_M4_v2_downstream_analysis.pdf")
    plt.savefig(f"{PLOT_OUT_PATH}/RECOMB_M4_v2_downstream_analysis.png", dpi=300)
    plt.close()


# ─────────────────────────────────────────────────────────────────────────────
# Correlation plots: utility loss vs. each downstream metric
# ─────────────────────────────────────────────────────────────────────────────
DOWNSTREAM_COLUMNS = {
    "ld_loss": "LD loss",
    "wd_af": "AF",
    "wd_af_snp_only": "AF (SNPs only)",
    "wd_af_snps_0_05": "AF (SNPs MAF $<$ 0.05)",
}


def _build_correlation_metrics(data_df):
    """Return {pretty_name: (utility_loss_array, downstream_array)} for each downstream column."""
    metrics = {}
    for column, pretty_name in DOWNSTREAM_COLUMNS.items():
        metrics[pretty_name] = (
            data_df["utility_loss"].to_numpy(),
            data_df[column].to_numpy(),
        )
    return metrics


def _annotate_r2(ax, x, y):
    slope, intercept = np.polyfit(x, y, 1)
    y_pred = slope * x + intercept
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - (ss_res / ss_tot)
    ax.text(0.05, 0.95, f"R² = {r2:.2f}", transform=ax.transAxes,
            verticalalignment='top', fontproperties=FONT)


def create_binned_line_plot(stat_name, x, y, analysis_type, x_label="Utility Loss", num_bins=20):
    df_temp = pd.DataFrame({'x': x, 'y': y})
    df_temp['bin'] = pd.cut(df_temp['x'], bins=num_bins)
    stats = df_temp.groupby('bin', observed=True)['y'].agg(['mean', 'min', 'max']).reset_index()
    stats['x_mid'] = stats['bin'].apply(lambda b: b.mid).astype(float)

    fig, ax = plt.subplots(figsize=(3, 3))
    ax.fill_between(stats['x_mid'], stats['min'], stats['max'],
                    color=SECONDARY_COLORS[6], alpha=0.2, label='Range')
    ax.plot(stats['x_mid'], stats['mean'],
            color=SECONDARY_COLORS[6], marker='o', markersize=4, linewidth=1.5, label='Mean')

    _annotate_r2(ax, x, y)
    ax.set_xlabel(x_label)
    ax.set_ylabel(stat_name)

    plt.tight_layout()
    plt.savefig(CORRELATIONS_OUT + f"RECOMB_{analysis_type}_{stat_name}.pdf")
    plt.savefig(CORRELATIONS_OUT + f"RECOMB_{analysis_type}_{stat_name}.png")
    plt.close()


def plot_correlations(data_df=None):
    if data_df is None:
        data_df = _prepare_main_df()
    metrics = _build_correlation_metrics(data_df)
    for stat_name, (x, y) in metrics.items():
        create_binned_line_plot(stat_name, x, y, "downstream_in", "Utility Loss")


# ─────────────────────────────────────────────────────────────────────────────
# Entry point: regenerate both sets of plots
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    apply_style()
    data_df = _prepare_main_df()
    privacy_baseline = _prepare_baseline_df()
    plot_M4(data_df=data_df, privacy_baseline=privacy_baseline)
    plot_correlations(data_df=data_df)
