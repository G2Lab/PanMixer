import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

from tools.common.utils import load_data_all
from plot_util.my_color_palette import apply_style, SECONDARY_COLORS, MAIN_COLORS

from constants import (
    EXPERIMENT_PATH,
    PLOT_OUT_PATH,
)

EXPERIMENT_NUMBER = 8

def _load_subject_linking_scores(experiment_number, i):
    """Sum reverse genotype scores across all 22 chromosomes for subject index i.
    Returns (scores_others, score_obfuscated) — the per-candidate linking score
    distribution against the database (excluding self) and the linking score
    between the obfuscated and original individual."""
    chr1_dir = f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/chr1/{i}"
    scores_others = np.load(f"{chr1_dir}/reverse_genotypes_scores_others.npy").astype(float)
    score_obfuscated = np.load(f"{chr1_dir}/reverse_genotypes_score_obfuscated.npy").astype(float)
    for chrom in range(2, 23):
        scores_others += np.load(f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/chr{chrom}/{i}/reverse_genotypes_scores_others.npy")
        score_obfuscated += np.load(f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/chr{chrom}/{i}/reverse_genotypes_score_obfuscated.npy")
    return scores_others, float(score_obfuscated[0])


def plot_gap_score_boxplots():
    data_df = load_data_all(EXPERIMENT_NUMBER)
    subjects = data_df["subject"].astype(str).values
    populations = data_df["population_code"].astype(str).values

    others_data = []
    obf_scores = []
    for i in range(len(subjects)):
        scores_others, score_obfuscated = _load_subject_linking_scores(EXPERIMENT_NUMBER, i)
        others_data.append(scores_others)
        obf_scores.append(score_obfuscated)

    order = np.argsort(obf_scores)
    subjects = subjects[order]
    populations = populations[order]
    others_data = [others_data[k] for k in order]
    obf_scores = [obf_scores[k] for k in order]

    labels = [f"{s} ({p})" for s, p in zip(subjects, populations)]

    n = len(subjects)
    x_min = min(d.min() for d in others_data + [np.array(obf_scores)])
    x_max = max(d.max() for d in others_data + [np.array(obf_scores)])
    pad = 0.02 * (x_max - x_min)
    x_grid = np.linspace(x_min - pad, x_max + pad, 512)

    row_height = 1.0
    overlap = 0.4
    spacing = row_height * (1.0 - overlap)

    fig, ax = plt.subplots(figsize=(8, max(4, 0.25 * n)))

    fill_color = MAIN_COLORS[1]
    line_color = MAIN_COLORS[3]

    for i, d in enumerate(others_data):
        baseline = i * spacing
        if d.size > 1 and np.std(d) > 0:
            kde = gaussian_kde(d)
            density = kde(x_grid)
            density = density / density.max() * row_height
        else:
            density = np.zeros_like(x_grid)

        ax.fill_between(x_grid, baseline, baseline + density,
                        color=fill_color, alpha=0.7, edgecolor="black", linewidth=0.5,
                        zorder=2 + i * 0.001)

        # Line marking the obfuscated-vs-original linking score
        obf = obf_scores[i]
        if d.size > 1 and np.std(d) > 0:
            line_top = baseline + float(kde(np.array([obf]))[0]) / kde(x_grid).max() * row_height
        else:
            line_top = baseline + row_height
        ax.vlines(obf, baseline, line_top, color=line_color, linewidth=1.5,
                  zorder=2 + i * 0.001 + 0.0005)

    ax.set_yticks([i * spacing + row_height * 0.2 for i in range(n)])
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_ylim(-0.2, (n - 1) * spacing + row_height + 0.2)
    ax.set_xlabel("Linking score")
    ax.set_ylabel("Individual")

    legend_handles = [
        plt.Line2D([0], [0], color=line_color, linewidth=1.5,
                   label="Target's obfuscated-vs-original linking score"),
    ]
    ax.legend(handles=legend_handles, loc="upper center",
              bbox_to_anchor=(0.5, -0.05), frameon=False)

    plt.tight_layout()
    plt.savefig(f"{PLOT_OUT_PATH}/gap_score_ridge_exp{EXPERIMENT_NUMBER}.pdf")
    plt.savefig(f"{PLOT_OUT_PATH}/gap_score_ridge_exp{EXPERIMENT_NUMBER}.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    apply_style()
    plot_gap_score_boxplots()
