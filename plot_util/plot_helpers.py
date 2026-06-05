"""
Shared utilities for plot_util scripts.

Centralises constants and helper functions that were copy-pasted across
multiple plotting files.
"""

import numpy as np
import pandas as pd

from plot_util.my_color_palette import FONT
from constants import DEMOGRAPHICS_CSV

# ── Read-mapping metric definitions ──────────────────────────────────────────
GIRAFFE_METRICS = {
    "perfect": ("total_perfect",                   "total_aligned"),
    "gapless": ("total_gapless_softclips_allowed", "total_aligned"),
    "mapq60":  ("mapping_quality_max_60_reads",    "total_aligned"),
}

# ── Subjects used for read-mapping experiments ────────────────────────────────
READ_SUBJECTS = {
    "HG00138": "EUR",
    "HG00635": "EAS",
    "HG01112": "AMR",
    "HG01600": "EAS",
    "HG02698": "SAS",
    "NA12778": "EUR",
    "NA18853": "AFR",
}

PANGENOME_SUBJECTS = [
    "HG00438",
    "HG00733",
    "HG02145",
    "HG03492",
]


def collect_metric_points(df, subjects, metrics=None, col_prefix=""):
    """Compute per-metric read-mapping ratios for each subject.

    Returns a dict mapping metric name -> list of ratio values (floats).
    """
    if metrics is None:
        metrics = GIRAFFE_METRICS
    out = {m: [] for m in metrics}
    for subj in subjects:
        for m, (num_sfx, den_sfx) in metrics.items():
            num_col = f"{col_prefix}{subj}_{num_sfx}"
            den_col = f"{col_prefix}{subj}_{den_sfx}"
            if num_col not in df.columns or den_col not in df.columns:
                continue
            num = pd.to_numeric(df[num_col], errors="coerce")
            den = pd.to_numeric(df[den_col], errors="coerce")
            ratio = num / den
            out[m].extend(ratio[ratio.notna() & den.notna() & (den != 0)].tolist())
    return out


def add_marker(ax, x_values, y_values, target_x, xytext=(-35, 0)):
    """Add an annotated interpolated point marker at target_x on ax."""
    x_arr = np.asarray(x_values)
    y_arr = np.asarray(y_values)
    i1, i2 = np.argsort(np.abs(x_arr - target_x))[:2]
    x1, x2 = x_arr[i1], x_arr[i2]
    y1, y2 = y_arr[i1], y_arr[i2]
    target_y = y1 + (y2 - y1) * (target_x - x1) / (x2 - x1)
    ax.scatter(target_x, target_y, color="blue", zorder=5)
    ax.vlines(target_x, 0, target_y, color="blue", linestyle="--", zorder=5)
    ax.annotate(
        f"{target_y:.3f}",
        (target_x, target_y),
        textcoords="offset points",
        xytext=xytext,
        ha="center",
        fontsize=15,
        zorder=5,
        fontproperties=FONT,
    )


def find_when_gapscore_is_neg(df, score_col):
    """Return the minimum privacy_loss where score_col first goes negative.

    Iterates over all subjects, finds rows where score_col < 0, and returns
    the minimum of each subject's maximum privacy_loss in that region.
    """
    lowest = []
    for subject in df["subject"].unique():
        neg = df[(df["subject"] == subject) & (df[score_col] < 0)]
        print(len(neg))
        lowest.append(neg["privacy_loss"].max())
    print(lowest)
    return min(lowest)


def load_demographic_data():
    return pd.read_csv(DEMOGRAPHICS_CSV)


def add_demographic_information(data_df, demographic_data):
    return data_df.merge(demographic_data, on="subject")
