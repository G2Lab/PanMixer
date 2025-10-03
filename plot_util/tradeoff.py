import numpy as np
import pandas as pd
from tools.utils import load_data_multichromosome, load_data_all
from my_color_palette import apply_style, MAIN_COLORS, FONT
import matplotlib.pyplot as plt
import json

from constants import (
    PLOT_OUT_PATH,
    STARTING_DATA_PATH,
    TOTAL_UTILITY_JSON,
    DEMOGRAPHICS_CSV,
    ALL_POPULATIONS
)

PRIVACY_OPTIMIZATION_EXPERIMENT = 2
UTILITY_OPTIMIZATION_EXPERIMENT = 3


def load_demographic_data():
    # Load demographic data
    return pd.read_csv(DEMOGRAPHICS_CSV)

def add_demographic_information(data_df, demographic_data):
    data_df = data_df.merge(demographic_data, on="subject")
    return data_df

def add_marker(ax, x_values, y_values, target_x=0.0018):
    # Find indices of the two closest points to target_x
    x_values = x_values.to_numpy()
    y_values = y_values.to_numpy()
    closest_points = np.argsort(np.abs(x_values - target_x))[:2]
    
    # Extract the x and y values of the two closest points
    x1, x2 = x_values[closest_points]
    y1, y2 = y_values[closest_points]
    
    # Calculate the y-value at target_x using linear interpolation
    target_y = y1 + (y2 - y1) * (target_x - x1) / (x2 - x1)
    
    # Plot the marker at the target point
    ax.scatter(target_x, target_y, color="blue", zorder=5)
    ax.vlines(target_x, 0, target_y, color="blue", linestyle="--", zorder=5)
    ax.annotate(f"{target_y:.2f}", 
                (target_x, target_y), 
                textcoords="offset points", 
                xytext=(20, 0), 
                ha='center', fontsize=8,
                zorder=5,
                label= "Privacy Gain when Linking Fails")
    
def plot_figure_A():
    data_df_utility = load_data_multichromosome(UTILITY_OPTIMIZATION_EXPERIMENT)
    data_df_privacy = load_data_multichromosome(PRIVACY_OPTIMIZATION_EXPERIMENT)

    plot_subjects = data_df_utility[1]["subject"].unique()
    print(plot_subjects)
    
    aggre_data_df_utility = load_data_all(UTILITY_OPTIMIZATION_EXPERIMENT)
    aggre_data_df_privacy = load_data_all(PRIVACY_OPTIMIZATION_EXPERIMENT)

    subplots = ["All", "East Asian Ancestry", "American Ancestry", "South Asian Ancestry", "African Ancestry"]

    # Adjust hspace here to add more space between rows
    fig, axes = plt.subplots(2, len(subplots), sharex=True, sharey=True, 
                             figsize=(15, 10),
                             gridspec_kw={'hspace': 0.3})

    axes = np.atleast_2d(axes)

    # Add global row titles using fig.text
    fig.text(0.5, 0.915, "Optimized for Fixed Utility ($\\eta$)", ha="center", va="bottom", fontproperties=FONT)
    fig.text(0.5, 0.48, "Optimized for Fixed Privacy ($\\epsilon$)", ha="center", va="bottom", fontproperties=FONT)

    for i, population in enumerate(subplots):
        if population != "All":
            aggre_data_df_utility_subset = aggre_data_df_utility[aggre_data_df_utility["population"] == population]
        else:
            aggre_data_df_utility_subset = aggre_data_df_utility
        
        aggre_data_df_utility_subset = aggre_data_df_utility_subset.select_dtypes(include=[np.number])
        
        ax = axes[0, i]

        aggre_data_df_utility_subset["utility_loss_norm"] = aggre_data_df_utility_subset["utility_loss"] / aggre_data_df_utility_subset["true_max_utility_loss"]
        aggre_data_df_utility_subset["pmi_gain_norm"] = aggre_data_df_utility_subset["pmi_gain"] / aggre_data_df_utility_subset["max_pmi_gain"]

        utility_mean = aggre_data_df_utility_subset.groupby("capacity").mean()
        utillity_max = aggre_data_df_utility_subset.groupby("capacity").max()
        utility_min = aggre_data_df_utility_subset.groupby("capacity").min()

        x_no_reconstruction = utility_mean["utility_loss_norm"]
        y_no_reconstruction_mean = 1 - utility_mean["pmi_gain_norm"]
        y_no_reconstruction_max = 1 - utillity_max["pmi_gain_norm"]
        y_no_reconstruction_min = 1 - utility_min["pmi_gain_norm"]

        ax.set_title(population)
        ax.set_xlim([0, 1])
        ax.set_xlabel("Utility Loss")
        ax.set_ylabel("Privacy Risk")

                # Find the point with the highest utility loss
        max_utility_idx = utility_mean["utility_loss_norm"].idxmax()
        max_utility_x = utility_mean.loc[max_utility_idx, "utility_loss_norm"]
        max_utility_y = 1 - utility_mean.loc[max_utility_idx, "pmi_gain_norm"]

        # Plot the main result and shaded area
        ax.plot(x_no_reconstruction, y_no_reconstruction_mean, label="No Reconstruction", color=MAIN_COLORS[1], marker=".")
        ax.fill_between(x_no_reconstruction, y_no_reconstruction_min, y_no_reconstruction_max, color=MAIN_COLORS[1], alpha=0.2)

        # Add dashed line from highest utility loss point to the right (x=1)
        ax.hlines(y=max_utility_y, xmin=max_utility_x, xmax=1, linestyles="dashed", colors=MAIN_COLORS[1])

    for i, population in enumerate(subplots):
        if population != "All":
            aggre_data_df_subset = aggre_data_df_privacy[aggre_data_df_privacy["population"] == population]
        else:
            aggre_data_df_subset = aggre_data_df_privacy
        
        aggre_data_df_subset = aggre_data_df_subset.select_dtypes(include=[np.number])
        
        ax = axes[1, i]

        aggre_data_df_subset["utility_loss_norm"] = aggre_data_df_subset["utility_loss"] / aggre_data_df_subset["true_max_utility_loss"]
        aggre_data_df_subset["pmi_gain_norm"] = aggre_data_df_subset["pmi_gain"] / aggre_data_df_subset["max_pmi_gain"]

        utility_mean = aggre_data_df_subset.groupby("capacity").mean()
        utillity_max = aggre_data_df_subset.groupby("capacity").max()
        utility_min = aggre_data_df_subset.groupby("capacity").min()

        x_no_reconstruction = utility_mean["utility_loss_norm"]
        y_no_reconstruction_mean = 1 - utility_mean["pmi_gain_norm"]
        y_no_reconstruction_max = 1 - utillity_max["pmi_gain_norm"]
        y_no_reconstruction_min = 1 - utility_min["pmi_gain_norm"]

        max_utility_idx = utility_mean["utility_loss_norm"].idxmax()
        max_utility_x = utility_mean.loc[max_utility_idx, "utility_loss_norm"]
        max_utility_y = 1 - utility_mean.loc[max_utility_idx, "pmi_gain_norm"]

        # Plot the main result and shaded area
        ax.plot(x_no_reconstruction, y_no_reconstruction_mean, label="No Reconstruction", color=MAIN_COLORS[1], marker=".")
        ax.fill_between(x_no_reconstruction, y_no_reconstruction_min, y_no_reconstruction_max, color=MAIN_COLORS[1], alpha=0.2)

        # Add dashed line from highest utility loss point to the right (x=1)
        ax.hlines(y=max_utility_y, xmin=max_utility_x, xmax=1, linestyles="dashed", colors=MAIN_COLORS[1])

        ax.set_xlim([0, 1])
        ax.set_title(population)
        ax.set_xlabel("Utility Loss")
        ax.set_ylabel("Privacy Risk")

    plt.tight_layout()
    plt.savefig(f"{PLOT_OUT_PATH}/M2A_privacy_vs_utility_all.pdf")
    plt.close()

if __name__ == "__main__":
    apply_style()
    plot_figure_A()
