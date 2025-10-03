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
)

GAPSCORE_RESULTS = 2

def plot_combined_figure():
    # Load and merge data
    data_df = load_data_all(GAPSCORE_RESULTS)
    print(data_df)
    
    # Compute privacy gain metrics for each condition:
    # "privacy_loss" already exists for no reconstruction.  
    data_df["privacy_gain"] = data_df["pmi_gain_normalized"]  
    conditions = {
        "PanMixer": {"link_column": "gap_score", "pg_column": "privacy_gain", "shift": 0},
    }
    
    # Set up a figure with two panels (rows)
    import matplotlib.gridspec as gridspec

    fig = plt.figure(figsize=(10, 7))  # overall figure size
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[1, 3])  # adjust as needed

    axA = fig.add_subplot(gs[0])  # This will be vertically squashed (shorter)
    axB = fig.add_subplot(gs[1])  # Taller subplot

    # Define colors for each condition
    cond_colors = {"PanMixer": MAIN_COLORS[1]}

    for cond_name, cond_info in conditions.items():
        link_col = cond_info["link_column"]
        pg_col = cond_info["pg_column"]
        
        # Drop missing values and compute linking outcome as binary:
        df_cond = data_df.dropna(subset=[pg_col, link_col]).copy()
        df_cond["link_outcome"] = (df_cond[link_col] > 0).astype(int)

        df_cond.to_csv(f"{PLOT_OUT_PATH}/M3_GapScores_{cond_name}.csv")
        
        # Scatter plot the data points with a slight random jitter and a vertical shift
        x = 1 - df_cond[pg_col]
        y = df_cond["link_outcome"]
        y_jitter = y + np.random.uniform(-0.05, 0.05, size=len(y)) + cond_info["shift"]
        axA.scatter(x, y_jitter, alpha=0.5, color=cond_colors[cond_name], label=cond_name)

                # Define the transformation function
        def logistic_squashing(p, k=60):
            return 1 / (1 + np.exp(-k * (p - 0.5)))
        
        if x.values.size > 0:
            # If we have both classes, fit logistic regression:
            if len(np.unique(y)) > 1:
                model = LogisticRegression()
                X = x.values.reshape(-1, 1)
                model.fit(X, y)
                
                # Compute F1 score on training data
                y_pred_class = model.predict(X)
                f1 = f1_score(y, y_pred_class)

                coef = model.coef_[0][0]    # β₁
                intercept = model.intercept_[0]  # β₀

                decision_boundary = -intercept / coef

                print(f"Decision boundary for {cond_name}: {decision_boundary}")
                
                # Generate a grid of privacy gain values for plotting over the range of x:
                x_grid = np.linspace(x.min(), x.max(), 200)
                y_pred_prob = model.predict_proba(x_grid.reshape(-1, 1))[:, 1]
                y_pred_prob = logistic_squashing(y_pred_prob)
                
                # Shift the predicted probability curve by the same amount as the datapoints
                axA.plot(x_grid, y_pred_prob + cond_info["shift"],
                        color=cond_colors[cond_name],
                        lw=2,
                       label=f"{cond_name} (F1 = {f1:.2f})")
            else:
                # If only one class is present, plot a horizontal line at the unique (shifted) y value
                unique_y = np.unique(y)[0]
                y_value = unique_y + cond_info["shift"]
                axA.hlines(y_value, x.min(), x.max(),
                        color=cond_colors[cond_name],
                        lw=2,
                        label=f"{cond_name} (F1 = N/A)")
                
    axA.set_xlabel("Privacy Risk")
    axA.set_yticks([0, 1])
    axA.set_yticklabels(["Cannot Link", "Can Link"])
    axA.set_ylabel("Linkability")
    axA.legend(loc="upper left")
    axA.set_xscale("log")
    
    ### Panel B: Grouped Box Plots by Linking Outcome for Each Condition

    box_data = []      # List to hold data for each box plot
    positions = []     # Positions for the boxes on the x-axis
    xtick_labels = []  # Labels for each population and condition

    pos_counter = 0
    condition_spacing = 1.8  # Space between different conditions
    population_spacing = 0.5  # Space between populations within a condition
    box_width = 0.3  # Width of box plots

    # Iterate through each population

    data_df = load_data_all(GAPSCORE_RESULTS)

    data_df["privacy_gain"] = 1 - data_df["pmi_gain_normalized"]  

    for subpop in SUBPOPS:
        if subpop == "All":
            data_df_subset = data_df
        else:
            data_df_subset = data_df[data_df["population"] == subpop]

        # Iterate through each condition
        for cond_name, cond_info in conditions.items():
            link_col = cond_info["link_column"]
            pg_col = cond_info["pg_column"]

            df_cond = data_df_subset.dropna(subset=[pg_col, link_col]).copy()
            df_cond["link_outcome"] = (df_cond[link_col] > 0).astype(int)

            # Separate privacy loss values for the two groups:
            can_link = df_cond[df_cond["link_outcome"] == 1][pg_col]
            cannot_link = df_cond[df_cond["link_outcome"] == 0][pg_col]

            # Assign positions dynamically
            positions.append(pos_counter)  # Can Link
            positions.append(pos_counter + population_spacing)  # Cannot Link
            box_data.append(can_link)
            box_data.append(cannot_link)

            # Label each condition once at the first population group

            pos_counter += 2 * population_spacing  # Move to next population within condition

        pos_counter += condition_spacing  # Space before next condition
        xtick_labels.append(subpop)


    # Create the box plot
    bp = axB.boxplot(box_data, positions=positions, patch_artist=True, widths=box_width, showfliers=False)

    # Color coding: "magenta" for "Can Link" and "teal" for "Cannot Link"
    box_colors = [BAD_COLOR, GOOD_COLOR] * (len(SUBPOPS) * len(conditions))
    for patch, color in zip(bp['boxes'], box_colors):
        patch.set_facecolor(color)

    # Adjust x-axis labels to be centered within each population group
    xticks = [p * (condition_spacing + 2 * population_spacing) + (population_spacing / 2) for p in np.arange(0, len(SUBPOPS))]
    axB.set_xticks(xticks, labels=xtick_labels, rotation=45, ha="right")
    axB.set_ylabel("Privacy Risk")
    #log scale
    axB.set_yscale("log")

    legend_handles = [
        plt.Line2D([0], [0], color=BAD_COLOR, lw=4, label="Can Link"),
        plt.Line2D([0], [0], color=GOOD_COLOR, lw=4, label="Cannot Link")
    ]

    axB.legend(
        handles=legend_handles,
        loc="upper left",                   # anchor point inside the legend box
        bbox_to_anchor=(1.05, 1),           # position legend just outside the plot (right side)
        borderaxespad=0.,
        frameon=False
    )


    plt.tight_layout()
    plt.savefig(f"{PLOT_OUT_PATH}/M3_GapScores_AllPopulations.pdf")
    plt.close()

def plot_logistic_scatter():
    data_df = load_data_all(GAPSCORE_RESULTS)
    data_df["privacy_gain"] = data_df["pmi_gain_normalized"]

    fig, ax = plt.subplots(figsize=(10, 4))  # Adjusted for single panel

    conditions = {
        "PanMixer": {"link_column": "gap_score", "pg_column": "privacy_gain", "shift": 0},
    }
    cond_colors = {"PanMixer": MAIN_COLORS[1]}

    for cond_name, cond_info in conditions.items():
        link_col = cond_info["link_column"]
        pg_col = cond_info["pg_column"]

        df_cond = data_df.dropna(subset=[pg_col, link_col]).copy()
        df_cond["link_outcome"] = (df_cond[link_col] > 0).astype(int)
        df_cond.to_csv(f"{PLOT_OUT_PATH}/M3_GapScores_{cond_name}.csv")

        x = 1 - df_cond[pg_col]
        y = df_cond["link_outcome"]
        y_jitter = y + np.random.uniform(-0.05, 0.05, size=len(y)) + cond_info["shift"]
        ax.scatter(x, y_jitter, alpha=0.5, color=cond_colors[cond_name], label=cond_name)

        def logistic_squashing(p, k=60):
            return 1 / (1 + np.exp(-k * (p - 0.5)))

        if x.values.size > 0:
            if len(np.unique(y)) > 1:
                model = LogisticRegression()
                X = x.values.reshape(-1, 1)
                model.fit(X, y)

                y_pred_class = model.predict(X)
                f1 = f1_score(y, y_pred_class)

                coef = model.coef_[0][0]
                intercept = model.intercept_[0]
                decision_boundary = -intercept / coef

                print(f"Decision boundary for {cond_name}: {decision_boundary}")

                x_grid = np.linspace(x.min(), x.max(), 200)
                y_pred_prob = model.predict_proba(x_grid.reshape(-1, 1))[:, 1]
                y_pred_prob = logistic_squashing(y_pred_prob)

                ax.plot(x_grid, y_pred_prob + cond_info["shift"],
                        color=cond_colors[cond_name], lw=2,
                        label=f"{cond_name} (F1 = {f1:.2f})")
            else:
                unique_y = np.unique(y)[0]
                y_value = unique_y + cond_info["shift"]
                ax.hlines(y_value, x.min(), x.max(),
                          color=cond_colors[cond_name], lw=2,
                          label=f"{cond_name} (F1 = N/A)")

    ax.set_xlabel("Privacy Risk")
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["Cannot Link", "Can Link"])
    ax.set_ylabel("Linkability")
    ax.legend(loc="upper left")
    ax.set_xscale("log")
    plt.tight_layout()
    plt.savefig(f"{PLOT_OUT_PATH}/M3_GapScores_LogisticScatter.pdf")
    plt.close()

def plot_boxplot_by_population():
    data_df = load_data_all(GAPSCORE_RESULTS)
    data_df["privacy_gain"] = 1 - data_df["pmi_gain_normalized"]

    fig, ax = plt.subplots(figsize=(10, 6))  # Adjusted for single panel

    conditions = {
        "PanMixer": {"link_column": "gap_score", "pg_column": "privacy_gain", "shift": 0},
    }

    box_data = []
    positions = []
    xtick_labels = []
    pos_counter = 0
    condition_spacing = 1.8
    population_spacing = 0.5
    box_width = 0.3

    for subpop in SUBPOPS:
        if subpop == "All":
            data_df_subset = data_df
        else:
            data_df_subset = data_df[data_df["population"] == subpop]

        for cond_name, cond_info in conditions.items():
            link_col = cond_info["link_column"]
            pg_col = cond_info["pg_column"]

            df_cond = data_df_subset.dropna(subset=[pg_col, link_col]).copy()
            df_cond["link_outcome"] = (df_cond[link_col] > 0).astype(int)

            can_link = df_cond[df_cond["link_outcome"] == 1][pg_col]
            cannot_link = df_cond[df_cond["link_outcome"] == 0][pg_col]

            positions.append(pos_counter)
            positions.append(pos_counter + population_spacing)
            box_data.append(can_link)
            box_data.append(cannot_link)

            pos_counter += 2 * population_spacing

        pos_counter += condition_spacing
        xtick_labels.append(subpop)

    bp = ax.boxplot(box_data, positions=positions, patch_artist=True, widths=box_width, showfliers=False)
    box_colors = [BAD_COLOR, GOOD_COLOR] * (len(SUBPOPS) * len(conditions))
    for patch, color in zip(bp['boxes'], box_colors):
        patch.set_facecolor(color)

    xticks = [p * (condition_spacing + 2 * population_spacing) + (population_spacing / 2) for p in np.arange(0, len(SUBPOPS))]
    ax.set_xticks(xticks, labels=xtick_labels, rotation=45, ha="right")
    ax.set_ylabel("Privacy Risk")
    ax.set_yscale("log")

    legend_handles = [
        plt.Line2D([0], [0], color=BAD_COLOR, lw=4, label="Can Link"),
        plt.Line2D([0], [0], color=GOOD_COLOR, lw=4, label="Cannot Link")
    ]

    ax.legend(handles=legend_handles,
              loc="upper left",
              bbox_to_anchor=(1.05, 1),
              borderaxespad=0.,
              frameon=False)

    plt.tight_layout()
    plt.savefig(f"{PLOT_OUT_PATH}/M3_GapScores_BoxPlot.pdf")
    plt.close()

def plot_utility_loss_boxplot():
    data_df = load_data_all(GAPSCORE_RESULTS)
    data_df["utility_loss_normalized"] = data_df["utility_loss"] / data_df["true_max_utility_loss"]

    fig, axB = plt.subplots(figsize=(12, 7))

    conditions = {
        "PanMixer": {"link_column": "gap_score", "ul_column": "utility_loss_normalized"},
    }

    box_data = []
    positions = []
    xtick_labels = []

    pos_counter = 0
    condition_spacing = 1.8
    population_spacing = 0.5
    box_width = 0.3

    for subpop in SUBPOPS:
        if subpop == "All":
            data_df_subset = data_df
        else:
            data_df_subset = data_df[data_df["population"] == subpop]

        for cond_name, cond_info in conditions.items():
            link_col = cond_info["link_column"]
            ul_col = cond_info["ul_column"]

            df_cond = data_df_subset.dropna(subset=[ul_col, link_col]).copy()
            df_cond["link_outcome"] = (df_cond[link_col] > 0).astype(int)

            can_link = df_cond[df_cond["link_outcome"] == 1][ul_col]
            cannot_link = df_cond[df_cond["link_outcome"] == 0][ul_col]

            positions.append(pos_counter)
            positions.append(pos_counter + population_spacing)
            box_data.append(can_link)
            box_data.append(cannot_link)

            pos_counter += 2 * population_spacing

        pos_counter += condition_spacing
        xtick_labels.append(subpop)

    bp = axB.boxplot(box_data, positions=positions, patch_artist=True, widths=box_width, showfliers=False)

    box_colors = [BAD_COLOR, GOOD_COLOR] * (len(SUBPOPS) * len(conditions))
    for patch, color in zip(bp['boxes'], box_colors):
        patch.set_facecolor(color)

    xticks = [p * (condition_spacing + 2 * population_spacing) + (population_spacing / 2) for p in np.arange(0, len(SUBPOPS))]
    axB.set_xticks(xticks, labels=xtick_labels, rotation=45, ha="right")
    axB.set_ylabel("Utility Loss")

    legend_handles = [
        plt.Line2D([0], [0], color=BAD_COLOR, lw=4, label="Can Link"),
        plt.Line2D([0], [0], color=GOOD_COLOR, lw=4, label="Cannot Link")
    ]
    axB.legend(handles=legend_handles, loc="upper right")

    plt.tight_layout()
    plt.savefig(f"{PLOT_OUT_PATH}/M3_UtilityLoss_Boxplot.pdf")
    plt.close()

def plot_combo_privacy_utility():
    data_df = load_data_all(GAPSCORE_RESULTS)
    data_df["privacy_gain"] = 1 - data_df["pmi_gain_normalized"]
    data_df["utility_loss_normalized"] = data_df["utility_loss"] / data_df["true_max_utility_loss"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 7), sharey=False)

    conditions = {
        "PanMixer": {
            "link_column": "gap_score",
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
    privacy_data, privacy_pos, xticks = prepare_boxplot_data("pg_column")
    bp1 = ax1.boxplot(privacy_data, positions=privacy_pos, patch_artist=True, widths=0.3, showfliers=False)
    for patch, color in zip(bp1['boxes'], [BAD_COLOR, GOOD_COLOR] * len(privacy_data)):
        patch.set_facecolor(color)

    ax1.set_xticks([np.mean(privacy_pos[i*2:i*2+2]) for i in range(len(SUBPOPS))])
    ax1.set_xticklabels(xticks, rotation=45, ha="right")
    ax1.set_ylabel("Privacy Risk")
    ax1.set_yscale("log")
    #ax1.set_title("Privacy Loss")

    # Plot Utility Loss
    utility_data, utility_pos, _ = prepare_boxplot_data("ul_column")
    bp2 = ax2.boxplot(utility_data, positions=utility_pos, patch_artist=True, widths=0.3, showfliers=False)
    for patch, color in zip(bp2['boxes'], [BAD_COLOR, GOOD_COLOR] * len(utility_data)):
        patch.set_facecolor(color)

    ax2.set_xticks([np.mean(utility_pos[i*2:i*2+2]) for i in range(len(SUBPOPS))])
    ax2.set_xticklabels(xticks, rotation=45, ha="right")
    ax2.set_ylabel("Utility Loss")
    #ax2.set_title("Utility Loss")

    # Shared Legend
    legend_handles = [
        plt.Line2D([0], [0], color=BAD_COLOR, lw=4, label="Can Link"),
        plt.Line2D([0], [0], color=GOOD_COLOR, lw=4, label="Cannot Link")
    ]
    fig.legend(handles=legend_handles,
           loc="lower center",
           bbox_to_anchor=(0.5, -0.04),
           ncol=2,
           frameon=False)

    plt.tight_layout()
    plt.savefig(f"{PLOT_OUT_PATH}/M3_Privacy_Utility_Boxplots_SideBySide.pdf")
    plt.close()


if __name__ == "__main__":
    apply_style()
    plot_combined_figure()
    plot_logistic_scatter()
    plot_boxplot_by_population()
    plot_utility_loss_boxplot()
    plot_combo_privacy_utility()