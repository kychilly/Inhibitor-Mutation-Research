#!/usr/bin/env python3
import json
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

OUTPUT_DIR = Path("results/validation")

def plot_cleaned_correlation():
    csv_path = OUTPUT_DIR / "literature_correlation_summary.csv"
    if not csv_path.exists():
        print(f"[-] Missing {csv_path}")
        return

    df = pd.read_csv(csv_path)

    pearson_r, pearson_p = stats.pearsonr(df["DeltaDeltaG"], df["Log_IC50_Fold"])
    spearman_r, spearman_p = stats.spearmanr(df["DeltaDeltaG"], df["Log_IC50_Fold"])

    sns.set_theme(style="whitegrid", font="sans-serif")
    fig, ax = plt.subplots(figsize=(9, 6), dpi=300)

    palette = {"Gen 1": "#d95f02", "Gen 2": "#7570b3", "Gen 3": "#1b9e77"}

    # Regression line without massive confidence interval blowing out Y-axis
    sns.regplot(
        data=df,
        x="DeltaDeltaG",
        y="Log_IC50_Fold",
        scatter=False,
        ax=ax,
        color="gray",
        ci=None,
        line_kws={"linestyle": "--", "linewidth": 1.5}
    )

    # Scatter points
    sns.scatterplot(
        data=df,
        x="DeltaDeltaG",
        y="Log_IC50_Fold",
        hue="Gen",
        palette=palette,
        s=160,
        edgecolor="black",
        linewidth=1.0,
        ax=ax,
        zorder=5
    )

    # Label positioning offsets to prevent overlap
    offsets = {
        "Erlotinib": (10, -5),
        "Gefitinib": (10, -5),
        "Afatinib": (-65, 5),
        "Dacomitinib": (10, 0),
        "Osimertinib": (10, -5)
    }

    for _, row in df.iterrows():
        drug = row["Drug"]
        dx, dy = offsets.get(drug, (5, 5))
        ax.annotate(
            drug,
            xy=(row["DeltaDeltaG"], row["Log_IC50_Fold"]),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=10,
            fontweight="bold"
        )

    # Adjust Y-limits cleanly around actual data
    ax.set_ylim(-2.5, 4.0)

    ax.set_title(r"Computational $\Delta\Delta G$ vs. Literature Cell-Line Resistance ($\log_{10}$ IC$_{50}$ Ratio)", pad=14, fontweight="bold")
    ax.set_xlabel(r"Predicted Binding Energy Shift $\Delta\Delta G$ (kcal/mol)", fontweight="bold")
    ax.set_ylabel(r"Literature Resistance Ratio $\log_{10}(\text{IC}_{50,\text{T790M}} / \text{IC}_{50,\text{WT}})$", fontweight="bold")

    # Place stats box in top-right corner to avoid covering Erlotinib
    stats_text = f"Pearson $r = {pearson_r:.2f}$ ($p={pearson_p:.2f}$)\nSpearman $\\rho = {spearman_r:.2f}$ ($p={spearman_p:.2f}$)"
    ax.text(0.95, 0.95, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='gray'))

    ax.legend(title="Drug Gen", frameon=True, facecolor="white", edgecolor="gray", loc="lower left")

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "literature_correlation_plot.png", dpi=300)
    plt.close()
    print("[+] Saved clean correlation plot.")

if __name__ == "__main__":
    plot_cleaned_correlation()