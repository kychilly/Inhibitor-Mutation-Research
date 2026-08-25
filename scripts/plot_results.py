#!/usr/bin/env python3
"""
Reads Vina docking results across all conditions, calculates Mean ΔG and Standard Deviation,
and generates an inverted grouped bar chart (0 at bottom, -10 at top) with 0.2 interval minor
gridlines comparing Wildtype vs T790M mutant binding affinities across drug generations.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from pathlib import Path

# Paths
RESULTS_DIR = Path("data/processed/docking_results")
OUTPUT_DIR = Path("results/figures")

DRUGS_INFO = [
    {"drug": "erlotinib", "label": "Erlotinib", "gen": "Gen 1"},
    {"drug": "gefitinib", "label": "Gefitinib", "gen": "Gen 1"},
    {"drug": "afatinib", "label": "Afatinib", "gen": "Gen 2"},
    {"drug": "dacomitinib", "label": "Dacomitinib", "gen": "Gen 2"},
    {"drug": "osimertinib", "label": "Osimertinib", "gen": "Gen 3"},
]


def collect_data():
    records = []

    for item in DRUGS_INFO:
        drug = item["drug"]
        label = item["label"]
        gen = item["gen"]

        for rec, variant_label in [("wildtype", "Wildtype (WT)"), ("t790m", "T790M Mutant")]:
            pdbqt_dir = RESULTS_DIR / rec / drug / "pdbqt"
            candidates = sorted(pdbqt_dir.glob(f"{rec}_{drug}_seed*.pdbqt"))

            energies = []
            for path in candidates:
                try:
                    with open(path, "r", encoding="utf-8") as f:
                        for line in f:
                            if line.startswith("REMARK VINA RESULT:"):
                                affinity = float(line.split()[3])
                                energies.append(affinity)
                                break
                except (OSError, IndexError, ValueError):
                    continue

            if energies:
                mean_dg = np.mean(energies)
                std_dg = np.std(energies)
            else:
                mean_dg, std_dg = 0.0, 0.0

            records.append({
                "Drug": f"{label}\n({gen})",
                "Drug_Name": label,
                "Generation": gen,
                "Variant": variant_label,
                "Mean_DG": mean_dg,
                "Std_DG": std_dg,
            })

    return pd.DataFrame(records)


def plot_grouped_bars(df):
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Base styling
    sns.set_theme(style="white", font="sans-serif")
    plt.rcParams.update({"font.size": 12, "axes.labelsize": 13, "xtick.labelsize": 11, "ytick.labelsize": 11})

    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

    palette = {"Wildtype (WT)": "#2b5c8f", "T790M Mutant": "#d95f02"}

    x = np.arange(len(DRUGS_INFO))
    width = 0.35

    wt_df = df[df["Variant"] == "Wildtype (WT)"].reset_index(drop=True)
    mut_df = df[df["Variant"] == "T790M Mutant"].reset_index(drop=True)

    rects1 = ax.bar(
        x - width / 2,
        wt_df["Mean_DG"],
        width,
        yerr=wt_df["Std_DG"],
        capsize=5,
        label="Wildtype (WT)",
        color=palette["Wildtype (WT)"],
        edgecolor="black",
        linewidth=0.8,
        alpha=0.9
    )

    rects2 = ax.bar(
        x + width / 2,
        mut_df["Mean_DG"],
        width,
        yerr=mut_df["Std_DG"],
        capsize=5,
        label="T790M Mutant",
        color=palette["T790M Mutant"],
        edgecolor="black",
        linewidth=0.8,
        alpha=0.9
    )

    # Invert Y-axis so 0 is at the bottom and -11 is at the top
    ax.set_ylim(0, -11)

    # Y-axis Ticks & Gridlines
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))

    # Major grid (1.0 kcal/mol intervals)
    ax.grid(visible=True, which="major", axis="y", color="#cccccc", linestyle="-", linewidth=0.8)
    # Minor grid (0.2 kcal/mol intervals)
    ax.grid(visible=True, which="minor", axis="y", color="#e6e6e6", linestyle=":", linewidth=0.5)

    # Formatting
    ax.set_ylabel(r"Predicted Binding Affinity $\Delta G$ (kcal/mol)", fontweight="bold")
    ax.set_title("EGFR Inhibitor Binding Affinities Across Generations (WT vs. T790M)", pad=15, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(wt_df["Drug"])
    ax.legend(title="Receptor Variant", frameon=True, facecolor="white", edgecolor="gray")

    # Direct numeric labels placed inside bars near top
    for rects in [rects1, rects2]:
        for rect in rects:
            height = rect.get_height()
            ax.annotate(
                f"{height:.2f}",
                xy=(rect.get_x() + rect.get_width() / 2, height),
                xytext=(0, -14),
                textcoords="offset points",
                ha="center",
                va="top",
                fontsize=9,
                fontweight="bold",
                color="white"
            )

    plt.tight_layout()

    png_path = OUTPUT_DIR / "binding_affinity_comparison.png"
    svg_path = OUTPUT_DIR / "binding_affinity_comparison.svg"

    plt.savefig(png_path, dpi=300)
    plt.savefig(svg_path, format="svg")
    plt.close()

    print(f"[+] Saved inverted figure with 0.2 minor grid -> {png_path}")


def main():
    df = collect_data()
    if df.empty or df["Mean_DG"].abs().sum() == 0:
        print("[-] Error: No valid docking data found to plot.")
        return

    plot_grouped_bars(df)


if __name__ == "__main__":
    main()