"""
Visualizations:
1. Distribution plot: raw scatter of tumor cells (red) vs T-cells (blue) for
   two example patients (thin Breslow depth vs thick) -- shows WHERE cells are.
2. Density heatmap: 2D kernel-density-style histogram of T-cell concentration
   overlaid on tumor boundary -- shows HOW CONCENTRATED T-cells are, which a
   scatter plot alone can hide (same footprint, different clustering).
3. Summary scatter: mean NN distance vs Breslow depth across all 60 patients,
   with regression line and Pearson r annotated.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"
OUT_DIR = Path(__file__).parent.parent / "website" / "assets"
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.size": 11, "figure.dpi": 150})


def plot_distribution_pair(coords: pd.DataFrame, clinical: pd.DataFrame):
    thin_pid = clinical.sort_values("breslow_depth_mm").iloc[2]["patient_id"]
    thick_pid = clinical.sort_values("breslow_depth_mm").iloc[-3]["patient_id"]

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, pid, label in zip(axes, [thin_pid, thick_pid], ["Thin (low Breslow depth)", "Thick (high Breslow depth)"]):
        g = coords[coords.patient_id == pid]
        tumor = g[g.cell_type == "melanocyte"]
        tcell = g[g.cell_type == "CD8_Tcell"]
        depth = clinical[clinical.patient_id == pid]["breslow_depth_mm"].values[0]
        ax.scatter(tumor.x_um, tumor.y_um, c="#c0392b", s=8, alpha=0.6, label="Melanocyte (tumor)")
        ax.scatter(tcell.x_um, tcell.y_um, c="#2980b9", s=8, alpha=0.6, label="CD8+ T-cell")
        ax.set_title(f"{label}\nBreslow depth = {depth:.2f} mm")
        ax.set_xlabel("x (μm)"); ax.set_ylabel("y (μm)")
        ax.set_aspect("equal")
        ax.legend(fontsize=8, loc="upper right")
    fig.suptitle("Spatial Distribution: Tumor vs. T-cell Point Patterns", y=1.02)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "distribution_comparison.png", bbox_inches="tight")
    plt.close(fig)


def plot_density_heatmap(coords: pd.DataFrame, clinical: pd.DataFrame):
    thin_pid = clinical.sort_values("breslow_depth_mm").iloc[2]["patient_id"]
    thick_pid = clinical.sort_values("breslow_depth_mm").iloc[-3]["patient_id"]

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, pid, label in zip(axes, [thin_pid, thick_pid], ["Thin tumor", "Thick tumor"]):
        g = coords[coords.patient_id == pid]
        tumor = g[g.cell_type == "melanocyte"]
        tcell = g[g.cell_type == "CD8_Tcell"]
        depth = clinical[clinical.patient_id == pid]["breslow_depth_mm"].values[0]

        hb = ax.hexbin(tcell.x_um, tcell.y_um, gridsize=25, cmap="viridis", mincnt=1)
        ax.scatter(tumor.x_um, tumor.y_um, c="white", s=3, alpha=0.4, edgecolors="none")
        ax.set_title(f"{label} (Breslow = {depth:.2f} mm)\nT-cell density (color) vs. tumor outline (white dots)")
        ax.set_xlabel("x (μm)"); ax.set_ylabel("y (μm)")
        ax.set_aspect("equal")
        fig.colorbar(hb, ax=ax, label="T-cell density")
    fig.suptitle("Density Heatmap: CD8+ T-cell Concentration", y=1.02)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "density_heatmap.png", bbox_inches="tight")
    plt.close(fig)


def plot_correlation_summary(results: pd.DataFrame):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    ax = axes[0]
    ax.scatter(results.breslow_depth_mm, results.nn_mean_um, alpha=0.7, c="#8e44ad")
    m, b = np.polyfit(results.breslow_depth_mm, results.nn_mean_um, 1)
    xs = np.linspace(results.breslow_depth_mm.min(), results.breslow_depth_mm.max(), 50)
    ax.plot(xs, m * xs + b, color="black", linestyle="--", linewidth=1.5)
    ax.set_xlabel("Breslow depth (mm)")
    ax.set_ylabel("Mean NN distance (μm)")
    ax.set_title("NN Distance vs. Tumor Thickness")

    ax = axes[1]
    for outcome, color, lab in [(1, "#27ae60", "Survived"), (0, "#c0392b", "Did not survive")]:
        sub = results[results.survived_5yr == outcome]
        ax.scatter([outcome] * len(sub) + np.random.normal(0, 0.03, len(sub)), sub.nn_mean_um,
                   alpha=0.6, color=color, label=lab)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Did not survive", "Survived"])
    ax.set_ylabel("Mean NN distance (μm)")
    ax.set_title("NN Distance vs. Simulated 5-yr Survival")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "correlation_summary.png", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    coords = pd.read_csv(DATA_DIR / "synthetic_cell_coordinates.csv")
    clinical = pd.read_csv(DATA_DIR / "synthetic_clinical.csv")
    results = pd.read_csv(Path(__file__).parent / "patient_level_results.csv")

    plot_distribution_pair(coords, clinical)
    plot_density_heatmap(coords, clinical)
    plot_correlation_summary(results)
    print("Saved 3 figures to", OUT_DIR)
