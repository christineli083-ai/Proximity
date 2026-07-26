import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"
OUT_DIR = Path(__file__).parent.parent / "website" / "assets"
plt.rcParams.update({"font.size": 11, "figure.dpi": 150})

coords = pd.read_csv(DATA_DIR / "synthetic_cell_coordinates.csv")
cases = {"48": "infiltrated", "1": "intermediate", "45": "excluded"}

for pid_str, tag in cases.items():
    pid = float(pid_str)
    g = coords[coords.patient_id == pid]
    tumor = g[g.cell_type == "melanocyte"]
    tcell = g[g.cell_type == "CD8_Tcell"]

    fig, ax = plt.subplots(figsize=(5, 5))
    hb = ax.hexbin(tcell.x_um, tcell.y_um, gridsize=22, cmap="magma", mincnt=1)
    ax.scatter(tumor.x_um, tumor.y_um, c="#3FE0D0", s=6, alpha=0.55, edgecolors="none", label="Melanocyte")
    ax.set_xlabel("x (μm)"); ax.set_ylabel("y (μm)")
    ax.set_aspect("equal")
    ax.set_facecolor("#0d1117")
    fig.patch.set_facecolor("#0d1117")
    ax.tick_params(colors="#9aa4b2")
    ax.xaxis.label.set_color("#9aa4b2"); ax.yaxis.label.set_color("#9aa4b2")
    fig.tight_layout()
    fig.savefig(OUT_DIR / f"case_{tag}.png", bbox_inches="tight", facecolor="#0d1117")
    plt.close(fig)

print("done")
