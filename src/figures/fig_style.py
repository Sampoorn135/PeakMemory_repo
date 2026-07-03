"""Shared Nature-style matplotlib settings + palette for the memory-peak figures."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

FIGDIR = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak/figures")
FIGDIR.mkdir(parents=True, exist_ok=True)

# NPG (Nature) palette for the four datasets
DCOLOR = {
    "IMQ":      "#E64B35",   # red
    "HFSC":     "#4DBBD5",   # cyan
    "distil":   "#00A087",   # teal
    "pancreas": "#3C5488",   # navy
}
DATASETS = ["IMQ", "HFSC", "distil", "pancreas"]
DLABEL = {
    "IMQ": "IMQ\n(skin inflammation)",
    "HFSC": "HFSC\n(hair-follicle wound)",
    "distil": "distil\n(distal wound)",
    "pancreas": "pancreas\n(injury)",
}
SEQ = "#2B6CB0"      # primary blue
ACCENT = "#DD6B20"   # orange
GREY = "#4A5568"

plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.facecolor": "white",
    "figure.facecolor": "white",
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.titleweight": "bold",
    "axes.labelsize": 8,
    "axes.linewidth": 0.8,
    "axes.edgecolor": "#222222",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "legend.frameon": False,
    "axes.titlepad": 6,
})


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(FIGDIR / f"{name}.{ext}")
    print("wrote", FIGDIR / f"{name}.png")
    plt.close(fig)
