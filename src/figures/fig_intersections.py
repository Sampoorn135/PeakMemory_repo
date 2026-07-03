#!/usr/bin/env python3
"""Venn intersection diagrams for the four datasets — genes and enriched TFs."""
import sys; sys.path.insert(0, "/sessions/happy-blissful-fermi/mnt/Memory_peak")
sys.path.insert(0, "/sessions/happy-blissful-fermi/mnt/Memory_peak/figures_src")
from fig_style import plt, DCOLOR, DATASETS, FIGDIR, DLABEL
from pathlib import Path
from venn import venn
from annotate_and_overlap import load_tss, build_gene_loci, load_bed, annotate_set
from matplotlib.colors import to_rgba

BASE = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak")
GPATH = {
    "IMQ": BASE / "memory_analysis/results/04_final_memory",
    "HFSC": BASE / "hfsc/hfsc_analysis/results/04_final_memory",
    "distil": BASE / "distil/distil_analysis/results/04_final_memory",
    "pancreas": BASE / "epithelial_pancreas/pancreas_analysis/results/04_final_memory",
}
MOT = BASE / "cross_dataset_comparison/motifs"
COND = "50_Single"

# ---- gene-locus sets ----
by_chrom = load_tss(BASE / "mm10/mm10.tss")
t2l, locus_info = build_gene_loci(BASE / "mm10/mm10.rna")
gene_sets = {}
for d in DATASETS:
    _, loci = annotate_set(load_bed(GPATH[d] / f"StrictMemory{COND}.bed"), by_chrom, t2l, locus_info)
    gene_sets[d] = loci

# ---- enriched TF sets (q<0.05) ----
def enriched(path):
    s = set()
    with open(path, encoding="utf-8", errors="replace") as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 5: continue
            try: q = float(f[4])
            except ValueError: continue
            if q < 0.05: s.add(f[0].split("(")[0].split("/")[0].strip())
    return s
tf_sets = {d: enriched(MOT / d / "knownResults.txt") for d in DATASETS}

CMAP = [to_rgba(DCOLOR[d], 0.45) for d in DATASETS]

def make_venn(sets, title, fname, core_label):
    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    venn({d: sets[d] for d in DATASETS}, cmap=CMAP, fontsize=8.5,
         legend_loc="upper left", ax=ax)
    # bold the central (four-way) number
    core = len(set.intersection(*[sets[d] for d in DATASETS]))
    ax.set_title(title, fontsize=11, pad=12)
    ax.text(0.5, -0.02, f"four-way core = {core} {core_label}", transform=ax.transAxes,
            ha="center", fontsize=9, fontweight="bold", color="#1A202C")
    # fix legend labels to friendly names
    leg = ax.get_legend()
    if leg:
        for t, d in zip(leg.get_texts(), DATASETS):
            t.set_text(d); t.set_fontsize(8.5)
    for ext in ("png", "pdf"):
        fig.savefig(FIGDIR / f"{fname}.{ext}", dpi=300, bbox_inches="tight")
    print("wrote", FIGDIR / f"{fname}.png", "| core:", core)
    plt.close(fig)

make_venn(gene_sets, f"Memory-gene overlap across four datasets ({COND})",
          "Fig3_gene_venn", "genes")
make_venn(tf_sets, "Enriched-TF overlap across four datasets (q<0.05)",
          "Fig6_TF_venn", "TFs")
