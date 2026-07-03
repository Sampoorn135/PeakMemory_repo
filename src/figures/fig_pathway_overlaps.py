#!/usr/bin/env python3
"""Pathway enrichment lollipops for all pairwise + three-way overlap cores."""
import sys; sys.path.insert(0, "/sessions/happy-blissful-fermi/mnt/Memory_peak/figures_src")
from fig_style import plt, save
import numpy as np
from pathlib import Path

BASE = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak/cross_dataset_comparison/overlap_pathways")
ONTO = [("kegg", "KEGG", "#E64B35"), ("reactome", "Reactome", "#4DBBD5"),
        ("wikipathways", "WikiPathways", "#00A087"),
        ("biological_process", "GO BP", "#3C5488"),
        ("molecular_function", "GO MF", "#8491B4")]
STOP = {"biological regulation", "regulation of biological process", "cellular process",
        "positive regulation of biological process", "negative regulation of biological process",
        "positive regulation of cellular process", "negative regulation of cellular process",
        "binding", "protein binding", "regulation of cellular process",
        "single-organism process", "metabolic process", "regulation of metabolic process",
        "cellular metabolic process", "primary metabolic process",
        "anatomical structure development", "developmental process",
        "multicellular organism development", "system development",
        "regulation of multicellular organismal process", "cell communication"}
GENES = {"IMQ_HFSC": 356, "IMQ_distil": 520, "IMQ_pancreas": 489, "HFSC_distil": 525,
         "HFSC_pancreas": 417, "distil_pancreas": 576, "IMQ_HFSC_distil": 147,
         "IMQ_HFSC_pancreas": 117, "IMQ_distil_pancreas": 181, "HFSC_distil_pancreas": 155}

def parse(path):
    rows = {}
    try:
        for i, ln in enumerate(open(path, encoding="utf-8", errors="replace")):
            if i == 0: continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 11: continue
            t = f[1]
            try: p = float(f[2]); k = int(f[5])
            except ValueError: continue
            if t not in rows or p < rows[t][0]: rows[t] = (p, k)
    except FileNotFoundError:
        return []
    items = sorted(rows.items(), key=lambda kv: kv[1][0]); n = len(items)
    out = [(t, k, min(1.0, p*n/(r+1)), ) for r, (t, (p, k)) in enumerate(items)]
    for i in range(len(out)-2, -1, -1):
        if out[i][2] > out[i+1][2]: out[i] = (out[i][0], out[i][1], out[i+1][2])
    return out

PATHDB = [("kegg", "KEGG", "#E64B35"), ("reactome", "Reactome", "#4DBBD5"),
          ("wikipathways", "WikiPathways", "#00A087")]

def top_terms(case, per_onto=2):
    """Top specific pathways (KEGG/Reactome/WikiPathways), generic GO excluded."""
    out = []
    for fn, lab, col in PATHDB:
        terms = [(t, k, fdr, col) for (t, k, fdr) in parse(BASE / case / f"{fn}.txt")
                 if t.lower() not in STOP and "podnet" not in t.lower()]
        terms.sort(key=lambda x: x[2])
        out += terms[:per_onto]
    out.sort(key=lambda x: x[2])
    return out

def short(t, n=34):
    return t if len(t) <= n else t[:n-1] + "…"

def panel(ax, case):
    terms = top_terms(case)[::-1]
    if not terms:
        ax.text(0.5, 0.5, "no enrichment", ha="center", transform=ax.transAxes); ax.axis("off"); return
    y = np.arange(len(terms))
    vals = [-np.log10(max(t[2], 1e-12)) for t in terms]
    ax.hlines(y, 0, vals, color=[t[3] for t in terms], lw=2.0, alpha=0.55)
    ax.scatter(vals, y, s=[max(22, t[1]*12) for t in terms], c=[t[3] for t in terms],
               edgecolor="#222", linewidth=0.4, zorder=3)
    ax.axvline(-np.log10(0.05), color="#A0AEC0", ls="--", lw=0.8)
    ax.set_yticks(y); ax.set_yticklabels([short(t[0]) for t in terms], fontsize=6.2)
    ax.set_xlim(0, max(vals) * 1.18 + 0.2)
    nm = case.replace("_", " ∩ ")
    ax.set_title(f"{nm}\n({GENES.get(case,'?')} genes)", fontsize=7.6)
    ax.tick_params(axis="x", labelsize=6)

from matplotlib.lines import Line2D
onto_h = [Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=6,
                 markeredgecolor="#222", markeredgewidth=0.4, label=l) for _, l, c in PATHDB]

# --- pairwise (6 panels) ---
PAIRS = ["IMQ_HFSC", "IMQ_distil", "IMQ_pancreas", "HFSC_distil", "HFSC_pancreas", "distil_pancreas"]
fig, axes = plt.subplots(2, 3, figsize=(12.0, 6.0))
for ax, c in zip(axes.flat, PAIRS): panel(ax, c)
fig.suptitle("Pathway / GO enrichment of the pairwise memory-gene overlaps   (−log$_{10}$ FDR; dashed = 0.05; dot ∝ k)", fontsize=10, y=1.0)
fig.legend(handles=onto_h, title="Ontology", loc="lower center", ncol=5,
           bbox_to_anchor=(0.5, -0.04), fontsize=7, title_fontsize=7.5)
fig.tight_layout(rect=[0, 0.02, 1, 0.97])
save(fig, "Fig7_pairwise_pathways")

# --- three-way (4 panels) ---
TRIP = ["IMQ_HFSC_distil", "IMQ_HFSC_pancreas", "IMQ_distil_pancreas", "HFSC_distil_pancreas"]
fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.0))
for ax, c in zip(axes.flat, TRIP): panel(ax, c)
fig.suptitle("Pathway / GO enrichment of the three-way memory-gene overlaps   (−log$_{10}$ FDR; dashed = 0.05; dot ∝ k)", fontsize=10, y=1.0)
fig.legend(handles=onto_h, title="Ontology", loc="lower center", ncol=5,
           bbox_to_anchor=(0.5, -0.05), fontsize=7, title_fontsize=7.5)
fig.tight_layout(rect=[0, 0.03, 1, 0.96])
save(fig, "Fig8_threeway_pathways")
