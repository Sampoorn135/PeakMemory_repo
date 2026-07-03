#!/usr/bin/env python3
"""Figure 5 — conserved TF-motif enrichment across datasets, with sequence logos."""
import sys; sys.path.insert(0, "/sessions/happy-blissful-fermi/mnt/Memory_peak/figures_src")
from fig_style import plt, save, DATASETS, FIGDIR
import numpy as np, math, glob, os
import pandas as pd, logomaker, seaborn as sns
ROCKET = sns.color_palette("rocket_r", as_cmap=True)
from pathlib import Path

MOT = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak/cross_dataset_comparison/motifs")

# ---- parse knownResults per dataset: tf -> (neglog10p, q, pct_t, pct_b) ----
def parse_known(path):
    best = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9: continue
            tf = f[0].split("(")[0].split("/")[0].strip()
            try:
                logp = float(f[3]); q = float(f[4])
                pt = float(f[6].replace("%","")); pb = float(f[8].replace("%",""))
            except ValueError: continue
            nlp = -logp / math.log(10)
            if tf not in best or nlp > best[tf][0]:
                best[tf] = (nlp, q, pt, pb)
    return best

known = {d: parse_known(MOT / d / "knownResults.txt") for d in DATASETS}

# ---- PWM library: tf -> probability matrix (from IMQ known*.motif) ----
def load_pwms(motif_dir):
    pwm = {}
    for f in glob.glob(str(motif_dir / "known*.motif")):
        rows = []
        name = None
        for ln in open(f):
            if ln.startswith(">"):
                if name is None:
                    parts = ln.rstrip("\n").split("\t")
                    name = parts[1].split("(")[0].split("/")[0].strip() if len(parts) > 1 else None
            else:
                v = ln.split()
                if len(v) >= 4:
                    rows.append([float(v[0]), float(v[1]), float(v[2]), float(v[3])])
        if name and name not in pwm and rows:
            pwm[name] = np.array(rows)
    return pwm

pwm = load_pwms(MOT / "IMQ" / "knownResults")

# ---- conserved TFs to display (all in the four-way core), grouped by family ----
SEL = ["Fos", "Fra1", "JunB", "BATF", "Atf3", "NFAT:AP1",   # AP-1 / bZIP
       "CEBP", "Nrf2", "MafB",                               # other bZIP
       "Klf4", "KLF5",                                       # KLF/Sp
       "TEAD",                                               # Hippo
       "Foxo1", "FOXK1",                                     # Forkhead
       "Stat3", "NF1", "CTCF", "Pdx1"]                       # STAT, NF1, CTCF, pancreas
SEL = [t for t in SEL if t in pwm and all(t in known[d] for d in DATASETS)]

# heatmap matrix of -log10 p (capped for color)
H = np.array([[known[d][t][0] for d in DATASETS] for t in SEL])
QQ = np.array([[known[d][t][1] for d in DATASETS] for t in SEL])
VMAX = 60

N = len(SEL)
fig = plt.figure(figsize=(7.0, 0.34 * N + 1.5))
gs = fig.add_gridspec(N, 2, width_ratios=[2.0, 1.5], wspace=0.04, hspace=0.22,
                      left=0.16, right=0.97, top=0.9, bottom=0.1)

heat = fig.add_subplot(gs[:, 0])
im = heat.imshow(np.clip(H, 0, VMAX), cmap=ROCKET, aspect="auto", vmin=0, vmax=VMAX)
heat.set_xticks(range(4)); heat.set_xticklabels(DATASETS, rotation=30, ha="right", fontsize=7.5)
heat.set_yticks(range(N)); heat.set_yticklabels(SEL, fontsize=7.3)
heat.set_xlabel("")
# mark significant cells (q<0.05) with a dot
for i in range(N):
    for j in range(4):
        if QQ[i, j] < 0.05:
            heat.text(j, i, "•", ha="center", va="center", color="white", fontsize=6)
heat.set_title("Motif enrichment\n(−log$_{10}$ p; • = q<0.05)", fontsize=8.3)

# logos on the right, one per row
for i, tf in enumerate(SEL):
    lax = fig.add_subplot(gs[i, 1])
    mat = pwm[tf]
    mat = mat / mat.sum(axis=1, keepdims=True)
    df = pd.DataFrame(mat, columns=["A", "C", "G", "T"])
    info = logomaker.transform_matrix(df, from_type="probability", to_type="information")
    logomaker.Logo(info, ax=lax, color_scheme="classic", show_spines=False)
    lax.set_xticks([]); lax.set_yticks([]); lax.set_ylim(0, 2)
    if i == 0:
        lax.set_title("Motif", fontsize=8.3)

cax = fig.add_axes([0.16, 0.045, 0.3, 0.018])
cb = fig.colorbar(im, cax=cax, orientation="horizontal")
cb.set_label("−log$_{10}$ p (capped 60)", fontsize=6.6); cb.ax.tick_params(labelsize=6)

fig.suptitle("Conserved transcription-factor motif enrichment in memory peaks\n(four datasets; AP-1 / KLF / TEAD core)",
             fontsize=10, y=0.985)
fig.text(0.02, 0.95, "a", fontsize=14, fontweight="bold")
for ext in ("png", "pdf"):
    fig.savefig(FIGDIR / f"Fig5_TF_motif_heatmap.{ext}", dpi=300, bbox_inches="tight")
print("wrote", FIGDIR / "Fig5_TF_motif_heatmap.png", "| TFs shown:", len(SEL))
plt.close(fig)
