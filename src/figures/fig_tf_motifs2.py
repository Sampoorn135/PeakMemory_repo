#!/usr/bin/env python3
"""Improved conserved TF-motif figure: family-grouped heatmap + large logos."""
import sys; sys.path.insert(0, "/sessions/happy-blissful-fermi/mnt/Memory_peak/figures_src")
from fig_style import plt, DATASETS, FIGDIR
import numpy as np, math, glob
import pandas as pd, logomaker, seaborn as sns
from pathlib import Path
from matplotlib.colors import LinearSegmentedColormap
# clean white -> warm red gradient, tops out at a vivid red (no black/dark)
ROCKET = LinearSegmentedColormap.from_list(
    "memheat", ["#FFFFFF", "#FFE8CC", "#FDB87E", "#F47A55", "#E64B35", "#C0392B"])

MOT = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak/cross_dataset_comparison/motifs")

def parse_known(path):
    best = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9: continue
            tf = f[0].split("(")[0].split("/")[0].strip()
            try: logp = float(f[3]); q = float(f[4])
            except ValueError: continue
            nlp = -logp / math.log(10)
            if tf not in best or nlp > best[tf][0]:
                best[tf] = (nlp, q)
    return best

known = {d: parse_known(MOT / d / "knownResults.txt") for d in DATASETS}

def load_pwms(md):
    pwm = {}
    for f in glob.glob(str(md / "known*.motif")):
        rows = []; name = None
        for ln in open(f):
            if ln.startswith(">"):
                if name is None:
                    p = ln.rstrip("\n").split("\t")
                    name = p[1].split("(")[0].split("/")[0].strip() if len(p) > 1 else None
            else:
                v = ln.split()
                if len(v) >= 4: rows.append([float(x) for x in v[:4]])
        if name and name not in pwm and rows: pwm[name] = np.array(rows)
    return pwm
pwm = load_pwms(MOT / "IMQ" / "knownResults")

FAM = [
    ("AP-1 / bZIP", "#E64B35", ["Fos", "Fra1", "JunB", "BATF", "Atf3", "NFAT:AP1"]),
    ("bZIP (CEBP/Maf/Nrf2)", "#F39B7F", ["CEBP", "Nrf2", "MafB"]),
    ("KLF / Sp", "#00A087", ["Klf4", "KLF5"]),
    ("TEAD (Hippo)", "#3C5488", ["TEAD"]),
    ("Forkhead", "#8491B4", ["Foxo1", "FOXK1"]),
    ("STAT / NF1 / CTCF / Homeobox", "#7E6148", ["Stat3", "NF1", "CTCF", "Pdx1"]),
]
SEL, ROWCOL, FAMOF = [], [], []
for fam, col, tfs in FAM:
    for t in tfs:
        if t in pwm and all(t in known[d] for d in DATASETS):
            SEL.append(t); ROWCOL.append(col); FAMOF.append(fam)

H = np.array([[known[d][t][0] for d in DATASETS] for t in SEL])
QQ = np.array([[known[d][t][1] for d in DATASETS] for t in SEL])
VMAX = 60; N = len(SEL)

fig = plt.figure(figsize=(8.2, 0.46 * N + 1.6))
gs = fig.add_gridspec(N, 2, width_ratios=[1.5, 2.3], wspace=0.04, hspace=0.32,
                      left=0.17, right=0.97, top=0.88, bottom=0.13)
heat = fig.add_subplot(gs[:, 0])
im = heat.imshow(np.clip(H, 0, VMAX), cmap=ROCKET, aspect="auto", vmin=0, vmax=VMAX)
heat.set_xticks(range(4)); heat.set_xticklabels(DATASETS, rotation=30, ha="right", fontsize=8)
heat.set_yticks(range(N))
heat.set_yticklabels(SEL, fontsize=8)
for tick, c in zip(heat.get_yticklabels(), ROWCOL):
    tick.set_color(c); tick.set_fontweight("bold")
for i in range(N):
    for j in range(4):
        if QQ[i, j] < 0.05:
            dotc = "white" if min(H[i, j], VMAX) > VMAX * 0.45 else "#3A0A28"
            heat.text(j, i, "•", ha="center", va="center", color=dotc, fontsize=7)
# family separators + labels on far left
seen = {}
for i, fam in enumerate(FAMOF):
    seen.setdefault(fam, []).append(i)
for fam, idxs in seen.items():
    top = min(idxs) - 0.5
    if top > -0.5:
        heat.axhline(top, color="white", lw=1.6)
heat.set_title("Motif enrichment\n(−log$_{10}$ p; • = q<0.05)", fontsize=8.8)
# family colour legend (bottom)
from matplotlib.patches import Patch
fam_handles = [Patch(facecolor=col, edgecolor="none", label=fam) for fam, col, _ in FAM]
fig.legend(handles=fam_handles, title="TF family", loc="lower center",
           bbox_to_anchor=(0.55, -0.02), ncol=3, fontsize=6.5, title_fontsize=7,
           handlelength=1.0, columnspacing=1.2)

for i, tf in enumerate(SEL):
    lax = fig.add_subplot(gs[i, 1])
    mat = pwm[tf]; mat = mat / mat.sum(axis=1, keepdims=True)
    info = logomaker.transform_matrix(pd.DataFrame(mat, columns=list("ACGT")),
                                      from_type="probability", to_type="information")
    logomaker.Logo(info, ax=lax, color_scheme="classic", show_spines=False)
    lax.set_xticks([]); lax.set_yticks([]); lax.set_ylim(0, 2)
    if i == 0: lax.set_title("Motif (bits)", fontsize=8.8)

cax = fig.add_axes([0.30, 0.045, 0.28, 0.016])
cb = fig.colorbar(im, cax=cax, orientation="horizontal")
cb.set_label("−log$_{10}$ p (capped 60)", fontsize=6.6); cb.ax.tick_params(labelsize=6)
fig.suptitle("Conserved TF-motif enrichment in memory peaks across four datasets",
             fontsize=10.5, y=0.975)
for ext in ("png", "pdf"):
    fig.savefig(FIGDIR / f"Fig4_TF_motifs.{ext}", dpi=300, bbox_inches="tight")
print("wrote", FIGDIR / "Fig4_TF_motifs.png", "| TFs:", N)
plt.close(fig)
