#!/usr/bin/env python3
"""Figure 6 — overlap of enriched TFs across the four datasets (UpSet + counts)."""
import sys; sys.path.insert(0, "/sessions/happy-blissful-fermi/mnt/Memory_peak/figures_src")
from fig_style import plt, save, DATASETS, DCOLOR, FIGDIR
from pathlib import Path
from upsetplot import from_contents, UpSet

MOT = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak/cross_dataset_comparison/motifs")
QCUT = 0.05

def enriched(path):
    s = set()
    with open(path, encoding="utf-8", errors="replace") as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 5: continue
            try: q = float(f[4])
            except ValueError: continue
            if q < QCUT:
                s.add(f[0].split("(")[0].split("/")[0].strip())
    return s

sets = {d: enriched(MOT / d / "knownResults.txt") for d in DATASETS}
core = set.intersection(*[sets[d] for d in DATASETS])

data = from_contents({d: sets[d] for d in DATASETS})
fig = plt.figure(figsize=(8.4, 4.6))
up = UpSet(data, sort_by="cardinality", show_counts=True, facecolor="#2D3748",
           shading_color="#EDF2F7")
up.plot(fig=fig)
fig.suptitle(f"Enriched-TF overlap across four datasets (q<0.05)\nfour-way core = {len(core)} TFs (AP-1 / KLF / TEAD / CTCF)",
             fontsize=9.8, y=1.0)
for ext in ("png", "pdf"):
    fig.savefig(FIGDIR / f"Fig6_TF_upset.{ext}", dpi=300, bbox_inches="tight")
print("wrote", FIGDIR / "Fig6_TF_upset.png", "| four-way TF core:", len(core))
plt.close(fig)
