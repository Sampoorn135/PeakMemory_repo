#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Subset the aging-atlas SKIN scATAC (GSE288730) to interfollicular epidermal
BASAL cells (the IMQ-study EpdSC equivalent), pseudobulk by age (and sex), and
export SMALL files for downstream AP-1 / memory-peak analysis.

RUN THIS ON YOUR MAC (where the skin H5AD lives):
    pip install scanpy anndata numpy scipy pandas
    python3 analyze_epdsc_aging.py  /path/to/skin.h5ad  ./epdsc_out

It opens the H5AD in backed mode (low RAM), keeps only EpdSCs, sums counts per
age group, and writes:
    epdsc_out/peaks.bed                 - all peak coordinates (BED)
    epdsc_out/epdsc_pseudobulk.csv      - peak x (age[/sex]) count matrix (EpdSC only)
    epdsc_out/epdsc_cell_counts.csv     - #EpdSCs per age/sex
    epdsc_out/age_gained_21vs1.bed      - peaks up in 21mo vs 1mo (CPM log2FC>=1)
    epdsc_out/age_lost_21vs1.bed        - peaks down in 21mo vs 1mo
Drop those small files back to Claude to run AP-1 motif enrichment + overlap
with the IMQ memory peaks.
"""
import sys, re, gzip
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad

CELLTYPE_PATTERN = re.compile(r"interfollicular epidermal basal", re.I)
LOG2FC = 1.0          # |log2 fold-change| threshold for age DARs
MIN_CPM = 5.0         # require this CPM in the higher group

def pick_col(cols, *cands):
    low = {c.lower(): c for c in cols}
    for cand in cands:
        for lc, orig in low.items():
            if cand in lc:
                return orig
    return None

def main(h5ad, outdir):
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    print("opening (backed):", h5ad)
    A = ad.read_h5ad(h5ad, backed="r")
    obs = A.obs
    print("obs columns:", list(obs.columns))
    ct_col  = pick_col(obs.columns, "cell_type", "celltype", "cell.type", "main_cell", "annotation")
    age_col = pick_col(obs.columns, "age")
    sex_col = pick_col(obs.columns, "sex", "gender")
    print(f"using cell-type col = {ct_col} | age col = {age_col} | sex col = {sex_col}")
    if ct_col is None or age_col is None:
        sys.exit("Could not find cell-type/age columns — open A.obs and set them manually.")

    # mask EpdSCs
    ct = obs[ct_col].astype(str)
    mask = ct.str.contains(CELLTYPE_PATTERN)
    idx = np.where(mask.values)[0]
    print(f"EpdSC (interfollicular epidermal basal) cells: {len(idx):,} of {A.n_obs:,}")
    if len(idx) == 0:
        print("No match — unique cell types were:"); print(sorted(ct.unique())[:60]); sys.exit(1)

    sub = obs.iloc[idx]
    groups = sub[age_col].astype(str) if sex_col is None else \
             (sub[age_col].astype(str) + "|" + sub[sex_col].astype(str))
    print("cells per group:\n", groups.value_counts())
    groups.value_counts().rename("n_cells").to_csv(out / "epdsc_cell_counts.csv")

    # peak coordinates from var
    var = A.var
    if {"chr", "start", "end"}.issubset({c.lower() for c in var.columns}):
        cmap = {c.lower(): c for c in var.columns}
        chrom = var[cmap["chr"]].astype(str).values
        start = var[cmap["start"]].astype(int).values
        end   = var[cmap["end"]].astype(int).values
    else:  # parse from var_names like chr1-1000-2000 / chr1:1000-2000
        toks = [re.split(r"[-:_]", str(v)) for v in A.var_names]
        chrom = np.array([t[0] for t in toks])
        start = np.array([int(t[1]) for t in toks]); end = np.array([int(t[2]) for t in toks])
    with open(out / "peaks.bed", "w") as fh:
        for c, s, e in zip(chrom, start, end):
            fh.write(f"{c}\t{s}\t{e}\n")

    # pseudobulk: sum counts per group, streaming in chunks of cells
    glabels = sorted(groups.unique())
    gpos = {g: i for i, g in enumerate(glabels)}
    P = np.zeros((A.n_vars, len(glabels)), dtype=np.float64)
    CH = 20000
    gidx = groups.values
    for a in range(0, len(idx), CH):
        b = min(a + CH, len(idx))
        rows = idx[a:b]
        X = A[rows, :].X
        X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
        for j, g in enumerate(gidx[a:b]):
            P[:, gpos[g]] += np.asarray(X[j].sum(axis=0)).ravel()
        print(f"  pseudobulked {b:,}/{len(idx):,}")
    pdf = pd.DataFrame(P, columns=glabels)
    pdf.insert(0, "peak", [f"{c}:{s}-{e}" for c, s, e in zip(chrom, start, end)])
    pdf.to_csv(out / "epdsc_pseudobulk.csv", index=False)

    # simple age DAR: 21mo vs 1mo (CPM, log2FC), summing across sexes if present
    def age_digits(g):
        return re.sub(r"\D", "", g.split("|")[0])   # "21m"/"21mo"->"21", "1m"->"1"
    def colsum(age_token):
        cols = [g for g in glabels if age_digits(g) == age_token]
        return P[:, [gpos[c] for c in cols]].sum(axis=1) if cols else None
    old = colsum("21"); young = colsum("1")
    if old is not None and young is not None:
        cpm_o = old / max(old.sum(), 1) * 1e6
        cpm_y = young / max(young.sum(), 1) * 1e6
        l2 = np.log2((cpm_o + 1) / (cpm_y + 1))
        gained = (l2 >= LOG2FC) & (cpm_o >= MIN_CPM)
        lost   = (l2 <= -LOG2FC) & (cpm_y >= MIN_CPM)
        for name, m in [("age_gained_21vs1", gained), ("age_lost_21vs1", lost)]:
            with open(out / f"{name}.bed", "w") as fh:
                for c, s, e in zip(chrom[m], start[m], end[m]):
                    fh.write(f"{c}\t{s}\t{e}\n")
            print(f"{name}: {int(m.sum()):,} peaks")
    print("\nDONE ->", out, "\nDrop epdsc_pseudobulk.csv + age_gained/lost BEDs + peaks.bed back to Claude.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit("usage: python3 analyze_epdsc_aging.py <skin.h5ad> <outdir>")
    main(sys.argv[1], sys.argv[2])
