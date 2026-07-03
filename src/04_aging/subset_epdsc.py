#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Subset the aging-atlas SKIN scATAC peak matrix to ONLY the interfollicular
epidermal BASAL cells (the IMQ-study EpdSC population) and save a small H5AD
that keeps ALL those cells (every age + sex) and ALL peaks + metadata, so it can
be analysed flexibly later.

USAGE (on the machine that downloaded the file):
    pip install anndata scipy numpy pandas
    python3 subset_epdsc.py  GSM8774016_Skin_peak_count.h5ad  EpdSC_skin_peak.h5ad
"""
import sys, re
import anndata as ad

# the cell type to keep (matches "Skin.Interfollicular epidermal basal cells"
# or "Interfollicular epidermal basal cells" — case-insensitive, substring)
PATTERN = re.compile(r"interfollicular epidermal basal", re.I)

def main(inp, outp):
    print("opening (backed, low-RAM):", inp)
    A = ad.read_h5ad(inp, backed="r")
    print(f"total: {A.n_obs:,} cells x {A.n_vars:,} peaks")
    print("obs columns:", list(A.obs.columns))

    # auto-find the cell-type column (the categorical/text column that contains the pattern)
    ct_col = None
    for c in A.obs.columns:
        if str(A.obs[c].dtype) in ("category", "object"):
            if A.obs[c].astype(str).str.contains(PATTERN).any():
                ct_col = c
                break
    if ct_col is None:
        print("\n!! Could not auto-detect the cell-type column. Here are the text columns + values:")
        for c in A.obs.columns:
            if str(A.obs[c].dtype) in ("category", "object"):
                print("  ", c, "->", list(A.obs[c].astype(str).unique())[:40])
        sys.exit("Set ct_col manually in the script to the correct column name, then re-run.")
    print("cell-type column detected:", ct_col)

    mask = A.obs[ct_col].astype(str).str.contains(PATTERN).values
    n = int(mask.sum())
    print(f"interfollicular epidermal basal (EpdSC) cells: {n:,} of {A.n_obs:,}")
    if n == 0:
        sys.exit("No matching cells found — check the label spelling printed above.")

    # load only the subset into memory, then write
    try:
        sub = A[mask].to_memory()
    except Exception:
        sub = A[mask].copy()

    # report breakdown so you know it's balanced across ages
    for c in sub.obs.columns:
        lc = c.lower()
        if "age" in lc or lc in ("sex", "gender"):
            print(f"\n{c} breakdown:")
            print(sub.obs[c].value_counts())

    sub.write(outp)
    print(f"\nWROTE {outp}  ({sub.n_obs:,} cells x {sub.n_vars:,} peaks)")
    print("Send this file back. It keeps all EpdSCs across every age/sex + all peaks.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit("usage: python3 subset_epdsc.py <Skin_peak_count.h5ad> <output_EpdSC.h5ad>")
    main(sys.argv[1], sys.argv[2])
