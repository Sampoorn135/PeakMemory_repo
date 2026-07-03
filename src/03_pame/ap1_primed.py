#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AP-1-primed regions, layer 1: epidermal-SC (IMQ) memory peaks that contain an
AP-1(bZIP) motif. Scans peak sequences with HOMER's AP-1 PWM + detection
threshold against the mm10 genome FASTA.
"""
import numpy as np
from pathlib import Path

BASE = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak")
PEAKS = BASE / "memory_analysis/results/04_final_memory/StrictMemory50_Single.bed"
MOTIF = BASE / "cross_dataset_comparison/motifs/IMQ/knownResults/known10.motif"  # AP-1(bZIP)
FASTA = BASE / "mm10"
OUT = BASE / "cross_dataset_comparison/ap1_primed"; OUT.mkdir(parents=True, exist_ok=True)
B2I = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}

# ---- load motif (probabilities) + threshold ----
lines = MOTIF.read_text().splitlines()
hdr = lines[0].split("\t")
thr = float(hdr[2]); mname = hdr[1]
pwm = np.array([[float(x) for x in ln.split()] for ln in lines[1:] if ln.strip()])
m = pwm.shape[0]
# log-odds vs uniform background (0.25), pseudocount; col 4 (N) = very negative
llr = np.full((m, 5), -1e9)
llr[:, :4] = np.log((pwm + 1e-3) / 0.25)
# reverse-complement PWM (for the - strand): reverse positions, swap A<->T, C<->G
llr_rc = llr[::-1][:, [3, 2, 1, 0, 4]]

def seq_to_int(s):
    return np.frombuffer(s.translate(TRANS).encode("latin1"), dtype=np.uint8).astype(np.int64)

# translation table: bases -> 0..4 as bytes
import string
TRANS = {ord(k): v for k, v in {"A":0,"C":1,"G":2,"T":3,"a":0,"c":1,"g":2,"t":3}.items()}
TRANS = str.maketrans({chr(i): chr(B2I.get(chr(i).upper(), 4)) for i in range(256)})

def best_score(arr):
    L = len(arr)
    if L < m: return -1e9
    nw = L - m + 1
    f = np.zeros(nw); r = np.zeros(nw)
    for j in range(m):
        f += llr[j, arr[j:j+nw]]
        r += llr_rc[j, arr[j:j+nw]]
    return max(f.max(), r.max())

# ---- load peaks grouped by chrom ----
peaks = {}
order = []
for ln in open(PEAKS):
    c, s, e = ln.split("\t")[:3]; peaks.setdefault(c, []).append((int(s), int(e))); order.append((c, int(s), int(e)))

def load_chrom(c):
    seq = []
    with open(FASTA / f"{c}.fa") as fh:
        next(fh)
        for ln in fh: seq.append(ln.strip())
    return "".join(seq)

ap1 = set(); scores = {}
for c, ivs in peaks.items():
    fa = (FASTA / f"{c}.fa")
    if not fa.exists():
        continue
    chrom = load_chrom(c)
    chrom_i = np.frombuffer(chrom.translate(TRANS).encode("latin1"), dtype=np.uint8).astype(np.int64)
    for s, e in ivs:
        sc = best_score(chrom_i[s:e])
        scores[(c, s, e)] = sc
        if sc >= thr:
            ap1.add((c, s, e))
    print(f"  {c}: {len(ivs)} peaks scanned")

# ---- write outputs ----
with open(OUT / "IMQ_memory_AP1motif.bed", "w") as fh:
    for c, s, e in order:
        if (c, s, e) in ap1:
            fh.write(f"{c}\t{s}\t{e}\t{scores[(c,s,e)]:.2f}\n")
ntot = len(order); nap1 = len(ap1)
print(f"\nMotif: {mname}  (length {m}, HOMER threshold {thr:.2f})")
print(f"IMQ epidermal-SC memory peaks: {ntot}")
print(f"  AP-1 motif-positive (primed, layer 1): {nap1}  ({100*nap1/ntot:.1f}%)")
print("wrote", OUT / "IMQ_memory_AP1motif.bed")
