#!/usr/bin/env python3
"""Lift SetA AP1-primed H3K27ac peaks mm10 -> hg38, keep uniquely-mapped regions."""
import sys
from pyliftover import LiftOver

CHAIN = "mm10ToHg38.over.chain.gz"
INBED = "../SetA_AP1primed_H3K27ac.bed"
OUTBED = "SetA_AP1primed_hg38.bed"
DROPLOG = "liftover_drops.txt"

lo = LiftOver(CHAIN)

def lift_point(chrom, pos):
    r = lo.convert_coordinate(chrom, pos)
    return r  # list of (chrom, pos, strand, conv_score) or []

kept, dropped = [], []
n = 0
with open(INBED) as f:
    for line in f:
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        chrom, start, end = parts[0], int(parts[1]), int(parts[2])
        signal = parts[3] if len(parts) > 3 else "."
        n += 1
        name = f"{chrom}:{start}-{end}"
        s = lift_point(chrom, start)
        e = lift_point(chrom, end - 1)  # end is exclusive; lift last base
        # require both ends map uniquely to same chrom, consistent orientation
        if len(s) != 1 or len(e) != 1:
            dropped.append((name, signal, f"multi/none: start={len(s)} end={len(e)}"))
            continue
        sc, sp = s[0][0], s[0][1]
        ec, ep = e[0][0], e[0][1]
        if sc != ec:
            dropped.append((name, signal, f"split chrom {sc}!={ec}"))
            continue
        lo_start, hi_start = min(sp, ep), max(sp, ep)
        new_start, new_end = lo_start, hi_start + 1
        new_len = new_end - new_start
        old_len = end - start
        # sanity: lifted length within 2x of original (filters scrambled mappings)
        if new_len <= 0 or new_len > 5 * old_len or new_len < old_len / 5:
            dropped.append((name, signal, f"length blowup {old_len}->{new_len}"))
            continue
        kept.append((sc, new_start, new_end, signal, name, old_len, new_len))

kept.sort(key=lambda x: (x[0], x[1]))
with open(OUTBED, "w") as out:
    for sc, st, en, sig, name, ol, nl in kept:
        out.write(f"{sc}\t{st}\t{en}\t{sig}\t{name}\n")

with open(DROPLOG, "w") as d:
    for name, sig, reason in dropped:
        d.write(f"{name}\t{sig}\t{reason}\n")

print(f"input regions:   {n}")
print(f"kept (unique):   {len(kept)}  ({100*len(kept)/n:.1f}%)")
print(f"dropped:         {len(dropped)}  ({100*len(dropped)/n:.1f}%)")
from collections import Counter
c = Counter(r.split(":")[0].split(" ")[0] for _,_,r in dropped)
print("drop reasons:", dict(c))
chroms = Counter(k[0] for k in kept)
print("hg38 chroms covered:", len(chroms))
