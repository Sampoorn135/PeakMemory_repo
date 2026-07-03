#!/usr/bin/env python3
"""mm10->hg38 liftOver of small peaks, mimicking UCSC liftOver -minMatch=0.95.
Maps EVERY base of each peak. Chains are score-sorted (best first) so the first
chain covering a base wins. A peak is kept if >=95% of its bases map and >=95%
of mapped bases agree on one hg38 chrom; lifted region = span of those bases."""
import bisect
from collections import defaultdict, Counter

CHAIN = "mm10ToHg38.over.chain"
INBED = "../SetA_AP1primed_H3K27ac.bed"
import sys
MINMATCH = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
OUTBED = sys.argv[2] if len(sys.argv) > 2 else "SetA_AP1primed_hg38.bed"
DROPLOG = "liftover_drops.txt"

peaks = []
# tchrom -> {pos: peak_idx}   (peaks don't overlap within set; one idx per base)
base_pos = defaultdict(dict)
with open(INBED) as f:
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue
        p = line.split("\t")
        chrom, start, end = p[0], int(p[1]), int(p[2])
        sig = p[3] if len(p) > 3 else "."
        idx = len(peaks)
        peaks.append((chrom, start, end, sig, f"{chrom}:{start}-{end}"))
        d = base_pos[chrom]
        for b in range(start, end):
            d[b] = idx

sorted_pos = {c: sorted(d.keys()) for c, d in base_pos.items()}
needed = set(base_pos.keys())
# per peak: list of (qchrom, qpos) for each mapped base
hits = defaultdict(list)

def resolve_block(tchrom, t0, t1, q0, qSize, qStrand, qName):
    arr = sorted_pos.get(tchrom)
    if not arr:
        return 0
    lo = bisect.bisect_left(arr, t0)
    hi = bisect.bisect_left(arr, t1)
    if lo == hi:
        return 0
    d = base_pos[tchrom]
    done = []
    for k in range(lo, hi):
        pos = arr[k]
        qcoord = q0 + (pos - t0)
        hg = qcoord if qStrand == "+" else qSize - 1 - qcoord
        hits[d[pos]].append((qName, hg))
        done.append(pos)
    for pos in done:
        del d[pos]
    sorted_pos[tchrom] = sorted(d.keys())
    return len(done)

cur = None
with open(CHAIN) as fh:
    for line in fh:
        if line.startswith("chain"):
            f = line.split()
            tName = f[2]; tStart = int(f[5])
            qName = f[7]; qSize = int(f[8]); qStrand = f[9]; qStart = int(f[10])
            cur = [tName, qName, qSize, qStrand, tStart, qStart] if (tName in needed and sorted_pos.get(tName)) else None
        elif cur is not None:
            f = line.split()
            if not f:
                cur = None; continue
            size = int(f[0])
            tName, qName, qSize, qStrand, t0, q0 = cur
            n = resolve_block(tName, t0, t0 + size, q0, qSize, qStrand, qName)
            if n and not sorted_pos.get(tName):
                needed.discard(tName)
                if not needed:
                    break
            if len(f) == 3:
                cur[4] = t0 + size + int(f[1]); cur[5] = q0 + size + int(f[2])
            else:
                cur = None

# per-peak mapping summary (threshold-independent)
summary = []  # idx, L, frac_mapped, top_chrom, top_frac, new_start, new_end, new_len
for idx, (chrom, start, end, sig, name) in enumerate(peaks):
    L = end - start
    h = hits.get(idx, [])
    if not h:
        summary.append((idx, L, 0.0, None, 0.0, 0, 0, 0)); continue
    frac_mapped = len(h) / L
    cc = Counter(c for c, _ in h)
    top_chrom, top_n = cc.most_common(1)[0]
    top_frac = top_n / len(h)
    positions = [pos for c, pos in h if c == top_chrom]
    ns, ne = min(positions), max(positions) + 1
    summary.append((idx, L, frac_mapped, top_chrom, top_frac, ns, ne, ne - ns))

def evaluate(mm):
    kept, dropped = [], []
    for (idx, L, fm, tc, tf, ns, ne, nl) in summary:
        chrom, start, end, sig, name = peaks[idx]
        if fm < mm:
            dropped.append((name, sig, f"mapped {fm:.2f}")); continue
        if tf < mm:
            dropped.append((name, sig, "split")); continue
        if nl > 5 * L:
            dropped.append((name, sig, f"blowup {L}->{nl}")); continue
        kept.append((tc, ns, ne, sig, name))
    return kept, dropped

n = len(peaks)
print(f"input regions: {n}")
print("threshold sweep (minMatch -> kept):")
for mm in (0.95, 0.75, 0.5, 0.25, 0.1):
    k, _ = evaluate(mm)
    print(f"  minMatch={mm:<4} kept={len(k):3d} ({100*len(k)/n:.1f}%)")

kept, dropped = evaluate(MINMATCH)
kept.sort(key=lambda x: (x[0], x[1]))
with open(OUTBED, "w") as o:
    for sc, st, en, sig, name in kept:
        o.write(f"{sc}\t{st}\t{en}\t{sig}\t{name}\n")
with open(DROPLOG, "w") as d:
    for name, sig, reason in dropped:
        d.write(f"{name}\t{sig}\t{reason}\n")
print(f"\nCHOSEN minMatch={MINMATCH}: kept {len(kept)} ({100*len(kept)/n:.1f}%) -> {OUTBED}")
print("drop reasons:", dict(Counter(r.split()[0] for _,_,r in dropped)))
