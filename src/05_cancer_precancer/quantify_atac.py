#!/usr/bin/env python3
"""Quantify ATAC signal (mean0 = sum_coverage/region_len) for each lifted hg38
SetA region across the 13 GSE277274 samples, via remote range-read bigWig.
Parallelizes region queries with a thread pool (each thread = own keep-alive conn).
Checkpoints per sample."""
import os, csv, time, threading
from concurrent.futures import ThreadPoolExecutor
from bwremote import RemoteBigWig

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8517nnn/{gsm}/suppl/{fn}"
BED = "SetA_AP1primed_hg38.bed"
OUT = "atac_signal_long.csv"
WORKERS = 8

SAMPLES = [
    ("GSM8517614", "GSM8517614_rmdup_mapq30_AK_01.bwa.bw",     "AK_P6",     "AK"),
    ("GSM8517615", "GSM8517615_rmdup_mapq30_AK_02.bwa.bw",     "AK_P7",     "AK"),
    ("GSM8517616", "GSM8517616_rmdup_mapq30_AK_03.bwa.bw",     "AK_P8",     "AK"),
    ("GSM8517617", "GSM8517617_rmdup_mapq30_AK_05.bwa.bw",     "AK_P10",    "AK"),
    ("GSM8517618", "GSM8517618_rmdup_mapq30_AK_C01.bwa.bw",    "Normal_P6", "Normal"),
    ("GSM8517619", "GSM8517619_rmdup_mapq30_AK_C02.bwa.bw",    "Normal_P7", "Normal"),
    ("GSM8517620", "GSM8517620_rmdup_mapq30_AK_C03.bwa.bw",    "Normal_P8", "Normal"),
    ("GSM8517621", "GSM8517621_rmdup_mapq30_SCC_01_new.bwa.bw","cSCC_P1",   "cSCC"),
    ("GSM8517622", "GSM8517622_rmdup_mapq30_SCC_04_new.bwa.bw","cSCC_P4",   "cSCC"),
    ("GSM8517623", "GSM8517623_rmdup_mapq30_SCC_02.bwa.bw",    "cSCC_P2",   "cSCC"),
    ("GSM8517624", "GSM8517624_rmdup_mapq30_SCC_05.bwa.bw",    "cSCC_P5",   "cSCC"),
    ("GSM8517625", "GSM8517625_rmdup_mapq30_SCC_C01.bwa.bw",   "Normal_P1", "Normal"),
    ("GSM8517626", "GSM8517626_rmdup_mapq30_SCC_C02.bwa.bw",   "Normal_P2", "Normal"),
]

regions = []
with open(BED) as f:
    for line in f:
        c, s, e, sig, name = line.rstrip("\n").split("\t")
        regions.append((c, int(s), int(e), sig, name))

done = set()   # (gsm, region_name) already written
if os.path.exists(OUT):
    with open(OUT) as f:
        r = csv.reader(f); next(r, None)
        for row in r:
            done.add((row[0], row[3]))

write_header = not os.path.exists(OUT)
out = open(OUT, "a", newline="")
w = csv.writer(out)
if write_header:
    w.writerow(["gsm", "sample", "state", "region", "chrom", "start", "end", "signal", "mean0", "covfrac"])
wlock = threading.Lock()

def process_sample(gsm, fn, title, state):
    url = BASE.format(gsm=gsm, fn=fn)
    tls = threading.local()
    todo = [reg for reg in regions if (gsm, reg[4]) not in done]
    if not todo:
        return 0
    def get_rb():
        rb = getattr(tls, "rb", None)
        if rb is None:
            rb = RemoteBigWig(url, timeout=30, retries=8)
            tls.rb = rb
        return rb
    def query(reg):
        c, s, e, sig, name = reg
        rb = get_rb()
        res = rb.region_stats(c, s, e)
        if res is None:
            row = [gsm, title, state, name, c, s, e, sig, "", ""]
        else:
            ssum, cov, L = res
            row = [gsm, title, state, name, c, s, e, sig, ssum / L, cov / L]
        with wlock:
            w.writerow(row); out.flush()
        return 1
    n = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        for _ in ex.map(query, todo):
            n += 1
    return n

for gsm, fn, title, state in SAMPLES:
    t0 = time.time()
    try:
        n = process_sample(gsm, fn, title, state)
    except Exception as e:
        print(f"PARTIAL {gsm}: {e}", flush=True); continue
    if n:
        print(f"DONE {gsm} {title} ({state}) +{n} regions in {time.time()-t0:.0f}s", flush=True)
    else:
        print(f"skip {gsm} (already complete)", flush=True)

out.close()
