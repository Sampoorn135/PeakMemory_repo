#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pure-Python (standard-library only) H3K27ac bigWig reader  ->  Set A.

NO pip installs of any kind. Runs on the system python3, including 3.14.
Reads the mean H3K27ac signal over the AP-1-primed Set B regions (467) and a
random background (2000), calls H3K27ac-active regions, and writes Set A.

Driven by RunSetA.app.  Manual use:
    python3 extract_h3k27ac.py <bigWig path or https URL>
Inputs in CWD : SetB_named.bed, background_named.bed
Outputs in CWD: setB_H3K27ac_signal.tab, background_H3K27ac_signal.tab,
                SetA_AP1primed_H3K27ac.bed
"""
import sys, struct, zlib
try:
    sys.stdout.reconfigure(line_buffering=True)   # live progress in run.log
except Exception:
    pass

BW   = sys.argv[1] if len(sys.argv) > 1 else "GSE86900_Mouse_H3K27ac.bigWig"
SETB = "SetB_named.bed"
BG   = "background_named.bed"

BIGWIG_MAGIC = 0x888FFC26

# ----------------------------------------------------------------- byte readers
class LocalReader:
    is_local = True
    def __init__(self, path):
        self.f = open(path, "rb")
        self.f.seek(0, 2); self.size = self.f.tell()
    def preload(self, a, b):  # no-op for local (file seeks are cheap)
        pass
    def read(self, off, size):
        self.f.seek(off); return self.f.read(size)

class RemoteReader:
    """HTTP(S) range reader that REUSES one keep-alive connection (fast)."""
    is_local = False
    def __init__(self, url):
        import http.client
        from urllib.parse import urlparse
        self._httplib = http.client
        u = urlparse(url)
        self.scheme = u.scheme
        self.host = u.netloc
        self.path = u.path + (("?" + u.query) if u.query else "")
        self._conn = None
        self._total = None
        self._blobs = []
        self._fetch(0, 1)                 # prime connection + learn file size
        self.size = self._total
    def _connect(self):
        self._conn = (self._httplib.HTTPSConnection if self.scheme == "https"
                      else self._httplib.HTTPConnection)(self.host, timeout=120)
    def _fetch(self, off, size):
        last = off + size - 1
        for _ in range(4):                # retry/reconnect on a dropped keep-alive
            try:
                if self._conn is None:
                    self._connect()
                self._conn.request("GET", self.path, headers={"Range": f"bytes={off}-{last}"})
                resp = self._conn.getresponse()
                data = resp.read()         # fully read so the connection can be reused
                cr = resp.getheader("Content-Range")
                if cr and "/" in cr:
                    self._total = int(cr.split("/")[1])
                if resp.status in (200, 206):
                    return data
            except Exception:
                pass
            try:
                if self._conn:
                    self._conn.close()
            except Exception:
                pass
            self._conn = None
        raise IOError(f"range fetch failed at {off}-{last}")
    def preload(self, a, b):
        self._blobs.append((a, self._fetch(a, b - a)))
    def read(self, off, size):
        for s, b in self._blobs:
            if s <= off and off + size <= s + len(b):
                return b[off - s: off - s + size]
        return self._fetch(off, size)

# ----------------------------------------------------------------- BBI parsing
def read_header(r):
    d = r.read(0, 64)
    if struct.unpack("<I", d[:4])[0] == BIGWIG_MAGIC:
        e = "<"
    elif struct.unpack(">I", d[:4])[0] == BIGWIG_MAGIC:
        e = ">"
    else:
        sys.exit("ERROR: not a bigWig file (bad magic).")
    h = {"e": e}
    h["version"], h["zoomLevels"] = struct.unpack(e + "HH", d[4:8])
    h["chromTreeOffset"], h["fullDataOffset"], h["fullIndexOffset"] = struct.unpack(e + "QQQ", d[8:32])
    h["uncompressBufSize"] = struct.unpack(e + "I", d[52:56])[0]
    # bound the main R-tree index: it ends where the first zoom-level data begins
    h["indexEnd"] = None
    zl = h["zoomLevels"]
    if zl > 0:
        zd = r.read(64, zl * 24)
        zoom_data = []
        for i in range(zl):
            _rl, _res, dOff, _iOff = struct.unpack(e + "IIQQ", zd[i * 24:(i + 1) * 24])
            if dOff > h["fullIndexOffset"]:
                zoom_data.append(dOff)
        if zoom_data:
            h["indexEnd"] = min(zoom_data)
    return h

def read_chrom_tree(r, h):
    e = h["e"]; off = h["chromTreeOffset"]
    head = r.read(off, 32)
    _magic, _blockSize, keySize, _valSize = struct.unpack(e + "IIII", head[:16])
    chroms = {}
    def node(node_off):
        isLeaf, _res, count = struct.unpack(e + "BBH", r.read(node_off, 4))
        if isLeaf:
            isz = keySize + 8
            data = r.read(node_off + 4, count * isz)
            for i in range(count):
                it = data[i * isz:(i + 1) * isz]
                key = it[:keySize].split(b"\x00", 1)[0].decode()
                cid, csize = struct.unpack(e + "II", it[keySize:keySize + 8])
                chroms[key] = (cid, csize)
        else:
            isz = keySize + 8
            data = r.read(node_off + 4, count * isz)
            kids = [struct.unpack(e + "Q", data[i * isz + keySize:i * isz + keySize + 8])[0]
                    for i in range(count)]
            for k in kids:
                node(k)
    node(off + 32)
    return chroms

def find_blocks(r, h, cid, qs, qe):
    """Return list of (dataOffset,dataSize) leaf blocks overlapping cid:[qs,qe)."""
    e = h["e"]; off = h["fullIndexOffset"]
    blocks = []
    def before(eCix, eB): return (eCix < cid) or (eCix == cid and eB <= qs)
    def after(sCix, sB):  return (sCix > cid) or (sCix == cid and sB >= qe)
    def node(node_off):
        isLeaf, _res, count = struct.unpack(e + "BBH", r.read(node_off, 4))
        if isLeaf:
            isz = 32
            data = r.read(node_off + 4, count * isz)
            for i in range(count):
                it = data[i * isz:(i + 1) * isz]
                sC, sB, eC, eB = struct.unpack(e + "IIII", it[:16])
                if not before(eC, eB) and not after(sC, sB):
                    dOff, dSz = struct.unpack(e + "QQ", it[16:32])
                    blocks.append((dOff, dSz))
        else:
            isz = 24
            data = r.read(node_off + 4, count * isz)
            kids = []
            for i in range(count):
                it = data[i * isz:(i + 1) * isz]
                sC, sB, eC, eB = struct.unpack(e + "IIII", it[:16])
                if not before(eC, eB) and not after(sC, sB):
                    kids.append(struct.unpack(e + "Q", it[16:24])[0])
            for k in kids:
                node(k)
    node(off + 48)
    return blocks

def section_sum(buf, e, cid, qs, qe):
    """Sum value*overlap_bp over all wig sections in a decompressed block."""
    p, n, s = 0, len(buf), 0.0
    while p + 24 <= n:
        cId, cStart, _cEnd, itemStep, itemSpan, typ, _res, itemCount = struct.unpack(e + "IIIIIBBH", buf[p:p + 24])
        p += 24
        if typ == 1:      # bedGraph
            if cId == cid:
                for i in range(itemCount):
                    st, en, val = struct.unpack(e + "IIf", buf[p + i * 12:p + i * 12 + 12])
                    ov = min(en, qe) - max(st, qs)
                    if ov > 0: s += val * ov
            p += itemCount * 12
        elif typ == 2:    # varStep
            if cId == cid:
                for i in range(itemCount):
                    st, val = struct.unpack(e + "If", buf[p + i * 8:p + i * 8 + 8])
                    en = st + itemSpan
                    ov = min(en, qe) - max(st, qs)
                    if ov > 0: s += val * ov
            p += itemCount * 8
        elif typ == 3:    # fixedStep
            if cId == cid:
                for i in range(itemCount):
                    (val,) = struct.unpack(e + "f", buf[p + i * 4:p + i * 4 + 4])
                    st = cStart + i * itemStep; en = st + itemSpan
                    ov = min(en, qe) - max(st, qs)
                    if ov > 0: s += val * ov
            p += itemCount * 4
        else:
            break
    return s

def fetch_raw(reader, blocks):
    """blocks: set of (off,size). Returns {off: compressed_bytes}. Merges remote ranges."""
    out = {}
    if reader.is_local:
        for off, size in blocks:
            out[off] = reader.read(off, size)
        return out
    bl = sorted(blocks)
    merged = []
    for off, size in bl:
        if merged and off <= merged[-1][1] + 65536:
            merged[-1][1] = max(merged[-1][1], off + size)
        else:
            merged.append([off, off + size])
    print(f"  fetching {len(blocks)} data blocks in {len(merged)} ranges (keep-alive)...", flush=True)
    buffers = []
    for i, (s, e2) in enumerate(merged, 1):
        buffers.append((s, reader.read(s, e2 - s)))
        if i % 200 == 0:
            print(f"    {i}/{len(merged)} ranges fetched", flush=True)
    for off, size in blocks:
        for s, buf in buffers:
            if s <= off and off + size <= s + len(buf):
                out[off] = buf[off - s:off - s + size]; break
    return out

def pct(vals, p):
    xs = sorted(vals)
    if not xs: return 0.0
    k = (len(xs) - 1) * p / 100.0
    f = int(k); c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)

# ----------------------------------------------------------------- main
def load_bed(path):
    rows = []
    for ln in open(path):
        q = ln.split()
        if len(q) >= 4:
            rows.append((q[3], q[0], int(q[1]), int(q[2])))
    return rows

def main():
    reader = RemoteReader(BW) if BW.lower().startswith(("http://", "https://")) else LocalReader(BW)
    print("opening:", BW, "(remote)" if not reader.is_local else "(local)")
    h = read_header(reader)
    if not reader.is_local:
        reader.preload(h["chromTreeOffset"], h["fullDataOffset"])             # chrom tree
        reader.preload(h["fullIndexOffset"], h["indexEnd"] or reader.size)    # main R-tree index only
    chroms = read_chrom_tree(reader, h)
    chr1 = chroms.get("chr1", (None, None))[1]
    print(f"GENOME-BUILD CHECK -> chr1 length = {chr1}  (mm10=195471971, mm9=197195432)")
    if chr1 not in (195471971, None):
        print("   !! NOT mm10 -- tell Claude this number.")

    setb = load_bed(SETB)
    bg   = load_bed(BG)

    # find blocks for every region, collect unique blocks
    allb = set(); reg = {}
    for name, c, s, e2 in setb + bg:
        if c not in chroms:
            reg[(name)] = None; continue
        cid, csize = chroms[c]; e2 = min(e2, csize)
        b = find_blocks(reader, h, cid, s, e2) if e2 > s else []
        reg[name] = (cid, s, e2, b)
        allb.update(b)
    print(f"regions: {len(setb)} Set B + {len(bg)} background ; unique data blocks = {len(allb)}")

    raw = fetch_raw(reader, allb)
    dec = {}
    def get(off):
        if off not in dec:
            d = raw[off]
            dec[off] = zlib.decompress(d) if h["uncompressBufSize"] > 0 else d
        return dec[off]

    def mean_of(name):
        v = reg.get(name)
        if not v: return float("nan")
        cid, s, e2, blks = v
        if e2 <= s: return float("nan")
        tot = 0.0
        for off, _sz in blks:
            tot += section_sum(get(off), h["e"], cid, s, e2)
        return tot / (e2 - s)

    setb_sig = [(name, c, s, e2, mean_of(name)) for name, c, s, e2 in setb]
    bg_sig   = [(name, c, s, e2, mean_of(name)) for name, c, s, e2 in bg]

    bgv = [x[4] for x in bg_sig if x[4] == x[4]]   # drop NaN
    thr = pct(bgv, 95)
    print(f"background H3K27ac 95th-percentile threshold = {thr:.4f}")

    nA = 0
    with open("setB_H3K27ac_signal.tab", "w") as fh:
        fh.write("name\tchrom\tstart\tend\tH3K27ac_mean\tactive\n")
        for name, c, s, e2, val in setb_sig:
            a = int(val == val and val >= thr); nA += a
            fh.write(f"{name}\t{c}\t{s}\t{e2}\t{0.0 if val!=val else val:.4f}\t{a}\n")
    with open("background_H3K27ac_signal.tab", "w") as fh:
        fh.write("name\tchrom\tstart\tend\tH3K27ac_mean\n")
        for name, c, s, e2, val in bg_sig:
            fh.write(f"{name}\t{c}\t{s}\t{e2}\t{0.0 if val!=val else val:.4f}\n")
    with open("SetA_AP1primed_H3K27ac.bed", "w") as fh:
        for name, c, s, e2, val in setb_sig:
            if val == val and val >= thr:
                fh.write(f"{c}\t{s}\t{e2}\t{val:.4f}\n")

    print(f"\nSet A  (AP-1-primed AND H3K27ac-active) = {nA} of {len(setb)} Set B regions")
    print("DONE. Outputs written next to this script.")

if __name__ == "__main__":
    main()
