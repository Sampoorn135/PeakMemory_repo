#!/usr/bin/env python3
"""Minimal remote bigWig reader using HTTP range requests (little-endian bigWig).
Fetches only header + chrom B-tree + R-tree index + data blocks overlapping the
queried regions. Computes per-region coverage sum and covered-base count so callers
can form mean0 (sum/region_len) or mean (sum/covered)."""
import struct, zlib, http.client
from urllib.parse import urlparse

BIGWIG_MAGIC = 0x888FFC26

class RemoteBigWig:
    def __init__(self, url, timeout=30, retries=4):
        self.url = url; self.timeout = timeout; self.retries = retries
        u = urlparse(url)
        self.host = u.netloc; self.path = u.path
        self.https = (u.scheme == 'https')
        self.conn = None
        self._node_cache = {}   # cache R-tree/index node bytes (shared across regions)
        self._rtree_hdr = None
        self._read_header()
        self._read_chrom_tree()

    def _connect(self):
        if self.https:
            self.conn = http.client.HTTPSConnection(self.host, timeout=self.timeout)
        else:
            self.conn = http.client.HTTPConnection(self.host, timeout=self.timeout)

    def _get(self, start, length):
        import time as _t
        end = start + length - 1
        last = None
        for attempt in range(self.retries):
            try:
                if self.conn is None:
                    self._connect()
                self.conn.request('GET', self.path, headers={'Range': f'bytes={start}-{end}'})
                resp = self.conn.getresponse()
                data = resp.read()  # must fully read to reuse connection
                if resp.status not in (200, 206):
                    raise IOError(f'HTTP {resp.status}')
                return data
            except Exception as e:
                last = e
                try:
                    if self.conn: self.conn.close()
                except Exception:
                    pass
                self.conn = None
                _t.sleep(0.5 * (attempt + 1))  # backoff on 503/transient
        raise last

    def close(self):
        try:
            if self.conn: self.conn.close()
        except Exception:
            pass
        self.conn = None

    def _read_header(self):
        d = self._get(0, 64)
        magic = struct.unpack('<I', d[:4])[0]
        if magic != BIGWIG_MAGIC:
            raise ValueError(f'not little-endian bigWig: magic={magic:#x}')
        (self.magic, self.version, self.zoomLevels, self.chromTreeOffset,
         self.fullDataOffset, self.fullIndexOffset, self.fieldCount,
         self.definedFieldCount, self.autoSqlOffset, self.totalSummaryOffset,
         self.uncompressBufSize, _res) = struct.unpack('<IHHQQQHHQQIQ', d)

    def _read_chrom_tree(self):
        # B-tree header (32 bytes) at chromTreeOffset
        hdr = self._get(self.chromTreeOffset, 32)
        magic, blockSize, keySize, valSize, itemCount, _r = struct.unpack('<IIIIQQ', hdr)
        self.keySize = keySize
        # fetch whole tree region (chromTreeOffset -> fullDataOffset) and parse
        span = self.fullDataOffset - self.chromTreeOffset
        blob = self._get(self.chromTreeOffset, span)
        self.chroms = {}     # name -> (id, size)
        self.id2name = {}
        self._parse_btree_node(blob, 32, keySize)  # start after header

    def _parse_btree_node(self, blob, off, keySize):
        isLeaf, _res, count = struct.unpack_from('<BBH', blob, off)
        off += 4
        if isLeaf:
            for _ in range(count):
                key = blob[off:off+keySize].split(b'\x00', 1)[0].decode()
                off += keySize
                chromId, chromSize = struct.unpack_from('<II', blob, off)
                off += 8
                self.chroms[key] = (chromId, chromSize)
                self.id2name[chromId] = key
        else:
            child_offsets = []
            for _ in range(count):
                off += keySize  # key
                childOffset, = struct.unpack_from('<Q', blob, off)
                off += 8
                child_offsets.append(childOffset)
            for co in child_offsets:
                self._parse_btree_node(blob, co - self.chromTreeOffset, keySize)

    # ---- R-tree traversal ----
    def _find_blocks(self, chromId, start, end):
        """Return list of (dataOffset, dataSize) overlapping (chromId,start,end)."""
        if self._rtree_hdr is None:
            self._rtree_hdr = self._get(self.fullIndexOffset, 48)
        blocks = []
        self._rtree_node(self.fullIndexOffset + 48, chromId, start, end, blocks)
        return blocks

    def _overlap(self, sCi, sB, eCi, eB, chromId, start, end):
        # item range [ (sCi,sB), (eCi,eB) ) vs query (chromId,start,end)
        if (eCi < chromId) or (eCi == chromId and eB <= start):
            return False
        if (sCi > chromId) or (sCi == chromId and sB >= end):
            return False
        return True

    def _rtree_node(self, nodeOffset, chromId, start, end, out):
        d = self._node_cache.get(nodeOffset)
        if d is None:
            d = self._get(nodeOffset, 65536)  # node fits comfortably
            self._node_cache[nodeOffset] = d
        isLeaf, _res, count = struct.unpack_from('<BBH', d, 0)
        off = 4
        if isLeaf:
            for _ in range(count):
                sCi, sB, eCi, eB, dataOffset, dataSize = struct.unpack_from('<IIIIQQ', d, off)
                off += 32
                if self._overlap(sCi, sB, eCi, eB, chromId, start, end):
                    out.append((dataOffset, dataSize))
        else:
            children = []
            for _ in range(count):
                sCi, sB, eCi, eB, childOffset = struct.unpack_from('<IIIIQ', d, off)
                off += 24
                if self._overlap(sCi, sB, eCi, eB, chromId, start, end):
                    children.append(childOffset)
            for co in children:
                self._rtree_node(co, chromId, start, end, out)

    def _parse_block(self, raw, chromId, start, end):
        """Sum coverage*overlap and covered bp within (chromId,start,end) from one data block."""
        if self.uncompressBufSize > 0:
            raw = zlib.decompress(raw)
        ssum = 0.0; cov = 0
        pos = 0; n = len(raw)
        while pos + 24 <= n:
            cId, cStart, cEnd, itemStep, itemSpan, typ, _r, itemCount = \
                struct.unpack_from('<IIIIIBBH', raw, pos)
            pos += 24
            if typ == 1:      # bedGraph
                for _ in range(itemCount):
                    s, e, v = struct.unpack_from('<IIf', raw, pos); pos += 12
                    if cId == chromId:
                        o = min(e, end) - max(s, start)
                        if o > 0:
                            ssum += v * o; cov += o
            elif typ == 2:    # varStep
                for _ in range(itemCount):
                    s, v = struct.unpack_from('<If', raw, pos); pos += 8
                    e = s + itemSpan
                    if cId == chromId:
                        o = min(e, end) - max(s, start)
                        if o > 0:
                            ssum += v * o; cov += o
            elif typ == 3:    # fixedStep
                for i in range(itemCount):
                    v, = struct.unpack_from('<f', raw, pos); pos += 4
                    s = cStart + i * itemStep; e = s + itemSpan
                    if cId == chromId:
                        o = min(e, end) - max(s, start)
                        if o > 0:
                            ssum += v * o; cov += o
            else:
                break
        return ssum, cov

    def region_stats(self, chrom, start, end):
        """Return (sum, covered_bp, length). mean0=sum/length, mean=sum/covered."""
        if chrom not in self.chroms:
            return None
        chromId = self.chroms[chrom][0]
        blocks = self._find_blocks(chromId, start, end)
        ssum = 0.0; cov = 0
        seen = set()
        for off, size in blocks:
            if (off, size) in seen:
                continue
            seen.add((off, size))
            raw = self._get(off, size)
            s, c = self._parse_block(raw, chromId, start, end)
            ssum += s; cov += c
        return ssum, cov, end - start
