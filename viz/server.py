#!/usr/bin/env python3
"""Local server for the Rule 30 viewer: static files plus a small binary/JSON API.

    python server.py [port]

Data comes from $R30_DATA (default /Volumes/BLIPS/rule30/viz):
    tiles/<run>/   r30gen output (checkpoints, dense block levels, chain list, histograms)
    centre/        centre_prep output for the 10^9-bit centre column

Triangle tiles are T x T pixels at b cells per pixel (b a power of two). Tiles finer than a run's base
level, and pattern-search tiles up to PAT_MAX, are computed on demand by libr30 from the nearest
checkpoint; coarser tiles are read from the precomputed levels.
"""
import ctypes
import json
import math
import os
import sys
import threading
from collections import OrderedDict
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get("R30_DATA", "/Volumes/BLIPS/rule30/viz")
T = 256
PAT_MAX = 512
HBINS = 64

lib = ctypes.CDLL(os.path.join(HERE, "build", "libr30.dylib"))
P = np.ctypeslib.ndpointer
i64 = ctypes.c_int64
lib.r30_window.restype = ctypes.c_int
lib.r30_window.argtypes = [P(np.uint64, flags="C"), i64, i64, i64, i64, i64, i64, i64, i64, i64, ctypes.c_char_p,
                           ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, P(np.float64, flags="C")]
lib.r30c_find.restype = i64
lib.r30c_find.argtypes = [P(np.uint8, flags="C"), i64, i64, i64, ctypes.c_char_p, P(np.int64, flags="C"), i64,
                          P(np.int64, flags="C"), i64]
lib.r30c_runs.restype = None
lib.r30c_runs.argtypes = [P(np.uint8, flags="C"), i64, i64, i64, P(np.int64, flags="C"), P(np.int64, flags="C")]


class LRU:
    def __init__(self, cap):
        self.cap, self.d, self.lock = cap, OrderedDict(), threading.Lock()

    def get(self, k):
        with self.lock:
            if k in self.d:
                self.d.move_to_end(k)
                return self.d[k]
        return None

    def put(self, k, v):
        with self.lock:
            self.d[k] = v
            self.d.move_to_end(k)
            while len(self.d) > self.cap:
                self.d.popitem(last=False)


def ptr(a):
    return a.ctypes.data_as(ctypes.c_void_p) if a is not None else None


class Run:
    def __init__(self, name, path):
        m = json.load(open(os.path.join(path, "meta.json")))
        self.name, self.path, self.meta = name, path, m
        self.N, self.K, self.MID = m["N"], m["K"], m["MID"]
        self.levels = {}
        for L in m["levels"]:
            mm = np.memmap(os.path.join(path, f"level_{L['B']}.bin"), dtype=np.uint32, mode="r",
                           shape=(L["NR"], L["NC"], 5))
            self.levels[L["B"]] = (L, mm)
        self.base = min(self.levels)
        self.HB = m["HB"]
        self.hist = np.memmap(os.path.join(path, f"hist_{self.HB}.u16"), dtype=np.uint16, mode="r",
                              shape=(m["HR"], m["HC"], HBINS))
        self.chains = np.fromfile(os.path.join(path, "chains.i32"), dtype=np.int32).reshape(-1, 3)
        self.ck = LRU(64)
        self.tiles = LRU(256)  # block tiles are 1.5 MB each
        self._diag = None

    def info(self):
        return {"name": self.name, "N": self.N, "K": self.K, "base": self.base, "LMIN": self.meta["LMIN"],
                "levels": sorted(self.levels), "chains_logged": int(len(self.chains)), "hist": self.meta["hist"],
                "pat_max": PAT_MAX}

    def checkpoint(self, c):
        a = self.ck.get(c)
        if a is None:
            a = np.fromfile(os.path.join(self.path, "ck", f"ck_{c}.u64"), dtype=np.uint64)
            self.ck.put(c, a)
        return a

    def window(self, n0, n1, s0, s1, b, pat):
        c = (max(0, n0 - 128) // self.K) * self.K
        ck = self.checkpoint(c)
        H, W = n1 - n0, s1 - s0
        hist = np.zeros(HBINS, np.float64)
        if b == 1:
            cell, aux, blk = np.zeros(H * W, np.uint8), np.zeros(H * W, np.uint8), None
        else:
            cell = aux = None
            blk = np.zeros(((H + b - 1) // b) * ((W + b - 1) // b) * 6, np.uint32)
        rc = lib.r30_window(ck, (self.MID - c) >> 6, len(ck), self.MID, c, n0, n1, s0, s1, b,
                            pat.encode() if pat else None, ptr(cell), ptr(aux), ptr(blk), hist)
        if rc:
            raise RuntimeError(f"r30_window rc={rc}")
        return cell, aux, blk, hist

    def tile(self, b, ty, tx, pat):
        key = (b, ty, tx, pat)
        hit = self.tiles.get(key)
        if hit is not None:
            return hit
        n0, s0 = ty * T * b, tx * T * b
        n1, s1 = min(n0 + T * b, self.N), s0 + T * b
        hist = np.zeros(HBINS, np.float64)
        flags = 1 if pat else 0
        if n0 >= self.N or ty < 0 or s0 > n1 or s1 <= -n1:
            kind, body = (1, np.zeros(2 * T * T, np.uint8)) if b == 1 else (2, np.zeros(T * T * 6, np.uint32))
        elif b == 1:
            cell, aux, _, hist = self.window(n0, n1, s0, s1, 1, pat)
            full_c, full_a = np.zeros((T, T), np.uint8), np.zeros((T, T), np.uint8)
            full_c[: n1 - n0] = cell.reshape(n1 - n0, T)
            full_a[: n1 - n0] = aux.reshape(n1 - n0, T)
            kind, body = 1, np.concatenate([full_c.ravel(), full_a.ravel()])
        elif b < self.base or (pat and b <= PAT_MAX):
            _, _, blk, hist = self.window(n0, n1, s0, s1, b, pat)
            full = np.zeros((T, T, 6), np.uint32)
            rows = (n1 - n0 + b - 1) // b
            full[:rows] = blk.reshape(rows, T, 6)
            kind, body = 2, full.ravel()
        else:
            flags = 0
            L, mm = self.levels[b]
            full = np.zeros((T, T, 6), np.uint32)
            R0, C0 = ty * T, tx * T - L["CMIN"]
            r0, r1 = max(R0, 0), min(R0 + T, L["NR"])
            c0, c1 = max(C0, 0), min(C0 + T, L["NC"])
            if r1 > r0 and c1 > c0:
                full[r0 - R0:r1 - R0, c0 - C0:c1 - C0, :5] = mm[r0:r1, c0:c1]
            hist = self.hist_window(n0, n1, s0, s1)
            kind, body = 2, full.ravel()
        hdr = np.array([kind, T, b, flags], np.uint32).tobytes()
        out = hdr + hist.astype(np.float64).tobytes() + body.tobytes()
        self.tiles.put(key, out)
        return out

    def hist_window(self, n0, n1, s0, s1):
        HB, m = self.HB, self.meta
        r0, r1 = max(n0 // HB, 0), min(-(-n1 // HB), m["HR"])
        c0 = max(s0 // HB - m["HCMIN"], 0)
        c1 = min(-(-s1 // HB) - m["HCMIN"], m["HC"])
        if r1 <= r0 or c1 <= c0:
            return np.zeros(HBINS)
        return self.hist[r0:r1, c0:c1].reshape(-1, HBINS).sum(axis=0, dtype=np.float64)

    def chains_in(self, n0, n1, s0, s1, lmin, lmax, limit):
        ch = self.chains
        a, z = np.searchsorted(ch[:, 0], n0), np.searchsorted(ch[:, 0], n1)
        sub = ch[a:z]
        sel = sub[(sub[:, 1] >= s0) & (sub[:, 1] < s1) & (sub[:, 2] >= lmin) & (sub[:, 2] <= lmax)]
        total = len(sel)
        if total > limit:  # keep the longest
            sel = sel[np.argsort(-sel[:, 2], kind="stable")[:limit]]
        return total, sel

    def diagonals(self):
        """Persistent left-wedge structures: chains of length >= 18 with s < -n/4 that recur at the same
        distance dl = s + n from the left edge every 32 rows."""
        if self._diag is None:
            ch = self.chains[self.chains[:, 2] >= 18].astype(np.int64)
            w = ch[4 * ch[:, 1] < -ch[:, 0]]
            dl = w[:, 1] + w[:, 0]
            out = []
            order = np.lexsort((w[:, 0], dl))
            w, dl = w[order], dl[order]
            cuts = np.flatnonzero(np.diff(dl)) + 1
            for g in np.split(np.arange(len(w)), cuts):
                if len(g) < 3:
                    continue
                rows = w[g, 0]
                gaps = np.diff(rows)
                if np.all(gaps % 32 == 0):
                    out.append({"dl": int(dl[g[0]]), "n0": int(rows[0]), "n1": int(rows[-1]), "count": int(len(g)),
                                "phase": int(rows[0] % 32), "lengths": sorted(set(int(x) for x in w[g, 2])),
                                "min_gap": int(gaps.min()) if len(gaps) else 0})
            self._diag = out
        return self._diag


class Centre:
    def __init__(self, path):
        self.path = path
        self.summary = json.load(open(os.path.join(path, "summary.json")))
        self.nbits = self.summary["nbits"]
        self.bytes = np.fromfile(os.path.join(path, "centre.bits"), dtype=np.uint8)
        self.k0, self.k1 = self.summary["k0"], self.summary["k1"]
        self.lv = {k: np.fromfile(os.path.join(path, f"lvl_{k}.i32"), dtype=np.int32).reshape(-1, 3)
                   for k in range(self.k0, self.k1 + 1)}
        self.RB = self.summary["run_block"]
        self.runs = np.fromfile(os.path.join(path, "runs16.u32"), dtype=np.uint32).reshape(-1, 2, 64)
        self.runs_cum = np.concatenate([np.zeros((1, 2, 64), np.int64), np.cumsum(self.runs, axis=0, dtype=np.int64)])
        self.ties = np.fromfile(os.path.join(path, "ties.i64"), dtype=np.int64)

    def W_at(self, p):
        """W(p) = 2 * ones(first p bits) - p."""
        a = (p >> self.k0) << self.k0
        w = int(self.lv[self.k0][(a >> self.k0) - 1, 2]) if a else 0
        if p > a:
            bits = self.bits(a, p)
            w += 2 * int(bits.sum()) - (p - a)
        return w

    def bits(self, lo, hi):
        b = np.unpackbits(self.bytes[lo >> 3: (hi + 7) >> 3])
        return b[lo - ((lo >> 3) << 3): hi - ((lo >> 3) << 3)]

    def series(self, lo, hi, px):
        """Per-pixel min, max and final excess E = W/2 over p in (lo, hi]."""
        lo, hi = max(0, lo), min(self.nbits, hi)
        span = hi - lo
        px = max(1, min(px, 4000, span))
        per = span / px
        if per >= (1 << self.k0):
            k = min(self.k1, int(math.floor(math.log2(per))))
            B = 1 << k
            lv = self.lv[k]
            j0, j1 = lo // B, min(-(-hi // B), len(lv))
            n = min(px, j1 - j0)
            idx = (np.arange(n) * (j1 - j0)) // n
            seg = lv[j0:j1]
            mn = np.minimum.reduceat(seg[:, 0], idx)
            mx = np.maximum.reduceat(seg[:, 1], idx)
            ends = np.append(idx[1:], j1 - j0) - 1
            en = seg[ends, 2]
            p0 = (j0 + idx) * B
            p1 = np.minimum((j0 + ends + 1) * B, self.nbits)
            mode = f"blocks of 2^{k}"
        else:
            w0 = self.W_at(lo)
            b = self.bits(lo, hi).astype(np.int32)
            W = w0 + np.cumsum(2 * b - 1)
            n = px
            idx = (np.arange(n) * span) // n
            mn = np.minimum.reduceat(W, idx)
            mx = np.maximum.reduceat(W, idx)
            ends = np.append(idx[1:], span) - 1
            en = W[ends]
            p0, p1 = lo + idx, lo + ends + 1
            mode = "bits"
        return {"lo": lo, "hi": hi, "mode": mode, "p0": p0.tolist(), "p1": p1.tolist(),
                "emin": (mn / 2).tolist(), "emax": (mx / 2).tolist(), "eend": (en / 2).tolist()}

    def run_hist(self, lo, hi):
        lo, hi = max(0, lo), min(self.nbits, hi)
        RB = 1 << self.RB
        a, z = -(-lo // RB), hi // RB
        ones, zeros = np.zeros(64, np.int64), np.zeros(64, np.int64)
        o, zz = np.zeros(64, np.int64), np.zeros(64, np.int64)
        if a < z:
            full = self.runs_cum[z] - self.runs_cum[a]
            ones += full[0]
            zeros += full[1]
            parts = [(lo, a * RB), (z * RB, hi)]
        else:
            parts = [(lo, hi)]
        for x, y in parts:
            if y > x:
                lib.r30c_runs(self.bytes, self.nbits, x, y, o, zz)
                ones += o
                zeros += zz
        return {"lo": lo, "hi": hi, "black": ones.tolist(), "white": zeros.tolist()}

    def find(self, lo, hi, pat, px, maxpos=2000):
        lo, hi = max(0, lo), min(self.nbits, hi)
        pos, bins = np.zeros(maxpos, np.int64), np.zeros(max(1, px), np.int64)
        cnt = lib.r30c_find(self.bytes, self.nbits, lo, hi, pat.encode(), pos, maxpos, bins, len(bins))
        k = len(pat)
        expect = max(0, hi - lo) / 2 ** k
        return {"lo": lo, "hi": hi, "pattern": pat, "count": int(cnt), "expected": expect,
                "z": (cnt - expect) / math.sqrt(expect) if expect > 0 else 0,
                "positions": pos[: min(cnt, maxpos)].tolist(), "bins": bins.tolist()}


RUNS = {}
CENTRE = None


def load():
    global CENTRE
    tdir = os.path.join(DATA, "tiles")
    for name in sorted(os.listdir(tdir)) if os.path.isdir(tdir) else []:
        p = os.path.join(tdir, name)
        if os.path.exists(os.path.join(p, "meta.json")) and os.path.isdir(os.path.join(p, "ck")):
            try:
                RUNS[name] = Run(name, p)
            except Exception as e:  # a run still being generated has no levels yet
                print(f"skip run {name}: {e}", file=sys.stderr)
    c = os.path.join(DATA, "centre")
    if os.path.exists(os.path.join(c, "summary.json")):
        CENTRE = Centre(c)


def clean_pat(s):
    s = (s or "").strip()
    if s and (len(s) > 63 or any(ch not in "01" for ch in s)):
        raise ValueError("pattern must be 1-63 characters of 0/1")
    return s or None


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *a, **kw):
        super().__init__(*a, directory=os.path.join(HERE, "static"), **kw)

    def log_message(self, fmt, *args):
        pass

    def end_headers(self):
        self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def send(self, body, ctype="application/json", code=200):
        if isinstance(body, (dict, list)):
            body = json.dumps(body).encode()
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        u = urlparse(self.path)
        if not u.path.startswith("/api/"):
            return super().do_GET()
        q = {k: v[0] for k, v in parse_qs(u.query).items()}
        I = lambda k, d=0: int(float(q.get(k, d)))
        try:
            r = RUNS.get(q.get("run", ""))
            if u.path == "/api/meta":
                return self.send({"runs": {k: v.info() for k, v in RUNS.items()},
                                  "centre": CENTRE.summary if CENTRE else None, "tile": T})
            if u.path == "/api/tile":
                return self.send(r.tile(I("b", 1), I("ty"), I("tx"), clean_pat(q.get("pat"))),
                                 "application/octet-stream")
            if u.path == "/api/chains":
                total, sel = r.chains_in(I("n0"), I("n1"), I("s0"), I("s1"), I("lmin", 0), I("lmax", 1 << 30),
                                         I("limit", 20000))
                hdr = np.array([total, len(sel)], np.int32).tobytes()
                return self.send(hdr + sel.astype(np.int32).tobytes(), "application/octet-stream")
            if u.path == "/api/hist":
                return self.send({"hist": r.hist_window(I("n0"), I("n1"), I("s0"), I("s1")).tolist(), "HB": r.HB})
            if u.path == "/api/diagonals":
                return self.send(r.diagonals())
            if u.path == "/api/rows":  # per-row records n, centre bit, ones, new chains, active pyramids
                n0, n1 = max(0, I("n0")), min(r.N, I("n1"))
                a = np.memmap(os.path.join(r.path, "rows.i64"), dtype=np.int64, mode="r").reshape(-1, 5)
                return self.send(np.ascontiguousarray(a[n0:n1]).tobytes(), "application/octet-stream")
            if u.path == "/api/centre/series":
                return self.send(CENTRE.series(I("lo"), I("hi"), I("px", 1000)))
            if u.path == "/api/centre/bits":
                lo, hi = I("lo"), I("hi")
                hi = min(hi, lo + (1 << 20), CENTRE.nbits)
                return self.send(CENTRE.bits(max(lo, 0), hi).astype(np.uint8).tobytes(), "application/octet-stream")
            if u.path == "/api/centre/runs":
                return self.send(CENTRE.run_hist(I("lo"), I("hi")))
            if u.path == "/api/centre/find":
                pat = clean_pat(q.get("pat"))
                if not pat:
                    raise ValueError("empty pattern")
                return self.send(CENTRE.find(I("lo"), I("hi"), pat, I("px", 1000)))
            if u.path == "/api/centre/ties":
                return self.send(CENTRE.ties.tolist())
            return self.send({"error": "unknown endpoint"}, code=404)
        except Exception as e:
            return self.send({"error": f"{type(e).__name__}: {e}"}, code=400)


class Server(ThreadingHTTPServer):
    daemon_threads = True

    def handle_error(self, request, client_address):
        if isinstance(sys.exc_info()[1], (BrokenPipeError, ConnectionResetError)):
            return  # the browser dropped a tile it no longer needs
        super().handle_error(request, client_address)


def main():
    port = int(sys.argv[1]) if len(sys.argv) > 1 else int(os.environ.get("PORT", 8730))
    load()
    print(f"runs: {', '.join(RUNS) or 'none'}; centre: {'yes' if CENTRE else 'no'}", flush=True)
    srv = Server(("127.0.0.1", port), Handler)
    print(f"Rule 30 viewer on http://127.0.0.1:{port}/", flush=True)
    srv.serve_forever()


if __name__ == "__main__":
    main()
