#!/usr/bin/env python3
"""Checks for the viewer's engines and data. Prints one PASS/FAIL line per check.

    python validate.py [run]          (run defaults to 1m)

1. r30gen per-row records vs ../out_1m-era rows.i64 (rule30.c) and chains >= 18 vs trace_1m/chains.i64 (xfiles.c).
2. libr30 r30_window vs a plain numpy byte-array evolution (cells, chain cells, pyramid labels, pattern
   matches, block counts) on random windows, including windows that straddle the triangle's edges.
3. On-demand blocks vs the precomputed level grid, and level totals vs the per-row records.
4. Centre column: the Data Repository bits vs rule30.c's centre bits for the first 10^6 rows, plus the
   known totals.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import server  # noqa: E402

STATS = "/Volumes/BLIPS/rule30/2026-10-22_stats_1m/rows.i64"
TRACE = "/Volumes/BLIPS/rule30/2026-10-23_xfiles/trace_1m/chains.i64"
fails = 0


def check(name, ok, detail=""):
    global fails
    fails += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}", flush=True)


def brute(N, pat):
    """Byte-per-cell reference: rows, chain cells, pyramid membership with labels, pattern starts/cover."""
    W = 2 * N + 8
    mid = N + 4
    row = np.zeros(W, np.uint8)
    row[mid] = 1
    lab = np.zeros(W, np.uint8)
    Pm = np.zeros(W, bool)
    cells, chainc, labs, starts, cover, chains = [], [], [], [], [], []
    k = len(pat) if pat else 0
    pv = np.array([c == "1" for c in pat], np.uint8) if pat else None
    for n in range(N):
        L, C, R = np.roll(row, 1), row, np.roll(row, -1)
        is111 = (L & C & R).astype(bool)
        is000 = ~(L | C | R).astype(bool)
        in2 = is111 & (np.roll(is111, 1) | np.roll(is111, -1))
        cc = np.zeros(W, np.uint8)
        d = np.diff(np.concatenate(([0], in2.astype(np.int8), [0])))
        st, en = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
        newlab = np.where(is000, lab, 0).astype(np.uint8)
        for a, e in zip(st, en):
            ln = e - a
            chains.append((n, a - mid, ln))
            cc[a:e] = 1
            newlab[a:e + 1] = min(ln, 63)
        cells.append(row.copy())
        chainc.append(cc)
        labs.append(np.where(Pm, lab, 0).astype(np.uint8))
        if k:
            s = np.ones(W, bool)
            for j in range(k):
                s &= np.roll(row, -j) == pv[j]
            idx = np.arange(W) - mid
            s &= (idx >= -n) & (idx <= n - k + 1)
            cv = np.zeros(W, bool)
            for j in range(k):
                cv |= np.roll(s, j)
            starts.append(s)
            cover.append(cv)
        seed = in2 | np.roll(in2, 1)
        Pm = seed | (is000 & Pm)
        lab = newlab
        row = (L ^ (C | R)).astype(np.uint8)
    out = {"mid": mid, "cells": np.array(cells), "chain": np.array(chainc), "lab": np.array(labs), "chains": chains}
    if k:
        out["start"], out["cover"] = np.array(starts), np.array(cover)
    return out


def main():
    run_name = sys.argv[1] if len(sys.argv) > 1 else "1m"
    server.load()
    run = server.RUNS[run_name]
    N = run.N
    rows = np.fromfile(os.path.join(run.path, "rows.i64"), dtype=np.int64).reshape(-1, 5)
    ref = np.fromfile(STATS, dtype=np.int64).reshape(-1, 20)[:N]
    m = min(len(rows), len(ref))
    for name, a, b in [("row", 0, 0), ("centre bit", 1, 1), ("ones", 2, 9), ("new chains", 3, 6), ("active", 4, 7)]:
        bad = np.flatnonzero(rows[:m, a] != ref[:m, b])
        check(f"r30gen vs rule30.c rows.i64: {name}", len(bad) == 0, f"{m} rows, {len(bad)} differ")
    tr = np.fromfile(TRACE, dtype=np.int64).reshape(-1, 3)
    tr = tr[tr[:, 0] < N]
    mine = run.chains[run.chains[:, 2] >= 18].astype(np.int64)
    check("r30gen vs xfiles.c chains.i64 (length >= 18)", np.array_equal(mine, tr), f"{len(mine)} vs {len(tr)} chains")
    hist = np.fromfile(os.path.join(os.path.dirname(STATS), "hist.i64"), dtype=np.int64)
    if N == 1000000:
        g = np.array(run.meta["hist"], np.int64)
        check("r30gen chain-length histogram vs rule30.c hist.i64", np.array_equal(g[:63], hist[:63]),
              f"total {g.sum()} vs {hist.sum()}")

    # 2. window renderer vs brute force
    BN, pat = 1500, "11011"
    br = brute(BN, pat)
    mid = br["mid"]
    rng = np.random.default_rng(30)
    probe = run if run.N >= BN else None
    if probe is None:
        check("window renderer", False, "run too small")
    else:
        wins = [(0, 256, -128, 128), (300, 556, -700, -444), (1000, 1256, 1100, 1356), (1244, 1500, -1500, -1244)]
        for _ in range(8):
            n0 = int(rng.integers(0, BN - 256))
            s0 = int(rng.integers(-n0 - 300, n0 + 50))
            wins.append((n0, n0 + 256, s0, s0 + 256))
        bad = 0
        for n0, n1, s0, s1 in wins:
            cell, aux, _, h = probe.window(n0, n1, s0, s1, 1, pat)
            cell, aux = cell.reshape(n1 - n0, s1 - s0), aux.reshape(n1 - n0, s1 - s0)
            sl = (slice(n0, n1), slice(mid + s0, mid + s1))
            exp = br["cells"][sl] | (br["chain"][sl] << 1) | (br["lab"][sl] << 2)
            expa = br["start"][sl].astype(np.uint8) | (br["cover"][sl].astype(np.uint8) << 1)
            eh = np.zeros(64)
            for n, s, ln in br["chains"]:
                if n0 <= n < n1 and s0 <= s < s1:
                    eh[min(ln, 63)] += 1
            ok = np.array_equal(cell, exp) and np.array_equal(aux, expa) and np.array_equal(h, eh)
            bad += not ok
            if not ok:
                print("   mismatch", (n0, n1, s0, s1), int((cell != exp).sum()), int((aux != expa).sum()),
                      np.flatnonzero(h != eh))
        check("r30_window b=1 vs numpy (cells, chain cells, pyramid labels, '11011' matches, histogram)", bad == 0,
              f"{len(wins)} windows of 256x256")
        bad = 0
        for b in (2, 8, 64):
            for _ in range(3):
                n0 = int(rng.integers(0, (BN - 4 * b) // b)) * b
                n1 = min(n0 + 16 * b, BN)
                n1 = n0 + ((n1 - n0) // b) * b
                s0 = int(rng.integers(-n0 // b - 2, n0 // b + 2)) * b - 4 * b
                s1 = s0 + 16 * b
                _, _, blk, h = probe.window(n0, n1, s0, s1, b, pat)
                blk = blk.reshape((n1 - n0) // b, 16, 6)
                sl = (slice(n0, n1), slice(mid + s0, mid + s1))
                idx = np.arange(s0, s1)[None, :]
                nn = np.arange(n0, n1)[:, None]
                inside = (idx >= -nn) & (idx <= nn)
                def red(x):
                    return x.reshape((n1 - n0) // b, b, 16, b).sum(axis=(1, 3))
                e = [red(inside), red(br["cells"][sl] & inside), red((br["lab"][sl] > 0) & inside),
                     None, None, red(br["start"][sl] & inside)]
                heads = np.zeros(((n1 - n0) // b, 16), np.int64)
                mx = np.zeros_like(heads)
                for n, s, ln in br["chains"]:
                    if n0 <= n < n1 and s0 <= s < s1:
                        heads[(n - n0) // b, (s - s0) // b] += 1
                        mx[(n - n0) // b, (s - s0) // b] = max(mx[(n - n0) // b, (s - s0) // b], ln)
                e[3], e[4] = heads, mx
                ok = all(np.array_equal(blk[:, :, i], e[i]) for i in range(6))
                bad += not ok
                if not ok:
                    print("   block mismatch", b, (n0, n1, s0, s1), [int((blk[:, :, i] != e[i]).sum()) for i in range(6)])
        check("r30_window b=2,8,64 blocks vs numpy (cells, ones, pyramid cells, heads, max length, matches)", bad == 0)

    # 3. on-demand vs precomputed levels; level totals vs per-row records
    L, mm = run.levels[run.base]
    B = run.base
    bad = 0
    for ty, tx in [(0, -1), (0, 0), (1, -1), (3, -4), (3, 2)]:
        if (ty + 1) * 256 * (B // 2) > N:
            continue
        t = np.frombuffer(run.tile(B // 2, ty, tx, None), dtype=np.uint8)
        blk = t[528:].view(np.uint32).reshape(256, 256, 6)
        agg = blk.reshape(128, 2, 128, 2, 6).sum(axis=(1, 3))
        R0, C0 = ty * 128, tx * 128 - L["CMIN"]
        if C0 < 0 or C0 + 128 > L["NC"]:
            continue
        lv = np.asarray(mm[R0:R0 + 128, C0:C0 + 128])
        bad += not np.array_equal(agg[:, :, :4], lv[:, :, :4])
    check(f"on-demand b={B // 2} tiles (2x2 summed) vs precomputed level_{B}", bad == 0)
    nr = N // B
    lvl = np.asarray(mm[:nr]).sum(axis=1, dtype=np.int64)
    rs = rows[: nr * B].reshape(nr, B, 5).sum(axis=1)
    check(f"level_{B} ones per block row vs rows.i64", np.array_equal(lvl[:, 1], rs[:, 2]))
    check(f"level_{B} chain heads per block row vs rows.i64 new chains", np.array_equal(lvl[:, 3], rs[:, 3]))
    for Bc in sorted(run.levels)[1:]:
        tot = np.asarray(run.levels[Bc][1]).reshape(-1, 5).sum(axis=0, dtype=np.int64)
        if not np.array_equal(tot[:4], np.asarray(mm).reshape(-1, 5).sum(axis=0, dtype=np.int64)[:4]):
            check(f"level_{Bc} totals", False)
            break
    else:
        check("all coarser levels preserve totals", True, f"{len(run.levels)} levels")

    # 4. centre column
    c = server.Centre(os.path.join(server.DATA, "centre"))
    ours = np.fromfile(STATS, dtype=np.int64).reshape(-1, 20)[:, 1]
    theirs = c.bits(0, len(ours))
    check("Data Repository centre bits vs rule30.c centre column", np.array_equal(ours, theirs), f"{len(ours)} bits")
    s = c.summary
    check("centre totals", s["ones"] == 500025038 and s["last_tie"] == 195112,
          f"ones {s['ones']}, last tie {s['last_tie']}, ties {s['ties']}")
    walk = np.loadtxt("/Volumes/BLIPS/rule30/2026-10-23_xfiles/xfile3_billion/walk_every_1e6.csv", delimiter=",")
    mine = np.array([c.lv[c.k0][(p >> c.k0) - 1, 2] / 2 if p % (1 << c.k0) == 0 else c.W_at(p) / 2
                     for p in range(10**6, 10**9 + 1, 10**6)])
    check("excess at every 10^6 vs walk_every_1e6.csv", np.array_equal(mine, walk), f"{len(walk)} checkpoints")
    print("ALL PASS" if not fails else f"{fails} FAILED")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
