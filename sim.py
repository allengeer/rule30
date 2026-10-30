#!/usr/bin/env python3
"""Rule 30 per-row statistics, run until the result set reaches a byte budget.

Rule 30: new = left XOR (centre OR right), from a single black cell.
Row n spans columns -n..n (offset from the centre column).

Per row n this records:
  * row period     -- smallest p with row[i] == row[i+p] across the row's light cone
  * centre period  -- smallest period of the centre-column prefix c[0..n] (KMP, incremental);
                      equals n+1 while the column has shown no repetition
  * E(cell)        -- per column, fraction of black cells seen in that column so far
                      (only rows where the column is inside the light cone count)
  * pyramids       -- per column, how many pyramids have covered it so far

A pyramid is a maximal chain of >= 2 adjacent 111 neighbourhoods (Rule 30 case 8,
black-black-black) in one row. That chain whitens the cells under it plus the 110 cell
to its right, and the white triangle then shrinks down through 000 neighbourhoods.
The first one sits on row 3, cols +1..+2, and its apex lands on row 5 (row 6 counting
from 1), col +2.

Outputs (in --out):
  summary.csv  one line of scalars per row
  e.f32        per row, float32 E for columns -n..n   (row n starts at offset 4*n*n)
  pyr.u16      per row, uint16 pyramid counts for columns -n..n (offset 2*n*n)
"""
import argparse
import math
import os
import sys
import time

import numpy as np

RULE30 = np.array([(30 >> k) & 1 for k in range(8)], dtype=np.uint8)
BYTES_PER_CELL = 4 + 2  # float32 E + uint16 pyramid count


def rows_for_budget(budget, summary_bytes_per_row=96):
    # sum_{n<N} (2n+1) * BYTES_PER_CELL = BYTES_PER_CELL * N^2, plus one summary line per row
    a, b = BYTES_PER_CELL, summary_bytes_per_row
    return int((-b + math.sqrt(b * b + 4 * a * budget)) / (2 * a))


def smallest_period(bits):
    """Smallest p such that bits is p-periodic (KMP failure function)."""
    s = bits.tobytes()
    n = len(s)
    fail = [0] * n
    k = 0
    for i in range(1, n):
        while k and s[i] != s[k]:
            k = fail[k - 1]
        if s[i] == s[k]:
            k += 1
        fail[i] = k
    return n - fail[-1]


def run(out, budget):
    os.makedirs(out, exist_ok=True)
    N = rows_for_budget(budget)
    W = 2 * N + 3  # one dead cell of padding each side
    mid = N + 1

    row = np.zeros(W, dtype=np.uint8)
    row[mid] = 1
    ones = np.zeros(W, dtype=np.uint64)
    pyr_count = np.zeros(W, dtype=np.uint16)
    pyr_id = np.zeros(W, dtype=np.int64)  # pyramid owning each white cell, 0 = none
    next_id = 1

    # incremental KMP over the centre column
    centre = bytearray()
    cfail = []
    ck = 0

    f_e = open(os.path.join(out, "e.f32"), "wb")
    f_p = open(os.path.join(out, "pyr.u16"), "wb")
    f_s = open(os.path.join(out, "summary.csv"), "w")
    f_s.write("row,row_period,centre_period,centre_bit,centre_E,centre_pyramids,"
              "new_pyramids,active_pyramids,total_pyramids,mean_E_left,mean_E_right,"
              "max_abs_dev_centre_band\n")

    t0 = time.time()
    for n in range(N):
        lo, hi = mid - n, mid + n + 1
        cone = row[lo:hi]
        ones[lo:hi] += cone

        # E per column: column c entered the cone at row |c|
        offs = np.abs(np.arange(-n, n + 1))
        E = (ones[lo:hi] / (n - offs + 1)).astype(np.float32)
        E.tofile(f_e)
        pyr_count[lo:hi].tofile(f_p)

        # centre column prefix period
        b = row[mid]
        centre.append(b)
        i = len(centre) - 1
        if i == 0:
            cfail.append(0)
        else:
            while ck and centre[i] != centre[ck]:
                ck = cfail[ck - 1]
            if centre[i] == centre[ck]:
                ck += 1
            cfail.append(ck)
        centre_period = len(centre) - cfail[-1]

        row_period = smallest_period(cone)

        # next row + pyramid bookkeeping
        nb = (row[:-2] << 2) | (row[1:-1] << 1) | row[2:]  # neighbourhood of cells 1..W-2
        new = np.zeros_like(row)
        new[1:-1] = RULE30[nb]
        new_id = np.zeros_like(pyr_id)
        inherit = (nb == 0) & (pyr_id[1:-1] > 0)
        new_id[1:-1][inherit] = pyr_id[1:-1][inherit]

        is111 = np.concatenate(([0], (nb == 7).astype(np.int8), [0]))
        d = np.diff(is111)
        starts = np.flatnonzero(d == 1)  # index into nb == cell index - 1
        ends = np.flatnonzero(d == -1)   # exclusive
        centre_pyr = int(pyr_count[mid])  # before this row's seeds, matching pyr.u16
        new_pyr = 0
        for s, e in zip(starts, ends):
            if e - s >= 2:
                a, z = s + 1, e + 1  # cells a..z-1 under the chain, plus cell z (the 110)
                new_id[a:z + 1] = next_id
                pyr_count[a:z + 1] += 1
                next_id += 1
                new_pyr += 1
        active = len(np.unique(new_id[new_id > 0]))

        band = min(n, 64)
        dev = float(np.max(np.abs(E[n - band:n + band + 1] - 0.5))) if n else 0.5
        f_s.write(f"{n},{row_period},{centre_period},{b},{E[n]:.6f},{centre_pyr},"
                  f"{new_pyr},{active},{next_id - 1},"
                  f"{float(E[:n].mean()) if n else 0:.6f},{float(E[n + 1:].mean()) if n else 0:.6f},"
                  f"{dev:.6f}\n")

        row, pyr_id = new, new_id
        if n % 250 == 0:
            f_s.flush()
            print(f"row {n}/{N}  {time.time() - t0:.1f}s  centre_E={E[n]:.4f} "
                  f"pyramids={next_id - 1}", flush=True)

    for f in (f_e, f_p, f_s):
        f.close()
    total = sum(os.path.getsize(os.path.join(out, x)) for x in os.listdir(out))
    print(f"done: {N} rows, {total / 1e6:.1f} MB in {out}, {time.time() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "out"))
    ap.add_argument("--budget-mb", type=float, default=100)
    args = ap.parse_args()
    run(args.out, int(args.budget_mb * 1e6))
    sys.exit(0)
