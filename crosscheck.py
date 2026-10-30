#!/usr/bin/env python3
"""Independent check of xfiles.c chain logging: plain numpy byte arithmetic, no bit packing.

Evolves Rule 30 to row hi, and for rows lo..hi-1 lists every chain of >= minlen adjacent 111
neighbourhoods (a black run of length B holds a chain of B - 2), then diffs against chains.i64.
"""
import sys

import numpy as np

lo, hi, chains_path = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
minlen = 18
W = 2 * hi + 3
mid = hi + 1
row = np.zeros(W, dtype=np.uint8)
row[mid] = 1
mine = []
for n in range(hi):
    if n >= lo:
        b = np.concatenate(([0], row, [0])).astype(np.int8)
        d = np.diff(b)
        starts, ends = np.flatnonzero(d == 1), np.flatnonzero(d == -1)  # black runs [starts, ends)
        for a, e in zip(starts, ends):
            B = e - a
            if B - 2 >= minlen:
                mine.append((n, a + 1 - mid, B - 2))  # first 111 centre is the run's second cell
    nxt = np.zeros_like(row)
    nxt[1:-1] = row[:-2] ^ (row[1:-1] | row[2:])
    row = nxt
    if n % 50000 == 0:
        print("row", n, flush=True)

ch = np.fromfile(chains_path, dtype=np.int64).reshape(-1, 3)
ch = ch[(ch[:, 0] >= lo) & (ch[:, 0] < hi)]
theirs = set(map(tuple, ch.tolist()))
mine = set(mine)
print(f"rows {lo}..{hi - 1}: python {len(mine)} chains, xfiles.c {len(theirs)}, "
      f"only python {len(mine - theirs)}, only C {len(theirs - mine)}")
for L in (19, 20, 21):
    print(L, sum(1 for x in mine if x[2] == L), sum(1 for x in theirs if x[2] == L))
