#!/usr/bin/env python3
"""Independent numpy check of `xfiles frontier` and `xfiles damage`.
usage: frontier_check.py dir N T seed"""
import sys

import numpy as np

d, N, T, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])

# Frontier: keep the triangle as rows of width 2N+1, left edge of row m at index N-m.
L = np.fromfile(f"{d}/frontier.i64", np.int64)[:N]
rows = np.zeros((N, 2 * N + 3), np.uint8)
rows[0, N + 1] = 1
for m in range(1, N):
    p = rows[m - 1]
    rows[m, 1:-1] = p[:-2] ^ (p[1:-1] | p[2:])
bad = 0
for m in range(32, N):
    dd = np.arange(0, 2 * m + 1)
    now = rows[m, N + 1 - m + dd]
    old = np.zeros_like(now)                        # diagonal d of row m-32, zero past its right edge
    k = 2 * (m - 32) + 1
    old[:k] = rows[m - 32, N + 1 - (m - 32) + dd[:k]]
    viol = np.flatnonzero(now != old)
    want = viol[0] if len(viol) else -1
    bad += want != L[m]
print(f"frontier: {N - 32} rows compared, {bad} mismatches; L/m at end = {L[N-1]/(N-1):.4f}")

# Damage: the same splitmix64 stream, plain numpy evolution on the ring.
M = (1 << 64) - 1
def splitmix(s):
    s = (s + 0x9e3779b97f4a7c15) & M
    z = s
    z = ((z ^ (z >> 30)) * 0xbf58476d1ce4e5b9) & M
    z = ((z ^ (z >> 27)) * 0x94d049bb133111eb) & M
    return s, z ^ (z >> 31)
nw = (2 * T + 128) // 64
W = nw * 64
words, s = [], seed
for _ in range(nw):
    s, z = splitmix(s)
    words.append(z)
a = np.array([(w >> i) & 1 for w in words for i in range(64)], np.uint8)
b = a.copy()
c = 3 * W // 10
b[c] ^= 1
front = np.fromfile(f"{d}/damage_{seed}.i32", np.int32)
mis = 0
for t in range(T):
    a = np.roll(a, 1) ^ (a | np.roll(a, -1))
    b = np.roll(b, 1) ^ (b | np.roll(b, -1))
    mis += (c - np.flatnonzero(a != b)[0]) != front[t]
print(f"damage: {T} steps compared, {mis} mismatches; front/t at end = {front[T-1]/T:.4f}")
