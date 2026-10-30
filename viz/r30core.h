// Shared bit-parallel Rule 30 helpers for the viewer. Definitions match ../rule30.c and ../xfiles.c:
// rows start at n = 0 with one black cell, columns s = -n..n live at absolute bit MID + s, bit i's left
// neighbour is bit i-1. A chain is a maximal run of >= 2 adjacent 111 neighbourhoods; its seed cells are
// the chain's cells plus the 110 cell to its right, and membership persists through 000 neighbourhoods.
#ifndef R30CORE_H
#define R30CORE_H
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef uint64_t u64;

static inline u64 lsh(const u64 *x, int64_t w) { return (x[w] << 1) | (x[w - 1] >> 63); }  // left neighbour's value
static inline u64 rsh(const u64 *x, int64_t w) { return (x[w] >> 1) | (x[w + 1] << 63); }  // right neighbour's value
static inline int bit(const u64 *x, int64_t i) { return (x[i >> 6] >> (i & 63)) & 1; }
// value of cell i+j at bit i, 0 <= j < 64
static inline u64 ahead(const u64 *x, int64_t w, int j) { return j ? (x[w] >> j) | (x[w + 1] << (64 - j)) : x[w]; }
// value of cell i-j at bit i, 0 <= j < 64
static inline u64 behind(const u64 *x, int64_t w, int j) { return j ? (x[w] << j) | (x[w - 1] >> (64 - j)) : x[w]; }

// cells lo..hi (absolute, inclusive) within word w
static inline u64 span(int64_t w, int64_t lo, int64_t hi) {
    int64_t a = w * 64, b = a + 63;
    if (hi < a || lo > b) return 0;
    if (lo < a) lo = a;
    if (hi > b) hi = b;
    int64_t len = hi - lo + 1;
    return (len == 64 ? ~0ULL : ((1ULL << len) - 1)) << (lo - a);
}

static inline int64_t floordiv(int64_t a, int64_t b) { return a >= 0 ? a / b : -((-a + b - 1) / b); }

static inline u64 *alloc64(int64_t n) {
    u64 *p = calloc(n, sizeof(u64));
    if (!p) { perror("calloc"); exit(1); }
    return p;
}

// One generation over words [wlo, whi]: next row, 111/000 masks, chain mask in2, seeds and next membership.
static inline void gen_masks(const u64 *row, u64 *nxt, u64 *is111, u64 *is000, int64_t wlo, int64_t whi) {
    for (int64_t w = wlo; w <= whi; w++) {
        u64 L = lsh(row, w), C = row[w], R = rsh(row, w);
        is111[w] = L & C & R;
        is000[w] = ~L & ~C & ~R;
        nxt[w] = L ^ (C | R);
    }
}
static inline void gen_in2(const u64 *is111, u64 *in2, int64_t wlo, int64_t whi) {
    for (int64_t w = wlo; w <= whi; w++) in2[w] = is111[w] & (lsh(is111, w) | rsh(is111, w));
}
#endif
