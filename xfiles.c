// Instrumented Rule 30 runs for the X-File investigations. Pyramid definitions match rule30.c.
//
//   xfiles trace N dir
//       From a single black cell for N rows. Writes:
//         chains.i64   (row, start offset, length) for every chain of length >= LOGMIN
//         hist2.i64    [HIST lengths][BUCKETS] chain counts by length and by distance to the
//                      nearer edge, bucket = floor(log2(distance + 1))
//         edge.i64     [3][EDGE]: ones at distance d from the left edge, from the right
//                      edge, and how many rows reached distance d
//   xfiles random W T warm seed dir
//       Rule 30 on a ring of W cells from iid fair random bits, for T steps. Writes
//         random_rows.f64  (step, density, new/width, active/width)
//         random_hist.i64  chain-length histogram over steps >= warm
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef uint64_t u64;
#define HIST 256
#define BUCKETS 24
#define EDGE 8192
#define LOGMIN 18

static u64 *alloc(int64_t n) { u64 *p = calloc(n, sizeof(u64)); if (!p) { perror("calloc"); exit(1); } return p; }
static inline u64 lsh(const u64 *x, int64_t w) { return (x[w] << 1) | (x[w - 1] >> 63); }
static inline u64 rsh(const u64 *x, int64_t w) { return (x[w] >> 1) | (x[w + 1] << 63); }
static inline int bit(const u64 *x, int64_t i) { return (x[i >> 6] >> (i & 63)) & 1; }
static int lg(int64_t d) { int b = 0; while (d > 1) { d >>= 1; b++; } return b; }

static FILE *open_or_die(const char *dir, const char *name) {
    char path[4096];
    snprintf(path, sizeof path, "%s/%s", dir, name);
    FILE *f = fopen(path, "wb");
    if (!f) { perror(path); exit(1); }
    return f;
}

// One generation over words [wlo, whi]. Fills next row, membership Pn; returns new chains via callback.
typedef void (*chain_cb)(int64_t start_cell, int64_t len, void *ctx);
static void generation(u64 *row, u64 *nxt, u64 *is111, u64 *is000, u64 *in2, u64 *P, u64 *Pn,
                       int64_t wlo, int64_t whi, void (*ghost)(u64 *), int64_t *newp, int64_t *active,
                       chain_cb cb, void *ctx) {
    if (ghost) ghost(row);
    for (int64_t w = wlo; w <= whi; w++) {
        u64 L = lsh(row, w), C = row[w], R = rsh(row, w);
        is111[w] = L & C & R;
        is000[w] = ~L & ~C & ~R;
        nxt[w] = L ^ (C | R);
    }
    if (ghost) ghost(is111);
    for (int64_t w = wlo; w <= whi; w++) in2[w] = is111[w] & (lsh(is111, w) | rsh(is111, w));
    if (ghost) ghost(in2);
    int64_t cur = -1, np = 0, ac = 0;
    for (int64_t w = wlo; w <= whi; w++) {
        u64 l = lsh(in2, w), seed = in2[w] | l;
        Pn[w] = seed | (is000[w] & P[w]);
        u64 starts = in2[w] & ~l, ends = in2[w] & ~rsh(in2, w), ev = starts | ends;
        np += __builtin_popcountll(starts);
        while (ev) {
            int j = __builtin_ctzll(ev);
            ev &= ev - 1;
            int64_t pos = w * 64 + j;
            if ((starts >> j) & 1) cur = pos;
            if (((ends >> j) & 1) && cur >= 0) cb(cur, pos - cur + 1, ctx);
        }
    }
    if (ghost) ghost(Pn);
    for (int64_t w = wlo; w <= whi; w++) ac += __builtin_popcountll(Pn[w] & ~lsh(Pn, w));
    *newp = np;
    *active = ac;
}

// ---- trace -------------------------------------------------------------------------------
typedef struct { int64_t n, mid; int64_t *hist2; FILE *log; } trace_ctx;

static void trace_chain(int64_t start, int64_t len, void *vp) {
    trace_ctx *t = vp;
    int64_t s = start - t->mid, e = s + len - 1;     // offsets of the chain's 111 centres
    int64_t dl = s + t->n, dr = t->n - e;            // distance to left / right edge
    int64_t d = dl < dr ? dl : dr;
    int64_t L = len < HIST ? len : HIST - 1;
    t->hist2[L * BUCKETS + lg(d + 1)]++;
    if (len >= LOGMIN) {
        int64_t rec[3] = {t->n, s, len};
        fwrite(rec, sizeof rec, 1, t->log);
    }
}

static int run_trace(int64_t N, const char *dir) {
    int64_t MID = 128 + N, NW = (2 * N + 256) / 64 + 4;
    u64 *row = alloc(NW), *nxt = alloc(NW), *is111 = alloc(NW), *is000 = alloc(NW), *in2 = alloc(NW),
        *P = alloc(NW), *Pn = alloc(NW);
    int64_t *hist2 = calloc(HIST * BUCKETS, sizeof(int64_t));
    int64_t *edge = calloc(3 * EDGE, sizeof(int64_t));
    row[MID >> 6] |= 1ULL << (MID & 63);
    trace_ctx t = {0, MID, hist2, open_or_die(dir, "chains.i64")};
    for (int64_t n = 0; n < N; n++) {
        int64_t lim = 2 * n + 1 < EDGE ? 2 * n + 1 : EDGE;
        for (int64_t d = 0; d < lim; d++) {
            edge[d] += bit(row, MID - n + d);
            edge[EDGE + d] += bit(row, MID + n - d);
            edge[2 * EDGE + d]++;
        }
        t.n = n;
        int64_t np, ac;
        generation(row, nxt, is111, is000, in2, P, Pn, (MID - n - 1) >> 6, (MID + n + 1) >> 6, NULL,
                   &np, &ac, trace_chain, &t);
        u64 *x = row; row = nxt; nxt = x;
        x = P; P = Pn; Pn = x;
        if (n % 100000 == 0) { fflush(t.log); fprintf(stderr, "trace row %lld/%lld\n", (long long)n, (long long)N); }
    }
    fclose(t.log);
    FILE *f = open_or_die(dir, "hist2.i64"); fwrite(hist2, sizeof(int64_t), HIST * BUCKETS, f); fclose(f);
    f = open_or_die(dir, "edge.i64"); fwrite(edge, sizeof(int64_t), 3 * EDGE, f); fclose(f);
    fprintf(stderr, "trace done\n");
    return 0;
}

// ---- random ring -------------------------------------------------------------------------
static int64_t RW;  // real words, at indices 1..RW
static void ring_ghost(u64 *x) { x[0] = x[RW]; x[RW + 1] = x[1]; }
typedef struct { int64_t *hist; int on; } rand_ctx;
static void rand_chain(int64_t start, int64_t len, void *vp) {
    rand_ctx *r = vp;
    (void)start;
    if (r->on) r->hist[len < HIST ? len : HIST - 1]++;
}

static int run_random(int64_t W, int64_t T, int64_t warm, u64 seed, const char *dir) {
    RW = W / 64;
    int64_t NW = RW + 2;
    u64 *row = alloc(NW), *nxt = alloc(NW), *is111 = alloc(NW), *is000 = alloc(NW), *in2 = alloc(NW),
        *P = alloc(NW), *Pn = alloc(NW);
    u64 s = seed ? seed : 0x9E3779B97F4A7C15ULL;
    for (int64_t w = 1; w <= RW; w++) {  // splitmix64
        u64 z = (s += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        row[w] = z ^ (z >> 31);
    }
    int64_t hist[HIST] = {0};
    rand_ctx r = {hist, 0};
    FILE *f = open_or_die(dir, "random_rows.f64");
    for (int64_t t = 0; t < T; t++) {
        int64_t ones = 0, np, ac;
        for (int64_t w = 1; w <= RW; w++) ones += __builtin_popcountll(row[w]);
        r.on = t >= warm;
        generation(row, nxt, is111, is000, in2, P, Pn, 1, RW, ring_ghost, &np, &ac, rand_chain, &r);
        double rec[4] = {(double)t, (double)ones / W, (double)np / W, (double)ac / W};
        fwrite(rec, sizeof rec, 1, f);
        u64 *x = row; row = nxt; nxt = x;
        x = P; P = Pn; Pn = x;
    }
    fclose(f);
    f = open_or_die(dir, "random_hist.i64"); fwrite(hist, sizeof hist, 1, f); fclose(f);
    fprintf(stderr, "random done\n");
    return 0;
}

// ---- edge reach: ones by distance from each edge, sampling every stride-th row ----------
static int run_edge(int64_t N, int64_t D, int64_t stride, const char *dir) {
    int64_t MID = 128 + N, NW = (2 * N + 256) / 64 + 4;
    u64 *row = alloc(NW), *nxt = alloc(NW);
    int64_t *edge = calloc(3 * D, sizeof(int64_t));
    row[MID >> 6] |= 1ULL << (MID & 63);
    for (int64_t n = 0; n < N; n++) {
        if (n % stride == 0) {
            int64_t lim = 2 * n + 1 < D ? 2 * n + 1 : D;
            for (int64_t d = 0; d < lim; d++) {
                edge[d] += bit(row, MID - n + d);
                edge[D + d] += bit(row, MID + n - d);
                edge[2 * D + d]++;
            }
        }
        for (int64_t w = (MID - n - 1) >> 6; w <= (MID + n + 1) >> 6; w++)
            nxt[w] = lsh(row, w) ^ (row[w] | rsh(row, w));
        u64 *x = row; row = nxt; nxt = x;
    }
    FILE *f = open_or_die(dir, "edge_wide.i64"); fwrite(edge, sizeof(int64_t), 3 * D, f); fclose(f);
    fprintf(stderr, "edge done\n");
    return 0;
}

// ---- wedge: which left diagonals are periodic --------------------------------------------
// At every checkpoint row n (multiple of stride) look back H rows along each diagonal at
// distance d from the left edge (cells (m, -m + d), m = n-H+1..n) and record its smallest
// period p <= H/2 (0 = none). Writes wedge.i64: per checkpoint, n then a period per d < 2n+1-?.
static int run_wedge(int64_t N, int64_t stride, int64_t H, const char *dir) {
    int64_t MID = 128 + N, NW = (2 * N + 256) / 64 + 4;
    u64 **hist = malloc(sizeof(u64 *) * H);
    for (int64_t k = 0; k < H; k++) hist[k] = alloc(NW);
    u64 *row = hist[0];
    row[MID >> 6] |= 1ULL << (MID & 63);
    uint8_t *seq = malloc(H);
    FILE *f = open_or_die(dir, "wedge.i64");
    for (int64_t n = 0; n < N; n++) {
        row = hist[n % H];
        if (n >= H && n % stride == 0) {
            int64_t D = 2 * (n - H + 1) + 1;  // diagonals present in all H rows
            int64_t hdr[2] = {n, D};
            fwrite(hdr, sizeof hdr, 1, f);
            int16_t *per = malloc(sizeof(int16_t) * D);
            for (int64_t d = 0; d < D; d++) {
                for (int64_t k = 0; k < H; k++) {
                    int64_t m = n - H + 1 + k;
                    seq[k] = (uint8_t)bit(hist[m % H], MID - m + d);
                }
                int16_t best = 0;
                for (int64_t p = 1; p <= H / 2 && !best; p++) {
                    int ok = 1;
                    for (int64_t k = p; k < H && ok; k++) ok = seq[k] == seq[k - p];
                    if (ok) best = (int16_t)p;
                }
                per[d] = best;
            }
            fwrite(per, sizeof(int16_t), D, f);
            free(per);
            fflush(f);
            fprintf(stderr, "wedge checkpoint %lld\n", (long long)n);
        }
        u64 *nxt = hist[(n + 1) % H];
        for (int64_t w = (MID - n - 1) >> 6; w <= (MID + n + 1) >> 6; w++)
            nxt[w] = lsh(row, w) ^ (row[w] | rsh(row, w));
    }
    fclose(f);
    return 0;
}

// ---- rightmost black run per row (OEIS A094603, offset 1 = our row 0) ---------------------
static int run_rightrun(int64_t N, const char *dir) {
    int64_t MID = 128 + N, NW = (2 * N + 256) / 64 + 4;
    u64 *row = alloc(NW), *nxt = alloc(NW);
    row[MID >> 6] |= 1ULL << (MID & 63);
    uint16_t *out = malloc(sizeof(uint16_t) * N);
    for (int64_t n = 0; n < N; n++) {
        int64_t k = 0;
        while (bit(row, MID + n - k)) k++;
        out[n] = (uint16_t)k;
        for (int64_t w = (MID - n - 1) >> 6; w <= (MID + n + 1) >> 6; w++)
            nxt[w] = lsh(row, w) ^ (row[w] | rsh(row, w));
        u64 *x = row; row = nxt; nxt = x;
    }
    FILE *f = open_or_die(dir, "rightrun.u16"); fwrite(out, sizeof(uint16_t), N, f); fclose(f);
    return 0;
}

// ---- frontier: window-free edge of the frozen wedge ---------------------------------------
// Every wedge diagonal has a period dividing 32, so diagonal d breaks at row m iff
// cell(m, d) != cell(m-32, d), with d the distance from the left edge. L(m) is the smallest
// such d: the leftmost non-periodic cell of row m. Writes frontier.i64: L(m) for m = 0..N-1
// (-1 for m < 32).
#define FP 32
static int run_frontier(int64_t N, const char *dir) {
    int64_t MID = 128 + N, NW = (2 * N + 256) / 64 + 4;
    u64 *hist[FP + 1];
    for (int k = 0; k <= FP; k++) hist[k] = alloc(NW);
    hist[0][MID >> 6] |= 1ULL << (MID & 63);
    int64_t *L = malloc(sizeof(int64_t) * N);
    for (int64_t m = 0; m < N; m++) {
        u64 *row = hist[m % (FP + 1)];
        L[m] = -1;
        if (m >= FP) {
            u64 *old = hist[(m - FP) % (FP + 1)];   // cell d of row m-32 sits 32 bits higher
            for (int64_t w = (MID - m) >> 6; w <= (MID + m) >> 6; w++) {
                u64 v = row[w] ^ ((old[w] >> FP) | (old[w + 1] << (64 - FP)));
                if (v) { L[m] = w * 64 + __builtin_ctzll(v) - (MID - m); break; }
            }
        }
        u64 *nxt = hist[(m + 1) % (FP + 1)];
        for (int64_t w = (MID - m - 1) >> 6; w <= (MID + m + 1) >> 6; w++)
            nxt[w] = lsh(row, w) ^ (row[w] | rsh(row, w));
        if (m % 1000000 == 0) fprintf(stderr, "frontier row %lld L=%lld\n", (long long)m, (long long)L[m]);
    }
    FILE *f = open_or_die(dir, "frontier.i64"); fwrite(L, sizeof(int64_t), N, f); fclose(f);
    return 0;
}

// ---- damage: left spreading speed of a one-cell flip in a random row -----------------------
// Ring of W = 2T+128 cells from fair random bits (splitmix64), flip the cell at 0.3 W (the left
// front needs ~T/4 of room, the right front, at speed 1, needs T, so it never wraps), evolve
// both copies T steps and record how far left of the flip the leftmost difference has moved.
// Writes damage_<seed>.i32: the front distance after each step.
static u64 splitmix(u64 *s) {
    u64 z = (*s += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static void ring_step(u64 *x, u64 *y, int64_t NW) {
    x[0] = x[NW - 2]; x[NW - 1] = x[1];             // ghost words: words 1..NW-2 form the ring
    for (int64_t w = 1; w < NW - 1; w++) y[w] = lsh(x, w) ^ (x[w] | rsh(x, w));
}
static int run_damage(int64_t T, u64 seed, const char *dir) {
    int64_t NW = (2 * T + 128) / 64 + 2, W = (NW - 2) * 64;
    u64 *a = alloc(NW), *b = alloc(NW), *a2 = alloc(NW), *b2 = alloc(NW), s = seed;
    for (int64_t w = 1; w < NW - 1; w++) a[w] = b[w] = splitmix(&s);
    int64_t c = 64 + 3 * W / 10;                    // bit index, word 1 holds ring cells 0..63
    b[c >> 6] ^= 1ULL << (c & 63);
    int32_t *front = malloc(sizeof(int32_t) * T);
    int64_t lo = c;                                 // leftmost difference can move left by <= 1 per step
    for (int64_t t = 0; t < T; t++) {
        ring_step(a, a2, NW); ring_step(b, b2, NW);
        u64 *x = a; a = a2; a2 = x; x = b; b = b2; b2 = x;
        int64_t w = (lo - 1) >> 6;
        u64 v = (a[w] ^ b[w]) & (~0ULL << ((lo - 1) & 63));
        while (!v) { w++; v = a[w] ^ b[w]; }
        lo = w * 64 + __builtin_ctzll(v);
        front[t] = (int32_t)(c - lo);
    }
    char name[64]; snprintf(name, sizeof name, "damage_%llu.i32", (unsigned long long)seed);
    FILE *f = open_or_die(dir, name); fwrite(front, sizeof(int32_t), T, f); fclose(f);
    return 0;
}

int main(int argc, char **argv) {
    if (argc >= 4 && !strcmp(argv[1], "frontier")) return run_frontier(atoll(argv[2]), argv[3]);
    if (argc >= 5 && !strcmp(argv[1], "damage")) return run_damage(atoll(argv[2]), strtoull(argv[3], 0, 10), argv[4]);
    if (argc >= 4 && !strcmp(argv[1], "trace")) return run_trace(atoll(argv[2]), argv[3]);
    if (argc >= 4 && !strcmp(argv[1], "rightrun")) return run_rightrun(atoll(argv[2]), argv[3]);
    if (argc >= 6 && !strcmp(argv[1], "wedge")) return run_wedge(atoll(argv[2]), atoll(argv[3]), atoll(argv[4]), argv[5]);
    if (argc >= 6 && !strcmp(argv[1], "edge")) return run_edge(atoll(argv[2]), atoll(argv[3]), atoll(argv[4]), argv[5]);
    if (argc >= 7 && !strcmp(argv[1], "random"))
        return run_random(atoll(argv[2]), atoll(argv[3]), atoll(argv[4]), strtoull(argv[5], 0, 10), argv[6]);
    fprintf(stderr, "usage: %s trace N dir | random W T warm seed dir\n", argv[0]);
    return 2;
}
