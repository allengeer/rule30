// Bit-parallel Rule 30 engine for deep runs (1M+ rows).
//
// Same definitions as sim.py: rows from a single black cell, columns -n..n, a pyramid is a
// maximal chain of >= 2 adjacent 111 neighbourhoods; its seed cells are the cells under the
// chain plus the 110 cell to its right, and it lives on through 000 neighbourhoods.
//
//   rule30 stats  N out_dir [snap_rows...]   per-row records + per-column snapshots
//   rule30 period N out_dir shard shards     smallest period of each row n with n % shards == shard
//
// stats writes out_dir/rows.i64: 20 int64 fields per row (see FIELDS below), plus
// out_dir/hist.i64 (pyramid chain-length histogram) and out_dir/snap_<n>.bin
// (float32 E[2n+1] then uint32 pyramids[2n+1]) for each requested snapshot row.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef uint64_t u64;
#define PLANES 24
#define HIST 256
enum { F_ROW, F_CBIT, F_CONES, F_CPERIOD, F_CPYR, F_CINPYR, F_NEWPYR, F_ACTIVE, F_TOTAL,
       F_ONES, F_ONESL, F_ONESR, F_TRIPLE0, NFIELDS = F_TRIPLE0 + 8 };

static int64_t N, NW, MID;

static u64 *alloc(int64_t n) { u64 *p = calloc(n, sizeof(u64)); if (!p) { perror("calloc"); exit(1); } return p; }
static inline u64 lsh(const u64 *x, int64_t w) { return (x[w] << 1) | (x[w - 1] >> 63); }  // left neighbour's value
static inline u64 rsh(const u64 *x, int64_t w) { return (x[w] >> 1) | (x[w + 1] << 63); }  // right neighbour's value
static inline int bit(const u64 *x, int64_t i) { return (x[i >> 6] >> (i & 63)) & 1; }

// cells lo..hi (inclusive) within word w
static inline u64 span(int64_t w, int64_t lo, int64_t hi) {
    int64_t a = w * 64, b = a + 63;
    if (hi < a || lo > b) return 0;
    if (lo < a) lo = a;
    if (hi > b) hi = b;
    int64_t len = hi - lo + 1;
    return (len == 64 ? ~0ULL : ((1ULL << len) - 1)) << (lo - a);
}

// bit-sliced per-column counters
static void counter_add(u64 **planes, int64_t w, u64 m) {
    for (int k = 0; m && k < PLANES; k++) {
        u64 t = planes[k][w] & m;
        planes[k][w] ^= m;
        m = t;
    }
}
static uint32_t counter_get(u64 **planes, int64_t i) {
    uint32_t v = 0;
    for (int k = 0; k < PLANES; k++) v |= (uint32_t)bit(planes[k], i) << k;
    return v;
}

static int64_t smallest_period(const uint8_t *s, int64_t n, int32_t *fail) {
    int64_t k = 0;
    fail[0] = 0;
    for (int64_t i = 1; i < n; i++) {
        while (k && s[i] != s[k]) k = fail[k - 1];
        if (s[i] == s[k]) k++;
        fail[i] = (int32_t)k;
    }
    return n - fail[n - 1];
}

static void step(u64 *row, u64 *nxt, int64_t wlo, int64_t whi) {
    for (int64_t w = wlo; w <= whi; w++) nxt[w] = lsh(row, w) ^ (row[w] | rsh(row, w));
}

static void open_or_die(FILE **f, const char *dir, const char *name, const char *mode) {
    char path[4096];
    snprintf(path, sizeof path, "%s/%s", dir, name);
    if (!(*f = fopen(path, mode))) { perror(path); exit(1); }
}

static int run_period(const char *dir, int64_t shard, int64_t shards) {
    u64 *row = alloc(NW), *nxt = alloc(NW);
    uint8_t *s = malloc(2 * N + 1);
    int32_t *fail = malloc(sizeof(int32_t) * (2 * N + 1));
    row[MID >> 6] |= 1ULL << (MID & 63);
    char name[64];
    snprintf(name, sizeof name, "period_%lld_of_%lld.i64", (long long)shard, (long long)shards);
    FILE *f;
    open_or_die(&f, dir, name, "wb");
    for (int64_t n = 0; n < N; n++) {
        if (n % shards == shard) {
            for (int64_t i = 0; i <= 2 * n; i++) s[i] = bit(row, MID - n + i);
            int64_t rec[2] = {n, smallest_period(s, 2 * n + 1, fail)};
            fwrite(rec, sizeof rec, 1, f);
        }
        int64_t wlo = (MID - n - 1) >> 6, whi = (MID + n + 1) >> 6;
        step(row, nxt, wlo, whi);
        u64 *t = row; row = nxt; nxt = t;
    }
    fclose(f);
    return 0;
}

static int cmp_i64(const void *a, const void *b) {
    int64_t x = *(const int64_t *)a, y = *(const int64_t *)b;
    return (x > y) - (x < y);
}

static int run_stats(const char *dir, int64_t *snaps, int nsnaps) {
    u64 *row = alloc(NW), *nxt = alloc(NW), *is111 = alloc(NW), *is000 = alloc(NW),
        *in2 = alloc(NW), *seed = alloc(NW), *P = alloc(NW), *Pn = alloc(NW);
    u64 *ones[PLANES], *pyr[PLANES];
    for (int k = 0; k < PLANES; k++) { ones[k] = alloc(NW); pyr[k] = alloc(NW); }
    uint8_t *cbits = malloc(N);
    int32_t *cfail = malloc(sizeof(int32_t) * N);
    int64_t ck = 0, cones = 0, cpyr = 0, total = 0, hist[HIST] = {0};
    row[MID >> 6] |= 1ULL << (MID & 63);
    qsort(snaps, nsnaps, sizeof *snaps, cmp_i64);
    int si = 0;

    FILE *fr;
    open_or_die(&fr, dir, "rows.i64", "wb");
    time_t t0 = time(NULL);

    for (int64_t n = 0; n < N; n++) {
        int64_t rec[NFIELDS] = {0};
        int64_t wlo = (MID - n - 1) >> 6, whi = (MID + n + 1) >> 6;
        int64_t clo = MID - n, chi = MID + n;

        // centre column
        int b = bit(row, MID);
        cbits[n] = (uint8_t)b;
        cones += b;
        if (n) {
            while (ck && cbits[n] != cbits[ck]) ck = cfail[ck - 1];
            if (cbits[n] == cbits[ck]) ck++;
        }
        cfail[n] = (int32_t)ck;
        rec[F_ROW] = n; rec[F_CBIT] = b; rec[F_CONES] = cones; rec[F_CPERIOD] = n + 1 - ck;
        rec[F_CPYR] = cpyr; rec[F_CINPYR] = bit(P, MID);

        // pass 1: neighbourhoods, next row, densities, triples, column ones
        for (int64_t w = wlo; w <= whi; w++) {
            u64 L = lsh(row, w), C = row[w], R = rsh(row, w), m = span(w, clo, chi);
            is111[w] = L & C & R;
            is000[w] = ~L & ~C & ~R;
            nxt[w] = L ^ (C | R);
            rec[F_ONES] += __builtin_popcountll(C);
            rec[F_ONESL] += __builtin_popcountll(C & span(w, clo, MID - 1));
            rec[F_ONESR] += __builtin_popcountll(C & span(w, MID + 1, chi));
            for (int t = 0; t < 8; t++) {
                u64 x = ((t & 4) ? L : ~L) & ((t & 2) ? C : ~C) & ((t & 1) ? R : ~R);
                rec[F_TRIPLE0 + t] += __builtin_popcountll(x & m);
            }
            counter_add(ones, w, C);
        }
        // pass 2: chains of >= 2 adjacent 111s
        for (int64_t w = wlo; w <= whi; w++) in2[w] = is111[w] & (lsh(is111, w) | rsh(is111, w));
        // pass 3: seeds, membership of next row, chain-length histogram
        int64_t cur = 0;
        for (int64_t w = wlo; w <= whi; w++) {
            u64 l = lsh(in2, w);
            seed[w] = in2[w] | l;
            Pn[w] = seed[w] | (is000[w] & P[w]);
            u64 starts = in2[w] & ~l, ends = in2[w] & ~rsh(in2, w), ev = starts | ends;
            rec[F_NEWPYR] += __builtin_popcountll(starts);
            while (ev) {
                int j = __builtin_ctzll(ev);
                ev &= ev - 1;
                int64_t pos = w * 64 + j;
                if ((starts >> j) & 1) cur = pos;
                if ((ends >> j) & 1) { int64_t len = pos - cur + 1; hist[len < HIST ? len : HIST - 1]++; }
            }
        }
        for (int64_t w = wlo; w <= whi; w++) rec[F_ACTIVE] += __builtin_popcountll(Pn[w] & ~lsh(Pn, w));

        // snapshot (ones include row n; pyramid counts exclude seeds formed from row n)
        while (si < nsnaps && snaps[si] < n) si++;
        if (si < nsnaps && snaps[si] == n) {
            char name[64];
            snprintf(name, sizeof name, "snap_%lld.bin", (long long)n);
            FILE *fs;
            open_or_die(&fs, dir, name, "wb");
            float *E = malloc(sizeof(float) * (2 * n + 1));
            uint32_t *pc = malloc(sizeof(uint32_t) * (2 * n + 1));
            for (int64_t c = -n; c <= n; c++) {
                int64_t m = n - (c < 0 ? -c : c) + 1;
                E[c + n] = (float)counter_get(ones, MID + c) / (float)m;
                pc[c + n] = counter_get(pyr, MID + c);
            }
            fwrite(E, sizeof(float), 2 * n + 1, fs);
            fwrite(pc, sizeof(uint32_t), 2 * n + 1, fs);
            fclose(fs);
            free(E); free(pc);
        }

        for (int64_t w = wlo; w <= whi; w++) counter_add(pyr, w, seed[w]);
        cpyr += bit(seed, MID);
        total += rec[F_NEWPYR];
        rec[F_TOTAL] = total;
        fwrite(rec, sizeof rec, 1, fr);

        u64 *t = row; row = nxt; nxt = t;
        t = P; P = Pn; Pn = t;
        if (n % 50000 == 0) {
            fflush(fr);
            fprintf(stderr, "row %lld/%lld  %lds  centre_E=%.5f pyramids=%lld\n", (long long)n,
                    (long long)N, (long)(time(NULL) - t0), (double)cones / (n + 1), (long long)total);
        }
    }
    fclose(fr);
    FILE *fh;
    open_or_die(&fh, dir, "hist.i64", "wb");
    fwrite(hist, sizeof hist, 1, fh);
    fclose(fh);
    fprintf(stderr, "done: %lld rows in %lds\n", (long long)N, (long)(time(NULL) - t0));
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "usage: %s stats N dir [snap...] | period N dir shard shards\n", argv[0]);
        return 2;
    }
    N = atoll(argv[2]);
    MID = 128 + N;
    NW = (2 * N + 256) / 64 + 4;
    if (!strcmp(argv[1], "period")) return run_period(argv[3], atoll(argv[4]), atoll(argv[5]));
    int ns = argc - 4;
    int64_t *snaps = malloc(sizeof(int64_t) * (ns ? ns : 1));
    for (int i = 0; i < ns; i++) snaps[i] = atoll(argv[4 + i]);
    return run_stats(argv[3], snaps, ns);
}
