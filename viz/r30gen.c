// Tile generator for the Rule 30 viewer.
//
//   r30gen N K B LMIN HB out_dir
//
// Evolves N rows from a single black cell and writes into out_dir:
//   ck/ck_<n>.u64      row n as raw words (absolute word (MID-n)>>6 .. (MID+n)>>6), every K rows
//   level_<B>.bin      dense block grid, B x B cells per block (B a power of two, multiple of 64):
//                      NR = ceil(N/B) block rows x NC block columns, column c covering s in [cB, cB+B),
//                      c from CMIN; record = u32 cells, ones, pcells, heads, maxlen (20 bytes)
//                      cells: cells inside the triangle; ones: black cells; pcells: cells inside a pyramid
//                      (membership P of that row); heads: chains whose first 111 centre is in the block;
//                      maxlen: longest such chain
//   level_<2B>.bin ... coarser levels, each a 2x2 reduction, down to one block row or B = 32768
//   hist_<HB>.u16      chain-length histogram (64 bins, last = 63+) per HB x HB block of chain heads, saturating
//   chains.i32         (row, start offset, length) for every chain of length >= LMIN, by row
//   rows.i64           per row: n, centre bit, ones, new chains, active pyramids (rule30.c F_ACTIVE)
//   meta.json
#include <fcntl.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#include "r30core.h"

typedef struct { uint32_t cells, ones, pcells, heads, maxlen; } blk_t;
#define HBINS 64

static void die(const char *m) { perror(m); exit(1); }

static int open_out(const char *dir, const char *name) {
    char p[4096];
    snprintf(p, sizeof p, "%s/%s", dir, name);
    int fd = open(p, O_RDWR | O_CREAT | O_TRUNC, 0644);
    if (fd < 0) die(p);
    return fd;
}

static void pwrite_all(int fd, const void *buf, size_t len, off_t off) {
    const char *b = buf;
    while (len) {
        ssize_t k = pwrite(fd, b, len, off);
        if (k <= 0) die("pwrite");
        b += k; len -= k; off += k;
    }
}

// 2x2 reduction of a dense level; returns the new geometry through the pointers
static void reduce(const char *dir, int64_t B, int64_t NR, int64_t NC, int64_t CMIN,
                   int64_t *NR2, int64_t *NC2, int64_t *CMIN2) {
    char name[64];
    snprintf(name, sizeof name, "level_%lld.bin", (long long)B);
    char p[4096];
    snprintf(p, sizeof p, "%s/%s", dir, name);
    int in = open(p, O_RDONLY);
    if (in < 0) die(p);
    int64_t cmin2 = floordiv(CMIN, 2), cmax2 = floordiv(CMIN + NC - 1, 2);
    int64_t nr2 = (NR + 1) / 2, nc2 = cmax2 - cmin2 + 1;
    snprintf(name, sizeof name, "level_%lld.bin", (long long)(2 * B));
    int out = open_out(dir, name);
    if (ftruncate(out, nr2 * nc2 * (off_t)sizeof(blk_t))) die("ftruncate");
    blk_t *a = malloc(NC * sizeof(blk_t)), *o = malloc(nc2 * sizeof(blk_t));
    for (int64_t R = 0; R < nr2; R++) {
        memset(o, 0, nc2 * sizeof(blk_t));
        for (int64_t h = 0; h < 2; h++) {
            int64_t r = 2 * R + h;
            if (r >= NR) break;
            if (pread(in, a, NC * sizeof(blk_t), r * NC * (off_t)sizeof(blk_t)) != (ssize_t)(NC * sizeof(blk_t))) die("pread");
            for (int64_t j = 0; j < NC; j++) {
                blk_t *d = &o[floordiv(CMIN + j, 2) - cmin2], *s = &a[j];
                d->cells += s->cells; d->ones += s->ones; d->pcells += s->pcells; d->heads += s->heads;
                if (s->maxlen > d->maxlen) d->maxlen = s->maxlen;
            }
        }
        pwrite_all(out, o, nc2 * sizeof(blk_t), R * nc2 * (off_t)sizeof(blk_t));
    }
    close(in); close(out); free(a); free(o);
    *NR2 = nr2; *NC2 = nc2; *CMIN2 = cmin2;
}

int main(int argc, char **argv) {
    if (argc < 7) { fprintf(stderr, "usage: %s N K B LMIN HB out_dir\n", argv[0]); return 2; }
    int64_t N = atoll(argv[1]), K = atoll(argv[2]), B = atoll(argv[3]), LMIN = atoll(argv[4]), HB = atoll(argv[5]);
    const char *dir = argv[6];
    if (B < 64 || (B & (B - 1))) { fprintf(stderr, "B must be a power of two >= 64\n"); return 2; }
    int64_t MID = ((N + 256) / 64 + 1) * 64, NW = (MID + N + 256) / 64 + 2;
    int64_t WPB = B / 64;  // words per block column

    char ckdir[4096];
    snprintf(ckdir, sizeof ckdir, "%s/ck", dir);
    mkdir(dir, 0755);
    mkdir(ckdir, 0755);

    u64 *row = alloc64(NW), *nxt = alloc64(NW), *is111 = alloc64(NW), *is000 = alloc64(NW), *in2 = alloc64(NW),
        *P = alloc64(NW), *Pn = alloc64(NW);
    row[MID >> 6] |= 1ULL << (MID & 63);

    int64_t NR = (N + B - 1) / B, CMIN = floordiv(-(N - 1), B), NC = floordiv(N - 1, B) - CMIN + 1;
    int lv = open_out(dir, "level_tmp");
    char lvname[64];
    snprintf(lvname, sizeof lvname, "level_%lld.bin", (long long)B);
    if (ftruncate(lv, NR * NC * (off_t)sizeof(blk_t))) die("ftruncate");
    blk_t *acc = calloc(NC, sizeof(blk_t));

    int64_t HR = (N + HB - 1) / HB, HCMIN = floordiv(-(N - 1), HB), HC = floordiv(N - 1, HB) - HCMIN + 1;
    int hf = open_out(dir, "hist_tmp");
    char hname[64];
    snprintf(hname, sizeof hname, "hist_%lld.u16", (long long)HB);
    if (ftruncate(hf, HR * HC * HBINS * (off_t)2)) die("ftruncate");
    uint32_t *hacc = calloc(HC * HBINS, sizeof(uint32_t));
    uint16_t *hrow = calloc(HC * HBINS, sizeof(uint16_t));

    char p[4096];
    snprintf(p, sizeof p, "%s/chains.i32", dir);
    FILE *fch = fopen(p, "wb");
    snprintf(p, sizeof p, "%s/rows.i64", dir);
    FILE *frows = fopen(p, "wb");
    if (!fch || !frows) die("fopen");
    int64_t total_hist[HBINS] = {0}, nchains_logged = 0;
    time_t t0 = time(NULL);

    for (int64_t n = 0; n < N; n++) {
        int64_t wlo = (MID - n - 1) >> 6, whi = (MID + n + 1) >> 6;
        int64_t ones = 0, newp = 0, active = 0;
        if (n % K == 0) {
            snprintf(p, sizeof p, "%s/ck_%lld.u64", ckdir, (long long)n);
            FILE *f = fopen(p, "wb");
            if (!f) die(p);
            int64_t a = (MID - n) >> 6, b = (MID + n) >> 6;
            fwrite(row + a, sizeof(u64), b - a + 1, f);
            fclose(f);
        }
        gen_masks(row, nxt, is111, is000, wlo, whi);
        gen_in2(is111, in2, wlo, whi);
        int64_t cur = 0;
        for (int64_t w = wlo; w <= whi; w++) {
            u64 l = lsh(in2, w), seed = in2[w] | l;
            Pn[w] = seed | (is000[w] & P[w]);
            int64_t c = floordiv(w * 64 - MID, B) - CMIN;
            u64 m = span(w, MID - n, MID + n);
            if (m) {
                blk_t *k = &acc[c];
                int pc = __builtin_popcountll(row[w]);
                k->cells += __builtin_popcountll(m);
                k->ones += pc;
                k->pcells += __builtin_popcountll(P[w]);
                ones += pc;
            }
            u64 starts = in2[w] & ~l, ends = in2[w] & ~rsh(in2, w), ev = starts | ends;
            newp += __builtin_popcountll(starts);
            while (ev) {
                int j = __builtin_ctzll(ev);
                ev &= ev - 1;
                int64_t pos = w * 64 + j;
                if ((starts >> j) & 1) cur = pos;
                if ((ends >> j) & 1) {
                    int64_t len = pos - cur + 1, s = cur - MID, L = len < HBINS ? len : HBINS - 1;
                    blk_t *k = &acc[floordiv(s, B) - CMIN];
                    k->heads++;
                    if (len > k->maxlen) k->maxlen = (uint32_t)len;
                    hacc[(floordiv(s, HB) - HCMIN) * HBINS + L]++;
                    total_hist[L]++;
                    if (len >= LMIN) {
                        int32_t rec[3] = {(int32_t)n, (int32_t)s, (int32_t)len};
                        fwrite(rec, sizeof rec, 1, fch);
                        nchains_logged++;
                    }
                }
            }
        }
        for (int64_t w = wlo; w <= whi; w++) active += __builtin_popcountll(Pn[w] & ~lsh(Pn, w));
        int64_t rec[5] = {n, bit(row, MID), ones, newp, active};
        fwrite(rec, sizeof rec, 1, frows);

        if (n % B == B - 1 || n == N - 1) {
            pwrite_all(lv, acc, NC * sizeof(blk_t), (n / B) * NC * (off_t)sizeof(blk_t));
            memset(acc, 0, NC * sizeof(blk_t));
        }
        if (n % HB == HB - 1 || n == N - 1) {
            for (int64_t i = 0; i < HC * HBINS; i++) hrow[i] = hacc[i] > 65535 ? 65535 : (uint16_t)hacc[i];
            pwrite_all(hf, hrow, HC * HBINS * 2, (n / HB) * HC * HBINS * (off_t)2);
            memset(hacc, 0, HC * HBINS * sizeof(uint32_t));
        }
        u64 *t = row; row = nxt; nxt = t;
        t = P; P = Pn; Pn = t;
        if (n % 50000 == 0) {
            fflush(fch); fflush(frows);
            fprintf(stderr, "row %lld/%lld  %lds\n", (long long)n, (long long)N, (long)(time(NULL) - t0));
        }
    }
    fclose(fch); fclose(frows); close(lv); close(hf);
    char a[4096], b2[4096];
    snprintf(a, sizeof a, "%s/level_tmp", dir); snprintf(b2, sizeof b2, "%s/%s", dir, lvname); rename(a, b2);
    snprintf(a, sizeof a, "%s/hist_tmp", dir); snprintf(b2, sizeof b2, "%s/%s", dir, hname); rename(a, b2);
    fprintf(stderr, "evolved %lld rows in %lds; reducing\n", (long long)N, (long)(time(NULL) - t0));

    snprintf(p, sizeof p, "%s/meta.json", dir);
    FILE *fm = fopen(p, "w");
    fprintf(fm, "{\"N\": %lld, \"K\": %lld, \"MID\": %lld, \"LMIN\": %lld, \"HB\": %lld, \"HR\": %lld, \"HC\": %lld, "
                "\"HCMIN\": %lld, \"chains_logged\": %lld,\n \"hist\": [",
            (long long)N, (long long)K, (long long)MID, (long long)LMIN, (long long)HB, (long long)HR, (long long)HC,
            (long long)HCMIN, (long long)nchains_logged);
    for (int i = 0; i < HBINS; i++) fprintf(fm, "%s%lld", i ? ", " : "", (long long)total_hist[i]);
    fprintf(fm, "],\n \"levels\": [");
    int first = 1;
    for (int64_t b = B, nr = NR, nc = NC, cm = CMIN;;) {
        fprintf(fm, "%s\n  {\"B\": %lld, \"NR\": %lld, \"NC\": %lld, \"CMIN\": %lld}", first ? "" : ",", (long long)b,
                (long long)nr, (long long)nc, (long long)cm);
        first = 0;
        // u32 counts: a 2B x 2B block must hold < 2^32 cells
        if ((nr <= 1 && nc <= 2) || 2 * b > 32768) break;
        int64_t nr2, nc2, cm2;
        reduce(dir, b, nr, nc, cm, &nr2, &nc2, &cm2);
        b *= 2; nr = nr2; nc = nc2; cm = cm2;
    }
    fprintf(fm, "]}\n");
    fclose(fm);
    fprintf(stderr, "done in %lds\n", (long)(time(NULL) - t0));
    (void)WPB;
    return 0;
}
