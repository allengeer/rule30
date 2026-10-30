// Prepares the 10^9-bit centre column for the strip view.
//
//   centre_prep billion.mx offset nbits out_dir
//
// The Wolfram Data Repository MX stores the ByteArray uncompressed, highest bit first, at byte `offset`.
// Writes, with W(p) = 2 * (ones among the first p bits) - p (so the running excess E = W/2):
//   centre.bits            the raw packed bytes
//   lvl_<k>.i32            for blocks of 2^k bits, k = 12..30: min W, max W over p in (start, end], W(end)
//   runs16.u32             per 2^16-bit block: [64] black-run and [64] white-run counts by length (63 = 63+),
//                          runs counted in the block where they start, at full length
//   ties.i64               every p >= 1 with W(p) = 0
//   summary.json
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define K0 12
#define K1 30
#define RB 16

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: %s mx offset nbits out_dir\n", argv[0]); return 2; }
    long off = atol(argv[2]);
    int64_t nbits = atoll(argv[3]);
    const char *dir = argv[4];
    FILE *f = fopen(argv[1], "rb");
    if (!f) { perror(argv[1]); return 1; }
    int64_t nbytes = (nbits + 7) / 8;
    uint8_t *by = malloc(nbytes);
    fseek(f, off, SEEK_SET);
    if ((int64_t)fread(by, 1, nbytes, f) != nbytes) { fprintf(stderr, "short read\n"); return 1; }
    fclose(f);
    char p[4096];
    snprintf(p, sizeof p, "%s/centre.bits", dir);
    f = fopen(p, "wb"); fwrite(by, 1, nbytes, f); fclose(f);

    int64_t nb0 = (nbits + (1LL << K0) - 1) >> K0;
    int32_t *lv = malloc(nb0 * 3 * sizeof(int32_t));
    int64_t nrb = (nbits + (1LL << RB) - 1) >> RB;
    uint32_t *runs = calloc(nrb * 128, sizeof(uint32_t));
    snprintf(p, sizeof p, "%s/ties.i64", dir);
    FILE *ft = fopen(p, "wb");
    int64_t W = 0, ones = 0, ties = 0, last_tie = 0, run_start = 0, maxW = 0, minW_after_tie = 0;
    int prev = -1;
    for (int64_t i = 0; i < nbits; i++) {
        int b = (by[i >> 3] >> (7 - (i & 7))) & 1;
        if (b != prev) {
            if (prev >= 0) {
                int64_t L = i - run_start < 63 ? i - run_start : 63;
                runs[(run_start >> RB) * 128 + (prev ? 0 : 64) + L]++;
            }
            run_start = i;
            prev = b;
        }
        ones += b;
        W += b ? 1 : -1;
        int64_t blk = i >> K0;
        int32_t *q = lv + blk * 3;
        if ((i & ((1LL << K0) - 1)) == 0) q[0] = q[1] = (int32_t)W;
        if (W < q[0]) q[0] = (int32_t)W;
        if (W > q[1]) q[1] = (int32_t)W;
        q[2] = (int32_t)W;
        if (W > maxW) maxW = W;
        if (W == 0) { int64_t pp = i + 1; fwrite(&pp, 8, 1, ft); ties++; last_tie = pp; }
    }
    {
        int64_t L = nbits - run_start < 63 ? nbits - run_start : 63;
        runs[(run_start >> RB) * 128 + (prev ? 0 : 64) + L]++;
    }
    fclose(ft);
    (void)minW_after_tie;
    snprintf(p, sizeof p, "%s/runs16.u32", dir);
    f = fopen(p, "wb"); fwrite(runs, sizeof(uint32_t), nrb * 128, f); fclose(f);

    int64_t nb = nb0;
    for (int k = K0; k <= K1; k++) {
        snprintf(p, sizeof p, "%s/lvl_%d.i32", dir, k);
        f = fopen(p, "wb"); fwrite(lv, sizeof(int32_t), nb * 3, f); fclose(f);
        int64_t nb2 = (nb + 1) / 2;
        for (int64_t j = 0; j < nb2; j++) {
            int32_t *a = lv + 2 * j * 3, *d = lv + j * 3;
            int32_t mn = a[0], mx = a[1], e = a[2];
            if (2 * j + 1 < nb) {
                int32_t *c = a + 3;
                if (c[0] < mn) mn = c[0];
                if (c[1] > mx) mx = c[1];
                e = c[2];
            }
            d[0] = mn; d[1] = mx; d[2] = e;
        }
        nb = nb2;
    }
    snprintf(p, sizeof p, "%s/summary.json", dir);
    f = fopen(p, "w");
    fprintf(f, "{\"nbits\": %lld, \"ones\": %lld, \"excess\": %.1f, \"ties\": %lld, \"last_tie\": %lld, "
               "\"max_excess\": %.1f, \"k0\": %d, \"k1\": %d, \"run_block\": %d}\n",
            (long long)nbits, (long long)ones, W / 2.0, (long long)ties, (long long)last_tie, maxW / 2.0, K0, K1, RB);
    fclose(f);
    fprintf(stderr, "ones %lld excess %.1f ties %lld last tie %lld\n", (long long)ones, W / 2.0, (long long)ties,
            (long long)last_tie);
    return 0;
}
