// On-demand window renderer and centre-column helpers, loaded by server.py through ctypes.
//
// r30_window evolves from a checkpoint row c (the full row, as written by r30gen) down to n1 over a column
// band wide enough to cover the light cone of the window rows [n0, n1) x columns [s0, s1), and reports:
//   b == 1: per cell  cell[]: bit0 black, bit1 chain cell (a 111 centre inside a chain), bits 2..7 pyramid
//                     label = min(length, 63) of the chain whose pyramid holds this cell (0 = none)
//                     aux[]:  bit0 a pattern match starts here, bit1 the cell is covered by a match
//   b  > 1: per b x b block  6 u32: cells, ones, pcells, heads, maxlen, pattern matches (starts)
//   hist[64]: lengths of chains whose first 111 centre lies in the window (bin 63 = 63+)
// Pyramid membership starts empty at row c, so c must be at least 128 rows above n0 (or 0): a pyramid
// shrinks by >= 1 cell per side per row, so this is exact for chains shorter than ~250.
#include <string.h>

#include "r30core.h"

static int64_t imax(int64_t a, int64_t b) { return a > b ? a : b; }
static int64_t imin(int64_t a, int64_t b) { return a < b ? a : b; }

int r30_window(const u64 *ck, int64_t ck_wlo, int64_t ck_nw, int64_t MID, int64_t c, int64_t n0, int64_t n1,
               int64_t s0, int64_t s1, int64_t b, const char *pat, uint8_t *cell, uint8_t *aux, uint32_t *blk,
               double *hist) {
    if (c > n0 || n1 <= n0 || s1 <= s0 || b < 1) return -1;
    int k = pat ? (int)strlen(pat) : 0;
    if (k > 63) return -2;
    u64 pbits = 0;
    for (int j = 0; j < k; j++) pbits |= (u64)(pat[j] == '1') << j;

    int64_t reach = n1 - c + 192;
    int64_t A = imax(MID + s0 - reach, MID - n1 - 128), Z = imin(MID + s1 + reach, MID + n1 + 128);
    int64_t g0 = (A >> 6) - 2, g1 = (Z >> 6) + 2, G = g1 - g0 + 1;  // local words, absolute g0..g1
    u64 *buf = calloc(7 * (G + 2), sizeof(u64));
    if (!buf) return -3;
    // offset so x[w] is absolute word w
    u64 *row = buf + 1 - g0, *nxt = row + (G + 2), *is111 = nxt + (G + 2), *is000 = is111 + (G + 2),
        *in2 = is000 + (G + 2), *P = in2 + (G + 2), *Pn = P + (G + 2);
    for (int64_t w = imax(g0, ck_wlo); w <= imin(g1, ck_wlo + ck_nw - 1); w++) row[w] = ck[w - ck_wlo];

    int64_t W = s1 - s0, H = n1 - n0;
    int64_t a0 = MID + s0, a1 = MID + s1 - 1;          // window, absolute cells
    int64_t wa = a0 >> 6, wb = a1 >> 6;                  // window words
    int64_t lab_lo = c < n0 - 128 ? n0 - 128 : c;        // rows from which labels are tracked
    uint8_t *lab = b == 1 ? calloc(W, 1) : NULL, *labn = b == 1 ? calloc(W, 1) : NULL;
    int64_t nbr = (H + b - 1) / b, nbc = (W + b - 1) / b;
    if (b > 1) memset(blk, 0, nbr * nbc * 6 * sizeof(uint32_t));
    else { memset(cell, 0, H * W); memset(aux, 0, H * W); }
    memset(hist, 0, 64 * sizeof(double));
    u64 *M = calloc(wb - wa + 3, sizeof(u64));          // pattern starts, words wa-1..wb+1

    for (int64_t r = c; r < n1; r++) {
        gen_masks(row, nxt, is111, is000, g0 + 1, g1 - 1);
        gen_in2(is111, in2, g0 + 2, g1 - 2);
        for (int64_t w = g0 + 2; w <= g1 - 2; w++) Pn[w] = in2[w] | lsh(in2, w) | (is000[w] & P[w]);
        int inwin = r >= n0;
        int64_t y = r - n0;

        if (b == 1 && r >= lab_lo) {
            // labels: carry through 000, then stamp seeds of every chain that touches the window
            for (int64_t i = 0; i < W; i++) labn[i] = bit(is000, a0 + i) ? lab[i] : 0;
        }
        // chain events over the window words plus a margin for chains reaching in from the left
        int64_t e0 = imax(g0 + 2, wa - 2), e1 = imin(g1 - 2, wb + 1), cur = -1;
        if (r >= lab_lo || inwin) {
            for (int64_t w = e0; w <= e1; w++) {
                u64 l = lsh(in2, w), starts = in2[w] & ~l, ends = in2[w] & ~rsh(in2, w), ev = starts | ends;
                while (ev) {
                    int j = __builtin_ctzll(ev);
                    ev &= ev - 1;
                    int64_t pos = w * 64 + j;
                    if ((starts >> j) & 1) cur = pos;
                    if (((ends >> j) & 1) && cur >= 0) {
                        int64_t len = pos - cur + 1, L = len < 63 ? len : 63;
                        if (b == 1 && r >= lab_lo)
                            for (int64_t i = imax(cur, a0); i <= imin(pos + 1, a1); i++) labn[i - a0] = (uint8_t)L;
                        if (inwin && b == 1)
                            for (int64_t i = imax(cur, a0); i <= imin(pos, a1); i++) cell[y * W + (i - a0)] |= 2;
                        if (inwin && cur >= a0 && cur <= a1) {
                            hist[L] += 1;
                            if (b > 1) {
                                uint32_t *q = blk + ((y / b) * nbc + (cur - a0) / b) * 6;
                                q[3]++;
                                if ((uint32_t)len > q[4]) q[4] = (uint32_t)len;
                            }
                        }
                    }
                }
            }
        }
        if (inwin) {
            // pattern starts: cells i..i+k-1 equal pat, the whole match inside row r
            if (k) {
                for (int64_t w = wa - 1; w <= wb + 1; w++) {
                    u64 m = span(w, MID - r, MID + r - k + 1);
                    for (int j = 0; j < k && m; j++) m &= ((pbits >> j) & 1) ? ahead(row, w, j) : ~ahead(row, w, j);
                    M[w - wa + 1] = m;
                }
            }
            if (b == 1) {
                for (int64_t i = 0; i < W; i++) {
                    int64_t a = a0 + i;
                    uint8_t v = (uint8_t)bit(row, a);
                    if (bit(P, a)) v |= (uint8_t)(lab[i] << 2);
                    cell[y * W + i] |= v;
                }
                if (k) {
                    for (int64_t w = wa; w <= wb; w++) {
                        u64 st = M[w - wa + 1], cov = 0;
                        for (int j = 0; j < k; j++) cov |= behind(M + 1 - wa, w, j);  // starts at i-j cover i
                        u64 mk = span(w, a0, a1);
                        st &= mk; cov &= mk;
                        while (st) { int j = __builtin_ctzll(st); st &= st - 1; aux[y * W + (w * 64 + j - a0)] |= 1; }
                        while (cov) { int j = __builtin_ctzll(cov); cov &= cov - 1; aux[y * W + (w * 64 + j - a0)] |= 2; }
                    }
                }
            } else {
                uint32_t *brow = blk + (y / b) * nbc * 6;
                for (int64_t w = wa; w <= wb; w++) {
                    u64 inside = span(w, MID - r, MID + r) & span(w, a0, a1);
                    if (!inside) continue;
                    u64 C = row[w], Pw = P[w], Mw = k ? M[w - wa + 1] : 0;
                    if (b >= 64) {
                        uint32_t *q = brow + ((w * 64 - a0) / b) * 6;
                        q[0] += __builtin_popcountll(inside);
                        q[1] += __builtin_popcountll(C & inside);
                        q[2] += __builtin_popcountll(Pw & inside);
                        q[5] += __builtin_popcountll(Mw & inside);
                    } else {
                        for (int64_t o = 0; o < 64; o += b) {
                            u64 sub = inside & (((b == 64 ? 0 : (1ULL << b)) - 1) << o);
                            if (!sub) continue;
                            uint32_t *q = brow + ((w * 64 + o - a0) / b) * 6;
                            q[0] += __builtin_popcountll(sub);
                            q[1] += __builtin_popcountll(C & sub);
                            q[2] += __builtin_popcountll(Pw & sub);
                            q[5] += __builtin_popcountll(Mw & sub);
                        }
                    }
                }
            }
        }
        if (b == 1 && r >= lab_lo) { uint8_t *t = lab; lab = labn; labn = t; }
        u64 *t = row; row = nxt; nxt = t;
        t = P; P = Pn; Pn = t;
        // the band edges go stale; keep them zero so the garbage stays outside the light cone
        row[g0] = row[g0 + 1] = row[g1] = row[g1 - 1] = 0;
    }
    free(buf); free(lab); free(labn); free(M);
    return 0;
}

// ---- centre column (bytes packed highest bit first, as in the Wolfram Data Repository) ------------------
static inline int cbit(const uint8_t *by, int64_t i) { return (by[i >> 3] >> (7 - (i & 7))) & 1; }

// Occurrences of pat starting at p in [lo, hi) (match may run past hi up to nbits). Fills up to maxpos
// positions and per-pixel counts over npx equal bins of [lo, hi). Returns the total count.
int64_t r30c_find(const uint8_t *by, int64_t nbits, int64_t lo, int64_t hi, const char *pat, int64_t *pos,
                  int64_t maxpos, int64_t *bins, int64_t npx) {
    int k = (int)strlen(pat);
    if (k < 1 || k > 63) return -1;
    u64 want = 0, mask = (1ULL << k) - 1, reg = 0;
    for (int j = 0; j < k; j++) want = (want << 1) | (u64)(pat[j] == '1');
    int64_t cnt = 0, end = imin(hi + k - 1, nbits);
    if (bins) memset(bins, 0, npx * sizeof(int64_t));
    for (int64_t i = lo; i < end; i++) {
        reg = ((reg << 1) | (u64)cbit(by, i)) & mask;
        int64_t st = i - k + 1;
        if (st >= lo && reg == want) {
            if (cnt < maxpos) pos[cnt] = st;
            if (bins) bins[(int64_t)((double)(st - lo) / (double)(hi - lo) * npx)]++;
            cnt++;
        }
    }
    return cnt;
}

// Run-length histograms of runs starting in [lo, hi): ones[L], zeros[L] for L = 1..63 (bin 63 = 63+).
void r30c_runs(const uint8_t *by, int64_t nbits, int64_t lo, int64_t hi, int64_t *ones, int64_t *zeros) {
    memset(ones, 0, 64 * sizeof(int64_t));
    memset(zeros, 0, 64 * sizeof(int64_t));
    int64_t i = lo;
    if (i > 0) while (i < hi && cbit(by, i) == cbit(by, i - 1)) i++;  // skip the run that started before lo
    while (i < hi) {
        int v = cbit(by, i);
        int64_t j = i + 1;
        while (j < nbits && cbit(by, j) == v) j++;
        int64_t L = j - i < 63 ? j - i : 63;
        (v ? ones : zeros)[L]++;
        i = j;
    }
}
