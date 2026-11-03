// Q01 escape bound by bit-tree search. For each odd Y of exact bit length W, count the G-steps
// until bit 0 is clear (G = F∘F, F(Y) = (Y>>1) ^ (Y | (Y<<1))). Bit 0 after s steps depends only
// on bits 0..2s of Y, so bits are fixed from the bottom and a branch is dropped once it has
// reached an even number. Y lives in 256 bits; bits above 2(S - j) are masked away at step j,
// where S is the step cap, so nothing can overflow.
//   tree W S   -> prints W, the record step count, its holders, the first holder, and the
//                number of Y still odd after S steps (a counterexample to the bound S).
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define NL 4
typedef struct { uint64_t w[NL]; } bits;
static inline void shr1(const bits *a, bits *r) { for (int i = 0; i < NL; i++) r->w[i] = (a->w[i] >> 1) | (i + 1 < NL ? a->w[i + 1] << 63 : 0); }
static inline void shl1(const bits *a, bits *r) { for (int i = NL - 1; i >= 0; i--) r->w[i] = (a->w[i] << 1) | (i ? a->w[i - 1] >> 63 : 0); }
static inline void F(bits *y) {
    bits l, r; shr1(y, &l); shl1(y, &r);
    for (int i = 0; i < NL; i++) y->w[i] = l.w[i] ^ (y->w[i] | r.w[i]);
}
static inline void maskto(bits *y, int nb) {   // keep bits 0..nb-1
    for (int i = 0; i < NL; i++) {
        int lo = 64 * i;
        if (nb <= lo) y->w[i] = 0;
        else if (nb < lo + 64) y->w[i] &= (1ULL << (nb - lo)) - 1;
    }
}
static int W, S, rec; static uint64_t holders, longc, leaves; static bits first;
// Steps survived (bit 0 stays 1) by the known low bits 0..i, checking only steps <= i/2.
static int alive(const bits *y, int i) {
    bits z = *y; int s = i / 2;
    for (int j = 1; j <= s; j++) { F(&z); F(&z); if (!(z.w[0] & 1)) return 0; }
    return 1;
}
static void leaf(const bits *y) {
    leaves++;
    bits z = *y; int s = 0;
    while (z.w[0] & 1) {
        if (s == S) { longc++; return; }
        F(&z); F(&z); s++;
        maskto(&z, 2 * (S - s) + 1);
    }
    if (s > rec) { rec = s; holders = 1; first = *y; } else if (s == rec) holders++;
}
static void dfs(bits *y, int i) {   // bits 0..i-1 fixed
    if (i == W) { leaf(y); return; }
    for (int b = 0; b < 2; b++) {
        if (i == W - 1 && b == 0) continue;    // exact bit length W
        if (i == 0 && b == 0) continue;        // odd
        if (b) y->w[i >> 6] |= 1ULL << (i & 63);
        if (alive(y, i)) dfs(y, i + 1);
        y->w[i >> 6] &= ~(1ULL << (i & 63));
    }
}
int main(int argc, char **argv) {
    W = atoi(argv[1]); S = atoi(argv[2]);
    if (2 * S + 1 > 64 * NL || W > 2 * S) { fprintf(stderr, "need 2S+1 <= %d and W <= 2S\n", 64 * NL); return 2; }
    bits y; memset(&y, 0, sizeof y);
    dfs(&y, 0);
    printf("W=%d rec=%d holders=%llu long=%llu leaves=%llu first=", W, rec, (unsigned long long)holders,
           (unsigned long long)longc, (unsigned long long)leaves);
    for (int i = NL - 1, z = 1; i >= 0; i--) { if (z && !first.w[i] && i) continue; printf(z ? "%llx" : "%016llx", (unsigned long long)first.w[i]); z = 0; }
    printf("\n");
}
