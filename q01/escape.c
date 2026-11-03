// Q01 integer form: F(Y) = (Y>>1) ^ (Y | (Y<<1)), G = F∘F. For every odd Y < 2^K, count the
// G-steps until bit 0 is clear. Prints, per bit-length W, the record step count and how many
// Y attain it, plus the first record holder.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static inline u128 F(u128 y) { return (y >> 1) ^ (y | (y << 1)); }
int main(int argc, char **argv) {
    int K = atoi(argv[1]); uint64_t part = atoll(argv[2]), parts = atoll(argv[3]);
    int rec[129] = {0}; uint64_t cnt[129] = {0}, who[129] = {0};
    for (uint64_t y0 = 2 * part + 1; y0 < (1ULL << K); y0 += 2 * parts) {
        u128 y = y0; int s = 0;
        while (y & 1) { y = F(F(y)); s++; if (s > 60) { printf("LONG %llu\n", (unsigned long long)y0); break; } }
        int W = 64 - __builtin_clzll(y0);
        if (s > rec[W]) { rec[W] = s; cnt[W] = 1; who[W] = y0; } else if (s == rec[W]) cnt[W]++;
    }
    for (int W = 1; W <= K; W++)
        printf("%d %d %llu %llx\n", W, rec[W], (unsigned long long)cnt[W], (unsigned long long)who[W]);
}
