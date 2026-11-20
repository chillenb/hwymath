/* Verification of the rounded-reciprocal model of SVML's log family (hwy/rcp14_thresholds.h) against the
 * CPU's vrcp14pd + vrndscalepd (AVX-512F; run under Intel SDE on other machines). Compiled as C++ (the
 * header is C++).
 *
 * k(m) = 32 * vrndscalepd(vrcp14pd(m), imm 0x58) for m in [1, 2). Checks, with W, G, Z from the command line:
 *   0. bisection (assuming monotonicity) finds exactly the 16 thresholds of the header;
 *   1. every double within 2^W ulps of each threshold;
 *   2. both ends of every interval of 2^(52-G) doubles (mantissa prefixes of G bits) within 2^-Z of each
 *      threshold (|vrcp14pd error| < 2^-14 relative, so k is only uncertain within ~2^-13 of a threshold);
 *   3. both ends of all 2^24 intervals of 2^28 doubles covering [1, 2);
 *   4. 2^24 random mantissas;
 * have k equal to the step function of the thresholds, and that the Highway emulation (the same integer
 * operations in scalar C) also equals it at all those points. First (also without AVX-512): the emulation
 * equals the step function for every high word of m (it depends on nothing else), i.e. for every m in [1, 2).
 * usage: verify_rcp14 [W [G [Z]]]   defaults 20 32 13; the README quotes a run with 22 36 13.
 * Exit status: 0 = all equal, 1 = mismatch, 77 = vrcp14pd checks skipped (no AVX-512F). */
#include <immintrin.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "functions/log/hwy/rcp14_thresholds.h"

using namespace hwymath_log;

static double bd(uint64_t u) { double d; memcpy(&d, &u, 8); return d; }
static uint64_t db(double d) { uint64_t u; memcpy(&u, &d, 8); return u; }
static const uint64_t ONE = 0x3ff0000000000000ULL, TWO = 0x4000000000000000ULL;

/* k(m) from the CPU for 8 mantissas */
__attribute__((target("avx512f"))) static void cpu_k8(const uint64_t *mb, int *k) {
    __m512d v = _mm512_rcp14_pd(_mm512_loadu_pd((const double *)mb));
    double r[8];
    _mm512_storeu_pd(r, _mm512_roundscale_round_pd(v, 0x58, _MM_FROUND_NO_EXC));
    for (int l = 0; l < 8; l++) k[l] = (int)(r[l] * 32);
}
static int cpu_k(uint64_t mb) { uint64_t a[8]; int k[8]; for (int l = 0; l < 8; l++) a[l] = mb; cpu_k8(a, k); return k[0]; }

/* the step function of the header */
static int model_k(double m) { int k = 32; for (int j = 15; j >= 0; j--) if (m >= kRcpThreshold[j]) k = 16 + j; return k; }

/* the Highway emulation (log_hwy-inl.h: RoundedRcp14), in scalar C: integer operations on the high word */
static int emul_k(double m) {
    const uint64_t mb = db(m);
    const uint32_t h = (uint32_t)(mb >> 32), e = kRcpBins.v[(mb >> 47) & 31];
    return 16 + (int)(e & 15) + ((h | 15) < e);
}

static long g_bad, g_n;
static void check8(const uint64_t *a, const char *what) {
    int k[8];
    cpu_k8(a, k);
    for (int l = 0; l < 8; l++) {
        const double m = bd(a[l]);
        const int want = model_k(m), e = emul_k(m);
        if ((k[l] != want || e != want) && g_bad++ < 10) printf("%s MISMATCH m=%a: cpu %d, emulation %d, model %d\n", what, m, k[l], e, want);
    }
    g_n += 8;
}

int main(int argc, char **argv) {
    /* emulation only: it depends on the high 32 bits of m only, so this covers every m in [1, 2) */
    long n5 = 0;
    for (uint64_t h = 0x3ff00000u; h < 0x40000000u; h++) {
        for (int lo = 0; lo < 2; lo++) {
            const double m = bd(h << 32 | (lo ? 0xffffffffu : 0u));
            if (emul_k(m) != model_k(m) && g_bad++ < 10) printf("EMULATION MISMATCH m=%a: %d vs %d\n", m, emul_k(m), model_k(m));
            n5++;
        }
    }
    printf("emulation for all %ld high words of m (low word 0 and all ones)\n", n5 / 2);
    if (!__builtin_cpu_supports("avx512f")) {
        printf("verify_rcp14: emulation: %ld mismatches; checks against vrcp14pd SKIPPED (this CPU has no AVX-512F)\n", g_bad);
        return g_bad ? 1 : 77;
    }
    const int W = argc > 1 ? atoi(argv[1]) : 20, G = argc > 2 ? atoi(argv[2]) : 32, Z = argc > 3 ? atoi(argv[3]) : 13;
    uint64_t a[8];
    /* 0. bisection */
    for (int j = 31; j >= 16; j--) {  /* smallest m with k(m) <= j */
        uint64_t lo = ONE, hi = TWO - 1;
        while (hi - lo > 1) { const uint64_t mid = lo + (hi - lo) / 2; if (cpu_k(mid) <= j) hi = mid; else lo = mid; }
        if (hi != db(kRcpThreshold[j - 16]) && g_bad++ < 10) printf("THRESHOLD MISMATCH k=%d: cpu %a, header %a\n", j, bd(hi), kRcpThreshold[j - 16]);
    }
    printf("0. bisection: the 16 thresholds\n");
    /* 1. windows */
    long n0 = g_n;
    for (int j = 0; j < 16; j++) {
        const uint64_t t = db(kRcpThreshold[j]);
        for (uint64_t b = t - (1ULL << W); b < t + (1ULL << W); b += 8) { for (int l = 0; l < 8; l++) a[l] = b + l; check8(a, "window"); }
    }
    printf("1. %ld doubles within 2^%d ulps of the thresholds\n", g_n - n0, W);
    /* 2. intervals of G-bit prefixes near the thresholds */
    n0 = g_n;
    const uint64_t sh = 52 - G, w = 1ULL << sh, half = 1ULL << (52 - Z);
    for (int j = 0; j < 16; j++) {
        const uint64_t t = db(kRcpThreshold[j]);
        uint64_t start = ((t - half) >> sh) << sh, end = t + half;
        if (start < ONE) start = ONE;
        if (end > TWO) end = TWO;
        for (uint64_t b = start; b < end; b += 4 * w) {
            for (int l = 0; l < 4; l++) { a[2 * l] = b + l * w; a[2 * l + 1] = b + l * w + w - 1; }
            check8(a, "interval");
        }
    }
    printf("2. both ends of %ld intervals of 2^%d doubles within 2^-%d of the thresholds\n", (g_n - n0) / 2, 52 - G, Z);
    /* 3. coarse cover of [1, 2) */
    n0 = g_n;
    const uint64_t w3 = 1ULL << 28;
    for (uint64_t b = ONE; b < TWO; b += 4 * w3) {
        for (int l = 0; l < 4; l++) { a[2 * l] = b + l * w3; a[2 * l + 1] = b + l * w3 + w3 - 1; }
        check8(a, "cover");
    }
    printf("3. both ends of all %ld intervals of 2^28 doubles in [1, 2)\n", (g_n - n0) / 2);
    /* 4. random */
    n0 = g_n;
    uint64_t s = 0x9e3779b97f4a7c15ULL;
    for (int i = 0; i < (1 << 24) / 8; i++) {
        for (int l = 0; l < 8; l++) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; a[l] = ONE | (s >> 12); }
        check8(a, "random");
    }
    printf("4. %ld random mantissas\n", g_n - n0);
    printf("verify_rcp14: %ld mismatches\n", g_bad);
    return g_bad != 0;
}
