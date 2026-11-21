// Checks of the VRCP14PD/VRSQRT14PD/VGETEXPPD emulations (hwy/avx512_emul-inl.h), which the Highway asinh,
// acosh and atanh rely on for bit-identity with Intel's AVX-512 code:
//   - the scalar reference (approx14_ref.h, after Intel's RECIP14.c) against the instructions themselves,
//     if the CPU (or Intel SDE) has AVX-512F;
//   - every Highway target with FMA against the scalar reference.
// The results depend on x only through its exponent and the top 16 (rcp14) or 15 (rsqrt14) mantissa
// bits, so all of these prefixes are tested (with low bits 0, all ones and random) at several exponents.
// usage: test_approx14 [N_EXPONENTS]   exit status 1 on a mismatch
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "hwy/targets.h"
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#include "approx14_ref.h"

extern "C" {
void hwy_rcp14_array(const double *x, double *y, size_t n);
void hwy_rsqrt14_array(const double *x, double *y, size_t n);
void hwy_getexp_array(const double *x, double *y, size_t n);
int hwy_approx14_emulated(void);
}

#if defined(__x86_64__)
__attribute__((target("avx512f"))) static void hw(const std::vector<double> &x, std::vector<double> &r1,
                                                  std::vector<double> &r2, std::vector<double> &r3) {
    for (size_t i = 0; i + 8 <= x.size(); i += 8) {
        const __m512d v = _mm512_loadu_pd(&x[i]);
        _mm512_storeu_pd(&r1[i], _mm512_rcp14_pd(v));
        _mm512_storeu_pd(&r2[i], _mm512_rsqrt14_pd(v));
        _mm512_storeu_pd(&r3[i], _mm512_getexp_pd(v));
    }
}
#endif

static size_t compare(const char *what, const std::vector<double> &x, const std::vector<double> &got,
                      double (*ref)(double), bool (*in_domain)(double)) {
    size_t bad = 0, n = 0;
    for (size_t i = 0; i < x.size(); i++) {
        if (!in_domain(x[i])) continue;
        n++;
        const double r = ref(x[i]);
        if (dbits(got[i]) != dbits(r) && bad++ < 5) printf("  MISMATCH %s(%a): %a, reference %a\n", what, x[i], got[i], r);
    }
    printf("  %-8s %zu inputs, %zu mismatches\n", what, n, bad);
    return bad;
}

static double ref_getexp(double x) { return (double)((int)((dbits(x) >> 52) & 0x7ff) - 1023); }
static bool rcp_domain(double x) { return std::isnormal(x) && x > 0 && x < 0x1p1022; }
static bool rsqrt_domain(double x) { return (std::isnormal(x) && x > 0) || x == 0 || x == INFINITY || std::isnan(x); }
static bool getexp_domain(double x) { return std::isnormal(x) && x > 0; }
static bool any(double x) { return rcp_domain(x) || rsqrt_domain(x); }

int main(int argc, char **argv) {
    const int n_exp = argc > 1 ? atoi(argv[1]) : 6;
    std::vector<double> x;
    for (const double s : {0.0, -0.0, (double)INFINITY, (double)NAN, 1.0, 2.0, 4.0, 0x1p-1022, 0x1p1021, 0x1.fffffffffffffp1021})
        x.push_back(s);
    std::vector<int> exps = {0, 1, -1, 2, -53, 1021, -1022, -1021};  // [1, 2), [2, 4), ... and the extremes
    while ((int)exps.size() < n_exp + 8) exps.push_back((int)(rng() % 2043) - 1021);
    for (const int e : exps) {
        for (uint64_t p = 0; p < (1u << 16); p++) {
            for (int v = 0; v < 3; v++) {
                const uint64_t low = v == 0 ? 0 : v == 1 ? 0xfffffffffULL : rng() & 0xfffffffffULL;
                x.push_back(bitsd(((uint64_t)(e + 1023) << 52) | (p << 36) | low));
            }
        }
    }
    while (x.size() % 8) x.push_back(1.0);
    std::vector<double> r1(x.size()), r2(x.size()), r3(x.size());
    size_t bad = 0;
#if defined(__x86_64__)
    if (__builtin_cpu_supports("avx512f")) {
        printf("reference vs AVX-512 instructions:\n");
        hw(x, r1, r2, r3);
        bad += compare("rcp14", x, r1, ref_rcp14, rcp_domain);
        bad += compare("rsqrt14", x, r2, ref_rsqrt14, rsqrt_domain);
        bad += compare("getexp", x, r3, ref_getexp, getexp_domain);
    } else {
        printf("reference vs AVX-512 instructions: SKIPPED (no AVX-512F; run under Intel SDE)\n");
    }
#endif
    (void)any;
    for (int64_t target : hwy::SupportedAndGeneratedTargets()) {
        hwy::SetSupportedTargetsForTest(target);
        if (!hwy_approx14_emulated()) {
            printf("%s: no FMA, not used\n", hwy::TargetName(target));
            continue;
        }
        printf("%s emulation vs reference:\n", hwy::TargetName(target));
        hwy_rcp14_array(x.data(), r1.data(), x.size());
        hwy_rsqrt14_array(x.data(), r2.data(), x.size());
        hwy_getexp_array(x.data(), r3.data(), x.size());
        bad += compare("rcp14", x, r1, ref_rcp14, rcp_domain);
        bad += compare("rsqrt14", x, r2, ref_rsqrt14, rsqrt_domain);
        bad += compare("getexp", x, r3, ref_getexp, getexp_domain);
    }
    hwy::SetSupportedTargetsForTest(0);
    printf("test_approx14: %s\n", bad ? "FAILED" : "all checks passed");
    return bad ? 1 : 0;
}
