// Test and benchmark of the Highway cosh on every target that is compiled in and supported by this CPU,
// plus Intel's libsvml and SLEEF when available (x86-64).
//
// usage: hwy_test [N_RANDOM] [--no-bench] [--hash FILE]
//   --hash FILE  append "arch cosh target variant n hash" lines (FNV-1a of the result bits; variant svml =
//                the SVML algorithm, whose results must be bit-identical on every machine and target), and
//                "arch cosh svml-c svml n hash" for the C/AVX-512 translation when the CPU has AVX-512.
// Exit status 1 if a check fails:
//   - SVML-algorithm targets: max error < 1 ulp (HA), special values exact;
//   - all SVML-algorithm targets produce identical results, also with vectors capped at 2 and 4 lanes;
//   - x86-64 with AVX-512F/DQ (e.g. under Intel SDE): every SVML-algorithm target is bit-identical to the
//     C translation of Intel's AVX-512 code (svml_cosh_ha_array).
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "hwy/targets.h"

extern "C" {
#include "../cosh_d.h"
void hwy_cosh_array(const double *x, double *y, size_t len);
void hwy_cosh_array_capped(const double *x, double *y, size_t len, size_t lanes);
int hwy_cosh_variant(void);
size_t hwy_cosh_lanes(void);
}
#include "hwymath/accuracy.h"
#include "inputs.h"

#define FNAME "cosh"
static long double reference(long double x) { return coshl(x); }
// exact expectations: {x, cosh(x)}
static const double kSpecial[][2] = {
    {0.0, 1.0}, {-0.0, 1.0}, {INFINITY, INFINITY}, {-INFINITY, INFINITY},
    {1e-300, 1.0}, {-1e-300, 1.0}, {5e-324, 1.0}, {-5e-324, 1.0}, {0x1p-30, 1.0}, {-0x1.8p-26, 0x1.0000000000001p+0},
    {710.5, INFINITY}, {-710.5, INFINITY}, {1e300, INFINITY}, {-1e300, INFINITY},
};

#if defined(__x86_64__)
static const char *kArch = "x86_64";
#include <immintrin.h>
#define AVX2 __attribute__((target("avx2,fma")))
#define AVX512 __attribute__((target("avx512f")))
extern "C" {
#ifdef HAVE_INTEL_SVML
AVX2 __m256d __svml_cosh4_ha_l9(__m256d);
AVX512 __m512d __svml_cosh8_ha_z0(__m512d);
#endif
#ifdef HAVE_SLEEF  // declared here: sleef.h only declares them when compiled with -mavx2 / -mavx512f
AVX2 __m256d Sleef_coshd4_u10avx2(__m256d);
AVX512 __m512d Sleef_coshd8_u10avx512f(__m512d);
#endif
}
#define LOOP(NAME, ATTR, VT, W, LOAD, STORE, CALL)                       \
    ATTR static void NAME(const double *x, double *y, size_t n) {        \
        size_t i = 0;                                                    \
        for (; i + W <= n; i += W) STORE(y + i, CALL(LOAD(x + i)));      \
        if (i < n) {                                                     \
            double b[W] = {0}, o[W];                                     \
            memcpy(b, x + i, (n - i) * 8);                               \
            STORE(o, CALL(LOAD(b)));                                     \
            memcpy(y + i, o, (n - i) * 8);                               \
        }                                                                \
    }
#ifdef HAVE_INTEL_SVML
LOOP(intel_avx2, AVX2, __m256d, 4, _mm256_loadu_pd, _mm256_storeu_pd, __svml_cosh4_ha_l9)
LOOP(intel_avx512, AVX512, __m512d, 8, _mm512_loadu_pd, _mm512_storeu_pd, __svml_cosh8_ha_z0)
#endif
#ifdef HAVE_SLEEF
LOOP(sleef_avx2, AVX2, __m256d, 4, _mm256_loadu_pd, _mm256_storeu_pd, Sleef_coshd4_u10avx2)
LOOP(sleef_avx512, AVX512, __m512d, 8, _mm512_loadu_pd, _mm512_storeu_pd, Sleef_coshd8_u10avx512f)
#endif
#elif defined(__aarch64__)
static const char *kArch = "aarch64";
#else
static const char *kArch = "other";
#endif

typedef void (*array_fn)(const double *, double *, size_t);
static void scalar_libm(const double *x, double *y, size_t n) {
    for (size_t i = 0; i < n; i++) y[i] = std::cosh(x[i]);
}

static std::vector<long double> g_ref;

static acc_stats accuracy(const std::vector<double> &x, const std::vector<double> &y) {
    acc_stats s = ACC_STATS_INIT;
    for (size_t i = 0; i < x.size(); i++) acc_add(&s, x[i], y[i], g_ref[i]);
    return s;
}

static double bench(array_fn f) {  // ns/value, 4096 values uniform in [-600, 600] (no callout lanes)
    static double in[4096], out[4096];
    for (size_t i = 0; i < 4096; i++) in[i] = rng_uniform(-600.0, 600.0);
    return acc_bench(f, in, out, 4096, 2000);
}

static void print_row(const char *name, const char *alg, size_t lanes, const acc_stats &s, double ns) {
    char l[24] = "", t[24] = "";
    if (lanes) snprintf(l, sizeof l, "%zu", lanes);
    if (ns > 0) snprintf(t, sizeof t, "%.2f", ns);
    printf("%-24s %-10s %5s %9.4f %9.4f%% %8zu %8zu %10s\n", name, alg, l, s.max_err, s.n ? 100.0 * s.not_cr / s.n : 0.0,
           s.over1, s.mismatched_special, t);
}

int main(int argc, char **argv) {
    size_t n_random = 20000000;
    bool do_bench = true;
    const char *hash_file = nullptr;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--no-bench")) do_bench = false;
        else if (!strcmp(argv[i], "--hash") && i + 1 < argc) hash_file = argv[++i];
        else n_random = strtoull(argv[i], nullptr, 10);
    }
    std::vector<double> x(n_random + 8192);
    x.resize(make_inputs(x.data(), n_random));
    g_ref.resize(x.size());
    for (size_t i = 0; i < x.size(); i++) g_ref[i] = reference((long double)x[i]);
    std::vector<double> y(x.size());
    FILE *hf = hash_file ? fopen(hash_file, "a") : nullptr;
    if (hash_file && !hf) { perror(hash_file); return 1; }
    if (hf) fprintf(hf, "%s %s inputs - %zu %016llx\n", kArch, FNAME, x.size(), (unsigned long long)acc_hash(x.data(), x.size()));

    printf("%s, %zu inputs, errors vs %sl (long double with %d-bit mantissa)\n\n", kArch, x.size(), FNAME, LDBL_MANT_DIG);
    printf("%-24s %-10s %5s %9s %10s %8s %8s %10s\n", "implementation", "algorithm", "lanes", "max ulp", "not CR",
           ">1 ulp", "special", do_bench ? "ns/value" : "");
    int failures = 0;

    std::vector<double> ref_c;  // C translation of Intel's AVX-512 code
#if defined(__x86_64__)
    if (__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq")) {
        ref_c.resize(x.size());
        svml_cosh_ha_array(x.data(), ref_c.data(), x.size());
        const acc_stats s = accuracy(x, ref_c);
        print_row("C translation AVX-512", "svml", 8, s, do_bench ? bench(svml_cosh_ha_array) : 0);
        if (hf) fprintf(hf, "%s %s svml-c svml %zu %016llx\n", kArch, FNAME, x.size(),
                        (unsigned long long)acc_hash(ref_c.data(), ref_c.size()));
    }
#endif

    bool have_first = false;
    uint64_t first_hash = 0;
    for (int64_t target : hwy::SupportedAndGeneratedTargets()) {
        hwy::SetSupportedTargetsForTest(target);
        const int variant = hwy_cosh_variant();
        const char *tname = hwy::TargetName(target);
        hwy_cosh_array(x.data(), y.data(), x.size());
        const acc_stats s = accuracy(x, y);
        char name[64];
        snprintf(name, sizeof name, "hwy %s", tname);
        print_row(name, variant ? "svml" : "std::" FNAME, hwy_cosh_lanes(), s, do_bench ? bench(hwy_cosh_array) : 0);

        // special values (every target)
        for (const auto &sv : kSpecial) {
            double r;
            hwy_cosh_array(&sv[0], &r, 1);
            if (!(r == sv[1] && std::signbit(r) == std::signbit(sv[1]))) {
                printf("  FAIL: %s(%a) = %a, expected %a\n", FNAME, sv[0], r, sv[1]);
                failures++;
            }
        }
        if (s.mismatched_special) { printf("  FAIL: NaN/inf/zero results differ from the reference\n"); failures++; }
        if (!variant) continue;

        const uint64_t h = acc_hash(y.data(), y.size());
        if (hf) fprintf(hf, "%s %s %s svml %zu %016llx\n", kArch, FNAME, tname, y.size(), (unsigned long long)h);
        if (!(s.max_err < 1.0)) { printf("  FAIL: error bound (1 ulp) exceeded\n"); failures++; }
        if (!have_first) { have_first = true; first_hash = h; }
        else if (h != first_hash) { printf("  FAIL: results differ from the first SVML-algorithm target\n"); failures++; }
        if (!ref_c.empty()) {
            size_t diff = 0;
            for (size_t i = 0; i < x.size(); i++) {
                if (dbits(y[i]) != dbits(ref_c[i]) && diff++ < 5) {
                    printf("  x = %a: hwy %a, C translation %a\n", x[i], y[i], ref_c[i]);
                }
            }
            if (diff) { printf("  FAIL: %zu results differ from the C translation\n", diff); failures++; }
            else printf("  bit-identical to the C translation of the AVX-512 code\n");
        }
        // the 2- and 4-lane code paths (lookups) with vectors capped at 2/4 lanes on this target
        for (size_t lanes : {size_t{2}, size_t{4}}) {
            if (lanes >= hwy_cosh_lanes()) continue;
            std::vector<double> yc(x.size());
            hwy_cosh_array_capped(x.data(), yc.data(), x.size(), lanes);
            if (acc_hash(yc.data(), yc.size()) != h) {
                printf("  FAIL: results with vectors capped at %zu lanes differ\n", lanes);
                failures++;
            }
        }
    }
    hwy::SetSupportedTargetsForTest(0);

    struct { const char *name; array_fn fn; size_t lanes; bool avx512; } extras[] = {
#if defined(__x86_64__) && defined(HAVE_INTEL_SVML)
        {"Intel libsvml HA AVX2", intel_avx2, 4, false},
        {"Intel libsvml HA AVX-512", intel_avx512, 8, true},
#endif
#if defined(__x86_64__) && defined(HAVE_SLEEF)
        {"SLEEF u10 AVX2", sleef_avx2, 4, false},
        {"SLEEF u10 AVX-512", sleef_avx512, 8, true},
#endif
        {"glibc " FNAME " (scalar)", scalar_libm, 1, false},
    };
    for (const auto &e : extras) {
#if defined(__x86_64__)
        if (e.avx512 && !__builtin_cpu_supports("avx512f")) continue;
        if (e.lanes == 4 && !(__builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma"))) continue;
#endif
        e.fn(x.data(), y.data(), x.size());
        print_row(e.name, "", e.lanes, accuracy(x, y), do_bench ? bench(e.fn) : 0);
    }

    if (hf) fclose(hf);
    printf("\n%s\n", failures ? "FAILED" : "all checks passed");
    return failures ? 1 : 0;
}
