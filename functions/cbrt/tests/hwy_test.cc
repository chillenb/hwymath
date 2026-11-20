// Test and benchmark of the Highway cbrt on every target that is compiled in and supported by this CPU,
// plus the other implementations available here (extras*.c) and glibc's scalar cbrt.
//
// usage: hwy_test [N_RANDOM] [--no-bench] [--hash FILE]
//   --hash FILE  append "arch cbrt target svml n hash" lines (FNV-1a of the result bits; variant svml = the
//                SVML algorithm, bit-identical to Intel's __svml_cbrt8_ha), plus the C/AVX-512 translation's
//                hash ("svml-c-avx512") when the CPU has AVX-512F; the CI cross-check compares them.
// Exit status 1 if a check fails, for every target running the SVML algorithm:
//   - max error < 1 ulp (HA), special values (0, inf, NaN) exact
//   - all such targets produce identical results, and odd symmetry cbrt(-x) = -cbrt(x) holds bit for bit
//   - x86-64 with AVX-512F (e.g. under Intel SDE): identical to the C translation of Intel's code
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "hwy/targets.h"

extern "C" {
#include "../cbrt_d.h"
#include "extras.h"
#include "hwymath/accuracy.h"
void hwy_cbrt_array(const double *x, double *y, size_t len);
int hwy_cbrt_variant(void);
size_t hwy_cbrt_lanes(void);
}
#include "inputs.h"

#if defined(__x86_64__)
static const char *kArch = "x86_64";
#elif defined(__aarch64__)
static const char *kArch = "aarch64";
#else
static const char *kArch = "other";
#endif

static std::vector<long double> g_ref;  // cbrtl(x), computed once

static acc_stats accuracy(const std::vector<double> &x, const std::vector<double> &y) {
    acc_stats s = ACC_STATS_INIT;
    for (size_t i = 0; i < x.size(); i++) acc_add(&s, x[i], y[i], g_ref[i]);
    return s;
}

static double bench(void (*f)(const double *, double *, size_t)) {  // ns/value, 4096 values in [-1000, 1000]
    const size_t N = 4096;
    static double in[N], out[N];
    for (size_t i = 0; i < N; i++) in[i] = rng_uniform(-1000.0, 1000.0);
    return acc_bench(f, in, out, N, 2000);
}

static void print_row(const char *name, const char *alg, size_t lanes, const acc_stats &s, double ns) {
    char l[24] = "", t[24] = "";
    if (lanes) snprintf(l, sizeof l, "%zu", lanes);
    if (ns > 0) snprintf(t, sizeof t, "%.2f", ns);
    printf("%-24s %-10s %5s %9.4f %9.4f%% %8zu %8zu %10s\n", name, alg, l, s.max_err, s.n ? 100.0 * s.not_cr / s.n : 0.0,
           s.over1, s.mismatched_special, t);
}

static size_t count_diff(const std::vector<double> &a, const std::vector<double> &b) {
    size_t n = 0;
    for (size_t i = 0; i < a.size(); i++) n += dbits(a[i]) != dbits(b[i]);
    return n;
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
    std::vector<double> x(n_random + CBRT_FIXED_INPUTS);
    x.resize(make_inputs(x.data(), n_random));
    std::vector<double> xneg(x.size());
    for (size_t i = 0; i < x.size(); i++) xneg[i] = -x[i];
    g_ref.resize(x.size());
    for (size_t i = 0; i < x.size(); i++) g_ref[i] = cbrtl((long double)x[i]);
    std::vector<double> y(x.size()), yneg(x.size());
    FILE *hf = hash_file ? fopen(hash_file, "a") : nullptr;
    if (hf) fprintf(hf, "%s cbrt inputs - %zu %016llx\n", kArch, x.size(), (unsigned long long)acc_hash(x.data(), x.size()));

    std::vector<double> ref_c;  // the C translation of Intel's AVX-512 code (bit-identical to it)
#if defined(__x86_64__)
    if (__builtin_cpu_supports("avx512f")) {
        ref_c.resize(x.size());
        svml_cbrt_ha_array(x.data(), ref_c.data(), x.size());
        if (hf) fprintf(hf, "%s cbrt svml-c-avx512 svml %zu %016llx\n", kArch, x.size(),
                        (unsigned long long)acc_hash(ref_c.data(), x.size()));
    }
#endif

    printf("%s, %zu inputs, errors vs cbrtl (long double with %d-bit mantissa)%s\n\n", kArch, x.size(), LDBL_MANT_DIG,
           ref_c.empty() ? "" : "; AVX-512F present: comparing with the C translation");
    printf("%-24s %-10s %5s %9s %10s %8s %8s %10s\n", "implementation", "algorithm", "lanes", "max ulp", "not CR",
           ">1 ulp", "special", do_bench ? "ns/value" : "");
    int failures = 0;
    bool have_first = false;
    std::vector<double> first;  // results of the first target with the SVML algorithm
    for (int64_t target : hwy::SupportedAndGeneratedTargets()) {
        hwy::SetSupportedTargetsForTest(target);
        const int variant = hwy_cbrt_variant();
        const char *tname = hwy::TargetName(target);
        hwy_cbrt_array(x.data(), y.data(), x.size());
        const acc_stats s = accuracy(x, y);
        char name[64];
        snprintf(name, sizeof name, "hwy %s", tname);
        print_row(name, variant ? "svml" : "std::cbrt", hwy_cbrt_lanes(), s, do_bench ? bench(hwy_cbrt_array) : 0);
        if (!variant) continue;

        if (hf) fprintf(hf, "%s cbrt %s svml %zu %016llx\n", kArch, tname, y.size(), (unsigned long long)acc_hash(y.data(), y.size()));
        if (!(s.max_err < 1.0) || s.over1) { printf("  FAIL: error bound exceeded (x = %a)\n", s.worst_x); failures++; }
        if (s.mismatched_special) { printf("  FAIL: %zu special values wrong\n", s.mismatched_special); failures++; }
        hwy_cbrt_array(xneg.data(), yneg.data(), x.size());
        size_t asym = 0;
        for (size_t i = 0; i < x.size(); i++) asym += dbits(yneg[i]) != (dbits(y[i]) ^ 0x8000000000000000ULL);
        if (asym) { printf("  FAIL: cbrt(-x) != -cbrt(x) for %zu inputs\n", asym); failures++; }
        if (!have_first) {  // informational: exactness at perfect cubes and powers of 8 (all exponents)
            size_t ncube = 0, cube_inexact = 0, npow = 0, pow_inexact = 0;
            for (size_t i = 0; i < x.size(); i++) {
                if (!std::isfinite(x[i]) || x[i] == 0) continue;
                // perfect cubes: the rounded root r has <= 21 significant bits, so r^3 is exact in long double
                const double r = std::fabs((double)g_ref[i]);
                if ((dbits(r) & 0xffffffffULL) != 0 || (long double)r * r * r != fabsl((long double)x[i])) continue;
                const bool pow8 = (dbits(r) & 0x000fffffffffffffULL) == 0;  // r a power of two
                (pow8 ? npow : ncube)++;
                if (dbits(std::fabs(y[i])) != dbits(r)) (pow8 ? pow_inexact : cube_inexact)++;
            }
            printf("  exact cube roots: %zu of %zu powers of 8, %zu of %zu other perfect cubes\n", npow - pow_inexact,
                   npow, ncube - cube_inexact, ncube);
            first = y;
            have_first = true;
        } else if (const size_t nd = count_diff(first, y)) { printf("  FAIL: %zu results differ from the first SVML target\n", nd); failures++; }
        if (!ref_c.empty()) {
            const size_t nd = count_diff(ref_c, y);
            printf("  %s the C translation of Intel's AVX-512 code: %zu of %zu results differ\n", nd ? "FAIL: vs" : "bit-identical to", nd, y.size());
            failures += nd != 0;
        }
    }
    hwy::SetSupportedTargetsForTest(0);

    struct extra_impl ex[16];
    const size_t nex = extra_impls(ex, 16);
    for (size_t k = 0; k < nex; k++) {
        ex[k].fn(x.data(), y.data(), x.size());
        print_row(ex[k].name, "", 0, accuracy(x, y), do_bench ? bench(ex[k].fn) : 0);
    }
    for (size_t i = 0; i < x.size(); i++) y[i] = std::cbrt(x[i]);
    print_row("glibc cbrt (scalar)", "", 1, accuracy(x, y),
              do_bench ? bench([](const double *a, double *b, size_t n) { for (size_t i = 0; i < n; i++) b[i] = std::cbrt(a[i]); }) : 0);

    if (hf) fclose(hf);
    printf("\n%s\n", failures ? "FAILED" : "all checks passed");
    return failures ? 1 : 0;
}
