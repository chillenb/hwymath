// Test and benchmark of the Highway acos on every target that is compiled in and supported by this CPU,
// plus the other implementations available here (extras.c) and glibc's scalar acos.
//
// usage: hwy_test [N_RANDOM] [--no-bench] [--hash FILE]
//   --hash FILE  append "arch acos target svml n hash" lines (FNV-1a of the result bits; every FMA target
//                must give the same hash on every machine, equal to the C translation's, which is
//                written as "svml-c" when the CPU has AVX-512) for the CI cross-machine check.
// Exit status 1 if a check fails:
//   - FMA targets: max error <= kMaxUlp, no error > 1 ulp, NaN/inf/zero/|x| > 1 results as expected
//     (equal to the scalar callout, NaN = quieted input with flipped sign, acos(1) = +0, acos(-1) = RN(pi))
//   - all FMA targets produce identical results
//   - with AVX-512 (also under Intel SDE): every FMA target matches the C translation svml_acos8_ha
//     (= Intel's code) bit for bit
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "hwy/targets.h"

extern "C" {
#include "../acos_d.h"
#include "extras.h"
#include "hwymath/accuracy.h"
void hwy_acos_array(const double *x, double *y, size_t len);
int hwy_acos_variant(void);
size_t hwy_acos_lanes(void);
}
#include "inputs.h"

#if defined(__x86_64__)
static const char *kArch = "x86_64";
#elif defined(__aarch64__)
static const char *kArch = "aarch64";
#else
static const char *kArch = "other";
#endif

static const double kMaxUlp = 0.9;   // documented bound of the FMA targets (measured: see Performance.md)

static std::vector<long double> g_ref;  // acosl(x), computed once

static acc_stats accuracy(const std::vector<double> &x, const std::vector<double> &y) {
    acc_stats s = ACC_STATS_INIT;
    for (size_t i = 0; i < x.size(); i++) acc_add(&s, x[i], y[i], g_ref[i]);
    return s;
}

// expected result for inputs the vector algorithm does not compute itself, or 0 if x is ordinary
static bool special_expected(double x, double *r) {
    if (std::isnan(x)) { *r = -(x + x); return true; }  // quieted input with the sign flipped (SVML)
    if (std::fabs(x) > 1.0) { svml_dacos_ha_cout_rare(&x, r); return true; }
    if (std::fabs(x) == 1.0) { *r = x > 0 ? 0.0 : (double)acosl(-1.0L); return true; }
    return false;
}

static void print_row(const char *name, const char *alg, size_t lanes, const acc_stats &s, double ns) {
    char l[24] = "", t[24] = "";
    if (lanes) snprintf(l, sizeof l, "%zu", lanes);
    if (ns > 0) snprintf(t, sizeof t, "%.2f", ns);
    printf("%-24s %-10s %5s %9.4f %9.4f%% %8zu %8zu %10s\n", name, alg, l, s.max_err, s.n ? 100.0 * s.not_cr / s.n : 0.0,
           s.over1, s.mismatched_special, t);
}

// number of results that differ bit-wise; prints the first few
static size_t count_diff(const std::vector<double> &a, const std::vector<double> &b, const std::vector<double> *x = nullptr) {
    size_t n = 0;
    for (size_t i = 0; i < a.size(); i++) {
        if (dbits(a[i]) != dbits(b[i]) && n++ < 3 && x) printf("    x = %a: %a vs %a\n", (*x)[i], a[i], b[i]);
    }
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
    std::vector<double> x(n_random + 8192);
    x.resize(make_inputs(x.data(), n_random));
    g_ref.resize(x.size());
    for (size_t i = 0; i < x.size(); i++) g_ref[i] = acosl((long double)x[i]);
    std::vector<double> y(x.size());
    FILE *hf = hash_file ? fopen(hash_file, "a") : nullptr;
    if (hf) fprintf(hf, "%s acos inputs - %zu %016llx\n", kArch, x.size(), (unsigned long long)acc_hash(x.data(), x.size()));

    // benchmark inputs: uniform in [-1, 1]
    const size_t NB = 4096;
    std::vector<double> bx(NB), by(NB);
    for (size_t i = 0; i < NB; i++) bx[i] = rng_uniform(-1.0, 1.0);
    auto bench = [&](acc_array_fn fn) { return do_bench ? acc_bench(fn, bx.data(), by.data(), NB, 2000) : 0.0; };

    int failures = 0;
    bool have_avx512 = false;
    std::vector<double> ref_intel;
#if defined(__x86_64__)
    have_avx512 = __builtin_cpu_supports("avx512f");
    if (have_avx512) {  // the C translation (bit-identical to Intel's code)
        ref_intel.resize(x.size());
        svml_acos_ha_array(x.data(), ref_intel.data(), x.size());
        if (hf) fprintf(hf, "%s acos svml-c svml %zu %016llx\n", kArch, x.size(), (unsigned long long)acc_hash(ref_intel.data(), x.size()));
    }
#endif

    printf("%s, %zu inputs, errors vs acosl (long double with %d-bit mantissa)\n\n", kArch, x.size(), LDBL_MANT_DIG);
    printf("%-24s %-10s %5s %9s %10s %8s %8s %10s\n", "implementation", "algorithm", "lanes", "max ulp", "not CR",
           ">1 ulp", "special", do_bench ? "ns/value" : "");
    bool have_hash = false;
    uint64_t first_hash = 0;
    for (int64_t target : hwy::SupportedAndGeneratedTargets()) {
        hwy::SetSupportedTargetsForTest(target);
        const int variant = hwy_acos_variant();
        const char *tname = hwy::TargetName(target);
        hwy_acos_array(x.data(), y.data(), x.size());
        const acc_stats s = accuracy(x, y);
        char name[64];
        snprintf(name, sizeof name, "hwy %s", tname);
        print_row(name, variant ? "SVML" : "std::acos", hwy_acos_lanes(), s, bench(hwy_acos_array));
        if (!variant) continue;

        const uint64_t h = acc_hash(y.data(), y.size());
        if (hf) fprintf(hf, "%s acos %s svml %zu %016llx\n", kArch, tname, y.size(), (unsigned long long)h);
        if (s.max_err > kMaxUlp || s.over1 || s.mismatched_special) { printf("  FAIL: error bound exceeded\n"); failures++; }
        size_t bad_special = 0;
        for (size_t i = 0; i < x.size(); i++) {
            double r;
            if (special_expected(x[i], &r) && dbits(r) != dbits(y[i])) bad_special++;
        }
        if (bad_special) { printf("  FAIL: %zu special inputs with unexpected results\n", bad_special); failures++; }
        if (!have_hash) { first_hash = h; have_hash = true; }
        else if (h != first_hash) { printf("  FAIL: results differ from the first FMA target\n"); failures++; }
        if (have_avx512) {
            const size_t nd = count_diff(y, ref_intel, &x);
            printf("  bit-identical to svml_acos8_ha (Intel's code): %s (%zu differences)\n", nd ? "NO" : "yes", nd);
            if (nd) failures++;
        }
    }
    hwy::SetSupportedTargetsForTest(0);

#if defined(__x86_64__)
    if (have_avx512) {
        print_row("C SVML AVX-512", "SVML", 8, accuracy(x, ref_intel), bench(svml_acos_ha_array));
    }
#endif
    struct extra_impl ex[16];
    const size_t nex = extra_impls(ex, 16);
    for (size_t k = 0; k < nex; k++) {
        ex[k].fn(x.data(), y.data(), x.size());
        print_row(ex[k].name, "", 0, accuracy(x, y), bench(ex[k].fn));
    }
    for (size_t i = 0; i < x.size(); i++) y[i] = std::acos(x[i]);
    print_row("glibc acos (scalar)", "", 1, accuracy(x, y),
              bench([](const double *a, double *b, size_t n) { for (size_t i = 0; i < n; i++) b[i] = std::acos(a[i]); }));

    if (hf) fclose(hf);
    printf("\n%s\n", failures ? "FAILED" : "all checks passed");
    return failures ? 1 : 0;
}
