// Test and benchmark of the Highway atan2 on every target that is compiled in and supported by this CPU,
// plus the other implementations available here (extras.c) and glibc's scalar atan2.
//
// usage: hwy_test [N_RANDOM] [--no-bench] [--hash FILE]
//   --hash FILE  append "arch atan2 target svml n hash" lines (FNV-1a of the result bits; every FMA target
//                must give the same hash on every machine, equal to the C translation's, which is
//                written as "svml-c" when the CPU has AVX-512) for the CI cross-machine check.
// Exit status 1 if a check fails:
//   - FMA targets: max error <= kMaxUlp, no error > 1 ulp, results for zero/inf/NaN references as
//     expected; inputs outside the vector path's range give the scalar callout's result, or for zero
//     arguments (no NaN) the exact IEEE result (+-0, +-pi/2, +-pi rounded, = glibc's atan2)
//   - all FMA targets produce identical results
//   - with AVX-512 (also under Intel SDE): every FMA target matches the C translation svml_atan28_ha
//     (= Intel's code) bit for bit
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "hwy/targets.h"

extern "C" {
#include "../atan2_d.h"
#include "extras.h"
#include "hwymath/accuracy.h"
void hwy_atan2_array(const double *y, const double *x, double *r, size_t len);
int hwy_atan2_variant(void);
size_t hwy_atan2_lanes(void);
}
#include "inputs.h"

#if defined(__x86_64__)
static const char *kArch = "x86_64";
#elif defined(__aarch64__)
static const char *kArch = "aarch64";
#else
static const char *kArch = "other";
#endif

static const double kMaxUlp = 0.9;  // bound for the FMA targets (measured: see Performance.md)

static std::vector<double> g_y, g_x;
static std::vector<long double> g_ref;  // atan2l(y, x), computed once

struct Stats {
    acc_stats s = ACC_STATS_INIT;
    size_t worst = 0;
};

static Stats accuracy(const std::vector<double> &r) {
    Stats st;
    for (size_t i = 0; i < r.size(); i++) {
        const double before = st.s.max_err;
        acc_add(&st.s, g_y[i], r[i], g_ref[i]);
        if (st.s.max_err > before) st.worst = i;
    }
    return st;
}

static bool in_range(double v) { return std::fabs(v) >= 0x1p-1020 && std::fabs(v) < 0x1p993; }

// expected result for inputs the vector algorithm does not compute itself; false if (y, x) is ordinary
static bool special_expected(double y, double x, double *r) {
    if (in_range(x) && in_range(y)) return false;
    if ((x == 0.0 || y == 0.0) && !std::isnan(x) && !std::isnan(y)) *r = std::atan2(y, x);
    else svml_datan2_ha_cout_rare(&y, &x, r);
    return true;
}

static void print_row(const char *name, const char *alg, size_t lanes, const Stats &st, double ns) {
    const acc_stats &s = st.s;
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
    g_y.resize(n_random + INPUTS_EXTRA);
    g_x.resize(n_random + INPUTS_EXTRA);
    const size_t n = make_inputs(g_y.data(), g_x.data(), n_random);
    g_y.resize(n);
    g_x.resize(n);
    g_ref.resize(n);
    for (size_t i = 0; i < n; i++) g_ref[i] = atan2l((long double)g_y[i], (long double)g_x[i]);
    std::vector<double> r(n);
    FILE *hf = hash_file ? fopen(hash_file, "a") : nullptr;
    if (hf) {
        std::vector<double> both(g_y);
        both.insert(both.end(), g_x.begin(), g_x.end());
        fprintf(hf, "%s atan2 inputs - %zu %016llx\n", kArch, n, (unsigned long long)acc_hash(both.data(), both.size()));
    }

    // benchmark inputs: y, x uniform in [-2, 2]
    const size_t NB = 4096;
    std::vector<double> by(NB), bx(NB), br(NB);
    for (size_t i = 0; i < NB; i++) { by[i] = rng_uniform(-2.0, 2.0); bx[i] = rng_uniform(-2.0, 2.0); }
    auto bench = [&](acc_array2_fn fn) { return do_bench ? acc_bench2(fn, by.data(), bx.data(), br.data(), NB, 1000) : 0.0; };

    int failures = 0;
    bool have_avx512 = false;
    std::vector<double> ref_intel;
#if defined(__x86_64__)
    have_avx512 = __builtin_cpu_supports("avx512f");
    if (have_avx512) {  // the C translation (bit-identical to Intel's code)
        ref_intel.resize(n);
        svml_atan2_ha_array(g_y.data(), g_x.data(), ref_intel.data(), n);
        if (hf) fprintf(hf, "%s atan2 svml-c svml %zu %016llx\n", kArch, n, (unsigned long long)acc_hash(ref_intel.data(), n));
    }
#endif

    printf("%s, %zu input pairs, errors vs atan2l (long double with %d-bit mantissa)\n\n", kArch, n, LDBL_MANT_DIG);
    printf("%-24s %-10s %5s %9s %10s %8s %8s %10s\n", "implementation", "algorithm", "lanes", "max ulp", "not CR",
           ">1 ulp", "special", do_bench ? "ns/value" : "");
    bool have_hash = false;
    uint64_t first_hash = 0;
    for (int64_t target : hwy::SupportedAndGeneratedTargets()) {
        hwy::SetSupportedTargetsForTest(target);
        const int variant = hwy_atan2_variant();
        const char *tname = hwy::TargetName(target);
        hwy_atan2_array(g_y.data(), g_x.data(), r.data(), n);
        const Stats st = accuracy(r);
        char name[64];
        snprintf(name, sizeof name, "hwy %s", tname);
        print_row(name, variant ? "SVML" : "std::atan2", hwy_atan2_lanes(), st, bench(hwy_atan2_array));
        if (!variant) continue;
        printf("  worst: y = %a, x = %a\n", g_y[st.worst], g_x[st.worst]);

        const uint64_t h = acc_hash(r.data(), n);
        if (hf) fprintf(hf, "%s atan2 %s svml %zu %016llx\n", kArch, tname, n, (unsigned long long)h);
        if (st.s.max_err > kMaxUlp || st.s.over1 || st.s.mismatched_special) {
            printf("  FAIL: error bound exceeded\n");
            failures++;
        }
        size_t bad_special = 0;
        for (size_t i = 0; i < n; i++) {
            double e;
            if (special_expected(g_y[i], g_x[i], &e) && dbits(e) != dbits(r[i])) {
                if (bad_special++ < 5) printf("  special y=%a x=%a: %a, expected %a\n", g_y[i], g_x[i], r[i], e);
            }
        }
        if (bad_special) { printf("  FAIL: %zu special inputs with unexpected results\n", bad_special); failures++; }
        if (!have_hash) { first_hash = h; have_hash = true; }
        else if (h != first_hash) { printf("  FAIL: results differ from the first FMA target\n"); failures++; }
        if (have_avx512) {
            const size_t nd = count_diff(r, ref_intel);
            printf("  bit-identical to svml_atan28_ha (Intel's code): %s (%zu differences)\n", nd ? "NO" : "yes", nd);
            if (nd) failures++;
        }
    }
    hwy::SetSupportedTargetsForTest(0);

#if defined(__x86_64__)
    if (have_avx512) {
        print_row("C SVML AVX-512", "SVML", 8, accuracy(ref_intel), bench(svml_atan2_ha_array));
    }
#endif
    struct extra_impl ex[16];
    const size_t nex = extra_impls(ex, 16);
    for (size_t k = 0; k < nex; k++) {
        ex[k].fn(g_y.data(), g_x.data(), r.data(), n);
        print_row(ex[k].name, "", 0, accuracy(r), bench(ex[k].fn));
    }
    for (size_t i = 0; i < n; i++) r[i] = std::atan2(g_y[i], g_x[i]);
    print_row("glibc atan2 (scalar)", "", 1, accuracy(r), bench([](const double *a, const double *b, double *c, size_t m) {
                  for (size_t i = 0; i < m; i++) c[i] = std::atan2(a[i], b[i]);
              }));

    if (hf) fclose(hf);
    printf("\n%s\n", failures ? "FAILED" : "all checks passed");
    return failures ? 1 : 0;
}
