// Test and benchmark of a Highway trig function (sin, cos, tan) on every target that is compiled in and
// supported by this CPU, plus the other implementations available here (extras_tmpl.h). Included by
// functions/{sin,cos,tan}/tests/hwy_test.cc, which define
//   TRIG_NAME "sin", TRIG_HWY_ARRAY hwy_sin_array, TRIG_HWY_VARIANT hwy_sin_variant,
//   TRIG_HWY_LANES hwy_sin_lanes, TRIG_REFL sinl, TRIG_STD std::sin, TRIG_THRESHOLD 0x1p24 (start of the
//   large-argument path), and on x86-64 TRIG_SVML_ARRAY svml_sin_ha_array (the C/AVX-512 translation);
//   optionally TRIG_TARGET_CHECK(target_name), returning a number of failures (extra per-target checks).
//
// usage: hwy_test [N_RANDOM] [--no-bench] [--hash FILE]
//   --hash FILE  append "arch <f> <target> svml <n> <hash>" lines (FNV-1a of the result bits) for the
//                targets running the SVML algorithm, "arch <f> C-AVX512 svml <n> <hash>" for the C
//                translation when the CPU has AVX-512, and "arch <f> inputs - <n> <hash>": all svml hashes
//                must be equal (on every machine), which the CI cross-check verifies.
// Exit status 1 if a check fails:
//   - SVML-algorithm targets: max error < 1 ulp (the HA bound), NaN exactly where the reference is NaN,
//     signed zeros and infinities as in the reference
//   - all SVML-algorithm targets produce identical results
//   - with AVX-512F/DQ (or Intel SDE): every SVML-algorithm target matches the C translation (and hence
//     Intel's assembly, see test_<f>8) bit for bit
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "hwy/targets.h"

extern "C" {
#include "functions/sin/tests/extras_tmpl_decl.h"
#include "hwymath/accuracy.h"
void TRIG_HWY_ARRAY(const double *x, double *y, size_t len);
int TRIG_HWY_VARIANT(void);
size_t TRIG_HWY_LANES(void);
}
#include "functions/sin/tests/inputs.h"

#if defined(__x86_64__)
static const char *kArch = "x86_64";
#elif defined(__aarch64__)
static const char *kArch = "aarch64";
#else
static const char *kArch = "other";
#endif

#ifndef TRIG_TARGET_CHECK
#define TRIG_TARGET_CHECK(name) 0
#endif

static std::vector<long double> g_ref;  // reference results, computed once

struct Acc {
    acc_stats main = ACC_STATS_INIT, huge = ACC_STATS_INIT;  // |x| <= TRIG_THRESHOLD (and inf/NaN), larger |x|
};

static Acc accuracy(const std::vector<double> &x, const std::vector<double> &y) {
    Acc a;
    for (size_t i = 0; i < x.size(); i++) {
        const bool huge = std::fabs(x[i]) > TRIG_THRESHOLD && std::isfinite(x[i]);
        acc_add(huge ? &a.huge : &a.main, x[i], y[i], g_ref[i]);
    }
    return a;
}

// ns/value over 4096 cache-resident inputs: uniform in [-100, 100] (main path), and log-uniform over
// [2^30, 2^1000] (large-argument path)
static void bench(trig_array_fn f, double *ns_main, double *ns_huge) {
    const size_t N = 4096;
    static double in[N], out[N];
    for (size_t i = 0; i < N; i++) in[i] = rng_uniform(-100.0, 100.0);
    *ns_main = acc_bench(f, in, out, N, 2000);
    for (size_t i = 0; i < N; i++) {  // one generator call per statement (evaluation order is unspecified)
        const double m = rng_uniform(1.0, 2.0);
        const int e = 30 + (int)(rng() % 970);
        in[i] = std::ldexp(m, e);
    }
    *ns_huge = acc_bench(f, in, out, N, 200);
}

static void print_row(const char *name, const char *alg, size_t lanes, const Acc &a, bool do_bench, trig_array_fn f) {
    char l[24] = "", t1[24] = "", t2[24] = "";
    if (lanes) snprintf(l, sizeof l, "%zu", lanes);
    if (do_bench && f) {
        double m, h;
        bench(f, &m, &h);
        snprintf(t1, sizeof t1, "%.2f", m);
        snprintf(t2, sizeof t2, "%.1f", h);
    }
    const acc_stats &s = a.main, &h = a.huge;
    printf("%-26s %-10s %5s %8.4f %8.4f%% %6zu %4zu | %8.4f %8.4f%% %6zu | %8s %8s\n", name, alg, l, s.max_err,
           s.n ? 100.0 * s.not_cr / s.n : 0.0, s.over1, s.mismatched_special, h.max_err,
           h.n ? 100.0 * h.not_cr / h.n : 0.0, h.over1, t1, t2);
}

static int compare_bits(const char *what, const std::vector<double> &x, const std::vector<double> &y,
                        const std::vector<double> &ref) {
    size_t diff = 0;
    for (size_t i = 0; i < x.size(); i++) {
        if (dbits(y[i]) != dbits(ref[i]) && diff++ < 5) {
            printf("  %s: x = %a: %a, expected %a\n", what, x[i], y[i], ref[i]);
        }
    }
    if (diff) printf("  FAIL: %zu results differ from %s\n", diff, what);
    return diff != 0;
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
    std::vector<double> x(n_random + 40000);
    x.resize(make_inputs(x.data(), n_random));
    g_ref.resize(x.size());
    for (size_t i = 0; i < x.size(); i++) g_ref[i] = TRIG_REFL((long double)x[i]);
    std::vector<double> y(x.size()), first(x.size());
    FILE *hf = hash_file ? fopen(hash_file, "a") : nullptr;
    if (hf) fprintf(hf, "%s %s inputs - %zu %016llx\n", kArch, TRIG_NAME, x.size(), (unsigned long long)acc_hash(x.data(), x.size()));

    int failures = 0;
    std::vector<double> ref_c;  // the C/AVX-512 translation, if this CPU can run it
#if defined(__x86_64__)
    if (__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq")) {
        ref_c.resize(x.size());
        TRIG_SVML_ARRAY(x.data(), ref_c.data(), x.size());
        if (hf) fprintf(hf, "%s %s C-AVX512 svml %zu %016llx\n", kArch, TRIG_NAME, x.size(),
                        (unsigned long long)acc_hash(ref_c.data(), ref_c.size()));
    }
#endif

    printf("%s %s, %zu inputs (%s), errors in ulps vs %s (long double with %d-bit mantissa)\n", kArch, TRIG_NAME, x.size(),
           ref_c.empty() ? "no AVX-512: no direct comparison with the C translation" : "with AVX-512",
           TRIG_NAME "l", LDBL_MANT_DIG);
    printf("columns: max ulp, %% not correctly rounded, > 1 ulp, special-value mismatches for |x| <= %g (main path) and "
           "+-inf/NaN | the same for larger |x| (large-argument path) | ns/value for |x| <= 100 and for huge |x|\n\n",
           TRIG_THRESHOLD);
    printf("%-26s %-10s %5s %8s %9s %6s %4s | %8s %9s %6s | %8s %8s\n", "implementation", "algorithm", "lanes", "max ulp",
           "not CR", ">1ulp", "spec", "max ulp", "not CR", ">1ulp", "ns", "ns huge");
    bool have_first = false;
    for (int64_t target : hwy::SupportedAndGeneratedTargets()) {
        hwy::SetSupportedTargetsForTest(target);
        const int variant = TRIG_HWY_VARIANT();
        const char *tname = hwy::TargetName(target);
        TRIG_HWY_ARRAY(x.data(), y.data(), x.size());
        const Acc a = accuracy(x, y);
        char name[64];
        snprintf(name, sizeof name, "hwy %s", tname);
        print_row(name, variant ? "SVML" : "std::" TRIG_NAME, TRIG_HWY_LANES(), a, do_bench, TRIG_HWY_ARRAY);
        if (variant) {
            printf("  %s: worst x = %a (main path), %a (large path)\n", tname, a.main.worst_x, a.huge.worst_x);
            if (hf) fprintf(hf, "%s %s %s svml %zu %016llx\n", kArch, TRIG_NAME, tname, y.size(),
                            (unsigned long long)acc_hash(y.data(), y.size()));
            if (!(a.main.max_err < 1.0) || !(a.huge.max_err < 1.0) || a.main.over1 || a.huge.over1 ||
                a.main.mismatched_special || a.huge.mismatched_special) {
                printf("  FAIL: error bound exceeded or wrong special values\n");
                failures++;
            }
            if (!have_first) {
                first = y;
                have_first = true;
            } else {
                failures += compare_bits("the first SVML-algorithm target", x, y, first);
            }
            if (!ref_c.empty()) failures += compare_bits("the C/AVX-512 translation", x, y, ref_c);
            failures += TRIG_TARGET_CHECK(tname);
        }
    }
    hwy::SetSupportedTargetsForTest(0);

    struct extra_impl ex[16];
    const size_t nex = extra_impls(ex, 16);
    for (size_t k = 0; k < nex; k++) {
        ex[k].fn(x.data(), y.data(), x.size());
        print_row(ex[k].name, "", 0, accuracy(x, y), do_bench, ex[k].fn);
        if (!ref_c.empty()) {  // information only
            size_t diff = 0;
            for (size_t i = 0; i < x.size(); i++) diff += dbits(y[i]) != dbits(ref_c[i]);
            printf("  %zu results differ from the C/AVX-512 translation\n", diff);
        }
    }
    trig_array_fn scalar = [](const double *a, double *b, size_t n) {
        for (size_t i = 0; i < n; i++) b[i] = TRIG_STD(a[i]);
    };
    scalar(x.data(), y.data(), x.size());
    print_row("glibc " TRIG_NAME " (scalar)", "", 1, accuracy(x, y), do_bench, scalar);

    if (hf) fclose(hf);
    printf("\n%s\n", failures ? "FAILED" : "all checks passed");
    return failures ? 1 : 0;
}
