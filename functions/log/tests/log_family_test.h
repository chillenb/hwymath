// Test and benchmark of the Highway version of one log-family function (log, log2, log10, log1p) on every
// target that is compiled in and supported by this CPU, plus Highway's hwy/contrib/math version and the
// other implementations available here (extras*.c). Shared by functions/{log,log2,log10,log1p}/tests/
// hwy_test.cc, which define before including it:
//   FN_NAME "log"            REF(x)     long double reference, e.g. logl(x)
//   LIBM(x) scalar libm      HWY_ARRAY / CONTRIB_ARRAY / VARIANT / LANES   the hwy_<f>.cc exports
//   CALLOUT  the C scalar callout,  SVML_ARRAY  the C/AVX-512 translation (x86-64), SVML_ISA_OK()
//   BENCH_LO, BENCH_HI   benchmark input range
//
// usage: hwy_test [N_RANDOM] [--no-bench] [--hash FILE]
//   --hash FILE  append "arch <f> target svml n hash" lines (FNV-1a of the result bits; variant svml = the SVML
//                algorithm, bit-identical to Intel's AVX-512 code), "arch <f> C-AVX512 svml n hash" for the
//                C/AVX-512 translation when the CPU has AVX-512, and "arch <f> inputs - n hash".
// Exit status 1 if a check fails:
//   - SVML-algorithm targets: max error < 1 ulp (HA) and no special-value mismatch vs the reference
//   - all SVML-algorithm targets produce identical results
//   - special inputs (x outside the vector path): identical to the C scalar callout
//   - with AVX-512 (or under Intel SDE): every SVML-algorithm target matches the C/AVX-512 translation,
//     i.e. Intel's __svml_<f>8_ha, bit for bit
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "hwy/targets.h"
#include "hwymath/accuracy.h"

#if defined(__x86_64__)
static const char *kArch = "x86_64";
#elif defined(__aarch64__)
static const char *kArch = "aarch64";
#else
static const char *kArch = "other";
#endif

static void print_row(const char *name, const char *alg, size_t lanes, const acc_stats &s, double ns) {
    char l[24] = "", t[24] = "";
    if (lanes) snprintf(l, sizeof l, "%zu", lanes);
    if (ns > 0) snprintf(t, sizeof t, "%.2f", ns);
    printf("%-24s %-10s %5s %9.4f %9.4f%% %8zu %8zu %10s\n", name, alg, l, s.max_err, s.n ? 100.0 * s.not_cr / s.n : 0.0,
           s.over1, s.mismatched_special, t);
}

static std::vector<long double> g_ref;
static std::vector<double> g_bench_in;

static acc_stats accuracy(const std::vector<double> &x, const std::vector<double> &y) {
    acc_stats s = ACC_STATS_INIT;
    for (size_t i = 0; i < x.size(); i++) acc_add(&s, x[i], y[i], g_ref[i]);
    return s;
}

static double bench(acc_array_fn f) {
    static double out[4096];
    return acc_bench(f, g_bench_in.data(), out, g_bench_in.size(), 2000);
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
    for (size_t i = 0; i < x.size(); i++) g_ref[i] = REF((long double)x[i]);
    for (int i = 0; i < 4096; i++) g_bench_in.push_back(rng_uniform(BENCH_LO, BENCH_HI));
    std::vector<double> y(x.size());
    FILE *hf = hash_file ? fopen(hash_file, "a") : nullptr;
    if (hash_file && !hf) { perror(hash_file); return 1; }
    if (hf) fprintf(hf, "%s %s inputs - %zu %016llx\n", kArch, FN_NAME, x.size(), (unsigned long long)acc_hash(x.data(), x.size()));

    // results of the scalar callout for the inputs outside the vector path
    std::vector<size_t> special;
    std::vector<double> special_ref;
    for (size_t i = 0; i < x.size(); i++) {
        if (IS_SPECIAL(x[i])) {
            double r;
            CALLOUT(&x[i], &r);
            special.push_back(i);
            special_ref.push_back(r);
        }
    }

    std::vector<double> ref_svml;  // the C/AVX-512 translation, bit-identical to Intel's assembly
#if defined(__x86_64__)
    if (SVML_ISA_OK()) {
        ref_svml.resize(x.size());
        SVML_ARRAY(x.data(), ref_svml.data(), x.size());
        if (hf) fprintf(hf, "%s %s C-AVX512 svml %zu %016llx\n", kArch, FN_NAME, x.size(),
                        (unsigned long long)acc_hash(ref_svml.data(), ref_svml.size()));
    }
#endif

    printf("%s %s, %zu inputs (%zu special), errors vs a long double reference (%d-bit mantissa)%s\n\n", kArch, FN_NAME,
           x.size(), special.size(), LDBL_MANT_DIG,
           ref_svml.empty() ? "" : "; AVX-512: every target compared bit for bit with the C/AVX-512 translation");
    printf("%-24s %-10s %5s %9s %10s %8s %8s %10s\n", "implementation", "algorithm", "lanes", "max ulp", "not CR",
           ">1 ulp", "special", do_bench ? "ns/value" : "");
    int failures = 0;
    bool have_first = false;
    uint64_t first_hash = 0;
    std::vector<double> first;  // results of the first SVML-algorithm target
    for (int64_t target : hwy::SupportedAndGeneratedTargets()) {
        hwy::SetSupportedTargetsForTest(target);
        const int variant = VARIANT();
        const char *tname = hwy::TargetName(target);
        HWY_ARRAY(x.data(), y.data(), x.size());
        const acc_stats s = accuracy(x, y);
        char name[64];
        snprintf(name, sizeof name, "hwy %s", tname);
        print_row(name, variant ? "SVML" : "libm", LANES(), s, do_bench ? bench(HWY_ARRAY) : 0);

        if (variant) {
            const uint64_t h = acc_hash(y.data(), y.size());
            if (hf) fprintf(hf, "%s %s %s svml %zu %016llx\n", kArch, FN_NAME, tname, y.size(), (unsigned long long)h);
            if (s.max_err >= 1.0 || s.over1 || s.mismatched_special) { printf("  FAIL: error bound exceeded\n"); failures++; }
            if (!have_first) { have_first = true; first_hash = h; first = y; }
            else if (h != first_hash) { printf("  FAIL: results differ from the first SVML target\n"); failures++; }
            size_t diff_special = 0;
            for (size_t k = 0; k < special.size(); k++) diff_special += dbits(y[special[k]]) != dbits(special_ref[k]);
            if (diff_special) { printf("  FAIL: %zu special inputs differ from the scalar callout\n", diff_special); failures++; }
            if (!ref_svml.empty()) {
                size_t diff = 0;
                for (size_t i = 0; i < x.size(); i++) {
                    if (dbits(y[i]) != dbits(ref_svml[i]) && diff++ < 5) {
                        printf("  x = %a: hwy %a, C/AVX-512 %a\n", x[i], y[i], ref_svml[i]);
                    }
                }
                if (diff) { printf("  FAIL: %zu results differ from the C/AVX-512 translation\n", diff); failures++; }
            }
        }
        CONTRIB_ARRAY(x.data(), y.data(), x.size());
        print_row("  hwy contrib", "", 0, accuracy(x, y), do_bench ? bench(CONTRIB_ARRAY) : 0);
    }
    hwy::SetSupportedTargetsForTest(0);

    struct extra_impl ex[16];
    const size_t nex = extra_impls(ex, 16);
    for (size_t k = 0; k < nex; k++) {
        ex[k].fn(x.data(), y.data(), x.size());
        print_row(ex[k].name, "", 0, accuracy(x, y), do_bench ? bench(ex[k].fn) : 0);
        if (!first.empty() && strstr(ex[k].name, "Intel SVML HA AVX-512")) {  // information only
            size_t diff = 0, diff_nan = 0;
            for (size_t i = 0; i < x.size(); i++) {
                if (dbits(y[i]) == dbits(first[i])) continue;
                if (std::isnan(y[i]) && std::isnan(first[i])) diff_nan++;
                else diff++;
            }
            printf("  (Intel's libsvml binary vs this version: %zu different results, %zu NaNs with other bits)\n", diff, diff_nan);
        }
    }
    for (size_t i = 0; i < x.size(); i++) y[i] = LIBM(x[i]);
    print_row("glibc (scalar)", "", 1, accuracy(x, y),
              do_bench ? bench([](const double *a, double *b, size_t n) { for (size_t i = 0; i < n; i++) b[i] = LIBM(a[i]); }) : 0);

    if (hf) fclose(hf);
    printf("\ncolumns: max error and share of results that are not correctly rounded, results > 1 ulp, special-value\n"
           "mismatches (NaN/inf/zero results); ns/value: best of 7 runs over 4096 values in [%g, %g]\n",
           (double)BENCH_LO, (double)BENCH_HI);
    printf("\n%s\n", failures ? "FAILED" : "all checks passed");
    return failures ? 1 : 0;
}
