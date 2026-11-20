// Test and benchmark of the Highway asinh on every target that is compiled in and supported by this CPU,
// plus the other implementations available here (extras*.c).
//
// usage: hwy_test [N_RANDOM] [--no-bench] [--hash FILE]
//   --hash FILE  append "arch asinh target svml n hash" lines (FNV-1a of the result bits) for every target
//                with the SVML algorithm, and for the C/AVX-512 translation if the CPU has AVX-512F: all of
//                them must be bit-identical on every machine, which the CI cross-check job verifies.
// Exit status 1 if a check fails:
//   - SVML-algorithm targets: max error < 1 ulp (vs asinhl), special results (inf, NaN, +-0) exact
//   - all SVML-algorithm targets produce identical results
//   - with AVX-512F (or under Intel SDE): every SVML-algorithm target is bit-identical to the C translation
//     of Intel's code (itself checked against the assembly by test_asinh8)
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "hwy/targets.h"

extern "C" {
#include "../asinh_d.h"
#include "extras.h"
#include "hwymath/accuracy.h"
void hwy_asinh_array(const double *x, double *y, size_t len);
int hwy_asinh_variant(void);
size_t hwy_asinh_lanes(void);
}
#include "inputs.h"

#define FNAME "asinh"
static long double ref_fn(long double x) { return asinhl(x); }
static const double kBenchLo = -20.0, kBenchHi = 20.0;  // benchmark inputs: uniform in this range

#if defined(__x86_64__)
static const char *kArch = "x86_64";
#elif defined(__aarch64__)
static const char *kArch = "aarch64";
#else
static const char *kArch = "other";
#endif

static std::vector<long double> g_ref;

static acc_stats accuracy(const std::vector<double> &x, const std::vector<double> &y) {
    acc_stats s = ACC_STATS_INIT;
    for (size_t i = 0; i < x.size(); i++) acc_add(&s, x[i], y[i], g_ref[i]);
    return s;
}

static double bench(void (*f)(const double *, double *, size_t)) {
    const size_t N = 4096;
    static double in[N], out[N];
    for (size_t i = 0; i < N; i++) in[i] = rng_uniform(kBenchLo, kBenchHi);
    return acc_bench(f, in, out, N, 1000);
}

static void print_row(const char *name, const char *alg, size_t lanes, const acc_stats &s, double ns) {
    char l[24] = "", t[24] = "";
    if (lanes) snprintf(l, sizeof l, "%zu", lanes);
    if (ns > 0) snprintf(t, sizeof t, "%.2f", ns);
    printf("%-24s %-12s %5s %9.4f %9.4f%% %8zu %8zu %10s\n", name, alg, l, s.max_err, s.n ? 100.0 * s.not_cr / s.n : 0.0,
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
    std::vector<double> x(n_random + 16384);
    x.resize(make_inputs(x.data(), n_random));
    g_ref.resize(x.size());
    for (size_t i = 0; i < x.size(); i++) g_ref[i] = ref_fn((long double)x[i]);
    std::vector<double> y(x.size());
    FILE *hf = hash_file ? fopen(hash_file, "a") : nullptr;
    if (hf) fprintf(hf, "%s %s inputs - %zu %016llx\n", kArch, FNAME, x.size(), (unsigned long long)acc_hash(x.data(), x.size()));

    // the C translation of Intel's AVX-512 code (bit-identical to the assembly), if this CPU can run it
    std::vector<double> ref_c;
#if defined(__x86_64__)
    if (__builtin_cpu_supports("avx512f")) {
        ref_c.resize(x.size());
        svml_asinh_ha_array(x.data(), ref_c.data(), x.size());
        if (hf) fprintf(hf, "%s %s C-AVX512 svml %zu %016llx\n", kArch, FNAME, x.size(),
                        (unsigned long long)acc_hash(ref_c.data(), x.size()));
    }
#endif

    printf("%s, %zu inputs, errors vs %sl (long double with %d-bit mantissa)%s\n\n", kArch, x.size(), FNAME,
           LDBL_MANT_DIG, ref_c.empty() ? "" : "; SVML targets compared with the C/AVX-512 translation");
    printf("%-24s %-12s %5s %9s %10s %8s %8s %10s\n", "implementation", "algorithm", "lanes", "max ulp", "not CR",
           ">1 ulp", "special", do_bench ? "ns/value" : "");
    int failures = 0;
    bool have_hash = false;
    uint64_t first_hash = 0;
    for (int64_t target : hwy::SupportedAndGeneratedTargets()) {
        hwy::SetSupportedTargetsForTest(target);
        const int variant = hwy_asinh_variant();
        const char *tname = hwy::TargetName(target);
        hwy_asinh_array(x.data(), y.data(), x.size());
        const acc_stats s = accuracy(x, y);
        char name[64];
        snprintf(name, sizeof name, "hwy %s", tname);
        print_row(name, variant ? "SVML" : "std::" FNAME, hwy_asinh_lanes(), s, do_bench ? bench(hwy_asinh_array) : 0);
        if (s.mismatched_special) { printf("  FAIL: special results differ from the reference\n"); failures++; }
        if (!variant) continue;
        const uint64_t h = acc_hash(y.data(), y.size());
        if (hf) fprintf(hf, "%s %s %s svml %zu %016llx\n", kArch, FNAME, tname, y.size(), (unsigned long long)h);
        printf("  max error at x = %a\n", s.worst_x);
        if (s.max_err >= 1.0) { printf("  FAIL: error bound (1 ulp) exceeded at x = %a\n", s.worst_x); failures++; }
        if (!have_hash) { have_hash = true; first_hash = h; }
        else if (h != first_hash) { printf("  FAIL: results differ from the first SVML target\n"); failures++; }
        if (!ref_c.empty()) {
            size_t diff = 0;
            for (size_t i = 0; i < x.size(); i++) {
                if (dbits(y[i]) != dbits(ref_c[i]) && diff++ < 5)
                    printf("  x = %a: %a, C translation %a\n", x[i], y[i], ref_c[i]);
            }
            if (diff) { printf("  FAIL: %zu results differ from the C/AVX-512 translation\n", diff); failures++; }
            else printf("  bit-identical to the C/AVX-512 translation\n");
        }
    }
    hwy::SetSupportedTargetsForTest(0);

    struct extra_impl ex[16];
    const size_t nex = extra_impls(ex, 16);
    for (size_t k = 0; k < nex; k++) {
        ex[k].fn(x.data(), y.data(), x.size());
        print_row(ex[k].name, "", 0, accuracy(x, y), do_bench ? bench(ex[k].fn) : 0);
        if (!ref_c.empty() && strstr(ex[k].name, "Intel")) {
            size_t diff = 0;
            for (size_t i = 0; i < x.size(); i++) diff += dbits(y[i]) != dbits(ref_c[i]) && !(std::isnan(y[i]) && std::isnan(ref_c[i]));
            printf("  %zu results differ from the C/AVX-512 translation (NaN payloads ignored)\n", diff);
        }
    }
    for (size_t i = 0; i < x.size(); i++) y[i] = std::asinh(x[i]);
    print_row("glibc " FNAME " (scalar)", "", 1, accuracy(x, y),
              do_bench ? bench([](const double *a, double *b, size_t n) { for (size_t i = 0; i < n; i++) b[i] = std::asinh(a[i]); }) : 0);

    if (hf) fclose(hf);
    printf("\n%s\n", failures ? "FAILED" : "all checks passed");
    return failures ? 1 : 0;
}
