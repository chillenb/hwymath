// Test and benchmark of the Highway exp2 on every target that is compiled in and supported by this CPU,
// plus Highway's hwy/contrib/math Exp2 and the other implementations available here (extras.c).
//
// usage: hwy_test [N_RANDOM] [--no-bench] [--hash FILE]
//   --hash FILE  append "arch exp2 target variant n hash" lines (FNV-1a of the result bits; variant svml =
//                the SVML algorithm, also written for the C/AVX-512 translation when the CPU has AVX-512):
//                results must be bit-identical on every machine and target (CI cross-check).
// Exit status 1 if a check fails:
//   - SVML-algorithm targets: max error < 1 ulp, no special-value mismatch (NaN, inf, 0)
//   - all SVML-algorithm targets produce identical results
//   - x86-64 with AVX-512F/DQ: every SVML-algorithm target matches the C/AVX-512 translation bit for bit
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "hwy/targets.h"

extern "C" {
#include "../exp2_d.h"
#include "extras.h"
void hwy_exp2_array(const double *x, double *y, size_t len);
void hwy_contrib_exp2_array(const double *x, double *y, size_t len);
int hwy_exp2_variant(void);
size_t hwy_exp2_lanes(void);
}
#include "hwymath/accuracy.h"
#include "inputs.h"

#if defined(__x86_64__)
static const char *kArch = "x86_64";
#elif defined(__aarch64__)
static const char *kArch = "aarch64";
#else
static const char *kArch = "other";
#endif

static std::vector<long double> g_ref;  // exp2l(x), computed once

static acc_stats accuracy(const std::vector<double> &x, const std::vector<double> &y, double xmax) {
    acc_stats s = ACC_STATS_INIT;
    for (size_t i = 0; i < x.size(); i++) {
        if (!(x[i] <= xmax) && !std::isnan(x[i])) continue;
        acc_add(&s, x[i], y[i], g_ref[i]);
    }
    return s;
}

static double bench(acc_array_fn f) {  // ns/value, 4096 values in [-1000, 1000]
    const size_t N = 4096;
    static double in[N], out[N];
    for (size_t i = 0; i < N; i++) in[i] = rng_uniform(-1000.0, 1000.0);
    return acc_bench(f, in, out, N, 2000);
}

static void print_row(const char *name, const char *alg, size_t lanes, const acc_stats &s, double ns) {
    char l[24] = "", t[24] = "";
    if (lanes) snprintf(l, sizeof l, "%zu", lanes);
    if (ns > 0) snprintf(t, sizeof t, "%.2f", ns);
    printf("%-24s %-10s %5s %9.4f %9.4f%% %8zu %8zu %10s\n", name, alg, l, s.max_err,
           s.n ? 100.0 * s.not_cr / s.n : 0.0, s.over1, s.mismatched_special, t);
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
    std::vector<double> x(n_random + INPUTS_EXTRA);
    x.resize(make_inputs(x.data(), n_random));
    g_ref.resize(x.size());
    for (size_t i = 0; i < x.size(); i++) g_ref[i] = exp2l((long double)x[i]);
    std::vector<double> y(x.size());
    FILE *hf = hash_file ? fopen(hash_file, "a") : nullptr;
    if (hf) fprintf(hf, "%s exp2 inputs - %zu %016llx\n", kArch, x.size(), (unsigned long long)acc_hash(x.data(), x.size()));

    std::vector<double> ref_c;  // C/AVX-512 translation (bit-identical to Intel's assembly)
#if defined(__x86_64__)
    if (__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq")) {
        ref_c.resize(x.size());
        svml_exp2_ha_array(x.data(), ref_c.data(), x.size());
        if (hf) fprintf(hf, "%s exp2 C-AVX512 svml %zu %016llx\n", kArch, x.size(),
                        (unsigned long long)acc_hash(ref_c.data(), ref_c.size()));
    }
#endif

    printf("%s, %zu inputs, errors vs exp2l (long double with %d-bit mantissa)%s\n\n", kArch, x.size(), LDBL_MANT_DIG,
           ref_c.empty() ? "" : "; C/AVX-512 translation available for bit-identity checks");
    printf("%-24s %-10s %5s %9s %10s %8s %8s %10s\n", "implementation", "algorithm", "lanes", "max ulp", "not CR",
           ">1 ulp", "special", do_bench ? "ns/value" : "");
    int failures = 0;
    bool have_first = false;
    uint64_t first_hash = 0;
    for (int64_t target : hwy::SupportedAndGeneratedTargets()) {
        hwy::SetSupportedTargetsForTest(target);
        const int variant = hwy_exp2_variant();
        const char *tname = hwy::TargetName(target);
        hwy_exp2_array(x.data(), y.data(), x.size());
        acc_stats s = accuracy(x, y, INFINITY);
        char name[64];
        snprintf(name, sizeof name, "hwy %s", tname);
        print_row(name, variant ? "svml" : "std::exp2", hwy_exp2_lanes(), s, do_bench ? bench(hwy_exp2_array) : 0);

        if (variant) {
            const uint64_t h = acc_hash(y.data(), y.size());
            if (hf) fprintf(hf, "%s exp2 %s svml %zu %016llx\n", kArch, tname, y.size(), (unsigned long long)h);
            if (s.max_err >= 1.0 || s.mismatched_special) { printf("  FAIL: error bound exceeded or special mismatch\n"); failures++; }
            if (!have_first) { first_hash = h; have_first = true; }
            else if (h != first_hash) { printf("  FAIL: results differ from the first SVML-algorithm target\n"); failures++; }
            if (!ref_c.empty()) {
                size_t diff = 0;
                for (size_t i = 0; i < x.size(); i++) {
                    if (dbits(y[i]) != dbits(ref_c[i]) && diff++ < 5)
                        printf("  x=%a: hwy %a  C/AVX-512 %a\n", x[i], y[i], ref_c[i]);
                }
                if (diff) { printf("  FAIL: %zu results differ from the C/AVX-512 translation\n", diff); failures++; }
                else printf("  bit-identical to the C/AVX-512 translation\n");
            }
        }
        hwy_contrib_exp2_array(x.data(), y.data(), x.size());
        print_row("  hwy contrib Exp2", "", 0, accuracy(x, y, 1024.0), do_bench ? bench(hwy_contrib_exp2_array) : 0);
    }
    hwy::SetSupportedTargetsForTest(0);

    struct extra_impl ex[16];
    const size_t nex = extra_impls(ex, 16);
    for (size_t k = 0; k < nex; k++) {
        ex[k].fn(x.data(), y.data(), x.size());
        print_row(ex[k].name, "", 0, accuracy(x, y, INFINITY), do_bench ? bench(ex[k].fn) : 0);
    }
    for (size_t i = 0; i < x.size(); i++) y[i] = std::exp2(x[i]);
    print_row("glibc exp2 (scalar)", "", 1, accuracy(x, y, INFINITY),
              do_bench ? bench([](const double *a, double *b, size_t n) { for (size_t i = 0; i < n; i++) b[i] = std::exp2(a[i]); }) : 0);

    if (hf) fclose(hf);
    printf("\n%s\n", failures ? "FAILED" : "all checks passed");
    return failures ? 1 : 0;
}
