// Test and benchmark of the two Highway exps (hwy_exp.cc: hp algorithm; hwy_exp_svml.cc: Intel's SVML
// algorithm) on every target that is compiled in and supported by this CPU, plus Highway's
// hwy/contrib/math Exp and the other implementations available here (extras*.c).
//
// usage: hwy_test [N_RANDOM] [--no-bench] [--hash FILE]
//   --hash FILE  append "arch exp target variant n hash" lines (FNV-1a of the result bits; variant hp4/hp8 =
//                hp algorithm with a 4/8-entry table, svml = Intel's algorithm, also written for the C
//                AVX2 port / AVX-512 translation and Intel's assembly when the CPU can run them): results
//                of the same variant must be bit-identical on every machine and target, which the CI
//                cross-check job verifies.
// Exit status 1 if a check fails:
//   - hp targets: max error <= 0.53 ulp and no error > 1 ulp
//   - all targets with the same hp variant produce identical results
//   - x86-64: the AVX2 target matches the C hp variant bit for bit, and all hp targets match the
//     bit-identical SVML port on special inputs/results (NaN, inf, 0, denormal)
//   - svml targets: max error < 1 ulp, no special-value mismatch (NaN, inf, 0), all identical; on x86-64
//     bit-identical to the C AVX2 port and, with AVX-512F, to the C/AVX-512 translation and to Intel's
//     __svml_exp8_ha assembly (which are also compared with each other)
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "hwy/targets.h"

extern "C" {
#include "../exp_d.h"
#include "extras.h"
void hwy_exp_array(const double *x, double *y, size_t len);
void hwy_contrib_exp_array(const double *x, double *y, size_t len);
int hwy_exp_variant(void);
size_t hwy_exp_lanes(void);
void hwy_exp_svml_array(const double *x, double *y, size_t len);
int hwy_exp_svml_variant(void);
size_t hwy_exp_svml_lanes(void);
}
#include "inputs.h"

#if defined(__x86_64__)
static const char *kArch = "x86_64";
#elif defined(__aarch64__)
static const char *kArch = "aarch64";
#else
static const char *kArch = "other";
#endif

#if defined(__x86_64__)
extern "C" __attribute__((target("avx512f"))) __m512d __svml_exp8_ha(__m512d x);  // Intel's assembly
__attribute__((target("avx512f"))) static void intel_asm_array(const double *x, double *y, size_t len) {
    size_t i = 0;
    for (; i + 8 <= len; i += 8) _mm512_storeu_pd(y + i, __svml_exp8_ha(_mm512_loadu_pd(x + i)));
    if (i < len) {
        const __mmask8 m = (__mmask8)((1u << (len - i)) - 1);
        _mm512_mask_storeu_pd(y + i, m, __svml_exp8_ha(_mm512_maskz_loadu_pd(m, x + i)));
    }
}
#endif

static double ulp_err(double y, long double ref) {
    double rd = (double)ref;
    long double ulp = (long double)std::nextafter(std::fabs(rd), INFINITY) - fabsl((long double)rd);
    if (std::fabs(rd) == DBL_MAX) ulp = 0x1p971L;
    return (double)(fabsl((long double)y - ref) / ulp);
}

struct Stats { double max_err = 0, worst_x = 0; size_t n = 0, not_cr = 0, over1 = 0; };

static std::vector<long double> g_ref;  // expl(x), computed once

static Stats accuracy(const std::vector<double> &x, const std::vector<double> &y, double xmax) {
    Stats s;
    for (size_t i = 0; i < x.size(); i++) {
        const long double ref = g_ref[i];
        if (!std::isfinite(x[i]) || x[i] > xmax || !std::isfinite((double)ref) || (double)ref == 0.0) continue;
        double e = ulp_err(y[i], ref);
        s.n++;
        s.not_cr += e > 0.5;
        s.over1 += e > 1.0;
        if (e > s.max_err) { s.max_err = e; s.worst_x = x[i]; }
    }
    return s;
}

// results whose reference is not a finite nonzero double must match exactly (NaN, same inf, same zero)
static size_t special_mismatches(const std::vector<double> &y) {
    size_t bad = 0;
    for (size_t i = 0; i < y.size(); i++) {
        const double rd = (double)g_ref[i];
        if (std::isfinite(rd) && rd != 0.0) continue;
        bad += std::isnan(rd) ? !std::isnan(y[i]) : !(y[i] == rd && std::signbit(y[i]) == std::signbit(rd));
    }
    return bad;
}

// number of results that differ from ref bit for bit (printing the first few)
static size_t bit_diff(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &ref,
                       const char *name) {
    size_t diff = 0;
    for (size_t i = 0; i < x.size(); i++) {
        if (dbits(y[i]) != dbits(ref[i]) && diff++ < 5) printf("  x=%a: %a, %s %a\n", x[i], y[i], name, ref[i]);
    }
    return diff;
}

static uint64_t fnv1a(const std::vector<double> &v) {
    uint64_t h = 0xcbf29ce484222325ULL;
    const unsigned char *p = reinterpret_cast<const unsigned char *>(v.data());
    for (size_t i = 0; i < v.size() * sizeof(double); i++) h = (h ^ p[i]) * 0x100000001b3ULL;
    return h;
}

static double now() {
    timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static double bench(void (*f)(const double *, double *, size_t)) {  // ns/value, 4096 values in [-700, 700]
    const size_t N = 4096;
    static double in[N], out[N];
    for (size_t i = 0; i < N; i++) in[i] = ((double)(rng() >> 11) * 0x1p-53) * 1400.0 - 700.0;
    double best = 1e9;
    for (int trial = 0; trial < 7; trial++) {
        double t0 = now();
        for (int r = 0; r < 2000; r++) f(in, out, N);
        best = std::fmin(best, (now() - t0) / (2000.0 * N));
    }
    return best * 1e9;
}

static void print_row(const char *name, const char *alg, size_t lanes, const Stats &s, double ns) {
    char l[24] = "", t[24] = "";
    if (lanes) snprintf(l, sizeof l, "%zu", lanes);
    if (ns > 0) snprintf(t, sizeof t, "%.2f", ns);
    printf("%-22s %-12s %5s %9.4f %9.4f%% %8zu %10s\n", name, alg, l, s.max_err, 100.0 * s.not_cr / s.n, s.over1, t);
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
    std::vector<double> x(n_random + 4096);
    x.resize(make_inputs(x.data(), n_random));
    g_ref.resize(x.size());
    for (size_t i = 0; i < x.size(); i++) g_ref[i] = expl((long double)x[i]);
    std::vector<double> y(x.size());
    FILE *hf = hash_file ? fopen(hash_file, "a") : nullptr;
    if (hf) fprintf(hf, "%s exp inputs - %zu %016llx\n", kArch, x.size(), (unsigned long long)fnv1a(x));

    int failures = 0;
    // references for the svml variant (x86-64): C AVX2 port, C/AVX-512 translation, Intel's assembly
    std::vector<std::pair<const char *, std::vector<double>>> svml_refs;
#if defined(__x86_64__)
    const bool have_avx2 = __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
    std::vector<double> ref_hp(x.size()), ref_ha(x.size());
    if (have_avx2) {
        svml_exp_hp_array_avx2(x.data(), ref_hp.data(), x.size());  // C hp variant
        svml_exp_ha_array_avx2(x.data(), ref_ha.data(), x.size());  // bit-identical SVML port
        svml_refs.emplace_back("C-AVX2", ref_ha);
    }
    if (__builtin_cpu_supports("avx512f")) {
        std::vector<double> r(x.size());
        svml_exp_ha_array(x.data(), r.data(), x.size());
        svml_refs.emplace_back("C-AVX512", r);
        intel_asm_array(x.data(), r.data(), x.size());
        svml_refs.emplace_back("asm-AVX512", r);
    }
#endif
    for (size_t k = 0; k < svml_refs.size(); k++) {  // the references must agree with each other
        const auto &ref = svml_refs[k].second;
        if (hf) fprintf(hf, "%s exp %s svml %zu %016llx\n", kArch, svml_refs[k].first, ref.size(),
                        (unsigned long long)fnv1a(ref));
        if (k > 0) {
            const size_t diff = bit_diff(x, ref, svml_refs[0].second, svml_refs[0].first);
            if (diff) { printf("FAIL: %s differs from %s in %zu results\n", svml_refs[k].first, svml_refs[0].first, diff); failures++; }
        }
    }

    printf("%s, %zu inputs, errors vs expl (long double with %d-bit mantissa)%s\n\n", kArch, x.size(), LDBL_MANT_DIG,
           svml_refs.size() > 1 ? "; C/AVX-512 translation and Intel's assembly available for bit-identity checks" : "");
    printf("%-22s %-12s %5s %9s %10s %8s %10s\n", "implementation", "algorithm", "lanes", "max ulp", "not CR", ">1 ulp",
           do_bench ? "ns/value" : "");
    std::map<int, uint64_t> variant_hash;
    bool have_svml_hash = false;
    uint64_t svml_hash = 0;  // hp variant -> result hash of the first target with it
    for (int64_t target : hwy::SupportedAndGeneratedTargets()) {
        hwy::SetSupportedTargetsForTest(target);
        const int variant = hwy_exp_variant();
        const char *tname = hwy::TargetName(target);
        hwy_exp_array(x.data(), y.data(), x.size());
        Stats s = accuracy(x, y, INFINITY);
        char name[64], alg[16];
        snprintf(name, sizeof name, "hwy %s", tname);
        snprintf(alg, sizeof alg, variant ? "hp, %d-entry" : "std::exp", 1 << variant);
        print_row(name, alg, hwy_exp_lanes(), s, do_bench ? bench(hwy_exp_array) : 0);

        if (variant) {
            const uint64_t h = fnv1a(y);
            if (hf) fprintf(hf, "%s exp %s hp%d %zu %016llx\n", kArch, tname, 1 << variant, y.size(), (unsigned long long)h);
            if (s.max_err > 0.53 || s.over1) { printf("  FAIL: error bound exceeded\n"); failures++; }
            auto it = variant_hash.find(variant);
            if (it == variant_hash.end()) variant_hash[variant] = h;
            else if (it->second != h) { printf("  FAIL: results differ from the first %d-entry target\n", 1 << variant); failures++; }
#if defined(__x86_64__)
            if (have_avx2) {
                size_t diff_special = 0, diff_hp = 0;
                for (size_t i = 0; i < x.size(); i++) {
                    const double r = ref_ha[i];
                    if ((!std::isfinite(x[i]) || !std::isfinite(r) || r == 0.0 || std::fabs(r) < 0x1p-1022) &&
                        dbits(y[i]) != dbits(r)) diff_special++;
                    diff_hp += dbits(y[i]) != dbits(ref_hp[i]);
                }
                if (diff_special) { printf("  FAIL: %zu special inputs differ from the SVML port\n", diff_special); failures++; }
                if (target == HWY_AVX2 && diff_hp) { printf("  FAIL: %zu results differ from the C hp variant\n", diff_hp); failures++; }
            }
#endif
        }

        // Intel's algorithm (hwy_exp_svml.cc)
        const int svariant = hwy_exp_svml_variant();
        hwy_exp_svml_array(x.data(), y.data(), x.size());
        s = accuracy(x, y, INFINITY);
        print_row(name, svariant ? "svml" : "std::exp", hwy_exp_svml_lanes(), s,
                  do_bench ? bench(hwy_exp_svml_array) : 0);
        if (svariant) {
            const uint64_t h = fnv1a(y);
            if (hf) fprintf(hf, "%s exp %s svml %zu %016llx\n", kArch, tname, y.size(), (unsigned long long)h);
            const size_t sm = special_mismatches(y);
            if (s.max_err >= 1.0 || s.over1 || sm) {
                printf("  FAIL: error bound exceeded or %zu special mismatches\n", sm);
                failures++;
            }
            if (!have_svml_hash) { svml_hash = h; have_svml_hash = true; }
            else if (h != svml_hash) { printf("  FAIL: results differ from the first svml target\n"); failures++; }
            for (const auto &ref : svml_refs) {
                const size_t diff = bit_diff(x, y, ref.second, ref.first);
                if (diff) { printf("  FAIL: %zu results differ from %s\n", diff, ref.first); failures++; }
                else printf("  svml: bit-identical to %s\n", ref.first);
            }
        }

        hwy_contrib_exp_array(x.data(), y.data(), x.size());
        Stats c = accuracy(x, y, 706.0);  // documented valid range of hwy contrib Exp: x <= 706
        snprintf(name, sizeof name, "  hwy contrib Exp");
        print_row(name, "", 0, c, do_bench ? bench(hwy_contrib_exp_array) : 0);
    }
    hwy::SetSupportedTargetsForTest(0);

    struct extra_impl ex[16];
    const size_t nex = extra_impls(ex, 16);
    for (size_t k = 0; k < nex; k++) {
        ex[k].fn(x.data(), y.data(), x.size());
        print_row(ex[k].name, "", 0, accuracy(x, y, INFINITY), do_bench ? bench(ex[k].fn) : 0);
    }
    for (size_t i = 0; i < x.size(); i++) y[i] = std::exp(x[i]);
    print_row("glibc exp (scalar)", "", 1, accuracy(x, y, INFINITY),
              do_bench ? bench([](const double *a, double *b, size_t n) { for (size_t i = 0; i < n; i++) b[i] = std::exp(a[i]); }) : 0);

    if (hf) fclose(hf);
    printf("\n%s\n", failures ? "FAILED" : "all checks passed");
    return failures ? 1 : 0;
}
