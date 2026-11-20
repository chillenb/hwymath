// Test and benchmark of the Highway pow on every target that is compiled in and supported by this CPU, plus
// the other implementations available here (extras*.c) and glibc's scalar pow.
//
// usage: hwy_test [N_RANDOM] [--no-bench] [--hash FILE]
//   --hash FILE  append "arch pow target variant n hash" lines (FNV-1a of the result bits, NaNs made
//                canonical because the default NaN's sign differs between x86-64 and AArch64; variant svml =
//                the SVML algorithm, bit-identical to Intel's AVX-512 code) and "arch pow inputs - n hash";
//                with AVX-512, also the C translation's hash as target C-AVX512: results of the same variant
//                must be bit-identical on every machine and target, which the CI cross-check verifies.
// Exit status 1 if a check fails:
//   - SVML-algorithm targets: no error > 1 ulp against powl, special values (NaN, inf, 0) as powl
//   - all SVML-algorithm targets produce identical results
//   - with AVX-512 (or under Intel SDE): every SVML-algorithm target matches the C translation of the
//     assembly (svml_pow_ha_array) bit for bit
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "hwy/targets.h"

extern "C" {
#include "../pow_d.h"
#include "extras.h"
#include "hwymath/accuracy.h"
void hwy_pow_array(const double *x, const double *y, double *z, size_t len);
int hwy_pow_variant(void);
void hwy_pow_primitive(int op, const double *a, const double *b, const double *c, double *out, size_t len);
size_t hwy_pow_lanes(void);
}
#include "inputs.h"

#if defined(__x86_64__)
static const char *kArch = "x86_64";
#elif defined(__aarch64__)
static const char *kArch = "aarch64";
#else
static const char *kArch = "other";
#endif

static std::vector<double> g_x, g_y;
static std::vector<long double> g_ref;  // powl(x, y), computed once

// hash with every NaN replaced by one canonical NaN: the default NaN of an invalid operation (0/0 for
// pow(x < 0, non-integer y)) is negative on x86-64 and positive on AArch64; all other bits must agree
static uint64_t result_hash(const std::vector<double> &z) {
    std::vector<double> c(z);
    for (double &v : c) if (std::isnan(v)) v = bitsd(0x7ff8000000000000ULL);
    return acc_hash(c.data(), c.size());
}

static bool is_snan(double v) { return std::isnan(v) && !((dbits(v) >> 51) & 1); }

static double g_worst_y;  // y of the worst case of the last accuracy() call (acc_stats only keeps x)

static acc_stats accuracy(const std::vector<double> &z) {
    acc_stats s = ACC_STATS_INIT;
    g_worst_y = 0;
    for (size_t i = 0; i < z.size(); i++) {
        if (is_snan(g_x[i]) || is_snan(g_y[i])) continue;  // pow(sNaN, 0): NaN in glibc, 1 in SVML and powl
        const double before = s.max_err;
        acc_add(&s, g_x[i], z[i], g_ref[i]);
        if (s.max_err != before) g_worst_y = g_y[i];
    }
    return s;
}

static double bench(pow_array_fn f) {  // ns/value, 4096 pairs: x in [0.01, 100], y in [-50, 50]
    const size_t N = 4096;
    static double in1[N], in2[N], out[N];
    for (size_t i = 0; i < N; i++) {
        in1[i] = rng_uniform(0.01, 100);
        in2[i] = rng_uniform(-50, 50);
    }
    return acc_bench2(f, in1, in2, out, N, 300);
}

static void print_row(const char *name, size_t lanes, const acc_stats &s, double ns) {
    char l[24] = "", t[24] = "";
    if (lanes) snprintf(l, sizeof l, "%zu", lanes);
    if (ns > 0) snprintf(t, sizeof t, "%.2f", ns);
    printf("%-26s %5s %9.4f %9.4f%% %8zu %8zu %10s   worst: pow(%a, %a)\n", name, l, s.max_err,
           s.n ? 100.0 * s.not_cr / s.n : 0.0, s.over1, s.mismatched_special, t, s.worst_x, g_worst_y);
}

#if defined(__x86_64__)
// the AVX-512 instructions that the Highway version emulates (op numbers as in hwy_pow_primitive)
__attribute__((target("avx512f,avx512dq"))) static void ref_primitive(int op, const double *a, const double *b,
                                                                       const double *c, double *out, size_t n) {
    const int RZ = _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC, RD = _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC;
    for (size_t i = 0; i + 8 <= n; i += 8) {
        const __m512d va = _mm512_loadu_pd(a + i), vb = _mm512_loadu_pd(b + i), vc = _mm512_loadu_pd(c + i);
        __m512d r;
        switch (op) {
            case 0: r = _mm512_roundscale_round_pd(_mm512_rcp14_pd(va), 0x58, _MM_FROUND_NO_EXC); break;
            case 1: r = _mm512_mul_round_pd(va, vb, RZ); break;
            case 2: r = _mm512_fmadd_round_pd(va, vb, vc, RZ); break;
            case 3: r = _mm512_add_round_pd(va, vb, RZ); break;
            case 4: r = _mm512_reduce_round_pd(va, 0x41, _MM_FROUND_NO_EXC); break;
            default: r = _mm512_add_round_pd(va, _mm512_set1_pd(0x1.8000000003ff0p+48), RD); break;
        }
        _mm512_storeu_pd(out + i, r);
    }
}

static double rnd_ldexp(double m, int emax) {  // m * 2^-(random in [0, emax)); rng calls in sequence
    return std::ldexp(m, -(int)(rng() % (uint64_t)emax));
}

static double rnd_scaled(int emin, int emax) {  // random sign, mantissa and exponent in [emin, emax]
    const double m = rng_uniform(1.0, 2.0);
    const int e = emin + (int)(rng() % (uint64_t)(emax - emin + 1));
    return (rng() & 1) ? -std::ldexp(m, e) : std::ldexp(m, e);
}

// operands for each emulated operation, in the ranges where pow uses them (no underflow in the error terms;
// see the README for why pow's results do not depend on the underflowing cases)
static void primitive_operands(int op, std::vector<double> &a, std::vector<double> &b, std::vector<double> &c) {
    a.clear(); b.clear(); c.clear();
    auto add = [&](double x, double y, double z) { a.push_back(x); b.push_back(y); c.push_back(z); };
    const size_t nr = 400000;
    switch (op) {
        case 0:  // every 16-bit prefix of the mantissa (4 members each) and random mantissas
            for (uint64_t p = 0; p < 65536; p++) {
                const uint64_t lows[4] = {0, (1ULL << 36) - 1, 1ULL << 35, rng() & ((1ULL << 36) - 1)};
                for (uint64_t l : lows) add(bitsd(0x3fe0000000000000ULL | p << 36 | l), 0, 0);
            }
            for (size_t i = 0; i < nr; i++) add(bitsd(0x3fe0000000000000ULL | (rng() >> 12)), 0, 0);
            break;
        case 1:  // products: random, exact (short mantissas), near powers of two, overflow
            for (size_t i = 0; i < nr; i++) {
                const double x = rnd_scaled(-400, 400);
                switch (i % 4) {
                    case 0: add(x, rnd_scaled(-400, 400), 0); break;
                    case 1: add(x, rnd_ldexp((double)(rng() >> 40), 30), 0); break;
                    case 2: add(x, (1.0 + rng_uniform(-0x1p-40, 0x1p-40)) / x, 0); break;
                    default: { const double v = rnd_scaled(1000, 1023); add(v, rnd_scaled(-10, 20), 0); } break;
                }
            }
            break;
        case 2:  // a*b + c: cancellation, c small or large relative to a*b, c = -RN(a*b)
            for (size_t i = 0; i < nr; i++) {
                const double x = rnd_scaled(-200, 200), y = rnd_scaled(-200, 200), p = x * y;
                switch (i % 5) {
                    case 0: add(x, y, -p * (1.0 + rng_uniform(-0x1p-45, 0x1p-45))); break;
                    case 1: add(x, y, rnd_ldexp(rnd_scaled(0, 0), 120) * std::fabs(p)); break;
                    case 2: add(x, y, 0x1p120 * rnd_ldexp(rnd_scaled(0, 0), 120) * std::fabs(p)); break;
                    case 3: add(x, y, -p); break;
                    default: { const double v = rnd_ldexp((double)(rng() >> 30), 60); add(x, v, rnd_scaled(-300, 300)); } break;
                }
            }
            break;
        case 3:  // a + b: all relative scales, exact cancellation, large values. (*) The emulation does not
                 // handle an overflowing sum (inf instead of DBL_MAX); in pow that needs |P| ~ DBL_MAX, a
                 // callout lane either way
            for (size_t i = 0; i < nr; i++) {
                const double x = rnd_scaled(-500, 500);
                switch (i % 4) {
                    case 0: add(x, rnd_ldexp(rnd_scaled(0, 0), 120) * std::fabs(x), 0); break;
                    case 1: add(x, -x * (1.0 + rng_uniform(-0x1p-50, 0x1p-50)), 0); break;
                    case 2: add(x, -x, 0); break;
                    default: { const double v = rnd_scaled(1020, 1021); add(v, rnd_scaled(1020, 1021), 0); } break;  // (*)
                }
            }
            break;
        default:  // reduce / RD(s + shifter): s in [-1100, 1100], multiples of 1/16 +- ulps, tiny, (-1/16, 0)
            for (size_t i = 0; i < nr; i++) {
                switch (i % 5) {
                    case 0: add(rng_uniform(-1100, 1100), 0, 0); break;
                    case 1: { const double v = (double)((int64_t)(rng() % 35200) - 17600) / 16;
                              add(bitsd(dbits(v) + (rng() % 5) - 2), 0, 0); break; }
                    case 2: add(rnd_scaled(-1074, -4), 0, 0); break;
                    case 3: add(-rng_uniform(0, 0.0625), 0, 0); break;
                    default: add((double)((int64_t)(rng() % 35200) - 17600) / 16, 0, 0); break;
                }
            }
            add(0.0, 0, 0); add(-0.0, 0, 0); add(-0.0625, 0, 0); add(-0.03125, 0, 0); add(5e-324, 0, 0); add(-5e-324, 0, 0);
            break;
    }
    while (a.size() % 64) add(a[0], b[0], c[0]);
}
#endif

static void scalar_pow(const double *x, const double *y, double *z, size_t n) {
    for (size_t i = 0; i < n; i++) z[i] = std::pow(x[i], y[i]);
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
    g_x.resize(n_random + POW_INPUTS_FIXED);
    g_y.resize(n_random + POW_INPUTS_FIXED);
    const size_t n = make_inputs(g_x.data(), g_y.data(), n_random);
    g_x.resize(n);
    g_y.resize(n);
    g_ref.resize(n);
    for (size_t i = 0; i < n; i++) g_ref[i] = powl((long double)g_x[i], (long double)g_y[i]);
    std::vector<double> z(n);
    FILE *hf = hash_file ? fopen(hash_file, "a") : nullptr;
    if (hf) {
        std::vector<double> xy(2 * n);
        for (size_t i = 0; i < n; i++) { xy[2 * i] = g_x[i]; xy[2 * i + 1] = g_y[i]; }
        fprintf(hf, "%s pow inputs - %zu %016llx\n", kArch, n, (unsigned long long)acc_hash(xy.data(), xy.size()));
    }

    // reference: the C translation of the assembly (needs AVX-512; e.g. under Intel SDE)
    std::vector<double> ref_c;
#if defined(__x86_64__)
    if (__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq")) {
        ref_c.resize(n);
        svml_pow_ha_array(g_x.data(), g_y.data(), ref_c.data(), n);
        if (hf) fprintf(hf, "%s pow C-AVX512 svml %zu %016llx\n", kArch, n, (unsigned long long)result_hash(ref_c));
    }
#endif

    printf("%s, %zu input pairs, errors vs powl (long double with %d-bit mantissa)%s\n\n", kArch, n, LDBL_MANT_DIG,
           ref_c.empty() ? "" : "; bit-identity checked against the C translation of the assembly");
    printf("%-26s %5s %9s %10s %8s %8s %10s\n", "implementation", "lanes", "max ulp", "not CR", ">1 ulp", "special",
           do_bench ? "ns/value" : "");
    int failures = 0;
    bool have_hash = false;
    uint64_t svml_hash = 0;
    for (int64_t target : hwy::SupportedAndGeneratedTargets()) {
        hwy::SetSupportedTargetsForTest(target);
        const int variant = hwy_pow_variant();
        const char *tname = hwy::TargetName(target);
        hwy_pow_array(g_x.data(), g_y.data(), z.data(), n);
        const acc_stats s = accuracy(z);
        char name[64];
        snprintf(name, sizeof name, "hwy %s%s", tname, variant ? "" : " (std::pow)");
        print_row(name, hwy_pow_lanes(), s, do_bench ? bench(hwy_pow_array) : 0);
        if (!variant) continue;
        const uint64_t h = result_hash(z);
        if (hf) fprintf(hf, "%s pow %s svml %zu %016llx\n", kArch, tname, n, (unsigned long long)h);
        if (s.over1 || s.mismatched_special) {
            printf("  FAIL: %zu errors > 1 ulp, %zu special-value mismatches\n", s.over1, s.mismatched_special);
            failures++;
        }
        if (!have_hash) { have_hash = true; svml_hash = h; }
        else if (h != svml_hash) { printf("  FAIL: results differ from the first SVML-algorithm target\n"); failures++; }
        if (!ref_c.empty()) {
            size_t diff = 0;
            for (size_t i = 0; i < n; i++) {
                if (dbits(z[i]) != dbits(ref_c[i]) && diff++ < 5) {
                    printf("  differs from the C translation: pow(%a, %a) = %a, C %a\n", g_x[i], g_y[i], z[i], ref_c[i]);
                }
            }
            printf("  %zu of %zu results differ from the C translation of the assembly\n", diff, n);
            if (diff) failures++;
        }
    }
    hwy::SetSupportedTargetsForTest(0);

#if defined(__x86_64__)
    // the emulated AVX-512 operations against the instructions, on operands covering their use in pow
    if (!ref_c.empty()) {
        static const char *op_names[6] = {"vrndscalepd(vrcp14pd)", "vmulpd {rz-sae}", "vfmadd {rz-sae}",
                                          "vaddpd {rz-sae}", "vreducepd 0x41", "vaddpd {rd-sae}"};
        printf("\nemulated AVX-512 operations vs the instructions (operands, mismatches per target):\n");
        std::vector<double> a, b, c, r, e;
        const std::vector<int64_t> targets = hwy::SupportedAndGeneratedTargets();
        for (int op = 0; op < 6; op++) {
            primitive_operands(op, a, b, c);
            r.resize(a.size());
            e.resize(a.size());
            ref_primitive(op, a.data(), b.data(), c.data(), r.data(), a.size());
            printf("  %-22s %8zu:", op_names[op], a.size());
            for (int64_t target : targets) {
                hwy::SetSupportedTargetsForTest(target);
                if (!hwy_pow_variant()) continue;
                hwy_pow_primitive(op, a.data(), b.data(), c.data(), e.data(), a.size());
                size_t bad = 0;
                for (size_t i = 0; i < a.size(); i++) {
                    if (dbits(e[i]) != dbits(r[i]) && bad++ < 3) {
                        printf("\n    %s: op(%a, %a, %a) = %a, instruction %a", hwy::TargetName(target), a[i], b[i], c[i],
                               e[i], r[i]);
                    }
                }
                printf(" %s %zu", hwy::TargetName(target), bad);
                if (bad) failures++;
            }
            printf("\n");
        }
        hwy::SetSupportedTargetsForTest(0);
    }
#endif

    struct extra_impl ex[16];
    const size_t nex = extra_impls(ex, 16);
    for (size_t k = 0; k < nex; k++) {
        ex[k].fn(g_x.data(), g_y.data(), z.data(), n);
        print_row(ex[k].name, 0, accuracy(z), do_bench ? bench(ex[k].fn) : 0);
    }
    scalar_pow(g_x.data(), g_y.data(), z.data(), n);
    print_row("glibc pow (scalar)", 1, accuracy(z), do_bench ? bench(scalar_pow) : 0);

    if (hf) fclose(hf);
    printf("\n%s\n", failures ? "FAILED" : "all checks passed");
    return failures ? 1 : 0;
}
