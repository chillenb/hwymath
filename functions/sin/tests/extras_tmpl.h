/* Other implementations of sin/cos/tan available on this machine, as array functions, for comparison
 * with the Highway versions (hwy_test). Included by functions/{sin,cos,tan}/tests/extras.c, which define
 *   TRIG_F          the function name, e.g. sin
 *   TRIG_SVML_ARRAY the array version of the C/AVX-512 translation, e.g. svml_sin_ha_array
 *   TRIG_NO_LIBMVEC (optional) if libmvec has no vector version (tan: none in the glibc 2.28 sysroot)
 * x86-64: the C/AVX-512 translation (if the CPU has AVX-512F/DQ), Intel's libsvml (HAVE_INTEL_SVML;
 * __svml_<f>4_ha_l9 and __svml_<f>8_ha_z0, only called), SLEEF u10 (HAVE_SLEEF), glibc libmvec, each in
 * AVX2 and AVX-512 versions as supported by the CPU. AArch64: SLEEF u10 Advanced SIMD. */
#include <string.h>

#include "functions/sin/tests/extras_tmpl_decl.h"

#define TRIG_CAT_(a, b) a##b
#define TRIG_CAT(a, b) TRIG_CAT_(a, b)
#define TRIG_CAT3(a, b, c) TRIG_CAT(TRIG_CAT(a, b), c)

#define ADD(NAME, FN)                                                                                 \
    if (k < max) {                                                                                    \
        out[k].name = NAME;                                                                           \
        out[k].fn = FN;                                                                               \
        k++;                                                                                          \
    }

#if defined(__x86_64__)
#include <immintrin.h>

#define AVX2 __attribute__((target("avx2,fma")))
#define AVX512 __attribute__((target("avx512f,avx512dq")))

#define LIBMVEC4 TRIG_CAT(_ZGVdN4v_, TRIG_F)
#define LIBMVEC8 TRIG_CAT(_ZGVeN8v_, TRIG_F)
#define SLEEF4 TRIG_CAT3(Sleef_, TRIG_F, d4_u10avx2)
#define SLEEF8 TRIG_CAT3(Sleef_, TRIG_F, d8_u10avx512f)
#define INTEL4 TRIG_CAT3(__svml_, TRIG_F, 4_ha_l9)
#define INTEL8 TRIG_CAT3(__svml_, TRIG_F, 8_ha_z0)

#ifndef TRIG_NO_LIBMVEC
AVX2 __m256d LIBMVEC4(__m256d);
AVX512 __m512d LIBMVEC8(__m512d);
#endif
#ifdef HAVE_SLEEF /* declared here: sleef.h only declares them when compiled with -mavx2 / -mavx512f */
AVX2 __m256d SLEEF4(__m256d);
AVX512 __m512d SLEEF8(__m512d);
#endif
#ifdef HAVE_INTEL_SVML
AVX2 __m256d INTEL4(__m256d);
AVX512 __m512d INTEL8(__m512d);
#endif

/* array loops: full vectors, then the tail through a zero-padded buffer */
#define LOOP4(NAME, CALL)                                                                             \
    AVX2 static void NAME(const double *x, double *y, size_t n) {                                     \
        size_t i = 0;                                                                                 \
        for (; i + 4 <= n; i += 4) _mm256_storeu_pd(y + i, CALL(_mm256_loadu_pd(x + i)));             \
        if (i < n) {                                                                                  \
            double b[4] = {0, 0, 0, 0}, o[4];                                                         \
            memcpy(b, x + i, (n - i) * 8);                                                            \
            _mm256_storeu_pd(o, CALL(_mm256_loadu_pd(b)));                                            \
            memcpy(y + i, o, (n - i) * 8);                                                            \
        }                                                                                             \
    }
#define LOOP8(NAME, CALL)                                                                             \
    AVX512 static void NAME(const double *x, double *y, size_t n) {                                   \
        size_t i = 0;                                                                                 \
        for (; i + 8 <= n; i += 8) _mm512_storeu_pd(y + i, CALL(_mm512_loadu_pd(x + i)));             \
        if (i < n) {                                                                                  \
            double b[8] = {0, 0, 0, 0, 0, 0, 0, 0}, o[8];                                             \
            memcpy(b, x + i, (n - i) * 8);                                                            \
            _mm512_storeu_pd(o, CALL(_mm512_loadu_pd(b)));                                            \
            memcpy(y + i, o, (n - i) * 8);                                                            \
        }                                                                                             \
    }

#ifndef TRIG_NO_LIBMVEC
LOOP4(libmvec_avx2, LIBMVEC4)
LOOP8(libmvec_avx512, LIBMVEC8)
#endif
#ifdef HAVE_INTEL_SVML
LOOP4(intel_ha_avx2, INTEL4)
LOOP8(intel_ha_avx512, INTEL8)
#endif
#ifdef HAVE_SLEEF
LOOP4(sleef_avx2, SLEEF4)
LOOP8(sleef_avx512, SLEEF8)
#endif

size_t
extra_impls(struct extra_impl *out, size_t max)
{
    size_t k = 0;
    const int avx2 = __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
    const int avx512 = __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq");
    if (avx2) {
#ifdef HAVE_INTEL_SVML
        ADD("Intel SVML HA AVX2", intel_ha_avx2)
#endif
#ifdef HAVE_SLEEF
        ADD("SLEEF u10 AVX2", sleef_avx2)
#endif
#ifndef TRIG_NO_LIBMVEC
        ADD("libmvec AVX2", libmvec_avx2)
#endif
    }
    if (avx512) {
        ADD("C SVML translation AVX-512", TRIG_SVML_ARRAY)
#ifdef HAVE_INTEL_SVML
        ADD("Intel SVML HA AVX-512", intel_ha_avx512)
#endif
#ifdef HAVE_SLEEF
        ADD("SLEEF u10 AVX-512", sleef_avx512)
#endif
#ifndef TRIG_NO_LIBMVEC
        ADD("libmvec AVX-512", libmvec_avx512)
#endif
    }
    return k;
}

#elif defined(__aarch64__)

#ifdef HAVE_SLEEF
#include <arm_neon.h>
#include <sleef.h>

#define SLEEF2 TRIG_CAT3(Sleef_, TRIG_F, d2_u10advsimd)

static void
sleef_advsimd(const double *x, double *y, size_t n)
{
    size_t i = 0;
    for (; i + 2 <= n; i += 2) vst1q_f64(y + i, SLEEF2(vld1q_f64(x + i)));
    for (; i < n; i++) {
        double b[2] = {x[i], 0}, o[2];
        vst1q_f64(o, SLEEF2(vld1q_f64(b)));
        y[i] = o[0];
    }
}
#endif

size_t
extra_impls(struct extra_impl *out, size_t max)
{
    size_t k = 0;
#ifdef HAVE_SLEEF
    ADD("SLEEF u10 AdvSIMD", sleef_advsimd)
#else
    (void)out;
    (void)max;
#endif
    return k;
}

#else

size_t
extra_impls(struct extra_impl *out, size_t max)
{
    (void)out;
    (void)max;
    return 0;
}

#endif
#undef ADD
