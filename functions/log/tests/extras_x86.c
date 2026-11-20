/* x86-64 comparison implementations of FN (log, log2, log10 or log1p; set with -DFN=...): this project's
 * C/AVX-512 translation, glibc libmvec (log only in glibc 2.28), SLEEF (HAVE_SLEEF) and Intel's libsvml
 * (HAVE_INTEL_SVML), in AVX2 and, when the CPU supports it, AVX-512 versions. */
#include <immintrin.h>
#include <string.h>

#include "extras.h"

#define CAT_(a, b, c) a##b##c
#define CAT(a, b, c) CAT_(a, b, c)
#define STR_(a) #a
#define STR(a) STR_(a)
#define AVX2 __attribute__((target("avx2,fma")))
#define AVX512 __attribute__((target("avx512f,avx512dq")))

#define LOG_IS_log 1
#define IS_LOG CAT(LOG_IS_, FN, )

AVX512 void CAT(svml_, FN, _ha_array)(const double *x, double *y, size_t len);  /* the C translation */
#if IS_LOG == 1
AVX2 __m256d _ZGVdN4v_log(__m256d);
AVX512 __m512d _ZGVeN8v_log(__m512d);
#endif
#ifdef HAVE_SLEEF  /* declared here: sleef.h only declares them when compiled with -mavx2 / -mavx512f */
AVX2 __m256d CAT(Sleef_, FN, d4_u10avx2)(__m256d);
AVX512 __m512d CAT(Sleef_, FN, d8_u10avx512f)(__m512d);
#endif
#ifdef HAVE_INTEL_SVML
AVX2 __m256d CAT(__svml_, FN, 4_ha_l9)(__m256d);
AVX512 __m512d CAT(__svml_, FN, 8_ha_z0)(__m512d);
#endif

/* array loops: full vectors, then the tail through a buffer padded with 1.0 */
#define LOOP4(NAME, CALL)                                                                   \
    AVX2 static void NAME(const double *x, double *y, size_t n) {                           \
        size_t i = 0;                                                                       \
        for (; i + 4 <= n; i += 4) _mm256_storeu_pd(y + i, CALL(_mm256_loadu_pd(x + i)));   \
        if (i < n) {                                                                        \
            double b[4] = {1, 1, 1, 1}, o[4];                                               \
            memcpy(b, x + i, (n - i) * 8);                                                  \
            _mm256_storeu_pd(o, CALL(_mm256_loadu_pd(b)));                                  \
            memcpy(y + i, o, (n - i) * 8);                                                  \
        }                                                                                   \
    }
#define LOOP8(NAME, CALL)                                                                   \
    AVX512 static void NAME(const double *x, double *y, size_t n) {                         \
        size_t i = 0;                                                                       \
        for (; i + 8 <= n; i += 8) _mm512_storeu_pd(y + i, CALL(_mm512_loadu_pd(x + i)));   \
        if (i < n) {                                                                        \
            double b[8] = {1, 1, 1, 1, 1, 1, 1, 1}, o[8];                                   \
            memcpy(b, x + i, (n - i) * 8);                                                  \
            _mm512_storeu_pd(o, CALL(_mm512_loadu_pd(b)));                                  \
            memcpy(y + i, o, (n - i) * 8);                                                  \
        }                                                                                   \
    }

#if IS_LOG == 1
LOOP4(libmvec_avx2, _ZGVdN4v_log)
LOOP8(libmvec_avx512, _ZGVeN8v_log)
#endif
#ifdef HAVE_INTEL_SVML
LOOP4(intel_ha_avx2, CAT(__svml_, FN, 4_ha_l9))
LOOP8(intel_ha_avx512, CAT(__svml_, FN, 8_ha_z0))
#endif
#ifdef HAVE_SLEEF
LOOP4(sleef_avx2, CAT(Sleef_, FN, d4_u10avx2))
LOOP8(sleef_avx512, CAT(Sleef_, FN, d8_u10avx512f))
#endif

size_t
extra_impls(struct extra_impl *out, size_t max)
{
    size_t k = 0;
#define ADD(NAME, F) if (k < max) { out[k].name = NAME; out[k].fn = F; k++; }
    const int avx2 = __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
    const int avx512 = __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq");
    if (avx2) {
#ifdef HAVE_INTEL_SVML
        ADD("Intel SVML HA AVX2", intel_ha_avx2)
#endif
#ifdef HAVE_SLEEF
        ADD("SLEEF u10 AVX2", sleef_avx2)
#endif
#if IS_LOG == 1
        ADD("libmvec AVX2", libmvec_avx2)
#endif
    }
    if (avx512) {
        ADD("C SVML AVX-512", CAT(svml_, FN, _ha_array))
#ifdef HAVE_INTEL_SVML
        ADD("Intel SVML HA AVX-512", intel_ha_avx512)
#endif
#ifdef HAVE_SLEEF
        ADD("SLEEF u10 AVX-512", sleef_avx512)
#endif
#if IS_LOG == 1
        ADD("libmvec AVX-512", libmvec_avx512)
#endif
    }
    (void)avx2;
    return k;
#undef ADD
}
