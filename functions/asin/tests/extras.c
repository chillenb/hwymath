/* Comparison implementations of asin (see extras.h). x86-64 only; none elsewhere. */
#include <string.h>

#include "extras.h"

#if defined(__x86_64__)
#include <immintrin.h>

#define AVX2 __attribute__((target("avx2,fma")))
#define AVX512 __attribute__((target("avx512f")))

#ifdef HAVE_SLEEF  /* declared here: sleef.h only declares them when compiled with -mavx2 / -mavx512f */
AVX2 __m256d Sleef_asind4_u10avx2(__m256d);
AVX512 __m512d Sleef_asind8_u10avx512f(__m512d);
#endif
#ifdef HAVE_INTEL_SVML
AVX2 __m256d __svml_asin4_ha_l9(__m256d);
AVX512 __m512d __svml_asin8_ha_z0(__m512d);
#endif

/* array loops: full vectors, then the tail through a zero-padded buffer */
#define LOOP4(NAME, CALL)                                                                   \
    AVX2 static void NAME(const double *x, double *y, size_t n) {                           \
        size_t i = 0;                                                                       \
        for (; i + 4 <= n; i += 4) _mm256_storeu_pd(y + i, CALL(_mm256_loadu_pd(x + i)));   \
        if (i < n) {                                                                        \
            double b[4] = {0, 0, 0, 0}, o[4];                                               \
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
            double b[8] = {0, 0, 0, 0, 0, 0, 0, 0}, o[8];                                   \
            memcpy(b, x + i, (n - i) * 8);                                                  \
            _mm512_storeu_pd(o, CALL(_mm512_loadu_pd(b)));                                  \
            memcpy(y + i, o, (n - i) * 8);                                                  \
        }                                                                                   \
    }

#ifdef HAVE_INTEL_SVML
LOOP4(intel_ha_avx2, __svml_asin4_ha_l9)
LOOP8(intel_ha_avx512, __svml_asin8_ha_z0)
#endif
#ifdef HAVE_SLEEF
LOOP4(sleef_avx2, Sleef_asind4_u10avx2)
LOOP8(sleef_avx512, Sleef_asind8_u10avx512f)
#endif

size_t
extra_impls(struct extra_impl *out, size_t max)
{
    size_t k = 0;
#define ADD(NAME, FN) if (k < max) { out[k].name = NAME; out[k].fn = FN; k++; }
    const int avx2 = __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
    const int avx512 = __builtin_cpu_supports("avx512f");
    (void)avx2; (void)avx512;
    if (avx2) {
#ifdef HAVE_INTEL_SVML
        ADD("Intel SVML HA AVX2", intel_ha_avx2)
#endif
#ifdef HAVE_SLEEF
        ADD("SLEEF u10 AVX2", sleef_avx2)
#endif
    }
    if (avx512) {
#ifdef HAVE_INTEL_SVML
        ADD("Intel SVML HA AVX-512", intel_ha_avx512)
#endif
#ifdef HAVE_SLEEF
        ADD("SLEEF u10 AVX-512", sleef_avx512)
#endif
    }
    return k;
#undef ADD
}

#else

size_t
extra_impls(struct extra_impl *out, size_t max)
{
    (void)out; (void)max;
    return 0;
}

#endif
