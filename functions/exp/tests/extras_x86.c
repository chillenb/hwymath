/* x86-64 comparison implementations: this project's C kernels, glibc libmvec, SLEEF (HAVE_SLEEF) and
 * Intel's libsvml (HAVE_INTEL_SVML), in AVX2 and, when the CPU supports it, AVX-512 versions. */
#include <immintrin.h>
#include <string.h>

#include "../exp_d.h"
#include "extras.h"

#define AVX2 __attribute__((target("avx2,fma")))
#define AVX512 __attribute__((target("avx512f")))

AVX2 __m256d _ZGVdN4v_exp(__m256d);
AVX512 __m512d _ZGVeN8v_exp(__m512d);
#ifdef HAVE_SLEEF  /* declared here: sleef.h only declares them when compiled with -mavx2 / -mavx512f */
AVX2 __m256d Sleef_expd4_u10avx2(__m256d);
AVX512 __m512d Sleef_expd8_u10avx512f(__m512d);
#endif
#ifdef HAVE_INTEL_SVML
AVX2 __m256d __svml_exp4_ha_l9(__m256d);
AVX512 __m512d __svml_exp8_ha_z0(__m512d);
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

LOOP4(libmvec_avx2, _ZGVdN4v_exp)
LOOP8(libmvec_avx512, _ZGVeN8v_exp)
#ifdef HAVE_INTEL_SVML
LOOP4(intel_ha_avx2, __svml_exp4_ha_l9)
LOOP8(intel_ha_avx512, __svml_exp8_ha_z0)
#endif
#ifdef HAVE_SLEEF
LOOP4(sleef_avx2, Sleef_expd4_u10avx2)
LOOP8(sleef_avx512, Sleef_expd8_u10avx512f)
#endif

size_t
extra_impls(struct extra_impl *out, size_t max)
{
    size_t k = 0;
#define ADD(NAME, FN) if (k < max) { out[k].name = NAME; out[k].fn = FN; k++; }
    const int avx2 = __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
    const int avx512 = __builtin_cpu_supports("avx512f");
    if (avx2) {
        ADD("C hp AVX2", svml_exp_hp_array_avx2)
        ADD("C fast AVX2", svml_exp_fast_array_avx2)
        ADD("C SVML port AVX2", svml_exp_ha_array_avx2)
#ifdef HAVE_INTEL_SVML
        ADD("Intel SVML HA AVX2", intel_ha_avx2)
#endif
#ifdef HAVE_SLEEF
        ADD("SLEEF u10 AVX2", sleef_avx2)
#endif
        ADD("libmvec AVX2", libmvec_avx2)
    }
    if (avx512) {
        ADD("C SVML AVX-512", svml_exp_ha_array)
#ifdef HAVE_INTEL_SVML
        ADD("Intel SVML HA AVX-512", intel_ha_avx512)
#endif
#ifdef HAVE_SLEEF
        ADD("SLEEF u10 AVX-512", sleef_avx512)
#endif
        ADD("libmvec AVX-512", libmvec_avx512)
    }
    return k;
#undef ADD
}
