/* C/intrinsics translation of Intel SVML's svml_z0_atan2_d_ha.s (BSD-3-Clause, Copyright (C) 2023 Intel). */
#ifndef SVML_ATAN2_D_HA_H
#define SVML_ATAN2_D_HA_H

#include <stddef.h>

#if defined(__x86_64__)
#include <immintrin.h>

/* atan2(y, x) of 8 pairs; needs AVX-512F at run time. Equivalent of __svml_atan28_ha (y in zmm0, x in zmm1). */
__attribute__((target("avx512f"))) __m512d svml_atan28_ha(__m512d y, __m512d x);

/* r[i] = atan2(y[i], x[i]) for i < len; needs AVX-512F at run time. */
__attribute__((target("avx512f"))) void svml_atan2_ha_array(const double *y, const double *x, double *r, size_t len);

/* The same with vrcp14pd (implementation-defined bits) replaced by the portable approximation of the
 * Highway version, trunc16(RN(1/B)): reference for the bit-identity of the Highway targets. */
__attribute__((target("avx512f"))) __m512d svml_atan28_ha_portable(__m512d y, __m512d x);
__attribute__((target("avx512f"))) void svml_atan2_ha_portable_array(const double *y, const double *x, double *r,
                                                                      size_t len);

#endif /* __x86_64__ */

/* Scalar atan2 used for special lanes (equivalent of __svml_datan2_ha_cout_rare_internal), translated
 * mechanically from the assembly. Plain C, runs anywhere. Returns 0. */
int svml_datan2_ha_cout_rare(const double *py, const double *px, double *pr);

/* Result of the vector routine's branch for zero arguments (x or y = +-0, neither NaN): portable C,
 * bit-identical to the vector code. */
double svml_atan2_ha_zero_args(double y, double x);

#endif
