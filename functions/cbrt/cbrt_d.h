/* C/intrinsics translation of Intel SVML's svml_z0_cbrt_d_ha.s (BSD-3-Clause, Copyright (C) 2023 Intel). */
#ifndef SVML_CBRT_D_HA_H
#define SVML_CBRT_D_HA_H

#include <stddef.h>

#if defined(__x86_64__)
#include <immintrin.h>

/* cbrt of 8 doubles; needs AVX-512F at run time. Equivalent of __svml_cbrt8_ha. */
__attribute__((target("avx512f"))) __m512d svml_cbrt8_ha(__m512d x);

/* y[i] = cbrt(x[i]) for i < len; needs AVX-512F at run time. */
__attribute__((target("avx512f"))) void svml_cbrt_ha_array(const double *x, double *y, size_t len);

#endif /* __x86_64__ */

/* Scalar cbrt (equivalent of __svml_dcbrt_ha_cout_rare_internal). The vector routine never calls it
 * (it handles 0, inf, NaN and denormals itself); it is VML's scalar fallback. Plain C, runs anywhere.
 * Always returns 0. */
int svml_dcbrt_ha_cout_rare(const double *px, double *pr);

#endif
