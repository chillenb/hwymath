/* C/intrinsics translation of Intel SVML's svml_z0_log1p_d_ha.s (BSD-3-Clause, Copyright (C) 2023 Intel). */
#ifndef SVML_LOG1P_D_HA_H
#define SVML_LOG1P_D_HA_H

#include <stddef.h>

#if defined(__x86_64__)
#include <immintrin.h>

/* log1p of 8 doubles; needs AVX-512F at run time. Equivalent of __svml_log1p8_ha. */
__attribute__((target("avx512f"))) __m512d svml_log1p8_ha(__m512d x);

/* y[i] = log1p(x[i]) for i < len; needs AVX-512F at run time. */
__attribute__((target("avx512f"))) void svml_log1p_ha_array(const double *x, double *y, size_t len);

#endif /* __x86_64__ */

/* Scalar log1p used for special lanes (equivalent of __svml_dlog1p_ha_cout_rare_internal). Plain C, runs
 * anywhere. Returns 0, 1 (x < -1 or -inf: NaN, invalid) or 2 (x = -1: -inf, divide-by-zero). */
int svml_dlog1p_ha_cout_rare(const double *px, double *pr);

#endif
