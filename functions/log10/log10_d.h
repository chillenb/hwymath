/* C/intrinsics translation of Intel SVML's svml_z0_log10_d_ha.s (BSD-3-Clause, Copyright (C) 2023 Intel). */
#ifndef SVML_LOG10_D_HA_H
#define SVML_LOG10_D_HA_H

#include <stddef.h>

#if defined(__x86_64__)
#include <immintrin.h>

/* log10 of 8 doubles; needs AVX-512F and AVX-512DQ at run time. Equivalent of __svml_log108_ha. */
__attribute__((target("avx512f,avx512dq"))) __m512d svml_log108_ha(__m512d x);

/* y[i] = log10(x[i]) for i < len; needs AVX-512F/DQ at run time. */
__attribute__((target("avx512f,avx512dq"))) void svml_log10_ha_array(const double *x, double *y, size_t len);

#endif /* __x86_64__ */

/* Scalar log10 used for special lanes (equivalent of __svml_dlog10_ha_cout_rare_internal). Plain C, runs
 * anywhere. Returns 0, 1 (x < 0 or -inf: NaN, invalid) or 2 (x = +-0: -inf, divide-by-zero). */
int svml_dlog10_ha_cout_rare(const double *px, double *pr);

#endif
