/* C/intrinsics translation of Intel SVML's svml_z0_log_d_ha.s (BSD-3-Clause, Copyright (C) 2023 Intel). */
#ifndef SVML_LOG_D_HA_H
#define SVML_LOG_D_HA_H

#include <stddef.h>

#if defined(__x86_64__)
#include <immintrin.h>

/* log of 8 doubles; needs AVX-512F and AVX-512DQ at run time. Equivalent of __svml_log8_ha. */
__attribute__((target("avx512f,avx512dq"))) __m512d svml_log8_ha(__m512d x);

/* y[i] = log(x[i]) for i < len; needs AVX-512F/DQ at run time. */
__attribute__((target("avx512f,avx512dq"))) void svml_log_ha_array(const double *x, double *y, size_t len);

#endif /* __x86_64__ */

/* Scalar log used for special lanes (equivalent of __svml_dlog_ha_cout_rare_internal). Plain C, runs
 * anywhere. Returns 0, 1 (x < 0 or -inf: NaN, invalid) or 2 (x = +-0: -inf, divide-by-zero). */
int svml_dlog_ha_cout_rare(const double *px, double *pr);

#endif
