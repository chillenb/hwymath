/* C/intrinsics translation of Intel SVML's svml_z0_tanh_d_ha.s (BSD-3-Clause, Copyright (C) 2023 Intel). */
#ifndef SVML_TANH_D_HA_H
#define SVML_TANH_D_HA_H

#include <stddef.h>

#if defined(__x86_64__)
#include <immintrin.h>

/* tanh of 8 doubles; needs AVX-512F/DQ at run time. Equivalent of __svml_tanh8_ha. */
__attribute__((target("avx512f,avx512dq"))) __m512d svml_tanh8_ha(__m512d x);

/* y[i] = tanh(x[i]) for i < len; needs AVX-512F/DQ at run time. */
__attribute__((target("avx512f,avx512dq"))) void svml_tanh_ha_array(const double *x, double *y, size_t len);

#endif /* __x86_64__ */

/* Scalar tanh used for special lanes (equivalent of __svml_dtanh_ha_cout_rare_internal).
 * Plain C, runs anywhere. Always returns 0. */
int svml_dtanh_ha_cout_rare(const double *px, double *pr);

#endif
