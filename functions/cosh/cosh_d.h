/* C/intrinsics translation of Intel SVML's svml_z0_cosh_d_ha.s (BSD-3-Clause, Copyright (C) 2023 Intel). */
#ifndef SVML_COSH_D_HA_H
#define SVML_COSH_D_HA_H

#include <stddef.h>

#if defined(__x86_64__)
#include <immintrin.h>

/* cosh of 8 doubles; needs AVX-512F/DQ at run time. Equivalent of __svml_cosh8_ha. */
__attribute__((target("avx512f,avx512dq"))) __m512d svml_cosh8_ha(__m512d x);

/* y[i] = cosh(x[i]) for i < len; needs AVX-512F/DQ at run time. */
__attribute__((target("avx512f,avx512dq"))) void svml_cosh_ha_array(const double *x, double *y, size_t len);

#endif /* __x86_64__ */

/* Scalar cosh used for special lanes (equivalent of __svml_dcosh_ha_cout_rare_internal).
 * Plain C, runs anywhere. Returns 0, or 3 (overflow). */
int svml_dcosh_ha_cout_rare(const double *px, double *pr);

#endif
