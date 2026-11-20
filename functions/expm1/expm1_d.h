/* C/intrinsics translation of Intel SVML's svml_z0_expm1_d_ha.s (BSD-3-Clause, Copyright (C) 2023 Intel). */
#ifndef SVML_EXPM1_D_HA_H
#define SVML_EXPM1_D_HA_H

#include <stddef.h>

#if defined(__x86_64__)
#include <immintrin.h>

/* expm1 of 8 doubles; needs AVX-512F and AVX-512DQ at run time. Equivalent of __svml_expm18_ha. */
__attribute__((target("avx512f,avx512dq"))) __m512d svml_expm18_ha(__m512d x);

/* y[i] = expm1(x[i]) for i < len; needs AVX-512F and AVX-512DQ at run time. */
__attribute__((target("avx512f,avx512dq"))) void svml_expm1_ha_array(const double *x, double *y, size_t len);

#endif /* __x86_64__ */

/* Scalar expm1 used for special lanes (equivalent of __svml_dexpm1_ha_cout_rare_internal).
 * Plain C, runs anywhere. Returns 0, 3 (overflow) or 4 (underflow). */
int svml_dexpm1_ha_cout_rare(const double *px, double *pr);

#endif
