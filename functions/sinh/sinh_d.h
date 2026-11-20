/* C/intrinsics translation of Intel SVML's svml_z0_sinh_d_ha.s (BSD-3-Clause, Copyright (C) 2023 Intel). */
#ifndef SVML_SINH_D_HA_H
#define SVML_SINH_D_HA_H

#include <stddef.h>

#if defined(__x86_64__)
#include <immintrin.h>

/* sinh of 8 doubles; needs AVX-512F/DQ at run time. Equivalent of __svml_sinh8_ha. */
__attribute__((target("avx512f,avx512dq"))) __m512d svml_sinh8_ha(__m512d x);

/* y[i] = sinh(x[i]) for i < len; needs AVX-512F/DQ at run time. */
__attribute__((target("avx512f,avx512dq"))) void svml_sinh_ha_array(const double *x, double *y, size_t len);

#endif /* __x86_64__ */

/* Scalar sinh used for special lanes (equivalent of __svml_dsinh_ha_cout_rare_internal).
 * Plain C, runs anywhere. Returns 0, or 3 (overflow). */
int svml_dsinh_ha_cout_rare(const double *px, double *pr);

#endif
