/* C/intrinsics translation of Intel SVML's svml_z0_exp2_d_ha.s (BSD-3-Clause, Copyright (C) 2023 Intel). */
#ifndef SVML_EXP2_D_HA_H
#define SVML_EXP2_D_HA_H

#include <stddef.h>

#if defined(__x86_64__)
#include <immintrin.h>

/* exp2 of 8 doubles; needs AVX-512F and AVX-512DQ at run time. Equivalent of __svml_exp28_ha. */
__attribute__((target("avx512f,avx512dq"))) __m512d svml_exp28_ha(__m512d x);

/* y[i] = exp2(x[i]) for i < len; needs AVX-512F and AVX-512DQ at run time. */
__attribute__((target("avx512f,avx512dq"))) void svml_exp2_ha_array(const double *x, double *y, size_t len);

#endif /* __x86_64__ */

/* Scalar exp2 used for special lanes (equivalent of __svml_dexp2_ha_cout_rare_internal).
 * Plain C, runs anywhere. Returns 0, 3 (overflow) or 4 (underflow). */
int svml_dexp2_ha_cout_rare(const double *px, double *pr);

#endif
