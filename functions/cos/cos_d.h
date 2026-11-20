/* C/intrinsics translation of Intel SVML's svml_z0_cos_d_ha.s (BSD-3-Clause, Copyright (C) 2023 Intel). */
#ifndef SVML_COS_D_HA_H
#define SVML_COS_D_HA_H

#include <stddef.h>

#if defined(__x86_64__)
#include <immintrin.h>

/* cos of 8 doubles; needs AVX-512F and AVX-512DQ at run time. Equivalent of __svml_cos8_ha. */
__attribute__((target("avx512f,avx512dq"))) __m512d svml_cos8_ha(__m512d x);

/* y[i] = cos(x[i]) for i < len; needs AVX-512F and AVX-512DQ at run time. */
__attribute__((target("avx512f,avx512dq"))) void svml_cos_ha_array(const double *x, double *y, size_t len);

#endif /* __x86_64__ */

/* Scalar routine for the lanes with x = +-inf or NaN (equivalent of __svml_dcos_ha_cout_rare_internal).
 * Plain C, runs anywhere. Returns 1 for +-inf (invalid operation, result NaN), else 0; writes *pr only
 * for inf and NaN inputs. */
int svml_dcos_ha_cout_rare(const double *px, double *pr);

#endif
