/* C/intrinsics translation of Intel SVML's svml_z0_tan_d_ha.s (BSD-3-Clause, Copyright (C) 2023 Intel). */
#ifndef SVML_TAN_D_HA_H
#define SVML_TAN_D_HA_H

#include <stddef.h>

#if defined(__x86_64__)
#include <immintrin.h>

/* tan of 8 doubles; needs AVX-512F and AVX-512DQ at run time. Equivalent of __svml_tan8_ha. */
__attribute__((target("avx512f,avx512dq"))) __m512d svml_tan8_ha(__m512d x);

/* evaluation part of svml_tan8_ha (exported for make tan-insns): tan(j*pi/16 + R + Rlo), j = idx mod 16 */
__attribute__((target("avx512f,avx512dq"))) __m512d svml_tan8_ha_eval(__m512i idx, __m512d r, __m512d rlo);

/* y[i] = tan(x[i]) for i < len; needs AVX-512F and AVX-512DQ at run time. */
__attribute__((target("avx512f,avx512dq"))) void svml_tan_ha_array(const double *x, double *y, size_t len);

#endif /* __x86_64__ */

/* Scalar routine for the lanes with x = +-inf or NaN (equivalent of __svml_dtan_ha_cout_rare_internal).
 * Plain C, runs anywhere. Returns 1 for +-inf (invalid operation, result NaN), else 0; writes *pr only
 * for inf and NaN inputs. */
int svml_dtan_ha_cout_rare(const double *px, double *pr);

#endif
