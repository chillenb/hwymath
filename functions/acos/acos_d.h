/* C/intrinsics translation of Intel SVML's svml_z0_acos_d_ha.s (BSD-3-Clause, Copyright (C) 2023 Intel). */
#ifndef SVML_ACOS_D_HA_H
#define SVML_ACOS_D_HA_H

#include <stddef.h>

#if defined(__x86_64__)
#include <immintrin.h>

/* acos of 8 doubles; needs AVX-512F at run time. Equivalent of __svml_acos8_ha. */
__attribute__((target("avx512f"))) __m512d svml_acos8_ha(__m512d x);

/* y[i] = acos(x[i]) for i < len; needs AVX-512F at run time. */
__attribute__((target("avx512f"))) void svml_acos_ha_array(const double *x, double *y, size_t len);

#endif /* __x86_64__ */

/* Scalar acos used for special lanes (equivalent of __svml_dacos_ha_cout_rare_internal), translated
 * mechanically from the assembly. Plain C, runs anywhere. Returns 0, or 1 (domain error: |x| > 1 or inf). */
int svml_dacos_ha_cout_rare(const double *px, double *pr);

#endif
