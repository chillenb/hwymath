/* C/intrinsics translation of Intel SVML's svml_z0_asinh_d_ha.s (BSD-3-Clause, Copyright (C) 2023 Intel). */
#ifndef SVML_ASINH_D_HA_H
#define SVML_ASINH_D_HA_H

#include <stddef.h>

#if defined(__x86_64__)
#include <immintrin.h>

/* asinh of 8 doubles; needs AVX-512F at run time. Equivalent of __svml_asinh8_ha. */
__attribute__((target("avx512f"))) __m512d svml_asinh8_ha(__m512d x);

/* y[i] = asinh(x[i]) for i < len; needs AVX-512F at run time. */
__attribute__((target("avx512f"))) void svml_asinh_ha_array(const double *x, double *y, size_t len);
#endif

/* Scalar routine for the lanes with |x| > DBL_MAX or NaN (equivalent of __svml_dasinh_ha_cout_rare_internal):
 * returns x for +-inf, x*x (quiet NaN) for NaN; return code 0. Plain C, runs anywhere. */
int svml_dasinh_ha_cout_rare(const double *px, double *pr);

#endif
