/* C/intrinsics translation of Intel SVML's svml_z0_atanh_d_ha.s (BSD-3-Clause, Copyright (C) 2023 Intel). */
#ifndef SVML_ATANH_D_HA_H
#define SVML_ATANH_D_HA_H

#include <stddef.h>

#if defined(__x86_64__)
#include <immintrin.h>

/* atanh of 8 doubles; needs AVX-512F at run time. Equivalent of __svml_atanh8_ha. */
__attribute__((target("avx512f"))) __m512d svml_atanh8_ha(__m512d x);

/* y[i] = atanh(x[i]) for i < len; needs AVX-512F at run time. */
__attribute__((target("avx512f"))) void svml_atanh_ha_array(const double *x, double *y, size_t len);
#endif

/* Scalar routine for the lanes with |x| >= 1 or NaN (equivalent of __svml_datanh_ha_cout_rare_internal):
 * x*x (quiet NaN) for NaN (code 0), +-inf with the divide-by-zero flag for +-1 (code 2), and NaN with the
 * invalid flag for everything else (code 1). Plain C, runs anywhere. */
int svml_datanh_ha_cout_rare(const double *px, double *pr);

#endif
