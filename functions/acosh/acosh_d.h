/* C/intrinsics translation of Intel SVML's svml_z0_acosh_d_ha.s (BSD-3-Clause, Copyright (C) 2023 Intel). */
#ifndef SVML_ACOSH_D_HA_H
#define SVML_ACOSH_D_HA_H

#include <stddef.h>

#if defined(__x86_64__)
#include <immintrin.h>

/* acosh of 8 doubles; needs AVX-512F at run time. Equivalent of __svml_acosh8_ha. */
__attribute__((target("avx512f"))) __m512d svml_acosh8_ha(__m512d x);

/* y[i] = acosh(x[i]) for i < len; needs AVX-512F at run time. */
__attribute__((target("avx512f"))) void svml_acosh_ha_array(const double *x, double *y, size_t len);
#endif

/* Scalar routine for the lanes with x <= 1, x > DBL_MAX or NaN (equivalent of __svml_dacosh_ha_cout_rare_internal):
 * returns x*x (quiet NaN) for NaN, +inf for +inf and +0 for 1 (code 0), and NaN with the invalid flag for
 * everything else (code 1). Plain C, runs anywhere. */
int svml_dacosh_ha_cout_rare(const double *px, double *pr);

#endif
