/* C/intrinsics translation of Intel SVML's svml_z0_pow_d_ha.s (BSD-3-Clause, Copyright (C) 2023 Intel). */
#ifndef SVML_POW_D_HA_H
#define SVML_POW_D_HA_H

#include <stddef.h>

#if defined(__x86_64__)
#include <immintrin.h>

/* pow of 8 pairs of doubles; needs AVX-512F at run time. Equivalent of __svml_pow8_ha. */
__attribute__((target("avx512f,avx512dq"))) __m512d svml_pow8_ha(__m512d x, __m512d y);

/* z[i] = pow(x[i], y[i]) for i < len; needs AVX-512F at run time. */
__attribute__((target("avx512f,avx512dq"))) void svml_pow_ha_array(const double *x, const double *y, double *z,
                                                                   size_t len);
#endif /* __x86_64__ */

/* Scalar pow used for special lanes (equivalent of __svml_dpow_ha_cout_rare_internal). Plain C, runs
 * anywhere. Writes *pr; returns 0, or 1 for a domain error (x < 0 with non-integer y) or a pole
 * (x = 0 with y < 0). */
int svml_dpow_ha_cout_rare(const double *px, const double *py, double *pr);

#endif
