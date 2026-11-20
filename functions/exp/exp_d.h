/* C/intrinsics translation of Intel SVML's svml_z0_exp_d_ha.s (BSD-3-Clause, Copyright (C) 2023 Intel). */
#ifndef SVML_EXP_D_HA_H
#define SVML_EXP_D_HA_H

#include <stddef.h>

#if defined(__x86_64__)
#include <immintrin.h>

/* exp of 8 doubles; needs AVX-512F at run time. Equivalent of __svml_exp8_ha. */
__attribute__((target("avx512f"))) __m512d svml_exp8_ha(__m512d x);

/* y[i] = exp(x[i]) for i < len; needs AVX-512F at run time. */
__attribute__((target("avx512f"))) void svml_exp_ha_array(const double *x, double *y, size_t len);

/* AVX2/FMA3 port of the same algorithm (avx2/svml_exp_d_ha_avx2.c): exp of 4 doubles, and array version. */
__attribute__((target("avx2,fma"))) __m256d svml_exp4_ha_avx2(__m256d x);
__attribute__((target("avx2,fma"))) void svml_exp_ha_array_avx2(const double *x, double *y, size_t len);

/* Faster AVX2 variant (avx2/exp_d_avx2_fast.c): 8-entry tables, round-to-nearest reduction, own
 * degree-7 minimax polynomial; not bit-identical to the SVML algorithm. */
__attribute__((target("avx2,fma"))) __m256d svml_exp4_fast_avx2(__m256d x);
__attribute__((target("avx2,fma"))) void svml_exp_fast_array_avx2(const double *x, double *y, size_t len);

/* Extra-precision AVX2 variant (avx2/exp_d_avx2_hp.c): like the fast variant, with an exact
 * two-part reduction and an error-free reconstruction, aiming at ~0.5 ulp. */
__attribute__((target("avx2,fma"))) __m256d svml_exp4_hp_avx2(__m256d x);
__attribute__((target("avx2,fma"))) void svml_exp_hp_array_avx2(const double *x, double *y, size_t len);

#endif /* __x86_64__ */

/* Scalar exp used for special lanes (equivalent of __svml_dexp_ha_cout_rare_internal).
 * Plain C, runs anywhere. Returns 0, 3 (overflow) or 4 (underflow). */
int svml_dexp_ha_cout_rare(const double *px, double *pr);

#endif
