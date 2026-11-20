/*
 * Large-argument path shared by the C/AVX-512 translations of Intel SVML's high-accuracy double sin, cos
 * and tan (svml_z0_{sin,cos,tan}_d_ha.s): the three assemblies contain the same Payne-Hanek-style
 * reduction (instruction for instruction, up to register allocation and scheduling), and sin and cos also
 * share its evaluation tail. Included by svml_sin_d_ha.c, svml_cos_d_ha.c and svml_tan_d_ha.c.
 *
 * Original: Copyright (C) 2023 Intel Corporation, SPDX-License-Identifier: BSD-3-Clause.
 * This translation is distributed under the same license.
 *
 * Reduction (all lanes; the callers merge the lanes with |x| above their threshold):
 *   xs = |x| (sin) or x (cos, tan); e = biased exponent, M = 2^52 + mantissa (integer, < 2^53).
 *   T0, T1, T2 = trig_red_Tab[3e .. 3e+2]: 192 bits of 2/pi (in units of 2pi turns, i.e. of 1/(2pi))
 *   aligned for exponent e, so that the 128 bits hi:lo of (M * T) >> 64 (a truncated 53x192-bit product
 *   computed from ten 32x32->64-bit products) hold the fraction of x/(2pi) as a 128-bit fixed-point
 *   number (the integer part, a multiple of 2pi, has been dropped by the table).
 *   hi (bits 127..76 with sign and exponent of 1.0) becomes the double F = +-1.f, which is rounded to a
 *   multiple of 2^-9 by adding 1.5*2^43: N = that multiple mod 512 (x = N*2pi/512 + r), and the
 *   remainder F - N/512 (exact) plus the medium and low 52/24-bit pieces of the fraction, converted to
 *   doubles by the usual "or exponent bits, subtract the same power of 2" trick, give RHi + RLo.
 *   (RHi + RLo) * 2pi = P + low, with an FMA for the exact error of RHi*2piHi.
 * Inputs with |x| <= 2^-20 pass through unreduced (unreachable here: callers use it for |x| > 2^16).
 */
#ifndef SVML_TRIG_REDUCTION_H
#define SVML_TRIG_REDUCTION_H

#include <immintrin.h>

#include "svml_trig_reduction_data.inc"

#define TRIG_RN (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)
#define TRIG_V1(c) _mm512_set1_pd(c)
#define TRIG_U1(c) _mm512_set1_epi64((long long)(c))
#define TRIG_AS_D(v) _mm512_castsi512_pd(v)
#define TRIG_AS_I(v) _mm512_castpd_si512(v)

/* Payne-Hanek reduction of xs: x = N*2pi/512 + P + low (|P| <= pi/512 roughly), N in [0, 511].
 * nonfinite: lanes where x is inf or NaN (they get the scalar callout). */
__attribute__((target("avx512f,avx512dq"), always_inline)) static inline void
trig_reduce_avx512(__m512d xs, __m512d *p_hi, __m512d *p_lo, __m512i *p_n, __mmask8 *nonfinite)
{
    const __m512i xb = TRIG_AS_I(xs);
    const __m512i lo32 = TRIG_U1(trig_red_Low32);

    /* inf/NaN lanes */
    const __m512i ebits = _mm512_and_si512(xb, TRIG_U1(trig_red_ExpMask));
    const __mmask8 finite = _mm512_cmp_round_pd_mask(TRIG_AS_D(ebits), TRIG_AS_D(TRIG_U1(trig_red_ExpMask)),
                                                     _CMP_NEQ_UQ, _MM_FROUND_NO_EXC);
    const __m512i nanmask = _mm512_mask_andnot_epi64(TRIG_U1(trig_red_AllOnes), finite, ebits, ebits);
    *nonfinite = _mm512_cmp_round_pd_mask(TRIG_AS_D(nanmask), TRIG_AS_D(nanmask), _CMP_UNORD_Q, _MM_FROUND_NO_EXC);

    /* the three table words for this exponent (gathers, as in the assembly) */
    const __m512i e = _mm512_srli_epi64(_mm512_and_si512(xb, TRIG_U1(trig_red_ExpMask)), 52);
    const __m256i off = _mm512_cvtepi64_epi32(_mm512_slli_epi64(_mm512_add_epi64(e, _mm512_slli_epi64(e, 1)), 3));
    const __m512i t0 = TRIG_AS_I(_mm512_i32gather_pd(off, (const void *)trig_red_Tab, 1));
    const __m512i t1 = TRIG_AS_I(_mm512_i32gather_pd(off, (const void *)(trig_red_Tab + 1), 1));
    const __m512i t2 = TRIG_AS_I(_mm512_i32gather_pd(off, (const void *)(trig_red_Tab + 2), 1));
    const __m512i t0h = _mm512_srli_epi64(t0, 32), t1h = _mm512_srli_epi64(t1, 32), t2h = _mm512_srli_epi64(t2, 32);
    const __m512i t0l = _mm512_and_si512(t0, lo32), t1l = _mm512_and_si512(t1, lo32), t2l = _mm512_and_si512(t2, lo32);

    /* significand with the implicit bit (wrong for denormals, as noted in the assembly; harmless here) */
    const __m512i m = _mm512_add_epi64(_mm512_and_si512(xb, TRIG_U1(trig_red_MantMask)), TRIG_U1(trig_red_ImplicitBit));
    const __m512i mh = _mm512_srli_epi64(m, 32), ml = _mm512_and_si512(m, lo32);

    /* the big multiplication and carry propagation (vpmullq of 32-bit values: exact products) */
    const __m512i mh_t1h = _mm512_mullo_epi64(mh, t1h);
    const __m512i mh_t0l = _mm512_mullo_epi64(mh, t0l);
    const __m512i mh_t1l = _mm512_mullo_epi64(mh, t1l);
    const __m512i mh_t2h = _mm512_mullo_epi64(mh, t2h);
    const __m512i mh_t2l = _mm512_mullo_epi64(mh, t2l);
    const __m512i ml_t0l = _mm512_mullo_epi64(ml, t0l);
    const __m512i ml_t1l = _mm512_mullo_epi64(ml, t1l);
    const __m512i ml_t2h = _mm512_mullo_epi64(ml, t2h);
    const __m512i ml_t0h = _mm512_mullo_epi64(ml, t0h);
    const __m512i ml_t1h = _mm512_mullo_epi64(ml, t1h);
    const __m512i s1 = _mm512_add_epi64(mh_t2h, _mm512_srli_epi64(ml_t2h, 32));
    const __m512i s2 = _mm512_add_epi64(mh_t1l, _mm512_srli_epi64(ml_t1l, 32));
    const __m512i s3 = _mm512_add_epi64(mh_t0l, _mm512_srli_epi64(ml_t0l, 32));
    const __m512i s4 = _mm512_add_epi64(mh_t1h, _mm512_srli_epi64(ml_t1h, 32));
    const __m512i s5 = _mm512_add_epi64(_mm512_and_si512(ml_t0h, lo32), s3);
    const __m512i s6 = _mm512_add_epi64(_mm512_srli_epi64(mh_t2l, 32), _mm512_add_epi64(_mm512_and_si512(ml_t1l, lo32), s1));
    const __m512i s7 = _mm512_add_epi64(_mm512_srli_epi64(s6, 32), _mm512_add_epi64(_mm512_and_si512(ml_t1h, lo32), s2));
    const __m512i s8 = _mm512_add_epi64(_mm512_srli_epi64(s7, 32), _mm512_add_epi64(s4, _mm512_and_si512(ml_t0l, lo32)));
    const __m512i s9 = _mm512_add_epi64(s5, _mm512_srli_epi64(s8, 32));
    const __m512i wlo = _mm512_add_epi64(_mm512_slli_epi64(s7, 32), _mm512_and_si512(s6, lo32));  /* bits 63..0 */
    const __m512i whi = _mm512_add_epi64(_mm512_slli_epi64(s9, 32), _mm512_and_si512(s8, lo32));  /* bits 127..64 */

    /* sign of xs into all pieces */
    const __m512i sign = _mm512_and_si512(xb, TRIG_U1(trig_red_SignMask));

    /* high part +-1.f (bits 127..76), rounded to a multiple of 2^-9: N, and the exact remainder */
    const __m512d fhi = TRIG_AS_D(_mm512_or_si512(_mm512_xor_si512(sign, TRIG_U1(trig_red_OneBits)), _mm512_srli_epi64(whi, 12)));
    const __m512d sh = _mm512_add_round_pd(TRIG_V1(trig_red_Shifter), fhi, TRIG_RN);
    __m512d rem = _mm512_sub_round_pd(fhi, _mm512_sub_round_pd(sh, TRIG_V1(trig_red_Shifter), TRIG_RN), TRIG_RN);

    /* low piece (bits 23..0 of wlo, times 2^-104... as 1.xxx*2^-52 minus 2^-52) and medium piece
     * (bits 75..24: low 12 bits of whi, high 40 of wlo) */
    const __m512i two104 = _mm512_xor_si512(sign, TRIG_U1(trig_red_Two104Bits));
    const __m512i two52 = _mm512_xor_si512(sign, TRIG_U1(trig_red_Two52Bits));
    const __m512d plo = _mm512_sub_round_pd(
            TRIG_AS_D(_mm512_or_si512(_mm512_slli_epi64(_mm512_and_si512(wlo, TRIG_U1(trig_red_Low24)), 28), two104)),
            TRIG_AS_D(two104), TRIG_RN);
    const __m512d pmed = _mm512_sub_round_pd(
            TRIG_AS_D(_mm512_or_si512(_mm512_or_si512(_mm512_slli_epi64(_mm512_and_si512(whi, TRIG_U1(trig_red_Low12)), 40),
                                                       _mm512_srli_epi64(wlo, 24)), two52)),
            TRIG_AS_D(two52), TRIG_RN);

    /* two aligned pieces RHi + RLo, then times 2pi */
    const __m512d rhi = _mm512_add_round_pd(rem, pmed, TRIG_RN);
    const __m512d phi = _mm512_mul_round_pd(rhi, TRIG_V1(trig_red_TwoPiHi), TRIG_RN);
    rem = _mm512_sub_round_pd(rem, rhi, TRIG_RN);
    const __m512d perr = _mm512_fmsub_round_pd(rhi, TRIG_V1(trig_red_TwoPiHi), phi, TRIG_RN);
    const __m512d rlo = _mm512_add_round_pd(_mm512_add_round_pd(pmed, rem, TRIG_RN), plo, TRIG_RN);
    const __m512d t = _mm512_fmadd_round_pd(TRIG_V1(trig_red_TwoPiLo), rhi, perr, TRIG_RN);
    *p_lo = _mm512_fmadd_round_pd(TRIG_V1(trig_red_TwoPiHi), rlo, t, TRIG_RN);
    *p_hi = phi;
    *p_n = _mm512_and_si512(TRIG_AS_I(sh), TRIG_U1(trig_red_Mask9));
}

/* Reduced argument r + rlo of the large path: P + low, or xs itself (and 0) if |x| <= 2^-20. */
__attribute__((target("avx512f,avx512dq"), always_inline)) static inline void
trig_reduced_avx512(__m512d xs, __m512d phi, __m512d plo, __m512d *r, __m512d *rlo)
{
    const __m512d ax = TRIG_AS_D(_mm512_and_si512(TRIG_AS_I(xs), TRIG_U1(trig_red_AbsMask)));
    const __mmask8 big = _mm512_cmp_round_pd_mask(ax, TRIG_V1(trig_red_Tiny), _CMP_NLE_UQ, _MM_FROUND_NO_EXC);
    const __mmask8 small = _mm512_cmp_round_pd_mask(ax, TRIG_V1(trig_red_Tiny), _CMP_NGT_UQ, _MM_FROUND_NO_EXC);
    const __m512i ones = TRIG_U1(trig_red_AllOnes);
    const __m512i keep_x = _mm512_mask_andnot_epi64(ones, big, TRIG_AS_I(ax), TRIG_AS_I(ax));    /* ~0 if tiny */
    const __m512i keep_p = _mm512_mask_andnot_epi64(ones, small, TRIG_AS_I(ax), TRIG_AS_I(ax));  /* ~0 if not tiny */
    *r = TRIG_AS_D(_mm512_or_si512(_mm512_and_si512(TRIG_AS_I(phi), keep_p), _mm512_and_si512(keep_x, TRIG_AS_I(xs))));
    *rlo = TRIG_AS_D(_mm512_and_si512(TRIG_AS_I(plo), keep_p));
}

/*
 * Evaluation tail of the sin and cos large path: with B = N*2pi/512 and the 4-entry table row
 * {T0, T1, T2, T3} = T[N] (sin: T0 + T1 = cos B with T1 = 2^round(log2 cos B), T2 + T3 = sin B;
 * cos: the same for cos B = sin(B + pi/2)), f(B + r) = T2 + (T0 + T1)*r + T2*(cos r - 1) + (T0 + T1)*(sin r - r)
 * + T3 + rlo*(T0 + T1 - T2*r), evaluated with an exact split of T2 + T1*r + T0*r.
 */
__attribute__((target("avx512f,avx512dq"), always_inline)) static inline __m512d
trig_sincos_large_tail_avx512(__m512d r, __m512d rlo, __m512i n, const double *tab, double s1, double s2,
                              double c1, double c2, double c3)
{
    const __m256i off = _mm512_cvtepi64_epi32(_mm512_slli_epi64(n, 5));
    const __m512d t0 = _mm512_i32gather_pd(off, (const void *)tab, 1);
    const __m512d t1 = _mm512_i32gather_pd(off, (const void *)(tab + 1), 1);
    const __m512d t2 = _mm512_i32gather_pd(off, (const void *)(tab + 2), 1);
    const __m512d r2 = _mm512_mul_round_pd(r, r, TRIG_RN);
    const __m512d a = _mm512_fmadd_round_pd(t1, r, t2, TRIG_RN);          /* T2 + T1*r */
    const __m512d hi = _mm512_fmadd_round_pd(r, t0, a, TRIG_RN);          /* T2 + T1*r + T0*r */
    __m512d e1 = _mm512_sub_round_pd(t2, a, TRIG_RN);
    __m512d e2 = _mm512_sub_round_pd(a, hi, TRIG_RN);
    e1 = _mm512_fmadd_round_pd(t1, r, e1, TRIG_RN);                       /* error of a */
    e2 = _mm512_fmadd_round_pd(r, t0, e2, TRIG_RN);                       /* error of hi */
    const __m512d sp = _mm512_fmadd_round_pd(TRIG_V1(s2), r2, TRIG_V1(s1), TRIG_RN);
    const __m512d err = _mm512_add_round_pd(e2, e1, TRIG_RN);
    const __m512d ct = _mm512_add_round_pd(t0, t1, TRIG_RN);              /* T0 + T1 */
    const __m512d sp3 = _mm512_mul_round_pd(_mm512_mul_round_pd(sp, r2, TRIG_RN), r, TRIG_RN);  /* sin r - r */
    const __m512d corr = _mm512_fnmadd_round_pd(t2, r, ct, TRIG_RN);     /* T0 + T1 - T2*r */
    __m512d lo = _mm512_fmadd_round_pd(ct, sp3, err, TRIG_RN);
    __m512d cp = _mm512_fmadd_round_pd(TRIG_V1(c3), r2, TRIG_V1(c2), TRIG_RN);
    cp = _mm512_fmadd_round_pd(r2, cp, TRIG_V1(c1), TRIG_RN);
    cp = _mm512_mul_round_pd(cp, r2, TRIG_RN);                            /* cos r - 1 */
    lo = _mm512_fmadd_round_pd(t2, cp, lo, TRIG_RN);
    const __m512d t3 = _mm512_i32gather_pd(off, (const void *)(tab + 3), 1);
    const __m512d lo2 = _mm512_fmadd_round_pd(corr, rlo, t3, TRIG_RN);
    return _mm512_add_round_pd(hi, _mm512_add_round_pd(lo, lo2, TRIG_RN), TRIG_RN);
}

#endif
