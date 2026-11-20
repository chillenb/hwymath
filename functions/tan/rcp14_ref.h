/* Scalar reference of VRCP14PD (approximate reciprocal, relative error < 2^-14), following Intel's
 * reference implementation RCP14D in recip14.c ("Reference Implementations for IA Approximation
 * Instructions VRCP14, VRSQRT14, VRCP28, VRSQRT28, and VEXP2", Intel, 2015; BSD-3-Clause,
 * Copyright (c) 2015 Intel Corporation), for MXCSR.DAZ = MXCSR.FTZ = 0 (the defaults; SVML runs with them).
 * Rewritten with integer arithmetic; the coefficient table is Intel's RCP14_Coeff. Plain C.
 *
 * For 1 < m < 2 (m = significand, not a power of 2): i = top 6 fraction bits, l = next 10 bits;
 *   rcp14(m) = truncate_to_17_bits(a[i] - b[i]*(m17 - y[i])) / 2^18, m17 = m truncated to 16 fraction
 *   bits, y[i] = 1 + 1/128 + i/64, so that (a[i] - b[i]*(l - 512)/256) / 2^18 lies in [1/2, 1) and its
 *   truncation to 17 significant bits is floor((256*a[i] - b[i]*(l - 512)) / 512) * 2^-17.
 * Other inputs: exponent scaling (the result's biased exponent is 2045 - e), powers of 2 exactly,
 * results below the normal range denormalized by truncation, denormal inputs normalized first (and
 * |x| <= 2^-1024 gives +-inf), +-0 -> +-inf, +-inf -> +-0, NaN -> the quieted NaN. */
#ifndef HWYMATH_RCP14_REF_H
#define HWYMATH_RCP14_REF_H

#include <stdint.h>

#ifdef __cplusplus
#define RCP14_CONST static constexpr
#define RCP14_INLINE static constexpr inline
#else
#define RCP14_CONST static const __attribute__((unused))
#define RCP14_INLINE static inline
#endif

/* RCP14_Coeff from recip14.c: b[i] (10 bits) and a[i] (18 bits) alternate */
RCP14_CONST uint32_t rcp14_coeff[128] = {
    1009, 260119, 977, 256148, 949, 252296, 921, 248558, 893, 244929, 869, 241405, 843, 237981, 821, 234652,
    797,  231416, 777, 228266, 755, 225202, 735, 222220, 717, 219314, 699, 216485, 681, 213727, 663, 211038,
    647,  208417, 631, 205859, 617, 203364, 601, 200929, 587, 198551, 573, 196229, 561, 193960, 547, 191743,
    535,  189576, 523, 187458, 513, 185387, 501, 183360, 491, 181377, 479, 179439, 469, 177540, 459, 175681,
    451,  173860, 441, 172077, 433, 170330, 423, 168618, 415, 166940, 407, 165295, 399, 163682, 391, 162101,
    385,  160550, 377, 159027, 369, 157535, 363, 156069, 357, 154631, 349, 153219, 343, 151832, 337, 150470,
    331,  149133, 325, 147819, 319, 146528, 315, 145260, 309, 144012, 303, 142787, 299, 141582, 293, 140397,
    289,  139232, 285, 138085, 279, 136959, 275, 135853, 271, 134763, 267, 133689, 263, 132631, 259, 131589,
};

/* 17-bit significand q (2^16 <= q < 2^17) of rcp14(1.f) for a fraction f (52 bits, f != 0) */
RCP14_INLINE uint64_t
rcp14_sig17(uint64_t f)
{
    const uint32_t i = (uint32_t)(f >> 46), l = (uint32_t)(f >> 36) & 1023;
    const int64_t a = rcp14_coeff[2 * i + 1], b = rcp14_coeff[2 * i];
    return (uint64_t)((256 * a - b * ((int64_t)l - 512)) >> 9);
}

static inline uint64_t
rcp14_ref_bits(uint64_t x)
{
    const uint64_t sign = x & 0x8000000000000000ULL;
    uint64_t e = (x >> 52) & 0x7ff, f = x & 0x000fffffffffffffULL;
    if (e == 0x7ff) return f ? x | 0x0008000000000000ULL : sign;  /* NaN: quieted; inf: 0 */
    if (e == 0) {
        if (f == 0) return sign | 0x7ff0000000000000ULL;            /* +-0: +-inf */
        if (f <= 0x0004000000000000ULL) return sign | 0x7ff0000000000000ULL;  /* reciprocal overflows */
        /* normalize: 1.f' * 2^(-1022 - s); the reciprocal is normal */
        int s = 0;
        while (!(f & 0x0010000000000000ULL)) { f <<= 1; s++; }
        f &= 0x000fffffffffffffULL;
        uint64_t be = 0x7fc + (uint64_t)s;  /* = 0x7fd - (1 - s) */
        if (f == 0) return sign | ((be + 1) << 52);  /* power of 2: exact */
        return sign | (be << 52) | ((rcp14_sig17(f) - 0x10000) << 36);
    }
    if (f == 0) {  /* power of 2: exact */
        if (e == 0x7fe) return sign | 0x0008000000000000ULL;      /* 2^-1023, denormal */
        return sign | ((0x7fe - e) << 52);
    }
    const uint64_t q = rcp14_sig17(f) << 36;  /* 53-bit significand, 2^52 <= q < 2^53 */
    if (e == 0x7fe) return sign | (q >> 2);                       /* denormal result */
    if (e == 0x7fd) return sign | (q >> 1);                       /* denormal result */
    return sign | ((0x7fd - e) << 52) | (q & 0x000fffffffffffffULL);
}

#endif
