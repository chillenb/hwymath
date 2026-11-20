/* Exact model of the reciprocal step of SVML's cbrt: vrndscalepd(vrcp14pd(m), imm 0x48), i.e. vrcp14pd
 * rounded to nearest to 4 fractional bits, for m in +-[1, 2). The result is k/16 with k = 8..16.
 *
 * vrcp14pd's bits are implementation-defined (documented relative error < 2^-14). On Intel's
 * implementation (Intel SDE, which models it) its output for m in [1, 2) depends only on the top 16
 * mantissa bits of m, has 16 fractional bits, and is monotone; rounded to 4 fractional bits it switches
 * from k+1 to k at the thresholds below, all multiples of 2^-16. Exact rounding of 16/m would switch at
 * 32/(2k+1) instead: these differ by up to 1.4 prefixes (2.2e-5), so the result *does* depend on the
 * approximation, and the thresholds are needed. tests/test_rcp14.c checks this model against the
 * instructions for all 2^16 prefixes (both signs, several low-bit patterns) and random inputs.
 * C and C++. */
#ifndef CBRT_RCP14_MODEL_H
#define CBRT_RCP14_MODEL_H

/* |m| >= cbrt_rcp14_threshold[i] -> k <= 15 - i (i = 0..7) */
static const double cbrt_rcp14_threshold[8] = {
    0x1.0842p+0,  /* 16 -> 15 (32/31 = 0x1.08421084...) */
    0x1.1a7dp+0,  /* 15 -> 14 (32/29 = 0x1.1a7b9611...) */
    0x1.2f69p+0,  /* 14 -> 13 (32/27 = 0x1.2f684bda...) */
    0x1.47adp+0,  /* 13 -> 12 (32/25 = 0x1.47ae147a...) */
    0x1.642cp+0,  /* 12 -> 11 (32/23 = 0x1.642c8590...) */
    0x1.8619p+0,  /* 11 -> 10 (32/21 = 0x1.86186186...) */
    0x1.af29p+0,  /* 10 ->  9 (32/19 = 0x1.af286bca...) */
    0x1.e1e2p+0,  /*  9 ->  8 (32/17 = 0x1.e1e1e1e1...) */
};

/* k such that vrndscalepd(vrcp14pd(m), 0x48) = +-k/16, for |m| in [1, 2) */
static inline int
cbrt_rcp14_k(double m)
{
    const double a = m < 0 ? -m : m;
    int k = 16;
    for (int i = 0; i < 8; i++) k -= a >= cbrt_rcp14_threshold[i];
    return k;
}

#endif
