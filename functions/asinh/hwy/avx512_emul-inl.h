// Exact Highway emulations of the AVX-512 instructions used by Intel SVML's asinh/acosh/atanh (and other
// log-based routines): VRCP14PD, VRSQRT14PD, VGETEXPPD (positive normal inputs) and the 16-entry
// VPERMT2PD/VPERMI2PD table lookup. Shared by functions/{asinh,acosh,atanh}/hwy; per-target header
// (re-included by foreach_target.h). Needs HWY_NATIVE_FMA for the exactness of NegMulAdd below.
//
// VRCP14PD/VRSQRT14PD are "implementation-defined" approximations in the SDM, but Intel publishes a
// bit-exact reference (RECIP14.c, see tools/gen_approx14.py), which Intel SDE follows and which
// tests/test_approx14 checks against the hardware instructions. It evaluates a linear function on one of
// 64 segments of the mantissa, r = a - 256*b*(xt - y), with 17- or 18-bit integer coefficients a, b. Here
// the coefficients are computed (one division, a square root for rsqrt14) plus a 4-bit correction read
// from 64-bit words (approx14_tables.h), because a 64-entry lookup would need gathers.

// include guard for the static data (once), then the per-target toggle
#include "functions/asinh/hwy/approx14_tables.h"

#if defined(HWYMATH_AVX512_EMUL_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef HWYMATH_AVX512_EMUL_INL_H_
#undef HWYMATH_AVX512_EMUL_INL_H_
#else
#define HWYMATH_AVX512_EMUL_INL_H_
#endif

#include <stddef.h>
#include <stdint.h>

#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace hwymath {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

// small non-negative integers (< 2^52) to double, exactly: 2^52 + v has v in its low mantissa bits
template <class D, class VU = hn::VFromD<hn::RebindToUnsigned<D>>>
HWY_INLINE hn::VFromD<D> SmallU64ToF64(D d, VU v) {
    const hn::RebindToUnsigned<D> du;
    return hn::Sub(hn::BitCast(d, hn::Or(v, hn::Set(du, 0x4330000000000000ULL))), hn::Set(d, 0x1p52));
}

// 4-bit field k (0..63) of the four 64-bit words w[0..3], as a double
template <class D, class VU = hn::VFromD<hn::RebindToUnsigned<D>>>
HWY_INLINE hn::VFromD<D> Nibble(D d, VU k, const uint64_t *w) {
    const hn::RebindToUnsigned<D> du;
    const auto b4 = hn::TestBit(k, hn::Set(du, 16));
    const auto w01 = hn::IfThenElse(b4, hn::Set(du, w[1]), hn::Set(du, w[0]));
    const auto w23 = hn::IfThenElse(b4, hn::Set(du, w[3]), hn::Set(du, w[2]));
    const auto word = hn::IfThenElse(hn::TestBit(k, hn::Set(du, 32)), w23, w01);
    const auto sh = hn::ShiftLeft<2>(hn::And(k, hn::Set(du, 15)));
    return SmallU64ToF64(d, hn::And(hn::Shr(word, sh), hn::Set(du, 15)));
}

// VRCP14PD for positive normal x whose reciprocal is normal (2^-1022 < 1/x), bit-exact.
template <class D>
HWY_INLINE hn::VFromD<D> Rcp14(D d, hn::VFromD<D> x) {
#if HWY_TARGET <= HWY_AVX3  // the instruction itself (this function emulates it bit-exactly)
    (void)d;
    return hn::ApproximateReciprocal(x);
#else
    const hn::RebindToUnsigned<D> du;
    const auto u = hn::BitCast(du, x);
    const auto m = hn::And(u, hn::Set(du, 0x000fffffffffffffULL));
    const auto e = hn::ShiftRight<52>(u);
    const auto one = hn::Set(du, 0x3ff0000000000000ULL);
    // segment i = top 6 mantissa bits; its centre y = 1 + (2i + 1)/128; x truncated to 16 fraction bits
    const auto y = hn::BitCast(d, hn::Or(hn::Or(hn::And(m, hn::Set(du, 0x000fc00000000000ULL)),
                                                hn::Set(du, 1ULL << 45)), one));
    const auto xt = hn::BitCast(d, hn::Or(hn::And(m, hn::Set(du, 0x000ffff000000000ULL)), one));
    const auto q = hn::Div(hn::Set(d, 1.0), y);
    // a = floor(2^18/y) + correction, 256*b = 256*(2*floor(512/y^2) + 1)
    const auto a = hn::Add(hn::Floor(hn::Mul(q, hn::Set(d, 262144.0))),
                           Nibble(d, hn::ShiftRight<46>(m), kRcp14Res));
    const auto b256 = hn::MulAdd(hn::Floor(hn::Mul(hn::Mul(q, q), hn::Set(d, 512.0))), hn::Set(d, 512.0),
                                 hn::Set(d, 256.0));
    const auto r = hn::NegMulAdd(b256, hn::Sub(xt, y), a);  // exact: a - 256*b*(xt - y), in (2^17, 2^18)
    auto rb = hn::And(hn::BitCast(du, r), hn::Set(du, 0xfffffff000000000ULL));  // truncated to 17 bits
    // powers of two: the exact reciprocal (r = 2^18)
    rb = hn::IfThenElse(hn::Eq(m, hn::Zero(du)), hn::Set(du, 1041ULL << 52), rb);
    // r*2^-18 * 2^(1023 - e): add (1005 - e) to the exponent field (mod 2^64)
    return hn::BitCast(d, hn::Add(rb, hn::ShiftLeft<52>(hn::Sub(hn::Set(du, 1005), e))));
#endif
}

// VRSQRT14PD for x = +-0 (+-inf), +inf (+0), positive normal x (bit-exact) and NaN (x itself, not
// quieted). Negative and denormal inputs are not emulated.
template <class D>
HWY_INLINE hn::VFromD<D> Rsqrt14(D d, hn::VFromD<D> x) {
#if HWY_TARGET <= HWY_AVX3  // the instruction itself (this function emulates it bit-exactly)
    (void)d;
    return hn::ApproximateReciprocalSqrt(x);
#else
    const hn::RebindToUnsigned<D> du;
    const hn::RebindToSigned<D> di;
    const auto u = hn::BitCast(du, x);
    const auto m = hn::And(u, hn::Set(du, 0x000fffffffffffffULL));
    const auto e = hn::And(hn::ShiftRight<52>(u), hn::Set(du, 0x7ff));
    const auto one = hn::Set(du, 0x3ff0000000000000ULL);
    // odd unbiased exponent: the mantissa is scaled to [2, 4), else to [1, 2)
    const auto hi = hn::Eq(hn::And(e, hn::Set(du, 1)), hn::Zero(du));
    const auto hid = hn::RebindMask(d, hi);
    // segment i = top 5 mantissa bits, centre y1 = 1 + (2i + 1)/64 (y = 2*y1 on [2, 4)), x truncated to
    // 15 fraction bits; on [2, 4) 128*d*(xt - y) = 256*d*(xt1 - y1)
    const auto y1b = hn::Or(hn::Or(hn::And(m, hn::Set(du, 0x000f800000000000ULL)), hn::Set(du, 1ULL << 46)), one);
    const auto y1 = hn::BitCast(d, y1b);
    const auto y = hn::BitCast(d, hn::Add(y1b, hn::IfThenElseZero(hi, hn::Set(du, 1ULL << 52))));
    const auto xt1 = hn::BitCast(d, hn::Or(hn::And(m, hn::Set(du, 0x000fffe000000000ULL)), one));
    const auto q = hn::Div(hn::Set(d, 1.0), hn::Sqrt(y));
    const auto q2 = hn::Mul(q, q);
    const auto q3 = hn::Mul(q2, q);
    const auto q5 = hn::Mul(q3, q2);
    const auto alpha = hn::IfThenElse(hid, hn::Set(d, kRsqrt14Alpha[1]), hn::Set(d, kRsqrt14Alpha[0]));
    const auto k = hn::Or(hn::ShiftRight<47>(m), hn::IfThenElseZero(hi, hn::Set(du, 32)));
    // c = floor(2^19/sqrt(y) + alpha*y^-2.5 + beta) + correction, 256*d = 256*(2*floor(K/y^1.5) + 1)
    const auto c = hn::Add(hn::Floor(hn::MulAdd(q, hn::Set(d, 524288.0),
                                                hn::MulAdd(alpha, q5, hn::Set(d, kRsqrt14Beta)))),
                           Nibble(d, k, kRsqrt14Res));
    const auto K = hn::IfThenElse(hid, hn::Set(d, 1024.0), hn::Set(d, 512.0));
    const auto d256 = hn::MulAdd(hn::Floor(hn::Mul(q3, K)), hn::Set(d, 512.0), hn::Set(d, 256.0));
    const auto r = hn::NegMulAdd(d256, hn::Sub(xt1, y1), c);  // exact, in [2^18, 2^19)
    auto rb = hn::And(hn::BitCast(du, r), hn::Set(du, 0xfffffff000000000ULL));
    // x = 4^n exactly: 2^-n (r = 2^19)
    rb = hn::IfThenElse(hn::AndNot(hi, hn::Eq(m, hn::Zero(du))), hn::Set(du, 1042ULL << 52), rb);
    // r*2^-19 * 2^-n, n = floor((e - 1023)/2)
    const auto n = hn::ShiftRight<1>(hn::Sub(hn::BitCast(di, e), hn::Set(di, 1023)));
    auto res = hn::BitCast(d, hn::Sub(rb, hn::ShiftLeft<52>(hn::BitCast(du, hn::Add(n, hn::Set(di, 19))))));
    // +-0 -> +-inf, +inf -> +0, NaN -> x
    const auto abs = hn::And(u, hn::Set(du, 0x7fffffffffffffffULL));
    res = hn::IfThenElse(hn::RebindMask(d, hn::Eq(abs, hn::Zero(du))),
                         hn::BitCast(d, hn::Or(hn::Xor(u, abs), hn::Set(du, 0x7ff0000000000000ULL))), res);
    res = hn::IfThenElse(hn::RebindMask(d, hn::Eq(u, hn::Set(du, 0x7ff0000000000000ULL))), hn::Zero(d), res);
    return hn::IfThenElse(hn::IsNaN(x), x, res);
#endif
}

// VGETEXPPD for positive normal x: the unbiased exponent as a double
template <class D>
HWY_INLINE hn::VFromD<D> GetExpPositiveNormal(D d, hn::VFromD<D> x) {
    const hn::RebindToUnsigned<D> du;
    return hn::Sub(hn::BitCast(d, hn::Or(hn::ShiftRight<52>(hn::BitCast(du, x)), hn::Set(du, 0x4330000000000000ULL))),
                   hn::Set(d, 0x1p52 + 1023.0));
}

// VPERMT2PD with a 16-entry table: tab[idx & 15] (the instruction uses the low 4 index bits)
template <class D, class VU = hn::VFromD<hn::RebindToUnsigned<D>>>
HWY_INLINE hn::VFromD<D> Lookup16(D d, const double *tab, VU idx) {
    const hn::RebindToUnsigned<D> du;
    const hn::RebindToSigned<D> di;
    const size_t N = hn::Lanes(d);
    idx = hn::And(idx, hn::Set(du, 15));
    if (N >= 16 && N % 16 == 0) {  // the whole table in one vector
        return hn::TableLookupLanes(hn::LoadN(d, tab, 16), hn::IndicesFromVec(d, hn::BitCast(di, idx)));
    }
    if (N >= 2 && 16 % (2 * N) == 0) {  // 16/(2N) two-vector lookups (AVX-512: one vpermt2pd), then selects
        const auto local = hn::IndicesFromVec(d, hn::BitCast(di, hn::And(idx, hn::Set(du, 2 * N - 1))));
        auto r = hn::TwoTablesLookupLanes(d, hn::LoadU(d, tab), hn::LoadU(d, tab + N), local);
        for (size_t base = 2 * N; base < 16; base += 2 * N) {
            const auto t = hn::TwoTablesLookupLanes(d, hn::LoadU(d, tab + base), hn::LoadU(d, tab + base + N), local);
            r = hn::IfThenElse(hn::RebindMask(d, hn::Ge(idx, hn::Set(du, base))), t, r);
        }
        return r;
    }
    return hn::GatherIndex(d, tab, hn::BitCast(di, idx));  // single lanes, odd SVE vector lengths
}

}  // namespace HWY_NAMESPACE
}  // namespace hwymath
HWY_AFTER_NAMESPACE();

#endif  // HWYMATH_AVX512_EMUL_INL_H_
