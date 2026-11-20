// Exact emulation of AVX-512's vrcp14pd (approximate reciprocal, relative error < 2^-14) with Highway,
// following Intel's reference algorithm (functions/tan/rcp14_ref.h; MXCSR.DAZ = FTZ = 0): for a
// significand 1.f (f != 0), with i = the top 6 fraction bits and l = the next 10,
//   rcp14 = floor((256*a[i] + 512*b[i] - b[i]*l) / 512) * 2^-17   (a 17-bit significand in [1/2, 1)),
// the exponent mirrored (biased 2045 - e), exact powers of 2, denormal inputs scaled by 2^52 first,
// denormal results by shifting, +-0 and |x| <= 2^-1024 -> +-inf, +-inf -> +-0, NaN -> quiet NaN.
// The 64 coefficient pairs {a, b} (18 + 10 bits) are packed into 32-bit words and looked up with
// permutes (u32 lanes: 2 x vpermt2d + 1 select on AVX-512, 4 x (2 vpermd) + 3 selects on AVX2).

#ifndef HWYMATH_RCP14_TABLES_H
#define HWYMATH_RCP14_TABLES_H
#include <cstdint>

#include "functions/tan/rcp14_ref.h"

namespace hwymath_rcp14 {
struct alignas(64) PackedCoeff {
    uint32_t v[64];  // a[i] << 10 | b[i]
    constexpr PackedCoeff() : v{} {
        for (int i = 0; i < 64; i++) v[i] = rcp14_coeff[2 * i + 1] << 10 | rcp14_coeff[2 * i];
    }
};
constexpr PackedCoeff kCoeff{};
}  // namespace hwymath_rcp14
#endif

#if defined(HWYMATH_RCP14_HWY_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef HWYMATH_RCP14_HWY_INL_H_
#undef HWYMATH_RCP14_HWY_INL_H_
#else
#define HWYMATH_RCP14_HWY_INL_H_
#endif

#include "hwy/highway.h"
#include "functions/sin/hwy/trig_hwy-inl.h"

HWY_BEFORE_NAMESPACE();
namespace hwymath_rcp14 {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64

template <class D>
HWY_INLINE hn::VFromD<D> Rcp14(D d, hn::VFromD<D> x) {
    const hn::RebindToUnsigned<D> du;
    const hn::RebindToSigned<D> di;
    const hn::Repartition<uint32_t, decltype(du)> d32;
    using VU = hn::VFromD<decltype(du)>;
    const VU mant = hn::Set(du, 0x000fffffffffffffULL);
    const VU xb = hn::BitCast(du, x);
    const VU sign = hn::And(xb, hn::Set(du, 0x8000000000000000ULL));
    const VU ax = hn::And(xb, hn::Set(du, 0x7fffffffffffffffULL));

    // denormals: scale by 2^52 (exact), then correct the exponent
    const auto denorm = hn::Lt(ax, hn::Set(du, 0x0010000000000000ULL));
    const VU xn = hn::BitCast(du, hn::IfThenElse(hn::RebindMask(d, denorm), hn::Mul(x, hn::Set(d, 0x1p52)), x));
    const VU f = hn::And(xn, mant);
    const auto e = hn::BitCast(di, hn::ShiftRight<52>(hn::And(xn, hn::Set(du, 0x7ff0000000000000ULL))));
    const auto eadj = hn::Sub(e, hn::IfThenElseZero(hn::RebindMask(di, denorm), hn::Set(di, 52)));

    // 17-bit significand from the coefficient table
    const VU i = hn::ShiftRight<46>(f);
    const VU l = hn::And(hn::ShiftRight<36>(f), hn::Set(du, 1023));
    const auto w32 = hwymath_trig::HWY_NAMESPACE::Lookup<64>(d32, kCoeff.v, kCoeff.v, hn::BitCast(d32, i));
    const VU w = hn::And(hn::BitCast(du, w32), hn::Set(du, 0xffffffffULL));
    const VU a = hn::ShiftRight<10>(w), b = hn::And(w, hn::Set(du, 1023));
    const VU bl = hn::MulEven(hn::BitCast(d32, b), hn::BitCast(d32, l));
    const VU q = hn::ShiftRight<9>(hn::Sub(hn::Add(hn::ShiftLeft<8>(a), hn::ShiftLeft<9>(b)), bl));

    // significand with the implicit bit and biased exponent of the result (powers of 2: exact)
    const auto pow2 = hn::Eq(f, hn::Zero(du));
    const VU sig = hn::IfThenElse(pow2, hn::Set(du, 0x0010000000000000ULL), hn::ShiftLeft<36>(q));
    const auto rexp = hn::Sub(hn::IfThenElse(hn::RebindMask(di, pow2), hn::Set(di, 0x7fe), hn::Set(di, 0x7fd)), eadj);
    const auto normal = hn::Gt(rexp, hn::Zero(di));
    const VU res_n = hn::Or(sign, hn::Or(hn::ShiftLeft<52>(hn::BitCast(du, rexp)), hn::And(sig, mant)));
    const VU shift = hn::BitCast(du, hn::IfThenElseZero(hn::Not(normal), hn::Sub(hn::Set(di, 1), rexp)));  // 1 or 2
    const VU res_d = hn::Or(sign, hn::Shr(sig, shift));
    VU res = hn::IfThenElse(hn::RebindMask(du, normal), res_n, res_d);

    // +-0 and tiny denormals: +-inf; +-inf: +-0; NaN: quieted
    res = hn::IfThenElse(hn::Le(ax, hn::Set(du, 0x0004000000000000ULL)), hn::Or(sign, hn::Set(du, 0x7ff0000000000000ULL)),
                         res);
    const auto infnan = hn::Ge(ax, hn::Set(du, 0x7ff0000000000000ULL));
    const auto nan = hn::Gt(ax, hn::Set(du, 0x7ff0000000000000ULL));
    res = hn::IfThenElse(infnan, hn::IfThenElse(nan, hn::Or(xb, hn::Set(du, 0x0008000000000000ULL)), sign), res);
    return hn::BitCast(d, res);
}

// The same for inputs with biased exponent in [1, 0x7fc] (normal, with a normal reciprocal); other
// inputs give unspecified results. Used by the main path of tan, whose denominators always lie in this
// range (functions/tan/README.md).
template <class D>
HWY_INLINE hn::VFromD<D> Rcp14Normal(D d, hn::VFromD<D> x) {
    const hn::RebindToUnsigned<D> du;
    const hn::Repartition<uint32_t, decltype(du)> d32;
    using VU = hn::VFromD<decltype(du)>;
    const VU xb = hn::BitCast(du, x);
    const VU f = hn::And(xb, hn::Set(du, 0x000fffffffffffffULL));
    const VU se = hn::Xor(xb, f);  // sign and exponent bits
    const VU i = hn::ShiftRight<46>(f);
    const VU l = hn::And(hn::ShiftRight<36>(f), hn::Set(du, 1023));
    const auto w32 = hwymath_trig::HWY_NAMESPACE::Lookup<64>(d32, kCoeff.v, kCoeff.v, hn::BitCast(d32, i));
    const VU w = hn::And(hn::BitCast(du, w32), hn::Set(du, 0xffffffffULL));
    const VU a = hn::ShiftRight<10>(w), b = hn::And(w, hn::Set(du, 1023));
    const VU bl = hn::MulEven(hn::BitCast(d32, b), hn::BitCast(d32, l));
    const VU q = hn::ShiftRight<9>(hn::Sub(hn::Add(hn::ShiftLeft<8>(a), hn::ShiftLeft<9>(b)), bl));
    // |result| bits = ((0x7fd - e) << 52) + ((q << 36) - 2^52) = (0x7fc << 52) + (q << 36) - (e << 52),
    // and ((0x7fe - e) << 52) = (0x7fc << 52) + 2^53 - (e << 52) for powers of 2
    const auto pow2 = hn::Eq(f, hn::Zero(du));
    const VU sum = hn::Add(hn::IfThenElse(pow2, hn::Set(du, 0x0020000000000000ULL), hn::ShiftLeft<36>(q)),
                           hn::Set(du, 0x7fc0000000000000ULL));
    const VU sgn = hn::And(xb, hn::Set(du, 0x8000000000000000ULL));
    return hn::BitCast(d, hn::Or(sgn, hn::Sub(sum, hn::Xor(se, sgn))));
}

#endif  // HWY_NATIVE_FMA && HWY_HAVE_FLOAT64

}  // namespace HWY_NAMESPACE
}  // namespace hwymath_rcp14
HWY_AFTER_NAMESPACE();

#endif  // HWYMATH_RCP14_HWY_INL_H_
