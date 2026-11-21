// Portable (Google Highway) version of Intel SVML's AVX-512 high-accuracy cbrt (svml/svml_cbrt_d_ha.c),
// with run-time dispatch to the best target. Same algorithm, constants and operation order; the
// AVX-512-only instructions are replaced by exact equivalents, so that on every target with true FMA
// (HWY_NATIVE_FMA) the results are bit-identical to Intel's __svml_cbrt8_ha on Intel CPUs:
//   vgetexppd/vgetmantpd  exponent and mantissa bits; denormals scaled by 2^64 first
//   vrcp14pd + vrndscalepd(0x48)  rcp = k/16 with k = 16 - #{i : |m| >= T_i}: vrcp14pd's output depends
//                         on the top 16 mantissa bits only, and its rounding to 4 fractional bits switches at
//                         8 thresholds T_i (multiples of 2^-16, rcp14_model.h, checked exhaustively by
//                         tests/test_rcp14.c); exact rounding of 16/m would differ for ~1 in 13000 inputs
//   vrndscalepd(0x09)     Floor
//   vpermpd/vpermt2pd     TableLookupLanes/TwoTablesLookupLanes on 4- and 8-entry tables (k = 16 by a select)
//   vscalefpd             multiplication by 2^K built from exponent bits (the result is always normal)
//   0, inf, NaN           x + x, which is what the AVX-512 code returns for them (it needs no callout)
// On AVX-512 targets (HWYMATH_NATIVE_AVX512) none of this is emulated: CbrtNative is Intel's instruction
// sequence (vgetexppd, vgetmantpd, vrcp14pd, vrndscalepd, vpermpd, vpermt2pd, vscalefpd), with no
// special-value handling, exactly like the assembly.
// Targets without FMA (SSE2/SSSE3/SSE4, EMU128, SCALAR): std::cbrt.
// Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.

#include <cmath>
#include <cstddef>
#include <cstdint>

#include "functions/cbrt/hwy/hwy_cbrt_constants.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "functions/cbrt/hwy/hwy_cbrt.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"
#include "hwymath/hwy/avx512_native-inl.h"

HWY_BEFORE_NAMESPACE();
namespace hwycbrt {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

#if HWYMATH_NATIVE_AVX512

// vpermpd / vpermt2pd with the raw index vector: the instructions use only its low 3 / 4 bits, so no index
// masking is needed (Highway's TableLookupLanes requires in-range indices, which debug builds assert)
HWY_INLINE hn::Vec512<double> PermVar(hn::Vec512<double> table, hn::Vec512<uint64_t> idx) {
    return hn::Vec512<double>{_mm512_permutexvar_pd(idx.raw, table.raw)};
}
HWY_INLINE hn::Vec512<double> PermVar2(hn::Vec512<double> lo, hn::Vec512<uint64_t> idx, hn::Vec512<double> hi) {
    return hn::Vec512<double>{_mm512_permutex2var_pd(lo.raw, idx.raw, hi.raw)};
}

// Intel's __svml_cbrt8_ha instruction for instruction (operation order of svml/svml_cbrt_d_ha.c). Denormals
// are handled by vgetexppd/vgetmantpd, and 0, +-inf, NaN come out right by themselves: no branch, no select.
template <class D>
HWY_INLINE hn::VFromD<D> CbrtNative(D d, hn::VFromD<D> x) {
    namespace nat = hwymath::HWY_NAMESPACE;
    const hn::RebindToUnsigned<D> du;
    const auto E = nat::GetExp(x);                                             // floor(log2|x|)
    const auto m = nat::GetMant<0>(x);                                         // +-[1, 2), sign of x
    const auto rcp = hn::ApproximateReciprocal(m);                             // vrcp14pd
    const auto ES = hn::Add(E, hn::Set(d, kShifter));                          // 1.5*2^52 + E
    const auto sgn = hn::And(x, hn::BitCast(d, hn::Set(du, 0x8000000000000000ULL)));
    const auto rcpR = nat::RoundScale<0x48>(rcp);                              // to 4 fractional bits
    const auto t = hn::MulSub(hn::Set(d, kOneThird), ES, hn::Set(d, kTwo51));  // E/3 + 1/4 (+ tiny)
    const auto R = hn::MulSub(m, rcpR, hn::Set(d, kOne));                      // m*rcp - 1
    const auto K = nat::RoundScale<0x09>(t);                                   // floor(E/3)
    const auto idx = hn::ShiftRight<49>(hn::BitCast(du, rcpR));                // k - 8 in bits 0..3
    const auto jS = hn::NegMulAdd(hn::Set(d, kThree), K, ES);                  // 1.5*2^52 + (E mod 3)
    auto p = hn::MulAdd(hn::Set(d, kC9), R, hn::Set(d, kC8));
    const auto j = hn::BitCast(du, jS);
    const auto Sh = PermVar(hn::LoadU(d, vSh), j);
    const auto Th = PermVar2(hn::LoadU(d, vTh), idx, hn::LoadU(d, vTh + 8));
    auto tl_hi = hn::LoadU(d, vTl + 8);
    __asm__("" : "+v"(tl_hi.raw));  // all zero: GCC would rebuild it with a vxorpd in every iteration
    const auto Tl = PermVar2(hn::LoadU(d, vTl), idx, tl_hi);
    const auto Sl = PermVar(hn::LoadU(d, vSl), j);

    const auto ShTh = hn::Mul(Sh, Th);
    const auto R2 = hn::Mul(R, R);
    const auto q45 = hn::MulAdd(hn::Set(d, kC5), R, hn::Set(d, kC4));
    const auto q67 = hn::MulAdd(hn::Set(d, kC7), R, hn::Set(d, kC6));
    p = hn::MulAdd(R2, p, q67);
    auto lo = hn::MulSub(Th, Sh, ShTh);  // (Sh*Th)_low
    p = hn::MulAdd(R2, p, q45);
    lo = hn::MulAdd(Tl, Sh, lo);
    lo = hn::MulAdd(Th, Sl, lo);
    const auto q23 = hn::MulAdd(hn::Set(d, kC3), R, hn::Set(d, kC2));
    p = hn::MulAdd(R2, p, q23);
    const auto q01 = hn::MulAdd(hn::Set(d, kC1), R, hn::Set(d, kC0));
    const auto ShThR = hn::Mul(ShTh, R);
    p = hn::MulAdd(R2, p, q01);
    lo = hn::MulAdd(p, ShThR, lo);
    const auto y = hn::Add(lo, ShTh);
    return hn::Or(nat::Scalef(y, K), sgn);
}

void CbrtArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    size_t i = 0;
    for (; i + N <= len; i += N) hn::StoreU(CbrtNative(d, hn::LoadU(d, x + i)), d, y + i);
    if (i < len) hn::StoreN(CbrtNative(d, hn::LoadN(d, x + i, len - i)), d, y + i, len - i);
}

#elif HWY_NATIVE_FMA && HWY_HAVE_FLOAT64

// lookup of a 4-entry table (padded copy in pad), idx < 4
template <class D, class VI>
HWY_INLINE hn::VFromD<D> Lookup4(D d, const double *pad, VI idx) {
    const size_t N = hn::Lanes(d);
    if (N >= 4 && N <= 64) return hn::TableLookupLanes(hn::LoadU(d, pad), hn::IndicesFromVec(d, idx));
    if (N == 2) return hn::TwoTablesLookupLanes(d, hn::LoadU(d, pad), hn::LoadU(d, pad + 2), hn::IndicesFromVec(d, idx));
    return hn::GatherIndex(d, pad, idx);  // single-lane vectors (not on x86 FMA targets)
}

// lookup of an 8-entry table: pad = the table padded, lo/hi = its halves padded; idx < 8
template <class D, class VI>
HWY_INLINE hn::VFromD<D> Lookup8(D d, const double *pad, const double *lo, const double *hi, VI idx) {
    const size_t N = hn::Lanes(d);
    if (N >= 8 && N <= 64) return hn::TableLookupLanes(hn::LoadU(d, pad), hn::IndicesFromVec(d, idx));
    if (N == 4) return hn::TwoTablesLookupLanes(d, hn::LoadU(d, pad), hn::LoadU(d, pad + 4), hn::IndicesFromVec(d, idx));
    if (N == 1) return hn::GatherIndex(d, pad, idx);
    // N = 2, or 6 (SVE-384): two 4-entry lookups
    const hn::RebindToSigned<D> di;
    const auto i3 = hn::And(idx, hn::Set(di, 3));
    const auto upper = hn::RebindMask(d, hn::Ne(hn::And(idx, hn::Set(di, 4)), hn::Zero(di)));
    return hn::IfThenElse(upper, Lookup4(d, hi, i3), Lookup4(d, lo, i3));
}

// kAllNormal: all lanes are normal numbers (the common case: no denormal/0/inf/NaN handling needed)
template <bool kAllNormal, class D>
HWY_INLINE hn::VFromD<D> CbrtCore(D d, hn::VFromD<D> x, hn::VFromD<D> ax) {
    const hn::RebindToUnsigned<D> du;
    const hn::RebindToSigned<D> di;
    const auto sign = hn::And(hn::BitCast(du, x), hn::Set(du, 0x8000000000000000ULL));

    // getexp/getmant: denormals are scaled by 2^64 first (exact) and the exponent corrected
    hn::VFromD<decltype(du)> xb, bias;
    if constexpr (kAllNormal) {
        xb = hn::BitCast(du, ax);
        bias = hn::Set(du, kShifterBits - 1023);
    } else {
        const auto denorm = hn::Lt(ax, hn::Set(d, 0x1p-1022));
        xb = hn::BitCast(du, hn::IfThenElse(denorm, hn::Mul(ax, hn::Set(d, 0x1p64)), ax));
        bias = hn::IfThenElse(hn::RebindMask(du, denorm), hn::Set(du, kShifterBits - 1023 - 64),
                              hn::Set(du, kShifterBits - 1023));
    }
    const auto ES = hn::BitCast(d, hn::Add(hn::ShiftRight<52>(xb), bias));  // 1.5*2^52 + E (vaddpd, exact)
    const auto m = hn::BitCast(d, hn::Or(hn::And(xb, hn::Set(du, 0x000fffffffffffffULL)),
                                         hn::Set(du, 0x3ff0000000000000ULL)));  // |getmant(x)| in [1, 2)

    // vrndscalepd(vrcp14pd(m), 0x48) = k/16, k = 8..16; kb = bits of 2^52 + k
    // count the thresholds below m (a true mask is -1)
    auto kb = hn::Set(di, int64_t{0x4330000000000010});  // bits of 2^52 + 16
#define CBRT_STEP(i) kb = hn::Add(kb, hn::BitCast(di, hn::VecFromMask(d, hn::Ge(m, hn::Set(d, cbrt_rcp14_threshold[i])))))
    CBRT_STEP(0); CBRT_STEP(1); CBRT_STEP(2); CBRT_STEP(3); CBRT_STEP(4); CBRT_STEP(5); CBRT_STEP(6); CBRT_STEP(7);
#undef CBRT_STEP
    const auto rcp = hn::Mul(hn::Sub(hn::BitCast(d, kb), hn::Set(d, 0x1p52)), hn::Set(d, 0.0625));
    const auto idx = hn::And(kb, hn::Set(di, 7));  // k - 8 for k = 8..15 (k = 16: selected below)
    const auto k16 = hn::Lt(m, hn::Set(d, cbrt_rcp14_threshold[0]));

    const auto t = hn::MulSub(hn::Set(d, kOneThird), ES, hn::Set(d, kTwo51));  // E/3 + 1/4 (+ tiny)
    const auto R = hn::MulSub(m, rcp, hn::Set(d, kOne));                       // m*rcp - 1
    const auto K = hn::Floor(t);                                                 // floor(E/3)
    const auto jS = hn::NegMulAdd(hn::Set(d, kThree), K, ES);                    // 1.5*2^52 + (E mod 3)
    auto p = hn::MulAdd(hn::Set(d, kC9), R, hn::Set(d, kC8));
    const auto j = hn::And(hn::BitCast(di, jS), hn::Set(di, 3));
    const auto Sh = Lookup4(d, kShPad.v, j);
    const auto Th = hn::IfThenElse(k16, hn::Set(d, 1.0), Lookup8(d, kThPad.v, kThLo.v, kThHi.v, idx));
    const auto Tl = hn::IfThenZeroElse(k16, Lookup8(d, kTlPad.v, kTlLo.v, kTlHi.v, idx));
    const auto Sl = Lookup4(d, kSlPad.v, j);

    const auto ShTh = hn::Mul(Sh, Th);
    const auto R2 = hn::Mul(R, R);
    const auto q45 = hn::MulAdd(hn::Set(d, kC5), R, hn::Set(d, kC4));
    const auto q67 = hn::MulAdd(hn::Set(d, kC7), R, hn::Set(d, kC6));
    p = hn::MulAdd(R2, p, q67);
    auto lo = hn::MulSub(Th, Sh, ShTh);  // (Sh*Th)_low
    p = hn::MulAdd(R2, p, q45);
    lo = hn::MulAdd(Tl, Sh, lo);
    lo = hn::MulAdd(Th, Sl, lo);
    const auto q23 = hn::MulAdd(hn::Set(d, kC3), R, hn::Set(d, kC2));
    p = hn::MulAdd(R2, p, q23);
    const auto q01 = hn::MulAdd(hn::Set(d, kC1), R, hn::Set(d, kC0));
    const auto ShThR = hn::Mul(ShTh, R);
    p = hn::MulAdd(R2, p, q01);
    lo = hn::MulAdd(p, ShThR, lo);
    const auto y = hn::Add(lo, ShTh);

    // scalef(y, K) = y * 2^K, exact: K in [-358, 341] and y in [1, 2), so the result is normal
    const auto scale = hn::ShiftLeft<52>(hn::BitCast(du, hn::Add(K, hn::Set(d, kShifter + 1023))));
    const auto res = hn::BitCast(d, hn::Or(hn::BitCast(du, hn::Mul(y, hn::BitCast(d, scale))), sign));
    if constexpr (kAllNormal) return res;
    const auto finite_nonzero = hn::And(hn::Gt(ax, hn::Zero(d)), hn::Lt(ax, hn::Inf(d)));
    return hn::IfThenElse(finite_nonzero, res, hn::Add(x, x));
}

template <class D>
HWY_INLINE hn::VFromD<D> Cbrt(D d, hn::VFromD<D> x) {
    const auto ax = hn::Abs(x);
    const auto normal = hn::And(hn::Ge(ax, hn::Set(d, 0x1p-1022)), hn::Le(ax, hn::Set(d, 0x1.fffffffffffffp+1023)));
    if (HWY_LIKELY(hn::AllTrue(d, normal))) return CbrtCore<true>(d, x, ax);
    return CbrtCore<false>(d, x, ax);
}

void CbrtArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    size_t i = 0;
    for (; i + N <= len; i += N) hn::StoreU(Cbrt(d, hn::LoadU(d, x + i)), d, y + i);
    if (i < len) hn::StoreN(Cbrt(d, hn::LoadN(d, x + i, len - i)), d, y + i, len - i);
}

#else  // no fused multiply-add: the algorithm's FMA steps would not be reproduced

void CbrtArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    for (size_t i = 0; i < len; i++) y[i] = std::cbrt(x[i]);
}

#endif

size_t LanesF64() { return hn::Lanes(hn::ScalableTag<double>()); }

int Variant() {  // 1: SVML algorithm, 0: std::cbrt fallback
#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64
    return 1;
#else
    return 0;
#endif
}

}  // namespace HWY_NAMESPACE
}  // namespace hwycbrt
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace hwycbrt {
HWY_EXPORT(CbrtArray);
HWY_EXPORT(Variant);
HWY_EXPORT(LanesF64);
// dispatch inside the namespace: with a single (static) target, HWY_DYNAMIC_DISPATCH(F) becomes N_<TARGET>::F
static void DispatchCbrt(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(CbrtArray)(x, y, len); }
static int DispatchVariant() { return HWY_DYNAMIC_DISPATCH(Variant)(); }
static size_t DispatchLanes() { return HWY_DYNAMIC_DISPATCH(LanesF64)(); }
}  // namespace hwycbrt

extern "C" {
void hwy_cbrt_array(const double *x, double *y, size_t len) { hwycbrt::DispatchCbrt(x, y, len); }
int hwy_cbrt_variant(void) { return hwycbrt::DispatchVariant(); }
size_t hwy_cbrt_lanes(void) { return hwycbrt::DispatchLanes(); }
}
#endif  // HWY_ONCE
