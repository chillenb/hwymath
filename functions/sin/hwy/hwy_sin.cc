// Portable (Google Highway) port of SVML's high-accuracy double sin (svml/svml_sin_d_ha.c, the
// translation of svml_z0_sin_d_ha.s), with run-time dispatch to the best target. Per target:
//   - targets with true FMA (HWY_NATIVE_FMA): the SVML algorithm, operation for operation, so that the
//     results are bit-identical to Intel's AVX-512 code. The AVX-512 features it uses are emulated
//     exactly: vpermt2pd (16-entry tables) by TableLookupLanes/TwoTablesLookupLanes (trig_hwy-inl.h),
//     the {rn-sae} operations by ordinary round-to-nearest operations, masked operations by selects,
//     vpmullq and the gathers of the large-argument path by MulEven and scalar loads;
//   - targets without FMA (SSE2/SSSE3/SSE4, EMU128): std::sin.
// Vectors with a lane |x| > 2^24 (or inf/NaN) are recomputed with the large-argument path (out of line);
// lanes with x = +-inf or NaN then use the scalar SVML callout (svml_dsin_ha_cout_rare).
// Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.

#include <cmath>
#include <cstddef>
#include <cstdint>

extern "C" int svml_dsin_ha_cout_rare(const double *px, double *pr);

#include "functions/sin/hwy/trig_tables.h"
#include "functions/sin/svml/svml_sin_d_ha_constants.inc"

#ifndef HWYMATH_HWY_SIN_TABLES
#define HWYMATH_HWY_SIN_TABLES
namespace hwysin {
constexpr hwymath_trig::Padded64<double, 16> kShPad(sin_vSh);
constexpr hwymath_trig::Padded64<double, 16> kChPad(sin_vCh);
constexpr hwymath_trig::Padded64<uint64_t, 16> kClSlPad(sin_vClSl);
}  // namespace hwysin
#endif

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "functions/sin/hwy/hwy_sin.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"
#include "functions/sin/hwy/trig_hwy-inl.h"

HWY_BEFORE_NAMESPACE();
namespace hwysin {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64
using hwymath_trig::HWY_NAMESPACE::Lookup;

// sin(x); kLarge: also the large-argument path for the lanes with |x| > 2^24 (and inf/NaN)
template <bool kLarge, class D>
HWY_INLINE hn::VFromD<D> SinCore(D d, hn::VFromD<D> x) {
    namespace ht = hwymath_trig::HWY_NAMESPACE;
    const hn::RebindToUnsigned<D> du;
    const auto y = hn::MulAdd(x, hn::Set(d, sin_vInvPi), hn::Set(d, sin_vRShift));
    const auto yi = hn::BitCast(du, y);
    const auto mzero = hn::Eq(hn::BitCast(du, x), hn::Set(du, sin_vMinusZero));
    const auto dn = hn::Sub(y, hn::Set(d, sin_vRShift));
    const auto ch = Lookup<16>(d, sin_vCh, kChPad.v, yi);
    const auto sh = Lookup<16>(d, sin_vSh, kShPad.v, yi);
    const auto cl = hn::BitCast(d, Lookup<16>(du, sin_vClSl, kClSlPad.v, yi));
    const auto ys = hn::IfThenElse(mzero, hn::Xor(yi, hn::Set(du, sin_vZeroSignBit)), yi);  // sin(-0) = -0
    const auto r1 = hn::NegMulAdd(dn, hn::Set(d, sin_vPI1), x);
    const auto n4 = hn::ShiftRight<4>(ys);

    const auto r2 = hn::NegMulAdd(dn, hn::Set(d, sin_vPI2), r1);
    const auto r = hn::NegMulAdd(dn, hn::Set(d, sin_vPI3), r2);
    auto e2 = hn::Sub(r1, r2);
    const auto rr = hn::Mul(r, r);
    e2 = hn::NegMulAdd(dn, hn::Set(d, sin_vPI2), e2);
    const auto e3 = hn::MulAdd(hn::Set(d, sin_vPI3), dn, hn::Sub(r, r2));
    const auto rlo = hn::Sub(e2, e3);
    auto cp = hn::MulAdd(hn::Set(d, sin_vC4), rr, hn::Set(d, sin_vC3));
    auto sp = hn::MulAdd(hn::Set(d, sin_vS4), rr, hn::Set(d, sin_vS3));
    cp = hn::MulAdd(rr, cp, hn::Set(d, sin_vC2));
    const auto sl = hn::BitCast(d, hn::ShiftLeft<32>(hn::BitCast(du, cl)));
    cp = hn::MulAdd(rr, cp, hn::Set(d, sin_vC1));
    sp = hn::MulAdd(rr, sp, hn::Set(d, sin_vS2));
    const auto hi = hn::MulAdd(r, ch, sh);  // Sh + R*Ch
    sp = hn::MulAdd(rr, sp, hn::Set(d, sin_vS1));
    const auto rch_l = hn::MulSub(ch, r, hn::Sub(hi, sh));  // error of R*Ch
    const auto csr = hn::NegMulAdd(r, sh, ch);
    const auto rch = hn::Mul(ch, r);
    const auto shc = hn::Mul(sh, cp);
    const auto lo1 = hn::MulAdd(r, cl, rch_l);
    const auto lo2 = hn::MulAdd(rlo, csr, sl);
    const auto sign = hn::BitCast(d, hn::ShiftLeft<63>(n4));
    auto q = hn::MulAdd(rch, sp, shc);
    q = hn::MulAdd(rr, q, hn::Add(lo1, lo2));
    auto res = hn::Xor(hn::Add(hi, q), sign);

    if constexpr (kLarge) {
        const auto ax = hn::Abs(x);
        const auto large = hn::Not(hn::Le(ax, hn::Set(d, sin_vRangeVal)));
        hn::VFromD<D> phi, plo, z, zlo;
        hn::VFromD<decltype(du)> n;
        ht::TrigReduceLarge(d, ax, phi, plo, n);
        ht::TrigReduced(d, ax, phi, plo, z, zlo);
        const auto v = ht::SinCosLargeTail(d, z, zlo, n, sin_hT, sin_hS1, sin_hS2, sin_hC1, sin_hC2, sin_hC3);
        res = hn::IfThenElse(large, hn::Xor(v, hn::Xor(ax, x)), res);
    }
    return res;
}

// rare: a vector with |x| > 2^24, inf or NaN: recompute with the large path, then the scalar callout
static HWY_NOINLINE __attribute__((cold)) void SinSpecial(const double *x, double *y, size_t n) {
    const hn::ScalableTag<double> d;
    hn::StoreN(SinCore<true>(d, hn::LoadN(d, x, n)), d, y, n);
    for (size_t k = 0; k < n; k++) {
        if (!std::isfinite(x[k])) svml_dsin_ha_cout_rare(&x[k], &y[k]);
    }
}

void SinArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    const auto range = hn::Set(d, sin_vRangeVal);
    size_t i = 0;
    for (; i + N <= len; i += N) {
        const auto v = hn::LoadU(d, x + i);
        hn::StoreU(SinCore<false>(d, v), d, y + i);
        if (HWY_UNLIKELY(!hn::AllTrue(d, hn::Le(hn::Abs(v), range)))) SinSpecial(x + i, y + i, N);
    }
    if (i < len) {
        const size_t rem = len - i;
        const auto v = hn::LoadN(d, x + i, rem);
        hn::StoreN(SinCore<false>(d, v), d, y + i, rem);
        if (!hn::AllTrue(d, hn::Le(hn::Abs(v), range))) SinSpecial(x + i, y + i, rem);
    }
}

int Variant() { return 1; }  // the SVML algorithm

#else  // no fused multiply-add: the algorithm cannot be reproduced

void SinArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    for (size_t i = 0; i < len; i++) y[i] = std::sin(x[i]);
}

int Variant() { return 0; }  // std::sin

#endif

size_t LanesF64() {
#if HWY_HAVE_FLOAT64
    return hn::Lanes(hn::ScalableTag<double>());
#else
    return 1;
#endif
}

}  // namespace HWY_NAMESPACE
}  // namespace hwysin
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace hwysin {
HWY_EXPORT(SinArray);
HWY_EXPORT(Variant);
HWY_EXPORT(LanesF64);
// dispatch inside the namespace: with a single (static) target, HWY_DYNAMIC_DISPATCH(F) becomes N_<TARGET>::F
static void DispatchSin(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(SinArray)(x, y, len); }
static int DispatchVariant() { return HWY_DYNAMIC_DISPATCH(Variant)(); }
static size_t DispatchLanes() { return HWY_DYNAMIC_DISPATCH(LanesF64)(); }
}  // namespace hwysin

extern "C" {
void hwy_sin_array(const double *x, double *y, size_t len) { hwysin::DispatchSin(x, y, len); }
int hwy_sin_variant(void) { return hwysin::DispatchVariant(); }
size_t hwy_sin_lanes(void) { return hwysin::DispatchLanes(); }
}
#endif  // HWY_ONCE
