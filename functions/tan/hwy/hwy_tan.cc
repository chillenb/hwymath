// Portable (Google Highway) port of SVML's high-accuracy double tan (svml/svml_tan_d_ha.c, the
// translation of svml_z0_tan_d_ha.s), with run-time dispatch to the best target. Per target:
//   - targets with true FMA (HWY_NATIVE_FMA): the SVML algorithm, operation for operation, so that the
//     results are bit-identical to Intel's AVX-512 code. The AVX-512 features it uses are emulated
//     exactly: vpermt2pd (16-entry tables) by TableLookupLanes/TwoTablesLookupLanes, vrcp14pd by
//     Intel's published reference algorithm (rcp14_hwy-inl.h), the {rn-sae} operations by ordinary
//     round-to-nearest operations, masked operations by selects, vpmullq and the gathers of the
//     large-argument path by MulEven and scalar loads (functions/sin/hwy/trig_hwy-inl.h);
//   - targets without FMA (SSE2/SSSE3/SSE4, EMU128): std::tan.
// Vectors with a lane |x| > 65536.0035 (or inf/NaN) are recomputed with the large-argument path (out of
// line); lanes with x = +-inf or NaN then use the scalar SVML callout (svml_dtan_ha_cout_rare).
// Also exports the vrcp14pd emulations for testing (hwy_tan_rcp14_array, hwy_tan_rcp14_normal_array).
// Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.

#include <cmath>
#include <cstddef>
#include <cstdint>

extern "C" int svml_dtan_ha_cout_rare(const double *px, double *pr);

#include "functions/sin/hwy/trig_tables.h"
#include "functions/tan/svml/svml_tan_d_ha_constants.inc"

#ifndef HWYMATH_HWY_TAN_TABLES
#define HWYMATH_HWY_TAN_TABLES
namespace hwytan {
constexpr hwymath_trig::Padded64<double, 16> kThPad(tan_vTh);
constexpr hwymath_trig::Padded64<double, 16> kTlPad(tan_vTl);
}  // namespace hwytan
#endif

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "functions/tan/hwy/hwy_tan.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"
#include "functions/sin/hwy/trig_hwy-inl.h"
#include "functions/tan/hwy/rcp14_hwy-inl.h"

HWY_BEFORE_NAMESPACE();
namespace hwytan {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64
using hwymath_rcp14::HWY_NAMESPACE::Rcp14;
using hwymath_rcp14::HWY_NAMESPACE::Rcp14Normal;
using hwymath_trig::HWY_NAMESPACE::Lookup;

// tan(x); kLarge: also the large-argument path for the lanes with |x| > 65536.0035 (and inf/NaN)
template <bool kLarge, class D>
HWY_INLINE hn::VFromD<D> TanCore(D d, hn::VFromD<D> x) {
    namespace ht = hwymath_trig::HWY_NAMESPACE;
    const hn::RebindToUnsigned<D> du;
    const hn::RebindToSigned<D> di;
    // reduction, main path
    const auto y = hn::MulAdd(x, hn::Set(d, tan_vInvPi), hn::Set(d, tan_vRShift));
    const auto dn = hn::Sub(y, hn::Set(d, tan_vRShift));
    const auto r1 = hn::NegMulAdd(dn, hn::Set(d, tan_vPI1), x);
    const auto r2 = hn::NegMulAdd(dn, hn::Set(d, tan_vPI2), r1);
    auto r = hn::NegMulAdd(dn, hn::Set(d, tan_vPI3), r2);
    const auto e2 = hn::NegMulAdd(dn, hn::Set(d, tan_vPI2), hn::Sub(r1, r2));
    const auto e3 = hn::MulAdd(hn::Set(d, tan_vPI3), dn, hn::Sub(r, r2));
    auto rlo = hn::Sub(e2, e3);
    auto idx = hn::BitCast(du, y);

    if constexpr (kLarge) {
        const auto large = hn::Not(hn::Le(hn::Abs(x), hn::Set(d, tan_vRangeVal)));
        hn::VFromD<D> phi, plo, z, zlo;
        hn::VFromD<decltype(du)> m;
        ht::TrigReduceLarge(d, x, phi, plo, m);
        ht::TrigReduced(d, x, phi, plo, z, zlo);
        // M = 16*k + i, k = (M + 8 - sign(z)) >> 4
        const auto zsign = hn::ShiftRight<63>(hn::BitCast(du, z));
        const auto k = hn::ShiftRight<4>(hn::Sub(hn::Add(m, hn::Set(du, tan_lHalf)), zsign));
        const auto i = hn::ConvertTo(d, hn::Sub(hn::BitCast(di, m), hn::BitCast(di, hn::ShiftLeft<4>(k))));
        const auto t = hn::MulAdd(hn::Set(d, tan_lPi256Hi), i, z);
        const auto u = hn::NegMulAdd(i, hn::Set(d, tan_lPi256Hi), t);
        const auto w = hn::MulAdd(hn::Set(d, tan_lPi256Lo), i, hn::Sub(z, u));
        r = hn::IfThenElse(large, t, r);
        idx = hn::IfThenElse(hn::RebindMask(du, large), k, idx);
        rlo = hn::IfThenElse(large, hn::Add(w, zlo), rlo);
    }

    // evaluation
    const auto rr = hn::Mul(r, r);
    const auto tl = Lookup<16>(d, tan_vTl, kTlPad.v, idx);
    const auto th = Lookup<16>(d, tan_vTh, kThPad.v, idx);
    auto p = hn::MulAdd(hn::Set(d, tan_vP5), rr, hn::Set(d, tan_vP4));
    p = hn::MulAdd(rr, p, hn::Set(d, tan_vP3));
    p = hn::MulAdd(rr, p, hn::Set(d, tan_vP2));
    p = hn::MulAdd(rr, p, hn::Set(d, tan_vP1));
    const auto q = hn::NegMulSub(rr, hn::Mul(p, r), rlo);  // -(R^3*poly) - Rlo
    const auto pt = hn::Sub(r, q);
    const auto plo = hn::Sub(hn::Sub(r, pt), q);
    const auto one = hn::Set(d, tan_vOne);
    const auto num = hn::Add(pt, th);
    const auto nlo = hn::Add(hn::Add(hn::Sub(pt, hn::Sub(num, th)), tl), plo);
    const auto den = hn::NegMulAdd(pt, th, one);
    auto dlo = hn::MulAdd(pt, th, hn::Sub(den, one));
    dlo = hn::MulAdd(plo, th, dlo);
    dlo = hn::MulAdd(pt, tl, dlo);
    // the denominator of the main path is always in the range of Rcp14Normal; the rare path uses the
    // complete emulation anyway (and its inf/NaN lanes are replaced by the callout)
    auto rcp = kLarge ? Rcp14(d, den) : Rcp14Normal(d, den);
    const auto e = hn::MulAdd(dlo, rcp, hn::NegMulAdd(den, rcp, one));
    rcp = hn::MulAdd(e, rcp, rcp);
    const auto quo = hn::Mul(rcp, num);
    auto g = hn::MulSub(quo, den, num);
    g = hn::NegMulAdd(quo, dlo, g);
    g = hn::Sub(g, nlo);
    return hn::NegMulAdd(rcp, g, quo);
}

// rare: a vector with |x| > 65536.0035, inf or NaN: recompute with the large path, then the scalar callout
static HWY_NOINLINE __attribute__((cold)) void TanSpecial(const double *x, double *y, size_t n) {
    const hn::ScalableTag<double> d;
    hn::StoreN(TanCore<true>(d, hn::LoadN(d, x, n)), d, y, n);
    for (size_t k = 0; k < n; k++) {
        if (!std::isfinite(x[k])) svml_dtan_ha_cout_rare(&x[k], &y[k]);
    }
}

void TanArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    const auto range = hn::Set(d, tan_vRangeVal);
    size_t i = 0;
    for (; i + N <= len; i += N) {
        const auto v = hn::LoadU(d, x + i);
        hn::StoreU(TanCore<false>(d, v), d, y + i);
        if (HWY_UNLIKELY(!hn::AllTrue(d, hn::Le(hn::Abs(v), range)))) TanSpecial(x + i, y + i, N);
    }
    if (i < len) {
        const size_t rem = len - i;
        const auto v = hn::LoadN(d, x + i, rem);
        hn::StoreN(TanCore<false>(d, v), d, y + i, rem);
        if (!hn::AllTrue(d, hn::Le(hn::Abs(v), range))) TanSpecial(x + i, y + i, rem);
    }
}

// the vrcp14pd emulation alone, for testing
void Rcp14Array(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    size_t i = 0;
    for (; i + N <= len; i += N) hn::StoreU(Rcp14(d, hn::LoadU(d, x + i)), d, y + i);
    if (i < len) hn::StoreN(Rcp14(d, hn::LoadN(d, x + i, len - i)), d, y + i, len - i);
}

void Rcp14NormalArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    size_t i = 0;
    for (; i + N <= len; i += N) hn::StoreU(Rcp14Normal(d, hn::LoadU(d, x + i)), d, y + i);
    if (i < len) hn::StoreN(Rcp14Normal(d, hn::LoadN(d, x + i, len - i)), d, y + i, len - i);
}

int Variant() { return 1; }  // the SVML algorithm

#else  // no fused multiply-add: the algorithm cannot be reproduced

void TanArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    for (size_t i = 0; i < len; i++) y[i] = std::tan(x[i]);
}

void Rcp14Array(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    for (size_t i = 0; i < len; i++) y[i] = 1.0 / x[i];  // not an emulation (unused: Variant() == 0)
}

void Rcp14NormalArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) { Rcp14Array(x, y, len); }

int Variant() { return 0; }  // std::tan

#endif

size_t LanesF64() {
#if HWY_HAVE_FLOAT64
    return hn::Lanes(hn::ScalableTag<double>());
#else
    return 1;
#endif
}

}  // namespace HWY_NAMESPACE
}  // namespace hwytan
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace hwytan {
HWY_EXPORT(TanArray);
HWY_EXPORT(Rcp14Array);
HWY_EXPORT(Rcp14NormalArray);
HWY_EXPORT(Variant);
HWY_EXPORT(LanesF64);
// dispatch inside the namespace: with a single (static) target, HWY_DYNAMIC_DISPATCH(F) becomes N_<TARGET>::F
static void DispatchTan(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(TanArray)(x, y, len); }
static void DispatchRcp14(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(Rcp14Array)(x, y, len); }
static void DispatchRcp14Normal(const double *x, double *y, size_t len) {
    HWY_DYNAMIC_DISPATCH(Rcp14NormalArray)(x, y, len);
}
static int DispatchVariant() { return HWY_DYNAMIC_DISPATCH(Variant)(); }
static size_t DispatchLanes() { return HWY_DYNAMIC_DISPATCH(LanesF64)(); }
}  // namespace hwytan

extern "C" {
void hwy_tan_array(const double *x, double *y, size_t len) { hwytan::DispatchTan(x, y, len); }
void hwy_tan_rcp14_array(const double *x, double *y, size_t len) { hwytan::DispatchRcp14(x, y, len); }
void hwy_tan_rcp14_normal_array(const double *x, double *y, size_t len) { hwytan::DispatchRcp14Normal(x, y, len); }
int hwy_tan_variant(void) { return hwytan::DispatchVariant(); }
size_t hwy_tan_lanes(void) { return hwytan::DispatchLanes(); }
}
#endif  // HWY_ONCE
