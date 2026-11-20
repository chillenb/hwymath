// Portable (Google Highway) port of SVML's high-accuracy double cos (svml/svml_cos_d_ha.c, the
// translation of svml_z0_cos_d_ha.s; same structure as functions/sin/hwy/hwy_sin.cc), with run-time
// dispatch to the best target. Per target:
//   - targets with true FMA (HWY_NATIVE_FMA): the SVML algorithm, operation for operation, so that the
//     results are bit-identical to Intel's AVX-512 code. The AVX-512 features it uses are emulated
//     exactly: vpermt2pd (16-entry tables) by TableLookupLanes/TwoTablesLookupLanes (trig_hwy-inl.h),
//     the {rn-sae} operations by ordinary round-to-nearest operations, masked operations by selects,
//     vpmullq and the gathers of the large-argument path by MulEven and scalar loads;
//   - targets without FMA (SSE2/SSSE3/SSE4, EMU128): std::cos.
// Vectors with a lane |x| > 2^24 (or inf/NaN) are recomputed with the large-argument path (out of line);
// lanes with x = +-inf or NaN then use the scalar SVML callout (svml_dcos_ha_cout_rare).
// Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.

#include <cmath>
#include <cstddef>
#include <cstdint>

extern "C" int svml_dcos_ha_cout_rare(const double *px, double *pr);

#include "functions/sin/hwy/trig_tables.h"
#include "functions/cos/svml/svml_cos_d_ha_constants.inc"

#ifndef HWYMATH_HWY_COS_TABLES
#define HWYMATH_HWY_COS_TABLES
namespace hwycos {
constexpr hwymath_trig::Padded64<double, 16> kShPad(cos_vSh);
constexpr hwymath_trig::Padded64<double, 16> kChPad(cos_vCh);
constexpr hwymath_trig::Padded64<uint64_t, 16> kClSlPad(cos_vClSl);
}  // namespace hwycos
#endif

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "functions/cos/hwy/hwy_cos.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"
#include "functions/sin/hwy/trig_hwy-inl.h"

HWY_BEFORE_NAMESPACE();
namespace hwycos {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64
using hwymath_trig::HWY_NAMESPACE::Lookup;

// cos(x); kLarge: also the large-argument path for the lanes with |x| > 2^24 (and inf/NaN)
template <bool kLarge, class D>
HWY_INLINE hn::VFromD<D> CosCore(D d, hn::VFromD<D> x) {
    namespace ht = hwymath_trig::HWY_NAMESPACE;
    const hn::RebindToUnsigned<D> du;
    const auto y = hn::MulAdd(x, hn::Set(d, cos_vInvPi), hn::Set(d, cos_vRShift));
    const auto yi = hn::BitCast(du, y);
    const auto dn = hn::Sub(y, hn::Set(d, cos_vRShift));
    const auto sh = Lookup<16>(d, cos_vSh, kShPad.v, yi);
    const auto ch = Lookup<16>(d, cos_vCh, kChPad.v, yi);
    const auto cl = hn::BitCast(d, Lookup<16>(du, cos_vClSl, kClSlPad.v, yi));
    const auto r1 = hn::NegMulAdd(dn, hn::Set(d, cos_vPI1), x);
    const auto n4 = hn::ShiftRight<4>(yi);

    const auto r2 = hn::NegMulAdd(dn, hn::Set(d, cos_vPI2), r1);
    const auto r = hn::NegMulAdd(dn, hn::Set(d, cos_vPI3), r2);
    auto e2 = hn::Sub(r1, r2);
    const auto rr = hn::Mul(r, r);
    e2 = hn::NegMulAdd(dn, hn::Set(d, cos_vPI2), e2);
    const auto e3 = hn::MulAdd(hn::Set(d, cos_vPI3), dn, hn::Sub(r, r2));
    auto cp = hn::MulAdd(hn::Set(d, cos_vC4), rr, hn::Set(d, cos_vC3));
    const auto rlo = hn::Sub(e2, e3);
    auto sp = hn::MulAdd(hn::Set(d, cos_vS4), rr, hn::Set(d, cos_vS3));
    cp = hn::MulAdd(rr, cp, hn::Set(d, cos_vC2));
    const auto sl = hn::BitCast(d, hn::ShiftLeft<32>(hn::BitCast(du, cl)));
    cp = hn::MulAdd(rr, cp, hn::Set(d, cos_vC1));
    sp = hn::MulAdd(rr, sp, hn::Set(d, cos_vS2));
    const auto hi = hn::NegMulAdd(r, sh, ch);  // Ch - R*Sh
    sp = hn::MulAdd(rr, sp, hn::Set(d, cos_vS1));
    const auto rsh_l = hn::MulSub(r, sh, hn::Sub(ch, hi));  // error of R*Sh
    const auto sc = hn::MulAdd(r, ch, sh);
    const auto rsh = hn::Mul(sh, r);
    const auto chc = hn::Mul(ch, cp);
    const auto lo1 = hn::MulAdd(r, sl, rsh_l);
    const auto lo2 = hn::NegMulAdd(rlo, sc, cl);
    const auto sign = hn::BitCast(d, hn::ShiftLeft<63>(n4));
    auto q = hn::NegMulAdd(rsh, sp, chc);
    q = hn::MulAdd(rr, q, hn::Sub(lo2, lo1));
    auto res = hn::Xor(hn::Add(hi, q), sign);

    if constexpr (kLarge) {
        const auto large = hn::Not(hn::Le(hn::Abs(x), hn::Set(d, cos_vRangeVal)));
        hn::VFromD<D> phi, plo, z, zlo;
        hn::VFromD<decltype(du)> n;
        ht::TrigReduceLarge(d, x, phi, plo, n);
        ht::TrigReduced(d, x, phi, plo, z, zlo);
        const auto v = ht::SinCosLargeTail(d, z, zlo, n, cos_hT, cos_hS1, cos_hS2, cos_hC1, cos_hC2, cos_hC3);
        res = hn::IfThenElse(large, v, res);
    }
    return res;
}

// rare: a vector with |x| > 2^24, inf or NaN: recompute with the large path, then the scalar callout
static HWY_NOINLINE __attribute__((cold)) void CosSpecial(const double *x, double *y, size_t n) {
    const hn::ScalableTag<double> d;
    hn::StoreN(CosCore<true>(d, hn::LoadN(d, x, n)), d, y, n);
    for (size_t k = 0; k < n; k++) {
        if (!std::isfinite(x[k])) svml_dcos_ha_cout_rare(&x[k], &y[k]);
    }
}

void CosArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    const auto range = hn::Set(d, cos_vRangeVal);
    size_t i = 0;
    for (; i + N <= len; i += N) {
        const auto v = hn::LoadU(d, x + i);
        hn::StoreU(CosCore<false>(d, v), d, y + i);
        if (HWY_UNLIKELY(!hn::AllTrue(d, hn::Le(hn::Abs(v), range)))) CosSpecial(x + i, y + i, N);
    }
    if (i < len) {
        const size_t rem = len - i;
        const auto v = hn::LoadN(d, x + i, rem);
        hn::StoreN(CosCore<false>(d, v), d, y + i, rem);
        if (!hn::AllTrue(d, hn::Le(hn::Abs(v), range))) CosSpecial(x + i, y + i, rem);
    }
}

int Variant() { return 1; }  // the SVML algorithm

#else  // no fused multiply-add: the algorithm cannot be reproduced

void CosArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    for (size_t i = 0; i < len; i++) y[i] = std::cos(x[i]);
}

int Variant() { return 0; }  // std::cos

#endif

size_t LanesF64() {
#if HWY_HAVE_FLOAT64
    return hn::Lanes(hn::ScalableTag<double>());
#else
    return 1;
#endif
}

}  // namespace HWY_NAMESPACE
}  // namespace hwycos
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace hwycos {
HWY_EXPORT(CosArray);
HWY_EXPORT(Variant);
HWY_EXPORT(LanesF64);
// dispatch inside the namespace: with a cosgle (static) target, HWY_DYNAMIC_DISPATCH(F) becomes N_<TARGET>::F
static void DispatchCos(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(CosArray)(x, y, len); }
static int DispatchVariant() { return HWY_DYNAMIC_DISPATCH(Variant)(); }
static size_t DispatchLanes() { return HWY_DYNAMIC_DISPATCH(LanesF64)(); }
}  // namespace hwycos

extern "C" {
void hwy_cos_array(const double *x, double *y, size_t len) { hwycos::DispatchCos(x, y, len); }
int hwy_cos_variant(void) { return hwycos::DispatchVariant(); }
size_t hwy_cos_lanes(void) { return hwycos::DispatchLanes(); }
}
#endif  // HWY_ONCE
