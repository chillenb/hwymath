// Portable (Google Highway) version of Intel SVML's AVX-512 atanh (svml/svml_atanh_d_ha.c), with run-time
// dispatch to the best target. Per target:
//   - targets with true FMA (HWY_NATIVE_FMA): the SVML algorithm, operation for operation, with the
//     AVX-512-only instructions emulated exactly (common/include/hwymath/hwy/avx512_emul-inl.h: VRCP14PD,
//     VGETEXPPD, 16-entry VPERMT2PD/VPERMI2PD lookups), so results are bit-identical to Intel's
//     __svml_atanh8_ha on every such target;
//   - targets without FMA (SSE2/SSSE3/SSE4): std::atanh.
// Lanes with |x| >= 1 or NaN use the scalar SVML callout (svml_datanh_ha_cout_rare).
// Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.

#include <cmath>
#include <cstddef>
#include <cstdint>

extern "C" int svml_datanh_ha_cout_rare(const double *px, double *pr);

#ifndef HWYATANH_CONSTANTS
#define HWYATANH_CONSTANTS
namespace hwyatanh {
#include "functions/atanh/svml/svml_atanh_d_ha_constants.inc"
}  // namespace hwyatanh
#endif

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "functions/atanh/hwy/hwy_atanh.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"
#include "hwymath/hwy/avx512_emul-inl.h"

HWY_BEFORE_NAMESPACE();
namespace hwyatanh {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64
namespace em = hwymath::HWY_NAMESPACE;

// scalar fix-up for |x| >= 1 or NaN (rare: kept out of line)
static HWY_NOINLINE __attribute__((cold)) void FixSpecial(const double *x, double *y, size_t n) {
    for (size_t k = 0; k < n; k++) {
        if (!(std::fabs(x[k]) < kOne)) svml_datanh_ha_cout_rare(&x[k], &y[k]);
    }
}

// The operations of svml_atanh8_ha() (same names, same order where it matters)
template <class D>
HWY_INLINE hn::VFromD<D> AtanhCore(D d, hn::VFromD<D> x) {
    const hn::RebindToUnsigned<D> du;
    const auto one = hn::Set(d, kOne);
    const auto y = hn::Abs(x);
    const auto Yp = hn::Add(y, one);
    const auto Ym = hn::Sub(one, y);
    const auto sign = hn::Xor(x, y);
    const auto Yph = hn::Sub(Yp, one);
    const auto nYmh = hn::Sub(Ym, one);
    const auto rcpP = em::Rcp14(d, Yp);
    const auto rcpM = em::Rcp14(d, Ym);
    const auto add = hn::Set(du, kRcpRoundAdd), mask = hn::Set(du, kRcpRoundMask);
    const auto Ypl = hn::Sub(y, Yph);
    const auto RcpP = hn::BitCast(d, hn::And(hn::Add(hn::BitCast(du, rcpP), add), mask));
    const auto RcpM = hn::BitCast(d, hn::And(hn::Add(hn::BitCast(du, rcpM), add), mask));
    const auto Yml = hn::Add(y, nYmh);
    auto Rp = hn::MulSub(RcpP, Yp, one);
    const auto Kp = em::GetExpPositiveNormal(d, RcpP);
    Rp = hn::MulAdd(Ypl, RcpP, Rp);
    const auto idxM = hn::ShiftRight<48>(hn::BitCast(du, RcpM));
    const auto idxP = hn::ShiftRight<48>(hn::BitCast(du, RcpP));
    auto Rm = hn::MulSub(Ym, RcpM, one);
    Rm = hn::NegMulAdd(Yml, RcpM, Rm);
    const auto Km = em::GetExpPositiveNormal(d, RcpM);
    const auto dK = hn::Sub(Km, Kp);
    const auto TlP = em::Lookup16(d, kLogTl, idxP);
    const auto TlM = em::Lookup16(d, kLogTl, idxM);
    const auto ThP = em::Lookup16(d, kLogTh, idxP);
    const auto ThM = em::Lookup16(d, kLogTh, idxM);
    const auto dTl = hn::Sub(TlM, TlP);
    auto pP = hn::MulAdd(hn::Set(d, kLogC9), Rp, hn::Set(d, kLogC8));
    auto pM = hn::MulAdd(hn::Set(d, kLogC9), Rm, hn::Set(d, kLogC8));
    const auto dTh = hn::Sub(ThM, ThP);
    pP = hn::MulAdd(Rp, pP, hn::Set(d, kLogC7));
    pM = hn::MulAdd(Rm, pM, hn::Set(d, kLogC7));
    const auto Hi0 = hn::MulAdd(hn::Set(d, kL2H), dK, dTh);
    pP = hn::MulAdd(Rp, pP, hn::Set(d, kLogC6));
    pM = hn::MulAdd(Rm, pM, hn::Set(d, kLogC6));
    const auto Lo0 = hn::MulAdd(hn::Set(d, kL2L), dK, dTl);
    const auto Rp2 = hn::Mul(Rp, Rp);
    pP = hn::MulAdd(Rp, pP, hn::Set(d, kLogC5));
    pM = hn::MulAdd(Rm, pM, hn::Set(d, kLogC5));
    const auto H1 = hn::Add(Rp, Hi0);
    pP = hn::MulAdd(Rp, pP, hn::Set(d, kLogC4));
    pM = hn::MulAdd(Rm, pM, hn::Set(d, kLogC4));
    const auto Rm2 = hn::Mul(Rm, Rm);
    pP = hn::MulAdd(Rp, pP, hn::Set(d, kLogC3));
    pM = hn::MulAdd(Rm, pM, hn::Set(d, kLogC3));
    const auto nRph = hn::Sub(Hi0, H1);
    const auto H2 = hn::Sub(H1, Rm);
    pP = hn::MulAdd(Rp, pP, hn::Set(d, kLogC2));
    pM = hn::MulAdd(Rm, pM, hn::Set(d, kLogC2));
    pP = hn::MulAdd(Rp, pP, hn::Set(d, kLogC1));
    pM = hn::MulAdd(Rm, pM, hn::Set(d, kLogC1));
    const auto Rpl = hn::Add(Rp, nRph);
    const auto nRmh = hn::Sub(H2, H1);
    pP = hn::MulAdd(Rp2, pP, Lo0);
    pM = hn::NegMulAdd(Rm2, pM, Rpl);
    const auto Rml = hn::Add(Rm, nRmh);
    const auto half = hn::Or(sign, hn::Set(d, kHalf));
    auto lo = hn::Add(pP, pM);
    lo = hn::Sub(lo, Rml);
    return hn::Mul(hn::Add(H2, lo), half);
}

void AtanhArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    const auto one = hn::Set(d, kOne);
    size_t i = 0;
    // two independent vectors per iteration: one evaluation is a long dependency chain (latency-bound)
    for (; i + 2 * N <= len; i += 2 * N) {
        const auto v0 = hn::LoadU(d, x + i);
        const auto v1 = hn::LoadU(d, x + i + N);
        hn::StoreU(AtanhCore(d, v0), d, y + i);
        hn::StoreU(AtanhCore(d, v1), d, y + i + N);
        if (HWY_UNLIKELY(!hn::AllTrue(d, hn::And(hn::Lt(hn::Abs(v0), one), hn::Lt(hn::Abs(v1), one)))))
            FixSpecial(x + i, y + i, 2 * N);
    }
    for (; i + N <= len; i += N) {
        const auto v = hn::LoadU(d, x + i);
        hn::StoreU(AtanhCore(d, v), d, y + i);
        if (HWY_UNLIKELY(!hn::AllTrue(d, hn::Lt(hn::Abs(v), one)))) FixSpecial(x + i, y + i, N);
    }
    if (i < len) {
        const size_t rem = len - i;
        const auto v = hn::LoadN(d, x + i, rem);
        hn::StoreN(AtanhCore(d, v), d, y + i, rem);
        FixSpecial(x + i, y + i, rem);
    }
}

#else  // no fused multiply-add: the algorithm's error-free steps would not be exact

void AtanhArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    for (size_t i = 0; i < len; i++) y[i] = std::atanh(x[i]);
}

#endif

size_t LanesF64() { return hn::Lanes(hn::ScalableTag<double>()); }

int Variant() {  // 1: SVML algorithm, 0: std::atanh fallback
    return HWY_NATIVE_FMA && HWY_HAVE_FLOAT64;
}

}  // namespace HWY_NAMESPACE
}  // namespace hwyatanh
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace hwyatanh {
HWY_EXPORT(AtanhArray);
HWY_EXPORT(Variant);
HWY_EXPORT(LanesF64);
// dispatch inside the namespace: with a single (static) target, HWY_DYNAMIC_DISPATCH(F) becomes N_<TARGET>::F
static void DispatchAtanh(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(AtanhArray)(x, y, len); }
static int DispatchVariant() { return HWY_DYNAMIC_DISPATCH(Variant)(); }
static size_t DispatchLanes() { return HWY_DYNAMIC_DISPATCH(LanesF64)(); }
}  // namespace hwyatanh

extern "C" {
void hwy_atanh_array(const double *x, double *y, size_t len) { hwyatanh::DispatchAtanh(x, y, len); }
int hwy_atanh_variant(void) { return hwyatanh::DispatchVariant(); }
size_t hwy_atanh_lanes(void) { return hwyatanh::DispatchLanes(); }
}
#endif  // HWY_ONCE
