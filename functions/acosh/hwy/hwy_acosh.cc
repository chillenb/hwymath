// Portable (Google Highway) version of Intel SVML's AVX-512 acosh (svml/svml_acosh_d_ha.c), with run-time
// dispatch to the best target. Per target:
//   - targets with true FMA (HWY_NATIVE_FMA): the SVML algorithm, operation for operation, with the
//     AVX-512-only instructions emulated exactly (common/include/hwymath/hwy/avx512_emul-inl.h: VRCP14PD,
//     VRSQRT14PD, VGETEXPPD, 16-entry VPERMT2PD lookups; masked operations as selects), so results are
//     bit-identical to Intel's __svml_acosh8_ha on every such target;
//   - targets without FMA (SSE2/SSSE3/SSE4): std::acosh.
// Lanes with x <= 1, x > DBL_MAX (inf) or NaN use the scalar SVML callout (svml_dacosh_ha_cout_rare).
// Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.

#include <cmath>
#include <cstddef>
#include <cstdint>

extern "C" int svml_dacosh_ha_cout_rare(const double *px, double *pr);

#ifndef HWYACOSH_CONSTANTS
#define HWYACOSH_CONSTANTS
namespace hwyacosh {
#include "functions/acosh/svml/svml_acosh_d_ha_constants.inc"
}  // namespace hwyacosh
#endif

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "functions/acosh/hwy/hwy_acosh.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"
#include "hwymath/hwy/avx512_emul-inl.h"

HWY_BEFORE_NAMESPACE();
namespace hwyacosh {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64
namespace em = hwymath::HWY_NAMESPACE;

// scalar fix-up for x <= 1, x > DBL_MAX or NaN (rare: kept out of line)
static HWY_NOINLINE __attribute__((cold)) void FixSpecial(const double *x, double *y, size_t n) {
    for (size_t k = 0; k < n; k++) {
        if (!(x[k] > kOne && x[k] <= kMaxX)) svml_dacosh_ha_cout_rare(&x[k], &y[k]);
    }
}

// The operations of svml_acosh8_ha() (same names, same order where it matters)
template <class D>
HWY_INLINE hn::VFromD<D> AcoshCore(D d, hn::VFromD<D> x) {
    const hn::RebindToUnsigned<D> du;
    const auto one = hn::Set(d, kOne);
    const auto x2 = hn::Mul(x, x);
    const auto Y0 = hn::Sub(x2, one);
    const auto t = hn::Sub(x2, Y0);
    const auto x2l = hn::MulSub(x, x, x2);
    const auto nOnel = hn::Sub(t, one);
    const auto Yl0 = hn::Add(x2l, nOnel);
    const auto Y = hn::Add(Y0, Yl0);
    const auto large = hn::Not(hn::Lt(x, hn::Set(d, kLargeX)));
    const auto R0 = em::Rsqrt14(d, Y);
    const auto notsmall = hn::Not(hn::Lt(Y, hn::Set(d, kSmallX)));  // true for NaN (x^2 overflow)
    auto ps = hn::MulAdd(hn::Set(d, kSmallC3), Y, hn::Set(d, kSmallC2));
    const auto Sh = hn::Mul(Y, R0);
    ps = hn::MulAdd(Y, ps, hn::Set(d, kSmallC1));
    const auto Ylh = hn::Sub(Y, Y0);
    const auto Shx = hn::Add(Sh, x);
    const auto pY = hn::Mul(ps, Y);
    const auto Yl = hn::Sub(Yl0, Ylh);
    const auto YR0l = hn::MulSub(R0, Y, Sh);
    auto Eh = hn::NegMulAdd(R0, Sh, one);
    const auto Sl = hn::MulAdd(R0, Yl, YR0l);
    Eh = hn::NegMulAdd(R0, Sl, Eh);
    const auto ShEh = hn::Mul(Sh, Eh);
    auto pS = hn::MulAdd(hn::Set(d, kSqrtC4), Eh, hn::Set(d, kSqrtC3));
    pS = hn::MulAdd(Eh, pS, hn::Set(d, kSqrtC2));
    const auto Shh = hn::Sub(Shx, x);
    pS = hn::MulAdd(Eh, pS, hn::Set(d, kSqrtC1));
    const auto Shl = hn::Sub(Sh, Shh);
    const auto Sc = hn::MulAdd(pS, ShEh, Sl);
    auto s = hn::Add(Sh, Sc);
    auto Xin = hn::Add(Shx, Sc);
    s = hn::MulAdd(pY, s, Sc);
    const auto Slh = hn::Sub(Xin, Shx);
    auto res = hn::Add(s, Sh);
    const auto Sll = hn::Sub(Sc, Slh);
    Xin = hn::IfThenElse(large, hn::Mul(x, hn::Set(d, kLargeScale)), Xin);
    auto Xinl = hn::Add(Sll, Shl);
    const auto rcp = em::Rcp14(d, Xin);
    Xinl = hn::IfThenZeroElse(large, Xinl);
    const auto Rcp = hn::BitCast(d, hn::And(hn::Add(hn::BitCast(du, rcp), hn::Set(du, kRcpRoundAdd)),
                                            hn::Set(du, kRcpRoundMask)));
    const auto Rh0 = hn::MulSub(Xin, Rcp, one);
    auto K = em::GetExpPositiveNormal(d, Rcp);
    const auto idx = hn::ShiftRight<48>(hn::BitCast(du, Rcp));
    K = hn::IfThenElse(large, hn::Sub(K, hn::Set(d, kLargeK)), K);
    const auto Tl = em::Lookup16(d, kLogTl, idx);
    auto Th = em::Lookup16(d, kLogTh, idx);
    const auto R = hn::MulAdd(Xinl, Rcp, Rh0);
    const auto Tlx = hn::MulAdd(Xinl, Rcp, Tl);
    Th = hn::NegMulAdd(hn::Set(d, kL2H), K, Th);
    auto p = hn::MulAdd(hn::Set(d, kLogC9), R, hn::Set(d, kLogC8));
    p = hn::MulAdd(R, p, hn::Set(d, kLogC7));
    p = hn::MulAdd(R, p, hn::Set(d, kLogC6));
    const auto Lo = hn::NegMulAdd(hn::Set(d, kL2L), K, Tlx);
    p = hn::MulAdd(R, p, hn::Set(d, kLogC5));
    const auto Hi = hn::Add(Rh0, Th);
    p = hn::MulAdd(R, p, hn::Set(d, kLogC4));
    const auto nRh = hn::Sub(Th, Hi);
    p = hn::MulAdd(R, p, hn::Set(d, kLogC3));
    const auto R2 = hn::Mul(R, R);
    p = hn::MulAdd(R, p, hn::Set(d, kLogC2));
    const auto Rl = hn::Add(Rh0, nRh);
    p = hn::MulAdd(R, p, hn::Set(d, kLogC1));
    p = hn::MulAdd(R2, p, Lo);
    p = hn::Add(p, Rl);
    return hn::IfThenElse(notsmall, hn::Add(Hi, p), res);
}

void AcoshArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    const auto one = hn::Set(d, kOne);
    const auto maxx = hn::Set(d, kMaxX);
    size_t i = 0;
    // two independent vectors per iteration: one evaluation is a long dependency chain (latency-bound)
    for (; i + 2 * N <= len; i += 2 * N) {
        const auto v0 = hn::LoadU(d, x + i);
        const auto v1 = hn::LoadU(d, x + i + N);
        hn::StoreU(AcoshCore(d, v0), d, y + i);
        hn::StoreU(AcoshCore(d, v1), d, y + i + N);
        if (HWY_UNLIKELY(!hn::AllTrue(d, hn::And(hn::And(hn::Gt(v0, one), hn::Le(v0, maxx)), hn::And(hn::Gt(v1, one), hn::Le(v1, maxx))))))
            FixSpecial(x + i, y + i, 2 * N);
    }
    for (; i + N <= len; i += N) {
        const auto v = hn::LoadU(d, x + i);
        hn::StoreU(AcoshCore(d, v), d, y + i);
        if (HWY_UNLIKELY(!hn::AllTrue(d, hn::And(hn::Gt(v, one), hn::Le(v, maxx))))) FixSpecial(x + i, y + i, N);
    }
    if (i < len) {
        const size_t rem = len - i;
        const auto v = hn::LoadN(d, x + i, rem);
        hn::StoreN(AcoshCore(d, v), d, y + i, rem);
        FixSpecial(x + i, y + i, rem);
    }
}

#else  // no fused multiply-add: the algorithm's error-free steps would not be exact

void AcoshArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    for (size_t i = 0; i < len; i++) y[i] = std::acosh(x[i]);
}

#endif

size_t LanesF64() { return hn::Lanes(hn::ScalableTag<double>()); }

int Variant() {  // 1: SVML algorithm, 0: std::acosh fallback
    return HWY_NATIVE_FMA && HWY_HAVE_FLOAT64;
}

}  // namespace HWY_NAMESPACE
}  // namespace hwyacosh
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace hwyacosh {
HWY_EXPORT(AcoshArray);
HWY_EXPORT(Variant);
HWY_EXPORT(LanesF64);
// dispatch inside the namespace: with a single (static) target, HWY_DYNAMIC_DISPATCH(F) becomes N_<TARGET>::F
static void DispatchAcosh(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(AcoshArray)(x, y, len); }
static int DispatchVariant() { return HWY_DYNAMIC_DISPATCH(Variant)(); }
static size_t DispatchLanes() { return HWY_DYNAMIC_DISPATCH(LanesF64)(); }
}  // namespace hwyacosh

extern "C" {
void hwy_acosh_array(const double *x, double *y, size_t len) { hwyacosh::DispatchAcosh(x, y, len); }
int hwy_acosh_variant(void) { return hwyacosh::DispatchVariant(); }
size_t hwy_acosh_lanes(void) { return hwyacosh::DispatchLanes(); }
}
#endif  // HWY_ONCE
