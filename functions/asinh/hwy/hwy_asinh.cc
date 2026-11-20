// Portable (Google Highway) version of Intel SVML's AVX-512 asinh (svml/svml_asinh_d_ha.c), with run-time
// dispatch to the best target. Per target:
//   - targets with true FMA (HWY_NATIVE_FMA): the SVML algorithm, operation for operation, with the
//     AVX-512-only instructions emulated exactly (hwy/avx512_emul-inl.h: VRCP14PD, VRSQRT14PD, VGETEXPPD,
//     16-entry VPERMT2PD lookups; masked operations as selects), so results are bit-identical to Intel's
//     __svml_asinh8_ha on every such target;
//   - targets without FMA (SSE2/SSSE3/SSE4): std::asinh.
// Lanes with |x| > DBL_MAX (inf) or NaN use the scalar SVML callout (svml_dasinh_ha_cout_rare).
// Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.

#include <cmath>
#include <cstddef>
#include <cstdint>

extern "C" int svml_dasinh_ha_cout_rare(const double *px, double *pr);

#ifndef HWYASINH_CONSTANTS
#define HWYASINH_CONSTANTS
namespace hwyasinh {
#include "functions/asinh/svml/svml_asinh_d_ha_constants.inc"
}  // namespace hwyasinh
#endif

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "functions/asinh/hwy/hwy_asinh.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"
#include "functions/asinh/hwy/avx512_emul-inl.h"

HWY_BEFORE_NAMESPACE();
namespace hwyasinh {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64
namespace em = hwymath::HWY_NAMESPACE;

// scalar fix-up for |x| > DBL_MAX or NaN (rare: kept out of line)
static HWY_NOINLINE __attribute__((cold)) void FixSpecial(const double *x, double *y, size_t n) {
    for (size_t k = 0; k < n; k++) {
        if (!(std::fabs(x[k]) <= kMaxX)) svml_dasinh_ha_cout_rare(&x[k], &y[k]);
    }
}

// The operations of svml_asinh8_ha() (same names, same order where it matters)
template <class D>
HWY_INLINE hn::VFromD<D> AsinhCore(D d, hn::VFromD<D> x) {
    const hn::RebindToUnsigned<D> du;
    const auto one = hn::Set(d, kOne);
    const auto x2 = hn::Mul(x, x);
    const auto A = hn::Max(one, x2);
    auto B = hn::Min(one, x2);
    auto ps = hn::MulAdd(hn::Set(d, kSmallC3), x2, hn::Set(d, kSmallC2));
    const auto Y = hn::Add(x2, one);
    const auto ax = hn::Abs(x);
    ps = hn::MulAdd(x2, ps, hn::Set(d, kSmallC1));
    const auto R0 = em::Rsqrt14(d, Y);
    const auto notsmall = hn::Not(hn::Lt(ax, hn::Set(d, kSmallX)));
    const auto Bh = hn::Sub(Y, A);
    const auto sign = hn::Xor(ax, x);
    const auto Sh = hn::Mul(Y, R0);
    auto res = hn::Mul(ps, x2);
    B = hn::Sub(B, Bh);
    const auto YR0l = hn::MulSub(R0, Y, Sh);
    res = hn::MulAdd(ax, res, ax);
    const auto x2l = hn::MulSub(x, x, x2);
    auto Sl = hn::Add(B, x2l);
    Sl = hn::MulAdd(R0, Sl, YR0l);
    const auto large = hn::Not(hn::Lt(ax, hn::Set(d, kLargeX)));
    auto Eh = hn::NegMulAdd(R0, Sh, one);
    Eh = hn::NegMulAdd(R0, Sl, Eh);
    const auto ShEh = hn::Mul(Sh, Eh);
    auto pS = hn::MulAdd(hn::Set(d, kSqrtC4), Eh, hn::Set(d, kSqrtC3));
    pS = hn::MulAdd(Eh, pS, hn::Set(d, kSqrtC2));
    pS = hn::MulAdd(Eh, pS, hn::Set(d, kSqrtC1));
    const auto Shx = hn::Add(Sh, ax);
    auto Sc = hn::MulAdd(pS, ShEh, Sl);
    const auto Xh = hn::Sub(Shx, Sh);
    auto Xin = hn::Add(Shx, Sc);
    const auto Xl = hn::Sub(ax, Xh);
    const auto Slh = hn::Sub(Xin, Shx);
    Xin = hn::IfThenElse(large, hn::Mul(ax, hn::Set(d, kLargeScale)), Xin);
    Sc = hn::Sub(Sc, Slh);
    const auto rcp = em::Rcp14(d, Xin);
    auto Xinl = hn::Add(Sc, Xl);
    Xinl = hn::IfThenZeroElse(large, Xinl);
    const auto Rcp = hn::BitCast(d, hn::And(hn::Add(hn::BitCast(du, rcp), hn::Set(du, kRcpRoundAdd)),
                                            hn::Set(du, kRcpRoundMask)));
    const auto idx = hn::ShiftRight<48>(hn::BitCast(du, Rcp));
    const auto Rh0 = hn::MulSub(Xin, Rcp, one);
    auto K = em::GetExpPositiveNormal(d, Rcp);
    const auto Tl = em::Lookup16(d, kLogTl, idx);
    K = hn::IfThenElse(large, hn::Sub(K, hn::Set(d, kLargeK)), K);
    auto Th = em::Lookup16(d, kLogTh, idx);
    const auto R = hn::MulAdd(Xinl, Rcp, Rh0);
    const auto Tlx = hn::MulAdd(Xinl, Rcp, Tl);
    auto p = hn::MulAdd(hn::Set(d, kLogC9), R, hn::Set(d, kLogC8));
    Th = hn::NegMulAdd(hn::Set(d, kL2H), K, Th);
    p = hn::MulAdd(R, p, hn::Set(d, kLogC7));
    const auto Hi = hn::Add(Rh0, Th);
    p = hn::MulAdd(R, p, hn::Set(d, kLogC6));
    const auto Lo = hn::NegMulAdd(hn::Set(d, kL2L), K, Tlx);
    const auto nRh = hn::Sub(Th, Hi);
    p = hn::MulAdd(R, p, hn::Set(d, kLogC5));
    const auto Rl = hn::Add(Rh0, nRh);
    p = hn::MulAdd(R, p, hn::Set(d, kLogC4));
    const auto R2 = hn::Mul(R, R);
    p = hn::MulAdd(R, p, hn::Set(d, kLogC3));
    p = hn::MulAdd(R, p, hn::Set(d, kLogC2));
    p = hn::MulAdd(R, p, hn::Set(d, kLogC1));
    p = hn::MulAdd(R2, p, Lo);
    p = hn::Add(p, Rl);
    res = hn::IfThenElse(notsmall, hn::Add(Hi, p), res);
    return hn::Xor(res, sign);
}

void AsinhArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    const auto maxx = hn::Set(d, kMaxX);
    size_t i = 0;
    // two independent vectors per iteration: one evaluation is a long dependency chain (latency-bound)
    for (; i + 2 * N <= len; i += 2 * N) {
        const auto v0 = hn::LoadU(d, x + i);
        const auto v1 = hn::LoadU(d, x + i + N);
        hn::StoreU(AsinhCore(d, v0), d, y + i);
        hn::StoreU(AsinhCore(d, v1), d, y + i + N);
        if (HWY_UNLIKELY(!hn::AllTrue(d, hn::And(hn::Le(hn::Abs(v0), maxx), hn::Le(hn::Abs(v1), maxx)))))
            FixSpecial(x + i, y + i, 2 * N);
    }
    for (; i + N <= len; i += N) {
        const auto v = hn::LoadU(d, x + i);
        hn::StoreU(AsinhCore(d, v), d, y + i);
        if (HWY_UNLIKELY(!hn::AllTrue(d, hn::Le(hn::Abs(v), maxx)))) FixSpecial(x + i, y + i, N);
    }
    if (i < len) {
        const size_t rem = len - i;
        const auto v = hn::LoadN(d, x + i, rem);
        hn::StoreN(AsinhCore(d, v), d, y + i, rem);
        FixSpecial(x + i, y + i, rem);
    }
}

#else  // no fused multiply-add: the algorithm's error-free steps would not be exact

void AsinhArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    for (size_t i = 0; i < len; i++) y[i] = std::asinh(x[i]);
}

#endif

size_t LanesF64() { return hn::Lanes(hn::ScalableTag<double>()); }

int Variant() {  // 1: SVML algorithm, 0: std::asinh fallback
    return HWY_NATIVE_FMA && HWY_HAVE_FLOAT64;
}

}  // namespace HWY_NAMESPACE
}  // namespace hwyasinh
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace hwyasinh {
HWY_EXPORT(AsinhArray);
HWY_EXPORT(Variant);
HWY_EXPORT(LanesF64);
// dispatch inside the namespace: with a single (static) target, HWY_DYNAMIC_DISPATCH(F) becomes N_<TARGET>::F
static void DispatchAsinh(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(AsinhArray)(x, y, len); }
static int DispatchVariant() { return HWY_DYNAMIC_DISPATCH(Variant)(); }
static size_t DispatchLanes() { return HWY_DYNAMIC_DISPATCH(LanesF64)(); }
}  // namespace hwyasinh

extern "C" {
void hwy_asinh_array(const double *x, double *y, size_t len) { hwyasinh::DispatchAsinh(x, y, len); }
int hwy_asinh_variant(void) { return hwyasinh::DispatchVariant(); }
size_t hwy_asinh_lanes(void) { return hwyasinh::DispatchLanes(); }
}
#endif  // HWY_ONCE
