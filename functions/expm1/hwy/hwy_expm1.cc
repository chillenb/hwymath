// Portable (Google Highway) version of Intel SVML's AVX-512 expm1 (svml_z0_expm1_d_ha.s, the algorithm
// of svml/svml_expm1_d_ha.c), with run-time dispatch to the best target. On every target with true FMA
// (HWY_NATIVE_FMA) the results are bit-identical to Intel's AVX-512 code; the AVX-512-only operations
// are emulated exactly:
//   - vscalefpd(Th, max(N, -128)) = Th*2^floor(max(N, -128)): 2^floor(...) from the bits of
//     max(N, -128) + shifter (a multiple of 1/16 in the shifter's mantissa); the product is a normal
//     number, so the multiplication is exact like vscalefpd;
//   - vpermt2pd (16-entry tables): TableLookupLanes / TwoTablesLookupLanes (+ selects on 4/2-lane
//     vectors), no gathers;
//   - vblendmpd / vmaxpd / bitwise masks: IfThenElse / Max / And / Or.
// On AVX-512 targets (HWYMATH_NATIVE_AVX512) the instructions themselves are used instead: vpermt2pd on
// the raw bits of t (no index masking), vscalefpd, and one vcmppd (NLT_UQ) for the callout lanes.
// Everything else is the same sequence of round-to-nearest FMAs, multiplications and additions.
// Targets without FMA (SSE2/SSSE3/SSE4) use std::expm1. Lanes with x >= 707.703..., +inf or NaN use
// the scalar SVML callout (svml_dexpm1_ha_cout_rare). Also exports Highway's own hwy/contrib/math Expm1
// for comparison.
// Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.

#include <cmath>
#include <cstddef>
#include <cstdint>

extern "C" int svml_dexpm1_ha_cout_rare(const double *px, double *pr);

#include "functions/expm1/hwy/hwy_expm1_constants.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "functions/expm1/hwy/hwy_expm1.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"
#include "hwy/contrib/math/math-inl.h"
#include "hwymath/hwy/avx512_native-inl.h"

HWY_BEFORE_NAMESPACE();
namespace hwyexpm1 {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64

// scalar fix-up for lanes the vector path does not handle: x >= threshold or NaN (rare: kept out of line)
static HWY_NOINLINE __attribute__((cold)) void FixSpecial(const double *x, double *y, size_t n) {
    for (size_t k = 0; k < n; k++) {
        if (!(x[k] < kThreshold)) svml_dexpm1_ha_cout_rare(&x[k], &y[k]);
    }
}

// a[j], b[j] for j < 16 from two 16-entry tables (vpermt2pd on AVX-512; pa/pb: tables repeated to 64
// entries): one vector for >= 16 lanes, two for 8 lanes, two pairs of vectors + a select for 4 lanes,
// four pairs + two selects for 2 lanes.
template <class D, class VU>
HWY_INLINE void Lookup16(D d, VU j, const double *ta, const double *pa, const double *tb, const double *pb,
                         hn::VFromD<D> &a, hn::VFromD<D> &b) {
    const hn::RebindToSigned<D> di;
    const size_t N = hn::Lanes(d);
    const auto ji = hn::BitCast(di, j);
    if (N >= 16 && N <= 64) {
        const auto idx = hn::IndicesFromVec(d, ji);
        a = hn::TableLookupLanes(hn::LoadU(d, pa), idx);
        b = hn::TableLookupLanes(hn::LoadU(d, pb), idx);
    } else if (N == 8) {
        const auto idx = hn::IndicesFromVec(d, ji);
        a = hn::TwoTablesLookupLanes(d, hn::LoadU(d, ta), hn::LoadU(d, ta + 8), idx);
        b = hn::TwoTablesLookupLanes(d, hn::LoadU(d, tb), hn::LoadU(d, tb + 8), idx);
    } else if (N == 4) {
        const auto idx = hn::IndicesFromVec(d, hn::And(ji, hn::Set(di, 7)));
        const auto hi = hn::RebindMask(d, hn::TestBit(ji, hn::Set(di, 8)));
        auto oct = [&](const double *t) { return hn::TwoTablesLookupLanes(d, hn::LoadU(d, t), hn::LoadU(d, t + 4), idx); };
        a = hn::IfThenElse(hi, oct(ta + 8), oct(ta));
        b = hn::IfThenElse(hi, oct(tb + 8), oct(tb));
    } else if (N == 2) {
        const auto idx = hn::IndicesFromVec(d, hn::And(ji, hn::Set(di, 3)));
        const auto b2 = hn::RebindMask(d, hn::TestBit(ji, hn::Set(di, 4)));
        const auto b3 = hn::RebindMask(d, hn::TestBit(ji, hn::Set(di, 8)));
        auto quad = [&](const double *t) { return hn::TwoTablesLookupLanes(d, hn::LoadU(d, t), hn::LoadU(d, t + 2), idx); };
        a = hn::IfThenElse(b3, hn::IfThenElse(b2, quad(ta + 12), quad(ta + 8)), hn::IfThenElse(b2, quad(ta + 4), quad(ta)));
        b = hn::IfThenElse(b3, hn::IfThenElse(b2, quad(tb + 12), quad(tb + 8)), hn::IfThenElse(b2, quad(tb + 4), quad(tb)));
    } else {  // e.g. single-lane targets
        a = hn::GatherIndex(d, ta, ji);
        b = hn::GatherIndex(d, tb, ji);
    }
}

#if HWYMATH_NATIVE_AVX512
// vpermt2pd with the raw index vector (the instruction uses only its low 4 bits; Highway's
// TwoTablesLookupLanes requires in-range indices, which debug builds assert)
// lanes with !(x < thr) (x >= thr or NaN): one vcmppd NLT_UQ, like Intel's code
#endif

template <class D>
HWY_INLINE hn::VFromD<D> AndBits(D d, hn::VFromD<D> v, uint64_t bits) {
    const hn::RebindToUnsigned<D> du;
    return hn::BitCast(d, hn::And(hn::BitCast(du, v), hn::Set(du, bits)));
}

template <class D>
HWY_INLINE hn::VFromD<D> Expm1Core(D d, hn::VFromD<D> x) {
    const hn::RebindToUnsigned<D> du;
    const auto S = hn::Set(d, kShifter);
    const auto t = hn::MulAdd(x, hn::Set(d, kLog2e), S);  // N = t - S: x*log2(e) to a multiple of 1/16
    const auto n = hn::Sub(t, S);
#if HWYMATH_NATIVE_AVX512
    const auto tb = hn::BitCast(du, t);  // table index in the low 4 bits
    auto th = hwymath::HWY_NAMESPACE::Permute2Var(hn::LoadU(d, vTh), tb, hn::LoadU(d, vTh + 8));
    const auto tlr = hwymath::HWY_NAMESPACE::Permute2Var(hn::LoadU(d, vTlr), tb, hn::LoadU(d, vTlr + 8));
    const auto sgn = AndBits(d, x, kSignMaskBits);
    const auto r = hn::NegMulAdd(n, hn::Set(d, kLn2Hi), x);  // R = x - N*ln2_hi
    th = hwymath::HWY_NAMESPACE::Scalef(th, hn::Max(n, hn::Set(d, kMinN)));  // Th' = Th*2^floor(max(N, -128))
#else
    hn::VFromD<D> th, tlr;
    Lookup16(d, hn::And(hn::BitCast(du, t), hn::Set(du, 15)), vTh, ThPad(), vTlr, TlrPad(), th, tlr);
    const auto sgn = AndBits(d, x, kSignMaskBits);
    const auto r = hn::NegMulAdd(n, hn::Set(d, kLn2Hi), x);  // R = x - N*ln2_hi
    // vscalefpd(Th, max(N, -128)): max(N, -128) + S holds 16*(N + 1023) in its low mantissa bits
    const auto tc = hn::BitCast(du, hn::Add(hn::Max(n, hn::Set(d, kMinN)), S));
    th = hn::Mul(th, hn::BitCast(d, hn::ShiftLeft<52>(hn::ShiftRight<4>(tc))));
#endif
    auto rr = hn::NegMulAdd(n, hn::Set(d, kLn2Lo), r);  // R' = R - N*ln2_lo
    const auto small = hn::Lt(th, hn::Set(d, 1.0));
    auto rl = hn::Sub(r, rr);
    rr = AndBits(d, rr, kRMaskBits);  // |R'| < 2 even for special cases
    rl = hn::NegMulAdd(hn::Set(d, kLn2Lo), n, rl);  // Rl = -(N*ln2_lo)_low
    rl = AndBits(d, rl, kRMaskBits);
    const auto p56 = hn::MulAdd(hn::Set(d, kC6), rr, hn::Set(d, kC5));
    const auto p34 = hn::MulAdd(hn::Set(d, kC4), rr, hn::Set(d, kC3));
    const auto r2 = hn::Mul(rr, rr);
    const auto lo1 = hn::MulAdd(rr, tlr, hn::Add(rl, tlr));  // Rl + Tlr + Tlr*R
    auto p = hn::MulAdd(hn::Set(d, kC8), rr, hn::Set(d, kC7));
    p = hn::MulAdd(r2, p, p56);
    p = hn::MulAdd(r2, p, p34);
    p = hn::MulAdd(rr, p, hn::Set(d, kC2));
    p = hn::MulAdd(r2, p, lo1);  // R^2 P(R) + Tlr*(1+R) + Rl
    // Th' - 1 by Fast2Sum (larger magnitude first) and Th'*R as a double-double
    const auto m1 = hn::Set(d, -1.0);
    const auto tm1 = hn::Add(th, m1);
    const auto tl1 = hn::Add(hn::Sub(hn::IfThenElse(small, m1, th), tm1), hn::IfThenElse(small, th, m1));
    const auto hi = hn::MulAdd(rr, th, tm1);
    const auto lo = hn::MulSub(rr, th, hn::Sub(hi, tm1));
    const auto res = hn::Add(hi, hn::Add(hn::MulAdd(p, th, tl1), lo));
    return hn::Or(res, sgn);
}

void Expm1Array(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    const auto thr = hn::Set(d, kThreshold);
    size_t i = 0;
    for (; i + N <= len; i += N) {
        const auto v = hn::LoadU(d, x + i);
        hn::StoreU(Expm1Core(d, v), d, y + i);
#if HWYMATH_NATIVE_AVX512
        const auto special = hwymath::HWY_NAMESPACE::Cmp<_CMP_NLT_UQ>(v, thr);
#else
        const auto special = hn::Not(hn::Lt(v, thr));
#endif
        if (HWY_UNLIKELY(!hn::AllFalse(d, special))) FixSpecial(x + i, y + i, N);
    }
    if (i < len) {
        const size_t rem = len - i;
        const auto v = hn::LoadN(d, x + i, rem);
        hn::StoreN(Expm1Core(d, v), d, y + i, rem);
        FixSpecial(x + i, y + i, rem);
    }
}

#else  // no fused multiply-add: the algorithm would not be reproduced exactly

void Expm1Array(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    for (size_t i = 0; i < len; i++) y[i] = std::expm1(x[i]);
}

#endif

// Highway's own expm1 (hwy/contrib/math), for comparison
void ContribExpm1Array(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
#if HWY_HAVE_FLOAT64
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    size_t i = 0;
    for (; i + N <= len; i += N) hn::StoreU(hn::Expm1(d, hn::LoadU(d, x + i)), d, y + i);
    if (i < len) hn::StoreN(hn::Expm1(d, hn::LoadN(d, x + i, len - i)), d, y + i, len - i);
#else
    for (size_t i = 0; i < len; i++) y[i] = std::expm1(x[i]);
#endif
}

size_t LanesF64() { return hn::Lanes(hn::ScalableTag<double>()); }

int Variant() {  // 1: SVML algorithm (bit-identical), 0: std::expm1 fallback
#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64
    return 1;
#else
    return 0;
#endif
}

}  // namespace HWY_NAMESPACE
}  // namespace hwyexpm1
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace hwyexpm1 {
HWY_EXPORT(Expm1Array);
HWY_EXPORT(ContribExpm1Array);
HWY_EXPORT(Variant);
HWY_EXPORT(LanesF64);
// dispatch inside the namespace: with a single (static) target, HWY_DYNAMIC_DISPATCH(F) becomes N_<TARGET>::F
static void DispatchExpm1(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(Expm1Array)(x, y, len); }
static void DispatchContrib(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(ContribExpm1Array)(x, y, len); }
static int DispatchVariant() { return HWY_DYNAMIC_DISPATCH(Variant)(); }
static size_t DispatchLanes() { return HWY_DYNAMIC_DISPATCH(LanesF64)(); }
}  // namespace hwyexpm1

extern "C" {
void hwy_expm1_array(const double *x, double *y, size_t len) { hwyexpm1::DispatchExpm1(x, y, len); }
void hwy_contrib_expm1_array(const double *x, double *y, size_t len) { hwyexpm1::DispatchContrib(x, y, len); }
int hwy_expm1_variant(void) { return hwyexpm1::DispatchVariant(); }
size_t hwy_expm1_lanes(void) { return hwyexpm1::DispatchLanes(); }
}
#endif  // HWY_ONCE
