// Portable (Google Highway) version of Intel SVML's AVX-512 exp2 (svml_z0_exp2_d_ha.s, the algorithm
// of svml/svml_exp2_d_ha.c), with run-time dispatch to the best target. On every target with true FMA
// (HWY_NATIVE_FMA) the results are bit-identical to Intel's AVX-512 code; the AVX-512-only operations
// are emulated exactly:
//   - vreducepd(x, 0x41) = x - floor(16 x)/16, subtraction rounded toward -inf: Floor, then a
//     round-to-nearest subtraction whose exact error (Fast2Sum) decides a 1-ulp downward correction;
//   - vaddpd {rd-sae} x + 1.5*2^48 (only its low 4 bits, floor(16 x) mod 16, are used) and
//     vscalefpd(y, x) = y*2^floor(x): both from the integer floor(16 x), placed in the mantissa of a
//     shifter; 2^floor(x) is a normal number and y*2^floor(x) is exact for all non-callout lanes;
//   - vpermt2pd (16-entry tables): TableLookupLanes / TwoTablesLookupLanes (+ selects on 4/2-lane
//     vectors), no gathers.
// Targets without FMA (SSE2/SSSE3/SSE4) use std::exp2. Lanes with |x| >= 1021.9995 or NaN use the
// scalar SVML callout (svml_dexp2_ha_cout_rare; for NaN it returns x*x, the quieted NaN that the
// AVX-512 code returns). Also exports Highway's own hwy/contrib/math Exp2 for comparison.
// Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.

#include <cmath>
#include <cstddef>
#include <cstdint>

extern "C" int svml_dexp2_ha_cout_rare(const double *px, double *pr);

#include "functions/exp2/hwy/hwy_exp2_constants.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "functions/exp2/hwy/hwy_exp2.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"
#include "hwy/contrib/math/math-inl.h"
#include "hwymath/hwy/avx512_native-inl.h"

HWY_BEFORE_NAMESPACE();
namespace hwyexp2 {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64

// scalar fix-up for lanes the vector path does not handle: |x| >= threshold or NaN (rare: kept out of line)
static HWY_NOINLINE __attribute__((cold)) void FixSpecial(const double *x, double *y, size_t n) {
    for (size_t k = 0; k < n; k++) {
        if (!(std::fabs(x[k]) < kThreshold)) svml_dexp2_ha_cout_rare(&x[k], &y[k]);
    }
}

// Th[j], Tl[j] for j < 16 (vpermt2pd on AVX-512): one vector for >= 16 lanes, two for 8 lanes,
// two pairs of vectors + a select for 4 lanes, four pairs + two selects for 2 lanes.
template <class D, class VU>
HWY_INLINE void Lookup16(D d, VU j, hn::VFromD<D> &th, hn::VFromD<D> &tl) {
    const hn::RebindToSigned<D> di;
    const size_t N = hn::Lanes(d);
    const auto ji = hn::BitCast(di, j);
    if (N >= 16 && N <= 64) {
        const auto idx = hn::IndicesFromVec(d, ji);
        th = hn::TableLookupLanes(hn::LoadU(d, ThPad()), idx);
        tl = hn::TableLookupLanes(hn::LoadU(d, TlPad()), idx);
    } else if (N == 8) {
        const auto idx = hn::IndicesFromVec(d, ji);
        th = hn::TwoTablesLookupLanes(d, hn::LoadU(d, vTh), hn::LoadU(d, vTh + 8), idx);
        tl = hn::TwoTablesLookupLanes(d, hn::LoadU(d, vTl), hn::LoadU(d, vTl + 8), idx);
    } else if (N == 4) {
        const auto idx = hn::IndicesFromVec(d, hn::And(ji, hn::Set(di, 7)));
        const auto hi = hn::RebindMask(d, hn::TestBit(ji, hn::Set(di, 8)));
        th = hn::IfThenElse(hi, hn::TwoTablesLookupLanes(d, hn::LoadU(d, vTh + 8), hn::LoadU(d, vTh + 12), idx),
                            hn::TwoTablesLookupLanes(d, hn::LoadU(d, vTh), hn::LoadU(d, vTh + 4), idx));
        tl = hn::IfThenElse(hi, hn::TwoTablesLookupLanes(d, hn::LoadU(d, vTl + 8), hn::LoadU(d, vTl + 12), idx),
                            hn::TwoTablesLookupLanes(d, hn::LoadU(d, vTl), hn::LoadU(d, vTl + 4), idx));
    } else if (N == 2) {
        const auto idx = hn::IndicesFromVec(d, hn::And(ji, hn::Set(di, 3)));
        const auto b2 = hn::RebindMask(d, hn::TestBit(ji, hn::Set(di, 4)));
        const auto b3 = hn::RebindMask(d, hn::TestBit(ji, hn::Set(di, 8)));
        auto quad = [&](const double *t) {  // t[0..3]
            return hn::TwoTablesLookupLanes(d, hn::LoadU(d, t), hn::LoadU(d, t + 2), idx);
        };
        th = hn::IfThenElse(b3, hn::IfThenElse(b2, quad(vTh + 12), quad(vTh + 8)),
                            hn::IfThenElse(b2, quad(vTh + 4), quad(vTh)));
        tl = hn::IfThenElse(b3, hn::IfThenElse(b2, quad(vTl + 12), quad(vTl + 8)),
                            hn::IfThenElse(b2, quad(vTl + 4), quad(vTl)));
    } else {  // e.g. single-lane targets
        th = hn::GatherIndex(d, vTh, ji);
        tl = hn::GatherIndex(d, vTl, ji);
    }
}

template <class D>
HWY_INLINE hn::VFromD<D> Exp2Core(D d, hn::VFromD<D> x) {
    const hn::RebindToUnsigned<D> du;
#if HWYMATH_NATIVE_AVX512
    // Intel's instructions: vreducepd 0x41 {sae}, vaddpd {rd-sae}, vpermt2pd, vscalefpd
    namespace nat = hwymath::HWY_NAMESPACE;
    const auto r = nat::Reduce<0x41>(x);
    const auto j = hn::And(hn::BitCast(du, nat::AddR<nat::kRD>(x, hn::Set(d, kShifter))), hn::Set(du, 15));
#else
    // vreducepd(x, 0x41): r = x - floor(16 x)/16, rounded toward -inf. f <= x and the subtraction is
    // exact unless -1/16 < x < 0 (f = -1/16, |f| > |x|), so err = x - (r + f) is its exact error
    // (Fast2Sum); err < 0 means that RN rounded up: step one ulp down (r > 0 there).
    const auto f16 = hn::Floor(hn::Mul(x, hn::Set(d, 16.0)));
    const auto f = hn::Mul(f16, hn::Set(d, 0.0625));
    const auto rn = hn::Sub(x, f);
    const auto err = hn::Sub(x, hn::Add(rn, f));
    const auto r = hn::BitCast(d, hn::Add(hn::BitCast(du, rn),
                                          hn::BitCast(du, hn::VecFromMask(d, hn::Lt(err, hn::Zero(d))))));

    // index j = floor(16 x) mod 16 (the low bits of x + 1.5*2^48 rounded down) and 2^floor(x)
    const auto tb = hn::BitCast(du, hn::Add(f16, hn::Set(d, kShifter52)));
    const auto j = hn::And(tb, hn::Set(du, 15));
    const auto scale = hn::BitCast(d, hn::ShiftLeft<52>(hn::ShiftRight<4>(tb)));
#endif

    auto p = hn::MulAdd(hn::Set(d, kC6), r, hn::Set(d, kC5));
    p = hn::MulAdd(r, p, hn::Set(d, kC4));
    p = hn::MulAdd(r, p, hn::Set(d, kC3));
    p = hn::MulAdd(r, p, hn::Set(d, kC2));
    p = hn::MulAdd(r, p, hn::Set(d, kC1));
    hn::VFromD<D> th, tl;
    Lookup16(d, j, th, tl);
    p = hn::MulAdd(hn::Mul(th, r), p, tl);  // Tl + Th*r*P(r)
    const auto y = hn::Add(p, th);
#if HWYMATH_NATIVE_AVX512
    return hwymath::HWY_NAMESPACE::Scalef(y, x);
#else
    return hn::Mul(y, scale);  // vscalefpd(y, x): exact, the result is normal
#endif
}

void Exp2Array(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    const auto thr = hn::Set(d, kThreshold);
    size_t i = 0;
    for (; i + N <= len; i += N) {
        const auto v = hn::LoadU(d, x + i);
        hn::StoreU(Exp2Core(d, v), d, y + i);
        if (HWY_UNLIKELY(!hn::AllFalse(d, hn::Not(hn::Lt(hn::Abs(v), thr))))) FixSpecial(x + i, y + i, N);
    }
    if (i < len) {
        const size_t rem = len - i;
        const auto v = hn::LoadN(d, x + i, rem);
        hn::StoreN(Exp2Core(d, v), d, y + i, rem);
        FixSpecial(x + i, y + i, rem);
    }
}

#else  // no fused multiply-add: the algorithm would not be reproduced exactly

void Exp2Array(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    for (size_t i = 0; i < len; i++) y[i] = std::exp2(x[i]);
}

#endif

// Highway's own exp2 (hwy/contrib/math), for comparison
void ContribExp2Array(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
#if HWY_HAVE_FLOAT64
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    size_t i = 0;
    for (; i + N <= len; i += N) hn::StoreU(hn::Exp2(d, hn::LoadU(d, x + i)), d, y + i);
    if (i < len) hn::StoreN(hn::Exp2(d, hn::LoadN(d, x + i, len - i)), d, y + i, len - i);
#else
    for (size_t i = 0; i < len; i++) y[i] = std::exp2(x[i]);
#endif
}

size_t LanesF64() { return hn::Lanes(hn::ScalableTag<double>()); }

int Variant() {  // 1: SVML algorithm (bit-identical), 0: std::exp2 fallback
#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64
    return 1;
#else
    return 0;
#endif
}

}  // namespace HWY_NAMESPACE
}  // namespace hwyexp2
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace hwyexp2 {
HWY_EXPORT(Exp2Array);
HWY_EXPORT(ContribExp2Array);
HWY_EXPORT(Variant);
HWY_EXPORT(LanesF64);
// dispatch inside the namespace: with a single (static) target, HWY_DYNAMIC_DISPATCH(F) becomes N_<TARGET>::F
static void DispatchExp2(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(Exp2Array)(x, y, len); }
static void DispatchContrib(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(ContribExp2Array)(x, y, len); }
static int DispatchVariant() { return HWY_DYNAMIC_DISPATCH(Variant)(); }
static size_t DispatchLanes() { return HWY_DYNAMIC_DISPATCH(LanesF64)(); }
}  // namespace hwyexp2

extern "C" {
void hwy_exp2_array(const double *x, double *y, size_t len) { hwyexp2::DispatchExp2(x, y, len); }
void hwy_contrib_exp2_array(const double *x, double *y, size_t len) { hwyexp2::DispatchContrib(x, y, len); }
int hwy_exp2_variant(void) { return hwyexp2::DispatchVariant(); }
size_t hwy_exp2_lanes(void) { return hwyexp2::DispatchLanes(); }
}
#endif  // HWY_ONCE
