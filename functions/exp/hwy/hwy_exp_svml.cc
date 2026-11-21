// Portable (Google Highway) version of Intel SVML's AVX-512 high-accuracy exp (svml_z0_exp_d_ha.s, the
// algorithm of svml/svml_exp_d_ha.c), with run-time dispatch to the best target. Not to be confused with
// hwy_exp.cc, the Highway version of this repository's own, more accurate hp algorithm.
//
// On every target with true FMA (HWY_NATIVE_FMA) the results are bit-identical to Intel's __svml_exp8_ha;
// the AVX-512-only operations are emulated exactly by the small helpers below (each one can be replaced
// by the native instruction on AVX-512 targets):
//   - vfmadd {rz-sae}: t = x*log2(e) + S (S = 1.5*2^48 + 1023/16) rounded toward zero (ShifterRZ):
//     the round-to-nearest FMA plus an FMA that gives the sign of its exact rounding error (S - t is
//     exact by Sterbenz); where RN rounded up (t > 0), step one ulp down on the bit pattern;
//   - vpermt2pd (16-entry tables, index = low 4 bits of t; Lookup16): two-source permute on 8 lanes,
//     TableLookupLanes of a repeated table for >= 16 lanes, TwoTablesLookupLanes + selects on 4/2-lane
//     vectors; no gathers;
//   - vscalefpd(y, n) = y*2^floor(n) (Scalef): bits 4..15 of t are 1023 + floor(n), so 2^floor(n) is
//     built by two shifts; it is a normal number and y*2^floor(n) is exact for all non-callout lanes.
// The `and` with 0xbfff... that Intel applies to r (it keeps |r| < 2 for special inputs) is a no-op for
// every lane whose vector result is kept, so it is omitted. Lanes with |x| >= 707.703... and NaN lanes
// use the scalar SVML callout (svml_dexp_ha_cout_rare; for NaN it returns x*x, the quieted NaN that the
// AVX-512 code returns). Targets without FMA (SSE2/SSSE3/SSE4) use std::exp.
// Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.

#include <cmath>
#include <cstddef>
#include <cstdint>

extern "C" int svml_dexp_ha_cout_rare(const double *px, double *pr);

#include "functions/exp/hwy/hwy_exp_svml_constants.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "functions/exp/hwy/hwy_exp_svml.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace hwyexpsvml {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64

// scalar fix-up for lanes the vector path does not handle: |x| >= threshold or NaN (rare: kept out of line)
static HWY_NOINLINE __attribute__((cold)) void FixSpecial(const double *x, double *y, size_t n) {
    for (size_t k = 0; k < n; k++) {
        if (!(std::fabs(x[k]) < kThreshold)) svml_dexp_ha_cout_rare(&x[k], &y[k]);
    }
}

// vfmadd {rz-sae}(x, log2(e), S), returned as bits. Exact for all lanes with S/2 < t < 2S (every
// non-callout lane): err = RN(x*log2(e) + (S - t)) has the sign of the exact rounding error, and err < 0
// means that RN rounded up (away from zero, t > 0), so t - 1 ulp is the round-toward-zero result.
template <class D>
HWY_INLINE hn::VFromD<hn::RebindToUnsigned<D>> ShifterRZ(D d, hn::VFromD<D> x) {
    const hn::RebindToUnsigned<D> du;
    const auto L = hn::Set(d, kLog2e);
    const auto S = hn::Set(d, kShifter);
    const auto t = hn::MulAdd(x, L, S);
    const auto err = hn::MulAdd(x, L, hn::Sub(S, t));
    const auto up = hn::Lt(err, hn::Zero(d));
    const auto tb = hn::BitCast(du, t);
#if HWY_ARCH_X86 && HWY_TARGET <= HWY_AVX3
    return hn::MaskedSubOr(tb, hn::RebindMask(du, up), tb, hn::Set(du, 1));  // one masked vpsubq
#else
    return hn::Add(tb, hn::BitCast(du, hn::VecFromMask(d, up)));  // mask lanes are all ones = -1
#endif
}

// vpermt2pd: th = Th[j], tlr = Tlr[j] with j = tb mod 16
template <class D, class VU>
HWY_INLINE void Lookup16(D d, VU tb, hn::VFromD<D> &th, hn::VFromD<D> &tlr) {
    const hn::RebindToSigned<D> di;
    const size_t N = hn::Lanes(d);
#if HWY_ARCH_X86 && HWY_TARGET <= HWY_AVX3
    // 8 lanes: vpermt2pd itself, which only reads the low 4 bits of the index
    if (N == 8) {
        const auto idx = hn::IndicesFromVec(d, hn::BitCast(di, tb));
        th = hn::TwoTablesLookupLanes(d, hn::LoadU(d, vTh), hn::LoadU(d, vTh + 8), idx);
        tlr = hn::TwoTablesLookupLanes(d, hn::LoadU(d, vTlr), hn::LoadU(d, vTlr + 8), idx);
        return;
    }
#elif HWY_TARGET == HWY_AVX2
    // 4 lanes: two-table lookups (vpermd pairs + blendv on bit 2) of entries 0-7 and 8-15, then blendv on
    // bit 3 moved into the sign bit (AVX2 blendv only reads the sign bit; Highway's own AVX2
    // TwoTablesLookupLanes relies on the same). The 32-bit vpermd indices only use bits 0-2 of j, so no
    // masking is needed.
    if (N == 4) {
        const auto idx = hn::IndicesFromVec(d, hn::BitCast(di, tb));
        const auto hi = hn::MaskFromVec(hn::BitCast(d, hn::ShiftLeft<60>(tb)));
        th = hn::IfThenElse(hi, hn::TwoTablesLookupLanes(d, hn::LoadU(d, vTh + 8), hn::LoadU(d, vTh + 12), idx),
                            hn::TwoTablesLookupLanes(d, hn::LoadU(d, vTh), hn::LoadU(d, vTh + 4), idx));
        tlr = hn::IfThenElse(hi, hn::TwoTablesLookupLanes(d, hn::LoadU(d, vTlr + 8), hn::LoadU(d, vTlr + 12), idx),
                             hn::TwoTablesLookupLanes(d, hn::LoadU(d, vTlr), hn::LoadU(d, vTlr + 4), idx));
        return;
    }
#endif
    const auto ji = hn::BitCast(di, hn::And(tb, hn::Set(hn::RebindToUnsigned<D>(), 15)));
    if (N >= 16 && N <= 64) {  // the (repeated) table in one vector
        const auto idx = hn::IndicesFromVec(d, ji);
        th = hn::TableLookupLanes(hn::LoadU(d, ThPad()), idx);
        tlr = hn::TableLookupLanes(hn::LoadU(d, TlrPad()), idx);
    } else if (N == 8) {
        const auto idx = hn::IndicesFromVec(d, ji);
        th = hn::TwoTablesLookupLanes(d, hn::LoadU(d, vTh), hn::LoadU(d, vTh + 8), idx);
        tlr = hn::TwoTablesLookupLanes(d, hn::LoadU(d, vTlr), hn::LoadU(d, vTlr + 8), idx);
    } else if (N == 4) {
        const auto idx = hn::IndicesFromVec(d, hn::And(ji, hn::Set(di, 7)));
        const auto hi = hn::RebindMask(d, hn::TestBit(ji, hn::Set(di, 8)));
        th = hn::IfThenElse(hi, hn::TwoTablesLookupLanes(d, hn::LoadU(d, vTh + 8), hn::LoadU(d, vTh + 12), idx),
                            hn::TwoTablesLookupLanes(d, hn::LoadU(d, vTh), hn::LoadU(d, vTh + 4), idx));
        tlr = hn::IfThenElse(hi, hn::TwoTablesLookupLanes(d, hn::LoadU(d, vTlr + 8), hn::LoadU(d, vTlr + 12), idx),
                             hn::TwoTablesLookupLanes(d, hn::LoadU(d, vTlr), hn::LoadU(d, vTlr + 4), idx));
    } else if (N == 2) {
        const auto idx = hn::IndicesFromVec(d, hn::And(ji, hn::Set(di, 3)));
        const auto b2 = hn::RebindMask(d, hn::TestBit(ji, hn::Set(di, 4)));
        const auto b3 = hn::RebindMask(d, hn::TestBit(ji, hn::Set(di, 8)));
        auto quad = [&](const double *t) {  // t[0..3]
            return hn::TwoTablesLookupLanes(d, hn::LoadU(d, t), hn::LoadU(d, t + 2), idx);
        };
        th = hn::IfThenElse(b3, hn::IfThenElse(b2, quad(vTh + 12), quad(vTh + 8)),
                            hn::IfThenElse(b2, quad(vTh + 4), quad(vTh)));
        tlr = hn::IfThenElse(b3, hn::IfThenElse(b2, quad(vTlr + 12), quad(vTlr + 8)),
                             hn::IfThenElse(b2, quad(vTlr + 4), quad(vTlr)));
    } else {  // e.g. single-lane targets
        th = hn::GatherIndex(d, vTh, ji);
        tlr = hn::GatherIndex(d, vTlr, ji);
    }
}

// vscalefpd(y, n) = y*2^floor(n), for the lanes kept from the vector path: bits 4..15 of tb hold
// 1023 + floor(n) (in [2, 2044]), so 2^floor(n) is normal and the product is exact.
template <class D, class VU>
HWY_INLINE hn::VFromD<D> Scalef(D d, hn::VFromD<D> y, hn::VFromD<D> /*n*/, VU tb) {
    return hn::Mul(y, hn::BitCast(d, hn::ShiftLeft<52>(hn::ShiftRight<4>(tb))));
}

template <class D>
HWY_INLINE hn::VFromD<D> ExpCore(D d, hn::VFromD<D> x) {
    const auto tb = ShifterRZ(d, x);                         // its low 4 bits are the table index
    const auto n = hn::Sub(hn::BitCast(d, tb), hn::Set(d, kShifter));  // x*log2(e) truncated to 1/16
    hn::VFromD<D> th, tlr;
    Lookup16(d, tb, th, tlr);

    auto r = hn::NegMulAdd(n, hn::Set(d, kLn2Hi), x);        // x - N*ln2_hi
    r = hn::NegMulAdd(hn::Set(d, kLn2Lo), n, r);             // - N*ln2_lo

    const auto r2 = hn::Mul(r, r);
    const auto p56 = hn::MulAdd(hn::Set(d, kC6), r, hn::Set(d, kC5));
    const auto p34 = hn::MulAdd(hn::Set(d, kC4), r, hn::Set(d, kC3));
    const auto p12 = hn::MulAdd(hn::Set(d, kC2), r, hn::Set(d, kC1));
    const auto p36 = hn::MulAdd(p56, r2, p34);
    const auto p16 = hn::MulAdd(p36, r2, p12);
    const auto p = hn::MulAdd(p16, r, tlr);
    const auto y = hn::MulAdd(th, p, th);                    // Th*(1+p)
    return Scalef(d, y, n, tb);
}

void ExpSvmlArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    const auto thr = hn::Set(d, kThreshold);
    size_t i = 0;
    for (; i + N <= len; i += N) {
        const auto v = hn::LoadU(d, x + i);
        hn::StoreU(ExpCore(d, v), d, y + i);
        if (HWY_UNLIKELY(!hn::AllTrue(d, hn::Lt(hn::Abs(v), thr)))) FixSpecial(x + i, y + i, N);
    }
    if (i < len) {
        const size_t rem = len - i;
        const auto v = hn::LoadN(d, x + i, rem);
        hn::StoreN(ExpCore(d, v), d, y + i, rem);
        FixSpecial(x + i, y + i, rem);
    }
}

#else  // no fused multiply-add: the algorithm would not be reproduced exactly

void ExpSvmlArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    for (size_t i = 0; i < len; i++) y[i] = std::exp(x[i]);
}

#endif

size_t LanesF64() { return hn::Lanes(hn::ScalableTag<double>()); }

int Variant() {  // 1: SVML algorithm (bit-identical), 0: std::exp fallback
#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64
    return 1;
#else
    return 0;
#endif
}

}  // namespace HWY_NAMESPACE
}  // namespace hwyexpsvml
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace hwyexpsvml {
HWY_EXPORT(ExpSvmlArray);
HWY_EXPORT(Variant);
HWY_EXPORT(LanesF64);
// dispatch inside the namespace: with a single (static) target, HWY_DYNAMIC_DISPATCH(F) becomes N_<TARGET>::F
static void DispatchExp(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(ExpSvmlArray)(x, y, len); }
static int DispatchVariant() { return HWY_DYNAMIC_DISPATCH(Variant)(); }
static size_t DispatchLanes() { return HWY_DYNAMIC_DISPATCH(LanesF64)(); }
}  // namespace hwyexpsvml

extern "C" {
void hwy_exp_svml_array(const double *x, double *y, size_t len) { hwyexpsvml::DispatchExp(x, y, len); }
int hwy_exp_svml_variant(void) { return hwyexpsvml::DispatchVariant(); }
size_t hwy_exp_svml_lanes(void) { return hwyexpsvml::DispatchLanes(); }
}
#endif  // HWY_ONCE
