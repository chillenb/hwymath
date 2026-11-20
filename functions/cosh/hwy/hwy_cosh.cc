// Portable (Google Highway) version of Intel SVML's AVX-512 high-accuracy cosh (svml/svml_cosh_d_ha.c):
// the same algorithm, constants and operation order, with run-time dispatch to the best target.
//   - targets with true FMA (HWY_NATIVE_FMA): bit-identical to the AVX-512 code. The vector path uses
//     only round-to-nearest arithmetic, so the only AVX-512-specific operations are emulated exactly:
//       vpermt2pd (16-entry tables)   see Lookup below
//       vpmovqd + vpcmpgtd + vmovmskps (high words of |x| > 0x40861d99)  64-bit signed compare
//   - targets without FMA (SSE2/SSSE3/SSE4): std::cosh.
// Lanes with |x| >= 707.70 (and inf/NaN) use the scalar SVML callout (svml_dcosh_ha_cout_rare).
// Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

extern "C" int svml_dcosh_ha_cout_rare(const double *px, double *pr);

#include "functions/cosh/hwy/hwy_cosh_constants.h"

// Table lookups on targets with 2 or 4 lanes: 1 = by rows (one load per lane from the interleaved table,
// then a transpose), 0 = permute/blend tree (TwoTablesLookupLanes + IfThenElse). See README.md.
#ifndef HWYCOSH_ROWS
#define HWYCOSH_ROWS 1
#endif

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "functions/cosh/hwy/hwy_cosh.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace hwycosh {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64

static inline bool NeedsCallout(double x) {
    uint64_t b;
    std::memcpy(&b, &x, sizeof b);
    return (b & ~kSignMask) > kThresholdBits;
}

// scalar fix-up for the lanes Intel's code sends to the callout (rare: kept out of line)
static HWY_NOINLINE __attribute__((cold)) void FixSpecial(const double *x, double *y, size_t n) {
    for (size_t k = 0; k < n; k++) {
        if (NeedsCallout(x[k])) svml_dcosh_ha_cout_rare(&x[k], &y[k]);
    }
}

// 16-entry lookup tab[j] (j < 16, as 64-bit integers) with permutes: one vpermt2pd on AVX-512, a tree of
// two-table lookups and selects on bits 2/3 of j for 4 and 2 lanes.
template <class D, class VI>
HWY_INLINE hn::VFromD<D> Lookup16(D d, const double *tab, VI j) {
    const hn::RebindToSigned<D> di;
    const size_t N = hn::Lanes(d);
    if (N >= 16) {
        return hn::TableLookupLanes(hn::LoadN(d, tab, 16), hn::IndicesFromVec(d, j));
    }
    if (N == 8) {
        return hn::TwoTablesLookupLanes(d, hn::LoadU(d, tab), hn::LoadU(d, tab + 8), hn::IndicesFromVec(d, j));
    }
    const auto bit3 = hn::RebindMask(d, hn::TestBit(j, hn::Set(di, 8)));
    if (N == 4) {
        const auto idx = hn::IndicesFromVec(d, hn::And(j, hn::Set(di, 7)));
        const auto lo = hn::TwoTablesLookupLanes(d, hn::LoadU(d, tab), hn::LoadU(d, tab + 4), idx);
        const auto hi = hn::TwoTablesLookupLanes(d, hn::LoadU(d, tab + 8), hn::LoadU(d, tab + 12), idx);
        return hn::IfThenElse(bit3, hi, lo);
    }
    if (N == 2) {
        const auto bit2 = hn::RebindMask(d, hn::TestBit(j, hn::Set(di, 4)));
        const auto idx = hn::IndicesFromVec(d, hn::And(j, hn::Set(di, 3)));
        const auto t0 = hn::TwoTablesLookupLanes(d, hn::LoadU(d, tab), hn::LoadU(d, tab + 2), idx);
        const auto t1 = hn::TwoTablesLookupLanes(d, hn::LoadU(d, tab + 4), hn::LoadU(d, tab + 6), idx);
        const auto t2 = hn::TwoTablesLookupLanes(d, hn::LoadU(d, tab + 8), hn::LoadU(d, tab + 10), idx);
        const auto t3 = hn::TwoTablesLookupLanes(d, hn::LoadU(d, tab + 12), hn::LoadU(d, tab + 14), idx);
        return hn::IfThenElse(bit3, hn::IfThenElse(bit2, t3, t2), hn::IfThenElse(bit2, t1, t0));
    }
    return hn::GatherIndex(d, tab, j);  // not reached on FMA targets (N = 1)
}

// The three lookups TpH, TpL, TmH at index j.
template <class D, class VI, class V = hn::VFromD<D>>
HWY_INLINE void Lookup(D d, VI j, V &tph, V &tpl, V &tmh) {
    const size_t N = hn::Lanes(d);
#if HWYCOSH_ROWS
    // By rows: the {TpH, TpL, TmH, 0} row of each lane is loaded from the interleaved table and the
    // rows are transposed (4 lanes: 4 loads + 6 shuffles; 2 lanes: 4 loads + 3 zips).
    if (N == 4 && HWY_ARCH_X86) {  // 128-bit-block semantics of InterleaveLower/Upper: AVX2
        alignas(32) int64_t idx[4];
        hn::Store(j, hn::RebindToSigned<D>(), idx);
        const auto r0 = hn::Load(d, kRows.v + 4 * idx[0]);
        const auto r1 = hn::Load(d, kRows.v + 4 * idx[1]);
        const auto r2 = hn::Load(d, kRows.v + 4 * idx[2]);
        const auto r3 = hn::Load(d, kRows.v + 4 * idx[3]);
        const auto t0 = hn::InterleaveLower(d, r0, r1);  // r0[0] r1[0] | r0[2] r1[2]
        const auto t1 = hn::InterleaveUpper(d, r0, r1);  // r0[1] r1[1] | r0[3] r1[3]
        const auto t2 = hn::InterleaveLower(d, r2, r3);
        const auto t3 = hn::InterleaveUpper(d, r2, r3);
        tph = hn::ConcatLowerLower(d, t2, t0);
        tpl = hn::ConcatLowerLower(d, t3, t1);
        tmh = hn::ConcatUpperUpper(d, t2, t0);
        return;
    }
    if (N == 2) {
        alignas(16) int64_t idx[2];
        hn::Store(j, hn::RebindToSigned<D>(), idx);
        const double *p0 = kRows.v + 4 * idx[0];
        const double *p1 = kRows.v + 4 * idx[1];
        const auto a0 = hn::Load(d, p0), b0 = hn::Load(d, p0 + 2);
        const auto a1 = hn::Load(d, p1), b1 = hn::Load(d, p1 + 2);
        tph = hn::InterleaveLower(d, a0, a1);
        tpl = hn::InterleaveUpper(d, a0, a1);
        tmh = hn::InterleaveLower(d, b0, b1);
        return;
    }
#endif
    (void)N;
    tph = Lookup16(d, kTpH, j);
    tpl = Lookup16(d, kTpL, j);
    tmh = Lookup16(d, kTmH, j);
}

template <class D>
HWY_INLINE hn::VFromD<D> CoshCore(D d, hn::VFromD<D> x) {
    using V = hn::VFromD<D>;
    const hn::RebindToUnsigned<D> du;
    const hn::RebindToSigned<D> di;
    const V S = hn::Set(d, kShifter);
    const V ax = hn::AndNot(hn::BitCast(d, hn::Set(du, kSignMask)), x);
    const V m = hn::MulAdd(ax, hn::Set(d, kLog2e), S);  // a*log2(e) + 1.5*2^48
    const auto mb = hn::BitCast(du, m);
    V tph, tpl, tmh;
    Lookup(d, hn::BitCast(di, hn::And(mb, hn::Set(du, 15))), tph, tpl, tmh);

    const V n = hn::Sub(m, S);
    V r = hn::NegMulAdd(n, hn::Set(d, kLn2Hi), ax);  // a - N*ln2_hi
    r = hn::NegMulAdd(n, hn::Set(d, kLn2Lo), r);     // - N*ln2_lo
    const V r2 = hn::Mul(r, r);
    const auto e = hn::And(hn::ShiftLeft<48>(mb), hn::Set(du, kExpMask));  // M << 52
    const V tpl2 = hn::BitCast(d, hn::Add(e, hn::BitCast(du, tpl)));
    const V tp = hn::BitCast(d, hn::Add(e, hn::BitCast(du, tph)));
    const V tm = hn::BitCast(d, hn::Sub(hn::BitCast(du, tmh), e));

    const V sum = hn::Add(tp, tm);
    const V dif = hn::Sub(tp, tm);
    V po = hn::MulAdd(r2, hn::Set(d, kS7), hn::Set(d, kS5));
    po = hn::MulAdd(r2, po, hn::Set(d, kS3));
    V pe = hn::MulAdd(r2, hn::Set(d, kC6), hn::Set(d, kC4));
    po = hn::Mul(po, r2);
    pe = hn::MulAdd(r2, pe, hn::Set(d, kC2));
    po = hn::MulAdd(r, po, r);
    pe = hn::Mul(pe, r2);
    pe = hn::MulAdd(sum, pe, tpl2);
    V res = hn::MulAdd(dif, po, pe);
    res = hn::Add(res, tm);
    return hn::Add(res, tp);
}

template <class D>
HWY_INLINE bool AnyCallout(D d, hn::VFromD<D> x) {
    const hn::RebindToSigned<D> di;
    const auto ab = hn::And(hn::BitCast(di, x), hn::Set(di, static_cast<int64_t>(~kSignMask)));
    return !hn::AllFalse(di, hn::Gt(ab, hn::Set(di, static_cast<int64_t>(kThresholdBits))));
}

template <class D>
HWY_INLINE void CoshArrayD(D d, const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    const size_t N = hn::Lanes(d);
    size_t i = 0;
    for (; i + N <= len; i += N) {
        const auto v = hn::LoadU(d, x + i);
        hn::StoreU(CoshCore(d, v), d, y + i);
        if (HWY_UNLIKELY(AnyCallout(d, v))) FixSpecial(x + i, y + i, N);
    }
    if (i < len) {
        const size_t rem = len - i;
        const auto v = hn::LoadN(d, x + i, rem);
        hn::StoreN(CoshCore(d, v), d, y + i, rem);
        FixSpecial(x + i, y + i, rem);
    }
}

void CoshArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    CoshArrayD(hn::ScalableTag<double>(), x, y, len);
}

// Tests only: vectors of at most 2 or 4 lanes, to exercise the 2- and 4-lane lookups on wider targets.
void CoshArrayCapped(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len, size_t lanes) {
    if (lanes == 2) CoshArrayD(hn::CappedTag<double, 2>(), x, y, len);
    else CoshArrayD(hn::CappedTag<double, 4>(), x, y, len);
}

int Variant() { return 1; }  // 1: SVML algorithm, 0: std::cosh fallback

#else  // no fused multiply-add: the SVML algorithm cannot be reproduced

void CoshArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    for (size_t i = 0; i < len; i++) y[i] = std::cosh(x[i]);
}

void CoshArrayCapped(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len, size_t) {
    CoshArray(x, y, len);
}

int Variant() { return 0; }

#endif

size_t LanesF64() { return hn::Lanes(hn::ScalableTag<double>()); }

}  // namespace HWY_NAMESPACE
}  // namespace hwycosh
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace hwycosh {
HWY_EXPORT(CoshArray);
HWY_EXPORT(CoshArrayCapped);
HWY_EXPORT(Variant);
HWY_EXPORT(LanesF64);
// dispatch inside the namespace: with a single (static) target, HWY_DYNAMIC_DISPATCH(F) becomes N_<TARGET>::F
static void DispatchCosh(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(CoshArray)(x, y, len); }
static void DispatchCapped(const double *x, double *y, size_t len, size_t lanes) {
    HWY_DYNAMIC_DISPATCH(CoshArrayCapped)(x, y, len, lanes);
}
static int DispatchVariant() { return HWY_DYNAMIC_DISPATCH(Variant)(); }
static size_t DispatchLanes() { return HWY_DYNAMIC_DISPATCH(LanesF64)(); }
}  // namespace hwycosh

extern "C" {
void hwy_cosh_array(const double *x, double *y, size_t len) { hwycosh::DispatchCosh(x, y, len); }
void hwy_cosh_array_capped(const double *x, double *y, size_t len, size_t lanes) {
    hwycosh::DispatchCapped(x, y, len, lanes);
}
int hwy_cosh_variant(void) { return hwycosh::DispatchVariant(); }
size_t hwy_cosh_lanes(void) { return hwycosh::DispatchLanes(); }
}
#endif  // HWY_ONCE
