// Portable (Google Highway) version of Intel SVML's AVX-512 high-accuracy tanh (svml/svml_tanh_d_ha.c):
// the same algorithm, constants and operation order, with run-time dispatch to the best target.
//   - targets with true FMA (HWY_NATIVE_FMA): bit-identical to the AVX-512 code. The vector path uses
//     only round-to-nearest arithmetic, so the only AVX-512-specific operations are emulated exactly:
//       vpermt2pd (19 tables of 16 entries)   see Group below
//       vpmovqd + 32-bit and/sub/max/min/shift on the high words (subinterval index)
//                                             the same clamp on the exponent bits in 32-bit lanes
//       vpcmpgtd + vmovmskps (callout lanes)  64-bit signed compare
//   - targets without FMA (SSE2/SSSE3/SSE4): std::tanh.
// Lanes with |x| >= 1.5*2^1023 (and inf/NaN) use the scalar SVML callout (svml_dtanh_ha_cout_rare).
// Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

extern "C" int svml_dtanh_ha_cout_rare(const double *px, double *pr);

#include "functions/tanh/hwy/hwy_tanh_constants.h"

// Table lookups on targets with 2 or 4 lanes, for each group of 4 tables: 2 = by rows from the interleaved
// table (4 lanes: half-row loads combined into vectors, then interleaved; 2 lanes: row loads, then
// interleaved), 1 = the same with full-row loads and a 4x4 transpose on 4 lanes, 0 = permute/blend tree
// (TwoTablesLookupLanes + IfThenElse). Measurements: README.md.
#ifndef HWYTANH_ROWS
#define HWYTANH_ROWS 2
#endif

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "functions/tanh/hwy/hwy_tanh.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace hwytanh {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64

static inline bool NeedsCallout(double x) {
    uint64_t b;
    std::memcpy(&b, &x, sizeof b);
    return (b & kAbsMask) > kThresholdBits;
}

// scalar fix-up for the lanes Intel's code sends to the callout (rare: kept out of line)
static HWY_NOINLINE __attribute__((cold)) void FixSpecial(const double *x, double *y, size_t n) {
    for (size_t k = 0; k < n; k++) {
        if (NeedsCallout(x[k])) svml_dtanh_ha_cout_rare(&x[k], &y[k]);
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

// Coefficient tables 4g .. 4g+3 at index j (table 19 is zero). idx: j stored to memory (rows lookups).
template <int g, class D, class VI, class V = hn::VFromD<D>>
HWY_INLINE void Group(D d, VI j, const int64_t *idx, V &c0, V &c1, V &c2, V &c3) {
    const size_t N = hn::Lanes(d);
#if HWYTANH_ROWS
    // By rows, per group of 4 tables: 4 lanes: 8 half-row loads + 4 inserts + 4 unpacks (HWYTANH_ROWS 2)
    // or 4 row loads + 8 shuffles (1); 2 lanes: 4 loads + 4 zips
    if (N == 4 && HWY_ARCH_X86) {  // 128-bit-block semantics of InterleaveLower/Upper: AVX2
#if HWYTANH_ROWS == 2
        const hn::Half<D> dh;
        const double *p0 = kRows.v + kRowLen * idx[0] + 4 * g, *p1 = kRows.v + kRowLen * idx[1] + 4 * g;
        const double *p2 = kRows.v + kRowLen * idx[2] + 4 * g, *p3 = kRows.v + kRowLen * idx[3] + 4 * g;
        const auto a = hn::Combine(d, hn::Load(dh, p2), hn::Load(dh, p0));  // r0[0] r0[1] | r2[0] r2[1]
        const auto b = hn::Combine(d, hn::Load(dh, p3), hn::Load(dh, p1));
        const auto e = hn::Combine(d, hn::Load(dh, p2 + 2), hn::Load(dh, p0 + 2));
        const auto f = hn::Combine(d, hn::Load(dh, p3 + 2), hn::Load(dh, p1 + 2));
        c0 = hn::InterleaveLower(d, a, b);
        c1 = hn::InterleaveUpper(d, a, b);
        c2 = hn::InterleaveLower(d, e, f);
        c3 = hn::InterleaveUpper(d, e, f);
#else
        const auto r0 = hn::Load(d, kRows.v + kRowLen * idx[0] + 4 * g);
        const auto r1 = hn::Load(d, kRows.v + kRowLen * idx[1] + 4 * g);
        const auto r2 = hn::Load(d, kRows.v + kRowLen * idx[2] + 4 * g);
        const auto r3 = hn::Load(d, kRows.v + kRowLen * idx[3] + 4 * g);
        const auto t0 = hn::InterleaveLower(d, r0, r1);  // r0[0] r1[0] | r0[2] r1[2]
        const auto t1 = hn::InterleaveUpper(d, r0, r1);  // r0[1] r1[1] | r0[3] r1[3]
        const auto t2 = hn::InterleaveLower(d, r2, r3);
        const auto t3 = hn::InterleaveUpper(d, r2, r3);
        c0 = hn::ConcatLowerLower(d, t2, t0);
        c1 = hn::ConcatLowerLower(d, t3, t1);
        c2 = hn::ConcatUpperUpper(d, t2, t0);
        c3 = hn::ConcatUpperUpper(d, t3, t1);
#endif
        return;
    }
    if (N == 2) {
        const double *p0 = kRows.v + kRowLen * idx[0] + 4 * g;
        const double *p1 = kRows.v + kRowLen * idx[1] + 4 * g;
        const auto a0 = hn::Load(d, p0), b0 = hn::Load(d, p0 + 2);
        const auto a1 = hn::Load(d, p1), b1 = hn::Load(d, p1 + 2);
        c0 = hn::InterleaveLower(d, a0, a1);
        c1 = hn::InterleaveUpper(d, a0, a1);
        c2 = hn::InterleaveLower(d, b0, b1);
        c3 = hn::InterleaveUpper(d, b0, b1);
        return;
    }
#endif
    (void)N;
    (void)idx;
    c0 = Lookup16(d, kTab[4 * g], j);
    c1 = Lookup16(d, kTab[4 * g + 1], j);
    c2 = Lookup16(d, kTab[4 * g + 2], j);
    if constexpr (4 * g + 3 < kTables) {
        c3 = Lookup16(d, kTab[4 * g + 3], j);
    } else {
        c3 = hn::Zero(d);
    }
}

template <class D>
HWY_INLINE hn::VFromD<D> TanhCore(D d, hn::VFromD<D> x) {
    using V = hn::VFromD<D>;
    const hn::RebindToUnsigned<D> du;
    const hn::RebindToSigned<D> di;
    const hn::Repartition<int32_t, D> d32;
    const V ax = hn::And(x, hn::BitCast(d, hn::Set(du, kAbsMask)));
    const V sign = hn::And(x, hn::BitCast(d, hn::Set(du, kSignMask)));

    // j = clamp(u - kIndexBase, 0, 15), u = (bits >> 51) & 0xfff, in 32-bit lanes (the upper halves of the
    // 64-bit lanes are 0 and stay 0: max(0, base) - base)
    const auto u = hn::And(hn::ShiftRight<51>(hn::BitCast(du, x)), hn::Set(du, 0xfff));
    auto w = hn::Max(hn::BitCast(d32, u), hn::Set(d32, static_cast<int32_t>(kIndexBase)));
    w = hn::Min(w, hn::Set(d32, static_cast<int32_t>(kIndexMax)));
    w = hn::Sub(w, hn::Set(d32, static_cast<int32_t>(kIndexBase)));
    const auto j = hn::BitCast(di, w);
    alignas(64) int64_t idx[HWY_MAX_LANES_D(D) < 4 ? 4 : HWY_MAX_LANES_D(D)];
    if (HWYTANH_ROWS && hn::Lanes(d) <= 4) hn::Store(j, di, idx);

    V b, p0h, p0l, p1;
    Group<0>(d, j, idx, b, p0h, p0l, p1);
    const V y = hn::Sub(ax, b);  // |x| - B
    V c0, c1, c2, c3;
    Group<4>(d, j, idx, c0, c1, c2, c3);  // P14 P15 P16 -
    V p = hn::MulAdd(y, c2, c1);
    p = hn::MulAdd(y, p, c0);
    Group<3>(d, j, idx, c0, c1, c2, c3);  // P10 .. P13
    p = hn::MulAdd(y, p, c3);
    p = hn::MulAdd(y, p, c2);
    p = hn::MulAdd(y, p, c1);
    p = hn::MulAdd(y, p, c0);
    Group<2>(d, j, idx, c0, c1, c2, c3);  // P6 .. P9
    p = hn::MulAdd(y, p, c3);
    p = hn::MulAdd(y, p, c2);
    p = hn::MulAdd(y, p, c1);
    p = hn::MulAdd(y, p, c0);
    Group<1>(d, j, idx, c0, c1, c2, c3);  // P2 .. P5
    p = hn::MulAdd(y, p, c3);
    p = hn::MulAdd(y, p, c2);
    p = hn::MulAdd(y, p, c1);
    p = hn::MulAdd(y, p, c0);
    V q = hn::Mul(p, y);
    q = hn::MulAdd(y, q, p0l);
    p = hn::MulAdd(y, p1, q);
    const V res = hn::Add(p, p0h);
    return hn::Or(res, sign);
}

template <class D>
HWY_INLINE bool AnyCallout(D d, hn::VFromD<D> x) {
    const hn::RebindToSigned<D> di;
    const auto ab = hn::And(hn::BitCast(di, x), hn::Set(di, static_cast<int64_t>(kAbsMask)));
    return !hn::AllFalse(di, hn::Gt(ab, hn::Set(di, static_cast<int64_t>(kThresholdBits))));
}

template <class D>
HWY_INLINE void TanhArrayD(D d, const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    const size_t N = hn::Lanes(d);
    size_t i = 0;
    for (; i + N <= len; i += N) {
        const auto v = hn::LoadU(d, x + i);
        hn::StoreU(TanhCore(d, v), d, y + i);
        if (HWY_UNLIKELY(AnyCallout(d, v))) FixSpecial(x + i, y + i, N);
    }
    if (i < len) {
        const size_t rem = len - i;
        const auto v = hn::LoadN(d, x + i, rem);
        hn::StoreN(TanhCore(d, v), d, y + i, rem);
        FixSpecial(x + i, y + i, rem);
    }
}

void TanhArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    TanhArrayD(hn::ScalableTag<double>(), x, y, len);
}

// Tests only: vectors of at most 2 or 4 lanes, to exercise the 2- and 4-lane lookups on wider targets.
void TanhArrayCapped(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len, size_t lanes) {
    if (lanes == 2) TanhArrayD(hn::CappedTag<double, 2>(), x, y, len);
    else TanhArrayD(hn::CappedTag<double, 4>(), x, y, len);
}

int Variant() { return 1; }  // 1: SVML algorithm, 0: std::tanh fallback

#else  // no fused multiply-add: the SVML algorithm cannot be reproduced

void TanhArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    for (size_t i = 0; i < len; i++) y[i] = std::tanh(x[i]);
}

void TanhArrayCapped(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len, size_t) {
    TanhArray(x, y, len);
}

int Variant() { return 0; }

#endif

size_t LanesF64() { return hn::Lanes(hn::ScalableTag<double>()); }

}  // namespace HWY_NAMESPACE
}  // namespace hwytanh
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace hwytanh {
HWY_EXPORT(TanhArray);
HWY_EXPORT(TanhArrayCapped);
HWY_EXPORT(Variant);
HWY_EXPORT(LanesF64);
// dispatch inside the namespace: with a single (static) target, HWY_DYNAMIC_DISPATCH(F) becomes N_<TARGET>::F
static void DispatchTanh(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(TanhArray)(x, y, len); }
static void DispatchCapped(const double *x, double *y, size_t len, size_t lanes) {
    HWY_DYNAMIC_DISPATCH(TanhArrayCapped)(x, y, len, lanes);
}
static int DispatchVariant() { return HWY_DYNAMIC_DISPATCH(Variant)(); }
static size_t DispatchLanes() { return HWY_DYNAMIC_DISPATCH(LanesF64)(); }
}  // namespace hwytanh

extern "C" {
void hwy_tanh_array(const double *x, double *y, size_t len) { hwytanh::DispatchTanh(x, y, len); }
void hwy_tanh_array_capped(const double *x, double *y, size_t len, size_t lanes) {
    hwytanh::DispatchCapped(x, y, len, lanes);
}
int hwy_tanh_variant(void) { return hwytanh::DispatchVariant(); }
size_t hwy_tanh_lanes(void) { return hwytanh::DispatchLanes(); }
}
#endif  // HWY_ONCE
