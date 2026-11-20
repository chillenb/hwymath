// Per-target helpers shared by the Highway versions of SVML's log, log2, log10 and log1p
// (functions/{log,log2,log10,log1p}/hwy/hwy_*.cc): exact emulations of the AVX-512-only operations of the
// SVML algorithm. Re-included once per target by foreach_target.h (include guard toggled per target).
// Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.

// clang-format off
#if defined(HWYMATH_LOG_HWY_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef HWYMATH_LOG_HWY_INL_H_
#undef HWYMATH_LOG_HWY_INL_H_
#else
#define HWYMATH_LOG_HWY_INL_H_
#endif
// clang-format on

#include <cstddef>
#include <cstdint>

#include "functions/log/hwy/log_tables.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace hwymath_log {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

// v[i] = tab[idx[i]] for a 16-entry table (vpermt2pd on AVX-512), idx < 16. pad: the table repeated to
// 64 entries, for vectors of 16 or more lanes. The lane count may be a run-time value (SVE, RVV).
template <class D, class VU>
HWY_INLINE hn::VFromD<D> Lookup16(D d, const double *HWY_RESTRICT tab, const double *HWY_RESTRICT pad, VU idx) {
    const hn::RebindToSigned<D> di;
    const hn::RebindToUnsigned<D> du;
    const size_t N = hn::Lanes(d);
    if (N >= 16 && N <= 64) {  // table (repeated) in one vector
        return hn::TableLookupLanes(hn::LoadU(d, pad), hn::IndicesFromVec(d, hn::BitCast(di, idx)));
    }
    if (N == 8) {  // two vectors: vpermt2pd on AVX-512
        return hn::TwoTablesLookupLanes(d, hn::LoadU(d, tab), hn::LoadU(d, tab + 8),
                                        hn::IndicesFromVec(d, hn::BitCast(di, idx)));
    }
    if (N == 4 || N == 2) {  // 16/(2N) two-vector lookups, then selects on the high index bits
        const auto lane = hn::IndicesFromVec(d, hn::BitCast(di, hn::And(idx, hn::Set(du, 2 * N - 1))));
        const auto t = [&](size_t c) {  // entries 2N*c .. 2N*c + 2N-1
            return hn::TwoTablesLookupLanes(d, hn::LoadU(d, tab + 2 * N * c), hn::LoadU(d, tab + 2 * N * c + N), lane);
        };
        const auto bit = [&](uint64_t b) { return hn::RebindMask(d, hn::Ne(hn::And(idx, hn::Set(du, b)), hn::Zero(du))); };
        const auto hi8 = bit(8);
        if (N == 4) return hn::IfThenElse(hi8, t(1), t(0));
        const auto hi4 = bit(4);
        return hn::IfThenElse(hi8, hn::IfThenElse(hi4, t(3), t(2)), hn::IfThenElse(hi4, t(1), t(0)));
    }
    return hn::GatherIndex(d, tab, hn::BitCast(di, idx));  // single-lane targets
}

// v[i] = tab[idx[i]] for a 32-entry table of uint32 (vpermt2d on AVX-512), idx < 32. pad: the table
// repeated to 64 entries. The lane count may be a run-time value (SVE, RVV).
template <class D32, class VU32>
HWY_INLINE VU32 Lookup32U32(D32 d32, const uint32_t *HWY_RESTRICT tab, const uint32_t *HWY_RESTRICT pad, VU32 idx) {
    const hn::RebindToSigned<D32> di;
    const size_t N = hn::Lanes(d32);
    if (N >= 32 && N <= 64) {  // table (repeated) in one vector
        return hn::TableLookupLanes(hn::LoadU(d32, pad), hn::IndicesFromVec(d32, hn::BitCast(di, idx)));
    }
    if (N == 16) {  // two vectors
        return hn::TwoTablesLookupLanes(d32, hn::LoadU(d32, tab), hn::LoadU(d32, tab + 16),
                                        hn::IndicesFromVec(d32, hn::BitCast(di, idx)));
    }
    if (N == 8 || N == 4) {  // 32/(2N) two-vector lookups, then selects on the high index bits
        const auto lane = hn::IndicesFromVec(d32, hn::BitCast(di, hn::And(idx, hn::Set(d32, uint32_t(2 * N - 1)))));
        const auto t = [&](size_t c) {  // entries 2N*c .. 2N*c + 2N-1
            return hn::TwoTablesLookupLanes(d32, hn::LoadU(d32, tab + 2 * N * c), hn::LoadU(d32, tab + 2 * N * c + N), lane);
        };
        const auto bit = [&](uint32_t b) { return hn::Ne(hn::And(idx, hn::Set(d32, b)), hn::Zero(d32)); };
        const auto hi16 = bit(16);
        if (N == 8) return hn::IfThenElse(hi16, t(1), t(0));
        const auto hi8 = bit(8);
        return hn::IfThenElse(hi16, hn::IfThenElse(hi8, t(3), t(2)), hn::IfThenElse(hi8, t(1), t(0)));
    }
    return hn::GatherIndex(d32, tab, hn::BitCast(di, idx));  // other lane counts (not used by the f64 targets)
}

// rcp = vrndscalepd(vrcp14pd(m), RN to 1/32) = k/32 for m in [1, 2), computed exactly from the thresholds
// in rcp14_thresholds.h, and the table index j = top 4 mantissa bits of rcp = k mod 16. Integer operations
// on the high 32 bits of m (see rcp14_thresholds.h: kRcpBins). For any other bit pattern of m (special lanes)
// the results are finite garbage with an index < 16.
template <class D, class V = hn::VFromD<D>, class VU = hn::VFromD<hn::RebindToUnsigned<D>>>
HWY_INLINE void RoundedRcp14(D d, V m, V &rcp, VU &idx) {
    const hn::RebindToUnsigned<D> du;
    const hn::Repartition<uint32_t, D> d32;
    const auto mb = hn::BitCast(du, m);
    // the bin b = top 5 mantissa bits, in the high 32-bit half of each 64-bit lane (low half: 0)
    const auto b = hn::And(hn::ShiftRight<15>(mb), hn::Set(du, uint64_t{31} << 32));
    const auto e = hn::BitCast(du, Lookup32U32(d32, kRcpBins.v, Tables::BinsPad(), hn::BitCast(d32, b)));
    // high halves: (hi32(m) | 15) - entry < 0 <=> m below the threshold of the bin (32-bit subtraction)
    const auto diff = hn::BitCast(du, hn::Sub(hn::BitCast(d32, hn::Or(mb, hn::Set(du, uint64_t{15} << 32))), hn::BitCast(d32, e)));
    const auto jk = hn::Add(hn::And(hn::ShiftRight<32>(e), hn::Set(du, 15)), hn::ShiftRight<63>(diff));  // k - 16, 0..16
    rcp = hn::BitCast(d, hn::Add(hn::Set(du, 0x3fe0000000000000ULL), hn::ShiftLeft<48>(jk)));  // k/32 (k = 32: 1.0)
    idx = hn::And(jk, hn::Set(du, 15));
}

// x = 2^e * m with m in [1, 2) for finite x > 0 (vgetexppd/vgetmantpd), including denormals.
template <class D, class V = hn::VFromD<D>>
HWY_INLINE void ExpMant(D d, V x, V &e, V &m) {
    const hn::RebindToUnsigned<D> du;
    const auto denorm = hn::Lt(x, hn::Set(d, 0x1p-1022));
    const auto bits = hn::BitCast(du, hn::IfThenElse(denorm, hn::Mul(x, hn::Set(d, 0x1p52)), x));
    m = hn::BitCast(d, hn::Or(hn::And(bits, hn::Set(du, 0x000fffffffffffffULL)), hn::Set(du, 0x3ff0000000000000ULL)));
    // 2^52 + biased exponent, exactly, from the exponent bits
    const auto eb = hn::BitCast(d, hn::Or(hn::ShiftRight<52>(bits), hn::Set(du, 0x4330000000000000ULL)));
    e = hn::Sub(eb, hn::IfThenElse(denorm, hn::Set(d, 0x1p52 + 1075.0), hn::Set(d, 0x1p52 + 1023.0)));
}

// Same for normal x > 0 only (log1p: 1 + x >= 2^-53).
template <class D, class V = hn::VFromD<D>>
HWY_INLINE void ExpMantNormal(D d, V x, V &e, V &m) {
    const hn::RebindToUnsigned<D> du;
    const auto bits = hn::BitCast(du, x);
    m = hn::BitCast(d, hn::Or(hn::And(bits, hn::Set(du, 0x000fffffffffffffULL)), hn::Set(du, 0x3ff0000000000000ULL)));
    const auto eb = hn::BitCast(d, hn::Or(hn::ShiftRight<52>(bits), hn::Set(du, 0x4330000000000000ULL)));
    e = hn::Sub(eb, hn::Set(d, 0x1p52 + 1023.0));
}

}  // namespace HWY_NAMESPACE
}  // namespace hwymath_log
HWY_AFTER_NAMESPACE();

#endif  // HWYMATH_LOG_HWY_INL_H_
