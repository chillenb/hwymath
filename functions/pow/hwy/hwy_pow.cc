// Portable (Google Highway) version of Intel SVML's AVX-512 high-accuracy pow (svml/svml_pow_d_ha.c), with
// run-time dispatch to the best target. Same algorithm, constants and operation order. AVX-512 targets
// (HWYMATH_NATIVE_AVX512) use the assembly's instructions; elsewhere the AVX-512-only operations are emulated
// exactly, so on every target with true FMA (HWY_NATIVE_FMA) the results are bit-identical to Intel's
// __svml_pow8_ha:
//   vgetmantpd/vgetexppd   exponent and mantissa bits (denormals scaled by 2^52 first)
//   vrcp14pd + vrndscalepd the rounded reciprocal depends only on the top 16 fraction bits of the
//                          mantissa (Intel's reference implementation, verified exhaustively under SDE);
//                          round(32/q) on the truncated mantissa q, corrected near the 32 switch points
//                          (AVX-512 targets use the instructions themselves)
//                          with a table (rcp14_switch.inc)
//   {rz-sae} mul/fma/add   round-to-nearest result, stepped toward zero by one ulp when the exact error
//                          (TwoProduct / TwoSum / ErrFma) has the opposite sign
//   {rd-sae} add           shifter + floor(16 S)/16 (exact)
//   vreducepd (floor, RD)  S - floor(16 S)/16, rounded down when inexact (only for -1/32 < S < 0)
//   vpermt2pd + blend      TableLookupLanes / TwoTablesLookupLanes (+ selects for vectors < 16 lanes)
//   vfpclasspd, vcmppd     compares
// Lanes that the SVML vector path sends to its scalar callout (x <= 0, x or y inf/NaN, |y log2 x| >
// 1021.5) use the C translation of that callout (svml_dpow_ha_cout_rare). Targets without FMA use std::pow.
// Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.

#include <cmath>
#include <cstddef>
#include <cstdint>

extern "C" int svml_dpow_ha_cout_rare(const double *px, const double *py, double *pr);

#include "functions/pow/hwy/hwy_pow_constants.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "functions/pow/hwy/hwy_pow.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"
#include "hwymath/hwy/avx512_native-inl.h"

HWY_BEFORE_NAMESPACE();
namespace hwypow {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64

template <class D>
using V = hn::VFromD<D>;

// scalar callout for the lanes flagged in bits (rare: kept out of line)
static HWY_NOINLINE __attribute__((cold)) void FixSpecial(const double *x, const double *y, double *z, size_t n,
                                                          const uint8_t *bits) {
    for (size_t k = 0; k < n; k++) {
        if (bits[k / 8] >> (k % 8) & 1) svml_dpow_ha_cout_rare(&x[k], &y[k], &z[k]);
    }
}

// entries [lo, lo + kSize) of a table, kSize >= 2N, idx < T: TwoTablesLookupLanes on chunks of 2N entries
// (vpermt2pd on AVX-512, like the assembly) combined by a binary tree of selects on the higher index bits
// (blend on the sign bit). Recursion instead of arrays: SVE vectors are sizeless and cannot be array
// elements, also on the fixed-size SVE_256/SVE2_128 targets.
template <size_t kSize, class D, class VU, class VI>
HWY_INLINE V<D> LookupTree(D d, const double *tab, VI ind, VU idx) {
    const size_t N = hn::Lanes(d);
    if (kSize <= 2 * N) return hn::TwoTablesLookupLanes(d, hn::LoadU(d, tab), hn::LoadU(d, tab + N), ind);
    if constexpr (kSize > 2) {
        const hn::RebindToSigned<D> di;
        constexpr int kBit = hwy::FloorLog2(kSize / 2);  // index bit selecting the upper half
        const auto lo = LookupTree<kSize / 2>(d, tab, ind, idx);
        const auto hi = LookupTree<kSize / 2>(d, tab + kSize / 2, ind, idx);
        return hn::IfNegativeThenElse(hn::BitCast(d, hn::ShiftLeft<63 - kBit>(hn::BitCast(di, idx))), hi, lo);
    }
    return hn::Zero(d);  // not reached
}

// lookup in a T-entry table (T = 16 or 32), idx < T: one TableLookupLanes if a vector holds the table,
// otherwise LookupTree
template <size_t T, class D, class VU>
HWY_INLINE V<D> Lookup(D d, const double *tab, const double *pad, VU idx) {
    const hn::RebindToSigned<D> di;
    const hn::RebindToUnsigned<D> du;
    const size_t N = hn::Lanes(d);
    if (N >= T && N <= 64) return hn::TableLookupLanes(hn::LoadU(d, pad), hn::IndicesFromVec(d, hn::BitCast(di, idx)));
    if (N == 1) return hn::Set(d, tab[hn::GetLane(idx)]);
    // TwoTablesLookupLanes: indices must be < 2N (NEON's tbl returns 0 for larger ones)
    const auto ind = hn::IndicesFromVec(d, hn::BitCast(di, hn::And(idx, hn::Set(du, 2 * N - 1))));
    return LookupTree<T>(d, tab, ind, idx);
}

// round toward zero, given the round-to-nearest result r and (the sign of) the exact error of r: step r
// one ulp toward zero when the error has the opposite sign (callers ensure r != 0 when err != 0; for an
// RN product or sum that holds: a zero result is exact or has an error of its own sign)
template <class D>
HWY_INLINE V<D> TowardZero(D d, V<D> r, V<D> err) {
    const hn::RebindToUnsigned<D> du;
    const hn::RebindToSigned<D> di;
    const auto opposite = hn::RebindMask(d, hn::Lt(hn::BitCast(di, hn::Xor(hn::BitCast(du, r), hn::BitCast(du, err))),
                                                    hn::Zero(di)));
    const auto away = hn::And(opposite, hn::Ne(err, hn::Zero(d)));
    return hn::BitCast(d, hn::Add(hn::BitCast(du, r), hn::BitCast(du, hn::VecFromMask(d, away))));  // bits - 1
}

// error of the rounded sum s = a + b (TwoSum, exact)
template <class D>
HWY_INLINE V<D> SumErr(V<D> a, V<D> b, V<D> s) {
    const auto bb = hn::Sub(s, a);
    return hn::Add(hn::Sub(a, hn::Sub(s, bb)), hn::Sub(b, bb));
}

// RZ(a*b + c): sign of the exact error of the RN fma from ErrFma (Boldo & Muller 2011)
template <class D>
HWY_INLINE V<D> FmaRZ(D d, V<D> a, V<D> b, V<D> c) {
#if HWYMATH_NATIVE_AVX512
    (void)d;
    return hwymath::HWY_NAMESPACE::MulAddR<hwymath::HWY_NAMESPACE::kRZ>(a, b, c);  // vfmadd {rz-sae}
#else
    const auto r1 = hn::MulAdd(a, b, c);
    const auto u1 = hn::Mul(a, b);
    const auto u2 = hn::MulSub(a, b, u1);
    const auto a1 = hn::Add(c, u2);
    const auto zz = SumErr<D>(c, u2, a1);
    const auto b1 = hn::Add(u1, a1);
    const auto b2 = SumErr<D>(u1, a1, b1);
    const auto g = hn::Add(hn::Sub(b1, r1), b2);
    // ErrFma's Fast2Sum(g, zz) = (r2, r3) gives a*b + c = r1 + r2 + r3 exactly with |r3| <= ulp(r2)/2, and
    // r2 = 0 only if g + zz = 0 (then r3 = 0): the sign of the error is the sign of r2 = RN(g + zz)
    // (r1 = 0 with a nonzero error only happens if the error underflows; leave such results alone)
    return TowardZero(d, r1, hn::IfThenZeroElse(hn::Eq(r1, hn::Zero(d)), hn::Add(g, zz)));
#endif
}

// vrndscalepd(vrcp14pd(m), 1/32, RN) for m in [1/2, 1): depends only on q = m truncated to 16 fraction bits
// and equals round(32/q)/32 except within a few prefixes of the 32 switch points (32/q within 2^-7 of a
// half-integer; the switch points are at most 3.5 prefixes = 0.0034 in 32/q from the exact midpoints)
template <class D>
HWY_INLINE V<D> RcpRounded(D d, V<D> mant) {
#if HWYMATH_NATIVE_AVX512
    // the instructions themselves: vrcp14pd, vrndscalepd 0x58 (to 1/32, RN, {sae})
    (void)d;
    return hwymath::HWY_NAMESPACE::RoundScale<0x58>(hn::ApproximateReciprocal(mant));
#else
    const hn::RebindToUnsigned<D> du;
    const auto q = hn::BitCast(d, hn::And(hn::BitCast(du, mant), hn::Set(du, ~((uint64_t{1} << 36) - 1))));
    const auto t = hn::Div(hn::Set(d, 32.0), q);  // floor() and round() of it are exact
    auto k = hn::Round(t);
    const auto ft = hn::Floor(t);
    const auto near = hn::Lt(hn::Abs(hn::Sub(hn::Sub(t, ft), hn::Set(d, 0.5))), hn::Set(d, 0x1p-7));
    if (HWY_UNLIKELY(!hn::AllFalse(d, near))) {
        const auto kmid = hn::Min(ft, hn::Set(d, 63.0));  // switch between kmid and kmid + 1
        const auto ki = hn::And(hn::BitCast(du, hn::Add(kmid, hn::Set(d, 0x1p52 - 32))), hn::Set(du, 31));
        const auto thr = Lookup<32>(d, kRcpSwitch.v, kRcpSwitchPad.v, ki);
        const auto kc = hn::Add(kmid, hn::IfThenElseZero(hn::Lt(q, thr), hn::Set(d, 1.0)));
        k = hn::IfThenElse(near, kc, k);
    }
    return hn::Mul(k, hn::Set(d, 1.0 / 32));
#endif
}

// RZ(a*b): the error of the RN product is exact (fma)
template <class D>
HWY_INLINE V<D> MulRZ(D d, V<D> a, V<D> b) {
#if HWYMATH_NATIVE_AVX512
    (void)d;
    return hwymath::HWY_NAMESPACE::MulR<hwymath::HWY_NAMESPACE::kRZ>(a, b);  // vmulpd {rz-sae}
#else
    const auto p = hn::Mul(a, b);
    return TowardZero(d, p, hn::MulSub(a, b, p));
#endif
}

// RZ(a + b) (TwoSum)
template <class D>
HWY_INLINE V<D> AddRZ(D d, V<D> a, V<D> b) {
#if HWYMATH_NATIVE_AVX512
    (void)d;
    return hwymath::HWY_NAMESPACE::AddR<hwymath::HWY_NAMESPACE::kRZ>(a, b);  // vaddpd {rz-sae}
#else
    const auto s = hn::Add(a, b);
    return TowardZero(d, s, SumErr<D>(a, b, s));
#endif
}

// vreducepd(s, imm 0x41) = s - floor(16 s)/16 with the subtraction rounded down, and sh = vaddpd {rd-sae}
// s + kShifter (= kShifter + floor(16 s)/16: the table index and exponent in its low bits)
template <class D>
HWY_INLINE V<D> Reduce16(D d, V<D> s, V<D> &sh) {
#if HWYMATH_NATIVE_AVX512
    namespace nat = hwymath::HWY_NAMESPACE;
    sh = nat::AddR<nat::kRD>(s, hn::Set(d, kShifter));
    return nat::Reduce<0x41>(s);
#else
    const hn::RebindToUnsigned<D> du;
    const auto f16 = hn::Floor(hn::Mul(s, hn::Set(d, 16.0)));
    const auto f = hn::Mul(f16, hn::Set(d, 1.0 / 16));
    auto red = hn::Sub(s, f);  // inexact only for -1/32 < s < 0 (then f = -1/16: Fast2Sum error below)
    const auto rerr = hn::Sub(s, hn::Add(red, f));
    red = hn::IfThenElse(hn::Lt(rerr, hn::Zero(d)), hn::BitCast(d, hn::Sub(hn::BitCast(du, red), hn::Set(du, 1))), red);
    red = hn::IfThenElse(hn::Eq(red, hn::Zero(d)), hn::Set(d, -0.0), red);  // x - x = -0 when rounding down
    sh = hn::MulAdd(f16, hn::Set(d, 1.0 / 16), hn::Set(d, kShifter));  // = RD(s + shifter), exact
    return red;
#endif
}

#if HWYMATH_NATIVE_AVX512
// vpermt2pd lookup of entries tab[idx & 15] (the instruction ignores the higher index bits, so no masking)
HWY_INLINE hn::Vec512<double> Permute16(const double *tab, hn::Vec512<uint64_t> idx) {
    return hn::Vec512<double>{_mm512_permutex2var_pd(_mm512_loadu_pd(tab), idx.raw, _mm512_loadu_pd(tab + 8))};
}
// 32-entry lookup as in the assembly: vpermt2pd on index bits 0-3 of both halves, blend on bit 4
HWY_INLINE hn::Vec512<double> Permute32(const double *tab, hn::Vec512<uint64_t> idx, hn::Mask512<double> hi) {
    return hn::IfThenElse(hi, Permute16(tab + 16, idx), Permute16(tab, idx));
}
// vcmppd _CMP_NLE_UQ: !(a <= b), true also for NaN (one compare, where Highway needs Le and Not)
HWY_INLINE hn::Mask512<double> NotLe(hn::Vec512<double> a, hn::Vec512<double> b) {
    return hn::Mask512<double>{_mm512_cmp_pd_mask(a.raw, b.raw, _CMP_NLE_UQ)};
}
#endif

template <class D>
HWY_INLINE V<D> PowCore(D d, V<D> x, V<D> y, hn::MFromD<D> &callout) {
    const hn::RebindToUnsigned<D> du;
#if HWYMATH_NATIVE_AVX512
    // Intel's instructions: vfpclasspd, vgetmantpd/vgetexppd {sae}
    namespace nat = hwymath::HWY_NAMESPACE;
    const auto special = hn::Or(nat::FpClass<0xdf>(x), nat::FpClass<0x99>(y));  // x <= 0, x or y inf/NaN
    const auto mant = nat::GetMant<0x0a>(x);  // to [1/2, 1) (NaN for negative x: callout lanes)
    auto expon = nat::GetExp(x);
#else
    const hn::RebindToSigned<D> di;
    const auto inf = hn::Inf(d);

    // vfpclasspd: x <= 0 (incl. -0 and negative denormals), x or y inf/NaN
    const auto x_ok = hn::And(hn::Gt(x, hn::Zero(d)), hn::Lt(x, inf));
    const auto y_ok = hn::Lt(hn::Abs(y), inf);

    // vgetmantpd (to [1/2, 1)) and vgetexppd, exact for positive normal and denormal x
    const auto den = hn::RebindMask(d, hn::Lt(hn::BitCast(di, x), hn::Set(di, int64_t{0x0010000000000000})));
    const auto xn = hn::BitCast(du, hn::IfThenElse(den, hn::Mul(x, hn::Set(d, 0x1p52)), x));
    const auto eb = hn::And(hn::ShiftRight<52>(xn), hn::Set(du, 0x7ff));
    const auto mant = hn::BitCast(d, hn::Or(hn::And(xn, hn::Set(du, 0x000fffffffffffffULL)), hn::Set(du, 0x3fe0000000000000ULL)));
    auto expon = hn::Sub(hn::BitCast(d, hn::Or(eb, hn::Set(du, 0x4330000000000000ULL))),
                         hn::IfThenElse(den, hn::Set(d, 0x1p52 + 1075), hn::Set(d, 0x1p52 + 1023)));
#endif

    const auto rcp = RcpRounded(d, mant);

    // log2(x)
#if HWYMATH_NATIVE_AVX512
    // vpermt2pd ignores the index bits above bit 3: no masking; blend on bit 4 (vpsllq + vpmovq2m)
    const auto idx = hn::ShiftRight<47>(hn::BitCast(du, rcp));
    const auto hi16 = hn::MaskFromVec(hn::BitCast(d, hn::ShiftLeft<59>(idx)));
    const auto th = Permute32(vLogTh, idx, hi16);
    expon = hn::IfThenElse(hn::Lt(rcp, hn::Set(d, kRcpLimit)), hn::Add(expon, hn::Set(d, kOne)), expon);
    const auto tl = Permute32(vLogTl, idx, hi16);
#else
    const auto idx = hn::And(hn::ShiftRight<47>(hn::BitCast(du, rcp)), hn::Set(du, 31));
    const auto th = Lookup<32>(d, vLogTh, kLogThPad.v, idx);
    expon = hn::IfThenElse(hn::Lt(rcp, hn::Set(d, kRcpLimit)), hn::Add(expon, hn::Set(d, kOne)), expon);
    const auto tl = Lookup<32>(d, vLogTl, kLogTlPad.v, idx);
#endif
    const auto mh = hn::Mul(mant, hn::Set(d, kHalf));
    const auto R = hn::MulSub(rcp, mh, hn::Set(d, kHalf));  // (rcp*m - 1)/2, exact

    const auto R2 = hn::Mul(R, R);
    const auto p910 = hn::MulAdd(R, hn::Set(d, kA9), hn::Set(d, kA8));
    const auto p78 = hn::MulAdd(R, hn::Set(d, kA7), hn::Set(d, kA6));
    const auto eth = hn::Add(th, expon);
    const auto p710 = hn::MulAdd(p910, R2, p78);
    const auto p12 = hn::MulAdd(R, hn::Set(d, kA1), hn::Set(d, kA0));
    const auto p56 = hn::MulAdd(R, hn::Set(d, kA5), hn::Set(d, kA4));
    const auto p34 = hn::MulAdd(R, hn::Set(d, kA3), hn::Set(d, kA2));
    const auto R4 = hn::Mul(R2, R2);
    const auto p36 = hn::MulAdd(p56, R2, p34);
    const auto p310 = hn::MulAdd(p710, R4, p36);
    const auto p110 = hn::MulAdd(p310, R2, p12);
    const auto tlp = hn::MulAdd(p110, R, tl);

    const auto c1h = hn::Set(d, kC1H);
    const auto s1 = hn::MulAdd(R, c1h, eth);
    const auto rh = hn::Sub(s1, eth);
    const auto s2 = hn::NegMulAdd(R, rh, s1);
    const auto rl = hn::MulSub(R, c1h, rh);
    const auto rrh = hn::Sub(s1, s2);
    const auto rl2 = hn::NegMulAdd(rl, R, rl);
    const auto rrl = hn::MulSub(rh, R, rrh);
    const auto low1 = hn::Sub(rl2, rrl);
    const auto low = hn::Add(tlp, low1);
    const auto high = hn::Add(s2, low);

    // y*log2(x), rounded toward zero ({rz-sae} on AVX-512, emulated elsewhere)
    const auto P = MulRZ(d, y, high);
    const auto tlh = hn::Sub(high, s2);
    const auto zl = hn::MulSub(y, high, P);  // exact (the RZ is a no-op)
    const auto tll = hn::Sub(low, tlh);
    const auto T = FmaRZ(d, y, tll, zl);
    const auto S = AddRZ(d, P, T);
    const auto sp = hn::Sub(S, P);

    // 2^(S + low)
    V<D> sh;
    const auto red = Reduce16(d, S, sh);
    const auto lw = hn::Sub(T, sp);
    auto re = hn::Add(red, lw);
    re = hn::BitCast(d, hn::And(hn::BitCast(du, re), hn::Set(du, kRMaskBits)));
    const auto shb = hn::BitCast(du, sh);
#if HWYMATH_NATIVE_AVX512
    const auto tle = Permute16(vExpTle, shb);
    const auto the = Permute16(vExpThe, shb);
#else
    const auto jdx = hn::And(shb, hn::Set(du, 15));
    const auto tle = Lookup<16>(d, vExpTle, kExpTlePad.v, jdx);
    const auto the = Lookup<16>(d, vExpThe, kExpThePad.v, jdx);
#endif
    const auto re2 = hn::Mul(re, re);
#if HWYMATH_NATIVE_AVX512
    const auto range = NotLe(hn::Abs(S), hn::Set(d, kBound));
#else
    const auto range = hn::Not(hn::Le(hn::Abs(S), hn::Set(d, kBound)));
#endif
    auto pq = hn::MulAdd(re, hn::Set(d, kB6), hn::Set(d, kB5));
    const auto q34 = hn::MulAdd(re, hn::Set(d, kB4), hn::Set(d, kB3));
    const auto q12 = hn::MulAdd(re, hn::Set(d, kB2), hn::Set(d, kB1));
    const auto sc = hn::BitCast(d, hn::And(hn::ShiftLeft<48>(shb), hn::Set(du, kExpMaskBits)));  // 2^N
    pq = hn::MulAdd(pq, re2, q34);
    pq = hn::MulAdd(pq, re2, q12);
    pq = hn::MulAdd(pq, re, tle);
    pq = hn::MulAdd(pq, the, the);

#if HWYMATH_NATIVE_AVX512
    callout = hn::Or(special, range);
#else
    callout = hn::Or(hn::Not(hn::And(x_ok, y_ok)), range);
#endif
    return hn::Mul(pq, sc);
}

void PowArray(const double *HWY_RESTRICT x, const double *HWY_RESTRICT y, double *HWY_RESTRICT z, size_t len) {
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    uint8_t bits[16];
    size_t i = 0;
    for (; i + N <= len; i += N) {
        hn::MFromD<decltype(d)> m;
        hn::StoreU(PowCore(d, hn::LoadU(d, x + i), hn::LoadU(d, y + i), m), d, z + i);
        if (HWY_UNLIKELY(!hn::AllFalse(d, m))) {
            hn::StoreMaskBits(d, m, bits);
            FixSpecial(x + i, y + i, z + i, N, bits);
        }
    }
    if (i < len) {
        const size_t rem = len - i;
        const auto one = hn::Set(d, 1.0);
        hn::MFromD<decltype(d)> m;
        hn::StoreN(PowCore(d, hn::LoadNOr(one, d, x + i, rem), hn::LoadNOr(one, d, y + i, rem), m), d, z + i, rem);
        if (!hn::AllFalse(d, m)) {
            hn::StoreMaskBits(d, m, bits);
            FixSpecial(x + i, y + i, z + i, rem, bits);
        }
    }
}

// the emulated AVX-512 operations on arrays, for the tests: op 0 = rcp, 1 = mul RZ, 2 = fma RZ, 3 = add RZ,
// 4 = reduce, 5 = RD(a + shifter). On AVX-512 targets these are the instructions themselves.
void Primitive(int op, const double *a, const double *b, const double *c, double *out, size_t len) {
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    for (size_t i = 0; i + N <= len; i += N) {
        const auto va = hn::LoadU(d, a + i), vb = hn::LoadU(d, b + i), vc = hn::LoadU(d, c + i);
        V<decltype(d)> r, sh;
        switch (op) {
            case 0: r = RcpRounded(d, va); break;
            case 1: r = MulRZ(d, va, vb); break;
            case 2: r = FmaRZ(d, va, vb, vc); break;
            case 3: r = AddRZ(d, va, vb); break;
            case 4: r = Reduce16(d, va, sh); break;
            default: Reduce16(d, va, sh); r = sh; break;
        }
        hn::StoreU(r, d, out + i);
    }
}

#else  // no fused multiply-add: the error-free transformations would not be exact

void PowArray(const double *HWY_RESTRICT x, const double *HWY_RESTRICT y, double *HWY_RESTRICT z, size_t len) {
    for (size_t i = 0; i < len; i++) z[i] = std::pow(x[i], y[i]);
}

void Primitive(int, const double *, const double *, const double *, double *, size_t) {}

#endif

size_t LanesF64() {
#if HWY_HAVE_FLOAT64
    return hn::Lanes(hn::ScalableTag<double>());
#else
    return 1;
#endif
}

int Variant() {  // 1: SVML algorithm (bit-identical), 0: std::pow fallback
#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64
    return 1;
#else
    return 0;
#endif
}

}  // namespace HWY_NAMESPACE
}  // namespace hwypow
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace hwypow {
HWY_EXPORT(PowArray);
HWY_EXPORT(Variant);
HWY_EXPORT(Primitive);
HWY_EXPORT(LanesF64);
// dispatch inside the namespace: with a single (static) target, HWY_DYNAMIC_DISPATCH(F) becomes N_<TARGET>::F
static void DispatchPow(const double *x, const double *y, double *z, size_t len) {
    HWY_DYNAMIC_DISPATCH(PowArray)(x, y, z, len);
}
static void DispatchPrimitive(int op, const double *a, const double *b, const double *c, double *out, size_t len) {
    HWY_DYNAMIC_DISPATCH(Primitive)(op, a, b, c, out, len);
}
static int DispatchVariant() { return HWY_DYNAMIC_DISPATCH(Variant)(); }
static size_t DispatchLanes() { return HWY_DYNAMIC_DISPATCH(LanesF64)(); }
}  // namespace hwypow

extern "C" {
void hwy_pow_array(const double *x, const double *y, double *z, size_t len) { hwypow::DispatchPow(x, y, z, len); }
// for the tests: the emulated AVX-512 operations (len a multiple of the vector length)
void hwy_pow_primitive(int op, const double *a, const double *b, const double *c, double *out, size_t len) {
    hwypow::DispatchPrimitive(op, a, b, c, out, len);
}
int hwy_pow_variant(void) { return hwypow::DispatchVariant(); }
size_t hwy_pow_lanes(void) { return hwypow::DispatchLanes(); }
}
#endif  // HWY_ONCE
