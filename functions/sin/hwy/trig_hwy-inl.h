// Highway building blocks shared by the ports of SVML's high-accuracy double sin, cos and tan
// (functions/{sin,cos,tan}/hwy/hwy_*.cc): exact emulations of the AVX-512 operations their algorithms
// use, so that the results are bit-identical to Intel's AVX-512 code on every target with a true FMA.
//   - Lookup<T>: vpermt2pd/vpermi2d-style lookup of a T-entry table (T = 16 or 64) with in-register
//     permutes: one TableLookupLanes if a vector holds the whole table, otherwise a tree of
//     TwoTablesLookupLanes and selects on the index bits (no gathers: they are slow with the Gather Data
//     Sampling mitigation on x86); scalable vectors shorter than the table use GatherIndex (SVE).
//   - TrigReduceLarge: the Payne-Hanek reduction of the large-argument path (vpmullq of 32-bit values =
//     MulEven of the 32-bit halves, exact; the table lookups by exponent, done with vgatherdpd in the
//     assembly, become scalar loads: the path is rare), TrigReduced (|x| <= 2^-20 pass-through),
//     SinCosLargeTail (evaluation of the sin/cos large path).
// All floating-point operations round to nearest, like the {rn-sae} instructions of the assembly; with
// -ffp-contract=off the compiler fuses nothing else. See functions/sin/svml/svml_trig_reduction.h for
// the algorithm. Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.

#include "functions/sin/hwy/trig_tables.h"

// Per-target code (re-included by foreach_target.h for every target).
#if defined(HWYMATH_TRIG_HWY_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef HWYMATH_TRIG_HWY_INL_H_
#undef HWYMATH_TRIG_HWY_INL_H_
#else
#define HWYMATH_TRIG_HWY_INL_H_
#endif

#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace hwymath_trig {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64

// ---------------------------------------------------------------------------------------------------
// Table lookups with permutes

// tab[idx] for kCount entries starting at tab, kN lanes per vector (compile-time), idx < kCount
template <size_t kCount, size_t kN, class D, class VI>
HWY_INLINE hn::VFromD<D> LookupTree(D d, const hn::TFromD<D> *tab, VI idx) {
    const hn::RebindToSigned<D> di;
    if constexpr (kCount <= 2 * kN) {
        static_assert(kCount == 2 * kN, "table smaller than two vectors");
        const auto i = hn::And(hn::BitCast(di, idx), hn::Set(di, static_cast<hn::TFromD<decltype(di)>>(kCount - 1)));
        return hn::TwoTablesLookupLanes(d, hn::LoadU(d, tab), hn::LoadU(d, tab + kN), hn::IndicesFromVec(d, i));
    } else {
        const auto lo = LookupTree<kCount / 2, kN>(d, tab, idx);
        const auto hi = LookupTree<kCount / 2, kN>(d, tab + kCount / 2, idx);
        // select on index bit log2(kCount/2), moved into the sign bit (a single blendv on AVX2)
        constexpr int kShift = static_cast<int>(sizeof(hn::TFromD<D>) * 8 - 1) - hwy::FloorLog2(kCount / 2);
        const auto sel = hn::ShiftLeft<kShift>(hn::BitCast(di, idx));
        return hn::BitCast(d, hn::IfNegativeThenElse(sel, hn::BitCast(di, hi), hn::BitCast(di, lo)));
    }
}

// tab[idx mod kT] (kT a power of 2); pad = the table repeated to 64 entries
template <size_t kT, class D, class VI>
HWY_INLINE hn::VFromD<D> Lookup(D d, const hn::TFromD<D> *tab, const hn::TFromD<D> *pad, VI idx) {
    const hn::RebindToSigned<D> di;
    using TI = hn::TFromD<decltype(di)>;
    const auto i = hn::And(hn::BitCast(di, idx), hn::Set(di, static_cast<TI>(kT - 1)));
#if HWY_HAVE_SCALABLE
    const size_t N = hn::Lanes(d);
    if (N >= kT && N <= 64) return hn::TableLookupLanes(hn::LoadU(d, pad), hn::IndicesFromVec(d, i));
    return hn::GatherIndex(d, tab, i);
#else
    constexpr size_t kN = hn::MaxLanes(d);
    if constexpr (kN >= kT) {
        static_assert(kN <= 64, "vector longer than the padded table");
        return hn::TableLookupLanes(hn::LoadU(d, pad), hn::IndicesFromVec(d, i));
    } else if constexpr (kN >= 2) {
        (void)pad;
        return LookupTree<kT, kN>(d, tab, i);
    } else {
        (void)pad;
        return hn::GatherIndex(d, tab, i);
    }
#endif
}

// base[idx[i]] for every lane: scalar loads (for the rare large-argument path)
template <class D, class VU>
HWY_INLINE hn::VFromD<D> LoadIndexed(D d, const hn::TFromD<D> *base, VU idx) {
    const hn::RebindToUnsigned<D> du;
    HWY_ALIGN uint64_t ix[hn::MaxLanes(du)];
    HWY_ALIGN hn::TFromD<D> buf[hn::MaxLanes(d)];
    hn::Store(idx, du, ix);
    for (size_t i = 0; i < hn::Lanes(d); i++) buf[i] = base[ix[i]];
    return hn::Load(d, buf);
}

// -(mul*x) - sub with one rounding, as vfnmsub*pd computes it. Highway's NegMulSub is Neg(MulAdd(mul, x,
// sub)) on NEON, which gives -0 instead of +0 when mul*x + sub is exactly 0 (e.g. tan(-0) became +0);
// -(mul*x) + (-sub) is the same operation in IEEE arithmetic on every target.
template <class V>
HWY_INLINE V NegMulSubExact(V mul, V x, V sub) {
    return hn::NegMulAdd(mul, x, hn::Neg(sub));
}

// ---------------------------------------------------------------------------------------------------
// Large-argument path

// 32x32 -> 64-bit product of the low halves of 64-bit lanes (vpmullq of values < 2^32)
template <class VU>
HWY_INLINE VU Mul32(VU a, VU b) {
    const hn::DFromV<VU> du;
    const hn::Repartition<uint32_t, decltype(du)> d32;
    return hn::MulEven(hn::BitCast(d32, a), hn::BitCast(d32, b));
}

// Payne-Hanek reduction of xs: x = N*2pi/512 + phi + plo, N in [0, 511] (trig_reduce_avx512)
template <class D>
HWY_INLINE void TrigReduceLarge(D d, hn::VFromD<D> xs, hn::VFromD<D> &phi, hn::VFromD<D> &plo,
                                hn::VFromD<hn::RebindToUnsigned<D>> &n) {
    const hn::RebindToUnsigned<D> du;
    using VU = hn::VFromD<decltype(du)>;
    const VU xb = hn::BitCast(du, xs);
    const VU lo32 = hn::Set(du, trig_red_Low32);

    const VU e = hn::ShiftRight<52>(hn::And(xb, hn::Set(du, trig_red_ExpMask)));
    const VU e3 = hn::Add(e, hn::ShiftLeft<1>(e));
    const VU t0 = LoadIndexed(du, trig_red_Tab, e3);
    const VU t1 = LoadIndexed(du, trig_red_Tab + 1, e3);
    const VU t2 = LoadIndexed(du, trig_red_Tab + 2, e3);
    const VU t0h = hn::ShiftRight<32>(t0), t1h = hn::ShiftRight<32>(t1), t2h = hn::ShiftRight<32>(t2);
    const VU t0l = hn::And(t0, lo32), t1l = hn::And(t1, lo32), t2l = hn::And(t2, lo32);

    const VU m = hn::Add(hn::And(xb, hn::Set(du, trig_red_MantMask)), hn::Set(du, trig_red_ImplicitBit));
    const VU mh = hn::ShiftRight<32>(m), ml = hn::And(m, lo32);

    const VU mh_t1h = Mul32(mh, t1h), mh_t0l = Mul32(mh, t0l), mh_t1l = Mul32(mh, t1l);
    const VU mh_t2h = Mul32(mh, t2h), mh_t2l = Mul32(mh, t2l), ml_t0l = Mul32(ml, t0l);
    const VU ml_t1l = Mul32(ml, t1l), ml_t2h = Mul32(ml, t2h), ml_t0h = Mul32(ml, t0h), ml_t1h = Mul32(ml, t1h);
    const VU s1 = hn::Add(mh_t2h, hn::ShiftRight<32>(ml_t2h));
    const VU s2 = hn::Add(mh_t1l, hn::ShiftRight<32>(ml_t1l));
    const VU s3 = hn::Add(mh_t0l, hn::ShiftRight<32>(ml_t0l));
    const VU s4 = hn::Add(mh_t1h, hn::ShiftRight<32>(ml_t1h));
    const VU s5 = hn::Add(hn::And(ml_t0h, lo32), s3);
    const VU s6 = hn::Add(hn::ShiftRight<32>(mh_t2l), hn::Add(hn::And(ml_t1l, lo32), s1));
    const VU s7 = hn::Add(hn::ShiftRight<32>(s6), hn::Add(hn::And(ml_t1h, lo32), s2));
    const VU s8 = hn::Add(hn::ShiftRight<32>(s7), hn::Add(s4, hn::And(ml_t0l, lo32)));
    const VU s9 = hn::Add(s5, hn::ShiftRight<32>(s8));
    const VU wlo = hn::Add(hn::ShiftLeft<32>(s7), hn::And(s6, lo32));
    const VU whi = hn::Add(hn::ShiftLeft<32>(s9), hn::And(s8, lo32));

    const VU sign = hn::And(xb, hn::Set(du, trig_red_SignMask));
    const auto fhi = hn::BitCast(d, hn::Or(hn::Xor(sign, hn::Set(du, trig_red_OneBits)), hn::ShiftRight<12>(whi)));
    const auto shifter = hn::Set(d, trig_red_Shifter);
    const auto sh = hn::Add(shifter, fhi);
    auto rem = hn::Sub(fhi, hn::Sub(sh, shifter));

    const VU two104 = hn::Xor(sign, hn::Set(du, trig_red_Two104Bits));
    const VU two52 = hn::Xor(sign, hn::Set(du, trig_red_Two52Bits));
    const auto pl = hn::Sub(hn::BitCast(d, hn::Or(hn::ShiftLeft<28>(hn::And(wlo, hn::Set(du, trig_red_Low24))), two104)),
                            hn::BitCast(d, two104));
    const auto pm = hn::Sub(
            hn::BitCast(d, hn::Or(hn::Or(hn::ShiftLeft<40>(hn::And(whi, hn::Set(du, trig_red_Low12))),
                                         hn::ShiftRight<24>(wlo)),
                                  two52)),
            hn::BitCast(d, two52));

    const auto twopi_hi = hn::Set(d, trig_red_TwoPiHi);
    const auto rhi = hn::Add(rem, pm);
    phi = hn::Mul(rhi, twopi_hi);
    rem = hn::Sub(rem, rhi);
    const auto perr = hn::MulSub(rhi, twopi_hi, phi);
    const auto rlo = hn::Add(hn::Add(pm, rem), pl);
    const auto t = hn::MulAdd(hn::Set(d, trig_red_TwoPiLo), rhi, perr);
    plo = hn::MulAdd(twopi_hi, rlo, t);
    n = hn::And(hn::BitCast(du, sh), hn::Set(du, trig_red_Mask9));
}

// reduced argument r + rlo: phi + plo, or xs itself (and 0) if |x| <= 2^-20 (trig_reduced_avx512)
template <class D>
HWY_INLINE void TrigReduced(D d, hn::VFromD<D> xs, hn::VFromD<D> phi, hn::VFromD<D> plo, hn::VFromD<D> &r,
                            hn::VFromD<D> &rlo) {
    const auto tiny = hn::Le(hn::Abs(xs), hn::Set(d, trig_red_Tiny));  // ordered: false for NaN
    const auto not_big = hn::Not(hn::Gt(hn::Abs(xs), hn::Set(d, trig_red_Tiny)));
    r = hn::Or(hn::IfThenZeroElse(not_big, phi), hn::IfThenElseZero(tiny, xs));
    rlo = hn::IfThenZeroElse(not_big, plo);
}

// evaluation tail of the sin/cos large path (trig_sincos_large_tail_avx512); tab: 512 rows {T0..T3}
template <class D, class VU>
HWY_INLINE hn::VFromD<D> SinCosLargeTail(D d, hn::VFromD<D> r, hn::VFromD<D> rlo, VU n, const double *tab, double s1,
                                         double s2, double c1, double c2, double c3) {
    const VU row = hn::ShiftLeft<2>(n);
    const auto t0 = LoadIndexed(d, tab, row);
    const auto t1 = LoadIndexed(d, tab + 1, row);
    const auto t2 = LoadIndexed(d, tab + 2, row);
    const auto t3 = LoadIndexed(d, tab + 3, row);
    const auto r2 = hn::Mul(r, r);
    const auto a = hn::MulAdd(t1, r, t2);
    const auto hi = hn::MulAdd(r, t0, a);
    auto e1 = hn::Sub(t2, a);
    auto e2 = hn::Sub(a, hi);
    e1 = hn::MulAdd(t1, r, e1);
    e2 = hn::MulAdd(r, t0, e2);
    const auto sp = hn::MulAdd(hn::Set(d, s2), r2, hn::Set(d, s1));
    const auto err = hn::Add(e2, e1);
    const auto ct = hn::Add(t0, t1);
    const auto sp3 = hn::Mul(hn::Mul(sp, r2), r);
    const auto corr = hn::NegMulAdd(t2, r, ct);
    auto lo = hn::MulAdd(ct, sp3, err);
    auto cp = hn::MulAdd(hn::Set(d, c3), r2, hn::Set(d, c2));
    cp = hn::MulAdd(r2, cp, hn::Set(d, c1));
    cp = hn::Mul(cp, r2);
    lo = hn::MulAdd(t2, cp, lo);
    const auto lo2 = hn::MulAdd(corr, rlo, t3);
    return hn::Add(hi, hn::Add(lo, lo2));
}

#endif  // HWY_NATIVE_FMA && HWY_HAVE_FLOAT64

}  // namespace HWY_NAMESPACE
}  // namespace hwymath_trig
HWY_AFTER_NAMESPACE();

#endif  // HWYMATH_TRIG_HWY_INL_H_
