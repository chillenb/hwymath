// AVX-512 instructions that Highway does not wrap, for the AVX-512 targets of the Highway versions: on those
// targets the functions use the instructions Intel's code uses instead of their exact emulations (which run on
// every other target and cost several instructions each: AVX-512 code is usually limited by its two vector
// ports, so instruction count is time). Per-target header (re-included by foreach_target.h).
//
// Usage, inside HWY_NAMESPACE:
//   #if HWYMATH_NATIVE_AVX512
//       y = nat::Scalef(y, n);            // nat = hwymath::HWY_NAMESPACE (or hwymath::HWY_NAMESPACE::)
//   #else
//       y = <the exact emulation>;
//   #endif
// HWYMATH_NATIVE_AVX512 is 1 on AVX-512 targets (HWY_TARGET <= HWY_AVX3: AVX-512 F, VL, DQ, BW), where
// ScalableTag<double> vectors are Vec512<double>; the functions take and return Vec512<double> (and
// Mask512<double>) only. Rounding-mode and {sae} operations suppress exceptions, as Intel's code does.

#if defined(HWYMATH_AVX512_NATIVE_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef HWYMATH_AVX512_NATIVE_INL_H_
#undef HWYMATH_AVX512_NATIVE_INL_H_
#else
#define HWYMATH_AVX512_NATIVE_INL_H_
#endif

#include "hwy/highway.h"

#undef HWYMATH_NATIVE_AVX512
#if HWY_TARGET <= HWY_AVX3
#define HWYMATH_NATIVE_AVX512 1
#else
#define HWYMATH_NATIVE_AVX512 0
#endif

HWY_BEFORE_NAMESPACE();
namespace hwymath {
namespace HWY_NAMESPACE {

#if HWYMATH_NATIVE_AVX512
namespace hn = hwy::HWY_NAMESPACE;
using V512d = hn::Vec512<double>;
using M512d = hn::Mask512<double>;

// rounding of the *_round operations (exceptions suppressed, as {rn,rd,ru,rz}-sae)
constexpr int kRN = _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC;
constexpr int kRD = _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC;
constexpr int kRU = _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC;
constexpr int kRZ = _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC;

template <int kR> HWY_INLINE V512d AddR(V512d a, V512d b) { return V512d{_mm512_add_round_pd(a.raw, b.raw, kR)}; }
template <int kR> HWY_INLINE V512d SubR(V512d a, V512d b) { return V512d{_mm512_sub_round_pd(a.raw, b.raw, kR)}; }
template <int kR> HWY_INLINE V512d MulR(V512d a, V512d b) { return V512d{_mm512_mul_round_pd(a.raw, b.raw, kR)}; }
template <int kR> HWY_INLINE V512d DivR(V512d a, V512d b) { return V512d{_mm512_div_round_pd(a.raw, b.raw, kR)}; }
// a*b + c, a*b - c, -(a*b) + c rounded once with kR
template <int kR> HWY_INLINE V512d MulAddR(V512d a, V512d b, V512d c) {
    return V512d{_mm512_fmadd_round_pd(a.raw, b.raw, c.raw, kR)};
}
template <int kR> HWY_INLINE V512d MulSubR(V512d a, V512d b, V512d c) {
    return V512d{_mm512_fmsub_round_pd(a.raw, b.raw, c.raw, kR)};
}
template <int kR> HWY_INLINE V512d NegMulAddR(V512d a, V512d b, V512d c) {
    return V512d{_mm512_fnmadd_round_pd(a.raw, b.raw, c.raw, kR)};
}

// vreducepd x, imm: x - round_to_2^-M(x) (imm = M << 4 | rounding)
template <int kImm> HWY_INLINE V512d Reduce(V512d x) { return V512d{_mm512_reduce_pd(x.raw, kImm)}; }
// vrndscalepd x, imm: round_to_2^-M(x)
template <int kImm> HWY_INLINE V512d RoundScale(V512d x) { return V512d{_mm512_roundscale_pd(x.raw, kImm)}; }
// vscalefpd y, n: y * 2^floor(n)
HWY_INLINE V512d Scalef(V512d y, V512d n) { return V512d{_mm512_scalef_pd(y.raw, n.raw)}; }
// vgetexppd {sae}: floor(log2|x|) as a double
HWY_INLINE V512d GetExp(V512d x) { return V512d{_mm512_getexp_round_pd(x.raw, _MM_FROUND_NO_EXC)}; }
// vgetmantpd {sae}, imm = sign control << 2 | interval
template <int kImm> HWY_INLINE V512d GetMant(V512d x) {
    return V512d{_mm512_getmant_round_pd(x.raw, static_cast<_MM_MANTISSA_NORM_ENUM>(kImm & 3),
                                         static_cast<_MM_MANTISSA_SIGN_ENUM>(kImm >> 2), _MM_FROUND_NO_EXC)};
}
// vfpclasspd x, imm: the categories selected by imm
template <int kImm> HWY_INLINE M512d FpClass(V512d x) { return M512d{_mm512_fpclass_pd_mask(x.raw, kImm)}; }
// vcmppd a, b, predicate (_CMP_*): e.g. _CMP_NLT_UQ / _CMP_NLE_UQ in one instruction (Highway: compare + knot)
template <int kPred> HWY_INLINE M512d Cmp(V512d a, V512d b) { return M512d{_mm512_cmp_pd_mask(a.raw, b.raw, kPred)}; }
// vpermpd: table[idx & 7]; vpermt2pd: (lo, hi)[idx & 15]. The instructions read only the low index bits, so
// no masking is needed (Highway's TableLookupLanes/TwoTablesLookupLanes require in-range indices). VI: any
// Vec512 of 64-bit integers.
template <class VI> HWY_INLINE V512d PermuteVar(V512d table, VI idx) {
    return V512d{_mm512_permutexvar_pd(idx.raw, table.raw)};
}
template <class VI> HWY_INLINE V512d Permute2Var(V512d lo, VI idx, V512d hi) {
    return V512d{_mm512_permutex2var_pd(lo.raw, idx.raw, hi.raw)};
}
// 16-entry table in memory: tab[idx & 15]
template <class VI> HWY_INLINE V512d Lookup16Raw(const double *tab, VI idx) {
    return V512d{_mm512_permutex2var_pd(_mm512_loadu_pd(tab), idx.raw, _mm512_loadu_pd(tab + 8))};
}
#endif  // HWYMATH_NATIVE_AVX512

}  // namespace HWY_NAMESPACE
}  // namespace hwymath
HWY_AFTER_NAMESPACE();

#endif  // HWYMATH_AVX512_NATIVE_INL_H_
