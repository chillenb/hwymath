// Array versions of the Highway emulations in hwymath/hwy/avx512_emul-inl.h, for test_approx14.
#include <cstddef>

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "common/tests/approx14/hwy_approx14.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"
#include "hwymath/hwy/avx512_emul-inl.h"

HWY_BEFORE_NAMESPACE();
namespace hwyapprox14 {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

template <int kOp>
void Apply(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    for (size_t i = 0; i < len; i += N) {
        const size_t n = len - i < N ? len - i : N;
        const auto v = hn::LoadN(d, x + i, n);
        hn::VFromD<decltype(d)> r;
        if constexpr (kOp == 0) r = hwymath::HWY_NAMESPACE::Rcp14(d, v);
        else if constexpr (kOp == 1) r = hwymath::HWY_NAMESPACE::Rsqrt14(d, v);
        else r = hwymath::HWY_NAMESPACE::GetExpPositiveNormal(d, v);
        hn::StoreN(r, d, y + i, n);
    }
#else
    (void)x, (void)y, (void)len;
#endif
}
void Rcp14Array(const double *x, double *y, size_t n) { Apply<0>(x, y, n); }
void Rsqrt14Array(const double *x, double *y, size_t n) { Apply<1>(x, y, n); }
void GetExpArray(const double *x, double *y, size_t n) { Apply<2>(x, y, n); }
int Emulated() { return HWY_NATIVE_FMA && HWY_HAVE_FLOAT64; }

}  // namespace HWY_NAMESPACE
}  // namespace hwyapprox14
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace hwyapprox14 {
HWY_EXPORT(Rcp14Array);
HWY_EXPORT(Rsqrt14Array);
HWY_EXPORT(GetExpArray);
HWY_EXPORT(Emulated);
static void D0(const double *x, double *y, size_t n) { HWY_DYNAMIC_DISPATCH(Rcp14Array)(x, y, n); }
static void D1(const double *x, double *y, size_t n) { HWY_DYNAMIC_DISPATCH(Rsqrt14Array)(x, y, n); }
static void D2(const double *x, double *y, size_t n) { HWY_DYNAMIC_DISPATCH(GetExpArray)(x, y, n); }
static int D3() { return HWY_DYNAMIC_DISPATCH(Emulated)(); }
}  // namespace hwyapprox14
extern "C" {
void hwy_rcp14_array(const double *x, double *y, size_t n) { hwyapprox14::D0(x, y, n); }
void hwy_rsqrt14_array(const double *x, double *y, size_t n) { hwyapprox14::D1(x, y, n); }
void hwy_getexp_array(const double *x, double *y, size_t n) { hwyapprox14::D2(x, y, n); }
int hwy_approx14_emulated(void) { return hwyapprox14::D3(); }
}
#endif
