// Cross-platform check of the Highway asin/acos/atan/atan2, built for one static target
// (-DHWY_COMPILE_ONLY_STATIC) and one function (-DXC_F=asin|acos|atan|atan2, -I functions/<f>/tests):
// writes the inputs and results to files (compared bit for bit between targets and architectures by
// functions/asin/tools/xcheck_invtrig.sh) and prints the accuracy vs the long double function (binary128
// on AArch64, x87 extended on x86-64).
// usage: xcheck OUTPUT_PREFIX [N_RANDOM]       Build with -ffp-contract=off (identical inputs).
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "hwy/targets.h"

#define XC_CAT2(a, b) a##b
#define XC_CAT(a, b) XC_CAT2(a, b)
#define XC_STR2(a) #a
#define XC_STR(a) XC_STR2(a)
#define XC_HWY(suffix) XC_CAT(XC_CAT(hwy_, XC_F), suffix)

extern "C" {
#if XC_ARGS == 2
void XC_HWY(_array)(const double *y, const double *x, double *r, size_t len);
#else
void XC_HWY(_array)(const double *x, double *y, size_t len);
#endif
int XC_HWY(_variant)(void);
size_t XC_HWY(_lanes)(void);
}
#include "inputs.h"

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s OUTPUT_PREFIX [N_RANDOM]\n", argv[0]); return 2; }
    const size_t n_random = argc > 2 ? strtoull(argv[2], nullptr, 10) : 200000;
    std::vector<double> x(n_random + 16384), r;
#if XC_ARGS == 2
    std::vector<double> y(n_random + 16384);
    const size_t n = make_inputs(y.data(), x.data(), n_random);
    x.resize(n);
    y.resize(n);
    r.resize(n);
    XC_HWY(_array)(y.data(), x.data(), r.data(), n);
#else
    const size_t n = make_inputs(x.data(), n_random);
    x.resize(n);
    r.resize(n);
    XC_HWY(_array)(x.data(), r.data(), n);
#endif

    const std::string pre = argv[1];
    FILE *f = fopen((pre + ".in").c_str(), "wb");
#if XC_ARGS == 2
    fwrite(y.data(), 8, n, f);
#endif
    fwrite(x.data(), 8, n, f);
    fclose(f);
    f = fopen((pre + ".out").c_str(), "wb");
    fwrite(r.data(), 8, n, f);
    fclose(f);

    double max_err = 0;
    size_t cnt = 0, not_cr = 0;
    for (size_t i = 0; i < n; i++) {
#if XC_ARGS == 2
        const long double ref = XC_CAT(XC_F, l)((long double)y[i], (long double)x[i]);
#else
        const long double ref = XC_CAT(XC_F, l)((long double)x[i]);
#endif
        const double rd = (double)ref;
        if (!std::isfinite(rd) || rd == 0.0) continue;
        long double ulp = (long double)std::nextafter(std::fabs(rd), INFINITY) - fabsl((long double)rd);
        const double e = (double)(fabsl((long double)r[i] - ref) / ulp);
        cnt++;
        not_cr += e > 0.5;
        if (e > max_err) max_err = e;
    }
    printf("%-6s %-10s %2zu lanes  %s  %zu inputs: max %.4f ulp, not correctly rounded %.4f%%  (long double: %d-bit mantissa)\n",
           XC_STR(XC_F), hwy::TargetName(HWY_STATIC_TARGET), XC_HWY(_lanes)(), XC_HWY(_variant)() ? "SVML port" : "std lib  ",
           n, max_err, 100.0 * not_cr / cnt, LDBL_MANT_DIG);
    return 0;
}
