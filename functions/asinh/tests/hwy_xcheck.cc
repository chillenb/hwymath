// Cross-platform check of the Highway asinh/acosh/atanh, built for one static target
// (-DHWY_COMPILE_ONLY_STATIC) with -DXFN=asinh|acosh|atanh: writes the inputs and results to files
// (compare them bit for bit between targets/architectures) and prints the accuracy vs the long double
// function (binary128 on AArch64, x87 extended on x86-64).
// usage: hwy_xcheck OUTPUT_PREFIX [N_RANDOM]       Build with -ffp-contract=off (identical inputs).
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "hwy/targets.h"

#define CAT3(a, b, c) a##b##c
#define XCAT3(a, b, c) CAT3(a, b, c)
#define STR(a) #a
#define XSTR(a) STR(a)
#define INPUTS_H XSTR(functions/XFN/tests/inputs.h)

extern "C" {
void XCAT3(hwy_, XFN, _array)(const double *x, double *y, size_t len);
int XCAT3(hwy_, XFN, _variant)(void);
size_t XCAT3(hwy_, XFN, _lanes)(void);
#include "hwymath/accuracy.h"
}
#include INPUTS_H

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s OUTPUT_PREFIX [N_RANDOM]\n", argv[0]); return 2; }
    const size_t n_random = argc > 2 ? strtoull(argv[2], nullptr, 10) : 2000000;
    std::vector<double> x(n_random + 16384);
    x.resize(make_inputs(x.data(), n_random));
    std::vector<double> y(x.size());
    XCAT3(hwy_, XFN, _array)(x.data(), y.data(), x.size());

    const std::string pre = argv[1];
    FILE *f = fopen((pre + ".in").c_str(), "wb"); fwrite(x.data(), 8, x.size(), f); fclose(f);
    f = fopen((pre + ".out").c_str(), "wb"); fwrite(y.data(), 8, y.size(), f); fclose(f);

    acc_stats s = ACC_STATS_INIT;
    for (size_t i = 0; i < x.size(); i++) acc_add(&s, x[i], y[i], XCAT3(, XFN, l)((long double)x[i]));
    printf("%-10s %2zu lanes  %-12s ", hwy::TargetName(HWY_STATIC_TARGET), XCAT3(hwy_, XFN, _lanes)(),
           XCAT3(hwy_, XFN, _variant)() ? "SVML" : "std::" XSTR(XFN));
    acc_print(stdout, XSTR(XFN), &s);
    return 0;
}
