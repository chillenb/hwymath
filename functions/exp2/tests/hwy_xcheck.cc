// Cross-platform check of the Highway exp2, built for one static target (-DHWY_COMPILE_ONLY_STATIC):
// writes the inputs and results to files (compare them bit for bit between targets/architectures)
// and prints accuracy vs exp2l (binary128 long double on AArch64, x87 extended on x86-64).
// usage: hwy_xcheck OUTPUT_PREFIX [N_RANDOM]       Build with -ffp-contract=off (identical inputs).
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "hwy/targets.h"

extern "C" {
void hwy_exp2_array(const double *x, double *y, size_t len);
int hwy_exp2_variant(void);
size_t hwy_exp2_lanes(void);
}
#include "hwymath/accuracy.h"
#include "inputs.h"

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s OUTPUT_PREFIX [N_RANDOM]\n", argv[0]); return 2; }
    size_t n_random = argc > 2 ? strtoull(argv[2], nullptr, 10) : 2000000;
    std::vector<double> x(n_random + INPUTS_EXTRA);
    x.resize(make_inputs(x.data(), n_random));
    std::vector<double> y(x.size());
    hwy_exp2_array(x.data(), y.data(), x.size());

    std::string pre = argv[1];
    FILE *f = fopen((pre + ".in").c_str(), "wb"); fwrite(x.data(), 8, x.size(), f); fclose(f);
    f = fopen((pre + ".out").c_str(), "wb"); fwrite(y.data(), 8, y.size(), f); fclose(f);

    acc_stats s = ACC_STATS_INIT;
    for (size_t i = 0; i < x.size(); i++) acc_add(&s, x[i], y[i], exp2l((long double)x[i]));
    char name[64];
    snprintf(name, sizeof name, "%s %zu lanes %s", hwy::TargetName(HWY_STATIC_TARGET), hwy_exp2_lanes(),
             hwy_exp2_variant() ? "svml" : "std::exp2");
    acc_print(stdout, name, &s);
    return 0;
}
