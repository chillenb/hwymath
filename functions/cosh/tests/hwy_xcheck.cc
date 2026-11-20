// Cross-platform check of the Highway cosh, built for one static target (-DHWY_COMPILE_ONLY_STATIC): writes
// the inputs and results to files (compared bit for bit between targets/architectures by
// tools/xcheck.sh) and prints the accuracy vs coshl (binary128 on AArch64, x87 extended on x86-64).
// usage: hwy_xcheck OUTPUT_PREFIX [N_RANDOM]       Build with -ffp-contract=off (identical inputs).
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "hwy/targets.h"

extern "C" {
void hwy_cosh_array(const double *x, double *y, size_t len);
int hwy_cosh_variant(void);
size_t hwy_cosh_lanes(void);
}
#include "hwymath/accuracy.h"
#include "inputs.h"

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s OUTPUT_PREFIX [N_RANDOM]\n", argv[0]); return 2; }
    const size_t n_random = argc > 2 ? strtoull(argv[2], nullptr, 10) : 2000000;
    std::vector<double> x(n_random + 8192);
    x.resize(make_inputs(x.data(), n_random));
    std::vector<double> y(x.size());
    hwy_cosh_array(x.data(), y.data(), x.size());

    const std::string pre = argv[1];
    FILE *f = fopen((pre + ".in").c_str(), "wb");
    fwrite(x.data(), 8, x.size(), f);
    fclose(f);
    f = fopen((pre + ".out").c_str(), "wb");
    fwrite(y.data(), 8, y.size(), f);
    fclose(f);

    acc_stats s = ACC_STATS_INIT;
    for (size_t i = 0; i < x.size(); i++) acc_add(&s, x[i], y[i], coshl((long double)x[i]));
    char name[64];
    snprintf(name, sizeof name, "%s, %zu lanes, %s", hwy::TargetName(HWY_STATIC_TARGET), hwy_cosh_lanes(),
             hwy_cosh_variant() ? "svml" : "std::cosh");
    acc_print(stdout, name, &s);
    return 0;
}
