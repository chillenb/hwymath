// Cross-platform check of the Highway version of one log-family function, built for one static target
// (-DHWY_COMPILE_ONLY_STATIC): writes the inputs and results to files (compared bit for bit between
// targets/architectures by tools/xcheck.sh) and prints the accuracy vs the long double reference (binary128
// on AArch64, x87 extended on x86-64). Included by functions/<f>/tests/hwy_xcheck.cc, which defines FN_NAME,
// REF(x), HWY_ARRAY, VARIANT, LANES.
// usage: hwy_xcheck OUTPUT_PREFIX [N_RANDOM]       Build with -ffp-contract=off (identical inputs).
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "hwy/targets.h"
#include "hwymath/accuracy.h"

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s OUTPUT_PREFIX [N_RANDOM]\n", argv[0]); return 2; }
    size_t n_random = argc > 2 ? strtoull(argv[2], nullptr, 10) : 2000000;
    std::vector<double> x(n_random + 40000);
    x.resize(make_inputs(x.data(), n_random));
    std::vector<double> y(x.size());
    HWY_ARRAY(x.data(), y.data(), x.size());

    std::string pre = argv[1];
    FILE *f = fopen((pre + ".in").c_str(), "wb"); fwrite(x.data(), 8, x.size(), f); fclose(f);
    f = fopen((pre + ".out").c_str(), "wb"); fwrite(y.data(), 8, y.size(), f); fclose(f);

    acc_stats s = ACC_STATS_INIT;
    for (size_t i = 0; i < x.size(); i++) acc_add(&s, x[i], y[i], REF((long double)x[i]));
    printf("%-6s %-10s %2zu lanes  %-5s %zu inputs: ", FN_NAME, hwy::TargetName(HWY_STATIC_TARGET), LANES(),
           VARIANT() ? "SVML" : "libm", x.size());
    acc_print(stdout, "", &s);
    return 0;
}
