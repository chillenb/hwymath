// Cross-platform check of the Highway pow, built for one static target (-DHWY_COMPILE_ONLY_STATIC): writes
// the inputs and results to files (compare them bit for bit between targets/architectures) and prints the
// accuracy vs powl (binary128 long double on AArch64, x87 extended on x86-64). NaN results are written as
// one canonical NaN (the default NaN of 0/0 is negative on x86-64, positive on AArch64).
// usage: hwy_xcheck OUTPUT_PREFIX [N_RANDOM]       Build with -ffp-contract=off (identical inputs).
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "hwy/targets.h"

extern "C" {
#include "hwymath/accuracy.h"
void hwy_pow_array(const double *x, const double *y, double *z, size_t len);
int hwy_pow_variant(void);
size_t hwy_pow_lanes(void);
}
#include "inputs.h"

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s OUTPUT_PREFIX [N_RANDOM]\n", argv[0]); return 2; }
    const size_t n_random = argc > 2 ? strtoull(argv[2], nullptr, 10) : 2000000;
    std::vector<double> x(n_random + POW_INPUTS_FIXED), y(x.size());
    const size_t n = make_inputs(x.data(), y.data(), n_random);
    x.resize(n);
    y.resize(n);
    std::vector<double> z(n);
    hwy_pow_array(x.data(), y.data(), z.data(), n);

    const std::string pre = argv[1];
    FILE *f = fopen((pre + ".in").c_str(), "wb");
    fwrite(x.data(), 8, n, f);
    fwrite(y.data(), 8, n, f);
    fclose(f);
    std::vector<double> zc(z);  // canonical NaNs: the default NaN's sign differs between x86-64 and AArch64
    for (double &v : zc) if (std::isnan(v)) v = bitsd(0x7ff8000000000000ULL);
    f = fopen((pre + ".out").c_str(), "wb");
    fwrite(zc.data(), 8, n, f);
    fclose(f);

    acc_stats s = ACC_STATS_INIT;
    for (size_t i = 0; i < n; i++) {
        const bool snan = (std::isnan(x[i]) && !((dbits(x[i]) >> 51) & 1)) || (std::isnan(y[i]) && !((dbits(y[i]) >> 51) & 1));
        if (!snan) acc_add(&s, x[i], z[i], powl((long double)x[i], (long double)y[i]));
    }
    printf("%-10s %2zu lanes  %-9s %zu pairs: max %.4f ulp, not correctly rounded %.4f%%, > 1 ulp: %zu, special "
           "mismatches: %zu (long double: %d-bit mantissa)\n",
           hwy::TargetName(HWY_STATIC_TARGET), hwy_pow_lanes(), hwy_pow_variant() ? "SVML" : "std::pow", n, s.max_err,
           100.0 * s.not_cr / s.n, s.over1, s.mismatched_special, LDBL_MANT_DIG);
    return 0;
}
