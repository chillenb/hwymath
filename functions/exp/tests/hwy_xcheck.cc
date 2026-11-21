// Cross-platform check of the two Highway exps (hp: hwy_exp.cc; svml, Intel's algorithm: hwy_exp_svml.cc),
// built for one static target (-DHWY_COMPILE_ONLY_STATIC): writes the inputs and results to files
// (PREFIX.in, PREFIX.out for hp, PREFIX.svml.out for svml; compare them bit for bit between
// targets/architectures) and prints accuracy vs expl (binary128 long double on AArch64, x87 extended on x86-64).
// usage: hwy_xcheck OUTPUT_PREFIX [N_RANDOM]       Build with -ffp-contract=off (identical inputs).
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "hwy/targets.h"

extern "C" {
void hwy_exp_array(const double *x, double *y, size_t len);
int hwy_exp_variant(void);
size_t hwy_exp_lanes(void);
void hwy_exp_svml_array(const double *x, double *y, size_t len);
int hwy_exp_svml_variant(void);
}
#include "inputs.h"

static void report(const char *alg, const std::vector<double> &x, const std::vector<double> &y) {
    double max_err = 0, worst = 0;
    size_t n = 0, not_cr = 0;
    for (size_t i = 0; i < x.size(); i++) {
        long double ref = expl((long double)x[i]);
        double rd = (double)ref;
        if (!std::isfinite(x[i]) || !std::isfinite(rd) || rd == 0.0) continue;
        long double ulp = (long double)std::nextafter(std::fabs(rd), INFINITY) - fabsl((long double)rd);
        if (std::fabs(rd) == DBL_MAX) ulp = 0x1p971L;
        double e = (double)(fabsl((long double)y[i] - ref) / ulp);
        n++;
        not_cr += e > 0.5;
        if (e > max_err) { max_err = e; worst = x[i]; }
    }
    printf("%-10s %2zu lanes  %-14s %zu inputs: max %.4f ulp (x = %.17g), not correctly rounded %.4f%%  (long double: %d-bit mantissa)\n",
           hwy::TargetName(HWY_STATIC_TARGET), hwy_exp_lanes(), alg, x.size(), max_err, worst, 100.0 * not_cr / n,
           LDBL_MANT_DIG);
}

static void write(const std::string &path, const std::vector<double> &v) {
    FILE *f = fopen(path.c_str(), "wb");
    fwrite(v.data(), 8, v.size(), f);
    fclose(f);
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s OUTPUT_PREFIX [N_RANDOM]\n", argv[0]); return 2; }
    size_t n_random = argc > 2 ? strtoull(argv[2], nullptr, 10) : 2000000;
    std::vector<double> x(n_random + 4096);
    x.resize(make_inputs(x.data(), n_random));
    std::vector<double> y(x.size()), ys(x.size());
    hwy_exp_array(x.data(), y.data(), x.size());
    hwy_exp_svml_array(x.data(), ys.data(), x.size());

    std::string pre = argv[1];
    write(pre + ".in", x);
    write(pre + ".out", y);
    write(pre + ".svml.out", ys);

    report(hwy_exp_variant() ? (hwy_exp_variant() == 2 ? "hp, 4-entry" : "hp, 8-entry") : "std::exp", x, y);
    report(hwy_exp_svml_variant() ? "svml" : "std::exp", x, ys);
    return 0;
}
