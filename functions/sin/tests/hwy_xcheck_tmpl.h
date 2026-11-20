// Cross-platform check of a Highway trig function (sin, cos, tan) built for one static target
// (-DHWY_COMPILE_ONLY_STATIC): writes the inputs and results to files (compared bit for bit between
// targets and architectures by functions/sin/tools/xcheck.sh) and checks accuracy and special values
// against the long double reference (binary128 on AArch64, x87 extended on x86-64). Included by
// functions/{sin,cos,tan}/tests/hwy_xcheck.cc, which define TRIG_NAME, TRIG_HWY_ARRAY, TRIG_HWY_VARIANT,
// TRIG_HWY_LANES, TRIG_REFL. Exit status 1 if the error reaches 1 ulp or a special value (NaN, signed
// zero, e.g. f(-0)) differs from the reference.
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
void TRIG_HWY_ARRAY(const double *x, double *y, size_t len);
int TRIG_HWY_VARIANT(void);
size_t TRIG_HWY_LANES(void);
}
#include "functions/sin/tests/inputs.h"

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s OUTPUT_PREFIX [N_RANDOM]\n", argv[0]);
        return 2;
    }
    const size_t n_random = argc > 2 ? strtoull(argv[2], nullptr, 10) : 2000000;
    std::vector<double> x(n_random + 40000);
    x.resize(make_inputs(x.data(), n_random));
    std::vector<double> y(x.size());
    TRIG_HWY_ARRAY(x.data(), y.data(), x.size());

    const std::string pre = argv[1];
    FILE *f = fopen((pre + ".in").c_str(), "wb");
    fwrite(x.data(), 8, x.size(), f);
    fclose(f);
    f = fopen((pre + ".out").c_str(), "wb");
    fwrite(y.data(), 8, y.size(), f);
    fclose(f);

    struct acc_stats s = ACC_STATS_INIT;
    size_t shown = 0;
    for (size_t i = 0; i < x.size(); i++) {
        const size_t before = s.mismatched_special;
        const long double ref = TRIG_REFL((long double)x[i]);
        acc_add(&s, x[i], y[i], ref);
        if (s.mismatched_special != before && shown++ < 5) {
            printf("  special value: x = %a: %a, expected %a\n", x[i], y[i], (double)ref);
        }
    }
    char name[64];
    snprintf(name, sizeof name, "%s %s %zu lanes %s", TRIG_NAME, hwy::TargetName(HWY_STATIC_TARGET), TRIG_HWY_LANES(),
             TRIG_HWY_VARIANT() ? "SVML" : "std");
    acc_print(stdout, name, &s);
    const bool ok = s.max_err < 1.0 && s.mismatched_special == 0;
    if (!ok) printf("  FAIL\n");
    return ok ? 0 : 1;
}
