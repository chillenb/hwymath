// Highway tan on every compiled and supported target: accuracy, bit-identity, comparisons (hwy_test_tmpl.h).
#include <cmath>
extern "C" {
#include "functions/tan/tan_d.h"
}
#define TRIG_NAME "tan"
#define TRIG_HWY_ARRAY hwy_tan_array
#define TRIG_HWY_VARIANT hwy_tan_variant
#define TRIG_HWY_LANES hwy_tan_lanes
#define TRIG_REFL tanl
#define TRIG_STD std::tan
#define TRIG_THRESHOLD 0x1.000000e4db24cp+16
#define TRIG_SVML_ARRAY svml_tan_ha_array

// extra check per SVML-algorithm target: the vrcp14pd emulations vs the reference algorithm (the
// complete one on all inputs, the one for the main path on inputs with biased exponent in [1, 0x7fc])
#include <cstdio>
#include <vector>
#include "functions/tan/rcp14_ref.h"
#include "functions/tan/tests/rcp14_inputs.h"
extern "C" void hwy_tan_rcp14_array(const double *x, double *y, size_t len);
extern "C" void hwy_tan_rcp14_normal_array(const double *x, double *y, size_t len);
static int rcp14_check(const char *target) {
    static std::vector<uint64_t> in;
    if (in.empty()) {
        const uint64_t saved = rng_state;  // keep the main inputs' random sequence independent of this
        in.resize(rcp14_max_inputs(1000000));
        in.resize(rcp14_inputs(in.data(), 1000000));
        rng_state = saved;
    }
    std::vector<uint64_t> out(in.size()), out_n(in.size());
    hwy_tan_rcp14_array(reinterpret_cast<const double *>(in.data()), reinterpret_cast<double *>(out.data()), in.size());
    hwy_tan_rcp14_normal_array(reinterpret_cast<const double *>(in.data()), reinterpret_cast<double *>(out_n.data()),
                               in.size());
    size_t bad = 0, n_normal = 0, bad_normal = 0;
    for (size_t i = 0; i < in.size(); i++) {
        const uint64_t ref = rcp14_ref_bits(in[i]), e = (in[i] >> 52) & 0x7ff;
        if (out[i] != ref && bad++ < 5) {
            printf("  rcp14 %s: x = %016llx: %016llx, reference %016llx\n", target, (unsigned long long)in[i],
                   (unsigned long long)out[i], (unsigned long long)ref);
        }
        if (e >= 1 && e <= 0x7fc) {
            n_normal++;
            if (out_n[i] != ref && bad_normal++ < 5) {
                printf("  rcp14 (normal range) %s: x = %016llx: %016llx, reference %016llx\n", target,
                       (unsigned long long)in[i], (unsigned long long)out_n[i], (unsigned long long)ref);
            }
        }
    }
    printf("  %s: vrcp14pd emulation vs reference: %zu inputs, %zu mismatches; normal-range version: %zu inputs, "
           "%zu mismatches\n", target, in.size(), bad, n_normal, bad_normal);
    return bad != 0 || bad_normal != 0;
}
#define TRIG_TARGET_CHECK(name) rcp14_check(name)
#include "functions/sin/tests/hwy_test_tmpl.h"
