/* Checks the AVX2 port against a scalar model of the AVX-512 instruction semantics of
 * __svml_exp8_ha, built with different mechanisms than the port's AVX2 replacements:
 *   {rz-sae} FMA  -> fesetround(FE_TOWARDZERO) + fma()        (port: RN fma + error-sign correction)
 *   vpermt2pd     -> array indexing with the low 4 bits of t (port: gathers)
 *   vscalefpd     -> ldexp(y, floor(N)), NaN N propagates      (port: exponent bits from t)
 * Must be compiled with -frounding-math. Runs on any AVX2 machine. */
#include <fenv.h>
#include <stdio.h>

#include "../svml_exp_d_ha.h"
#include "common.h"
#include "../constants.inc"

static double
model_lane(double x)
{
    fesetround(FE_TOWARDZERO);
    const double t = fma(x, vLog2e, vShifter);
    fesetround(FE_TONEAREST);
    const double ax = bitsd(dbits(x) & vAbsMaskBits);
    const double n = t - vShifter;
    const unsigned j = (unsigned)(dbits(t) & 15);
    const double tlr = vTlr[j], th = vTh[j];
    double r = fma(-n, vLn2Hi, x);
    const int callout = ax >= vThreshold;  /* ordered compare: false for NaN */
    r = fma(-vLn2Lo, n, r);
    r = bitsd(dbits(r) & vRMaskBits);
    const double r2 = r * r;
    const double p56 = fma(vC6, r, vC5), p34 = fma(vC4, r, vC3), p12 = fma(vC2, r, vC1);
    const double p36 = fma(p56, r2, p34);
    const double p16 = fma(p36, r2, p12);
    const double p = fma(p16, r, tlr);
    const double y = fma(th, p, th);
    double res;
    if (isnan(n)) res = n;            /* vscalefpd returns the (quiet) NaN of the scale operand */
    else if (isnan(y)) res = y;
    else res = ldexp(y, (int)fmax(fmin(floor(n), 4096.0), -4096.0));
    if (callout) svml_dexp_ha_cout_rare(&x, &res);
    return res;
}

__attribute__((target("avx2,fma"))) static void
port4(const double *x, double *y)
{
    _mm256_storeu_pd(y, svml_exp4_ha_avx2(_mm256_loadu_pd(x)));
}

int
main(int argc, char **argv)
{
    size_t n_random = argc > 1 ? strtoull(argv[1], NULL, 10) : 20000000;
    double *in = malloc((n_random + 4096 + 4) * sizeof *in);
    size_t n = make_inputs(in, n_random);
    while (n % 4) in[n++] = 0.0;
    size_t bad = 0;
    for (size_t i = 0; i < n; i += 4) {
        double y[4];
        port4(in + i, y);
        for (int l = 0; l < 4; l++) {
            double m = model_lane(in[i + l]);
            if (dbits(m) != dbits(y[l]) && bad++ < 10) {
                printf("MISMATCH x=%a: model %a  AVX2 port %a\n", in[i + l], m, y[l]);
            }
        }
    }
    printf("test_model: %zu inputs, %zu mismatches between the AVX2 port and the AVX-512 semantics model\n", n, bad);
    free(in);
    return bad != 0;
}
