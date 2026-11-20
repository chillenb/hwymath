/* Machine model used by the scalar callouts that functions/asin/tools/asm2c.py translates mechanically
 * from the SVML assembly (see the description there): memory accesses, double <-> bits, and the
 * evaluation of x86 condition codes from the operands of the last flag-setting instruction.
 * Portable C99 (+ C11 _Alignas); assumes little-endian byte order. */
#ifndef HWYMATH_ASM_MODEL_H
#define HWYMATH_ASM_MODEL_H

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#if defined(__GNUC__)
#define HWYMATH_HIDDEN __attribute__((visibility("hidden")))
#else
#define HWYMATH_HIDDEN
#endif

/* reads of a volatile zero: "pxor %xmm0, %xmm0" must not let the compiler fold e.g. 0.0/0.0 */
static volatile double asm_model_zero = 0.0;
#define fzero (asm_model_zero)

static inline double am_bitsd(uint64_t u) { double d; memcpy(&d, &u, sizeof d); return d; }
static inline uint64_t am_dbits(double d) { uint64_t u; memcpy(&u, &d, sizeof u); return u; }
#define bitsd am_bitsd
#define dbits am_dbits

static inline uint8_t ld8(uint64_t a) { uint8_t v; memcpy(&v, (const void *)(uintptr_t)a, 1); return v; }
static inline uint16_t ld16(uint64_t a) { uint16_t v; memcpy(&v, (const void *)(uintptr_t)a, 2); return v; }
static inline uint32_t ld32(uint64_t a) { uint32_t v; memcpy(&v, (const void *)(uintptr_t)a, 4); return v; }
static inline uint64_t ld64(uint64_t a) { uint64_t v; memcpy(&v, (const void *)(uintptr_t)a, 8); return v; }
static inline double ldd(uint64_t a) { double v; memcpy(&v, (const void *)(uintptr_t)a, 8); return v; }
static inline void st8(uint64_t a, uint8_t v) { memcpy((void *)(uintptr_t)a, &v, 1); }
static inline void st16(uint64_t a, uint16_t v) { memcpy((void *)(uintptr_t)a, &v, 2); }
static inline void st32(uint64_t a, uint32_t v) { memcpy((void *)(uintptr_t)a, &v, 4); }
static inline void st64(uint64_t a, uint64_t v) { memcpy((void *)(uintptr_t)a, &v, 8); }
static inline void std_(uint64_t a, double v) { memcpy((void *)(uintptr_t)a, &v, 8); }

/* Scalar SSE2 arithmetic (addsd, subsd, mulsd, divsd: IEEE binary64, round to nearest). Elsewhere the
 * NaN results are made those of x86 SSE: the first NaN operand (the destination register), quieted,
 * else the second, else (invalid operation) the default NaN, which is negative on x86 (positive on
 * AArch64). This keeps the callouts' results bit-identical across architectures. */
#if defined(__x86_64__) || defined(__i386__)
static inline double am_add(double a, double b) { return a + b; }
static inline double am_sub(double a, double b) { return a - b; }
static inline double am_mul(double a, double b) { return a * b; }
static inline double am_div(double a, double b) { return a / b; }
#else
static inline double am_nan(double a, double b) {
    if (isnan(a)) return am_bitsd(am_dbits(a) | 0x0008000000000000ULL);
    if (isnan(b)) return am_bitsd(am_dbits(b) | 0x0008000000000000ULL);
    return am_bitsd(0xfff8000000000000ULL);
}
static inline double am_add(double a, double b) { const double r = a + b; return isnan(r) ? am_nan(a, b) : r; }
static inline double am_sub(double a, double b) { const double r = a - b; return isnan(r) ? am_nan(a, b) : r; }
static inline double am_mul(double a, double b) { const double r = a * b; return isnan(r) ? am_nan(a, b) : r; }
static inline double am_div(double a, double b) { const double r = a / b; return isnan(r) ? am_nan(a, b) : r; }
#endif

/* kind of the last flag-setting instruction */
enum { FK_BAD = 0, FK_INT = 1, FK_COMI = 2, FK_UCOMI = 3 };

static inline uint64_t am_u(int w, uint64_t v) { return w == 64 ? v : v & ((1ULL << w) - 1); }
static inline int64_t am_s(int w, uint64_t v) {
    return w == 64 ? (int64_t)v : w == 32 ? (int64_t)(int32_t)v : w == 16 ? (int64_t)(int16_t)v : (int64_t)(int8_t)v;
}

/* comisd: signaling comparisons (invalid on any NaN), ucomisd: quiet ones (invalid only on sNaN), as in C */
#define AM_CC(name, int_expr, comi_expr, ucomi_expr)                                                 \
    static inline int cc_##name(int fk, int fw, uint64_t fa, uint64_t fb, double a, double b) {     \
        (void)fw; (void)fa; (void)fb; (void)a; (void)b;                                              \
        switch (fk) {                                                                                \
        case FK_INT: return (int_expr);                                                              \
        case FK_COMI: return (comi_expr);                                                            \
        case FK_UCOMI: return (ucomi_expr);                                                          \
        default: abort();                                                                            \
        }                                                                                            \
    }
AM_CC(e, am_u(fw, fa) == am_u(fw, fb), !(a < b) && !(a > b), !islessgreater(a, b))
AM_CC(ne, am_u(fw, fa) != am_u(fw, fb), (a < b) || (a > b), islessgreater(a, b))
AM_CC(b, am_u(fw, fa) < am_u(fw, fb), !(a >= b), !isgreaterequal(a, b))
AM_CC(be, am_u(fw, fa) <= am_u(fw, fb), !(a > b), !isgreater(a, b))
AM_CC(a, am_u(fw, fa) > am_u(fw, fb), a > b, isgreater(a, b))
AM_CC(ae, am_u(fw, fa) >= am_u(fw, fb), a >= b, isgreaterequal(a, b))
AM_CC(l, am_s(fw, fa) < am_s(fw, fb), (abort(), 0), (abort(), 0))
AM_CC(le, am_s(fw, fa) <= am_s(fw, fb), (abort(), 0), (abort(), 0))
AM_CC(g, am_s(fw, fa) > am_s(fw, fb), (abort(), 0), (abort(), 0))
AM_CC(ge, am_s(fw, fa) >= am_s(fw, fb), (abort(), 0), (abort(), 0))
AM_CC(p, (abort(), 0), !(a <= b) && !(a > b), isunordered(a, b))
AM_CC(np, (abort(), 0), (a <= b) || (a > b), !isunordered(a, b))
#undef AM_CC

/* cmovl: signed less (integer flags only) */
static inline int cc_l4(int fk, int fw, uint64_t fa, uint64_t fb) {
    if (fk != FK_INT) abort();
    return am_s(fw, fa) < am_s(fw, fb);
}

#endif
