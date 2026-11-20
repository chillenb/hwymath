// Compile-time tables for the Highway versions of the log family (included once, outside the per-target code).
#ifndef HWYMATH_LOG_TABLES_H
#define HWYMATH_LOG_TABLES_H

#include <cstddef>
#include <cstdint>

#include "functions/log/hwy/rcp14_thresholds.h"

namespace hwymath_log {

// 16-entry table repeated to fill vectors of up to 64 doubles, so that TableLookupLanes can use one vector.
// Must be built at compile time: function-local statics with run-time initialization put thread-safe init
// guards into the hot loop.
struct alignas(64) Padded16 {
    double v[64];
    explicit constexpr Padded16(const double *t) : v{} {
        for (size_t i = 0; i < 64; i++) v[i] = t[i % 16];
    }
};

struct alignas(64) Padded32U {
    uint32_t v[64];
    explicit constexpr Padded32U(const uint32_t *t) : v{} {
        for (size_t i = 0; i < 64; i++) v[i] = t[i % 32];
    }
};

struct Tables {
    static const uint32_t *BinsPad() { static constexpr Padded32U p(kRcpBins.v); return p.v; }
};

}  // namespace hwymath_log

#endif
