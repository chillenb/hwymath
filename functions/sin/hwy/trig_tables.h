// Compile-time tables for the Highway ports of SVML's sin, cos and tan (functions/{sin,cos,tan}/hwy):
// the large-argument reduction data and padded copies of small tables. Included once per translation
// unit (outside the per-target code).
// Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.
#ifndef HWYMATH_TRIG_TABLES_H
#define HWYMATH_TRIG_TABLES_H
#include <cstddef>
#include <cstdint>

#include "functions/sin/svml/svml_trig_reduction_data.inc"

namespace hwymath_trig {
// A table of kT entries repeated to 64 entries, built at compile time (padded tables with run-time
// initialization would put initialization guards into the hot loop).
template <typename T, size_t kT>
struct alignas(64) Padded64 {
    T v[64];
    explicit constexpr Padded64(const T *t) : v{} {
        for (size_t i = 0; i < 64; i++) v[i] = t[i % kT];
    }
};
}  // namespace hwymath_trig
#endif
