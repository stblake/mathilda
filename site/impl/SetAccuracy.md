---
source: src/precision.c
---
`builtin_set_accuracy` re-expresses a value to a target *accuracy* (digits after the decimal point) by converting accuracy to precision. It extracts the numeric accuracy `n` (integer/real/rational, or `MachinePrecision` which short-circuits to a machine-spec `numericalize`), then computes the required significant digits as `digits = n + log10(|x|)` using `expr_log10_abs`, floored at 1. It builds a `NumericSpec` (MPFR bits via `numeric_digits_to_bits(digits)`, or machine spec without MPFR) and calls `numericalize`. This is the standard "accuracy = digits past the point" approximation, not full significance-arithmetic semantics. Non-positive accuracy or unrecognised argument types return `NULL`.
