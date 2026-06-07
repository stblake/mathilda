# SetAccuracy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SetAccuracy[x, n]
    Returns an expression equivalent to x with numeric values
    re-rounded or promoted to n decimal digits of accuracy.
    Requires a USE_MPFR build for high-accuracy outputs.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_set_accuracy` re-expresses a value to a target *accuracy* (digits after the decimal point) by converting accuracy to precision. It extracts the numeric accuracy `n` (integer/real/rational, or `MachinePrecision` which short-circuits to a machine-spec `numericalize`), then computes the required significant digits as `digits = n + log10(|x|)` using `expr_log10_abs`, floored at 1. It builds a `NumericSpec` (MPFR bits via `numeric_digits_to_bits(digits)`, or machine spec without MPFR) and calls `numericalize`. This is the standard "accuracy = digits past the point" approximation, not Mathematica's full significance-arithmetic semantics. Non-positive accuracy or unrecognised argument types return `NULL`.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/precision.c`](https://github.com/stblake/mathilda/blob/main/src/precision.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
