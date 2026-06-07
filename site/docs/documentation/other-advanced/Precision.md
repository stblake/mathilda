# Precision

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Precision[x]
    Returns the number of decimal digits of precision in x.
    Exact numbers return Infinity; machine-precision reals return
    the symbol MachinePrecision; MPFR values return their declared
    precision in decimal digits.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_precision` (1-arg) delegates to the recursive `precision_of`. Exact quantities — `EXPR_INTEGER`, `EXPR_BIGINT`, exact `Rational`, exact numeric symbols — return `Infinity`. `EXPR_REAL` returns the symbol `MachinePrecision`. `EXPR_MPFR` returns its decimal precision computed from the stored bit precision (`mpfr_get_prec / log2(10)`) as an `EXPR_REAL`. For composite expressions, `Complex[re, im]` and general function heads take the minimum precision across parts via `precision_min` (which treats `MachinePrecision` as the constant `NUMERIC_MACHINE_PRECISION_DIGITS ≈ 15.95` when comparing against explicit MPFR digit counts). This mirrors Mathematica's rule that an expression is only as precise as its least-precise inexact part; the precision-tracking unit conversion (`LOG2_10`) is shared with `numeric.c`.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/precision.c`](https://github.com/stblake/mathilda/blob/main/src/precision.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
