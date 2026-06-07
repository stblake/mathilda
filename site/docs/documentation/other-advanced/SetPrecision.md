# SetPrecision

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SetPrecision[x, n]
    Returns an expression equivalent to x with numeric values
    re-rounded or promoted to n decimal digits of precision.
    Requires a USE_MPFR build for n > MachinePrecision.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_set_precision` is a two-argument wrapper: it parses the precision argument into a `NumericSpec` via `parse_prec_arg` (accepting an integer/real digit count or `MachinePrecision`) and drives `numericalize(value, spec)` — the same engine `N` uses — to re-represent the value to the requested number of significant digits (an `EXPR_MPFR` at the corresponding bit width when MPFR is built, otherwise machine `double`). Returns `NULL` if the precision argument is not a valid spec.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/precision.c`](https://github.com/stblake/mathilda/blob/main/src/precision.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
