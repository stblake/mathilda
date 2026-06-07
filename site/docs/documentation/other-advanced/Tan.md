# Tan

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Tan[z]
    gives the tangent of z. Equivalent to Sin[z] / Cos[z].
Tan is Listable. Singularities at z = Pi/2 + k Pi yield ComplexInfinity.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_tan` mirrors the `Sin` cascade in `src/trig.c`: `strip_inverse_call(arg, "ArcTan")` for `Tan[ArcTan[x]] -> x`; `try_simp_forward_of_inverse` for `Tan` of the other inverse trig functions (`Tan[ArcSin[x]] -> x/Sqrt[1-x^2]`, `Tan[ArcCos[x]] -> Sqrt[1-x^2]/x`, `Tan[ArcCot[x]] -> 1/x`); `odd_fold` for the odd symmetry `Tan[-x] -> -Tan[x]`; `trig_i_fold` for `Tan[I y] -> I Tanh[y]`; and `Tan[0] = 0`. Exact rational-multiple-of-`Pi` values are detected by `extract_pi_multiplier` and computed by `exact_tan` (a denominator-switch table analogous to `exact_sin`).

**Numeric.** MPFR-valued arguments go through `numeric_mpfr_apply_unary(..., mpfr_tan)` with an `mpfr_complex_tan` complex fallback; otherwise `get_approx` plus `ctan` produces a machine-precision real or `Complex` result, only when the argument is inexact. Symbolic arguments return `NULL`. Attributes: `ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/trig.c`](https://github.com/stblake/mathilda/blob/main/src/trig.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
