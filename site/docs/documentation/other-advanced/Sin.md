# Sin

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Sin[z]
    gives the sine of z (argument in radians).
Sin is Listable. Numeric inputs are evaluated via libm (Real) or MPFR
(arbitrary precision); rational multiples of Pi reduce to exact values.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_sin` is a single-argument cascade tried in order: (1) `strip_inverse_call(arg, "ArcSin")` collapses `Sin[ArcSin[x]] -> x`; (2) `try_simp_forward_of_inverse` rewrites `Sin` of the *other* inverse trig functions to radical forms (`Sin[ArcCos[x]] -> Sqrt[1-x^2]`, `Sin[ArcTan[x]] -> x/Sqrt[1+x^2]`); (3) `odd_fold` uses the odd symmetry `Sin[-x] -> -Sin[x]` whenever `expr_is_superficially_negative(arg)`; (4) `trig_i_fold` extracts an imaginary unit, `Sin[I y] -> I Sinh[y]`; (5) `Sin[0] = 0`. Exact special values come from `extract_pi_multiplier`, which detects `Pi` or `Times[Rational[n,d], Pi]`, handing `(n,d)` to `exact_sin`.

`exact_sin` reduces the angle into `[0, Pi/2]` (tracking a sign through the `[0,2Pi)` and `[0,Pi]` foldings), reduces the fraction by `gcd`, and switches on the denominator: closed radical forms are tabulated for `d` in {1,2,3,4,5,6,10,12}, e.g. `d==3 -> Sqrt[3]/2`, `d==5`/`d==10` give the golden-ratio nested-radical values, `d==12` gives `(Sqrt[6]±Sqrt[2])/4`. Anything outside the table returns `NULL` and stays symbolic.

**Numeric.** With MPFR built, an MPFR-valued argument is evaluated via `numeric_mpfr_apply_unary(..., mpfr_sin)`, falling back to `mpfr_complex_sin` for complex MPFR values, preserving the input precision. Otherwise `get_approx` extracts a `double complex` and, only if the value is genuinely inexact, applies `csin`; a real input yields `EXPR_REAL`, a complex one a `Complex[Real, Real]`.

`Sin` is registered with `ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED`; element-wise threading over lists is handled generically by the evaluator.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/trig.c`](https://github.com/stblake/mathilda/blob/main/src/trig.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
