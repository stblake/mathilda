# HornerForm

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
HornerForm[poly]
    rewrites the univariate polynomial poly in nested (Horner) form,
    which evaluates in n multiplications and n additions instead of
    the naive 2n.
HornerForm[poly, var]
    uses var as the recursion variable for multivariate poly.
HornerForm[poly1 / poly2, vars1, vars2]
    puts a rational function in Horner form, nested with respect to
    vars1 in the numerator and vars2 in the denominator.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_hornerform` (in `src/poly/poly.c`) rewrites a polynomial in nested multiply-add form. It first splits a rational input into numerator/denominator (scanning a top-level `Times`/`Power` for negative-exponent factors), then applies the recursive worker `horner_form_rec` to the numerator. That worker `expr_expand`s, verifies the expression is a polynomial in the lead variable via `PolynomialQ`, obtains its dense coefficients with `CoefficientList`, and folds from the highest-degree coefficient downward using the Horner recurrence `H ← c_i + v·H` (with zero-coefficient short-circuits). For multiple variables it recurses on each coefficient with the remaining variable list, producing a nested Horner form.

**Data structures.** Reuses the `CoefficientList` array as the coefficient sequence; nested `Times`/`Plus` nodes (built via `internal_times`/`internal_plus`) form the output.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- W. G. Horner, "A new method of solving numerical equations of all orders, by continuous approximation", Phil. Trans. R. Soc. 1819.
- Source: [`src/poly/poly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/poly.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
