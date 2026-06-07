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

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
