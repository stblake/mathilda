# Interpolation

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Interpolation[data]
    constructs an InterpolatingFunction that interpolates data, given
    as {f1, f2, ...} (values at x = 1, 2, ...), {{x1, f1}, ...} (values
    at given abscissae), or {{{x1, y1, ...}, f1}, ...} (an m-D tensor
    grid).
Interpolation[data, x]
    builds the interpolating function and evaluates it at x (a number,
    or a coordinate list in m-D).
Interpolation[{{{x1,...}, f1, df1, ddf1, ...}, ...}]
    reproduces supplied derivatives at the nodes (df = gradient, ddf =
    Hessian, ...) by tensor-product Hermite interpolation.
Interpolation[data, InterpolationOrder -> n]
    uses piecewise-polynomial pieces of degree n (default 3; 0 gives a
    piecewise-constant and 1 a piecewise-linear interpolant).
Interpolation[data, Method -> m]
    selects "Spline" (natural/cyclic cubic spline) or "Hermite"
    (piecewise cubic Hermite with estimated slopes).
Interpolation[data, PeriodicInterpolation -> True]
    builds a periodic interpolant (period = the data span; the data must
    repeat its first sample at the last). A per-dimension {True, False}
    list selects periodicity per axis.
Vector- or array-valued samples (f_i a list) are interpolated
component-wise and return an array of the same shape.
Works at machine or arbitrary (MPFR) precision, matching the data.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
