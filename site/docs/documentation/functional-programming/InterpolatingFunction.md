# InterpolatingFunction

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
InterpolatingFunction[domain, table]
    represents an approximate function whose values are found by
    interpolation. domain is {{x1min, x1max}, ...} with one interval
    per dimension; table is a list of {coord, value} data points on a
    full tensor grid (coord is a scalar for 1-D, an {x1, ..., xm} list
    for m-D).
InterpolatingFunction[...][x1, ..., xm]
    gives the interpolated value using tensor-product piecewise-
    polynomial (default order 3) interpolation. Arguments outside the
    domain are extrapolated with a warning.
Derivative[d1, ..., dm][InterpolatingFunction[...]]
    gives an InterpolatingFunction for the mixed partial derivative.
In standard output only the domain is shown; the rest is <>.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
