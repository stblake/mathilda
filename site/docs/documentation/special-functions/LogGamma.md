# LogGamma

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
LogGamma[z]
    gives the log-gamma function log(Gamma(z)), analytic except for a branch
cut on the negative reals. Exact at integer and half-integer z (with the
negative-axis branch term), divergent (Infinity) at non-positive integers,
and evaluated numerically for real and complex z at machine or arbitrary
(MPFR) precision. D[LogGamma[z], z] is PolyGamma[0, z]. Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- **Exact closed forms.** Integers reduce as `LogGamma[n] = Log[(n-1)!]`

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
