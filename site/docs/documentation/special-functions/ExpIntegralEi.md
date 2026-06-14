# ExpIntegralEi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ExpIntegralEi[z]
    gives the exponential integral Ei(z), the principal value of
    -Integral_{-z}^Infinity e^-t/t dt, with a branch cut on (-Infinity, 0).
ExpIntegralEi[0] = -Infinity, ExpIntegralEi[Infinity] = Infinity,
ExpIntegralEi[-Infinity] = 0, ExpIntegralEi[+-I Infinity] = +-I Pi. Real and
complex inputs evaluate numerically at machine or arbitrary (MPFR) precision;
D[ExpIntegralEi[z], z] = E^z/z. Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Exact special values: `ExpIntegralEi[0] = -Infinity`,

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
