# LogIntegral

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
LogIntegral[z]
    gives the logarithmic integral li(z), the principal value of
    Integral_0^z dt/ln t, equal to ExpIntegralEi[Log[z]], with a branch cut
    on (-Infinity, 1). LogIntegral[0] = 0, LogIntegral[1] = -Infinity,
LogIntegral[Infinity] = Infinity. Real and complex inputs evaluate
numerically at machine or arbitrary (MPFR) precision; D[LogIntegral[z], z] =
1/Log[z]. Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Exact special values: `LogIntegral[0] = 0`, `LogIntegral[1] = -Infinity`,

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
