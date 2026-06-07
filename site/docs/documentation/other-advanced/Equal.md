# Equal

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
lhs == rhs or Equal[lhs, rhs]
    tests mathematical equality. Numeric arguments decide directly
    (Integer / Rational exact comparison; Real / MPFR comparison with
    precision tolerance); structurally identical symbolic forms decide
    True; otherwise the call stays unevaluated as a symbolic equation.
Equal threads over Lists pairwise; chained Equal becomes Inequality.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
