# Floor

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Floor[x]
    gives the greatest integer less than or equal to x.
Floor[x, a]
    gives the greatest multiple of a less than or equal to x.
Floor is Listable. Exact (Integer / BigInt / Rational) inputs return
exact integers; Real / MPFR inputs are rounded toward -Infinity at
the input precision; symbolic inputs stay unevaluated.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
