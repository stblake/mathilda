# SubtractFrom

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SubtractFrom[x, dx] or x -= dx
    subtracts dx from x and returns the new value of x.
    x -= dx is equivalent to x = x - dx.

SubtractFrom has attribute HoldFirst. The first argument x can be a
symbol or a Part expression referring to an existing value; dx may be a
number, a symbolic expression, or a list. If x has no assigned value,
SubtractFrom::rvalue is emitted and the expression is left unevaluated.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `HoldFirst`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
