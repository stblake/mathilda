# PreIncrement

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PreIncrement[x] or ++x
    increases the value of x by 1, returning the new value of x.
    ++x is equivalent to x = x + 1.

PreIncrement has attribute HoldFirst. In PreIncrement[x], x can be a
symbol or a Part expression referring to an existing value.
If x has no assigned value, PreIncrement::rvalue is emitted and the
expression is left unevaluated.
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
