# Set

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
lhs = rhs or Set[lhs, rhs]
    evaluates rhs once and assigns the result to lhs.  When lhs is a
    symbol, the assignment is stored as an OwnValue; when lhs has the
    form f[args...] it is stored as a DownValue on f.  Set has attribute
    HoldFirst so lhs is not evaluated before the assignment.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `HoldFirst`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
