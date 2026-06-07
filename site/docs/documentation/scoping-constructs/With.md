# With

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
With[{x = x0, ...}, expr] specifies that x should be replaced by x0 throughout expr.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= x = 10; With[{x = 5}, x^2]
Out[1]= 25
```

## Implementation notes

- `HoldAll`, `Protected`.
- Replaces occurrences of symbols in the body before evaluation.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/scoping-constructs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/scoping-constructs.md)
