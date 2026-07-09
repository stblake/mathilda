# Scan

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Scan[f, expr]
    Applies f to each element of expr for its side effects
    and returns Null. Over an association, applies f to each value.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= s = 0; Scan[(s = s + #) &, <|"a" -> 1, "b" -> 2, "c" -> 3|>]; s
Out[1]= 6
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
