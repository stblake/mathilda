# TakeLargestBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
TakeLargestBy[list, f, n]
    Gives the n elements of list for which f is
    largest, in descending order of f. Over an association, ranks by f of
    each value.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= TakeLargest[{3, 1, 4, 1, 5, 9, 2, 6}, 3]
Out[1]= {9, 6, 5}

In[2]:= TakeLargest[<|"a" -> 3, "b" -> 9, "c" -> 1, "d" -> 6|>, 2]
Out[2]= <|"b" -> 9, "d" -> 6|>

In[3]:= TakeLargestBy[{-9, 2, -3, 5}, Abs, 2]
Out[3]= {-9, 5}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
