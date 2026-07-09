# DeleteDuplicatesBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
DeleteDuplicatesBy[expr, f]
    Keeps the first element for each distinct
    f[element], preserving order. Over an association, f is applied to the
    values and the surviving entries are kept (keys preserved).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= DeleteDuplicatesBy[{1, 12, 3, 14, 5}, EvenQ]
Out[1]= {1, 12}

In[2]:= DeleteDuplicatesBy[<|"a" -> 1, "b" -> 12, "c" -> 3|>, EvenQ]
Out[2]= <|"a" -> 1, "b" -> 12|>
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list/list_init.c`](https://github.com/stblake/mathilda/blob/main/src/list/list_init.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
