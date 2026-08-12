# ReverseSort

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ReverseSort[list]`**

Sorts into descending order (Reverse of Sort). Over an association, sorts the entries by value, descending.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= ReverseSort[{3, 1, 4, 1, 5, 9, 2}]
Out[1]= {9, 5, 4, 3, 2, 1, 1}

In[2]:= ReverseSort[<|"a" -> 3, "b" -> 1, "c" -> 2|>]
Out[2]= <|"a" -> 3, "c" -> 2, "b" -> 1|>
```

## Implementation notes

**Attributes:** `Protected`.

## See also

[ReverseSortBy](../../functional-programming/ReverseSortBy/), [Sort](../../data-structures/Sort/)

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_linalg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_linalg.c)
