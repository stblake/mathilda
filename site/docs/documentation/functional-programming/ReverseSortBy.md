# ReverseSortBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ReverseSortBy[list, f]`**

Sorts by f in descending order. Over an association, sorts by f of each value, descending.

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

## References

**See also:** [ReverseSort](../../functional-programming/ReverseSort/), [Sort](../../data-structures/Sort/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
