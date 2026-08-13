# TakeLargest

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`TakeLargest[list, n]`**

Gives the n largest elements of list, in descending order. Over an association, gives the n entries with the largest values (as an association).

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

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

## References

**See also:** [TakeSmallest](../../functional-programming/TakeSmallest/), [TakeLargestBy](../../functional-programming/TakeLargestBy/), [TakeSmallestBy](../../functional-programming/TakeSmallestBy/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
