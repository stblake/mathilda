# AnyTrue

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AnyTrue[list, test]`**

Gives True if test\[e\] is True for some element e (False for an empty list). Over an association, tests the values.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= AllTrue[{2, 4, 6}, EvenQ]
Out[1]= True

In[2]:= AnyTrue[{1, 3, 4}, EvenQ]
Out[2]= True

In[3]:= NoneTrue[{1, 3, 5}, EvenQ]
Out[3]= True
```

## Implementation notes

**Attributes:** `Protected`.

## See also

[AllTrue](../../functional-programming/AllTrue/), [NoneTrue](../../functional-programming/NoneTrue/)

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_ndarray.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray.c)
- Tests: [`tests/test_ndarray_selection.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_selection.c)
