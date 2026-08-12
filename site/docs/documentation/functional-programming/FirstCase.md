# FirstCase

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FirstCase[expr, patt]`**

Gives the first element of expr matching patt, or Missing\["NotFound"\]. FirstCase\[expr, patt, default\] uses default. Over an association, matches values and returns the first match.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= SelectFirst[{1, 3, 4, 5, 6}, EvenQ]
Out[1]= 4

In[2]:= FirstCase[{1, 2, 3, 4}, _?EvenQ]
Out[2]= 2

In[3]:= SelectFirst[{1, 3, 5}, EvenQ, None]
Out[3]= None
```

## Implementation notes

**Attributes:** `Protected`.

## See also

[SelectFirst](../../functional-programming/SelectFirst/)

## References

- Source: [`src/patterns.c`](https://github.com/stblake/mathilda/blob/main/src/patterns.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
