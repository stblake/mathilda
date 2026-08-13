# SelectFirst

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SelectFirst[list, pred]`**

Gives the first element e of list for which pred\[e\] is True, or Missing\["NotFound"\]. SelectFirst\[list, pred, default\] uses default. Over an association, tests values and returns the first match.

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

## References

**See also:** [FirstCase](../../functional-programming/FirstCase/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_ndarray_selection.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_selection.c)
- Tests: [`tests/test_pred_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_pred_compile.c)
