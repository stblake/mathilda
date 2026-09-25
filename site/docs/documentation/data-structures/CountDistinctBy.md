# CountDistinctBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CountDistinctBy[expr, f]`**

Gives the number of distinct values of f\[e\] over the elements e of expr (the values, for an association).

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= CountDistinctBy[{1, 2, 3, 4, 5}, EvenQ]
Out[1]= 2

In[2]:= CountDistinctBy[{"apple", "avocado", "banana"}, StringTake[#, 1] &]
Out[2]= 2
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
