# VectorQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VectorQ[expr]`**

gives True if expr is a list, none of whose elements are themselves lists, and gives False otherwise.

**`VectorQ[expr, test]`**

gives True only if test yields True when applied to each of the elements in expr.

**`VectorQ[expr, NumberQ] tests whether expr is a vector of numbers.`**

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

```mathematica
In[1]:= VectorQ[{1, 2, 3}]
Out[1]= True

In[2]:= VectorQ[{{1, 2}, {3, 4}}]
Out[2]= False
```

## Implementation notes

`builtin_vectorq` (`src/list.c`) returns `True` when the argument is a `List` (`is_listq`) none of whose elements is itself a `List`. With a second argument `test`, each element is instead required to satisfy `test[elem]` (evaluated, must yield `True`). Returns `False` as soon as a check fails.

**Attributes:** `Protected`.

## References

**See also:** [ListQ](../../expression-information/ListQ/), [MatrixQ](../../expression-information/MatrixQ/), [List](../../other-advanced/List/)

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
- Tests: [`tests/test_ml_pca.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_pca.c)
- Tests: [`tests/test_ndarray.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)

## Notes & additional examples

### Notes

`VectorQ` is `True` for a flat list none of whose elements are themselves lists; a list of lists (a matrix) is not a vector.
