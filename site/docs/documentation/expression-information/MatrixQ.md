# MatrixQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
MatrixQ[expr]
    gives True if expr is a list of lists that can represent a matrix, and gives False otherwise.
MatrixQ[expr, test]
    gives True only if test yields True when applied to each of the matrix elements in expr.

MatrixQ[expr] gives True only if expr is a list and each of its elements is a list of the same length,
containing no elements that are themselves lists.
MatrixQ[expr, NumberQ] tests whether expr is a numerical matrix.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_matrixq` (`src/list.c`) returns `True` when the argument is a non-empty `List` of equal-length `List` rows whose entries are non-lists; with an optional `test` argument each entry must instead satisfy `test[entry]` (evaluated to `True`). Any ragged row, empty outer list, or failing entry yields `False`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= MatrixQ[{{1, 2}, {3, 4}}]
Out[1]= True

In[2]:= MatrixQ[{1, 2, 3}]
Out[2]= False
```

### Notes

`MatrixQ` is `True` for a list of equal-length lists; a flat (rank-1) list is a vector, not a matrix.
