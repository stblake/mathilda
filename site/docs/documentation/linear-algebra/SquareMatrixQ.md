# SquareMatrixQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SquareMatrixQ[m]
    gives True if m is a square matrix (Dimensions[m] == {n, n}),
    and False otherwise.

Works for symbolic as well as numerical matrices.  Returns False on
non-list, ragged, rectangular, empty, or higher-rank tensor inputs.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= SquareMatrixQ[{{1, 2}, {3, 4}}]
Out[1]= True

In[2]:= SquareMatrixQ[{{1, 2, 3}, {4, 5, 6}}]
Out[2]= False

In[3]:= SquareMatrixQ[{1, 2, 3}]
Out[3]= False

In[4]:= SquareMatrixQ[{{1}, {2, 3}}]
Out[4]= False

In[5]:= SquareMatrixQ[{{a, b, c}, {d, e, f}, {g, h, i}}]
Out[5]= True
```

## Implementation notes

`builtin_square_matrix_q` is a pure shape test: it returns `True` iff the argument is a non-empty `List` of equal-length `List`s with row count equal to column count and no entry that is itself a `List` (rejecting ragged, rectangular, and higher-rank tensors). No element predicate is consulted, so `{{x}}` is square for any `x`. Exactly one argument is accepted; any other count emits `SquareMatrixQ::argx` and leaves the call unevaluated.

- `Protected`.
- Pure shape test: no element predicate or option is consulted.
- Works for symbolic as well as numerical matrices (`{{a,b},{c,d}}` is

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
