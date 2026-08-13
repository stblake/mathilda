# SquareMatrixQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SquareMatrixQ[m]`**

gives True if m is a square matrix (Dimensions\[m\] == {n, n}), and False otherwise.

<details>
<summary>Notes</summary>

Works for symbolic as well as numerical matrices.  Returns False on non-list, ragged, rectangular, empty, or higher-rank tensor inputs.

</details>

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

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

### Applications (6)

```mathematica
In[6]:= SquareMatrixQ[{{1, 2}, {3, 4}}]
Out[6]= True

In[7]:= SquareMatrixQ[{{1, 2, 3}, {4, 5, 6}}]
Out[7]= False

In[8]:= SquareMatrixQ[{{a, b}, {c, d}}]
Out[8]= True

In[9]:= SquareMatrixQ[IdentityMatrix[5]]
Out[9]= True

In[10]:= SquareMatrixQ[{1, 2, 3}]
Out[10]= False

In[11]:= SquareMatrixQ[{{1, 2}, {3}}]
Out[11]= False
```

## Implementation notes

`builtin_square_matrix_q` is a pure shape test: it returns `True` iff the argument is a non-empty `List` of equal-length `List`s with row count equal to column count and no entry that is itself a `List` (rejecting ragged, rectangular, and higher-rank tensors). No element predicate is consulted, so `{{x}}` is square for any `x`. Exactly one argument is accepted; any other count emits `SquareMatrixQ::argx` and leaves the call unevaluated.

- `Protected`.
- Pure shape test: no element predicate or option is consulted.
- Works for symbolic as well as numerical matrices (`{{a,b},{c,d}}` is
  square; entries are not evaluated).
- Returns `False` (rather than leaving unevaluated) on non-list, scalar,
  vector, empty (`{}`, `{{}}`), ragged, rectangular, or higher-rank
  tensor inputs.
- Exactly one argument is accepted; any other count emits a
  Mathematica-compatible

  ```
  SquareMatrixQ::argx: SquareMatrixQ called with N arguments; 1 argument is expected.
  ```

  to `stderr` and leaves the call unevaluated.

**Attributes:** `Protected`.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_square_matrix_q.c`](https://github.com/stblake/mathilda/blob/main/tests/test_square_matrix_q.c)

## Notes & additional examples

### Notes

`SquareMatrixQ[m]` is `True` exactly when `m` is a rank-2 tensor whose two
dimensions are equal, i.e. `Dimensions[m] == {n, n}`. It returns `False` on
non-list, ragged, rectangular, empty, or higher-rank inputs, and it works for
symbolic as well as numerical entries.
