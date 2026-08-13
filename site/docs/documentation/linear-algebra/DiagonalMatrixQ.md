# DiagonalMatrixQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DiagonalMatrixQ[m]`**

gives True if m is diagonal, and False otherwise.

**`DiagonalMatrixQ[m, k]`**

gives True if m has nonzero elements only on the k-th diagonal, and False otherwise.  Positive k refers to superdiagonals above the main diagonal; negative k refers to subdiagonals below it. Works for rectangular as well as square matrices.

<details>
<summary>Notes</summary>

Option: Tolerance -\> Automatic   numeric tolerance for approximate matrices. With Tolerance -\> t, off-diagonal entries e are taken to be zero when Abs\[e\] \<= t evaluates to True.  Without a tolerance the test is structural: only literal numeric zeros (Integer 0, Real 0.0, BigInt 0) count as zero.  Returns False on non-matrix, ragged, empty (i.e. {}), or higher-rank tensor inputs; an n-by-0 matrix (e.g. {{}, {}}) is vacuously diagonal.

</details>

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= DiagonalMatrixQ[{{a, 0, 0}, {0, b, 0}, {0, 0, c}}]
Out[1]= True

In[2]:= DiagonalMatrixQ[{{1, 0, 0}, {0, 0, 2}, {3, 0, 0}}]
Out[2]= False

In[3]:= DiagonalMatrixQ[{{0, a, 0}, {0, 0, b}, {0, 0, 0}}, 1]
Out[3]= True

In[4]:= DiagonalMatrixQ[{{0, 0, 0}, {a, 0, 0}, {0, b, 0}}, -1]
Out[4]= True

In[5]:= DiagonalMatrixQ[{{1, 0, 0}, {0, 2, 0}}]
Out[5]= True

In[6]:= DiagonalMatrixQ[{{1, 0}, {0, 2}, {0, 0}}]
Out[6]= True

In[7]:= DiagonalMatrixQ[IdentityMatrix[5]]
Out[7]= True
```

### Options (1)

```mathematica
In[8]:= DiagonalMatrixQ[{{1., 10^-12, 0}, {0, 2., 10^-13}, {0, 0, 3.}}, Tolerance -> 10^-12]
Out[8]= True
```

### Applications (5)

```mathematica
In[9]:= DiagonalMatrixQ[{{1, 0}, {0, 2}}]
Out[9]= True

In[10]:= DiagonalMatrixQ[{{1, 2}, {0, 3}}]
Out[10]= False

In[11]:= DiagonalMatrixQ[DiagonalMatrix[{a, b, c}]]
Out[11]= True

In[12]:= DiagonalMatrixQ[{{0, 5, 0}, {0, 0, 7}, {0, 0, 0}}, 1]
Out[12]= True

In[13]:= DiagonalMatrixQ[{{0.0, 1.0*10^-15}, {0, 0.0}}, Tolerance -> 10^-10]
Out[13]= True
```

## Implementation notes

`builtin_diagonal_matrix_q` tests whether a matrix has nonzero entries only on the `k`-th diagonal (default `k = 0`). It accepts an optional integer `k` at position 2 and a `Tolerance` option; an empty or malformed argument list yields a `DiagonalMatrixQ::argt`/`::nonopt` diagnostic, and missing args return `False`. After validating that the matrix is a rectangular `List` of `List`s, it returns `True`/`False` according to whether every off-`k`-diagonal entry is structurally (or within `Tolerance`) zero.

- `Protected`.
- Works for rectangular matrices, not only square -- only the entry-zero
  predicate and the shape constraints matter.
- Default test is structural: only literal numeric zeros (`Integer 0`,
  `Real 0.0`, `BigInt 0`) count as zero.  Symbolic off-diagonal entries
  fail the test, so the predicate is conservative.
- Returns `False` (rather than leaving unevaluated) on non-list, scalar,
  vector, ragged, or higher-rank tensor inputs.  `{}` is rejected; an
  `n`-by-`0` matrix (e.g. `{{}, {}}`) is vacuously diagonal and returns
  `True`.
- Zero positional arguments emits a Mathematica-compatible

  ```
  DiagonalMatrixQ::argt: DiagonalMatrixQ called with 0 arguments; 1 or 2 arguments are expected.
  ```

  to `stderr` and leaves the call unevaluated.

- More than two positional arguments (or any non-`Rule` junk in the
  option region) emits

  ```
  DiagonalMatrixQ::nonopt: Options expected (instead of <expr>) beyond position 2 in DiagonalMatrixQ[...]. An option must be a rule or a list of rules.
  ```

  to `stderr` and leaves the call unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [Rule](../../assignment-and-rules/Rule/)

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_diagonal_matrix_q.c`](https://github.com/stblake/mathilda/blob/main/tests/test_diagonal_matrix_q.c)

## Notes & additional examples

### Notes

Without a `Tolerance` option the test is structural: only literal numeric zeros off the main diagonal count as zero. Returns `False` on non-matrix, ragged, or higher-rank inputs.
