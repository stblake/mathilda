# SymmetricMatrixQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SymmetricMatrixQ[m]
    gives True if m is explicitly symmetric (m == Transpose[m]),
    and False otherwise.

Options:
    SameTest  -> Automatic   function used to test equality of entries.
    Tolerance -> Automatic   numeric tolerance for approximate matrices.

With SameTest -> f, entries m[i,j] and m[j,i] are taken to be equal
when f[m[i,j], m[j,i]] gives True.  With Tolerance -> t, entries are
accepted when Abs[m[i,j] - m[j,i]] <= t.  SymmetricMatrixQ uses the
definition m^T == m for both real- and complex-valued matrices, so a
complex symmetric matrix need not be Hermitian.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= SymmetricMatrixQ[{{1, 2}, {2, 3}}]
Out[1]= True

In[2]:= SymmetricMatrixQ[{{a, b, c}, {b, d, e}, {c, e, f}}]
Out[2]= True

In[3]:= SymmetricMatrixQ[{{1 + I, 2 - 3 I}, {2 - 3 I, 2 - 3 I}}]
Out[3]= True

In[4]:= SymmetricMatrixQ[{{1, 3 + 4 I}, {3 - 4 I, 2}}]   (* Hermitian, not symmetric *)
Out[4]= False
```

## Implementation notes

`builtin_symmetric_matrix_q` first applies the same square-matrix shape gate as `SquareMatrixQ`, then walks the strict upper triangle checking `m[i,j] == m[j,i]`. The comparison defaults to structural `expr_eq`, but a `SameTest -> f` option uses `symmetric_pair_sametest` and a `Tolerance -> t` option uses `symmetric_pair_tolerance`. Returns `False` on any shape rejection or mismatch; unrecognised options leave the call unevaluated.

- `Protected`.
- Default test is structural via `expr_eq`; the diagonal is exempt

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= SymmetricMatrixQ[{{1, 2}, {2, 1}}]
Out[1]= True

In[2]:= SymmetricMatrixQ[{{1, 2}, {3, 4}}]
Out[2]= False
```

A complex symmetric matrix is symmetric without being Hermitian:

```mathematica
In[1]:= SymmetricMatrixQ[{{1, I}, {I, 1}}]
Out[1]= True

In[2]:= HermitianMatrixQ[{{1, I}, {I, 1}}]
Out[2]= False
```

`Tolerance` accepts numerically near-symmetric matrices:

```mathematica
In[1]:= SymmetricMatrixQ[{{1.0, 2.0001}, {2.0, 1.0}}, Tolerance -> 0.001]
Out[1]= True
```

A custom `SameTest` relaxes equality of off-diagonal entries:

```mathematica
In[1]:= SymmetricMatrixQ[{{1, 2}, {3, 4}}, SameTest -> (Abs[#1 - #2] <= 1 &)]
Out[1]= True
```

### Notes

`SymmetricMatrixQ` uses the definition `m^T == m` for both real and complex matrices, so a complex symmetric matrix need not be Hermitian. Use the `SameTest` or `Tolerance` options for approximate or custom equality.
