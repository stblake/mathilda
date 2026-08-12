# HermitianMatrixQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HermitianMatrixQ[m]`**

gives True if m is explicitly Hermitian (m == ConjugateTranspose\[m\]), and False otherwise.

<details>
<summary>Notes</summary>

Options: SameTest  -\> Automatic   function used to test equality of entries. Tolerance -\> Automatic   numeric tolerance for approximate matrices. With SameTest -\> f, entries m\[i,j\] and Conjugate\[m\[j,i\]\] are taken to be equal when f\[m\[i,j\], Conjugate\[m\[j,i\]\]\] gives True.  With Tolerance -\> t, entries are accepted when Abs\[m\[i,j\] - Conjugate\[m\[j,i\]\]\] \<= t. Diagonal entries must satisfy the same test (i.e. be purely real for numeric matrices).

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= HermitianMatrixQ[{{1, 3 + 4 I}, {3 - 4 I, 2}}]
Out[1]= True

In[2]:= HermitianMatrixQ[{{0, a, b}, {Conjugate[a], 1, c}, {Conjugate[b], Conjugate[c], -1}}]
Out[2]= True

In[3]:= HermitianMatrixQ[{{1, 2 I}, {2 I, 3}}]
Out[3]= False
```

### Options (1)

```mathematica
In[4]:= HermitianMatrixQ[{{1.0, 2.0 + 0.01 I}, {2.0 - 0.02 I, 1.5}}, Tolerance -> 0.1]
Out[4]= True
```

### Applications (6)

```mathematica
In[1]:= HermitianMatrixQ[{{1, I}, {-I, 1}}]
Out[1]= True

In[2]:= HermitianMatrixQ[{{1, 2}, {3, 4}}]
Out[2]= False
```

A diagonal entry that is not real, or an off-diagonal pair that is not a
conjugate pair, breaks Hermiticity:

```mathematica
In[1]:= HermitianMatrixQ[{{1, 2 + I}, {2 + I, 1}}]
Out[1]= False

In[2]:= HermitianMatrixQ[{{2, 3 + I}, {3 - I, 5}}]
Out[2]= True
```

The predicate also handles inexact matrices, and a `Tolerance` option absorbs
floating-point noise on the diagonal:

```mathematica
In[1]:= HermitianMatrixQ[N[{{1, I}, {-I, 1}}]]
Out[1]= True

In[2]:= HermitianMatrixQ[{{1, I}, {-I, 2.0000001}}, Tolerance -> 0.001]
Out[2]= True
```

## Implementation notes

`builtin_hermitian_matrix_q` tests whether a matrix equals its conjugate transpose, i.e. `m[i,j] == Conjugate[m[j,i]]`. After validating that the argument is a non-empty square `List` of `List`s with no deeper nesting (returning `False` otherwise), it walks the upper triangle including the diagonal (the pair test is symmetric under `(i,j)↔(j,i)`) and checks each pair with one of three predicates: the default structural test (`hermitian_pair_structural`, exact for symbolic/exact-numeric entries), a user `SameTest -> f`, or `Tolerance -> t` (accepting pairs with `Abs[a - Conjugate[b]] <= t`). `SameTest`/`Tolerance` of `Automatic` fall through to the structural test; any unrecognised option leaves the call unevaluated. Returns `True`/`False`.

- `Protected`.
- Default test is structural: it accepts (Conjugate[a], a) / (a, Conjugate[a])
  symbolic pairs without requiring our `Conjugate` builtin to fold
  `Conjugate[Conjugate[x]] -> x`.
- Returns `False` (rather than leaving unevaluated) on non-matrix, non-square,
  ragged, empty, or higher-rank tensor inputs.
- Unknown options and non-Rule trailing arguments leave the call unevaluated.

**Attributes:** `Protected`.

## See also

[Conjugate](../../arithmetic/Conjugate/)

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_hermitian_matrix_q.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hermitian_matrix_q.c)
- Tests: [`tests/test_symmetric_matrix_q.c`](https://github.com/stblake/mathilda/blob/main/tests/test_symmetric_matrix_q.c)

## Notes & additional examples

### Notes

A matrix is Hermitian when `m == ConjugateTranspose[m]`; off-diagonal entries must be conjugates of their transpose partners and diagonal entries must be real. For real matrices this coincides with `SymmetricMatrixQ`.
