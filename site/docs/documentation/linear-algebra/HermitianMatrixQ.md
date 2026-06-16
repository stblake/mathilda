# HermitianMatrixQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
HermitianMatrixQ[m]
    gives True if m is explicitly Hermitian (m == ConjugateTranspose[m]),
    and False otherwise.

Options:
    SameTest  -> Automatic   function used to test equality of entries.
    Tolerance -> Automatic   numeric tolerance for approximate matrices.

With SameTest -> f, entries m[i,j] and Conjugate[m[j,i]] are taken to be
equal when f[m[i,j], Conjugate[m[j,i]]] gives True.  With Tolerance -> t,
entries are accepted when Abs[m[i,j] - Conjugate[m[j,i]]] <= t.
Diagonal entries must satisfy the same test (i.e. be purely real for
numeric matrices).
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_hermitian_matrix_q` tests whether a matrix equals its conjugate transpose, i.e. `m[i,j] == Conjugate[m[j,i]]`. After validating that the argument is a non-empty square `List` of `List`s with no deeper nesting (returning `False` otherwise), it walks the upper triangle including the diagonal (the pair test is symmetric under `(i,j)↔(j,i)`) and checks each pair with one of three predicates: the default structural test (`hermitian_pair_structural`, exact for symbolic/exact-numeric entries), a user `SameTest -> f`, or `Tolerance -> t` (accepting pairs with `Abs[a - Conjugate[b]] <= t`). `SameTest`/`Tolerance` of `Automatic` fall through to the structural test; any unrecognised option leaves the call unevaluated. Returns `True`/`False`.

- `Protected`.
- Default test is structural: it accepts (Conjugate[a], a) / (a, Conjugate[a])

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

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

### Notes

A matrix is Hermitian when `m == ConjugateTranspose[m]`; off-diagonal entries must be conjugates of their transpose partners and diagonal entries must be real. For real matrices this coincides with `SymmetricMatrixQ`.
