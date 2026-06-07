# LUDecomposition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
LUDecomposition[m]
    gives the LU decomposition of a square matrix m as a list
    {lu, p, c}.  The first element lu is the combined Doolittle
    factor matrix: its strictly-lower triangle is L (with an
    implicit unit diagonal) and its upper triangle is U.  The
    second element p is a 1-indexed row-permutation vector such
    that m[[p]] == l . u where l = LowerTriangularize[lu, -1] +
    IdentityMatrix[n] and u = UpperTriangularize[lu].  The third
    element c is an L-infinity condition-number estimate for
    approximate numerical matrices, or 0 for exact / symbolic m.

LUDecomposition works on every input family supported by the
rest of the linear-algebra builtins:
    - exact integer / rational matrices (output stays exact)
    - complex matrices
    - machine-precision Real matrices (LAPACK dgetrf / zgetrf
      with dgecon / zgecon for the condition estimate)
    - arbitrary-precision MPFR matrices (Householder-free
      Doolittle at the input precision; condition estimate via
      the explicit inverse)
    - free-symbolic matrices (output in closed symbolic form)

The algorithm is Doolittle's elimination with partial row
pivoting.  Numerical inputs use largest-|pivot| selection;
symbolic / exact inputs advance to the next non-zero pivot
only when the natural choice is provably zero.

A singular m emits LUDecomposition::sing and the factorisation
completes with a zero on U's diagonal at the singular step.

A non-square or empty matrix emits LUDecomposition::matsq and
the call is left unevaluated.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- `Protected`.
- Any non-empty rectangular `rows x cols` matrix is accepted.  Empty

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — LU factorisation with partial pivoting.
- L. N. Trefethen and D. Bau III, *Numerical Linear Algebra*, SIAM, 1997 — Gaussian elimination and LU factorisation.
- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= LUDecomposition[{{4, 3}, {6, 3}}]
Out[1]= {{{4, 3}, {3/2, -3/2}}, {1, 2}, 0}
```

```mathematica
In[1]:= LUDecomposition[{{2, 1}, {4, 1}}]
Out[1]= {{{2, 1}, {2, -1}}, {1, 2}, 0}
```

The combined factor reconstructs `m` via `L . U` (here `L = {{1, 0}, {3/2, 1}}`, `U = {{4, 3}, {0, -3/2}}`):

```mathematica
In[1]:= {{1, 0}, {3/2, 1}} . {{4, 3}, {0, -3/2}}
Out[1]= {{4, 3}, {6, 3}}
```

### Notes

`LUDecomposition` returns a list `{lu, p, c}`: `lu` is the combined Doolittle factor whose strictly-lower triangle is `L` (with an implicit unit diagonal) and whose upper triangle is `U`; `p` is the 1-indexed row-permutation vector (here `{1, 2}`, i.e. no swap was needed); and `c` is an `L`-infinity condition estimate that is `0` for exact or symbolic input. The relation is `m[[p]] == L . U`, as the manual reconstruction above confirms. The algorithm is Doolittle elimination with partial pivoting; exact integer inputs keep exact rational factors. A singular matrix emits `LUDecomposition::sing` and completes with a zero pivot on `U`'s diagonal; a non-square or empty matrix emits `LUDecomposition::matsq`.
