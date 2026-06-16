# SingularValueDecomposition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SingularValueDecomposition[m]
    gives the singular value decomposition of a matrix m as a
    list {u, sigma, v}, where sigma is a diagonal matrix and
    m == u . sigma . ConjugateTranspose[v].  u and v have
    orthonormal columns.

SingularValueDecomposition[m, k]
    gives the SVD associated with the k largest singular values
    (or |k| smallest when k is negative).
SingularValueDecomposition[m, UpTo[k]]
    gives the SVD for as many of the k largest singular values
    as are available (up to MatrixRank[m]).

SingularValueDecomposition[{m, a}]
    gives the generalized singular value decomposition of m
    with respect to a as {{u, ua}, {sigma, sigma_a}, v} such
    that m == u . sigma . ConjugateTranspose[v] and
    a == ua . sigma_a . ConjugateTranspose[v].  Uses LAPACK
    dggsvd3 / zggsvd3 for real / complex machine-precision
    inputs; exact-numeric input is numericalised to 53 bits and
    routed through the same path.  High-precision MPFR input is
    currently downgraded to machine precision and emits the
    ::gmpdwn warning.  Free-symbolic input emits ::nogsymb and
    leaves the call unevaluated.

SingularValueDecomposition works on every input family supported
by the rest of the linear-algebra builtins:
    - exact integer / rational matrices (output stays exact;
      singular values are Sqrt[...] forms when irrational)
    - complex matrices (u and v are unitary in the Hermitian
      inner product)
    - machine-precision Real matrices (LAPACK dgesdd / zgesdd,
      or dggsvd3 / zggsvd3 for the generalized form)
    - arbitrary-precision MPFR matrices (one-sided Jacobi SVD
      at the input precision, real and complex)
    - free-symbolic matrices (eigendecomposition of m^H . m;
      the call is left unevaluated when no closed form exists)

Options:
    Tolerance -> t       :  zero out singular values below t
    TargetStructure ->
      "Dense"            :  u, sigma, v all dense (default)
      "Structured"       :  sigma returned as DiagonalMatrix[{..}]

A non-rank-2 or empty matrix emits
SingularValueDecomposition::matrix and the call is left
unevaluated.  Generalized SVD with mismatched column counts
emits ::matdims.  An out-of-range k or UpTo[k] emits ::sval.
Unknown option keys / values emit ::opts and leave the call
unevaluated.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= SingularValueDecomposition[{{1, 2}, {1, 2}}]
Out[1]= {{{1/Sqrt[2], 1/Sqrt[2]}, {1/Sqrt[2], -1/Sqrt[2]}}, {{Sqrt[10], 0}, {0, 0}}, {{1/Sqrt[5], -2/Sqrt[5]}, {2/Sqrt[5], 1/Sqrt[5]}}}

In[2]:= SingularValueDecomposition[{{1.1, 2, 5}, {3, -11, 4.2}}, 1]
Out[2]= {{{0.0195749}, {0.999808}}, {{12.1526}}, {{0.248586}, {-0.901763}, {0.353593}}}

In[3]:= SingularValueDecomposition[N[{{1, 2}, {3, 4}}, 30]]
Out[3]= {{{0.4045535848337569316424487226274, -0.9145142956773044526791769738123}, {0.9145142956773044526791769738107, 0.4045535848337569316424487226246}}, {{5.464985704219042650451188493292, 0.0}, {0.0, 0.3659661906262578204229643842617}}, {{0.5760484367663207913310985819436, 0.8174155604703632730886523884647}, {0.8174155604703632730886523884647, -0.5760484367663207913310985819436}}}

In[4]:= SingularValueDecomposition[{{1.0, 0}, {1.0, 10^-14}}, Tolerance -> 10^-15]
Out[4]= {{{-0.707107, -0.707107}, {-0.707107, 0.707107}}, {{1.41421, 0.0}, {0.0, 7.07107e-15}}, {{-1.0, -5e-15}, {-5e-15, 1.0}}}
```

## Implementation notes

**Algorithm.** `builtin_singularvaluedecomposition` returns `{u, sigma, v}` with `m == u.sigma.ConjugateTranspose[v]`, supporting `SingularValueDecomposition[m, k]`, `UpTo[k]`, the generalized `{m, a}` form, and `Tolerance`/`TargetStructure` options. `svd_dispatch` selects a kernel per numeric domain; all three feed `svd_apply_postprocess`, which centralises truncation, tolerance, and `TargetStructure` handling.

- **Exact / symbolic** (`svd_symbolic_dispatch` → `svd_symbolic_core`): forms the smaller Gram matrix `B = mᴴm` (p × p) when `n >= p`, else `B = m mᴴ` (n × n), and eigendecomposes it through the evaluator's `Eigenvalues`/`Eigenvectors`, with a 2×2 closed-form fast path. The singular values are `Sqrt[lambda]`; the eigenvectors give the primary factor (`v` if `mᴴm`, else `u`); the other factor is recovered as `m.v.Sigma⁻¹` (resp. `mᴴ.u/sigma_i`) for the non-zero singular values, with the remaining columns filled by `qr_symbolic_core`'s orthogonal completion to span the null space. If the eigendecomposition has no closed form, this path returns NULL.
- **Inexact, `min_bits <= 53`** (`svd_machine_dispatch`): LAPACK `dgesdd`/`zgesdd` (divide-and-conquer) for the standard form, `dggsvd3`/`zggsvd3` for the generalized `{m, a}` form.
- **Inexact, `min_bits > 53`** (`svd_mpfr_dispatch`): one-sided Jacobi SVD over MPFR arrays (Demmel-Veselić), preceded by a QR/Paige-Van Loan reduction.

**Data structures.** `Expr*` `List`-of-`List` matrices for the symbolic path (Gram matrix, eigenpairs, `Sqrt`-valued `sigma`); dense `double` / interleaved-complex LAPACK buffers for the machine path; arbitrary-precision MPFR arrays for the high-precision path. The exact path picks the smaller of `mᴴm` and `m mᴴ` to keep the eigenproblem small. Inexact input uses the rationalise / numericalise round-trip at the minimum input precision.

**Complexity / limits.** Exact path cost is dominated by the symbolic eigendecomposition of a min(n,p)-square matrix and only succeeds when that closes; the generalized `{m, a}` form has no symbolic kernel and requires the machine path. Machine path is LAPACK-bound (~O(n p · min(n,p))). The `TargetStructure -> "Structured"` head is not fully realised (results are returned dense).

- `Protected`.
- Options: `Tolerance -> t` zeroes out singular values with

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed. (Johns Hopkins, 2013).
- J. Demmel and K. Veselić, "Jacobi's Method is More Accurate than QR", SIAM J. Matrix Anal. Appl. 13 (1992).
- Source: [`src/linalg/svdecomp.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/svdecomp.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= SingularValueDecomposition[{{2, 0}, {0, 3}}]
Out[1]= {{{0, 1}, {1, 0}}, {{3, 0}, {0, 2}}, {{0, 1}, {1, 0}}}
```

The middle factor is the diagonal matrix of singular values (here `3` and `2`, in descending order); the outer factors permute the axes so the larger value comes first.

```mathematica
In[1]:= With[{r = SingularValueDecomposition[N[{{1, 2}, {3, 4}}]]}, Chop[r[[1]] . r[[2]] . Transpose[r[[3]]]]]
Out[1]= {{1.0, 2.0}, {3.0, 4.0}}
```

The defining identity `m == u . sigma . ConjugateTranspose[v]` reconstructs the original matrix exactly (up to floating-point noise removed by `Chop`).

```mathematica
In[1]:= SingularValueDecomposition[N[{{1, 2}, {3, 4}}]]
Out[1]= {{{-0.404554, -0.914514}, {-0.914514, 0.404554}}, {{5.46499, 0.0}, {0.0, 0.365966}}, {{-0.576048, 0.817416}, {-0.817416, -0.576048}}}
```

### Notes

`SingularValueDecomposition[m]` returns `{u, sigma, v}` with `m == u . sigma . ConjugateTranspose[v]`, where `u` and `v` have orthonormal columns and `sigma` is diagonal with the singular values in descending order. Exact integer/rational matrices stay exact (singular values appear as `Sqrt[...]` forms when irrational); machine-precision real and complex inputs route through LAPACK, and high-precision MPFR input uses a one-sided Jacobi SVD. A two-argument form `SingularValueDecomposition[m, k]` keeps only the `k` largest singular values, and `SingularValueDecomposition[{m, a}]` computes the generalized SVD.
