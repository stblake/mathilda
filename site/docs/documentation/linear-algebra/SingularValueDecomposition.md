# SingularValueDecomposition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SingularValueDecomposition[m]`**

gives the singular value decomposition of a matrix m as a list {u, sigma, v}, where sigma is a diagonal matrix and m == u . sigma . ConjugateTranspose\[v\].  u and v have orthonormal columns.

**`SingularValueDecomposition[m, k]`**

gives the SVD associated with the k largest singular values (or |k| smallest when k is negative).

**`SingularValueDecomposition[m, UpTo[k]]`**

gives the SVD for as many of the k largest singular values as are available (up to MatrixRank\[m\]).

**`SingularValueDecomposition[{m, a}]`**

gives the generalized singular value decomposition of m with respect to a as {{u, ua}, {sigma, sigma\_a}, v} such that m == u . sigma . ConjugateTranspose\[v\] and a == ua . sigma\_a . ConjugateTranspose\[v\].  Uses LAPACK dggsvd3 / zggsvd3 for real / complex machine-precision inputs; exact-numeric input is numericalised to 53 bits and routed through the same path.  High-precision MPFR input is currently downgraded to machine precision and emits the ::gmpdwn warning.  Free-symbolic input emits ::nogsymb and leaves the call unevaluated.

<details>
<summary>Notes</summary>

SingularValueDecomposition works on every input family supported by the rest of the linear-algebra builtins: - exact integer / rational matrices (output stays exact; singular values are Sqrt\[...\] forms when irrational) - complex matrices (u and v are unitary in the Hermitian inner product) - machine-precision Real matrices (LAPACK dgesdd / zgesdd, or dggsvd3 / zggsvd3 for the generalized form) - arbitrary-precision MPFR matrices (one-sided Jacobi SVD at the input precision, real and complex) - free-symbolic matrices (eigendecomposition of m^H . m; the call is left unevaluated when no closed form exists) Options: Tolerance -\> t       :  zero out singular values below t TargetStructure -\> "Dense"            :  u, sigma, v all dense (default) "Structured"       :  sigma returned as DiagonalMatrix\[{..}\] A non-rank-2 or empty matrix emits SingularValueDecomposition::matrix and the call is left unevaluated.  Generalized SVD with mismatched column counts emits ::matdims.  An out-of-range k or UpTo\[k\] emits ::sval. Unknown option keys / values emit ::opts and leave the call unevaluated.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= sv = SingularValueDecomposition[{{1.2, 3.4}, {5.6, 7.8}, {9.0, 1.2}}]; sv[[2]]
Out[1]= {{12.4778, 0.0}, {0.0, 5.65202}, {0.0, 0.0}}
```

"thin" SVD with only the top singular value

```mathematica
In[2]:= SingularValueDecomposition[{{1.1, 2, 5}, {3, -11, 4.2}}, 1]
Out[2]= {{{0.0195749}, {0.999808}}, {{12.1526}}, {{0.248586}, {-0.901763}, {0.353593}}}
```

MPFR Jacobi -- output at the input precision

```mathematica
In[3]:= SingularValueDecomposition[N[{{1, 2}, {3, 4}}, 30]]
Out[3]= {{{0.4045535848337569316424487226274, -0.9145142956773044526791769738123}, {0.9145142956773044526791769738107, 0.4045535848337569316424487226246}}, {{5.464985704219042650451188493292, 0.0}, {0.0, 0.3659661906262578204229643842617}}, {{0.5760484367663207913310985819436, 0.8174155604703632730886523884647}, {0.8174155604703632730886523884647, -0.5760484367663207913310985819436}}}
```

### Options (1)

Explicit Tolerance preserves the tiny singular value

```mathematica
In[4]:= SingularValueDecomposition[{{1.0, 0}, {1.0, 10^-14}}, Tolerance -> 10^-15]
Out[4]= {{{-0.707107, -0.707107}, {-0.707107, 0.707107}}, {{1.41421, 0.0}, {0.0, 7.07107e-15}}, {{-1.0, -5e-15}, {-5e-15, 1.0}}}
```

### Applications (3)

```mathematica
In[5]:= SingularValueDecomposition[{{2, 0}, {0, 3}}]
Out[5]= {{{0, 1}, {1, 0}}, {{3, 0}, {0, 2}}, {{0, 1}, {1, 0}}}

In[6]:= With[{r = SingularValueDecomposition[N[{{1, 2}, {3, 4}}]]}, Chop[r[[1]] . r[[2]] . Transpose[r[[3]]]]]
Out[6]= {{1.0, 2.0}, {3.0, 4.0}}

In[7]:= SingularValueDecomposition[N[{{1, 2}, {3, 4}}]]
Out[7]= {{{-0.404554, -0.914514}, {-0.914514, 0.404554}}, {{5.46499, 0.0}, {0.0, 0.365966}}, {{-0.576048, 0.817416}, {-0.817416, -0.576048}}}
```

## Options & behaviour

**Generalized SVD (`{m, a}`) algorithm**:
- Real or complex machine-precision input is forwarded to LAPACK
  `dggsvd3` (real) or `zggsvd3` (complex).  The output of LAPACK,
  `A = U . D1 . (0 R) . Q^H` and `B = V . D2 . (0 R) . Q^H`, maps to
  Mathematica's form with `u = U`, `ua = V`, `v = Q`, and
  `sigma = D1 . (0 R)` and `sigma_a = D2 . (0 R)` materialised as
  dense rectangular matrices.  Mathilda handles both LAPACK
  layouts (`M >= K+L` and `M < K+L`, the latter splits `R` across
  the destroyed `A` and `B` buffers).
- Exact-numeric symbolic input (Integer / Rational / Complex of
  numeric) is numericalised to 53-bit Reals and re-dispatched
  through the LAPACK path.  Inputs with free symbolic content emit
  `::nogsymb` and leave the call unevaluated.
- High-precision MPFR input is also routed through the LAPACK path
  (with the `::gmpdwn` one-shot precision-downgrade warning) until a
  native MPFR generalized kernel lands.

> Implementation is split across:
> - `src/linalg/svdecomp.c` -- builtin entry, positional + option
>   parser (`svd_parse_args`), top-level dispatcher (`svd_dispatch`),
>   symbolic dispatcher + core (`svd_symbolic_core`), and the shared
>   post-processing helper `svd_apply_postprocess` (truncation /
>   tolerance / TargetStructure).
> - `src/linalg/svdecomp_machine.c` -- LAPACK fast path
>   (`svd_machine_dispatch`).  Loads to a column-major double buffer,
>   calls `mat_lapack_dgesdd` / `mat_lapack_zgesdd`, builds `u`,
>   rectangular `sigma`, and `v = ConjugateTranspose[V^H]` as Mathilda
>   lists.
> - `src/linalg/svdecomp_mpfr.c` -- MPFR one-sided Jacobi kernel
>   (`svd_mpfr_dispatch`).  Column-major MPFR arrays at the input
>   precision; orthonormal completion of `u` via Gram-Schmidt; emits
>   MPFR zeros (not Real 0.0) for off-diagonal sigma entries so
>   downstream Dot products preserve the working precision.

## Algorithm

svdecomp.c

```text
SingularValueDecomposition[m]                       -- {u, sigma, v}
SingularValueDecomposition[m, k]                    -- k largest (k<0 -> |k| smallest)
SingularValueDecomposition[m, UpTo[k]]              -- up to k largest
SingularValueDecomposition[{m, a}]                  -- generalized form
Options:  Tolerance -> t,  TargetStructure -> "Dense" | "Structured"

Strategy.  One algorithmic core per numeric domain:

  - Exact / symbolic   -> eigendecomposition of m^H . m (or m . m^H,
                          whichever is smaller), with a 2x2 closed-form
                          fast path.  Outputs Sqrt[lambda] for the
                          singular values; columns of v are orthonormal
                          eigenvectors; u = m . v . Sigma^-1 with the
                          null space completed via the existing
                          qr_symbolic_core orthogonal completion.

  - Inexact, min_bits <= 53  -> LAPACK dgesdd / zgesdd (standard) or
                                dggsvd3 / zggsvd3 (generalized).

  - Inexact, min_bits > 53   -> one-sided Jacobi SVD over MPFR arrays.
```

All three paths feed the result through svd_apply_postprocess so the truncation, tolerance, and TargetStructure logic lives in one place.

```text
Memory contract.  Standard builtin contract.  This file does NOT call
```

expr_free(res) - the evaluator owns `res` and frees it on a non-NULL return (MEMORY.md / SPEC.md Sec. 4.1).

## Implementation notes

**Algorithm.** `builtin_singularvaluedecomposition` returns `{u, sigma, v}` with `m == u.sigma.ConjugateTranspose[v]`, supporting `SingularValueDecomposition[m, k]`, `UpTo[k]`, the generalized `{m, a}` form, and `Tolerance`/`TargetStructure` options. `svd_dispatch` selects a kernel per numeric domain; all three feed `svd_apply_postprocess`, which centralises truncation, tolerance, and `TargetStructure` handling.

- **Exact / symbolic** (`svd_symbolic_dispatch` → `svd_symbolic_core`): forms the smaller Gram matrix `B = mᴴm` (p × p) when `n >= p`, else `B = m mᴴ` (n × n), and eigendecomposes it through the evaluator's `Eigenvalues`/`Eigenvectors`, with a 2×2 closed-form fast path. The singular values are `Sqrt[lambda]`; the eigenvectors give the primary factor (`v` if `mᴴm`, else `u`); the other factor is recovered as `m.v.Sigma⁻¹` (resp. `mᴴ.u/sigma_i`) for the non-zero singular values, with the remaining columns filled by `qr_symbolic_core`'s orthogonal completion to span the null space. If the eigendecomposition has no closed form, this path returns NULL.
- **Inexact, `min_bits <= 53`** (`svd_machine_dispatch`): LAPACK `dgesdd`/`zgesdd` (divide-and-conquer) for the standard form, `dggsvd3`/`zggsvd3` for the generalized `{m, a}` form.
- **Inexact, `min_bits > 53`** (`svd_mpfr_dispatch`): one-sided Jacobi SVD over MPFR arrays (Demmel-Veselić), preceded by a QR/Paige-Van Loan reduction.

**Data structures.** `Expr*` `List`-of-`List` matrices for the symbolic path (Gram matrix, eigenpairs, `Sqrt`-valued `sigma`); dense `double` / interleaved-complex LAPACK buffers for the machine path; arbitrary-precision MPFR arrays for the high-precision path. The exact path picks the smaller of `mᴴm` and `m mᴴ` to keep the eigenproblem small. Inexact input uses the rationalise / numericalise round-trip at the minimum input precision.

**Complexity / limits.** Exact path cost is dominated by the symbolic eigendecomposition of a min(n,p)-square matrix and only succeeds when that closes; the generalized `{m, a}` form has no symbolic kernel and requires the machine path. Machine path is LAPACK-bound (~O(n p · min(n,p))). The `TargetStructure -> "Structured"` head is not fully realised (results are returned dense).

- `Protected`.
- Options: `Tolerance -> t` zeroes out singular values with
  `|sigma_i| < t`; `TargetStructure -> "Dense" | "Structured"`
  (`"Structured"` wraps `sigma` as `DiagonalMatrix[{...}]` -- which
  currently evaluates back to a dense matrix because Mathilda has no
  unevaluated structured-matrix head yet).
- Works on every input family:
  - exact integer / rational matrices (output stays exact, with
    `Sqrt[lambda]` forms for the singular values when the
    characteristic polynomial of `m^H . m` doesn't simplify further);
  - complex matrices (`u` and `v` are unitary in the Hermitian inner
    product);
  - machine-precision Real matrices (output Real at machine precision
    via LAPACK `dgesdd` / `zgesdd`);
  - arbitrary-precision MPFR matrices (output at the input precision
    via a one-sided Jacobi SVD).
- Algorithm: dispatched on leaf precision.
  - **MachinePrecision inputs (`min_bits <= 53`)** use LAPACK
    divide-and-conquer SVD (`dgesdd` for real, `zgesdd` for complex)
    with `jobz = 'A'` so the full square `u` and `V^H` are returned.
    `V^H` is conjugate-transposed in place to produce `v`.
  - **MPFR inputs (`min_bits > 53`)** use a one-sided Jacobi SVD: sweep
    over column pairs of `A`, apply the rotation that diagonalises
    `[<a_i,a_i>, <a_i,a_j>; <a_i,a_j>, <a_j,a_j>]` to columns `i` and
    `j` of both `A` and the accumulating `V`.  After convergence (off-
    diagonal norm below `2^(-bits/2)`), column norms of `A` are the
    singular values and `A[:, i] / sigma_i` are the left singular
    vectors.  Orthonormal completion of the `n - rank` remaining
    columns of `u` uses Gram-Schmidt at the working precision.  Real
    and complex MPFR inputs are both supported: the complex path uses
    paired (re, im) MPFR buffers throughout and factors each 2x2
    rotation as a phase rotation (that makes the column inner product
    real-positive) followed by a real Jacobi rotation.
  - **Exact / symbolic inputs** compute the eigendecomposition of
    `m^H . m` (or `m . m^H`, whichever is larger so the null-space
    eigenvectors come back automatically), normalise the eigenvectors
    individually, derive the secondary singular vectors via
    `secondary = m . v_i / sigma_i`, and orthonormal-complete the
    smaller side via the existing `qr_symbolic_core` Gram-Schmidt
    pipeline.

    Because the kernel is an eigendecomposition, it inherits
    `Eigenvectors`' handling of irrational algebraic eigenvalues. A
    full-rank exact matrix whose gram has an **irreducible cubic**
    characteristic polynomial — `{{1,2,3},{4,5,6},{7,8,10}}`, for
    instance — used to come back unevaluated: the gram's eigenvalues are
    *casus irreducibilis*, `Eigenvectors` returned zero vectors, and
    normalising divided by a zero norm. Fixed 2026-08-03 by the adjugate
    route documented under `Eigenvectors`; it now reconstructs to
    `1.7e-34` in 0.81 s. **Radical depth, not rank deficiency, was the
    trigger** — that is a separate failure, fixed 2026-08-02.
  - When BLAS/LAPACK is unavailable at build time (`USE_LAPACK=0`),
    machine-precision inputs transparently route to the symbolic
    kernel; similarly `USE_MPFR=0` routes MPFR inputs to symbolic.
    Correctness is preserved across every combination; only
    performance changes.
- Issues `SingularValueDecomposition::matrix` and returns unevaluated
  if the argument is not a non-empty rank-2 tensor.
- Issues `SingularValueDecomposition::sval` and returns unevaluated
  for an out-of-range `k` or `UpTo[k]`.
- Issues `SingularValueDecomposition::opts` and returns unevaluated
  for an unknown option key or value.
- Issues `SingularValueDecomposition::matdims` and returns unevaluated
  when the generalized form `{m, a}` is called with mismatched column
  counts.
- Issues `SingularValueDecomposition::nogsymb` and returns unevaluated
  when the generalized form is called with free symbolic content (no
  closed-form symbolic kernel exists for the generalized GSVD).
- Issues `SingularValueDecomposition::gmpdwn` (one-shot) when the
  generalized form is called on high-precision MPFR input: the result
  is computed at machine precision via LAPACK because Mathilda does
  not yet have a native MPFR Paige/Van Loan kernel.

**Attributes:** `Protected`.

## References

**See also:** [Eigenvectors](../../linear-algebra/Eigenvectors/)

- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed. (Johns Hopkins, 2013).
- J. Demmel and K. Veselić, "Jacobi's Method is More Accurate than QR", SIAM J. Matrix Anal. Appl. 13 (1992).
- Source: [`src/linalg/svdecomp.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/svdecomp.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_ndarray_linalg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_linalg.c)
- Tests: [`tests/test_singularvaluedecomposition.c`](https://github.com/stblake/mathilda/blob/main/tests/test_singularvaluedecomposition.c)
- Tests: [`tests/test_singularvaluedecomposition_machine.c`](https://github.com/stblake/mathilda/blob/main/tests/test_singularvaluedecomposition_machine.c)
- Tests: [`tests/test_singularvaluedecomposition_mpfr.c`](https://github.com/stblake/mathilda/blob/main/tests/test_singularvaluedecomposition_mpfr.c)

## Notes & additional examples

### Notes

`SingularValueDecomposition[m]` returns `{u, sigma, v}` with `m == u . sigma . ConjugateTranspose[v]`, where `u` and `v` have orthonormal columns and `sigma` is diagonal with the singular values in descending order. Exact integer/rational matrices stay exact (singular values appear as `Sqrt[...]` forms when irrational); machine-precision real and complex inputs route through LAPACK, and high-precision MPFR input uses a one-sided Jacobi SVD. A two-argument form `SingularValueDecomposition[m, k]` keeps only the `k` largest singular values, and `SingularValueDecomposition[{m, a}]` computes the generalized SVD.
