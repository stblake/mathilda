# SchurDecomposition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SchurDecomposition[m]`**

gives the Schur decomposition of a numerical square matrix m as a list {q, t}, where q is orthonormal (unitary) and t is block upper-triangular (the Schur form), so that m == q . t . ConjugateTranspose\[q\].

**`SchurDecomposition[{m, a}]`**

gives the generalized (QZ) Schur decomposition as {q, s, p, t} with q, p orthonormal and s, t upper-triangular, so that m == q . s . ConjugateTranspose\[p\] and a == q . t . ConjugateTranspose\[p\].

<details>
<summary>Notes</summary>

Options: Pivoting -\> True also returns a scaling/permutation matrix d with m . d == d . q . t . ConjugateTranspose\[q\]. RealBlockDiagonalForm -\> True (default) keeps t real with 2x2 blocks for complex-conjugate eigenvalue pairs; -\> False makes t complex upper-triangular. TargetStructure -\> "Dense" (default) | "Structured"; both return dense matrices. SchurDecomposition is numerical: machine-precision Real/complex matrices (and NDArray / packed arrays) use LAPACK, and an arbitrary-precision real matrix uses the in-house MPFR Hessenberg + Francis QR at the input precision.  Arbitrary precision is native only for the standard real form with the default options; the complex, generalized, RealBlockDiagonalForm -\> False, and Pivoting cases are computed at machine precision even for a high-precision input.  A non-numeric (symbolic) matrix, for which there is no closed-form Schur decomposition, is left unevaluated, as is a non-square argument.

</details>

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= {q, t} = SchurDecomposition[{{2.7, 4.8, 8.1}, {-.6, 0, 0}, {.1, 0, .3}}]; Chop[{{2.7, 4.8, 8.1}, {-.6, 0, 0}, {.1, 0, .3}} - q . t . ConjugateTranspose[q]]
Out[1]= {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}

In[2]:= {q, t} = SchurDecomposition[{{1.81066, 0.31066, 1.5}, {-0.53033, 2.03033, 0.43934}, {-0.96967, -0.53033, 2.56066}}]; {UpperTriangularMatrixQ[t, -1], UpperTriangularMatrixQ[t]}
Out[2]= {True, False}

In[3]:= {q, s, p, t} = SchurDecomposition[{{{.5, 1}, {1.5, 2}}, {{2.5, 3}, {3.5, 4}}}]; Chop[{{.5, 1}, {1.5, 2}} - q . s . ConjugateTranspose[p]]
Out[3]= {{0, 0}, {0, 0}}
```

### Applications (3)

Gives {q, t} with q orthonormal and t upper-triangular

```mathematica
In[4]:= SchurDecomposition[{{1, 2}, {3, 4}}]
Out[4]= {{{-0.824565, -0.565767}, {0.565767, -0.824565}}, {{-0.372281, -1.0}, {0.0, 5.37228}}}
```

An already-diagonal matrix, computed numerically

```mathematica
In[5]:= SchurDecomposition[{{2, 0}, {0, 3}}]
Out[5]= {{{1.0, 0.0}, {0.0, 1.0}}, {{2.0, 0.0}, {0.0, 3.0}}}
```

Complex t, eigenvalues 1 +- 2 I on the diagonal

```mathematica
In[6]:= SchurDecomposition[N[{{3, -2}, {4, -1}}], RealBlockDiagonalForm -> False]
Out[6]= {{{-0.408248 - 0.408248*I, -0.689898 - 0.436701*I}, {-0.816497, 0.563299 - 0.126599*I}}, {{1.0 + 2.0*I, 4.44949 + 0.44949*I}, {0.0, 1.0 - 2.0*I}}}
```

## Algorithm

schurdecomp.c -- SchurDecomposition dispatcher.

Parses options, distinguishes the standard form SchurDecomposition[m] from the generalized form SchurDecomposition[{m, a}], classifies the numeric precision of the input, and routes to the right kernel:

```text
  - non-numeric (symbolic) matrix        -> NULL (call left unevaluated;
                                             a generic matrix has no
                                             closed-form Schur decomposition)
  - standard, real, arbitrary precision  -> schur_mpfr_standard_real
    (RealBlockDiagonalForm -> True, no Pivoting; also the no-LAPACK fallback)
  - everything else numeric              -> schur_machine_standard /
                                             schur_machine_generalized (LAPACK)
```

An NDArray / packed matrix is handled transparently: the machine kernel loads it via numarray.c's na_load_matrix (which accepts both an NDArray and a boxed

```text
List-of-Lists), and schur_matrix_order reads its rank-2 shape directly.  So no
```

separate ndla_* guard is needed -- but the head IS on src/pack.c's AWARE list, so the packing gate hands a packed argument straight through as an NDArray rather than materialising it (see docs/design/packed_arrays.md).

```text
Memory contract: standard builtin ownership (SPEC.md §4).  Never frees `res`.
```

## Implementation notes

**Algorithm.** `builtin_schurdecomposition` is a dispatcher. It parses the options (`Pivoting`, `RealBlockDiagonalForm` default `True`, `TargetStructure -> "Dense" | "Structured"`) with `schur_parse_options`, distinguishes the standard form `SchurDecomposition[m] -> {q, t}` (with `m == q . t . ConjugateTranspose[q]`, `q` orthonormal/unitary and `t` block upper-triangular) from the generalized QZ form `SchurDecomposition[{m, a}] -> {q, s, p, t}`, classifies the input's numeric precision, and routes: a non-numeric (symbolic) matrix returns `NULL` — a generic matrix has no closed-form Schur decomposition, so the call is left unevaluated; a standard, real, arbitrary-precision matrix goes to `schur_mpfr_standard_real` (which is also the no-LAPACK fallback); everything else numeric goes to `schur_machine_standard` / `schur_machine_generalized` via LAPACK. `RealBlockDiagonalForm -> True` keeps `t` real with 2×2 blocks for complex-conjugate eigenvalue pairs; `-> False` makes `t` complex upper-triangular; `Pivoting -> True` additionally returns a scaling/permutation matrix.

**Data structures.** The result is a nested `List` of machine-precision (or MPFR) reals or `Complex` values. The machine kernel loads the matrix through `numarray.c`'s `na_load_matrix`, which accepts both an `NDArray` and a boxed `List`-of-`List`s, and `schur_matrix_order` reads the rank-2 shape directly — so no separate `ndla_*` guard is needed, but the head is on `src/pack.c`'s `AWARE` list so the gate hands a packed argument straight through as an `NDArray` rather than materialising it. The LAPACK machine kernel (`schurdecomp_machine.c`) and the MPFR twin (`schurdecomp_mpfr.c`) are the two numeric backends.

**Complexity / limits.** `O(n^3)` for the LAPACK QR-iteration backend (Hessenberg reduction followed by the Francis double-shift QR / QZ sweep). The decomposition is defined only for numerical square matrices — a symbolic matrix is left unevaluated. The generalized form requires the two-matrix argument `{m, a}`; the arbitrary-precision MPFR path covers the standard real case at the input precision.

- `Protected`.
- **Numerical only** — a generic matrix has no closed-form symbolic Schur
  form, so a non-numeric (symbolic) matrix is left unevaluated. Input families:
  - machine-precision Real / complex matrices, and `NDArray` / packed arrays
    (read straight off the buffer);
  - exact integer / rational matrices (computed at machine precision, since
    there is no exact Schur form);
  - arbitrary-precision MPFR matrices — supported for the **standard real**
    form (output at the input precision).
- Options:
  - `Pivoting -> True` returns `{q, t, d}` with an extra scaling / permutation
    matrix `d`, such that `m . d == d . q . t . ConjugateTranspose[q]` (the
    balancing transform, applies to the single-matrix form).
  - `RealBlockDiagonalForm -> True` (default) keeps `t` real with 2×2 blocks
    for complex-conjugate eigenvalue pairs (so `t` is block — not strictly —
    upper-triangular); `-> False` makes `t` complex upper-triangular.
  - `TargetStructure -> "Dense"` (default) `| "Structured"` — both return
    dense matrices (Mathilda has no distinct structured-matrix type; `t` is
    already triangular and `q` orthonormal, so nothing is lost).
- Algorithm:
  - **Machine.** LAPACK `dgees` (real) / `zgees` (complex) for the standard
    form and `dgges` / `zgges` (the QZ algorithm) for the generalized form.
    `RealBlockDiagonalForm -> False` on a real matrix loads it complex and
    uses the `z` driver. `Pivoting -> True` balances with `dgebal` and
    reconstructs `d` by back-transforming the identity with `dgebak`. A
    100×100 real matrix takes a few milliseconds.
  - **Arbitrary precision.** The standard real form reuses the in-house
    Hessenberg + Francis QR (`eigen_schur_real_mpfr`), which accumulates the
    orthogonal Schur vectors alongside the quasi-triangular Schur form.
- Result-fidelity notes: the Schur form is not unique (`q`, `t` are fixed only
  up to the eigenvalue ordering and unitary freedom within degenerate blocks),
  so the contract is the residual `m ≈ q . t . ConjugateTranspose[q]` with
  `q` orthonormal and `t` (quasi-)triangular, not bit-reproducible factors.
- Limitations: at arbitrary precision only the standard real form is native;
  a complex or generalized arbitrary-precision matrix falls back to machine
  precision (a complex MPFR QR / an MPFR QZ would be needed for full
  precision).
- Issues `SchurDecomposition::argx` for zero arguments and
  `SchurDecomposition::sqma` for an argument that is neither a non-empty
  square matrix nor a pair of square matrices; both leave the call
  unevaluated.
- Not lowered by `Compile[]` (the result is a heterogeneous tuple of
  matrices, like `QRDecomposition` / `JordanDecomposition`); a packed /
  `NDArray[…]` argument is accepted (on `pack.c`'s AWARE list).

**Attributes:** `Protected`.

## References

**See also:** [NDArray](../../linear-algebra/NDArray/), [QRDecomposition](../../linear-algebra/QRDecomposition/), [JordanDecomposition](../../linear-algebra/JordanDecomposition/)

- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed. (Johns Hopkins, 2013), ch. 7 — The Unsymmetric Eigenvalue Problem.
- Source: [`src/linalg/schurdecomp.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/schurdecomp.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_schurdecomp.c`](https://github.com/stblake/mathilda/blob/main/tests/test_schurdecomp.c)

## Notes & additional examples

### Notes

`SchurDecomposition[m]` returns `{q, t}` with `q` orthonormal (unitary) and `t`
block upper-triangular, so that `m == q . t . ConjugateTranspose[q]`. By default
`RealBlockDiagonalForm -> True` keeps `t` real, using a 2×2 block for each
complex-conjugate eigenvalue pair; `-> False` makes `t` complex upper-triangular
with the eigenvalues on its diagonal. `Pivoting -> True` additionally returns a
scaling/permutation matrix.

The decomposition is numerical: a symbolic matrix has no closed-form Schur form
and is left unevaluated. `SchurDecomposition[{m, a}]` gives the generalized (QZ)
decomposition `{q, s, p, t}`. A machine or packed/`NDArray` matrix is read
straight off its buffer by the LAPACK kernel; arbitrary-precision real inputs use
an MPFR backend.
