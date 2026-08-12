# LinearSolve

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LinearSolve[m, b]`**

finds an x that solves the matrix equation m . x == b.

**`LinearSolve[m, b, Method -> "<name>"]`**

runs a specific elimination algorithm.

<details>
<summary>Notes</summary>

LinearSolve works on both numerical and symbolic matrices. The matrix m may be square or rectangular. The argument b may be a vector or a matrix; when b is a matrix (one column per RHS) LinearSolve returns a matrix of solutions. Higher-rank tensor inputs are also supported: when m has dimensions {d1, ..., d(N-1), n}, b may have dimensions {d1, ..., d(N-1), e1, ..., ep} and the result has dimensions {n, e1, ..., ep}. For under-determined systems LinearSolve returns a particular solution in which the free (non-pivot) variables are taken to be 0; Solve returns the general solution.  When the equation has no solution LinearSolve emits LinearSolve::nosol and returns unevaluated. Accepted method names: "Automatic"                 — alias for "DivisionFreeRowReduction" (default) "DivisionFreeRowReduction"  — Bareiss-like fraction-free Gauss-Jordan on \[m | b\] "OneStepRowReduction"       — classical Gauss-Jordan with division per pivot "CofactorExpansion"         — Cramer's rule via Laplace cofactor expansion (square non-singular m only; LinearSolve::cofnsq / ::cofsng on shape / singularity errors) Default implementation: fraction-free Gauss-Jordan elimination on the augmented matrix \[m | b\] (the Bareiss-like algorithm shared with RowReduce and Inverse), so exact integer / rational / symbolic inputs flow through without any spurious denominator blow-up.

</details>

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= LinearSolve[{{r, s}, {t, u}}, {y, z}]
Out[1]= {(u y)/(-s t + r u) - (s z)/(-s t + r u), -(t y)/(-s t + r u) + (r z)/(-s t + r u)}

In[2]:= LinearSolve[{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}]
Out[2]= {{-3, -4}, {4, 5}}

In[3]:= LinearSolve[{{1, 5}, {2, 6}, {3, 7}, {4, 8}}, {9, 10, 11, 12}]
Out[3]= {-1, 2}

In[4]:= LinearSolve[{{1, 2, 3}, {4, 5, 6}}, {6, 15}]
Out[4]= {0, 3, 0}

In[5]:= a = RandomInteger[{-3, 3}, {2, 3, 6}]; b = RandomInteger[{-3, 3}, {2, 3, 4, 5}]; Dimensions[LinearSolve[a, b]]
Out[5]= {6, 4, 5}
```

### Options (1)

```mathematica
In[6]:= LinearSolve[{{1, 2}, {3, 4}}, {5, 6}, Method -> "CofactorExpansion"]
Out[6]= {-4, 9/2}
```

### Applications (6)

```mathematica
In[1]:= LinearSolve[{{1, 2}, {3, 4}}, {5, 6}]
Out[1]= {-4, 9/2}
```

```mathematica
In[1]:= LinearSolve[{{2, 0}, {0, 4}}, {6, 8}]
Out[1]= {3, 2}
```

```mathematica
In[1]:= LinearSolve[{{1, 1, 1}, {1, 2, 4}, {1, 3, 9}}, {a, b, c}]
Out[1]= {3 a - 3 b + c, -5/2 a + 4 b - 3/2 c, 1/2 a - b + 1/2 c}
```

```mathematica
In[1]:= LinearSolve[{{1, c}, {c, 1}}, {1, 0}]
Out[1]= {1/(1 - c^2), -c/(1 - c^2)}
```

```mathematica
In[1]:= LinearSolve[{{1, 2}, {3, 4}}, {{1, 0}, {0, 1}}]
Out[1]= {{-2, 1}, {3/2, -1/2}}
```

```mathematica
In[1]:= LinearSolve[{{1, 2}, {2, 4}}, {1, 3}]
Out[1]= LinearSolve[{{1, 2}, {2, 4}}, {1, 3}]
```

## Algorithm

matsol.c

Method-aware dispatcher for `RowReduce[m]` and `LinearSolve[m, b]`.

Both builtins accept an optional `Method -> "<name>"` argument (RowReduce as arg 2, LinearSolve as arg 3) and route to one of three explicit algorithms:

```text
  "Automatic"                  -- alias for "DivisionFreeRowReduction"
  "DivisionFreeRowReduction"   -- Bareiss-like fraction-free Gauss-Jordan
  "OneStepRowReduction"        -- classical Gauss-Jordan, one division per
                                  pivot per element of the pivot row
  "CofactorExpansion"          -- LinearSolve: Cramer's rule via Laplace
                                  cofactor expansion;
                                  RowReduce: identity-if-invertible (for a
                                  non-singular square matrix) with fallback
                                  to DivisionFreeRowReduction on
                                  singular / rectangular / empty input.
```

The DivisionFreeRowReduction workers are direct lifts of the previous bodies of `builtin_rowreduce` and `builtin_linearsolve` in src/linalg.c, so existing behaviour is preserved bit-for-bit when the user does not supply a Method option.

Algorithm choice for the four supported matrix families:

```text
  - Machine precision (Real)            -> OneStep is fastest; DivFree fine.
  - Bignum integer                      -> DivFree avoids GCD blow-up.
  - MPFR / arbitrary precision          -> Same as machine precision.
  - Symbolic                            -> DivFree avoids algebraic growth;
                                           CofactorExpansion gives the
                                           textbook Cramer form for small n.
```

## Implementation notes

**Algorithm.** `builtin_linearsolve` solves `m . x = b` and shares its method dispatcher (`MatsolMethod`, `matsol_parse_method_option`) with `RowReduce` and `Inverse`. An optional `Method` argument selects:

- `"DivisionFreeRowReduction"` (default / `Automatic`): `linearsolve_divfree` performs Bareiss-like **fraction-free Gauss-Jordan** elimination on the augmented matrix `[m | b]` with a running pivot product to keep arithmetic exact (no GCD blow-up). After reduction it reads off a particular solution; the implementation handles rectangular/under-determined systems (pivot column bookkeeping per row) and inconsistent systems (emitting the appropriate diagnostic), and supports vector, matrix, and higher-rank RHS by tracking the trailing dimensions of `b`.
- `"OneStepRowReduction"`: classical Gauss-Jordan, one division per pivot.
- `"CofactorExpansion"`: Cramer's rule via Laplace cofactor expansion (small square non-singular systems).

The RHS shape is normalised: a rank-1 `b` against a rank-2 `m` returns a flat solution vector; a matrix RHS returns shape `{c, k}`; higher-rank `b` returns `{c, trail_dims…}`. Non-rectangular `m` emits a `LinearSolve` shape diagnostic.

**Data structures.** Dense flat `Expr**` augmented matrix `r × (c + k)`, row-major; `int* pivot_col_for_row` records pivot columns; the pivot product is a single `Expr*`. All arithmetic is symbolic/exact via the evaluator, so integer, rational, and symbolic systems solve exactly.

**Complexity / limits.** Fraction-free elimination is `O(r·c·(c+k))` arithmetic ops with bounded intermediate-coefficient growth (Bareiss). Inexact input is handled by the same exact pipeline used elsewhere in `linalg` (rationalise → solve → numericalise) rather than a floating-point LU.

- `Protected`.
- The matrix `m` may be square or rectangular.
- The argument `b` may be a vector (in which case the result is a
  vector) or a matrix (in which case the result is a matrix whose
  `k`-th column solves `m . x == b[[All, k]]`).
- **FLINT acceleration** (when built with FLINT): a square, nonsingular,
  all-integer/rational system is solved exactly via `fmpq_mat_solve` in
  polynomial time; the unique solution matches the classical division-free
  result. Non-square, singular, or symbolic systems fall through to the
  classical solver (which handles rectangular / underdetermined / inconsistent
  cases). Exposed directly as `` FLINT`LinearSolve `` (see *Structural
  Manipulation*).
- Higher-rank tensor inputs are supported. A rank-N `m` with
  dimensions `{d_1, ..., d_{N-1}, n}` is interpreted as a
  `(d_1 * ... * d_{N-1}) x n` linear system whose leading dimensions
  combine into rows; `b` then has dimensions
  `{d_1, ..., d_{N-1}, e_1, ..., e_p}` and the result has dimensions
  `{n, e_1, ..., e_p}`. When `p == 0` the result is a flat vector of
  length `n`.
- For under-determined systems LinearSolve returns one particular
  solution, with every free (non-pivot) variable set to 0; `Solve`
  returns the general solution.
- Issues `LinearSolve::nosol` and returns unevaluated when no solution
  exists.
- Issues `LinearSolve::matrix` / `::lvec` / `::lvec1` and returns
  unevaluated for shape errors.
- Default method is fraction-free Gauss-Jordan elimination on the
  augmented matrix `[m | b]` (the same Bareiss-like routine used by
  `RowReduce` and `Inverse`), so exact integer, rational, and
  symbolic inputs flow through with no spurious denominator blow-up.
- Lives in `src/linalg/linsolve.c`.
- Accepts an optional `Method -> "<name>"` argument:
  - `Method -> Automatic` or `Method -> "Automatic"` — default (alias for `"DivisionFreeRowReduction"`).
  - `Method -> "DivisionFreeRowReduction"` — Bareiss-like fraction-free Gauss-Jordan on `[m | b]`. Recommended for exact / symbolic inputs.
  - `Method -> "OneStepRowReduction"` — classical Gauss-Jordan with one division per pivot per element on `[m | b]`. Each per-cell update is canonicalised via `Together` so symbolic rationals reduce. Fast on numeric matrices.
  - `Method -> "CofactorExpansion"` — Cramer's rule. Requires square non-singular `m`; emits `LinearSolve::cofnsq` on a non-square `m`, and `LinearSolve::cofsng` on a structurally singular `m`. For matrix `b` the rule is applied column-by-column.
- Unknown method names emit `LinearSolve::method` and the call remains unevaluated.

**Attributes:** `Protected`.

## See also

[Solve](../../solutions-of-equations/Solve/), [RowReduce](../../linear-algebra/RowReduce/), [Inverse](../../linear-algebra/Inverse/), [Together](../../algebra/Together/)

## References

- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — solving linear systems by Gaussian elimination.
- L. N. Trefethen and D. Bau III, *Numerical Linear Algebra*, SIAM, 1997 — LU factorisation and triangular solves.
- Erwin H. Bareiss, "Sylvester's Identity and Multistep Integer-Preserving Gaussian Elimination", Mathematics of Computation 22 (1968).
- Source: [`src/linalg/linsolve.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/linsolve.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_compile_linalg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_linalg.c)
- Tests: [`tests/test_flint_bridge.c`](https://github.com/stblake/mathilda/blob/main/tests/test_flint_bridge.c)
- Tests: [`tests/test_lapack_builtin.c`](https://github.com/stblake/mathilda/blob/main/tests/test_lapack_builtin.c)
- Tests: [`tests/test_linearsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_linearsolve.c)

## Notes & additional examples

### Notes

`LinearSolve[m, b]` returns an `x` satisfying `m . x == b`; over exact integer inputs the solution is exact rational, as in the first example. The default algorithm is fraction-free (Bareiss-like) Gauss-Jordan elimination on the augmented matrix `[m | b]`, so no spurious denominators appear during elimination. When the system is inconsistent (last example) `LinearSolve::nosol` is emitted on stderr and the call is returned unevaluated; for under-determined systems a particular solution with free variables set to 0 is returned. The right-hand side may also be a matrix (one column per system), giving a matrix of solutions.
