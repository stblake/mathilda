# Inverse

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Inverse[m]
    gives the inverse of a square matrix m.
Inverse[m, Method -> "<name>"]
    runs a specific inversion algorithm.

Inverse works on both symbolic and numerical matrices.
For matrices with approximate real or complex numbers, the
inverse is generated to the maximum possible precision given the
input.  Inverse::sing is issued for singular matrices and
Inverse::matsq for non-square / empty input; in either case the
call is returned unevaluated.

Accepted method names:
  "Automatic"                 — alias for "DivisionFreeRowReduction" (default)
  "DivisionFreeRowReduction"  — Bareiss-like fraction-free Gauss-Jordan on [m | I]
  "OneStepRowReduction"       — classical Gauss-Jordan with division per pivot
  "CofactorExpansion"         — adjugate / determinant formula via Laplace expansion

An unknown method name emits Inverse::method and leaves the call
unevaluated.  Method -> Automatic (the symbol) is also accepted.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Inverse[{{1.4,2},{3,-6.7}}]
Out[1]= {{0.435631, 0.130039}, {0.195059, -0.0910273}}

In[2]:= Inverse[{{1,2,3},{4,2,2},{5,1,7}}]
Out[2]= {{-2/7, 11/42, 1/21}, {3/7, 4/21, -5/21}, {1/7, -3/14, 1/7}}

In[3]:= Inverse[{{u,v},{v,u}}]
Out[3]= {{u/(u^2 - v^2), -v/(u^2 - v^2)}, {-v/(u^2 - v^2), u/(u^2 - v^2)}}

In[4]:= Inverse[{{1.2,2.5,-3.2},{0.7,-9.4,5.8},{-0.2,0.3,6.4}}]
Out[4]= {{0.74546, 0.204249, 0.187629}, {0.0679223, -0.0847825, 0.110795}, {0.0201118, 0.010357, 0.15692}}

In[5]:= Inverse[{{2,3,2},{4,9,2},{7,2,4}}]
Out[5]= {{-8/13, 2/13, 3/13}, {1/26, 3/26, -1/13}, {55/52, -17/52, -3/26}}

In[6]:= Inverse[{{a,b},{c,d}}]
Out[6]= {{d/(-b c + a d), -b/(-b c + a d)}, {-c/(-b c + a d), a/(-b c + a d)}}

In[7]:= Inverse[{{1,2},{1,2}}]
Out[7]= Inverse[{{1, 2}, {1, 2}}]

In[8]:= a = {{1,2},{3,4}}; a . Inverse[a] == IdentityMatrix[2]
Out[8]= True
```

## Implementation notes

**Algorithm.** `builtin_inverse` validates that the argument is a non-empty square matrix, parses an optional `Method`, and dispatches (`MatsolMethod`) to one of three exact workers:

- `inverse_divfree` (default, `Method -> Automatic` / `"DivisionFreeRowReduction"`): Bareiss-like **fraction-free Gauss-Jordan** elimination on the augmented matrix `[A | I]`. A running pivot product `P` divides each updated entry so the elimination stays division-free (exact, no rational/GCD blow-up); the right half becomes `A^{-1}` once the left half is reduced to (a scalar multiple of) the identity. Singular matrices emit `Inverse::sing` and return `NULL`.
- `inverse_onestep` (`"OneStepRowReduction"`): classical Gauss-Jordan with one division per pivot, each entry canonicalised via `Together` so symbolic cancellations are still detected.
- `inverse_cofactor` (`"CofactorExpansion"`): the adjugate/determinant formula `A^{-1}[i,j] = (-1)^{i+j} det(M_{j,i})/det(A)`, with each minor computed by the same Laplace expansion used by `Det` — `O(n!·n^2)`, for tiny `n` only.

Inexact (`Real`/MPFR) matrices are routed through the standard `common_scan_inexact` → `common_rationalize_input` pipeline: the matrix is rationalised at the minimum precision present, inverted exactly, then numericalised back to that precision, so the rank/pivot decisions are exact.

**Data structures.** Dense flat `Expr**` augmented matrix of `n × 2n` element pointers, row-major; the pivot product is a single shared `Expr*`. The same `inv.c` module also implements `PseudoInverse` via a full-rank `B·C` decomposition from `RowReduce`.

**Complexity / limits.** Fraction-free Gauss-Jordan is `O(n^3)` arithmetic ops with controlled intermediate growth; cofactor expansion is `O(n!)`. There is no machine-precision LU `dgetrf`-style kernel — machine matrices are handled by exact rationalisation rather than floating-point LU.

- `Protected`.
- Works on both symbolic and numerical matrices.
- For matrices with approximate real or complex numbers, the inverse is generated to the maximum possible precision given the input.
- Issues `Inverse::sing` warning and returns unevaluated if the matrix is singular.
- Issues `Inverse::matsq` warning and returns unevaluated if the argument is not a non-empty square matrix.
- **FLINT acceleration** (when built with FLINT): an all-integer/rational matrix is inverted exactly via `fmpq_mat_inv` in polynomial time (default/division-free method); the inverse is unique so the result matches the classical path. Singular or symbolic matrices fall through unchanged. Exposed directly as `` FLINT`Inverse `` (see the FLINT` context section in *Structural Manipulation*).
- Satisfies the relation `a . Inverse[a] == Inverse[a] . a == IdentityMatrix[n]`.
- Satisfies the relation `Inverse[a . b] == Inverse[b] . Inverse[a]`.
- Accepts an optional `Method -> "<name>"` argument that selects the inversion algorithm. Shares the same method-name grammar as `RowReduce` and `LinearSolve`.
  - `Method -> Automatic` or `Method -> "Automatic"` (default) — alias for `"DivisionFreeRowReduction"`.
  - `Method -> "DivisionFreeRowReduction"` — Bareiss-like fraction-free Gauss-Jordan elimination on the augmented matrix `[A | I]`. Best choice for exact integer / rational / symbolic input — never produces a denominator larger than necessary.
  - `Method -> "OneStepRowReduction"` — classical Gauss-Jordan on `[A | I]` with one division per pivot per row entry. Each entry is canonicalised via `Together` so symbolic cancellations are still detected. Fast on numeric matrices.
  - `Method -> "CofactorExpansion"` — adjugate / determinant formula `A^-1[i,j] = (-1)^(i+j) Det[M[j,i]] / Det[A]`, with `Det` computed via Laplace cofactor expansion. Time complexity is `O(n! n^2)`; intended for small `n` or closed-form symbolic inverses of small matrices. Emits `Inverse::sing` if `Det[A]` is structurally zero.
- Unknown method names emit `Inverse::method` and the call remains unevaluated.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — Gauss-Jordan elimination and matrix inversion.
- Erwin H. Bareiss, "Sylvester's Identity and Multistep Integer-Preserving Gaussian Elimination", Mathematics of Computation 22 (1968).
- Source: [`src/linalg/inv.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/inv.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Inverse[{{2, 0}, {0, 4}}]
Out[1]= {{1/2, 0}, {0, 1/4}}
```

```mathematica
In[1]:= Inverse[{{1, 2}, {3, 4}}]
Out[1]= {{-2, 1}, {3/2, -1/2}}
```

```mathematica
In[1]:= Inverse[{{1, 1, 1}, {0, 1, 1}, {0, 0, 1}}]
Out[1]= {{1, -1, 0}, {0, 1, -1}, {0, 0, 1}}
```

```mathematica
In[1]:= Inverse[{{a, b}, {c, d}}]
Out[1]= {{d/(-b c + a d), -b/(-b c + a d)}, {-c/(-b c + a d), a/(-b c + a d)}}
```

```mathematica
In[1]:= Inverse[{{1, x}, {x, 1}}]
Out[1]= {{1/(1 - x^2), -x/(1 - x^2)}, {-x/(1 - x^2), 1/(1 - x^2)}}
```

```mathematica
In[1]:= Inverse[{{2.0, 1.0}, {1.0, 3.0}}]
Out[1]= {{0.6, -0.2}, {-0.2, 0.4}}
```

```mathematica
In[1]:= Inverse[{{1, 2}, {2, 4}}]
Out[1]= Inverse[{{1, 2}, {2, 4}}]
```

### Notes

Exact integer inputs yield exact rational inverses via fraction-free Gauss-Jordan elimination on the augmented matrix `[m | I]`. A singular matrix (last example) emits the `Inverse::sing` diagnostic on stderr and returns the call unevaluated; a non-square or empty argument emits `Inverse::matsq`. The `Method` option selects among `"DivisionFreeRowReduction"` (the default), `"OneStepRowReduction"`, and `"CofactorExpansion"`; an unrecognised method emits `Inverse::method`.
