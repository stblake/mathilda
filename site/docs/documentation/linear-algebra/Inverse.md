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

- `Protected`.
- Works on both symbolic and numerical matrices.
- For matrices with approximate real or complex numbers, the inverse is generated to the maximum possible precision given the input.
- Issues `Inverse::sing` warning and returns unevaluated if the matrix is singular.
- Issues `Inverse::matsq` warning and returns unevaluated if the argument is not a non-empty square matrix.
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
- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
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
In[1]:= Inverse[{{1, 2}, {2, 4}}]
Out[1]= Inverse[{{1, 2}, {2, 4}}]
```

### Notes

Exact integer inputs yield exact rational inverses via fraction-free Gauss-Jordan elimination on the augmented matrix `[m | I]`. A singular matrix (last example) emits the `Inverse::sing` diagnostic on stderr and returns the call unevaluated; a non-square or empty argument emits `Inverse::matsq`. The `Method` option selects among `"DivisionFreeRowReduction"` (the default), `"OneStepRowReduction"`, and `"CofactorExpansion"`; an unrecognised method emits `Inverse::method`.
