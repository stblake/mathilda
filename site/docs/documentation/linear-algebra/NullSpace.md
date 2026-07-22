# NullSpace

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
NullSpace[m]
    gives a list of vectors that forms a basis for the null
    space of the matrix m (i.e. vectors v such that m . v == 0).
NullSpace[m, Method -> "<name>"]
    runs a specific elimination algorithm.  Accepted method
    names are the same as RowReduce / LinearSolve / Inverse:
      "Automatic"                 — alias for "DivisionFreeRowReduction" (default)
      "DivisionFreeRowReduction"  — Bareiss-like fraction-free Gauss-Jordan
      "OneStepRowReduction"       — classical Gauss-Jordan with division per pivot
      "CofactorExpansion"         — identity-if-invertible (falls back to
                                     DivisionFreeRowReduction on singular /
                                     rectangular m)

NullSpace works on both numerical and symbolic matrices.  The
matrix m may be square or rectangular.  When m has full column
rank the result is the empty list `{}`.

Basis vectors are returned with the rightmost free column
first.  For exact integer / rational input each basis vector
is scaled to clear integer denominators, so the result is
integer-valued whenever the input is integer-valued.  For
symbolic input the basis vectors are left in their natural
rational form.

An unknown method name emits NullSpace::method and leaves the
call unevaluated.  A non-rank-2 / empty matrix emits
NullSpace::matrix and the call is left unevaluated.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= NullSpace[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= {{1, -2, 1}}

In[2]:= NullSpace[{{a, b}, {2 a, 2 b}}]
Out[2]= {{-b/a, 1}}

In[3]:= NullSpace[{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}}]
Out[3]= {}

In[4]:= NullSpace[{{a, b, c}, {c, b, a}, {0, 0, 0}}]
Out[4]= {{1, -(a/b + c/b), 1}}

In[5]:= NullSpace[{{3, 2, 2, 4}, {2, 3, -2, 7}, {3, 2, 5, 7}}]
Out[5]= {{12, -23, -5, 5}}

In[6]:= NullSpace[IdentityMatrix[5]]
Out[6]= {}

In[7]:= m = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}; m . First[NullSpace[m]]
Out[7]= {0, 0, 0}
```

## Implementation notes

**Algorithm.** `builtin_nullspace` computes a basis for `{v : m.v == 0}` by reduction to RREF rather than by an orthogonal-completion (QR/SVD) route. `nullspace_core` calls `RowReduce[m, Method -> ...]` through the evaluator (`call_rowreduce`, forwarding the optional `Method` option), flattens the RREF, and locates each row's pivot column as its leftmost structurally non-zero entry (`is_zero_poly`). For every *free* (non-pivot) column `f`, iterated rightmost-to-leftmost to match the standard ordering, it builds a length-`cols` basis vector with `v[f] = 1`, `v[p] = -RREF[row_of_p, f]` for each pivot column `p`, and `0` elsewhere.

**Data structures.** A flat `Expr**` of the RREF plus an `int* pivot_for_col` map; basis vectors are accumulated into a growable `Expr**`. For exact rational input each vector is scaled by the LCM of its entries' integer denominators (`clear_int_denominators`, via GMP `mpz_lcm` and the `Denominator` builtin) so an integer-valued nullspace comes out integer-valued; symbolic/inexact entries are left in natural form. Full-column-rank input returns `List[]`.

**Limits.** Inherits whatever the `RowReduce` dispatcher (default `DivisionFreeRowReduction`, Bareiss-like fraction-free Gauss-Jordan) can reduce; correctness for symbolic matrices depends on `is_zero_poly`/`Together` detecting cancellations. Bad `Method` values emit `NullSpace::method`; non-rectangular input emits `NullSpace::matrix`.

- `Protected`.
- Returns a list of linearly-independent vectors whose span equals

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/linalg/nullspace.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/nullspace.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

The classic rank-deficient magic-like matrix has a one-dimensional kernel:

```mathematica
In[1]:= NullSpace[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= {{1, -2, 1}}
```

A full-column-rank matrix has only the trivial null space, returned as the empty
list:

```mathematica
In[1]:= NullSpace[{{1, 0}, {0, 1}}]
Out[1]= {}
```

A wide matrix with a two-dimensional kernel — the basis vectors are scaled to
clear denominators and ordered with the rightmost free column first:

```mathematica
In[1]:= NullSpace[{{1, 2, 3, 4}, {2, 4, 6, 8}}]
Out[1]= {{-4, 0, 0, 1}, {-3, 0, 1, 0}, {-2, 1, 0, 0}}
```

NullSpace works symbolically — here the kernel is parametrized by `a`:

```mathematica
In[1]:= NullSpace[{{1, a}, {1, a}}]
Out[1]= {{-a, 1}}
```

The defining property `m . v == 0` can be checked directly:

```mathematica
In[1]:= m = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}; m . First[NullSpace[m]]
Out[1]= {0, 0, 0}
```

### Notes

`NullSpace[m]` returns a basis for the null space of `m` (the vectors `v` with
`m . v == 0`). The matrix may be square or rectangular; full column rank yields
the empty list `{}`. For exact integer or rational input each basis vector is
scaled to clear denominators, so the result stays integer-valued whenever the
input is; symbolic input keeps its natural rational form. A `Method` option
selects the elimination algorithm (default `"DivisionFreeRowReduction"`, a
Bareiss-like fraction-free Gauss-Jordan).
