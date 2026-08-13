# NullSpace

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NullSpace[m]`**

gives a list of vectors that forms a basis for the null space of the matrix m (i.e. vectors v such that m . v == 0).

**`NullSpace[m, Method -> "<name>"]`**

runs a specific elimination algorithm.  Accepted method names are the same as RowReduce / LinearSolve / Inverse: "Automatic"                 — alias for "DivisionFreeRowReduction" (default) "DivisionFreeRowReduction"  — Bareiss-like fraction-free Gauss-Jordan "OneStepRowReduction"       — classical Gauss-Jordan with division per pivot "CofactorExpansion"         — identity-if-invertible (falls back to DivisionFreeRowReduction on singular / rectangular m)

<details>
<summary>Notes</summary>

NullSpace works on both numerical and symbolic matrices.  The matrix m may be square or rectangular.  When m has full column rank the result is the empty list \`{}\`. Basis vectors are returned with the rightmost free column first.  For exact integer / rational input each basis vector is scaled to clear integer denominators, so the result is integer-valued whenever the input is integer-valued.  For symbolic input the basis vectors are left in their natural rational form. An unknown method name emits NullSpace::method and leaves the call unevaluated.  A non-rank-2 / empty matrix emits NullSpace::matrix and the call is left unevaluated.

</details>

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

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

### Applications (5)

```mathematica
In[8]:= NullSpace[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[8]= {{1, -2, 1}}

In[9]:= NullSpace[{{1, 0}, {0, 1}}]
Out[9]= {}

In[10]:= NullSpace[{{1, 2, 3, 4}, {2, 4, 6, 8}}]
Out[10]= {{-4, 0, 0, 1}, {-3, 0, 1, 0}, {-2, 1, 0, 0}}

In[11]:= NullSpace[{{1, a}, {1, a}}]
Out[11]= {{-a, 1}}

In[12]:= m = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}; m . First[NullSpace[m]]
Out[12]= {0, 0, 0}
```

## Options & behaviour

> Implementation lives in `src/linalg/nullspace.c` (registered by
> `matnull_init`).  The basis extraction calls back into RowReduce via
> the evaluator, so any future improvement to RowReduce automatically
> propagates to NullSpace.

## Algorithm

nullspace.c

```text
NullSpace[m]                       -- list of basis vectors v such that
                                      m . v == 0.
NullSpace[m, Method -> "<name>"]   -- explicit RREF method dispatch.
```

Implementation strategy:

```text
  1. Validate that `m` is a rank-2 non-empty matrix.  Wrong shape
     returns NULL so the call remains unevaluated (matching the
     RowReduce / LinearSolve convention).

  2. Reduce `m` to reduced row echelon form (RREF) by calling
     RowReduce as a builtin.  Routing through the existing
     `builtin_rowreduce` dispatcher avoids duplicating the
     Bareiss-like / OneStep / Cofactor algorithms and means
     NullSpace inherits any future RREF improvements automatically.

  3. Identify pivot columns from the RREF: the leftmost non-zero
     entry in each row marks one pivot column.

  4. For each free (non-pivot) column f, build a basis vector v of
     length cols by setting:
       - v[f] = 1
       - v[p] = -RREF[row_of_p, f]  for each pivot column p
       - v[other free] = 0
     Iterate `f` from the rightmost free column to the leftmost so
     the basis ordering matches Mathematica.

  5. For exact rational input, scale each vector by the LCM of its
     entries' integer denominators so the result is integer-valued
     whenever the input is integer-valued.  For symbolic / inexact
     input the vectors are left in their natural form.

Memory ownership: standard builtin contract.  On success this file
```

owns `res` and frees it; on NULL return the caller retains ownership of `res`.

## Implementation notes

**Algorithm.** `builtin_nullspace` computes a basis for `{v : m.v == 0}` by reduction to RREF rather than by an orthogonal-completion (QR/SVD) route. `nullspace_core` calls `RowReduce[m, Method -> ...]` through the evaluator (`call_rowreduce`, forwarding the optional `Method` option), flattens the RREF, and locates each row's pivot column as its leftmost structurally non-zero entry (`is_zero_poly`). For every *free* (non-pivot) column `f`, iterated rightmost-to-leftmost to match the standard ordering, it builds a length-`cols` basis vector with `v[f] = 1`, `v[p] = -RREF[row_of_p, f]` for each pivot column `p`, and `0` elsewhere.

**Data structures.** A flat `Expr**` of the RREF plus an `int* pivot_for_col` map; basis vectors are accumulated into a growable `Expr**`. For exact rational input each vector is scaled by the LCM of its entries' integer denominators (`clear_int_denominators`, via GMP `mpz_lcm` and the `Denominator` builtin) so an integer-valued nullspace comes out integer-valued; symbolic/inexact entries are left in natural form. Full-column-rank input returns `List[]`.

**Limits.** Inherits whatever the `RowReduce` dispatcher (default `DivisionFreeRowReduction`, Bareiss-like fraction-free Gauss-Jordan) can reduce; correctness for symbolic matrices depends on `is_zero_poly`/`Together` detecting cancellations. Bad `Method` values emit `NullSpace::method`; non-rectangular input emits `NullSpace::matrix`.

- `Protected`.
- Returns a list of linearly-independent vectors whose span equals
  `{v : m . v == 0}`. If `m` has full column rank the result is the
  empty list `{}`.
- Works on numerical (Integer / Rational / Real / MPFR / Complex),
  big-integer, and symbolic matrices.
- The matrix `m` may be square or rectangular.
- Basis vectors are returned with the **rightmost free column first**,
  matching the standard ordering.
- For exact integer / rational input each basis vector is scaled to
  clear integer denominators, so the result is integer-valued whenever
  the input is integer-valued. For symbolic input the basis vectors
  are left in their natural rational form.
- Internally calls `RowReduce[m, Method -> "<name>"]` and extracts the
  basis from the resulting RREF. The `Method` option therefore shares
  the same grammar as `RowReduce` / `LinearSolve` / `Inverse`:
  - `Method -> Automatic` (default) — alias for
    `"DivisionFreeRowReduction"`.
  - `Method -> "DivisionFreeRowReduction"` — Bareiss-like fraction-free
    Gauss-Jordan; best for exact integer / rational / symbolic input.
  - `Method -> "OneStepRowReduction"` — classical Gauss-Jordan with one
    division per pivot per element. Each entry is canonicalised via
    `Together`.
  - `Method -> "CofactorExpansion"` — identity-if-invertible path
    inside RowReduce; falls back to `"DivisionFreeRowReduction"` on
    singular / rectangular input.
- Issues `NullSpace::matrix` and returns unevaluated if the argument
  is not a non-empty rank-2 tensor.
- Issues `NullSpace::method` and returns unevaluated for unknown
  method names.

**Attributes:** `Protected`.

## References

**See also:** [RowReduce](../../linear-algebra/RowReduce/), [LinearSolve](../../linear-algebra/LinearSolve/), [Inverse](../../linear-algebra/Inverse/), [Together](../../algebra/Together/)

- Source: [`src/linalg/nullspace.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/nullspace.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_nullspace.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nullspace.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)

## Notes & additional examples

### Notes

`NullSpace[m]` returns a basis for the null space of `m` (the vectors `v` with
`m . v == 0`). The matrix may be square or rectangular; full column rank yields
the empty list `{}`. For exact integer or rational input each basis vector is
scaled to clear denominators, so the result stays integer-valued whenever the
input is; symbolic input keeps its natural rational form. A `Method` option
selects the elimination algorithm (default `"DivisionFreeRowReduction"`, a
Bareiss-like fraction-free Gauss-Jordan).
