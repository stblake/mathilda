# HermiteDecomposition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HermiteDecomposition[m]`**

Gives the Hermite normal form decomposition {u, r} of the integer matrix m: u is unimodular (Abs\[Det\[u\]\] == 1), r is the row Hermite normal form, and u . m == r.  r is in echelon shape with positive pivots and entries above each pivot reduced into \[0, pivot).

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= HermiteDecomposition[{{2, 3, 1}, {4, 1, 5}, {6, 2, 0}}]
Out[1]= {{{0, 2, -1}, {1, 4, -3}, {1, 7, -5}}, {{2, 0, 10}, {0, 1, 21}, {0, 0, 36}}}

In[2]:= h = HermiteDecomposition[{{1, 2, 3}, {4, 5, 6}}]; h[[1]] . {{1, 2, 3}, {4, 5, 6}} == h[[2]]
Out[2]= True
```

### Applications (4)

Gives {u, r} with u unimodular and u . m == r

```mathematica
In[3]:= HermiteDecomposition[{{1, 2}, {3, 4}}]
Out[3]= {{{-2, 1}, {3, -1}}, {{1, 0}, {0, 2}}}
```

A non-square 2x3 integer matrix

```mathematica
In[4]:= HermiteDecomposition[{{2, 3, 4}, {5, 6, 7}}]
Out[4]= {{{-2, 1}, {5, -2}}, {{1, 0, -1}, {0, 3, 6}}}
```

Already in Hermite form, so u is the identity

```mathematica
In[5]:= HermiteDecomposition[{{2, 0}, {0, 3}}]
Out[5]= {{{1, 0}, {0, 1}}, {{2, 0}, {0, 3}}}
```

The transform u has determinant +-1

```mathematica
In[6]:= Det[HermiteDecomposition[{{1, 2}, {3, 4}}][[1]]]
Out[6]= -1
```

## Algorithm

hnf.c -- Hermite Normal Form over Z, and the HermiteDecomposition builtin.

```text
`linalg_hnf` computes a unimodular P and row-HNF R with P*A == R, tracking
every integer row operation on P.  The elimination in each column uses the
```

extended-gcd 2x2 unimodular transform

```text
    [ s   t ] [row_r]      [ g*... ]           s*a_r + t*a_i = g
    [-b   a ] [row_i]  ->  [   0   ]  in col c, a = a_r/g, b = a_i/g,
```

whose determinant is s*a + t*b = (s*a_r + t*a_i)/g = 1, so P stays

```text
unimodular.  Pivots are then made positive and entries above each pivot are
reduced into [0, pivot).  This is the reusable integer primitive behind
```

exact linear Diophantine system solving (src/solve/solveint_linear.c).

## Implementation notes

**Algorithm.** `builtin_hermite_decomposition` gives `{u, r}` for an integer matrix `m`, where `u` is unimodular (`Abs[Det[u]] == 1`), `r` is the row Hermite normal form, and `u . m == r`. The core `linalg_hnf` works column by column over `Z`: it finds a pivot at or below the current pivot row, swaps it up, then eliminates every entry below it with the extended-gcd 2×2 unimodular transform `[[s, t], [-b, a]]`, where `mpz_gcdext` gives `s·a_r + t·a_i = g` and `a = a_r/g`, `b = a_i/g`. That transform has determinant `s·a + t·b = 1`, so the tracked `P` stays unimodular. Each pivot is then normalised positive and the entries above it are reduced into `[0, pivot)` by floor division. `P` starts as the identity and records every row operation; `R` starts as a copy of `A`; the routine returns the rank.

**Data structures.** The matrices are flat row-major `mpz_t` arrays (GMP), allocated by `hnf_mat_alloc` and released by `linalg_hnf_free`. The input `Expr` matrix is `flatten_tensor`'d into an `Expr**` and each cell coerced with `expr_to_mpz` after an `expr_is_integer_like` check; a non-integer entry aborts with the `HermiteDecomposition::intm` message. The result is rebuilt by `hnf_mat_to_expr` into a nested `List` of `expr_bigint_normalize`'d integers — `u` is `m×m`, `r` is `m×n`. A packed/`NDArray` argument is delisted via `linalg_delist_and_reeval` first (the head is on `src/pack.c`'s `AWARE` list for that reason), though the entries must still be exact integers.

**Complexity / limits.** `O(m·n·rank)` GMP row operations plus the coefficient growth characteristic of HNF over `Z`. The decomposition is defined only for a non-empty rectangular matrix with integer entries: a wrong argument count emits `HermiteDecomposition::argx`, a non-matrix `HermiteDecomposition::matrix`, and rational/real/symbolic entries `HermiteDecomposition::intm` — each routed through `mth_message`, each leaving the call unevaluated. The same `linalg_hnf` primitive backs exact linear Diophantine system solving (`src/solve/solveint_linear.c`).

**Attributes:** `Protected`.

## References

- H. Cohen, *A Course in Computational Algebraic Number Theory* (Springer, 1993), §2.4 — Hermite Normal Forms.
- Source: [`src/linalg/hnf.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/hnf.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_latticereduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_latticereduce.c)

## Notes & additional examples

### Notes

`HermiteDecomposition[m]` returns `{u, r}` where `u` is unimodular
(`Abs[Det[u]] == 1`), `r` is the row Hermite normal form, and `u . m == r`. The
form `r` is in echelon shape with positive pivots and every entry above a pivot
reduced into `[0, pivot)`.

The decomposition is defined only over the integers. A non-integer matrix is left
unevaluated with a message, as is a non-rectangular or empty argument. The same
integer HNF primitive drives Mathilda's exact linear Diophantine solving.
