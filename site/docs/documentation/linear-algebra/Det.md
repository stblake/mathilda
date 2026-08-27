# Det

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Det[m]`**

gives the determinant of the square matrix m.

<details>
<summary>Notes</summary>

Exact integer / rational / symbolic inputs use Bareiss-style fraction-free Gaussian elimination; machine-precision Real / Complex inputs dispatch to LAPACK LU (dgetrf / zgetrf) and accumulate the pivot-signed product of diagonal entries; arbitrary-precision MPFR inputs run a Doolittle LU at the input precision.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Det[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= 0

In[2]:= Det[{{1.7, 7.1, -2.7}, {2.2, 8.7, 3.2}, {3.2, -9.2, 1.2}}]
Out[2]= 251.572

In[3]:= Det[{{a, b, c}, {d, e, f}, {g, h, i}}]
Out[3]= c (-e g + d h) - b (-f g + d i) + a (-f h + e i)
```

### Applications (7)

```mathematica
In[4]:= Det[{{1, 2}, {3, 4}}]
Out[4]= -2

In[5]:= Det[{{2, 0, 0}, {0, 3, 0}, {0, 0, 4}}]
Out[5]= 24

In[6]:= Det[{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}}]
Out[6]= -3

In[7]:= Det[{{a, b}, {c, d}}]
Out[7]= -b c + a d

In[8]:= Det[{{1, 1, 1}, {a, b, c}, {a^2, b^2, c^2}}]
Out[8]= -a^2 b + a b^2 + a^2 c - b^2 c - a c^2 + b c^2

In[9]:= Det[Table[1/(i + j - 1), {i, 4}, {j, 4}]]
Out[9]= 1/6048000

In[10]:= Det[{{N[Pi, 40], 1}, {1, N[E, 40]}}]
Out[10]= 7.5397342226735670654635508695465744950351
```

## Implementation notes

**Algorithm.** `builtin_det` validates that the argument is a non-empty square rank-2 tensor (via `get_tensor_dims`; otherwise it emits `Det::matsq` and returns `NULL`), flattens it row-major into an `Expr**`, and computes the determinant by full **Laplace cofactor expansion** along the first row, recursively (`laplace_det`). Each cofactor term is built as `Times[±1, element, minor]` and accumulated with `Plus`, every product/sum being reduced through `eval_and_free`, so cancellation and symbolic simplification happen as the expansion unwinds. This keeps results exact and symbolic for integer, rational, and symbolic matrices.

**Data structures.** The matrix is a flat `Expr**` of `n*n` element pointers; recursion carries an explicit `int* cols` index set and a fixed `row` cursor, deleting one column per level rather than copying submatrices. `laplace_det` is exported via `linalg.h` and reused by `Cross`.

**Complexity / limits.** Cofactor expansion is `O(n!)`, so it is only practical for small `n`. There is no fraction-free Bareiss or LU fast path in this handler — the larger fraction-free Gauss-Jordan machinery lives in `inv.c`/`linsolve.c` for inversion and solving, not in `Det`.

- `Protected`.
- Evaluates the determinant of a square matrix symbolically or numerically using Laplace expansion.
- Returns a warning `Det::matsq` if `m` is not a non-empty square matrix.
- **FLINT acceleration** (when built with FLINT): a matrix whose entries are all integer or rational is computed exactly via `fmpq_mat_det` in polynomial time, avoiding the `O(n!)` Laplace expansion (e.g. a 12×12 Hilbert determinant is instant and exact). Symbolic matrices fall through to Laplace. The same kernel is exposed directly as `` FLINT`Det `` (see the FLINT` context section in *Structural Manipulation*).
- **Packed/NDArray fast path** (machine reals): the LU factorisation runs through LAPACK `dgetrf` / `zgetrf` (in-house partial-pivot `double` LU as fallback).
- **Overflow/underflow → arbitrary precision**: the determinant is the product of the LU pivots, which can exceed the IEEE-double range even when the matrix is ordinary — a 200×200 `RandomReal[{-10,10}]` matrix has `|Det| ≈ 10^340`. Rather than return `inf` (or `0` from a mid-product underflow), the pivot product is re-accumulated in a 53-bit-mantissa MPFR value whose exponent range is effectively unbounded, so the answer is finite and correct to machine precision (`Det[RandomReal[{-10,10},{200,200}]]` returns a `≈ -1.08×10^340` real). A genuinely singular matrix still returns `0`.
- **Arbitrary-precision matrices**: a genuine MPFR matrix (precision > machine) uses an `O(n^3)` MPFR LU determinant (`mpfr_det_dispatch`) rather than the `O(n!)` Laplace expansion, which previously hung for `n ≳ 12`.

**Attributes:** `Protected`.

## References

- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — Gaussian elimination and the LU view of the determinant.
- Source: [`src/linalg/det.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/det.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_characteristicpolynomial.c`](https://github.com/stblake/mathilda/blob/main/tests/test_characteristicpolynomial.c)
- Tests: [`tests/test_diagonal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_diagonal.c)
- Tests: [`tests/test_eigen.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eigen.c)
- Tests: [`tests/test_flint_bridge.c`](https://github.com/stblake/mathilda/blob/main/tests/test_flint_bridge.c)

## Notes & additional examples

### Notes

For a diagonal or triangular matrix the determinant is simply the product of the diagonal entries, as the third and fourth examples illustrate. Exact integer, rational, and symbolic inputs are handled by Bareiss-style fraction-free Gaussian elimination, so intermediate results never introduce spurious denominators and symbolic determinants come back fully factored where possible. A singular matrix returns 0 exactly. The argument must be a square matrix.
