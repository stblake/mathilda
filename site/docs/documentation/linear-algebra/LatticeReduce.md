# LatticeReduce

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LatticeReduce[m]`**

gives an LLL-reduced basis for the lattice spanned by the rows (vectors) of m.  The entries of m may be integers, Gaussian integers, rationals, or Gaussian rationals.  Reduction is exact (GMP rational arithmetic, so it is correct for both machine-size and arbitrary-precision entries) and preserves the lattice, its determinant, and every linear relation among the rows.  The rows must be linearly independent.

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= LatticeReduce[{{1, 0, 0, 1345}, {0, 1, 0, 35}, {0, 0, 1, 154}}]
Out[1]= {{0, 9, -2, 7}, {1, 1, -9, -6}, {1, -3, -8, 8}}

In[2]:= {w1, w2} = LatticeReduce[{{12, 2}, {13, 4}}]
Out[2]= {{1, 2}, {9, -4}}

In[3]:= a = {{1, 0, 0, -1}, {0, 1, 0, -2}, {0, 0, 1, -3}}; b = LatticeReduce[a]
Out[3]= {{1, 0, 0, -1}, {-1, 1, 0, -1}, {-1, -1, 1, 0}}
```

Relations preserved

```mathematica
In[4]:= b . {1, 2, 3, 1}
Out[4]= {0, 0, 0}
```

### Applications (5)

```mathematica
In[5]:= LatticeReduce[{{1, 1, 1}, {-1, 0, 2}, {3, 5, 6}}]
Out[5]= {{0, 1, 0}, {1, 0, 1}, {-2, 0, 1}}

In[6]:= Det[LatticeReduce[{{201, 37}, {1648, 297}}]]
Out[6]= -1279

In[7]:= Det[{{201, 37}, {1648, 297}}]
Out[7]= -1279

In[8]:= LatticeReduce[{{1, 0, 0, 31415927}, {0, 1, 0, 27182818}, {0, 0, 1, 16180340}}]
Out[8]= {{61, -183, 189, 113}, {-198, 108, 203, -182}, {-235, 146, 211, 323}}

In[9]:= LatticeReduce[{{1/2, 1}, {1, 1/3}}]
Out[9]= {{1/2, -2/3}, {1, 1/3}}
```

## Options & behaviour

> Implementation lives in `src/linalg/latticereduce.c` (registered by
> `linalg_init`). Each Gaussian-rational scalar is a pair of GMP `mpq_t`
> (`GRat`); the Gram–Schmidt data (`μ`, `|b*|²`) is computed once and
> maintained incrementally — updated in place on size reduction and on a
> Lovász swap via the conjugate-aware Cohen swap formulas — so no full
> recomputation is needed.

## Algorithm

latticereduce.c

LatticeReduce[{v1, v2, ...}] -- an LLL-reduced basis for the lattice spanned by the row vectors v_i.

```text
  LatticeReduce[m]   m an n x d matrix (List of n equal-length Lists).
                     Returns an n x d matrix whose rows form a reduced
                     basis of the same lattice (same Z- / Z[i]-module).
```

Entries may be:

```text
  - integers (machine int64 or GMP bigint),
  - rationals (Rational[p, q]),
  - Gaussian integers / Gaussian rationals (Complex[a, b] with a, b
    integer or rational).
```

Algorithm: the classical Lenstra-Lenstra-Lovas (LLL) reduction with

```text
Lovas parameter delta = 3/4, run in EXACT arithmetic.  Every quantity
```

is an exact Gaussian rational stored as a pair of GMP `mpq_t` (`GRat`), so the routine is correct for both machine-size and arbitrary-precision (bignum) inputs -- floating point is never used, which matters because LatticeReduce is most often used to discover integer relations where a rounding error would yield a wrong relation.

The Gram-Schmidt orthogonalisation is generalised to the complex

```text
(Hermitian) inner product  <x, y> = sum_k x_k conj(y_k),  so the same
code path handles real and Gaussian lattices.  Size reduction rounds
```

the mu coefficients to the nearest Gaussian integer; the Gram-Schmidt data (mu, |b*|^2) is maintained incrementally -- computed once up front, updated in place on size reduction, and updated on a Lovas swap via the conjugate-aware Cohen swap formulas (no full recomputation).

Because every basis transformation is an integer (Z, or Z[i] in the Gaussian case) row operation or a row swap, the lattice -- and hence Abs[Det] and every linear relation in the right null space -- is preserved exactly.

Linearly independent rows are required (every documented use of LatticeReduce -- integer-relation finding, basis reduction -- supplies

```text
a full-rank generating set).  A dependent generating set is detected
```

during Gram-Schmidt and the call is left unevaluated with a diagnostic.

```text
Memory ownership: standard builtin contract.  This file does NOT free
`res` on success or failure -- the evaluator owns it (see MEMORY.md /
```

SPEC.md S4.1).

## Implementation notes

**Algorithm.** `builtin_latticereduce` returns an **LLL-reduced** basis for the lattice spanned by the row vectors of the input matrix, using the classical Lenstra–Lenstra–Lovász reduction with Lovász parameter `δ = 3/4`, run entirely in **exact arithmetic**. Gram–Schmidt orthogonalisation is generalised to the Hermitian inner product `⟨x,y⟩ = Σ x_k conj(y_k)`, so the same code reduces real lattices and Gaussian (complex) lattices. The Gram–Schmidt data — the `μ` coefficients and the squared norms `|b*|²` — is maintained incrementally: computed once, updated in place on each size-reduction step (rounding `μ` to the nearest Gaussian integer), and updated on each Lovász swap via Cohen's conjugate-aware swap formulas (no full recomputation). Because every basis transformation is an integer (`Z`, or `Z[i]`) row operation or row swap, the lattice — and hence `Abs[Det]` and every relation in the right null space — is preserved exactly.

**Data structures.** Every scalar is an exact Gaussian rational `GRat` = a pair of GMP `mpq_t` (`re`, `im`); floating point is never used, which is essential when the reduction is used to discover integer relations where a rounding error would give a wrong relation. Inputs may be machine/bignum integers, rationals, or Gaussian integers/rationals. The basis is a dense array of `GRat` vectors.

**Complexity / limits.** Linearly independent rows are required; a rank-deficient generating set is detected during Gram–Schmidt and the call is left unevaluated with a diagnostic. LLL is polynomial in the dimension and the bit-size of the entries; the exact `mpq_t` arithmetic trades speed for guaranteed correctness.

- `Protected`.
- Returns an `n × d` matrix whose rows form a reduced basis of the same
  lattice (same `Z`- / `Z[i]`-module).
- Entries may be integers, rationals, Gaussian integers, or Gaussian
  rationals (`Complex[a, b]` with exact `a`, `b`).
- Exact Lenstra–Lenstra–Lovász reduction (Lovász parameter `δ = 3/4`)
  carried out entirely in GMP rational arithmetic — no floating point —
  so it is correct for both machine-size and arbitrary-precision (bignum)
  inputs. This matters for integer-relation finding, where a rounding
  error would yield a wrong relation.
- The Gram–Schmidt orthogonalisation uses the Hermitian inner product
  `⟨x, y⟩ = Σ x_k conj(y_k)`, so real and Gaussian lattices share one
  code path; size reduction rounds to the nearest Gaussian integer.
- The lattice — and hence `Abs[Det]` and every linear relation in the
  right null space — is preserved exactly.
- The rows must be linearly independent; a dependent generating set is
  reported with `LatticeReduce::dep`.
- Diagnostics: `LatticeReduce::argx` (wrong argument count),
  `LatticeReduce::matrix` (not a non-empty rectangular matrix),
  `LatticeReduce::latm` (an entry is not rational).

**Attributes:** `Protected`.

## References

- A. K. Lenstra, H. W. Lenstra, L. Lovász, "Factoring Polynomials with Rational Coefficients", Mathematische Annalen 261 (1982).
- Henri Cohen, *A Course in Computational Algebraic Number Theory* (Springer, 1993).
- Source: [`src/linalg/latticereduce.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/latticereduce.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_latticereduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_latticereduce.c)

## Notes & additional examples

### Notes

`LatticeReduce[m]` returns an LLL-reduced basis for the lattice spanned by the
rows of `m`. Entries may be integers, Gaussian integers, rationals, or Gaussian
rationals; arithmetic is exact (GMP rationals), so results are correct for both
machine-size and arbitrary-precision input. The reduction preserves the lattice,
its determinant, and every linear relation among the rows. The input rows must be
linearly independent.
