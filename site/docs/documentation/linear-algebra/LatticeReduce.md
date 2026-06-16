# LatticeReduce

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
LatticeReduce[m]
    gives an LLL-reduced basis for the lattice spanned by the rows
    (vectors) of m.  The entries of m may be integers, Gaussian
    integers, rationals, or Gaussian rationals.  Reduction is exact
    (GMP rational arithmetic, so it is correct for both machine-size
    and arbitrary-precision entries) and preserves the lattice, its
    determinant, and every linear relation among the rows.  The rows
    must be linearly independent.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= LatticeReduce[{{1, 0, 0, 1345}, {0, 1, 0, 35}, {0, 0, 1, 154}}]
Out[1]= {{0, 9, -2, 7}, {1, 1, -9, -6}, {1, -3, -8, 8}}

In[2]:= {w1, w2} = LatticeReduce[{{12, 2}, {13, 4}}]
Out[2]= {{1, 2}, {9, -4}}

In[3]:= b . {1, 2, 3, 1}    (* relations preserved *)
Out[3]= Dot[b, {1, 2, 3, 1}]

In[4]:= LatticeReduce[{{1, 2}, {3, 4.5}}]
Out[4]= LatticeReduce[{{1, 2}, {3, 4.5}}]
```

## Implementation notes

**Algorithm.** `builtin_latticereduce` returns an **LLL-reduced** basis for the lattice spanned by the row vectors of the input matrix, using the classical Lenstra–Lenstra–Lovász reduction with Lovász parameter `δ = 3/4`, run entirely in **exact arithmetic**. Gram–Schmidt orthogonalisation is generalised to the Hermitian inner product `⟨x,y⟩ = Σ x_k conj(y_k)`, so the same code reduces real lattices and Gaussian (complex) lattices. The Gram–Schmidt data — the `μ` coefficients and the squared norms `|b*|²` — is maintained incrementally: computed once, updated in place on each size-reduction step (rounding `μ` to the nearest Gaussian integer), and updated on each Lovász swap via Cohen's conjugate-aware swap formulas (no full recomputation). Because every basis transformation is an integer (`Z`, or `Z[i]`) row operation or row swap, the lattice — and hence `Abs[Det]` and every relation in the right null space — is preserved exactly.

**Data structures.** Every scalar is an exact Gaussian rational `GRat` = a pair of GMP `mpq_t` (`re`, `im`); floating point is never used, which is essential when the reduction is used to discover integer relations where a rounding error would give a wrong relation. Inputs may be machine/bignum integers, rationals, or Gaussian integers/rationals. The basis is a dense array of `GRat` vectors.

**Complexity / limits.** Linearly independent rows are required; a rank-deficient generating set is detected during Gram–Schmidt and the call is left unevaluated with a diagnostic. LLL is polynomial in the dimension and the bit-size of the entries; the exact `mpq_t` arithmetic trades speed for guaranteed correctness.

- `Protected`.
- Returns an `n × d` matrix whose rows form a reduced basis of the same

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- A. K. Lenstra, H. W. Lenstra, L. Lovász, "Factoring Polynomials with Rational Coefficients", Mathematische Annalen 261 (1982).
- Henri Cohen, *A Course in Computational Algebraic Number Theory* (Springer, 1993).
- Source: [`src/linalg/latticereduce.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/latticereduce.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= LatticeReduce[{{1, 1, 1}, {-1, 0, 2}, {3, 5, 6}}]
Out[1]= {{0, 1, 0}, {1, 0, 1}, {-2, 0, 1}}
```

LLL produces a short, nearly-orthogonal basis for the same lattice. Reduction preserves the lattice determinant exactly:

```mathematica
In[1]:= Det[LatticeReduce[{{201, 37}, {1648, 297}}]]
Out[1]= -1279

In[2]:= Det[{{201, 37}, {1648, 297}}]
Out[2]= -1279
```

Integer-relation detection: appending a scaled approximation column to the identity makes a short vector reveal a relation `61*pi - 183*e + 189*phi ~ 0` (the last coordinate is tiny):

```mathematica
In[1]:= LatticeReduce[{{1, 0, 0, 31415927}, {0, 1, 0, 27182818}, {0, 0, 1, 16180340}}]
Out[1]= {{61, -183, 189, 113}, {-198, 108, 203, -182}, {-235, 146, 211, 323}}
```

Reduction is exact, so rational and Gaussian-rational lattices are handled without rounding:

```mathematica
In[1]:= LatticeReduce[{{1/2, 1}, {1, 1/3}}]
Out[1]= {{1/2, -2/3}, {1, 1/3}}
```

### Notes

`LatticeReduce[m]` returns an LLL-reduced basis for the lattice spanned by the
rows of `m`. Entries may be integers, Gaussian integers, rationals, or Gaussian
rationals; arithmetic is exact (GMP rationals), so results are correct for both
machine-size and arbitrary-precision input. The reduction preserves the lattice,
its determinant, and every linear relation among the rows. The input rows must be
linearly independent.
