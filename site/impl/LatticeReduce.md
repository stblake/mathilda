---
references:
  - "A. K. Lenstra, H. W. Lenstra, L. Lovász, \"Factoring Polynomials with Rational Coefficients\", Mathematische Annalen 261 (1982)."
  - "Henri Cohen, *A Course in Computational Algebraic Number Theory* (Springer, 1993)."
source: src/linalg/latticereduce.c
---
**Algorithm.** `builtin_latticereduce` returns an **LLL-reduced** basis for the lattice spanned by the row vectors of the input matrix, using the classical Lenstra–Lenstra–Lovász reduction with Lovász parameter `δ = 3/4`, run entirely in **exact arithmetic**. Gram–Schmidt orthogonalisation is generalised to the Hermitian inner product `⟨x,y⟩ = Σ x_k conj(y_k)`, so the same code reduces real lattices and Gaussian (complex) lattices. The Gram–Schmidt data — the `μ` coefficients and the squared norms `|b*|²` — is maintained incrementally: computed once, updated in place on each size-reduction step (rounding `μ` to the nearest Gaussian integer), and updated on each Lovász swap via Cohen's conjugate-aware swap formulas (no full recomputation). Because every basis transformation is an integer (`Z`, or `Z[i]`) row operation or row swap, the lattice — and hence `Abs[Det]` and every relation in the right null space — is preserved exactly.

**Data structures.** Every scalar is an exact Gaussian rational `GRat` = a pair of GMP `mpq_t` (`re`, `im`); floating point is never used, which is essential when the reduction is used to discover integer relations where a rounding error would give a wrong relation. Inputs may be machine/bignum integers, rationals, or Gaussian integers/rationals. The basis is a dense array of `GRat` vectors.

**Complexity / limits.** Linearly independent rows are required; a rank-deficient generating set is detected during Gram–Schmidt and the call is left unevaluated with a diagnostic. LLL is polynomial in the dimension and the bit-size of the entries; the exact `mpq_t` arithmetic trades speed for guaranteed correctness.
