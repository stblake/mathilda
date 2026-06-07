---
references:
  - "A. K. Lenstra, H. W. Lenstra, L. Lovász, \"Factoring Polynomials with Rational Coefficients\", Mathematische Annalen 261 (1982)."
  - "Henri Cohen, *A Course in Computational Algebraic Number Theory* (Springer, 1993)."
source: src/linalg/latticereduce.c
---
**Algorithm.** `builtin_findintegernullvector` finds integers (or Gaussian integers, for complex inputs) `a = {a_1,…,a_n}`, not all zero, with `Σ a_i x_i = 0`, by **integer-relation detection via LLL** rather than PSLQ. It builds the relation lattice whose `i`-th row is `r_i = (e_i | round(2^b · x_i))` — the standard basis vector augmented with the scaled, rounded coordinate (Gaussian rounding when `x_i` is complex) — LLL-reduces it exactly using the same machinery as `LatticeReduce`, and reads the candidate relation off the leading components of the shortest reduced row. A rigorous certificate is computed from the LLL Gram–Schmidt bound `λ_1(L)² ≥ M2` (with `M2 = min_i |b*_i|²`) combined with the worst-case rounding error `|a·round(2^b x)| ≤ (√n/2)‖a‖`, giving a lower bound `B = √(M2 / (1 + (√n/2)²))` on the norm of any relation. `B` drives the no-relation / not-found diagnostics (`norel`, `lgrelb`, `rnfb`, `rnfu`).

**Data structures.** Exact Gaussian-rational `GRat` (pair of GMP `mpq_t`) scalars throughout, matching `LatticeReduce`. For inexact inputs the working precision `b` is derived from the inputs' Real/MPFR precision (`finv_max_prec_bits` under `USE_MPFR`).

**Complexity / limits.** Polynomial in `n` and the scaled bit-size; the exactness of the LLL pass plus the analytic norm bound let it both return a verified relation and *prove* none exists below a given norm (reported via the diagnostics).
