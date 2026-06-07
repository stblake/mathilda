---
source: src/poly/groebnerbasis.c
references:
  - "B. Buchberger, \"An Algorithm for Finding the Basis Elements of the Residue Class Ring of a Zero Dimensional Polynomial Ideal\" (PhD thesis, 1965)."
  - "R. Gebauer, H. M. Möller, \"On an installation of Buchberger's algorithm\", J. Symbolic Computation 1988."
  - "S. Collart, M. Kalkbrener, D. Mall, \"Converting bases with the Gröbner walk\", J. Symbolic Computation 1997."
  - "T. Becker, V. Weispfenning, *Gröbner Bases* (Springer, 1993)."
---
**Algorithm.** `builtin_groebner_basis` (in `src/poly/groebnerbasis.c`) is the front end; the engine is in `groebner.c` and the order conversion in `groebnerwalk.c`. The front end parses `GroebnerBasis[polys, vars]` (and the 3-arg elimination form `GroebnerBasis[polys, mainVars, elimVars]`), reads the `MonomialOrder` (Lexicographic default, DegreeReverseLexicographic, or EliminationOrder — forced by the 3-arg form), `CoefficientDomain` and `Method` options, converts each polynomial to the internal `GBPoly` via `gb_from_expr`, runs the core, and renders the basis back to expressions.

The core `gb_buchberger` is textbook **Buchberger** with **Gebauer–Möller** pair management: each new basis element is folded in with `gm_update`, which applies Buchberger criterion 1 (coprime leading monomials) and the Gebauer–Möller M/F/B criteria to discard provably useless pairs without re-walking exponents. The main loop picks the next pair by the normal strategy (`gm_pick_pair`), forms the S-polynomial (`gb_spoly`, using the pointwise-max lcm of the two leading exponent vectors), fully reduces it against the current basis (`gb_reduce`, multivariate division to a normal form), and adds any non-zero remainder (made monic). After saturation, `gb_finalize_basis` inter-reduces and normalises to the unique reduced Gröbner basis. A `tc_check_deadline` hook lets `TimeConstrained` abort between pairs.

For expensive target orders the front end can route through the **Gröbner walk** (`gb_groebner_walk`, the Collart–Kalkbrener–Mall algorithm): compute the basis under the cheap DegreeReverseLexicographic order, then walk along the weight path `w(t) = (1−t)·(1,…,1) + t·σ` across the Gröbner fan, rewriting through initial-form ideals at each cone wall, with a guaranteed fall-back to a direct `gb_buchberger` run under the target order on any overflow or degenerate lift.

**Data structures.** `GBPoly` holds rational coefficients as a parallel array of GMP `mpq_t` alongside a row-major `int* exps` (`n_terms × n_vars`) exponent matrix, with a borrowed `GBWeightMatrix` defining the order and an `elim_pivot` for elimination blocks. Pairs are `GMPair` records carrying a precomputed lcm and a `dead` flag.

**Complexity / limits.** Buchberger's algorithm is doubly-exponential in the worst case; lexicographic bases especially can blow up, which is why grevlex-then-walk is offered. Arithmetic is exact (`mpq_t`/`mpz`-guarded int64 weights). Unsupported option values emit a `GroebnerBasis::nimpl` diagnostic.
