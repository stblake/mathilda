---
source: src/poly/numberfieldintegralbasis.c
references:
  - "H. Cohen, *A Course in Computational Algebraic Number Theory*, GTM 138 (Springer, 1993), §6.1 (the Round 2 / Pohst–Zassenhaus maximal-order algorithm) and §4.4 (integral bases)."
  - "The FLINT library (https://flintlib.org), `qqbar` module, `fmpz_mat` Hermite normal form (`fmpz_mat_hnf`), and the number-field maximal-order layer."
---
**Algorithm.** `builtin_numberfieldintegralbasis` checks arity 1 and delegates
to `flint_qqbar_integral_basis`, which returns a `Z`-module basis of the ring of
integers `O_K` of `K = Q(a)`:

1. Build the algebraic-integer generator `phi = lc·alpha` of `Q(a)`
   (`algint_generator`), monic of degree `n`. If `n ≤ 1`, `O_K = Z` and the
   basis is `{1}`.
2. Hand `phi`'s monic integer defining polynomial to the number-field layer
   (`nf_field_create`), which runs **Dedekind's criterion** at every ramified
   prime and, where the equation order `Z[phi]` is not maximal, enlarges it to
   `O_K` by **Round 2 (Pohst–Zassenhaus)**. This is what makes the basis correct
   for non-monogenic fields, not just `{1, phi, …, phi^{n−1}}`.
3. Read back the numerator lattice `W` and common denominator `D`
   (`O_K = (1/D)·L`). The rows of `W` are arbitrarily ordered; the engine puts
   them in the standard presentation — the lower-triangular **Hermite normal
   form** in the `theta`-power basis, so `omega_0 = 1` and `omega_k` has degree
   exactly `k`. FLINT's `fmpz_mat_hnf` gives the upper-triangular HNF, so the
   lower-triangular one is obtained as a 180°-rotation applied to the
   column-reversed lattice.
4. Render each basis row `omega_i = (1/D) Σ_j H[i][j] phi^j` with
   `poly_to_algnum` — an `AlgebraicNumber[g, {..}]`, or a plain
   `Integer`/`Rational` when the row is rational (as `omega_0 = 1` always is).

A decline routes a `NumberFieldIntegralBasis::nintbas` message through
`mth_message`.

**Data structures.** FLINT `qqbar_t` for `phi`; the `NumberField*` maximal-order
object (its Round-2 basis `W` and denominator `D` as `mpz`/`fmpz`); `fmpz_mat`
for the HNF; `fmpq_poly` per row; and the `AlgebraicNumber` representation for
each returned element.

**Complexity / limits.** Dominated by the maximal-order computation —
factorisation of the field discriminant and the `p`-radical / Round-2 steps at
each ramified prime. `Listable` (each generator is its own field), `Protected`;
bounded by the degree cap `QQBAR_DEGREE_CAP = 120`. Declines when `O_K` cannot be
certified (the discriminant will not factor, or a prime exceeds the single-word
modulus) or FLINT is compiled out. The basis returned is correct but **not
unique** (any `Z`-basis of `O_K` is admissible; this one is the HNF
presentation).
