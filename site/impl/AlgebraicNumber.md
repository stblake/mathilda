---
source: src/poly/algebraicnumber.c
references:
  - "H. Cohen, *A Course in Computational Algebraic Number Theory*, GTM 138 (Springer, 1993), §4.2–4.3 (representation of number-field elements in a power basis)."
  - "The FLINT library (https://flintlib.org), `qqbar` module (exact algebraic numbers) and `nf`/`nf_elem` modules (number-field element arithmetic)."
---
**Algorithm.** `builtin_algebraicnumber` validates `AlgebraicNumber[theta,
{c0..cn}]` (arity 2) and delegates all field computation to
`flint_qqbar_algebraic_number`; the builtin itself only canonicalises and guards
against re-evaluation churn. The engine:

1. Builds the **algebraic-integer generator** `phi = lc·alpha` of `Q(theta)`,
   where `alpha = to_qqbar(theta)` and `lc` is the positive leading coefficient
   of `alpha`'s primitive integer minimal polynomial (`algint_generator`); `phi`
   is monic of degree `n`.
2. Forms `p(x) = Σ (c_i / lc^i) x^i` over `Q` and reduces it modulo `phi`'s monic
   minimal polynomial `M` (`fmpq_poly_rem`), giving a representative of degree
   `< n` in the power basis of `phi`.
3. Renders it with `poly_to_algnum`: a plain `Integer`/`Rational` when the value
   is rational (only the constant coefficient survives), otherwise the reduced
   `AlgebraicNumber[g, {d0..d_{n-1}}]` with the coefficient list padded to length
   `n` (the degree of the minimal polynomial).

A **fixpoint guard** (`expr_eq(cand, res)` → return `NULL`) leaves an input that
is already canonical untouched, so the evaluator does not loop — the
canonicalisation is idempotent. A **fast path** (via the `GENFD` generator-field
cache) short-circuits when `theta` is already the canonical integer generator
and the coefficients are already rational and reduced, returning the padded form
directly without re-running `to_qqbar`.

**Data structures.** The `AlgebraicNumber[theta, {c0..cn}]` representation
itself; FLINT `qqbar_t` for `alpha`/`phi`; `fmpq_poly` for `p` and `M`; and the
per-generator `GENFD` cache keyed by the generator expression (the minimal
polynomial is a pure function of it, so the cache never goes stale). The head
carries `NHoldAll` so `N` reaches the dedicated AlgebraicNumber branch rather
than threading into the generator and coefficient list.

**Complexity / limits.** Dominated by `to_qqbar` of the generator and the
`fmpq_poly` reduction modulo `M`; bounded by the degree cap
`QQBAR_DEGREE_CAP = 120`. Declines (stays unevaluated) for a non-algebraic
generator, malformed (non-integer/rational) coefficients, a degree-cap overflow,
or FLINT compiled out.
