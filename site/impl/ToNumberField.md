---
source: src/poly/tonumberfield.c
references:
  - "H. Cohen, *A Course in Computational Algebraic Number Theory*, GTM 138 (Springer, 1993), §4.2 (power-basis representation) and the primitive-element theorem for a compositum."
  - "The FLINT library (https://flintlib.org), `qqbar` module — field membership (`qqbar_express_in_field`) and primitive-element construction."
---
**Algorithm.** `builtin_tonumberfield` classifies the argument forms and defers
every field computation to `src/poly/flint_qqbar.c`; each produced
`AlgebraicNumber` is re-canonicalised by the evaluator.

- `ToNumberField[a, theta]` → `flint_qqbar_to_number_field`. Convert `a` and
  `theta` with `to_qqbar`, build the algebraic-integer generator
  `phi = algint_generator(theta)` of degree `n`, then express `a` as a
  polynomial in `phi` with `qqbar_express_in_field_esc_all` — which *escalates*
  the working precision for a named generator (fix A14), so a genuine relation
  that needs more than the default precision is not misread as non-membership.
  `poly_to_algnum` renders `AlgebraicNumber[g, {..}]` (or a rational); the result
  is `NULL` — leaving `ToNumberField` unevaluated — exactly when `a ∉ Q(theta)`.
- `ToNumberField[x]` → `flint_qqbar_to_number_field_self`, i.e. `x` expressed in
  `Q(x)` itself.
- `ToNumberField[{a1, …, ak}]`, `[…, Automatic]`, `[…, All]` →
  `flint_qqbar_to_number_field_common`. A **primitive element** of the
  compositum is built as `alpha = a1 + Σ c_i a_i`, with each `c` chosen by
  **degree** rather than by trial membership: `alpha + c·b ∈ Q(alpha, b)` for
  every `c`, so the largest `deg(alpha + c·b)` over `c ∈ 1..16` is the compositum
  degree (all but finitely many `c` give a primitive element). Choosing by degree
  sidesteps a membership test that could fail for want of precision — the defect
  that once made `ToNumberField[{2^(1/6), 2^(1/3)}]` (a degree-6 field resolvable
  only above 64 bits) come back unevaluated. Each `a_i` is then expressed in
  `Q(alpha)`.
- `ToNumberField[{a1, …, ak}, theta]` → each `a_i` expressed in `Q(theta)`;
  declines if any `a_i ∉ Q(theta)`.

**Data structures.** FLINT `qqbar_t` for the inputs, the primitive element and
`phi`; `fmpq_poly` for the field expression; the `AlgebraicNumber[theta,
{coeffs}]` representation for every output element. The common-field path uses a
`calloc`'d result array so an early express-in-field miss leaves the untouched
tail `NULL` for the error-path cleanup.

**Complexity / limits.** The single-generator form costs one express-in-field
solve; the common-field form adds a primitive-element search (16 candidate `c`
per adjoined generator). Bounded by the degree cap `QQBAR_DEGREE_CAP = 120`.
`Protected`, and deliberately **not** `Listable`: the `{a1, …}` form defines one
common field over the whole list, so threading element-wise would destroy its
meaning.
