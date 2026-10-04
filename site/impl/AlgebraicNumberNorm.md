---
source: src/poly/algebraicnumbernorm.c
references:
  - "H. Cohen, *A Course in Computational Algebraic Number Theory*, GTM 138 (Springer, 1993), §4.3 (norm and trace of algebraic numbers; transitivity in a tower)."
  - "The FLINT library (https://flintlib.org), `qqbar` module — minimal polynomial and field membership."
---
**Algorithm.** `builtin_algebraicnumbernorm` separates a trailing
`Extension -> theta` option from the positional argument
(`extract_extension_option`), checks arity, and delegates to
`flint_qqbar_algebraic_number_norm`. The **absolute norm** `N_{Q(a)/Q}(a)` is
read straight off `a`'s primitive integer minimal polynomial
`P(x) = c_n x^n + … + c_0` (content 1, `c_n > 0`): the product of its `n` roots
is `(−1)^n · c_0/c_n` (`qqbar_abs_norm`). With `Extension -> theta` the
**relative norm** `N_{Q(theta)/Q}(a)` is computed by transitivity of the norm in
the tower `Q ⊆ Q(a) ⊆ K = Q(theta)`:

```
N_{K/Q}(a) = N_{Q(a)/Q}(a)^{[K:Q(a)]} = (absolute norm)^{n/d}
```

with `n = deg minpoly(theta)` (via the algebraic-integer generator `phi`) and
`d = deg minpoly(a)`; membership `a ∈ Q(theta)` is decided by
`qqbar_express_in_field` (escalating working precision). The engine's return
codes map to: success (`Integer`/`Rational`); `AlgebraicNumberNorm::ext` (`a`
not an element of `Q(theta)`); `AlgebraicNumberNorm::nalg` (not a constant
algebraic number); silent decline (FLINT off). Both messages route through
`mth_message`.

**Data structures.** FLINT `qqbar_t` for `a` and `theta`, an `fmpq` accumulator
for the norm, and an `fmpq_poly` for the field-membership expression in the
relative case. No `AlgebraicNumber` object is constructed.

**Complexity / limits.** The absolute norm is `O(1)` once `a` is converted (one
minimal-polynomial read); the relative case adds a primitive-element build and an
express-in-field solve. `Listable` (a trailing `Extension` option is repeated
across the list), `Protected`; bounded by the degree cap
`QQBAR_DEGREE_CAP = 120`.
