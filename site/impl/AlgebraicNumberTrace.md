---
source: src/poly/algebraicnumbertrace.c
references:
  - "H. Cohen, *A Course in Computational Algebraic Number Theory*, GTM 138 (Springer, 1993), §4.3 (norm and trace of algebraic numbers; transitivity in a tower)."
  - "The FLINT library (https://flintlib.org), `qqbar` module — minimal polynomial and field membership."
---
**Algorithm.** `builtin_algebraicnumbertrace` mirrors `AlgebraicNumberNorm`
structurally: it separates a trailing `Extension -> theta` option, checks arity,
and delegates to `flint_qqbar_algebraic_number_trace`. The **absolute trace**
`Tr_{Q(a)/Q}(a)` is the sum of the roots of `a`'s primitive integer minimal
polynomial `P(x) = c_n x^n + … + c_0`, equal to `−c_{n−1}/c_n`
(`qqbar_abs_trace`). With `Extension -> theta` the **relative trace**
`Tr_{Q(theta)/Q}(a)` follows from transitivity of the trace in the tower
`Q ⊆ Q(a) ⊆ K = Q(theta)`:

```
Tr_{K/Q}(a) = [K:Q(a)] · Tr_{Q(a)/Q}(a) = (n/d) · (absolute trace)
```

with `n = deg minpoly(theta)`, `d = deg minpoly(a)`, and `a ∈ Q(theta)` decided
by `qqbar_express_in_field`. The contrast with the norm is exact: the trace is
**additive**, so it *scales* by the tower index `n/d`, whereas the multiplicative
norm is *raised to the power* `n/d`. Return codes map to the
`AlgebraicNumberTrace::ext` (`a` not in `Q(theta)`) and
`AlgebraicNumberTrace::nalg` (not a constant algebraic number) messages through
`mth_message`, a success `Integer`/`Rational`, or a silent decline (FLINT off).

**Data structures.** FLINT `qqbar_t` for `a` and `theta`, an `fmpq` accumulator
for the trace, and an `fmpq_poly` for membership in the relative case.

**Complexity / limits.** `O(1)` minimal-polynomial read for the absolute trace;
the relative case adds a primitive-element build and an express-in-field solve.
`Listable`, `Protected`; bounded by the degree cap `QQBAR_DEGREE_CAP = 120`.
