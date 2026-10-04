---
source: src/poly/algebraicintegerq.c
references:
  - "H. Cohen, *A Course in Computational Algebraic Number Theory*, GTM 138 (Springer, 1993), §4.1 (algebraic integers and minimal polynomials)."
  - "The FLINT library (https://flintlib.org), `qqbar` module — exact real and complex algebraic numbers via minimal polynomial plus isolating enclosure."
---
**Algorithm.** `builtin_algebraicintegerq` is a thin wrapper: it checks arity 1
and hands the argument to `flint_qqbar_algebraic_integer_q`, mapping the
tri-state result to `True` (1), `False` (0), or unevaluated (−1, FLINT compiled
out). The decision is exact. The engine converts `x` to a FLINT `qqbar_t` via
`to_qqbar` — the shared converter that accepts integers, rationals, radicals
`Power[base, p/q]`, roots of unity, the imaginary unit, `Root[]` objects and
`GoldenRatio`, combined by `+ - * / ^`. `x` is an algebraic integer **iff the
leading coefficient of its primitive integer minimal polynomial is 1** (monic):
`fmpz_poly_get_coeff_fmpz(lead, QQBAR_POLY(v), deg)` and test `fmpz_is_one`. A
rational `p/q` has minimal polynomial `q x − p`, so only true integers (`q = 1`)
qualify. Anything that is not a constant algebraic number — a free symbol, `Pi`,
`Log[2]` — fails the `to_qqbar` conversion and returns `0` (`False`), matching
WL.

**Data structures.** The FLINT `qqbar_t` (an exact algebraic number: its
primitive integer minimal polynomial as an `fmpz_poly` plus an isolating complex
enclosure). No `AlgebraicNumber[...]` object is built — this is a predicate that
only reads the minimal polynomial's leading coefficient.

**Complexity / limits.** Cost is dominated by the `to_qqbar` conversion of the
argument (field arithmetic over the compositum, bounded by the degree cap
`QQBAR_DEGREE_CAP = 120`). `Protected`, and deliberately **not** `Listable` —
`AlgebraicIntegerQ[list]` asks whether the list itself is an algebraic integer
(`False`), not whether its elements are. Declines (stays unevaluated) only when
FLINT is unavailable.
