---
source: src/rootreduce.c
references:
  - "H. Cohen, *A Course in Computational Algebraic Number Theory*, GTM 138 (Springer, 1993), ch. 4 (algebraic numbers, minimal polynomials, number fields)."
  - "The FLINT library (https://flintlib.org), `qqbar` module — exact real and complex algebraic numbers via minimal polynomial plus isolating enclosure, with no numeric zero oracle."
---
**Algorithm.** `builtin_rootreduce` is a dispatcher over three rigorous FLINT
engines, chosen by the shape of the argument. It first splits a trailing
`Method -> "Automatic" | "Recursive" | "NumberField"` option (recognised by
reading `Options[RootReduce]`, so a `Solve`-result rule `u -> value` is kept as
a positional argument, not mistaken for an option).

1. **Constant algebraic number** (no free symbol) → `flint_qqbar_canonical`
   (`src/poly/flint_qqbar.c`). `to_qqbar` converts an expression built from
   integers, rationals, radicals `Power[base, p/q]`, roots of unity `(−1)^(p/q)`,
   the imaginary unit, `Root[]` objects and `GoldenRatio`, combined by
   `+ − * / ^`, into a single `qqbar_t`. `qqbar_to_expr` then renders a canonical
   representative: a **rational** when `qqbar_is_rational`; `(a + b Sqrt[c])/q`
   for **degree 2** via `qqbar_get_quadratic`; and `Root[Function[minpoly &], k]`
   for **degree ≥ 3**, with the index `k` placed by WL's root ordering
   (`wl_root_index`) and the whole object memoised by `(minpoly, k)`.
   `Method -> "NumberField"` re-expresses the value through a single primitive
   element of the field (`number_field_value`).
2. **Parametric algebraic function** — a radical whose radicand carries a free
   variable — → `flint_algebraic_field_canonical` (`src/poly/flint_bridge.c`):
   the denominator is inverted in the field by an exact linear solve over the
   tower, with no numeric oracle.
3. **Polynomial / rational function in a free variable** with
   constant-algebraic coefficients → `flint_qqbar_reduce_coeffs`: each maximal
   constant-algebraic subexpression (a coefficient) is canonicalised via `qqbar`
   while the free-variable structure is left intact, so a vanishing radical
   coefficient reduces to `0` and its monomial drops out. (Plain polynomial
   cancellation is `Cancel`'s job, not done here.)

`RootReduce` also **threads** over `Equal`/`Unequal`/`Less`/`…`/`And`/`Or` and
over an immediate `Rule` (a `Solve` result entry). For a binary (in)equality of
constant algebraic numbers it is *decided exactly* — `flint_qqbar_equal`
(1/0/−1) for `Equal`/`Unequal`, `flint_qqbar_compare` (sign, or −2 undecided)
for the ordering relations — so `Sqrt[2] + Sqrt[3] == Sqrt[5 + 2 Sqrt[6]]`
returns `True` with no floating point. When the argument carries no algebraic
content it is returned unchanged (the positional arg is stolen out of `res`).

**Data structures.** FLINT `qqbar_t` (an exact algebraic number: primitive
integer minimal polynomial as an `fmpz_poly` plus an isolating complex
enclosure). Three session caches accelerate the hot paths: a `to_qqbar`
conversion cache, the shared `wl-roots` ordering cache, and a `(minpoly, k)` →
`Root[]`-object memo (`q2e_cache`). Equality and comparison are rigorous
decisions on the minimal polynomials, **not** a numeric zero test.

**Complexity / limits.** Dominated by `qqbar` arithmetic over the compositum of
the atoms; bounded by the degree cap `QQBAR_DEGREE_CAP = 120`, above which it
declines. A non-real ordering comparison is undecided (`−2`) and leaves the
relation unevaluated rather than guessing. `Listable`, `Protected`; option
`Method -> "Automatic" | "Recursive" | "NumberField"`.
