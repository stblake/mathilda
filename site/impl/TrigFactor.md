---
source: src/simp/trigsimp.c
---
**Algorithm.** `builtin_trigfactor_impl` factors trig expressions — broadly the
inverse of TrigExpand for the structural identities both support. It tries two
parallel paths and keeps the shorter result (by `trigfactor_leaf_count`):

- **Path A:** `trigfactor_run_pipeline` on the argument as-is. If this changes
  the expression at all it is trusted and Path B is skipped (Path B can be
  expensive on angle-sum arguments).
- **Path B:** only attempted when Path A was a no-op *and* the input has compound
  trig structure (`has_compound_trig_structure` — a `Power[trig, k≥2]` or a
  `Times` of two or more trig atoms). First `TrigExpand` the argument, then run
  the pipeline; kept only if strictly shorter than Path A. This catches
  cancellations that surface only after angle-sum expansion (e.g. `Cos[x+y] +
  Sin[x] Sin[y] -> Cos[x] Cos[y]`).

`trigfactor_run_pipeline` (trig canonicalizer suppressed throughout):
1. `ReplaceRepeated` with `trig_factor_to_sincos` — rewrite reciprocal heads as
   `Sin/Cos` ratios so `Factor` sees full polynomial structure.
2. `Together` over a common denominator.
3. `Factor` — Mathilda `Factor` treats trig atoms as polynomial variables.
   Skipped (left as the `Together` output) when the polynomial would stall:
   more than `TRIG_FACTOR_ATOM_THRESHOLD` distinct squared atoms, max trig-atom
   power above `TRIG_FACTOR_DEGREE_THRESHOLD`, or more than
   `TRIG_FACTOR_TOTAL_ATOM_THRESHOLD` distinct atoms.
4. `ReplaceRepeated` with `trig_factor_identities` — Pythagorean collapses (both
   signs, circular and hyperbolic), reverse angle-addition, reverse double-angle,
   factored-form `(Cos-Sin)(Cos+Sin) -> Cos[2x]` and `(Cosh±Sinh)` collapses, and
   the Weierstrass-style linear-combination factoring `a Sin[x] + b Cos[x] ->
   Sqrt[a^2+b^2] Sin[x + ArcTan[a, b]]` (gated by `NumberQ` and `Im == 0` on both
   coefficients so complex coefficients cannot collapse `Sqrt[a^2+b^2]` to zero).
5. `ReplaceRepeated` with `trig_factor_from_sincos` to restore `Tan`/`Sec`/...

**Data structures.** Static rule lists from `trigsimp_init`; tail patterns
(`r___`) on `Plus`/`Times` let each identity fire inside a larger sum/product,
with `ATTR_ORDERLESS` on `Plus`/`Times` driving the permutation search.
Memoized through the active `FactorMemo` via the `builtin_trigfactor` wrapper.
