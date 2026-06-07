---
source: src/simp/trigsimp.c
---
**Algorithm.** `builtin_trigexpand_impl` expands trig/hyperbolic functions of
sums and integer multiples into products and powers of single-argument trig
calls. The pipeline (with the trig canonicalizer suppressed throughout):

1. **Angle-addition + multiple-angle.** `ReplaceRepeated` with
   `trig_expand_rules`. The binary angle-addition forms
   `Sin[x_ + y__] :> Sin[x] Cos[Plus[y]] + Cos[x] Sin[Plus[y]]` (and `Cos`,
   `Sinh`, `Cosh`) recurse over the rest of the summands; the multiple-angle
   forms `Sin[n_Integer x_] /; n>1 :> Sin[(n-1)x] Cos[x] + Cos[(n-1)x] Sin[x]`
   reduce integer multiples down to `Sin[x]`/`Cos[x]`. Reciprocal heads
   (`Tan`/`Cot`/`Sec`/`Csc` and hyperbolic analogues) are rewritten as `Sin/Cos`
   ratios so the base rules apply. A large block of inverse-trig composition
   rules (`Cos[ArcSin[x]] :> Sqrt[1-x^2]`, etc.) is included.
2. **Expand** to distribute products of sums into a flat monomial sum.
3. **Pythagorean collapse.** If the expanded form has a denominator
   (`has_reciprocal_power`), no Pythagorean-eligible squared pair
   (`input_has_pythag_pair`), or more than `TRIG_FACTOR_ATOM_THRESHOLD` distinct
   squared trig atoms, only the direct-sum rules `trig_expand_pythag` are applied
   (`ReplaceRepeated`). Otherwise it first runs polynomial `Factor` — which turns
   `Sin[nx]^2 + Cos[nx]^2` into `(Sin[x]^2+Cos[x]^2)^n` — then applies
   `trig_expand_pythag` to collapse `(Sin^2+Cos^2)^n -> 1` (and the negated-sign
   and hyperbolic variants).
4. **Re-Expand** to restore the canonical monomial form.

**Data structures.** `trig_expand_rules` and `trig_expand_pythag` are static
`parse_expression`'d rule lists built in `trigsimp_init`. Threads over `List`
(via `ATTR_LISTABLE`) and over equations/inequalities/logic heads
(`trigexpand_threads_over`). Memoized through the active `FactorMemo` by the
`builtin_trigexpand` wrapper (`trig_memo_call`).
