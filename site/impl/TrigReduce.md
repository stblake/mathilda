---
source: src/simp/trigsimp.c
---
**Algorithm.** `builtin_trigreduce_impl` is the product/power-direction inverse of
TrigExpand: it rewrites products and integer powers of single-argument circular
and hyperbolic trig calls into single trig calls of compound (sum / multiple)
arguments, using the classical product-to-sum and power-reduction identities

```
Sin[a] Cos[b] = (Sin[a+b] + Sin[a-b]) / 2     Cos[a] Cos[b] = (Cos[a+b] + Cos[a-b]) / 2
Sin[a] Sin[b] = (Cos[a-b] - Cos[a+b]) / 2     Sin[x]^2 = (1 - Cos[2x]) / 2
Cos[x]^2 = (1 + Cos[2x]) / 2                   (and the Sinh/Cosh hyperbolic analogues)
```

The pipeline (trig canonicalizer suppressed throughout so the `Sin/Cos`
intermediate forms are not re-collapsed before the rules fire):

1. **To Sin/Cos.** `ReplaceRepeated` with `trig_factor_to_sincos` rewrites
   reciprocal heads (`Tan`/`Cot`/`Sec`/`Csc` and hyperbolic) as `Sin/Cos` ratios
   so the product-to-sum rules can see them.
2. **Iterate to a fixed point** (bounded at 16 iterations): alternate
   `ReplaceRepeated` with `trig_reduce_rules` (the power-reduction and
   product-to-sum identities, with each constructed compound argument wrapped in
   `Expand[...]` so e.g. `Sin[(x+y)-(x-y)]` canonicalizes to `Sin[2y]` before the
   surrounding trig head sees it) and an `Expand` step. The iteration is required
   because `Expand` re-exposes `Cos[2x]^2` terms hidden inside `(1 - Cos[2x])^2/4`
   after a power-reduction pass on `Sin[x]^4`, which the rule then reduces again.
   The 16-iteration cap covers exponents through `Sin[x]^65536`.
3. **Together** to combine over a common denominator so numerators appear as a
   single `Plus` for the collapse rules.
4. **Angle-addition collapse.** `ReplaceRepeated` with `trig_reduce_collapse`:
   coefficient-aware reverse angle-addition (`c. Sin[a]Cos[b] + c. Cos[a]Sin[b]
   :> c Sin[a+b]`, etc.) plus negative-argument cancellation rules guarded by
   `SameQ[Expand[a+b], 0]` that fold `Sin[a-b] + Sin[b-a]` (which the auto-evaluator
   leaves un-reduced) to zero.
5. **From Sin/Cos.** `ReplaceRepeated` with `trig_factor_from_sincos` restores
   `Tan`/`Sec`/`Csc` (and hyperbolic) where the ratio/reciprocal shape survives.
6. **Final canonicalisation:** `Expand` then `Together`, distributing outer scalars
   (`1/2 (2 Cos[a+b] + 2 Sin[a+b])` flattens) while keeping irreducible fractions
   like `(3 - 4 Cos[2x] + Cos[4x])/2` as a single rational.

**Data structures.** Four static rule lists (`trig_factor_to_sincos`,
`trig_reduce_rules`, `trig_reduce_collapse`, `trig_factor_from_sincos`) parsed in
`trigsimp_init`. `Times` is `ATTR_ORDERLESS`, so the matcher commutes factors and
only one direction of each product-to-sum pair needs to be written. Threads over
`List` (via `ATTR_LISTABLE`) and over equation/inequality/logic heads; memoized
through the active `FactorMemo` via the `builtin_trigreduce` wrapper.
