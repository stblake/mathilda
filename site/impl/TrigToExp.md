---
source: src/simp/trigsimp.c
---
**Algorithm.** `builtin_trigtoexp` rewrites circular, hyperbolic, and their
inverse functions into complex-exponential / logarithmic form. It applies the
`trig_to_exp_rules` rule list (a `ReplaceAll`) — e.g. `Sin[x] :> I/2 E^(-I x) -
I/2 E^(I x)`, `Cos[x] :> (E^(-I x) + E^(I x))/2`, `Tan`/`Cot`/`Sec`/`Csc`,
the `Sinh`/`Cosh`/`Tanh`/... hyperbolic analogues, and the `ArcSin :> -I Log[...]`
family of inverse identities — then runs `Expand` to distribute the result into
a flat sum. The `Times`/`Power`-level trig canonicalizer is suppressed for the
duration (`trig_canon_suppress_inc`/`dec`) so the intermediate `Sin/Cos` forms
are not prematurely re-collapsed back to `Tan` etc. before the rule fires.

**Data structures.** The rules are a single `parse_expression`'d `List` of
`RuleDelayed`s, built once in `trigsimp_init` and stored in the static
`trig_to_exp_rules`. Results are memoized through the active `FactorMemo`
(keyed under a `TrigToExp[arg]` head) so repeated Simplify candidate-set calls
on identical subexpressions hit the cache instead of re-running `ReplaceAll +
Expand`, which costs 30–130 ms on Tan-rich inputs.
