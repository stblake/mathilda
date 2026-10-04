---
source: src/solve/reduce_companions.c
---
**Algorithm.** `builtin_find_instance` peels trailing options (`Modulus` is
honoured; `Method` / `WorkingPrecision` / `RandomSeeding` are accepted and
ignored, since the search is exact and deterministic), reads the optional domain
and witness count `n`, then runs `fi_run_search` under a message mute — its
internal `Reduce` / `Solve` / `NMinimize` probes are speculative and, as in
Mathematica, must not leak diagnostics. The strategy is **soundness-first and
verify-gated**: every candidate point is checked against the original statement
(`expr /. point === True`, by `fi_verify`) before being accepted, so a returned
instance is never wrong, and a point that cannot be verified is discarded rather
than reported.

`fi_run_search` is a cascade of witness sources, tried cheapest/most-exact first.
`Reduce` is the satisfiability-and-solution-set oracle (step 1): its `False` is
taken as `{}` only for exactly-decidable systems — for transcendental or inexact
systems (`fi_is_transc_inexact`) that `False` is distrusted, and indexed variables
`c[i]` skip the oracle since `Reduce` rejects them — otherwise each top-level `Or`
clause is sampled into a point. Then `Solve` fallbacks with generated-parameter
instantiation reach parametric Diophantine families (the Pell equation); a bounded
integer box search runs over `Integers`; an equations-only retry, a 1-variable real
transcendental root bracket, a solve-one-then-sample grid, and Rabinowitsch ideal
saturation each cover systems the oracles decline; a Gröbner emptiness certificate
proves `{}` for a declined polynomial system; structured exact candidate sampling
finds branch-cut and open-region witnesses; and numerical feasibility
(`NMinimize` / least-infeasibility) is the last resort for transcendental/inexact
`Real` systems. Interval samples come from `rru_rational_between`. The `Booleans`
domain is handled separately by `fi_boolean`, reusing the `LogicalExpand` DNF
engine for satisfiability.

**Data structures.** Witnesses accumulate in a `FiWit` buffer and are returned in
`Solve`'s rule-list form `{{x -> v, ...}, ...}`. A "variable" is matched
structurally with `expr_eq`, so a plain symbol (`x`) and an indexed form (`c[i]`)
are both accepted. Everything is `Expr` trees driven through `evaluate`.

**Complexity / limits.** The default domain is `Complexes`, or `Reals` when the
statement carries an ordering (as in `Reduce`). `{}` is returned only when the set
is *provably* empty — a `Reduce` `False` on an exactly-decidable system, an
exhausted finite integer box, or a Gröbner certificate; when no witness is found
and emptiness is not proved, the call stays unevaluated. Cost is dominated by the
`Reduce`/`Solve`/Gröbner probes the cascade invokes; `Modulus -> p` searches over
`Z/pZ`.
