---
references:
  - "Joel S. Cohen, *Computer Algebra and Symbolic Computation: Mathematical Methods* (A K Peters, 2003)."
source: src/simp/simp.c
---
**Algorithm.** Simplify is a complexity-weighted, memoized candidate-set search
over the existing battery of algebraic transforms. `builtin_simplify`
(`simp_builtins.c`) parses options — a positional or `Assumptions ->`
assumption (combined into the `$Assumptions` default with `And`), a
`ComplexityFunction`, and `TransformationFunctions` — builds an `AssumeCtx` of
normalised facts, then drives the search. Inexact inputs are rationalised on
entry and numericalised on exit (`internal_rationalize_then_numericalize`).
`Equal`/`Less`/`And`/... heads are threaded manually (with a relational
rebalancing candidate); `List` is threaded via `ATTR_LISTABLE`.

The core driver is `simp_bottomup` (`simp_bottomup.c`), which descends into
every `Plus`/`Times`/`Power` child (memoizing each result in a `SimpMemo` hash
table) and dispatches each node through `simp_dispatch` → `simp_search`.
`simp_classify` routes pure-polynomial/rational shapes to dedicated pipelines
(`simp_pipeline_polynomial`/`_rational`/`_logexp`) and there are top-level fast
paths: a `SHAPE_RATIONAL` shortcut that runs `Together`/`Cancel`/`Factor` once
at the top, an algebraic-tower `Together[expr, Extension -> Automatic]` collapse,
and a `simp_trig_rational` substitution that maps trig/opaque subtrees to ground
symbols and works in the quotient ring.

`simp_search` is the heuristic engine. It seeds a `CandSet` with the input plus
the output of a long list of correctness-preserving rewriters (assumption rules,
log/exp identities, `SimpLogRules`, trig roundtrip, ExpToTrig, Pythagorean
square-completion / reduction / canonicalisation, trig-at-rational-Pi, tan-
addition, half-angle, radical product combine, sqrt/cube-root/algebraic
denesting, roots-of-unity, factorial decomposition, per-variable Collect). It
then runs `SIMP_ROUNDS` (= 2) rounds in which every seed is fed through the
`SIMP_TRANSFORMS` table — `Together, Cancel, Expand, ExpandNumerator,
ExpandDenominator, Factor, FactorSquareFree, FactorTerms, Apart, TrigExpand,
TrigFactor, TrigReduce, TrigToExp` — plus chained Pythagorean/radical/trig-Pi
passes. Each transform call is gated by `transform_can_fire` (cheap precondition
check) and the running best is updated by `update_best` against the complexity
score; candidates strictly worse than their parent are dropped from the next
round (with a loosened `2*parent + 8` bound for TrigExpand and the seed-phase
blow-up guard). A short-circuit (`simp_best_is_zero`) exits as soon as a literal
`0` is reached. Final polish passes apply `simp_lift_common_factor`,
`transform_pythag_reduce`, and `canon_negate_pairs`.

**Complexity scorer.** Candidates are ranked by `score_with_func`: with no user
`ComplexityFunction` this is `simp_default_complexity` (the SimplifyCount metric
— LeafCount with integers contributing their decimal-digit count) plus a
`nested_radical_penalty` (+3 per truly-nested radical). A user function is
evaluated as `f[candidate]`; an Integer result is used directly, a BigInt scores
`SIMP_SCORE_INF`, anything else falls back to the default. Lower wins; ties
favour the form that reached the score plateau first (force-take semantics let
assumption/log/factorial rewrites win even on a tie).

**Data structures.** `CandSet` (dynamic `Expr*` array, dedup via `expr_eq`, capped
at `SIMP_CAND_CAP` = 12); `SimpMemo` (256-bucket chained hash of input→best for
the bottom-up driver); `AssumeCtx` (flat fact array from `assume_ctx_from_expr`);
a per-Simplify `FactorMemo` shared by Factor/Trig* so duplicate subexpression
work hits the cache. Each transform is invoked as `f[candidate]` through the real
evaluator (`traced_call_unary`), so Simplify composes the same builtins exposed
in the REPL.

**Limits.** Not a complete decision procedure: no real inequality reasoner
(`Simplify[x>0, x>0]` folds only by literal fact match), and assumption-driven
wins depend on the structural provers in `simp_assume.c`.
