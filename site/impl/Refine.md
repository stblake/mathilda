---
source: src/refine.c
---
**Algorithm.** `builtin_refine` is a thin orchestrator over the assumption engine
`Simplify` already uses (`AssumeCtx`, `apply_assumption_rules`, declared in `simp.h`,
which reserves that type "so future modules (Refine, …) can share it"). It first splits
trailing options from positional arguments and forms the *effective* assumption with the
identical policy to `PossibleZeroQ`/`Simplify`: an `Assumptions -> X` option replaces
`$Assumptions`, a positional assumption is `And`-conjoined with `$Assumptions`
(`read_dollar_assumptions`), and the result is evaluated so `And[True, p]` canonicalises to
`p`. That expression is parsed into an `AssumeCtx` fact set by `assume_ctx_from_expr`. From
there two paths run. **(1) Predicate head** (`Element`, or a relational/logical head —
`Equal`, `Less`, `And`, …): `refine_decide_predicate` tries to settle it to `True`/`False`.
`Element[x, dom]` goes to the assumption-aware `element_decide`; `Equal`/`Unequal` first
run `zero_test_decide_assuming` on `lhs − rhs` (and an equality-substitution via
`apply_assumption_rules`), then fall back to `Reduce` over the complexes; inequalities and
logic go to a `Reduce`/CAD entailment over the reals, where `P` is `True` iff `Reduce[A &&
!P]` is unsatisfiable (literally `False`) and `False` iff `Reduce[A && P]` is. **(2) Rewrite
path** (everything else): `apply_assumption_rules` rewrites the expression under the facts
(the very rules `Simplify` applies — `Sqrt[x^2]->±x/Abs[x]`, `Log[x^p]->p Log[x]`, integer-`k`
trig, `Abs`/`Sign`/`Conjugate` under sign facts, …), then a bottom-up `deep_positivity_walk`
resolves `Sign[p]`, `Abs[p]` and `Sqrt[p^2]` for compound polynomials `p` the fast sign
prover could not settle, each via a `deep_sign` CAD run, and the tree is evaluated. With no
usable facts the call is the identity (it steals the positional expression out of `res`).

**Data structures.** Everything is `Expr` trees. Assumptions are held as an `AssumeCtx`
(flat `Expr*` fact array) borrowed from `simp`. Entailment is threaded through a
`RefineBudget { double deadline; int calls_left; }`; each `Reduce` query is built as
`Reduce[stmt, varlist, domain]`, where `collect_bare_vars` gathers the distinct bare-symbol
leaves (skipping protected real constants such as `Pi`, and *not* reusing the polynomial
`collect_variables`, which would atomise `x^p`/`Floor[x]` into pseudo-variables `Reduce`
rejects).

**Complexity / limits.** Dominated by the `Reduce`/CAD entailment, which is doubly
exponential in the number of variables; it is bounded two ways — the variable set is capped
at 6 (wider statements decline) and at most `REFINE_MAX_ENTAILMENT_CALLS` (24) `Reduce`
calls run per `Refine`. `TimeConstraint` (default 30 s) is a cooperative wall-clock budget
checked *before* each entailment via `simp_mono_seconds`; a running `Reduce` is never
preempted (CAD has no interior abort hook), and the async, malloc-lock-prone
`TimeConstrained[]` is deliberately avoided. Predicates it cannot decide, and expressions no
rule rewrites, pass through unchanged. `Refine` is purely symbolic/structural: it carries no
NDArray/packed kernel and no `Compile[]` lowering, by design, since it returns symbolic
expressions rather than machine numbers.
