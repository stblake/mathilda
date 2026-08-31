---
ticket: RG-1
created: 2026-08-30
type: plan
lifecycle: active
status: draft
full_plan: thoughts/shared/tickets/RG-1/plan.md
---

# RG-1 — `RandomGraph[{n, m}, k]`: one-page view

## Recommendation

Extract the single-graph body of `builtin_random_graph`
(`src/graph/generators.c:103-134`) into a static helper, make the builtin a thin dispatcher
on `arg_count`, and loop the helper `k` times into a `List`. `k` is validated with the
file's existing `as_count` (`:25-28`), which already gives the exact count semantics the
five sibling heads in `src/random.c` agree on: `k = 0` → `{}`, and negative, non-integer,
or symbolic `k` → unevaluated, silently.

Two phases: **(1)** the generator change, **(2)** tests, docstring, spec, changelog.
One source file of real change; the rest is docs and tests.

## Options considered

| Option | Verdict |
|---|---|
| Extract a helper, loop it `k` times (chosen) | Smallest diff that shares one code path between both arities. Keeps the current inline assembly — `make_graph` frees none of the arrays it is handed and takes an edge array, not a built `List`. |
| Generalize count-spec parsing into a shared helper | Rejected — no such helper exists; seven heads across three subsystems hand-roll their own. Would touch `src/random.c`, `src/ml/dist.c`, `src/list/constant_array.c`. |
| Replace the `evaluate(RandomSample[...])` round-trip with a C-level RNG call | Rejected — the round-trip is what makes `SeedRandom` work for free. Research resolved the O(n²) rework as out of scope. |
| Also fix the measured candidate-tree leak (Finding 4B) | Rejected — nearest call of the three. Research scoped it out; the human confirmed documenting over fixing. |

## Decision criteria

- Match the established convention rather than invent one: silent `NULL` for a bad count
  (zero `Message(` calls in `src/random.c`), `{Graph[...]}` for `k = 1`.
- Preserve `SeedRandom` reproducibility — non-negotiable, and free if the sampler call
  stays as-is.
- Do not let the refactor change the existing 1-arg form.

## Decisions

- **Flat `k` only.** The nested `RandomGraph[{n, m}, {k1, k2}]` form is out of scope, per
  the research's resolved question.
- **Extract a static helper, don't generalize.** There is no shared count-spec parser
  anywhere in the tree — `src/random.c` holds five hand-written copies. A sixth local check
  is the idiomatic choice here; a generalizing refactor would be scope creep.
- **Keep the `evaluate(RandomSample[...])` round-trip.** It is why `SeedRandom` works.
- **Fix the `n <= 1` case incidentally** (chosen by the human, 2026-08-30). The helper
  early-returns an edgeless graph when `maxe == 0` — and **only** then, not when `m == 0`,
  which works today and must keep consuming the RNG stream identically.
- **Document the leak, don't cap `k`** (chosen by the human, 2026-08-30). No sibling
  `Random*` head restricts its count, and a cap would be a non-Mathematica restriction.
- **Free the C arrays this diff allocates**, keeping the assembly inline rather than
  calling `make_graph`, which frees nothing. `expr_new_function` memcpys
  (`src/expr.c:257`), so the caller owns the array. Leak (A), fixed only where the diff
  lands.
- **Fail unevaluated instead of crashing on an absurd `n` or `k`.** Compute `maxe` in
  `unsigned long long`, bound it against `SIZE_MAX`, and check both `calloc`s — no
  arbitrary cap, just the head's existing `NULL` error channel.

## Non-goals

- **The candidate-tree leak, Finding 4B** (~364 KB/call at n=50). Scoped out by research;
  documented rather than fixed. It multiplies by `k`.
- **Leak (A) in the other four generators** (`CompleteGraph`, `CycleGraph`, both
  `PathGraph` branches). Same one-line shape, but out of this diff.
- **The nested `{k1, k2}` dimension form.**
- **Reworking the O(n²) candidate materialization.** Explicitly resolved as out of scope.
- **Any packed/NDArray/`Compile[]` surface.** `RandomGraph` returns a `Graph` object, not a
  machine array — the "non-machine object" exemption CLAUDE.md names. Verified absent from
  `src/pack.c`'s `AWARE` list, `src/ndkernels.c`, `src/ndinteger.c`, and `src/compile/`.

## Open Questions

### Unresolved

_None._

### Resolved

- [x] Should the `n <= 1` / `maxe == 0` case be fixed incidentally? — Yes. Confirmed with
  the human, 2026-08-30: the helper early-returns an edgeless graph, which the `k` form
  needs anyway and which removes a wrong answer rather than replicating it `k` times.
- [x] Cap `k` to bound the per-call leak, fix the leak, or document it? — Document only.
  Confirmed with the human, 2026-08-30. No sibling `Random*` head caps its count, and
  fixing the candidate-tree leak was already scoped out by the research.

## Requires Approval

The `n <= 1` fix (AC-14, AC-15) changes the observable behavior of the **existing** 1-arg
form: `RandomGraph[{0, 0}]` and `RandomGraph[{1, 0}]` go from unevaluated to returning an
edgeless graph. That is a bug fix, not a regression, but it is a behavior change outside
the ticket's literal scope and was approved in conversation on 2026-08-30.

## Architecture Impact

- New services introduced: none
- APIs changed: `RandomGraph` gains a second, optional argument — backward compatible yes
  (the 1-arg form is unchanged; the previously-unevaluated 2-arg form now evaluates)
- Data crossing a service boundary: none
- New external dependency: none
- Deployment topology change: none

## Subsystems & Dependencies

- Subsystems touched: graph (invocation: inline — no doc exists;
  `thoughts/shared/subsystems/` is absent from the tree)
- Interdependencies surfaced: graph ↔ random — `builtin_random_graph` reaches
  `RandomSample` through `evaluate()` rather than a C-level RNG call. This is load-bearing
  for `SeedRandom` and is also the cause of the `n <= 1` bug being fixed here.

## How we'll know it worked

19 acceptance criteria in the full plan, all REPL-observable, covering the `k` form, the
`k = 0` / `k = 1` boundary, invalid `k`, out-of-range `m` at both arities, seeded
determinism, draw independence, the unchanged 1-arg form, and the `n <= 1` fix. Plus a
build under `SDKROOT=$(xcrun --show-sdk-path) make -j8`, `make check-c99`, the graph unit
suite, and a `leaks --atExit` run.

## Review status

`plan-reviewer` ran on 2026-08-30 (scope-boundary + testability lenses): 2 BLOCKING,
2 WORTH FLAGGING, all four addressed in the full plan. The blocking pair were a
`make_graph`-vs-inline contradiction that would have reintroduced the leak the plan
commits to fixing, and an unapproved RNG-stream change to the existing 1-arg form at
`m == 0` — that early return is now `maxe == 0` only. Findings are transcribed verbatim
in the full plan's `## Plan Review`; `### Blocking` is clear.

## Appendix

Full detail, phase breakdown, code sketches, and the reviewer findings:
[`thoughts/shared/tickets/RG-1/plan.md`](plan.md).
Source research: [`thoughts/shared/tickets/RG-1/research.md`](research.md).
