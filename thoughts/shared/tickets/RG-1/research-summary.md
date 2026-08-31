---
created: 2026-08-30T21:18:27-0400
researcher: Michael Sollami
topic: "How does RandomGraph work today, and what would it take to add the RandomGraph[{n,m}, k] form that returns k random graphs?"
type: research
lifecycle: active
full_research: thoughts/shared/tickets/RG-1/research.md
---

# Research Summary: RandomGraph `k` form

**Full research (appendix)**: `thoughts/shared/tickets/RG-1/research.md`

## Recommendation

Add the flat `k` form directly in `src/graph/generators.c`: relax the `arg_count != 1`
gate at `:104`, extract the existing single-graph body into a static helper, and loop it
`k` times into a `List`. Parse `k` with the file's existing `as_count` (`:25-28`), which
already gives the codebase-standard behavior — `k = 0` → `{}`, negative/non-integer/
symbolic → unevaluated, no `Message[]`. `k = 1` must return `{Graph[...]}`, a one-element
list, not a bare graph.

## Options Considered

1. **Local loop over the existing path** *(recommended)* — ~15 lines, one file plus
   docstring/spec/test updates. Tradeoff: inherits the O(n²) candidate cost per graph and
   the per-call leak below, both multiplied by `k`.
2. **Rework sampling to O(m) first, then add `k`** — removes the n² floor and would let
   `k` scale. Tradeoff: touches the working 1-arg path and its `SeedRandom` guarantee;
   ruled out of scope.
3. **Generalize a shared count-spec parser across the `Random*` family** — would retire
   five duplicated validation loops in `src/random.c`. Tradeoff: a refactor of six
   unrelated heads riding along on a small feature; rejected.

## Decision Criteria

- Five sibling heads already establish the exact convention to copy, and there is **no**
  shared helper to reuse — so a sixth local check is idiomatic here, not sloppy.
- The `SeedRandom` guarantee comes from delegating to a real `RandomSample[...]` call
  through the evaluator; any rework risks it, and reproducibility is already tested.
- The packed/NDArray/`Compile[]` rule in CLAUDE.md does not bite: `RandomGraph` returns a
  `Graph`, a non-machine object, and appears in no `AWARE` list, ND kernel, or audit
  baseline.

## Open Questions

### Unresolved

_None._

### Resolved

- [x] Flat `k` only, or nested `{k1,k2}` too? — Flat only; nested is out of scope.
- [x] Rework the O(n²) candidate materialization? — No; reuse the existing path in a loop.
- [x] Prior attempt or deliberate deferral? — Neither; greenfield, single commit
  `56035303`.
- [x] Fix the memory leak found while measuring? — No; record as a finding only.

## Requires Approval

**A measured, pre-existing leak that scales with `k`.** `RandomGraph[{50,20}]` leaks
~364 KB *per call* — 2000 iterations leak 727 MB across 12.2M objects — because the whole
candidate-edge tree is never reclaimed (`generators.c:127`; `RandomSample` itself measures
0 leaks). All four generators additionally leak two `calloc`'d C arrays per call, because
`expr_new_function` memcpys rather than adopts (`src/expr.c:257`); `CompleteGraph[50]`
× 2000 = 21 MB.

Out of scope by decision — **finding only, not to be fixed here.** But the agreed
loop-`k`-times approach multiplies it by `k`: `RandomGraph[{50,20}, 1000]` would leak
~364 MB. Worth a changelog note or a docstring caveat rather than silence.

Two smaller pre-existing items, also not in scope: `RandomGraph[{0,0}]` and `{1,0}` are
unevaluated (empty candidate list hits `is_nonempty_list` in `RandomSample`), and
`RandomGraph` has no recorded `EXEMPT` entry in the packed-array audit tooling even though
it legitimately qualifies.
