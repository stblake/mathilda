---
ticket: RG-2
created: 2026-08-31
type: plan
lifecycle: active
status: draft
full_plan: thoughts/shared/tickets/RG-2/plan.md
---

# Plan Summary: `FindVertexColoring[g]` (RG-2)

**Full plan (appendix)**: `thoughts/shared/tickets/RG-2/plan.md`

## Recommendation

Add `FindVertexColoring[g]` as a new translation unit, `src/graph/vertexcoloring.c`, in four
phases. It returns a **minimal** colouring — integers in `VertexList` order — computed by
exact backtracking bracketed between a greedy-clique lower bound and a DSATUR upper bound.
Graphs above 128 vertices return unevaluated. The result is a packed int64 `List` with a
plain-`List` fallback. No existing head changes; `src/core.c` and the makefile need no edits.

This is the first combinatorial-search head in `src/graph/`: the subsystem's 27 existing
heads cover representation, traversal and connectivity only.

## Options Considered

1. **DSATUR only** — fast and simple, but returns valid-but-frequently-non-minimal
   colourings. Rejected: the failure is silent, a plausible integer list that contradicts
   the documented meaning of the function. A wrong answer is worse than a refusal.
2. **Exact search, ascending `k = 1..ub`** — what the first draft specified. Rejected after
   review: with no lower bound, `CompleteGraph[128]` sits *under* the cap and must refute
   k=1…127 before answering, which is a hang with no abort channel to escape it.
3. **Exact search over `lb..ub` with a clique lower bound** — chosen. `lb == ub` short-circuits
   the worst case to zero search steps.
4. **A time-based abort instead of a vertex cap** — rejected: no timeout facility exists,
   building one is larger than this feature, and a cap is deterministic where a timeout makes
   results machine-dependent.

## Decision Criteria

- **Minimality is the documented contract**, so an approximation that looks right is the
  worst outcome available.
- **The subsystem has no abort channel** — `return NULL` is the only error path — so
  tractability has to come from bounds, not from giving up mid-search.
- **The cap bounds size, not cost.** That distinction is why the lower bound matters more
  than the cap: size alone does not stop `CompleteGraph[128]`.
- **Exact backtracking seeded by DSATUR *is* Wolfram's `"BacktrackingDS"`**, so shipping
  only it is a documented subset rather than a divergence.
- **`AWARE` is consumer-only** (`src/pack.c:460-488`), so the packed producer path needs no
  audit registration — verified, not assumed.

## Decisions

- **Exact, not DSATUR-only.** A non-minimal result silently contradicts the documented
  semantics; a refusal does not.
- **DSATUR seeds the upper bound**, it is not the answer — cheap, and it prunes hard.
- **A greedy-clique lower bound, searching `lb..ub` not `1..ub`.** Without it
  `CompleteGraph[128]` — under the cap, so accepted — refutes k=1…127 before answering,
  which is a hang. With it `lb == ub` and the answer is immediate. This is what makes the
  cap safe rather than nominal.
- **`FVC_MAX_VERTICES = 128`, refusing above it** (human, 2026-08-31). No timeout or abort
  channel exists, so a hard cap is the only honest guard. Justified on *typical* not worst
  case. Named-cap precedent: `FM_MAX_CON` (`src/solve/reduce_fm.c:18`).
- **Form 1 only** (human). Forms 2–3 become a later mapping layer, keeping the packed
  obligation on one branch.
- **No `Method` option** (human). One algorithm that *is* a named Wolfram method beats a
  parameter with one legal value.
- **Phase 1 registers nothing** (human). The head appears only once the search is exact, so
  no build answers non-minimally.
- **Packed `List`, never visible `NDArray[...]`.** Unanimous producer precedent
  (`src/list/range.c:143-149`).

## Non-goals

- **`FindVertexColoring[g, {c1, ...}]` and `FindVertexColoring[g, l]`.** A later mapping
  layer over form 1: compute the integer index vector once, then substitute the caller's
  objects. Those objects are arbitrary expressions and cannot pack, which is exactly why
  keeping them out keeps the packed obligation on a single branch.
- **`Method -> "ILP"`** — no integer-programming solver exists in-tree.
- **`Method -> "HybridEA"`** — an evolutionary metaheuristic, a separate body of work.
- **`PerformanceGoal`.**
- **A timeout or abort mechanism for the subsystem.** The cap is a local guard, not the
  general facility `src/graph/` lacks.
- **Fixing the eight backwards-reading "frees res" header comments** (RG-2 research § 3).
- **`ChromaticNumber`**, though `Max` of this result is exactly that.

## The four phases

```
Phase 1  Algorithm, unregistered — cap, DSATUR bound, clique bound, exact search.
         Helpers non-static with graph.h prototypes (registration still deferred).
         Direct C unit tests; needs #include "graph.h" in test_graph.c.
Phase 2  Registration — the 7 sites incl. the 3 SYM_* sites; AC-1..AC-18 as
         assert_eval_eq rows.
Phase 3  Packed producer path — ndbuild_open_i64 open/fill + plain-List fallback;
         AC-19/AC-20 as C-level is_packed_list assertions.
Phase 4  Docs — graphs.md entry, changelog 2026-08-31.md.
```

## Open Questions

### Unresolved

_None._

### Resolved

- [x] Which call forms? — Form 1 only. _(stated by the human)_
- [x] Which `Method` values? — None; exact backtracking seeded by DSATUR is Wolfram's
  `"BacktrackingDS"`. _(stated by the human)_
- [x] NP-hardness without a timeout channel? — A hard vertex cap, named constant with its
  reasoning. _(stated by the human)_
- [x] Cap value? — 128, on typical not worst case. _(picked "128" over 64, 256, density-based)_
- [x] Does Phase 1 register a non-minimal head? — No; unregistered, registered in Phase 2.
  _(picked "Split, but gate Phase 1 behind no registration")_
- [x] Does a packed producer need an `AWARE` entry? — No; `AWARE` governs consumed
  arguments only. _(resolved by reading `src/pack.c:460-488`)_

## Requires Approval

Two items. **The cap is a user-visible refusal**: `FindVertexColoring` on a 129-vertex graph
returns unevaluated where Wolfram answers. That is a deliberate, documented divergence and
needs sign-off as such. **Phase 1 introduces direct C unit tests to `tests/test_graph.c`**,
which today tests exclusively through `assert_eval_eq`; an unregistered head cannot be
reached through the evaluator, so this follows from the phase split rather than being a free
choice.

## Architecture Impact

- New services introduced: none
- APIs changed: none — one new head, no existing signature touched
- Data crossing a service boundary: none
- New external dependency: none
- Deployment topology change: none

Routine tier. `architecture-guidance` resolves to `docs/design`, which holds no
vertex-colouring or graph-subsystem decision record, so there is nothing org-specific to
cite.

## Subsystems & Dependencies

- Subsystems touched: `graph` (invocation: inline — no subsystem doc exists at
  `thoughts/shared/subsystems/graph.md`)
- Interdependencies surfaced: `graph` ↔ `pack` — a new producer-direction use of
  `ndbuild_open_i64`, the subsystem's first

## Review outcome

`plan-reviewer` returned three BLOCKING and four WORTH FLAGGING findings; all three BLOCKING
are resolved and moved, with the full text in the long plan's `## Plan Review`. Two were real
design defects rather than wording:

- **AC-8 could not detect the bug it was written for.** On the path `c–a–b`, order `{c,a,b}`
  colours `{1,2,1}` and sorted `{a,b,c}` colours `{1,2,2}` — so "first two differ" is `True`
  under both. Strengthened to `col[[1]] === col[[3]] && col[[1]] =!= col[[2]]`.
- **The search had no lower bound**, so `CompleteGraph[128]` would hang under the cap; and
  the only timing witness, `CycleGraph[128]`, is an even cycle where DSATUR gives ub=2 and
  no search happens. Added the clique lower bound and three cap-boundary rows including a
  genuinely searched dense instance.
- **Phase 1's C unit tests were not mechanically achievable** — `static` helpers, no header
  prototype until Phase 2, and `test_graph.c` has no `#include "graph.h"`.

One residual stated honestly in the long plan: a dense graph where neither bound is tight can
still be slow, and AC-10b pins a budget for one instance rather than proving a general bound.
