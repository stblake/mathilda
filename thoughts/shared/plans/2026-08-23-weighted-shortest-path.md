---
created: 2026-08-23
source_sha: 3d872247
subsystems: [graph]
type: plan
lifecycle: active
status: implemented
---

# Weighted Shortest Path Implementation Plan

## TL;DR
Makes `FindShortestPath[g,s,t]` and `GraphDistance[g,s,t]` weight-aware: when `g` carries a
non-negative-numeric `EdgeWeight`, both dispatch to Dijkstra instead of BFS, matching real
Wolfram Language semantics and completing the follow-up ticket 1's own Non-goals named.
Falls back to the existing unweighted BFS whenever weights are absent, symbolic, or contain
a negative value — never regresses a previously-working call. Verified by extending
`tests/test_graph.c` and updating the one existing regression assertion (from ticket 1) that
this ticket intentionally supersedes.

## Overview
Ticket 1 added `EdgeWeight`/`WeightedAdjacencyMatrix` but explicitly deferred making
`FindShortestPath`/`GraphDistance` weight-aware (its own Non-goals: "real algorithmic scope
growth, not a few hours"). This ticket delivers that follow-up. `GraphAdj`
(`graph_util.c`), the structure both builtins currently use, has no weight storage and is
shared by 5 other builtins (ConnectedComponents, WeaklyConnectedComponents, FindSpanningTree, ConnectedGraphQ, VertexConnectivity) — so this plan builds a separate, call-scoped weighted adjacency
inside `shortestpath.c` rather than widening `GraphAdj` itself, directly applying the lesson
from ticket 1's `plan-reviewer`-caught defect (a shared choke point is a bigger blast radius
than it looks).

## Decisions
- **Dispatch is automatic, not a new builtin name** — real Wolfram Language uses edge weights
  automatically when present; matching that beats inventing `WeightedFindShortestPath`.
- **Dijkstra fires only for non-negative numeric weights; anything else falls back to BFS**
  unchanged — Dijkstra is incorrect on negative weights, and a previously-working call must
  not regress to unevaluated just because a weight is symbolic or negative.
- **No change to `GraphAdj`/`graph_build_adj`** — a local, per-call weighted structure lives
  only in the two changed builtins, reusing ticket 1's `graph_resolve_edge_weights`, keeping
  blast radius off the 5 other builtins sharing that structure.
- **Plain O(V²) Dijkstra, no priority queue** — matches this codebase's existing complexity
  tolerance for small-graph exact algorithms (`VertexConnectivity`'s own "brute-force ...
  intended for small graphs").
- **Ticket 1's `test_edge_weights` AC-11 lines for these two builtins must change** — they
  assert "weights ignored," the exact behavior this ticket removes. The other 5 builtins in
  that test are unaffected.

## Non-goals
- No Bellman-Ford or any negative-weight support — falls back to BFS instead.
- No priority-queue-based Dijkstra (O((V+E) log V)) — O(V²) is sufficient at this codebase's
  stated scale tolerance.
- No change to any other `graph_build_adj`-routed builtin (`ConnectedComponents`,
  `VertexConnectivity`, etc.) — weights are not meaningful for those algorithms.
- No A*, bidirectional search, or all-pairs shortest path.

## Acceptance Criteria

| ID | Given | When | Then | Input | Expected |
|---|---|---|---|---|---|
| AC-1 | A weighted graph where the shortest-hop-count path isn't the min-weight path | `FindShortestPath[g,s,t]` is called | it returns the min-weight path, not the min-hop path | `FindShortestPath[Graph[{1,2,3,4},{1->2,2->3,3->4,1->4},EdgeWeight->{1,1,1,10}],1,4]` | `{1, 2, 3, 4}` (weight 3) not `{1, 4}` (weight 10) |
| AC-2 | Same weighted graph | `GraphDistance[g,s,t]` is called | it returns the min total weight, not hop count | `GraphDistance[Graph[{1,2,3,4},{1->2,2->3,3->4,1->4},EdgeWeight->{1,1,1,10}],1,4]` | `3` |
| AC-3 | An unweighted graph | `FindShortestPath`/`GraphDistance` are called | behavior is byte-identical to before this ticket (BFS) | `FindShortestPath[CycleGraph[6],1,4]` | unchanged from pre-ticket output |
| AC-4 | A graph with a symbolic weight | `FindShortestPath[g,s,t]` is called | falls back to unweighted BFS, does not error | `FindShortestPath[Graph[{1,2,3},{1->2,2->3},EdgeWeight->{a,7}],1,3]` | `{1, 2, 3}` (BFS fallback, not unevaluated) |
| AC-5 | A graph with a negative weight | `GraphDistance[g,s,t]` is called | falls back to unweighted BFS | `GraphDistance[Graph[{1,2,3},{1->2,2->3},EdgeWeight->{-1,7}],1,3]` | `2` (BFS hop count, not a Dijkstra artifact) |
| AC-6 | An undirected weighted graph | `FindShortestPath[g,s,t]` is called | weights apply symmetrically | `FindShortestPath[Graph[{1,2,3},{1<->2,2<->3},EdgeWeight->{1,1}],1,3]` | `{1, 2, 3}` |
| AC-7 | An unreachable target on a weighted graph | `FindShortestPath`/`GraphDistance` are called | unreachable semantics unchanged | `FindShortestPath[Graph[{1,2,3},{1->2},EdgeWeight->{5}],1,3]` | `{}`; `GraphDistance[...]` → `Infinity` |

## Open Questions

### Unresolved
_None._

### Resolved
- [x] Dispatch mechanism, fallback rule, `GraphAdj` scope — see Decisions (from research).

## Plan Review

### Blocking
_None._

### Worth Flagging
_None._

### Resolved
**[BLOCKING] `double dist[]` has no stated path back to an exact `Expr`, so AC-2 is not
achievable as specified**
- Where: was in Phase 1 §2; `EXPR_REAL` prints distinctly from `EXPR_INTEGER`
  (`src/print.c`), and `assert_eval_eq` does exact string comparison, so a raw
  `expr_new_real(dist[it])` would print `12.` against AC-2's exact `3`/`12`.
- How addressed: `double dist[]` is now stated as internal-comparison-only; the actual
  returned value is reconstructed exactly via `evaluate(Plus[w1,...,wk])` over the real
  `Expr*` weights along the discovered path (Phase 1 §2, rewritten).

**[BLOCKING] The claim that both builtins' AC-11 lines in `test_edge_weights` need to change
was false for `FindShortestPath`**
- Where: was in Decisions/Components & Files Affected; verified against
  `tests/test_graph.c:350-359` directly — the AC-11 test graph
  (`Graph[{1,2,3},{1->2,2->3},EdgeWeight->{5,7}]`) has exactly one path from 1 to 3, so BFS
  and Dijkstra agree on `FindShortestPath`'s result; only `GraphDistance`'s hop-count `"2"`
  needs to become the weighted total `"12"`.
- How addressed: corrected everywhere this was claimed — only `GraphDistance`'s AC-11
  assertion changes.

**[WORTH FLAGGING] The numeric-type list for `graph_weights_usable` omitted `EXPR_MPFR` (a
live leaf type under the default `USE_MPFR ?= 1` build) and reinvented rather than reused
`expr_is_numeric_like` (`src/expr.c:412`), which already covers it**
- How addressed: `graph_weights_usable` now explicitly reuses `expr_is_numeric_like` (minus
  `Complex`) instead of a hand-rolled type list.

**[WORTH FLAGGING] "6 other builtins" / "7 call sites across 5 files" were both wrong —
`StronglyConnectedComponents` does not call `graph_build_adj` at all (own Tarjan structure);
the real numbers are 5 builtins / 6 call sites / 4 files. Same miscount was inherited from
ticket 1's (already-shipped) plan.**
- How addressed: corrected everywhere in this plan; ticket 1's already-merged plan/tests are
  left as historical record (not retroactively edited — see journal `KIT-FEEDBACK-GRAPH.md`
  for the discussion of whether to fix it there).

## Requires Approval
_None._ — this is the named follow-up from ticket 1's own Non-goals; no new scope call.

## Architecture Impact
- New services introduced: none
- APIs changed: none
- Data crossing a service boundary: none
- New external dependency: none
- Deployment topology change: none

(Behavior of two existing builtins changes for weighted graphs only, per ticket 1's own
Non-goals naming this as the intended follow-up — not an API contract break, since an
unweighted graph's behavior is provably unchanged, AC-3.)

## Subsystems & Dependencies
- Subsystems touched: graph (invocation: inline)
- Interdependencies surfaced: none

## Risks and Rollback
_None — standard tier, no architectural impact._

---

## Current State Analysis
- `src/graph/shortestpath.c` (full file, see research) — `bfs()`, `resolve()`, and both
  builtins, no weight awareness.
- `src/graph/graph.h:112-120` — `GraphAdj` has no weight field; `graph_build_adj` has 6 call
  sites across 4 files (`components.c` x2, `connectivity.c` x2, `spanningtree.c`,
  `shortestpath.c`), used by 5 other builtins (`ConnectedComponents`,
  `WeaklyConnectedComponents`, `FindSpanningTree`, `ConnectedGraphQ`, `VertexConnectivity`) —
  `StronglyConnectedComponents` builds its own Tarjan-specific structure and does not call
  `graph_build_adj`.
- `src/graph/graph_util.c` — `graph_resolve_edge_weights(g)` (ticket 1) already gives the
  per-edge weight list in `EdgeList` order.
- `tests/test_graph.c`'s `test_edge_weights` AC-11 block asserts the exact "ignore weights"
  behavior this ticket removes for `FindShortestPath`/`GraphDistance` specifically.

## Desired End State
Both builtins dispatch to Dijkstra on a non-negative-numeric-weighted graph, BFS otherwise;
`tests/test_graph.c` reflects the new behavior; `docs/spec/builtins/graphs.md` and this
week's changelog are updated. Verify via the Acceptance Criteria table above run against the
live REPL, plus `make check-c99`, `make check-packed-aware`, and the full `graph_tests`
suite.

### Key Discoveries:
- Weight resolution is already solved by ticket 1's `graph_resolve_edge_weights` — this
  ticket is purely the Dijkstra algorithm plus the numeric/non-negative gate, not new
  weight-plumbing.

## Components & Files Affected

| File | Change |
|---|---|
| `src/graph/shortestpath.c` | Add a local weighted-adjacency builder + O(V²) Dijkstra (`double dist[]` for internal vertex-selection comparisons only) + exact-value reconstruction (see below); both builtins check for a usable `EdgeWeight` (via a new `graph_weights_usable(g)` helper) and dispatch to Dijkstra or the existing BFS accordingly |
| `src/graph/graph_util.c`, `src/graph/graph.h` | New helper `graph_weights_usable(const Expr* g)`: `true` iff `g` has a 3-arg `EdgeWeight` and every weight satisfies `expr_is_numeric_like(w)` (the codebase's own existing generic numeric-type check — `src/expr.c:412`, already covers `EXPR_INTEGER`/`EXPR_BIGINT`/`EXPR_REAL`/`EXPR_MPFR`/`Rational`), is not `Complex`, and is `>= 0` |
| `tests/test_graph.c` | Update `test_edge_weights`'s `GraphDistance` AC-11 line only (`"2"` → `"12"`, the weighted total `5+7`) — `FindShortestPath`'s AC-11 assertion (`{1, 2, 3}`) is unaffected, since that specific test graph has only one path from vertex 1 to vertex 3, so BFS and Dijkstra necessarily agree on it. Add a new `test_weighted_shortest_path` for AC-1 through AC-7 (which do exercise multi-path graphs) |
| `docs/spec/builtins/graphs.md` | Update the `FindShortestPath`/`GraphDistance` bullets and remove the "remain unweighted BFS ... documented future extension" note added by ticket 1 |
| `docs/spec/changelog/2026-08-17.md` | New entry |

## Core Flow Diagram

```mermaid
flowchart TD
    A["FindShortestPath[g,s,t] / GraphDistance[g,s,t]"] --> B{graph_weights_usable(g)?}
    B -->|no: unweighted, symbolic, or has a negative weight| C[existing BFS path, unchanged]
    B -->|yes: EdgeWeight present, all non-negative numeric| D[build local weighted adjacency]
    D --> E[O(V^2) Dijkstra from s]
    E --> F[reconstruct path / distance to t]
    C --> G[return]
    F --> G[return]
```

## Alternatives Considered

### Extend `GraphAdj` with a weight array
**Rejected because:** `graph_build_adj` is shared by 6 call sites across 4 files
(`components.c`, `connectivity.c`, `spanningtree.c`, `shortestpath.c`) — widening it is
exactly the kind of shared-choke-point risk ticket 1's `plan-reviewer` pass caught as a real
defect. A local structure confined to `shortestpath.c` has zero blast radius on the other 6
builtins.

### A separate `WeightedFindShortestPath` builtin
**Rejected because:** real Wolfram Language dispatches on graph properties automatically,
not by function name — matching that is both more faithful and avoids two names for what a
user thinks of as one operation.

## Implementation Approach
Add `graph_weights_usable(g)` (graph_util.c) as the single gate both builtins check. Build
Dijkstra as a self-contained static function in `shortestpath.c`, parallel to the existing
`bfs()`, operating over a locally-built vertex-indexed weighted adjacency (reusing
`graph_resolve_edge_weights` + the same `GraphVIdx` index pattern already used elsewhere).
Both `builtin_find_shortest_path`/`builtin_graph_distance` branch once, at the top, on
`graph_weights_usable`.

## Phase 1: `graph_weights_usable` + Dijkstra + dispatch

### Overview
Implement the gate, the algorithm, and wire both builtins to it.

### Changes Required:

#### 1. Weight-usability gate
**File**: `src/graph/graph_util.c`, declared in `graph.h`
**Changes**: `int graph_weights_usable(const Expr* g)` — `graph_is_valid(g)` first;
`arg_count != 3` → `false`; otherwise walk `graph_resolve_edge_weights(g)` and require every
entry to satisfy `expr_is_numeric_like(w)` (`src/expr.c:412` — the codebase's own existing
generic numeric-type check, already covering `EXPR_INTEGER`/`EXPR_BIGINT`/`EXPR_REAL`/
`EXPR_MPFR`/`Rational`; do not hand-roll a narrower type list), be non-`Complex`, and convert
to a `double >= 0` via a small local `graph_weight_to_double(w)` helper (handles each of the
numeric leaf types `expr_is_numeric_like` accepts).

#### 2. Dijkstra, and how the exact output value is produced
**File**: `src/graph/shortestpath.c`
**Changes**: a static `dijkstra()` mirroring `bfs()`'s signature/shape (fills `parent[]`),
using a `double dist[]` **for internal vertex-selection comparisons only** — O(V²) array scan
for the minimum-unvisited-distance vertex each iteration (no heap, matching
`VertexConnectivity`'s existing complexity precedent). This resolves a real gap a
`plan-reviewer` pass caught in the previous draft: a raw `double` accumulator returned
directly as `GraphDistance`'s result would print as `EXPR_REAL` (e.g. `12.`, per
`src/print.c`'s real-vs-integer formatting), not the exact `12` AC-2 expects, and this
codebase treats exact arithmetic as load-bearing throughout. Fix: once `dijkstra()` finds the
parent chain to `t`, reconstruct the **exact** total by evaluating `Plus[w1, ..., wk]` (via
`evaluate()`) over the actual `Expr*` weights (from `graph_resolve_edge_weights`, matched to
the path's edges) — the same exact-arithmetic path every other numeric builtin in this
codebase already goes through, giving `GraphDistance` an exact `EXPR_INTEGER`/`Rational`
result whenever the inputs are exact, and only falling to `EXPR_REAL` if a weight genuinely
was (e.g. `EXPR_MPFR`). `FindShortestPath` needs no such reconstruction — it returns the
vertex path, not a distance value, and the `double`-based selection is only ever used to
choose *which* path, never printed itself.

#### 3. Dispatch
**File**: `src/graph/shortestpath.c`
**Changes**: both builtins check `graph_weights_usable(g)` once, at the top (before/instead
of building the plain `GraphAdj`), and call `dijkstra()` or the existing `bfs()` accordingly.

### Success Criteria:

#### Automated Verification:
- [x] Build succeeds: `make -j$(nproc)` (with `SDKROOT` set per GR-08)
- [x] Portability gate passes: `make check-c99`
- [x] Packed-array audit unaffected: `make check-packed-aware`
- [x] `tests/test_graph.c` passes: updated `test_edge_weights` + new
      `test_weighted_shortest_path` covering AC-1 through AC-7

#### Manual Verification:
- [x] Every Acceptance Criteria row run against the live `./Mathilda -file` REPL, matching
      Expected exactly
- [x] Pre-existing unweighted `FindShortestPath`/`GraphDistance` tests (from the original
      subsystem commit) still pass unmodified

**Implementation Note**: single-phase ticket; proceed straight to docs/changelog after this
phase's verification, per this session's "execute continuously once approved" convention.

---

## Testing Strategy
Extend `tests/test_graph.c` directly, same convention as ticket 1: update the superseded
AC-11 assertions in `test_edge_weights`, add `test_weighted_shortest_path` for the new
behavior.

### Edge Cases & Integration Scenarios:
- Multi-path graphs where hop-count and total-weight disagree (AC-1/AC-2) — the case that
  actually distinguishes Dijkstra from BFS
- Symbolic and negative weights (AC-4/AC-5) — the fallback path, not just the happy path
- Undirected weighted graphs (AC-6) — weight symmetry
- Unreachable target (AC-7) — unchanged semantics

### Manual Testing Steps:
1. Build and start `./Mathilda`
2. Run every Acceptance Criteria row, compare to Expected

## Performance Considerations
O(V²) per call — consistent with `VertexConnectivity`'s existing exact/small-graph
complexity tolerance in this subsystem; not intended for large graphs.

## Migration Notes
None — additive dispatch; unweighted graphs are provably unaffected (AC-3).

## References
- Research: `thoughts/shared/research/2026-08-23-weighted-shortest-path.md`
- Prior ticket (source of this follow-up): `thoughts/shared/plans/2026-08-22-graph-edge-weights.md`
- Direct template: `src/graph/shortestpath.c`'s existing `bfs()`

## Implementation Notes (post-hoc)

All Acceptance Criteria (AC-1 through AC-7) executed directly against the built `./Mathilda`
binary and matched exactly, including the exact-integer check (`Head[GraphDistance[...]]` →
`Integer`, not `Real`) that the plan-reviewer's first BLOCKING finding exists to guard.
`make check-c99`/`make check-packed-aware` both exit 0. `graph_tests` (17 tests, including
the new `test_weighted_shortest_path`) passes standalone. The verification ladder's `unit`
rung reports FAILED for the same pre-existing, unrelated reason documented in
`KIT-FEEDBACK-GRAPH.md` GR-13 (an unrelated flaky optimization test halts the
`for t in *_tests` loop alphabetically before `graph_tests` runs) — not a regression from
this change.
