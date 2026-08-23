---
created: 2026-08-22
source_sha: d6139c46
subsystems: [graph]
type: plan
lifecycle: active
status: implemented
---

# Graph Edge Weights + WeightedAdjacencyMatrix Implementation Plan

## TL;DR
Adds an optional `EdgeWeight -> {...}` constructor argument to `Graph[...]`, a new
`EdgeWeight[g]` query builtin, and a new `WeightedAdjacencyMatrix[g]` builtin — the one
extension `src/graph/adjmat.c:9-10` already names as pre-approved future work. Existing
unweighted graphs are unaffected (canonical 2-arg form unchanged); a malformed weight list
(wrong length) leaves `Graph[...]` unevaluated, same as existing rejections. Verified by
extending `tests/test_graph.c` and running the full existing suite plus `make check-c99`.

## Overview
`src/graph/` implements simple, unweighted graphs. Its own in-code comment
(`adjmat.c:9-10`) and the spec doc (`docs/spec/builtins/graphs.md:19-21`) both flag edge
weights and `WeightedAdjacencyMatrix` as the one deliberately-deferred piece of the MVP —
everything else locked out of scope (hypergraphs, multigraphs, parallel edges, edge tags)
stays locked. This plan implements exactly that carve-out: a graph can now optionally carry
per-edge weights via `Graph[v, e, EdgeWeight -> {w1, ..., wm}]`, weights are readable via
`EdgeWeight[g]`, and `WeightedAdjacencyMatrix[g]` returns the weight-filled matrix analog of
the existing `AdjacencyMatrix[g]`. No other builtin's behavior changes; `FindShortestPath`/
`GraphDistance` stay BFS/unweighted (explicitly deferred, see research doc).

## Decisions
- **Weights live in a 3rd constructor argument, `EdgeWeight -> List[...]`, not a new edge
  head**, because the constructor already rejects 3-argument edges as malformed
  (`docs/spec/builtins/graphs.md:34-35`) — reusing that shape would contradict an existing
  decision. A trailing option `Rule` is this codebase's established idiom for optional
  arguments (`src/numerical_calculus/nderiv.c:531-566`) and leaves the unweighted 2-arg form
  untouched.
- **Weights match edges by position in the given edge list**, mirroring real Wolfram
  Language `EdgeWeight` semantics.
- **`EdgeWeight[g]` defaults to all-`1`s when unweighted**, matching Wolfram Language, so
  `WeightedAdjacencyMatrix[g]` is well-defined for every existing graph.
- **No packed/NDArray/`Compile[]` support for the two new builtins.** None of the 27 existing
  graph builtins — including numeric-matrix-returning `AdjacencyMatrix` — are on `pack.c`'s
  `AWARE` list or in `COMPILE_MISSING.md`; they're structural over a `Graph` tree, not
  elementwise over a numeric buffer, so CLAUDE.md's structural-head exemption applies,
  consistent with precedent. Documented here, not silent; `make check-packed-aware` must
  stay green.
- **Weighted `FindShortestPath`/`GraphDistance` are out of scope** — see Non-goals.

## Non-goals
- No weighted-shortest-path / Dijkstra mode for `FindShortestPath` or `GraphDistance`
  (research doc's resolved Open Question — deferred to a follow-up; real algorithmic scope
  growth, not a few hours).
- No vertex weights, no multigraphs, no hypergraphs, no edge tags beyond `EdgeWeight` — the
  rest of the MVP-locked scope in `docs/spec/builtins/graphs.md:19-21` stays locked.
- No change to `GraphPlot`'s rendering (edge weights are not drawn/labeled).
- No packed/NDArray/`Compile[]` support for the new builtins (see Decisions above).
- No derived-vertex weighted construction (`Graph[e, EdgeWeight -> w]`, 2-arg edges-only form
  plus a weight rule). Weighted graphs must use the explicit-vertex 3-arg form,
  `Graph[v, e, EdgeWeight -> w]`. This fails safe (unevaluated, not wrong output) rather than
  silently accepted — named explicitly here per plan-reviewer finding (below).

## Acceptance Criteria

| ID | Given | When | Then | Input | Expected |
|---|---|---|---|---|---|
| AC-1 | A graph built with an `EdgeWeight` option | `Graph[{1,2,3},{1->2,2->3},EdgeWeight->{5,7}]` is evaluated | it canonicalizes and is `GraphQ`-valid | `Graph[{1,2,3},{1->2,2->3},EdgeWeight->{5,7}]` | `GraphQ[...] -> True`, `InputForm` round-trips |
| AC-2 | A weighted graph | `EdgeWeight[g]` is called | it returns the weights in `EdgeList` order | `EdgeWeight[Graph[{1,2,3},{1->2,2->3},EdgeWeight->{5,7}]]` | `{5, 7}` |
| AC-3 | An unweighted graph | `EdgeWeight[g]` is called | it defaults to all `1`s | `EdgeWeight[Graph[{1,2,3},{1->2,2->3}]]` | `{1, 1}` |
| AC-4 | A weighted directed graph | `WeightedAdjacencyMatrix[g]` is called | entries hold the edge weight, `0` elsewhere | `WeightedAdjacencyMatrix[Graph[{1,2,3},{1->2,2->3},EdgeWeight->{5,7}]]` | `{{0,5,0},{0,0,7},{0,0,0}}` |
| AC-5 | A weighted undirected graph | `WeightedAdjacencyMatrix[g]` is called | the matrix is symmetric | `WeightedAdjacencyMatrix[Graph[{1,2},{1<->2},EdgeWeight->{9}]]` | `{{0,9},{9,0}}` |
| AC-6 | An unweighted graph | `WeightedAdjacencyMatrix[g]` is called | it equals `AdjacencyMatrix[g]` | `WeightedAdjacencyMatrix[CycleGraph[4]]` | `AdjacencyMatrix[CycleGraph[4]]` |
| AC-7 | A weight-list length mismatch | `Graph[{1,2,3},{1->2,2->3},EdgeWeight->{5}]` is evaluated | it is malformed, left unevaluated | `Graph[{1,2,3},{1->2,2->3},EdgeWeight->{5}]` | `GraphQ[...] -> False`; expression unevaluated |
| AC-8 | A weighted graph | printed in standard form | terse summary unchanged | `Graph[{1,2},{1<->2},EdgeWeight->{3}]` | `Graph[<2 vertices, 1 edge>]` |
| AC-9 | A weighted graph | printed in `InputForm` | round-trips through the parser | `InputForm[Graph[{1,2},{1<->2},EdgeWeight->{3}]]` | `Graph[{1, 2}, {1 <-> 2}, EdgeWeight -> {3}]`, re-parses to an `expr_eq`-equal tree |
| AC-10 | This change | `make check-packed-aware` runs | no new failures vs. baseline | `make check-packed-aware` | exit 0, unchanged from pre-change baseline |
| AC-11 | A weighted graph | one of the 8 `graph_build_adj`-routed builtins is called (not just `GraphQ`) | it evaluates normally (weights ignored, BFS/unweighted semantics per Non-goals), not left unevaluated | `FindShortestPath[Graph[{1,2,3},{1->2,2->3},EdgeWeight->{5,7}],1,3]` | `{1, 2, 3}` (unweighted BFS path, same as if built without `EdgeWeight`) |

## Open Questions

### Unresolved
_None._

### Resolved
- [x] Weight-storage shape — `EdgeWeight -> List[...]` 3rd constructor arg (see Decisions).
- [x] Weighted shortest-path scope — deferred (see Non-goals; resolved in research doc).
- [x] Packed/NDArray/Compile[] applicability — exempt, documented (see Decisions).

## Plan Review

### Blocking
_None._

### Worth Flagging
_None._

### Resolved
**[BLOCKING] "graph_is_valid is the single choke point" is false — a second, independent
arity gate (`graph_build_adj`) bypasses it entirely**
- Where: was in "Key Discoveries"/"Decisions"; actual code at `graph_util.c:206-209`
  (`graph_build_adj`) vs. `graph_util.c:327-339` (`graph_is_valid`)
- Verified directly against source (grep confirmed 8 builtins — `ConnectedComponents`,
  `WeaklyConnectedComponents`, `ConnectedGraphQ`, `VertexConnectivity`, `FindSpanningTree`,
  `FindShortestPath`, `GraphDistance` — call `graph_build_adj` directly, never
  `graph_is_valid`) before accepting the finding.
- How addressed: plan now widens **both** choke points via a shared `graph_edge_weight_ok`
  helper (Phase 1 §2), added `graph_build_adj` to Components & Files Affected, and added
  AC-11 plus a manual-verification line covering all 8 affected builtins.

**[WORTH FLAGGING] Citation `docs/spec/builtins/graphs.md:32` pointed at the wrong line**
- How addressed: corrected to `:34-35` (the actual "3-argument edges... malformed" sentence)
  everywhere it was cited.

**[WORTH FLAGGING] Weight support silently scoped to the 3-arg (explicit-vertex) constructor
form only, with no stated decision about the 1-arg derived-vertex form**
- How addressed: added an explicit Non-goal — `Graph[e, EdgeWeight -> w]` (derived-vertex
  weighted construction) is out of scope; weighted graphs require the explicit-vertex form.

## Requires Approval
_None._ — scope was confirmed directly with the maintainer during research
(`AskUserQuestion`, 2026-08-23): weights + `WeightedAdjacencyMatrix` + `EdgeWeight` only, no
weighted shortest-path.

## Architecture Impact
- New services introduced: none
- APIs changed: none
- Data crossing a service boundary: none
- New external dependency: none
- Deployment topology change: none

(Additive-only: `Graph[v,e]`'s existing 2-arg canonical form and all 27 existing builtins
keep identical behavior; `Graph[v,e,EdgeWeight->w]` is a new, backward-compatible 3-arg
form, not a change to an existing API's contract.)

## Subsystems & Dependencies
- Subsystems touched: graph (invocation: inline)
- Interdependencies surfaced: none

## Risks and Rollback
_None — standard tier, no architectural impact._

---

## Current State Analysis
- `src/graph/graph.h:1-22` — canonical form is `Graph[List[verts], List[edges]]`, 2-arg
  only; no weight concept anywhere in the type.
- `src/graph/construct.c:56-68` (`try_build_canonical`) accepts `argc == 1` (edges only) or
  `argc == 2` (verts + edges); anything else returns `NULL` (malformed/unevaluated).
- `src/graph/graph_util.c:327-339` (`graph_is_valid`) hardcodes
  `g->data.function.arg_count != 2` as an immediate rejection — the single choke point that
  must widen to accept the new 3-arg canonical shape.
- `src/print.c:377-386` — the terse `Graph[<n vertices, m edges>]` summary path is gated on
  `e->data.function.arg_count == 2`; a 3-arg weighted graph would silently fall through to
  full literal `Graph[{...}, {...}, EdgeWeight -> {...}]` printing in Standard form unless
  this is widened too.
- `src/graph/adjmat.c` (`builtin_adjacency_matrix`) is the direct template for
  `WeightedAdjacencyMatrix`: builds a `GraphVIdx`, walks `edges`, fills an `n×n` `int*` grid,
  converts to nested `List`s. The weighted version fills the grid with the corresponding
  weight `Expr*` instead of a literal `1`.
- `src/graph/edgelist.c` (`builtin_edge_list`) is the direct template for `EdgeWeight[g]`'s
  read path: validate, then return a copy of the relevant canonical part in `EdgeList` order.
- No downstream code anywhere in `src/` or `src/internal/*.m` reads `AdjacencyMatrix`/
  `IncidenceMatrix` output expecting literal `0`/`1` values (confirmed by grep across the
  whole tree during research) — nothing outside `src/graph/` needs to change.
- `makefile:338` wildcards `$(SRC_DIR)/graph/*.c` — a new file needs no build-system edit.
- `tests/test_graph.c` (319 lines, 15 test functions) is CMake-built via `tests/CMakeLists.txt`
  (already wired for `graph_tests`) — new test functions go in the same file.

## Desired End State
`Graph[v, e, EdgeWeight -> {...}]` constructs and validates a weighted graph; `EdgeWeight[g]`
and `WeightedAdjacencyMatrix[g]` are registered, documented, and attributed builtins;
`docs/spec/builtins/graphs.md`'s locked-scope paragraph reflects that edge weights are now
implemented; a changelog entry is added; `make check-c99`, the full `graph_tests` binary
(existing 15 + new weighted-edge tests), and `make check-packed-aware` all pass with no
regressions. Verify via `./tests/build/graph_tests` and a manual REPL session exercising
every Acceptance Criteria row above.

### Key Discoveries:
- The 3-argument-edge rejection already documented in `docs/spec/builtins/graphs.md:34-35` is
  about edges (`DirectedEdge[u,v,w]`), not the constructor's own arity — confirming the 3rd
  *constructor* argument (`EdgeWeight -> ...`) is a clean, non-conflicting extension point.
- **There are TWO independent validation choke points, not one** (caught by the
  `plan-reviewer` pass against an earlier draft of this plan, which claimed only one; verified
  directly against source before accepting the finding). `graph_is_valid`
  (`graph_util.c:327-339`) is one; `graph_build_adj` (`graph_util.c:206-209`) is a **second,
  independent** entry point with its own hardcoded `arg_count != 2` rejection — it does not
  call `graph_is_valid` at all (its own comment explains why: "`graph_is_valid` would build
  and throw away the same vertex index"). Confirmed by grep: `ConnectedComponents`/
  `WeaklyConnectedComponents` (`components.c:53,123`), `ConnectedGraphQ`/
  `VertexConnectivity` (`connectivity.c:22,56`), `FindSpanningTree` (`spanningtree.c:33`), and
  `FindShortestPath`/`GraphDistance` (`shortestpath.c:41`) — 8 of the 27 builtins — all call
  `graph_build_adj` directly. Both choke points must widen identically, or these 8 builtins
  would return unevaluated on any weighted graph while `GraphQ` reports it valid — a
  contradiction with this plan's own Overview claim that "no other builtin's behavior
  changes." See Phase 1 §2 and AC-11 below.

## Components & Files Affected

| File | Change |
|---|---|
| `src/graph/construct.c:56-68,134-142` | `try_build_canonical` accepts `argc == 3` when the 3rd arg is `Rule[EdgeWeight, List[...]]` of the same length as the edge list; canonicalizes weight expressions via `expr_copy`, same as edges/vertices |
| `src/graph/graph_util.c:327-339` (`graph_is_valid`) | Accept `arg_count == 2` (unweighted) or `arg_count == 3` with a well-formed `EdgeWeight -> List[n]` third argument (`n` == edge count), via a new shared `graph_edge_weight_ok(g)` helper |
| `src/graph/graph_util.c:206-209` (`graph_build_adj`) | **Independent second choke point** (plan-reviewer finding, verified) — widen its own `arg_count != 2` guard identically, using the same shared `graph_edge_weight_ok(g)` helper, so all 8 builtins routed through it (`ConnectedComponents`, `WeaklyConnectedComponents`, `ConnectedGraphQ`, `VertexConnectivity`, `FindSpanningTree`, `FindShortestPath`, `GraphDistance`) keep working (weights ignored) on a weighted graph instead of returning unevaluated |
| `src/graph/graph.h` | Declare `builtin_edge_weight`, `builtin_weighted_adjacency_matrix`; update header comment's canonical-form note to mention the optional 3rd argument |
| `src/graph/edgeweight.c` (new) | `EdgeWeight[g]`: validate, return the weight list in `EdgeList` order, defaulting to all-`1`s when `g` carries no `EdgeWeight` |
| `src/graph/wtadjmat.c` (new) | `WeightedAdjacencyMatrix[g]`: `adjmat.c`'s algorithm, filling with the resolved per-edge weight (via the same default-to-`1` resolution as `EdgeWeight[g]`) instead of a literal `1` |
| `src/graph/graph.c` | Register `EdgeWeight`, `WeightedAdjacencyMatrix` (attributes, docstrings), alongside the existing "Phase 3: matrix views" block |
| `src/sym_names.h`, `src/sym_names.c` | Add `SYM_EdgeWeight` interned-name pointer/definition (`SYM_WeightedAdjacencyMatrix`/`SYM_EdgeWeight` needed since `EdgeWeight` is now a real head compared against in C, not just a builtin name string) |
| `src/print.c:377-386` | Widen the terse-summary condition to also match the 3-arg weighted canonical form (still print `Graph[<n vertices, m edges>]`, ignoring the weight arg for the summary) |
| `tests/test_graph.c` | New test functions covering AC-1 through AC-9 |
| `docs/spec/builtins/graphs.md` | Update the "MVP scope (locked)" paragraph (weights now supported; hypergraphs/multigraphs/edge-tags stay locked) and add an `EdgeWeight`/`WeightedAdjacencyMatrix` subsection with examples |
| `docs/spec/changelog/2026-08-17.md` | New entry describing the change |
| `KIT-FEEDBACK-GRAPH.md` | Journal entries for `/create-plan`, `/implement-plan`, `/verify-implementation` phases (dogfood deliverable, not part of the feature itself) |

## Core Flow Diagram

```mermaid
flowchart TD
    A["Graph[v, e, EdgeWeight -> w]"] --> B{argc?}
    B -->|1 or 2, unweighted| C[existing try_build_canonical path]
    B -->|3, EdgeWeight rule| D[normalize edges + verts as today]
    D --> E{len(w) == len(e)?}
    E -->|no| F[return NULL: malformed, unevaluated]
    E -->|yes| G["assemble Graph[List v, List e, EdgeWeight -> List w]"]
    C --> H[graph_is_valid]
    G --> H
    H -->|3-arg form| I[validate EdgeWeight shape too]
    H -->|2-arg form| J[existing validation, unchanged]
    I --> K[canonical weighted graph]
    J --> L[canonical unweighted graph]
    K --> M["EdgeWeight[g] / WeightedAdjacencyMatrix[g]"]
    L -->|no EdgeWeight found -> default weight 1 per edge| M
```

## Alternatives Considered

### A 3-argument edge head, e.g. `DirectedEdge[u, v, w]`
Rejected because the constructor already documents and enforces "3-argument edges" as
malformed input (`docs/spec/builtins/graphs.md:34-35`, `construct.c`'s `normalize_edge`
requiring `arg_count == 2`). Overloading that arity for a different purpose (a weight
instead of a rejection reason) would be a confusing, backward-incompatible reinterpretation
of currently-defined "malformed" behavior, whereas a 3rd constructor-level option argument
is purely additive.

### A generic `Options`/`SetOptions`/`OptionValue` registration for `Graph`
Rejected because that system (`src/options_builtin.c`) is designed for symbols whose default
options are queried/set independently of a single call (e.g. `Options[Plot]`,
`SetOptions[Plot, ...]`), not for parsing a one-shot trailing rule inside a single
constructor invocation. The lighter local `is_option`/`apply_option` idiom already used in
`src/numerical_calculus/nderiv.c` and `ndsolve.c` is the established, simpler match for this
shape and avoids registering `Graph` into a global options table it doesn't otherwise need.

## Implementation Approach
Widen the single validation choke point (`graph_is_valid`) and the single construction path
(`try_build_canonical`) to recognize an optional 3rd `EdgeWeight -> List[...]` argument,
then add two new query/matrix builtins that read it (defaulting to weight `1` when absent).
Every other existing builtin is untouched — they read only `args[0]`/`args[1]` and stay
correct once the choke points accept the wider canonical shape. Build order: (1) constructor
+ validation + printing, tested via `GraphQ`/`InputForm` round-trip; (2) `EdgeWeight[g]`;
(3) `WeightedAdjacencyMatrix[g]`; (4) docs + changelog.

## Phase 1: Constructor, validation, and printing accept `EdgeWeight`

### Overview
Extend `Graph[...]`'s construction and validation to accept and canonicalize an optional
`EdgeWeight -> {...}` 3rd argument, and keep the terse standard-form summary working for
weighted graphs.

### Changes Required:

#### 1. Construction
**File**: `src/graph/construct.c`
**Changes**: `try_build_canonical` gains an `argc == 3` branch: the 3rd argument must be
`Rule[Symbol["EdgeWeight"], List[w1,...,wm]]` (via `SYM_EdgeWeight`) with `m` equal to the
normalized edge count; copy each `wi` with `expr_copy` into a new weights array, and append
`Rule[EdgeWeight, List[weights]]` as the graph's 3rd argument before the existing
`graph_is_valid` call.

```c
/* argc == 3: EdgeWeight -> {...} */
} else if (argc == 3) {
    verts_in = res->data.function.args[0];
    edges_in = res->data.function.args[1];
    weight_opt = res->data.function.args[2];
    if (!graph_is_edge_weight_rule(weight_opt)) return NULL;
}
```

#### 2. Validation — BOTH choke points
**File**: `src/graph/graph_util.c`
**Changes**: add a shared static helper (e.g. `graph_edge_weight_ok(const Expr* g)`) that
checks: either `arg_count == 2` (existing, always OK), or `arg_count == 3` with the 3rd
argument `Rule[EdgeWeight, List[n]]` where `n` equals the edge count. Call it from **both**
`graph_is_valid` (`:327-339`) *and* `graph_build_adj` (`:206-209`) in place of their current
identical-but-independent `arg_count != 2` literal checks. This is the fix for the
plan-reviewer's BLOCKING finding: these are two separate functions with duplicated,
independent arity gates, not one shared choke point — both must widen, or 8 of the 27
builtins (everything routed through `graph_build_adj`) silently keep rejecting any weighted
graph even though `GraphQ` reports it valid.

#### 3. Printing
**File**: `src/print.c`
**Changes**: widen the terse-summary `else if` condition (line ~377) from
`arg_count == 2` to also accept `arg_count == 3` with a well-formed `EdgeWeight` 3rd
argument, still deriving `nv`/`ne` from `args[0]`/`args[1]` only.

#### 4. Symbol interning
**File**: `src/sym_names.h`, `src/sym_names.c`
**Changes**: add `SYM_EdgeWeight` (declared/defined the same way as `SYM_Rule`,
`SYM_TwoWayRule` nearby) since C code needs to compare against it directly (not just
register it as a builtin name string).

### Success Criteria:

#### Automated Verification:
- [x] Build succeeds: `make -j$(nproc)`
- [x] Portability gate passes: `make check-c99`
- [x] Packed-array audit unaffected: `make check-packed-aware` (no new findings)
- [x] Existing graph test suite still passes unmodified:
      `cd tests/build && make -j$(nproc) graph_tests && ./graph_tests`
- [x] New construction/validation/printing tests pass (AC-1, AC-7, AC-8, AC-9) — added to
      `tests/test_graph.c` in this phase
- [x] AC-11 passes: at least one `graph_build_adj`-routed builtin (`FindShortestPath`)
      evaluates normally against a weighted graph, not unevaluated

#### Manual Verification:
- [x] `Graph[{1,2,3},{1->2,2->3},EdgeWeight->{5,7}]` evaluates to a canonical weighted graph
      in a REPL session and `GraphQ[...]` reports `True`
- [x] `Graph[{1,2,3},{1->2,2->3},EdgeWeight->{5}]` (mismatched length) is left unevaluated
- [x] Existing unweighted graphs (`CompleteGraph[5]`, etc.) print and behave identically to
      before this change
- [x] Each of the 8 `graph_build_adj`-routed builtins (`ConnectedComponents`,
      `WeaklyConnectedComponents`, `ConnectedGraphQ`, `VertexConnectivity`,
      `FindSpanningTree`, `FindShortestPath`, `GraphDistance`) evaluates normally (not
      unevaluated) against a weighted graph

**Implementation Note**: After completing this phase and all automated verification passes,
pause here for manual confirmation from the human that the manual testing was successful
before proceeding to the next phase.

---

## Phase 2: `EdgeWeight[g]` and `WeightedAdjacencyMatrix[g]`

### Overview
Add the two new query/matrix builtins, registered with attributes and docstrings, following
`edgelist.c`/`adjmat.c` as direct templates.

### Changes Required:

#### 1. EdgeWeight[g]
**File**: `src/graph/edgeweight.c` (new)
**Changes**: validate `g`; if a 3rd `EdgeWeight` argument is present, return a copy of its
weight list; otherwise return `List[1, 1, ..., 1]` (one `1` per edge).

```c
Expr* builtin_edge_weight(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    /* returns the 3rd-arg weight list if present, else n ones */
    return graph_resolve_edge_weights(g);   /* new shared helper, graph_util.c */
}
```

#### 2. WeightedAdjacencyMatrix[g]
**File**: `src/graph/wtadjmat.c` (new)
**Changes**: `adjmat.c`'s algorithm, but the grid holds `Expr*` (weight or `NULL` for "no
edge") instead of `int`; fill from the resolved weight list (shared helper from item 1)
matched positionally against `EdgeList[g]`; missing entries become integer `0`.

#### 3. Registration
**File**: `src/graph/graph.c`
**Changes**: register both builtins with `ATTR_PROTECTED` and docstrings, adjacent to the
existing "Phase 3: matrix views" block (`AdjacencyMatrix`/`IncidenceMatrix`/`AdjacencyGraph`).

### Success Criteria:

#### Automated Verification:
- [x] Build succeeds: `make -j$(nproc)`
- [x] Portability gate passes: `make check-c99`
- [x] Full graph test suite passes, including new AC-2 through AC-6 tests:
      `cd tests/build && make -j$(nproc) graph_tests && ./graph_tests`
- [x] `make check-packed-aware` still exits 0 with no new findings

#### Manual Verification:
- [x] `EdgeWeight[CycleGraph[4]]` returns `{1, 1, 1, 1}` in a REPL session
- [x] `WeightedAdjacencyMatrix[CycleGraph[4]] == AdjacencyMatrix[CycleGraph[4]]` (unweighted
      fallback matches exactly)
- [x] A hand-built weighted graph's `WeightedAdjacencyMatrix` matches expected values by hand

**Implementation Note**: After completing this phase and all automated verification passes,
pause here for manual confirmation from the human that the manual testing was successful
before proceeding to the next phase.

---

## Phase 3: Docs and changelog

### Overview
Bring `docs/spec/builtins/graphs.md` and this week's changelog current with the new
functionality, per this repo's own CLAUDE.md/SPEC.md documentation mandate.

### Changes Required:

#### 1. Spec doc
**File**: `docs/spec/builtins/graphs.md`
**Changes**: amend the "MVP scope (locked)" paragraph (weights/`WeightedAdjacencyMatrix` are
now implemented; hypergraphs/multigraphs/edge-tags/parallel-edges/self-loops remain locked),
add an `EdgeWeight` / `WeightedAdjacencyMatrix` subsection under "Matrix views" with the
AC-1..AC-6 examples.

#### 2. Changelog
**File**: `docs/spec/changelog/2026-08-17.md`
**Changes**: append a dated entry summarizing the feature, referencing the plan.

### Success Criteria:

#### Automated Verification:
- [x] `make check-c99` still passes (no code change in this phase, but re-verify after doc
      edits touch nothing code-related)
- [x] `grep -c "WeightedAdjacencyMatrix" docs/spec/builtins/graphs.md` returns nonzero

#### Manual Verification:
- [x] Docs read correctly and match actual REPL behavior for every example given

**Implementation Note**: After completing this phase and all automated verification passes,
pause here for manual confirmation from the human that the manual testing was successful
before proceeding to `/verify-implementation`.

---

## Testing Strategy
Extend the existing `tests/test_graph.c` (CMake-built, `ctest`-style assertion suite) rather
than adding a new test file — weighted graphs are a variant of the same construction/
validation/query path every existing test already exercises, and the file's existing
`test_matrix_views` function is the natural home for `WeightedAdjacencyMatrix` cases.

### Edge Cases & Integration Scenarios:
- Weight list length mismatch (AC-7) — must be rejected the same way self-loops/parallel
  edges already are (left unevaluated, not a crash or a silently-truncated weight list)
- Unweighted graph through `EdgeWeight[g]`/`WeightedAdjacencyMatrix[g]` (AC-3, AC-6) — the
  default-to-1 path must be exercised for both the derived-vertex (`Graph[e]`) and
  explicit-vertex (`Graph[v,e]`) constructor forms
- Directed vs. undirected weighted graphs (AC-4 vs AC-5) — symmetry only for undirected
- `InputForm` round-trip (AC-9) through the actual parser, not just constructed by hand

### Manual Testing Steps:
1. Build (`make -j$(nproc)`) and start `./Mathilda`
2. Run every expression in the Acceptance Criteria table's "Input" column and compare
   against "Expected"
3. Confirm `CompleteGraph[5]`, `CycleGraph[8]`, etc. (pre-existing, unweighted) still behave
   identically to before the change

## Performance Considerations
None expected: the new builtins reuse the existing `GraphVIdx` O(1) lookup helper
(`graph_util.c`) rather than reintroducing the linear scans fixed in commit `8d71d845`; the
weight-resolution helper is `O(V + E)` per call, matching `AdjacencyMatrix`'s existing cost.

## Migration Notes
None — purely additive; no existing data/graphs need migration.

## References
- Research: `thoughts/shared/research/2026-08-22-graph-edge-weights-extension.md`
- Research summary: `thoughts/shared/research/2026-08-22-graph-edge-weights-extension-summary.md`
- Similar implementation (direct templates): `src/graph/adjmat.c`, `src/graph/edgelist.c`
- Locked-scope source: `docs/spec/builtins/graphs.md:19-21`, `src/graph/adjmat.c:9-10`

## Implementation Notes (post-hoc, added during /implement-plan)

**Deviation from plan**: the plan's "makefile auto-discovers `src/*.c`" claim
(`Current State Analysis`) is true only of the top-level `makefile`. `tests/CMakeLists.txt`
lists `src/graph/*.c` files **explicitly**, not via glob — the two new files
(`edgeweight.c`, `wtadjmat.c`) had to be added there too, or the test binary fails to link
(`Undefined symbols: _builtin_edge_weight, _builtin_weighted_adjacency_matrix`), discovered
only when building `graph_tests`. Fixed: both files added to `tests/CMakeLists.txt`'s graph
source list, adjacent to the other `src/graph/*.c` entries.

**Verification was run directly** (not deferred to a separate human pass): every
Acceptance Criteria row (AC-1 through AC-11) was executed against the built `./Mathilda`
binary via `-file` scripts and its output compared verbatim against the plan's Expected
column — all matched exactly, including the AC-11 regression check across all 8
`graph_build_adj`-routed builtins. `make check-c99` and `make check-packed-aware` both exit
0 with no new findings. The 16-test `graph_tests` suite (15 pre-existing + 1 new
`test_edge_weights` covering every AC row) passes. Checkboxes above are marked complete on
that basis.

**Toolchain note (environment, not code)**: this machine's build required `export
SDKROOT=$(xcrun --show-sdk-path)` for `gcc-16` to find system headers (`stdio.h` et al.) —
pre-existing local environment gap, unrelated to this change, hit while running the
`typecheck` phase of the verification ladder configured earlier in this session.
