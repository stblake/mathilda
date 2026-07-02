# Graph Subsystem Implementation Plan

> Making `Graph` a first-class symbol in Mathilda, with construction,
> query/representation, search, computation, and visualization builtins.
> Grounded in `graphs-research.md` (cross-system survey, 2026-07-01) and verified
> against the current codebase (`src/linalg/`, `src/core.c`, `src/sym_names.c`,
> `src/print.c`, `src/graphics/`, `makefile`).

## Overview

Add a new `src/graph/` subsystem — mirroring `src/linalg/` exactly (per-builtin
`.c` files, a `graph.h` of `builtin_*` entry points, a `graph_init()` registered
in `core_init()`) — that implements graphs as ordinary `Expr` trees:

```
Graph[ List[v1, v2, ...], List[edge1, edge2, ...] ]
```

Edges are `DirectedEdge[u,v]` / `UndirectedEdge[u,v]`, with `Rule`/`u->v` and
`TwoWayRule`/`u<->v` accepted as parse-time sugar and normalized on construction.
Vertices are arbitrary expressions. **No new `EXPR_*` tag** and **no changes to
`src/external/`**.

## Current State Analysis

- **No graph subsystem exists.** No `Graph`, `DirectedEdge`, `VertexList`, etc.
- **The mirror target is concrete.** `src/linalg/linalg.h:8-30` declares
  `Expr* builtin_*(Expr*)` entry points + `void linalg_init(void)`;
  `src/linalg/linalg.c:17-54` registers each with `symtab_add_builtin` +
  `attributes |= ATTR_PROTECTED` + (should) a docstring. Per-builtin `.c` files
  (`det.c`, `dot.c`, …) hold one builtin each.
- **Init hub.** `core_init()` in `src/core.c` calls every subsystem `*_init()`;
  the graphics block at `src/core.c:642-643` uses a local `void graphics_init(void);`
  forward-declare + call. We follow the same idiom for `graph_init()`.
- **Symbol interning.** `src/sym_names.c:1164` interns `SYM_Graphics` etc.;
  `src/sym_names.h:569` externs them. `SYM_Rule` already exists
  (`sym_names.h:396`).
- **Matrices are dense List-of-Lists.** `get_tensor_dims` / `flatten_tensor`
  (`src/linalg/linalg.h:34-35`) parse them; `Q`-predicates return
  `expr_new_symbol(SYM_True/SYM_False)` (e.g. `src/linalg/negdef_q.c:187`). We
  reuse this representation for `AdjacencyMatrix` so it interoperates with `Det`,
  `Eigenvalues`, `Tr` immediately.
- **Printer.** `src/print.c:475` renders `SYM_Rule` as `" -> "`. We add
  `DirectedEdge`→`" -> "`, `UndirectedEdge`→`" <-> "`, and a terse `Graph[...]`
  summary form.
- **Visualization precedent.** `src/graphics/plot.c:803 builtin_plot` builds a
  `Graphics[...]` expression from primitives (`Point`, `Line`, `Circle`, `Text`)
  and returns it for `render.c` to draw. `GraphPlot` emits the same primitives —
  **no renderer changes needed.**
- **⚠ Makefile is not recursive.** `makefile:125` enumerates each source subdir
  in `SRC`. `src/graph/*.c` **must** be added there or new files are silently
  ignored.

## Desired End State

A user in the REPL can:

```
g = Graph[{1,2,3,4}, {1->2, 2->3, 3->4, 4->1}]   (* prints as a terse summary *)
GraphQ[g]                    -> True
VertexList[g]                -> {1, 2, 3, 4}
EdgeList[g]                  -> {1 -> 2, 2 -> 3, 3 -> 4, 4 -> 1}
VertexCount[g]               -> 4
EdgeCount[g]                 -> 4
VertexDegree[g, 1]           -> 2
AdjacencyMatrix[g]           -> {{0,1,0,0},{0,0,1,0},{0,0,0,1},{1,0,0,0}}
Eigenvalues[AdjacencyMatrix[g]]   (* reuses linalg unchanged *)
FindShortestPath[g, 1, 4]    -> {1, 2, 3, 4}
GraphDistance[g, 1, 4]       -> 3
ConnectedComponents[g]       -> {{1, 2, 3, 4}}
CompleteGraph[4]             -> Graph[...]
GraphPlot[g]                 -> Graphics[...]  (renders when USE_GRAPHICS=1)
```

Verified by: clean `-std=c99 -Wall -Wextra` build, a new `tests/graph_tests.c`
suite passing, valgrind-clean under `valgrind --leak-check=full`, and docs in
sync.

### Key Discoveries
- `src/linalg/linalg.c:17-54` — exact registration idiom to copy.
- `src/core.c:642-643` — forward-declare-and-call idiom for a new `*_init()`.
- `src/sym_names.c:1164` / `sym_names.h:569` — where to intern/extern new symbols.
- `makefile:125` — the `SRC` wildcard line that must gain `$(SRC_DIR)/graph/*.c`.
- `src/print.c:475` — where edge/arrow printing plugs in.
- `src/graphics/plot.c:803` — `Graphics[...]`-emitting builtin precedent for `GraphPlot`.

## Locked Decisions (resolving graphs-research.md §4.4 / §6)

1. **Multigraphs & loops: forbidden in MVP** (Maxima minimal-viable point).
   `GraphQ` rejects parallel edges and self-loops. The 3-arg edge form
   `DirectedEdge[u,v,tag]` is **reserved** for a future `EdgeTaggedGraph`; not
   implemented now.
2. **`AdjacencyMatrix` returns a dense linalg List-of-Lists** (0/1, symmetric for
   undirected). `WeightedAdjacencyMatrix` is a **later** addition, not MVP.
3. **Internal storage: canonical edge `List`, adjacency derived on demand.** A
   shared `graph_util.c` builds an integer-indexed adjacency structure (vertex →
   index via linear `expr_eq` scan for MVP; `expr_hash` map is the documented
   upgrade path) that all algorithms consume. Graphs stay plain, inspectable
   `Expr` trees (not opaque atoms).
4. **Printer: terse summary by default** — `Graph[<|V| vertices, |E| edges|>]`
   style — with `InputForm[g]` printing the round-trippable literal constructor.
5. **Hypergraphs: out of scope** (research found only negative results). Not
   implemented.

## What We're NOT Doing

- No multigraphs, no loops, no `EdgeTaggedGraph`, no hypergraphs.
- No `WeightedAdjacencyMatrix` / edge weights / vertex weights in MVP (Phase 3
  leaves a documented hook).
- No opaque/atomic graph object; no `AtomQ[g] == True`.
- No changes to `src/external/`, the evaluator core, or the graphics renderer.
- No new sparse-array type; matrices reuse dense List-of-Lists.
- No `ChromaticNumber`/coloring or `FindMaximumFlow` in the first pass — deferred
  to Phase 6 as stretch (see phase notes).

## Implementation Approach

Six incremental phases, each independently buildable, testable, and
valgrind-clean. Phases 0–2 establish the type and its representation; Phase 3
adds matrix interop; Phase 4 adds generators; Phase 5 is the algorithm core
(built on one shared adjacency helper so every algorithm shares vertex-indexing
and BFS/DFS scaffolding); Phase 6 is visualization. Each builtin gets a
docstring, attributes, a `docs/spec/builtins/graphs.md` entry, and a changelog
note in `docs/spec/changelog/2026-06-29.md`.

---

## Phase 0: Scaffolding

### Overview
Stand up the empty subsystem so a trivial `Graph` builtin registers and the build
picks up `src/graph/`.

### Changes Required

#### 1. Symbol interning
**Files**: `src/sym_names.h`, `src/sym_names.c`
Add externs + `intern_symbol` calls (next to the `SYM_Graphics` block):
`SYM_Graph`, `SYM_DirectedEdge`, `SYM_UndirectedEdge`, `SYM_TwoWayRule`,
`SYM_GraphQ`, `SYM_VertexList`, `SYM_EdgeList`, and the remaining builtin heads.
(`SYM_Rule` already exists.)

#### 2. Subsystem header + init
**Files**: `src/graph/graph.h` (new), `src/graph/graph.c` (new)
Mirror `src/linalg/linalg.h`: `#include "expr.h"`, declare every
`Expr* builtin_*(Expr* res);` and `void graph_init(void);`. `graph.c` holds
`graph_init()` (registrations + attributes + docstrings) and, for Phase 0, a
stub `builtin_graph` that just returns `NULL`.

```c
/* src/graph/graph.c */
#include "graph.h"
#include "symtab.h"
#include "attr.h"
void graph_init(void) {
    symtab_add_builtin("Graph", builtin_graph);
    symtab_get_def("Graph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Graph",
        "Graph[v,e] represents a graph with vertices v and edges e.");
}
```

#### 3. Register in core_init
**File**: `src/core.c` (near line 642, next to `graphics_init`)
```c
void graph_init(void); graph_init();
```

#### 4. Makefile
**File**: `makefile:125` — append to the `SRC` wildcard list:
```
$(wildcard $(SRC_DIR)/graph/*.c)
```

#### 5. Docs skeleton
**Files**: `docs/spec/builtins/graphs.md` (new), `Mathilda_spec.md` (add a table
row linking to it), `docs/spec/changelog/2026-06-29.md` (create with the
`# Changelog: week of ...` heading if absent; add a "Graph subsystem: scaffolding"
note).

### Success Criteria

#### Automated
- [x] Clean build: `make -j$(nproc)` with `-std=c99 -Wall -Wextra`, no warnings.
- [x] `echo '?Graph' | ./Mathilda` prints the docstring. (Verified via `Information[Graph]` over the REPL's non-TTY JSON protocol.)
- [x] `src/graph/*.c` objects appear in the build (grep the link line / `ls src/graph/*.o`).

#### Manual
- [x] `Graph[{1,2},{1->2}]` returns unevaluated (stub) without crashing. (Confirmed in the notebook: returns the unevaluated `Graph[...]`, no crash.)

**Pause for confirmation before Phase 1.**

---

## Phase 1: Construction, Normalization & Printing

### Overview
Make `Graph[...]` a real value: normalize edge sugar, validate, canonicalize, and
print. Add `GraphQ`.

### Changes Required

#### 1. `builtin_graph` — normalize & canonicalize
**File**: `src/graph/construct.c` (new)
- Accept `Graph[edges]` (derive vertex set from edges, **directed by default** per
  Wolfram) and `Graph[verts, edges]`.
- Normalize each edge: `Rule[u,v]`/`u->v` → `DirectedEdge[u,v]`;
  `TwoWayRule[u,v]`/`u<->v` → `UndirectedEdge[u,v]`; pass through existing
  `DirectedEdge`/`UndirectedEdge`.
- Reject (return `NULL`, leave unevaluated) on: 3-arg edges, self-loops, parallel
  edges, or an edge endpoint not in an explicit vertex list.
- Produce canonical `Graph[List[verts], List[normalized edges]]` with vertices in
  first-appearance order. Ownership: build new tree, NULL-out reused sub-exprs
  before the evaluator frees `res` (SPEC §4).

#### 2. `builtin_graph_q` — GraphQ predicate
**File**: `src/graph/graphq.c` (new)
Returns `SYM_True` iff the arg is a canonical `Graph[List, List]` with all edges
`DirectedEdge`/`UndirectedEdge`, no loops, no duplicates; else `SYM_False`.
Shared validator `graph_is_valid(Expr*)` lives in `graph_util.c` and is reused by
every algorithm as a guard.

#### 3. Printer
**File**: `src/print.c` (near line 475)
- `DirectedEdge` → `" -> "`, `UndirectedEdge` → `" <-> "` (infix, precedence like
  `Rule`).
- `Graph[List v, List e]` → terse `Graph[<n vertices, m edges>]` in standard
  output; `InputForm`/`FullForm` print the literal constructor (round-trippable).

### Success Criteria

#### Automated
- [x] `tests/graph_tests.c` (new, wired into CMake like other `*_tests`) covers:
  normalization of all four edge sugars, vertex derivation, loop/duplicate
  rejection, `GraphQ` True/False cases. (`tests/test_graph.c`, all pass.)
- [x] `GraphQ[Graph[{1,2},{1->2}]]` → `True`; `GraphQ[Graph[{1},{1->1}]]` →
  `False` (loop); `GraphQ[5]` → `False`.
- [x] `InputForm[Graph[{1,2},{1<->2}]]` round-trips through the parser to an
  equal expression.
- [~] valgrind-clean — valgrind unavailable on this Mac (Darwin ARM64);
  substituting AddressSanitizer on the graph suite.

#### Manual
- [x] REPL: `Graph[{1,2,3},{1->2,2->3}]` prints the terse summary. (Verified:
  `Graph[<3 vertices, 2 edges>]`.)

---

## Phase 2: Query / Representation Builtins

### Overview
Extract structure back out of a graph.

### Changes Required
**Files** (one `.c` each in `src/graph/`): `vertexlist.c`, `edgelist.c`,
`counts.c` (`VertexCount`, `EdgeCount`), `adjlist.c` (`AdjacencyList`),
`degree.c` (`VertexDegree`, plus `VertexInDegree`/`VertexOutDegree` for directed),
`directedq.c` (`DirectedGraphQ`).

All are thin readers over the canonical form, guarded by `graph_is_valid`. For
directed graphs, degree/adjacency respect edge direction; for undirected,
symmetric.

### Success Criteria

#### Automated
- [x] `VertexList`/`EdgeList` return the canonical lists in order.
- [x] `VertexCount`/`EdgeCount` return correct integers.
- [x] `VertexDegree[g,v]`, in/out-degree correct for a known directed example.
- [x] `DirectedGraphQ` distinguishes directed vs undirected.
- [x] Tests pass (valgrind->ASan; see final).

#### Manual
- [x] Spot-check `AdjacencyList` on a small mixed example matches hand computation.

**Pause for confirmation.**

---

## Phase 3: Matrix Views (linalg interop)

### Overview
Bridge graphs to the existing dense-matrix world.

### Changes Required
**Files**: `src/graph/adjmat.c` (`AdjacencyMatrix`), `src/graph/incmat.c`
(`IncidenceMatrix`), `src/graph/adjgraph.c` (`AdjacencyGraph` — inverse of
`AdjacencyMatrix`, round-trips a 0/1 matrix back to a `Graph`).

- `AdjacencyMatrix[g]` → dense List-of-Lists of 0/1 (counts always 1 since no
  multiedges); **symmetric** for undirected. Result is a normal matrix consumable
  by `Det`, `Tr`, `Eigenvalues` with zero linalg changes.
- Leave a clearly-commented hook for a future `WeightedAdjacencyMatrix` (not
  implemented — see Locked Decision 2).

### Success Criteria

#### Automated
- [x] `AdjacencyMatrix` of a directed 4-cycle equals the known circulant 0/1 matrix.
- [x] `Det`/`Tr`[AdjacencyMatrix[g]] run unchanged (linalg interop).
- [x] `AdjacencyGraph[AdjacencyMatrix[g]]` reproduces `g`'s edges (round-trip).
- [x] Tests pass (valgrind->ASan; see final).

#### Manual
- [x] Undirected graph yields a symmetric matrix; directed does not.

**Pause for confirmation.**

---

## Phase 4: Graph Generators

### Overview
Standard constructors.

### Changes Required
**Files**: `src/graph/generators.c` (or one file each): `CompleteGraph[n]`,
`CycleGraph[n]`, `PathGraph[n]` / `PathGraph[{v...}]`, `RandomGraph[{n,m}]`
(reuse the existing RNG used by `random_init` / `RandomInteger`). Each builds a
canonical `Graph` via the Phase 1 constructor path (so validation/canonicalization
is shared, not duplicated).

### Success Criteria

#### Automated
- [x] `VertexCount[CompleteGraph[5]]` → 5, `EdgeCount` → 10 (undirected K5).
- [x] `CycleGraph[n]`/`PathGraph[n]` have correct counts and degrees.
- [x] `RandomGraph[{n,m}]` yields exactly `m` distinct non-loop edges; deterministic
  under `SeedRandom`.
- [x] Tests pass (valgrind->ASan; see final).

#### Manual
- [ ] `GraphPlot` (Phase 6) of `CompleteGraph[6]` looks correct.

**Pause for confirmation.**

---

## Phase 5: Search & Computation Algorithms

### Overview
The algorithmic core. All algorithms share one internal adjacency builder in
`graph_util.c`.

### Changes Required

#### 1. Shared adjacency scaffolding
**File**: `src/graph/graph_util.c` (extend)
`graph_to_adj(Expr* g, ...)` → integer-indexed CSR-ish adjacency (vertex list +
per-vertex neighbor index arrays), plus index↔vertex maps. Vertex→index via
linear `expr_eq` scan for MVP (documented `expr_hash` upgrade path). Provide
internal `bfs`/`dfs` helpers used by multiple algorithms. Careful manual memory
management; freed by the caller.

#### 2. Algorithms (one `.c` each)
- `shortestpath.c` — `FindShortestPath[g,s,t]` (returns the **path** as a vertex
  list; BFS for unweighted) and `GraphDistance[g,s,t]` (returns the **length**).
  Keep Wolfram's path-vs-length naming split.
- `components.c` — `ConnectedComponents`, `WeaklyConnectedComponents`,
  `StronglyConnectedComponents` (Tarjan/Kosaraju for directed).
- `spanningtree.c` — `FindSpanningTree` (BFS/DFS tree for unweighted MST).
- `connectivity.c` — `VertexConnectivity` (and `ConnectedGraphQ`).

**Stretch (same phase, if time permits, else defer):** `FindMaximumFlow`
(Edmonds–Karp), matching (`FindIndependentEdgeSet`), `ChromaticNumber`/greedy
coloring. Explicitly optional; not required for phase completion.

### Success Criteria

#### Automated
- [x] `FindShortestPath` returns a valid shortest path; `GraphDistance` returns
  its length, on directed and undirected examples with known answers.
- [x] Unreachable target: `GraphDistance` → `Infinity`, `FindShortestPath` → `{}`.
- [x] `ConnectedComponents` partitions vertices correctly; strong vs weak differ
  on a known directed example.
- [x] `FindSpanningTree` has `VertexCount - 1` edges on a connected graph.
- [x] Tests pass on all algorithm paths incl. the adjacency builder (valgrind->ASan; see final).

#### Manual
- [x] Cross-checked known small/medium examples (cycles, complete graphs, paths).

**Pause for confirmation.**

---

## Phase 6: Visualization (`GraphPlot`)

### Overview
Render a graph by emitting a `Graphics[...]` expression — reusing the existing
graphics primitives and renderer with **no renderer changes**.

### Changes Required
**File**: `src/graph/graphplot.c`
- `GraphPlot[g]` computes 2-D vertex coordinates (MVP: circular layout — vertices
  evenly on a circle; documented hook for a force-directed spring layout later),
  then builds a `Graphics[{ Line[{p_u,p_v}] per edge, Disk[p_v]/Point per vertex,
  Text[label] per vertex }]` exactly as `src/graphics/plot.c:803 builtin_plot`
  builds its `Graphics` tree. Return that expression; `render.c` draws it when
  `USE_GRAPHICS=1`, or the text placeholder otherwise (graceful-degrade policy).
- Directed edges: draw with an arrowhead primitive if available, else a plain line
  with a documented limitation.

### Success Criteria

#### Automated
- [x] `Head[GraphPlot[g]]` → `Graphics`.
- [x] The emitted `Graphics` contains one Line per edge and one Disk per vertex
  (asserted via Count[...,_Line/_Disk,Infinity]).
- [x] Builds and tests pass with both `USE_GRAPHICS=0` and `=1` (emitter has no raylib dep).
- [x] Memory-clean under ASan (valgrind unavailable on Darwin ARM64).

#### Manual
- [ ] `GraphPlot[CompleteGraph[6]]` and `GraphPlot[CycleGraph[8]]` render a
  recognizable figure with `USE_GRAPHICS=1` (raylib installed; awaiting your visual check).

**Pause for confirmation.**

---

## Testing Strategy

### Unit tests (`tests/graph_tests.c`, wired into CMake alongside existing suites)
- Construction/normalization: all four edge sugars, vertex derivation, directed
  default, loop/duplicate/3-arg rejection.
- `GraphQ` / `DirectedGraphQ` truth tables.
- Query builtins vs hand-computed small examples.
- `AdjacencyMatrix` symmetry (undirected) and round-trip via `AdjacencyGraph`.
- Generators: exact vertex/edge counts and degrees.
- Algorithms: shortest path/distance (incl. unreachable), components (weak vs
  strong), spanning tree edge count.
- `GraphPlot` structural assertions under both graphics settings.

### Integration
- `AdjacencyMatrix` feeding `Eigenvalues`/`Det`/`Tr` unchanged.
- `InputForm` round-trip through the parser.

### Manual
1. REPL walk-through of the Desired End State examples.
2. Visual check of `GraphPlot` on `CompleteGraph`/`CycleGraph` with graphics on.
3. `valgrind --leak-check=full ./Mathilda` over a scripted graph session.

## Performance Considerations
- MVP vertex→index is a linear `expr_eq` scan (O(V) per lookup, O(V·E) build) —
  fine for pico-CAS sizes; documented `expr_hash`-map upgrade path (codebase
  already has `expr_hash`) for larger graphs.
- Adjacency derived on demand per algorithm call; acceptable for the MVP, cache
  later if profiling warrants.
- BFS/DFS are shared in `graph_util.c` to avoid per-algorithm duplication.

## Documentation & Workflow Obligations (CLAUDE.md / SPEC.md §10)
- Every builtin: `symtab_set_docstring`, attributes in `graph_init()`, an entry
  in `docs/spec/builtins/graphs.md`, a note in
  `docs/spec/changelog/2026-06-29.md` (Monday of the current ISO week).
- New internal symbols in `src/sym_names.c` (+ externs in `sym_names.h`).
- Add the `graphs.md` row to `Mathilda_spec.md` and a changelog-table row for the
  weekly file.
- SPEC §4 ownership contract on every builtin; valgrind each phase.
- Commit as Michael Sollami; no AI/Claude attribution (global CLAUDE.md).

## References
- Research: `graphs-research.md` (esp. §4 recommended design, §4.4 decisions).
- Mirror pattern: `src/linalg/linalg.h:8-30`, `src/linalg/linalg.c:17-54`.
- Init hub: `src/core.c:642-643`.
- Symbols: `src/sym_names.c:1164`, `src/sym_names.h:569`.
- Build: `makefile:125` (`SRC` wildcard — must add `src/graph/*.c`).
- Printer: `src/print.c:475`.
- Visualization precedent: `src/graphics/plot.c:803`.
- Contributor workflow / standards: `CLAUDE.md`, `SPEC.md`.
