# Graphs, Hypergraphs & Networks — Cross-System Research

*Research report to inform adding a graph-theory subsystem to Mathilda (a pico
Mathematica clone in C99). Compiled 2026-07-01. All external claims below
survived 3-vote adversarial verification (25/25 confirmed, 0 refuted).*

---

## 1. Executive Summary

Across every surveyed system a graph is fundamentally **an adjacency structure
over a vertex set plus an edge collection**, but the systems split into two
design camps:

1. **Atomic / opaque graph objects** with rich attribute properties — Wolfram's
   `Graph` head (`AtomQ` returns `True`), SageMath's backend-selectable
   `DiGraph`.
2. **Explicitly exposed data structures** — NetworkX's dict-of-dicts-of-dicts,
   Maxima's Lisp adjacency lists, igraph's edge-list-plus-metadata.

Directed / undirected / mixed / weighted variants are universally handled either
through **constructor flags** (SageMath `loops`/`multiedges`/`weighted`; Wolfram
`DirectedEdges->False`) or **distinct classes** (NetworkX
`Graph`/`DiGraph`/`MultiGraph`/`MultiDiGraph`). Undirected edges are consistently
modeled as **two opposing directed edges**, yielding symmetric adjacency
matrices.

**Multigraph (parallel-edge) support is the key design differentiator.** Wolfram
needs `EdgeTaggedGraph` with unique per-edge tags because a plain `Graph` cannot
distinguish parallel edges; NetworkX adds a fourth dict level keyed by edge key;
Maxima supports neither multi-edges nor loops.

The expected algorithm surface is **broad and largely uniform** across systems
(shortest paths incl. named Dijkstra / Bellman-Ford / all-pairs, connectivity /
strong components, max flow, matching, coloring, MST) — signaling that a CAS
graph subsystem is expected to expose **distinct named algorithms**, not a single
generic dispatcher.

**For Mathilda:** the idiomatic representation is a `Graph[vertexList, edgeList]`
head where edges are `DirectedEdge[u,v]` / `UndirectedEdge[u,v]` (or
`Rule` / `TwoWayRule`), implemented as a new `src/graph/` subsystem mirroring
`src/linalg/` (per-builtin `.c` files, a `graph.h` of `builtin_*` entry points, a
`graph_init()` registered in `core_init()`), reusing the existing `Expr` tagged
union and List-of-Rules convention, with **no changes** to `src/external/`.

---

## 2. Per-System Findings

### 2.1 Wolfram Language (Mathilda's primary model)

**Representation — atomic, canonical, attribute-carrying.** A Wolfram `Graph` is
an *atomic raw object* (`AtomQ` → `True`), always normalized to a canonical
standard form `Graph[vertices, edges, ...]`. It carries vertex properties
(`VertexLabels`, `VertexCoordinates`, `VertexWeight`, `VertexStyle`) and edge
properties (`EdgeLabels`, `EdgeStyle`, `EdgeWeight`). Constructed via `Graph`,
`CompleteGraph`, `RandomGraph`, `AdjacencyGraph`, `GraphData`.
`[confidence: high · 3-0]`

**Edge syntax.** `Graph[{e1,...}]` builds from edges; `Graph[{v1,...},{e1,...}]`
adds explicit vertices.
- Directed: `u -> v` / `DirectedEdge[u,v]` / `Rule[u,v]`
- Undirected: `u <-> v` / `UndirectedEdge[u,v]` / `TwoWayRule[u,v]`

By default a **directed** graph is generated from a list of rules;
`DirectedEdges->False` reinterprets them as undirected. **Mixed graphs** (a
collection of both directed and undirected edges) are supported. `[3-0]`

**Multiple interchangeable views of one graph.** Matrix forms
(`AdjacencyMatrix`, `IncidenceMatrix`, `KirchhoffMatrix`,
`WeightedAdjacencyMatrix`) and list forms (`VertexList`, `EdgeList`,
`AdjacencyList`, `IncidenceList`, `EdgeRules`) coexist, and `AdjacencyGraph[amat]`
round-trips a matrix back to a graph. `[3-0]`

**Adjacency matrix semantics.** Entry `a_ij` is the *number* of directed edges
from vertex `i` to vertex `j`. An undirected edge = two opposing directed edges
⇒ **symmetric** matrix; directed graphs may be unsymmetric; multigraph
multiplicities appear as integer counts `> 1`. Wolfram returns a `SparseArray`
convertible via `Normal`. `AdjacencyMatrix` ignores weights —
`WeightedAdjacencyMatrix` handles those. `[3-0]`

**Weights & parallel edges.** `VertexWeight` / `EdgeWeight` accept any numeric or
symbolic expression. **Parallel edges are indistinguishable in a plain `Graph`**;
`EdgeTaggedGraph` assigns unique (auto-generated integer or user-supplied) tags,
making parallel edges individually identifiable and annotatable. This is the
central multigraph-identity decision. `[3-0]`

**Algorithms (core, current API).** `FindShortestPath` (returns the *path*),
`GraphDistance` (returns the *length*), `ConnectedComponents`,
`FindMaximumFlow`, `BreadthFirstScan`, `GraphUnion`. Note the useful naming split
between path-returning and length-returning functions. `[3-0]`

**Legacy `Combinatorica`.** As of Wolfram v10, most `Combinatorica` functionality
is built into the core system — the package is **obsolete**. Historically it
offered both a generic `ShortestPath` dispatcher (with an `Algorithm` option) and
distinct named implementations (`Dijkstra`, `BellmanFord`,
`AllPairsShortestPath`), plus `NetworkFlow`, `ResidualFlowGraph`,
`BipartiteMatching`, `MaximalMatching`, `StableMarriage`. This shows a generic
entry point and named algorithms can coexist. `[3-0]`

### 2.2 NetworkX (canonical exposed adjacency structure)

**Dict-of-dicts-of-dicts.** Graphs are an explicit nested dictionary: outer keyed
by node, inner by neighbor, innermost = the edge-attribute dict. `G[u][v]`
returns the edge attribute dictionary itself. Chosen over lists (for fast sparse
lookup) and sets (so data can be attached to an edge); edges are found/removed in
two dict look-ups. `[3-0]`

**Four classes via `Di` + `Multi` prefixes:** `Graph`, `DiGraph`, `MultiGraph`,
`MultiDiGraph`. `DiGraph` keeps **separate successor (`G.succ`) and predecessor
(`G.pred`)** structures (giving O(1) in-degree). `MultiGraph`/`MultiDiGraph` add a
**fourth dict level keyed by an edge key**. This is the class-per-variant
alternative to Wolfram's single-head-with-flags design. `[3-0]`

### 2.3 Maxima `graphs` package (minimal-viable design point)

Graphs are represented internally by **adjacency lists implemented as Lisp
structures**; edges/arcs are lists of length 2; vertices are integer ids. It
supports **simple undirected graphs and digraphs only** — explicitly **no
multiple edges, no loops, no hypergraphs**. A digraph may hold both `u->v` and
`v->u` without those counting as duplicates.

Algorithm suite: `shortest_path`, `shortest_weighted_path`,
`connected_components`, `strong_components`, `vertex_connectivity`,
`min_vertex_cut`, `max_flow`, `max_matching`, `chromatic_number`,
`minimum_spanning_tree`. `[3-0]`

This (**integer-id vertices + length-2 edge lists + simple-graph-only**) is the
closest analog to what a pico-CAS would ship first.

### 2.4 SageMath (mature, backend-selectable)

`DiGraph` defaults to **simple** (no loops, no multiedges) and **unweighted**,
tuned by separate `loops` / `multiedges` / `weighted` flags. Construction
**auto-detects** many input formats: integer vertex count, edge list,
dict-of-lists (out-neighbors), dict-of-dicts with labels, square adjacency
matrix, nonsquare incidence matrix, `[V, f]` rule form with a boolean adjacency
function, `dig6` string, and conversion from NetworkX/igraph. `[3-0]`

**Selectable backends** at construction: `'dense'`, `'sparse'` (default), and
`'static_sparse'`. `static_sparse` is faster, smaller in memory, and
**immutable** — so such graphs can be used as dictionary keys (`immutable=True`
is shorthand for it). Relevant to Mathilda later, since an immutable hashable
graph aligns with Mathilda's immutable-`Expr` convention. `[3-0]`

### 2.5 igraph (C library — closest architectural analog, partial coverage)

*(Extracted but not in the top-25 verified set; treat as strong-but-unverified.)*
igraph models a graph as **an edge list plus a metadata table** — a multiset of
integer-labeled vertex-ID pairs (ordered if directed, unordered if undirected)
rather than an adjacency matrix. In C it uses typed containers
(`igraph_vector_int_t` for vertex/edge ids), hand-rolled because C lacks
generics. It can export multiple views: adjacency matrix, adjacency list, edge
list, incidence list. This is the most directly transferable C99 design.

### 2.6 Hypergraphs (a gap — mostly negative results)

- **Maxima:** no hypergraph support.
- **SageMath:** `IncidenceStructure(ground_set, blocks)` models a hypergraph as
  an incidence structure / set system (verified as a source but not top-ranked).
- **Wolfram:** `HypergraphPlot` and the Wolfram Physics hypergraph model exist
  but were **not captured in any verified claim**.
- **HyperNetX** (PNNL) represents a hypergraph as a set system with an incidence
  store (source fetched, not top-ranked).

The natural Mathilda representation for a hypergraph edge is a **`List` of
vertices** (arbitrary arity) rather than a 2-element `DirectedEdge`/
`UndirectedEdge`, e.g. `Hypergraph[{v...}, {{a,b,c}, {c,d}, ...}]`.

---

## 3. Comparison Table

| System | Core structure | Directed/undirected | Weighted | Multigraph | Hypergraph | Vertex ids |
|--------|----------------|---------------------|----------|------------|------------|------------|
| **Wolfram** | Atomic `Graph[v,e,...]`, sparse-array views | flags (`DirectedEdges`), mixed OK | `EdgeWeight`/`VertexWeight` (numeric or symbolic) | via `EdgeTaggedGraph` tags | `HypergraphPlot` / Physics model (not verified) | any expression |
| **NetworkX** | dict-of-dicts-of-dicts | 4 classes (`Di`/`Multi`) | edge-attr dict | 4th dict level (edge key) | separate lib (HyperNetX) | any hashable |
| **Maxima** | Lisp adjacency lists | separate graph/digraph | weighted path fns | **none** | **none** | integers |
| **SageMath** | selectable dense/sparse/static_sparse backend | flags + classes | `weighted` flag | `multiedges` flag | `IncidenceStructure` | any |
| **igraph (C)** | edge list + metadata table | flag | attributes | supported | no | integers |

---

## 4. Recommended Mathilda Design

### 4.1 Expression representation

Reuse the existing `Expr` tagged union and the List-of-Rules convention. **No new
`EXPR_*` tag is needed** — a graph is a normal `EXPR_FUNCTION`:

```
Graph[ List[v1, v2, ...], List[edge1, edge2, ...] ]
```

with edge heads:
- `DirectedEdge[u, v]`   (accept `Rule[u,v]` / `u -> v` as sugar → normalize)
- `UndirectedEdge[u, v]` (accept `TwoWayRule[u,v]` / `u <-> v` as sugar)
- optional tag slot `DirectedEdge[u, v, tag]` for multigraphs (EdgeTaggedGraph
  style) — decide up front (see §4.4).

Vertices are arbitrary expressions (symbols, integers), matching Wolfram and
keeping uniformity with the rest of the system.

**Canonical form + `GraphQ`.** Following Wolfram, normalize any constructor input
to `Graph[List[verts], List[edges]]` inside `builtin_graph`, and gate structural
pattern-matching behind a `GraphQ` predicate rather than exposing internals.
Mathilda graphs stay plainly-inspectable `Expr` trees (unlike Wolfram's opaque
atom), which is simpler for a pico-CAS; hashability/immutability already come for
free from the immutable-`Expr` convention.

### 4.2 Subsystem structure (mirror `src/linalg/`)

Concrete pattern observed in the codebase:

- **`src/graph/graph.h`** — `#include "expr.h"`, declare each `Expr*
  builtin_*(Expr* res);` entry point and `void graph_init(void);` (mirrors
  `src/linalg/linalg.h:1-33`).
- **One `.c` per builtin** (e.g. `construct.c`, `adjmat.c`, `shortestpath.c`,
  `components.c`), plus a `graph.c` holding `graph_init()` and shared helpers —
  exactly the `src/linalg/` layout.
- **`graph_init()`** registers each builtin and sets attributes, following
  `src/linalg/linalg.c:17-54`:
  ```c
  symtab_add_builtin("Graph", builtin_graph);
  symtab_get_def("Graph")->attributes |= ATTR_PROTECTED;
  symtab_set_docstring("Graph", "Graph[v,e] represents a graph ...");
  ```
- **Register in `core_init()`** — add `graph_init();` alongside the other
  `*_init()` calls in `src/core.c` (see the block at `src/core.c:532-551`).
- **Intern internal symbols in `sym_names.c`** — `DirectedEdge`,
  `UndirectedEdge`, `Graph`, `EdgeTaggedGraph`, matching the
  `SYM_Graphics = intern_symbol("Graphics")` pattern at `src/sym_names.c:1164`.
- **⚠ Makefile is NOT recursive.** `makefile:125` lists every source subdir
  explicitly in the `SRC` wildcard. A new subsystem **must add**
  `$(wildcard $(SRC_DIR)/graph/*.c)` there — new files under `src/graph/` will
  otherwise be silently ignored.

### 4.3 Builtin surface (MVP → later)

**Construction:** `Graph`, `CompleteGraph`, `CycleGraph`, `PathGraph`,
`AdjacencyGraph` (inverse of `AdjacencyMatrix`), `RandomGraph`.

**Query/representation:** `GraphQ`, `VertexList`, `EdgeList`, `VertexCount`,
`EdgeCount`, `AdjacencyList`, `AdjacencyMatrix`, `IncidenceMatrix`,
`WeightedAdjacencyMatrix`, `DirectedGraphQ`, `VertexDegree`.

**Algorithms (named, per cross-system convention):** `FindShortestPath` (path)
and `GraphDistance` (length) — keep Wolfram's split; `ConnectedComponents`,
`StronglyConnectedComponents` / `WeaklyConnectedComponents`, `FindMaximumFlow`,
`FindSpanningTree` (MST), `VertexConnectivity`, graph coloring / `ChromaticNumber`,
matching. Expose distinct named functions; a generic dispatcher with an
`Algorithm` option is optional (Combinatorica precedent).

### 4.4 Key up-front decisions

1. **Multigraphs:** bake an optional tag slot into the edge
   (`DirectedEdge[u,v,tag]`, EdgeTaggedGraph-style) **now**, or forbid multi-edges
   (Maxima-style) and add tags later? Recommendation: **forbid multi-edges/loops
   in the MVP** (Maxima's minimal-viable point), reserve the 3-arg edge form for a
   later `EdgeTaggedGraph`.
2. **AdjacencyMatrix return type:** reuse the existing linalg **dense
   List-of-Lists** matrix (interoperates immediately with `Det`, `Eigenvalues`,
   `Tr`) for the MVP; add a sparse representation later as Wolfram does. Split
   unweighted `AdjacencyMatrix` (0/1 or counts) from `WeightedAdjacencyMatrix`.
3. **Internal storage:** simplest MVP stores the canonical edge `List` and derives
   adjacency on demand; a hash-map (the codebase already has `expr_hash`) from
   vertex → neighbor is the NetworkX-style upgrade path. igraph's edge-list +
   integer-id-vector model is the closest C99 precedent if performance matters.
4. **Printer:** teach `print.c` to render `DirectedEdge[u,v]` as `u -> v`,
   `UndirectedEdge[u,v]` as `u <-> v`, and `Graph[...]` either as a terse
   summary (`Graph[<n vertices, m edges>]`, Wolfram-atom style) or as the literal
   constructor for round-trippable `InputForm`.

### 4.5 Documentation & workflow obligations (per `CLAUDE.md` / `SPEC.md`)

Each new builtin needs: a `symtab_set_docstring`, appropriate attributes in
`graph_init()`, an entry in the relevant `docs/spec/builtins/` category file, and
a note in the current week's `docs/spec/changelog/<Monday>.md`. Add a
`docs/spec/builtins/graphs.md` category and reference it from `Mathilda_spec.md`.
Run `valgrind --leak-check=full` — the builtin ownership contract (§4 of SPEC)
applies: return a new `Expr*` on success (evaluator frees `res`), `NULL` when the
graph can't be evaluated, NULL-out reused sub-expressions before the wrapper is
freed.

---

## 5. Caveats & Coverage Gaps

- **Uneven system coverage.** Verified claims are dominated by Wolfram (the most
  directly relevant model), with solid primary-source coverage of NetworkX,
  Maxima, and SageMath. **Maple, GAP, and Magma produced no surviving verified
  claims.** **igraph** was characterized (edge-list + typed C vectors) but its
  claims did not make the top-25 verified set — treat §2.5 as strong-but-unverified.
- **Hypergraphs are essentially unaddressed.** The only verified hypergraph
  finding is negative (Maxima has none). Wolfram's `HypergraphPlot` / Physics
  hypergraph model was requested but not captured in a verified claim.
- **Vendor self-reporting.** All Wolfram claims cite Wolfram's own reference docs
  — authoritative for API facts, self-reported for design rationale.
- **`Combinatorica` is explicitly obsolete** (superseded in Wolfram v10); its
  algorithm inventory reflects legacy naming, not the current API.
- **Version sensitivity.** NetworkX (v3.6.1) and SageMath docs are current but
  version-specific.
- **Mathilda integration guidance is a recommendation, not an observed
  implementation** — inferred from inspecting `src/linalg/`, `src/core.c`,
  `sym_names.c`, and `makefile`. No graph subsystem exists yet. (The codebase
  facts cited in §4.2 — makefile non-recursion, the `*_init()` pattern, the
  `symtab_add_builtin` + attribute idiom — were verified directly against the
  repo.)

---

## 6. Open Questions

1. **Maple / GAP / Magma / igraph internals:** data structures, construction/query
   APIs, and algorithm suites. igraph is the closest architectural analog to a C99
   CAS and deserves a dedicated follow-up crawl of `igraph.org/c/doc/`.
2. **Wolfram hypergraphs:** how `HypergraphPlot`, the Wolfram Physics hypergraph
   model, and any `Hypergraph` object actually represent hyperedges — and whether
   any surveyed non-Wolfram system ships a first-class hypergraph structure.
3. **Multigraph policy for Mathilda:** tag slot from day one vs. deferred
   simple-graph-only (see §4.4.1).
4. **AdjacencyMatrix representation:** reuse linalg dense List-of-Lists vs. a
   SparseArray-style form; and how to split `AdjacencyMatrix` vs.
   `WeightedAdjacencyMatrix` (see §4.4.2).

---

## 7. Sources

**Primary (vendor / authoritative):**
- Wolfram — Graphs & Networks guide: <https://reference.wolfram.com/language/guide/GraphsAndNetworks.html>
- Wolfram — Graph Construction & Representation: <https://reference.wolfram.com/language/guide/GraphConstructionAndRepresentation.html>
- Wolfram — `Graph`: <https://reference.wolfram.com/language/ref/Graph.html>
- Wolfram — `AdjacencyMatrix`: <https://reference.wolfram.com/language/ref/AdjacencyMatrix.html>
- Wolfram — `EdgeTaggedGraph`: <https://reference.wolfram.com/language/ref/EdgeTaggedGraph.html>
- Wolfram — `IncidenceGraph`: <https://reference.wolfram.com/language/ref/IncidenceGraph.html>
- Wolfram (legacy) — Combinatorica Graph Algorithms: <https://reference.wolfram.com/language/Combinatorica/guide/GraphAlgorithms.html>
- NetworkX — Introduction / data structures: <https://networkx.org/documentation/stable/reference/introduction.html>
- Maxima — `graphs` package manual: <https://maxima.sourceforge.io/docs/manual/Package-graphs.html>
- SageMath — `DiGraph`: <https://doc.sagemath.org/html/en/reference/graphs/sage/graphs/digraph.html>
- SageMath — `IncidenceStructure` (hypergraphs): <https://doc.sagemath.org/html/en/reference/graphs/sage/combinat/designs/incidence_structures.html>
- igraph C — Basic interface: <https://igraph.org/c/doc/igraph-Basic.html>
- igraph C — Data structures: <https://igraph.org/c/doc/igraph-Data-structures.html>
- igraph (Python) — `Graph` API: <https://python.igraph.org/en/main/api/igraph.Graph.html>
- Graphviz — Cgraph library: <https://www.graphviz.org/pdf/cgraph.pdf>
- IGraph/M documentation: <https://szhorvat.net/mathematica/IGDocumentation/>
- HyperNetX — hypergraph class: <https://github.com/pnnl/HyperNetX/blob/master/hypernetx/classes/hypergraph.py>

**Secondary / reference:**
- SageMath Wiki — Graph Theory Software Survey: <https://wiki.sagemath.org/graph_survey>
- Wikipedia — Adjacency list: <https://en.wikipedia.org/wiki/Adjacency_list>
- GeeksforGeeks — adjacency list vs. matrix: <https://www.geeksforgeeks.org/dsa/comparison-between-adjacency-list-and-adjacency-matrix-representation-of-graph/>

*Research stats: 5 search angles · 22 sources fetched · 108 claims extracted ·
25 verified (25 confirmed, 0 refuted) · 104 agent calls.*
