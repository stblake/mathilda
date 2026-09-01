---
ticket: RG-2
created: 2026-08-31T09:07:35-0400
researcher: Michael Sollami
source_sha: 854419997e4b0d8b35fb401783a17b13db2495c9
branch: find-vertex-coloring
repository: mathilda
topic: "How are graph algorithms structured in this repo — location, shared substrate, representation, tests, and the conformance surface for a new algorithm?"
tags: [research, codebase, graph, architecture, conformance, packed-arrays]
subsystems: [graph]
type: research
lifecycle: active
status: complete
last_updated: 2026-08-31
last_updated_by: Michael Sollami
---

# Research: graph algorithm structure in `src/graph/`

**Date**: 2026-08-31T09:07:35-0400
**Researcher**: Michael Sollami
**Git Commit**: `854419997e4b0d8b35fb401783a17b13db2495c9`
**Branch**: `find-vertex-coloring`
**Repository**: mathilda

## TL;DR

A graph is a plain `Expr` tree — `Graph[List[verts], List[edges]]`, no C struct, no
`EXPR_*` tag — and every algorithm rebuilds an int-indexed `GraphAdj` from it per call.
27 heads exist, covering representation, traversal and connectivity; all 25 classic
combinatorial algorithms grepped for are absent, so a search/optimization head would be
the first of its kind. Unresolved: no edge-weight representation exists, and the ladder
trace covers only the consumer side — the producer side is unexamined and a measured
in-tree signal contradicts a blanket exemption.

## Summary

The subsystem is unusually consistent, which makes the conformance surface easy to state
and easy to satisfy. One builtin per translation unit; every prototype under a `Phase`
banner in `graph.h`; all registration centralized in `graph.c` as a rigid three-line
idiom; a genuine shared substrate in `graph_util.c` (vertex-index hash, adjacency
scaffold, component counter). Error handling is uniform to the point of being a rule:
zero `Message()` calls in the entire directory, `return NULL` for everything, and no
head takes any option at all.

Two structural facts dominate what a new algorithm would face. First, the subsystem
stops at traversal — there is no precedent for an option-bearing head, a head returning
a certificate or assignment, or one with an exponential worst case. `VertexConnectivity`
is the closest thing, and it is brute-force subset enumeration. Second, the packed/
NDArray/`Compile[]` ladder that `CLAUDE.md` mandates does not reach this subsystem on
the **consumer** side — no head accepts a packed argument — and not by documented
decision: 23 of 27 heads are absent from every audit's candidate universe. But the
**producer** side was not examined, and `src/ndarray.c:547` records a measured PageRank
speedup of 14.1 s → 486 ms attributed to graph kernels reading an index buffer directly.
A blanket "graphs are exempt" conclusion does not survive that.

## Open Questions

### Unresolved

_None._

### Resolved

- [x] Does a new graph algorithm returning machine-numeric vectors owe a packed producer
  path? — **Decide by return type.** A head returning a `Graph` object, or a vertex list
  of arbitrary expressions, carries no obligation; only a head returning a dense numeric
  vector or matrix does. This resolves the blanket-exemption problem without either
  asserting an exemption the evidence contradicts or committing every future head to
  packed work. Consequence for a plan: the return type must be settled *before* the
  ladder obligation can be stated, so it belongs in the plan's Decisions rather than
  discovered during implementation. _(picked "Not applicable — depends on return type"
  over three alternatives, including deferring to an empirical `--survival` run)_
- [x] `AdjacencyMatrix` — real gap or misclassification? — **Real gap, record it.** It
  returns a dense 0/1 matrix built element-by-element (`src/graph/adjmat.c:54`), so
  `NUMERIC` (`tools/numeric_coverage.py:252`) is the correct classification and the
  missing packed path is a genuine, previously unrecorded gap. Not to be fixed by
  removing it from the tool's `NUMERIC` set. _(picked "Real gap — record it" over
  correcting the tool or leaving it)_
- [x] Edge weights — how should a weighted algorithm handle their absence? — **A separate
  prior change.** The representation must be designed as its own ticket, touching
  `AdjacencyMatrix`, `graph_build_adj` and the canonical `Graph[...]` form. Any weighted
  algorithm is blocked until that lands; an unweighted algorithm is unaffected. _(picked
  "Separate prior change" over unweighted-only or designing it alongside)_
- [x] The eight backwards-reading "frees res" header comments — fix or leave? — **Leave
  them**, as a separate change. § 3 records which eight are wrong and points a new file at
  `graph.h:21`, which is enough to stop the template propagating the error. _(picked
  "Leave them — separate change" over a comment-only drive-by fix)_

- [x] Is there prior or abandoned graph-algorithm work not obvious from git? — RG-1
  (`RandomGraph[{n,m},k]`) is to be treated as the conformance template a new
  algorithm follows, not merely as history. _(chose "Treat RG-1 as the template to
  conform to" from the offered options)_
- [x] How deep should the map go into the mandatory numeric-surface ladder? — Trace it
  fully for graph heads. Done; see Finding 6. _(chose "Trace it fully for graph heads")_
- [x] Anything deliberately out of scope inside `src/graph/`? — `graphplot.c` and the
  rendering path. Its *implementation* was excluded from all analysis; `GraphPlot` is
  still counted in the head census (27 heads, the ladder arithmetic, test coverage),
  because excluding a registered head from a population count would misstate the
  subsystem. _(chose "graphplot.c / rendering")_
- [x] May the binary be built and probed to verify claims empirically? — Only if already
  current. It was **not** current (`src/facint.c`, `src/info.c`, `src/geometry.c` and
  others are newer than `./Mathilda`), so no probing was done and every finding below
  is source-derived. _(chose "Build only if it's already current")_

## Research Review

Reviewed 2026-08-31 by `plan-reviewer`, coverage-gap + contradiction lenses. Two BLOCKING
and five WORTH FLAGGING findings. Ten citation sites were independently sampled; nine said
what this document claimed. All findings were verified at source by the researcher before
being accepted — one was overstated by the reviewer and is recorded as corrected.

### Blocking

_None._

### Worth Flagging

**[WORTH FLAGGING] The `NDEBUG` finding said "silently disarms", which overstates what is lost**
- Where: § 7, the `NDEBUG` finding
- Why: `assert_eval_eq` prints `FAIL/Expected/Actual` to stderr *before* the `assert()`
  (`tests/test_utils.h:24-28`), so the failure is visible in the log; only the exit code
  is lost, which fools CTest. "Blind suite" and "green CTest on a failing suite" imply
  different fixes.
- Addressed: § 7 now states the distinction explicitly and names the fix (wrap in
  `ASSERT`), rather than calling the suite blind.

**[WORTH FLAGGING] "Probed" was used for grep results in a document that stakes its credibility on having probed nothing**
- Where: § TL;DR and § 8
- Why: the Resolved question and § 6 both state emphatically that the binary was stale and
  nothing was probed. Using "probed" for a source grep in the most-quoted lines invites
  the exact misreading the honesty note exists to prevent.
- Addressed: both now read "grepped for", with the source-only qualifier inline.

**[WORTH FLAGGING] `PageRank` is not absent from the tree, only from the symbol table**
- Where: § 8, "absent from the whole tree"
- Why: `PageRank` appears as prose at `src/ndarray.c:547` and in
  `docs/spec/changelog/2026-07-27.md`. Stated as a grep result, the claim was wrong for
  one of the 25.
- Addressed: § 8 carries an explicit qualifier distinguishing "no registered head" from
  "name appears nowhere".

**[WORTH FLAGGING] `graphplot.c` was declared out of all analysis yet `GraphPlot` is counted in every head-population claim**
- Where: § Open Questions Resolved vs. § 6's arithmetic, § 7's coverage claim, § 8's list
- Why: the scope note and the counts disagreed, and the ladder arithmetic that
  § Requires Approval rests on shifts by one depending on the reading.
- Addressed: the Resolved entry now scopes the exclusion to the rendering *implementation*
  and states that the head census deliberately includes `GraphPlot`.

**[WORTH FLAGGING] The header-comment boilerplate is inconsistent tree-wide and should not be copied verbatim — but not "wrong in all 16 files" as reported**
- Where: § 3, ownership discussion
- Why: RG-1's touch surface includes mirroring the header comment, so a false ownership
  claim would propagate into a new file. **Reviewer overstated the scope**: verified by
  grep, 8 files carry the backwards-reading bare "frees res" and 8 state it correctly as
  "the evaluator frees res", including `graph.h:21`. The code itself is uniform and
  correct in all 16.
- Addressed: § 3 now names both groups file-by-file and points a new file at the
  `graph.h:21` wording.

### Resolved

**[BLOCKING] Finding 6 claimed the numeric ladder was "traced fully" but never asked the producer-packing question, and the tree contains a measured counter-signal**
- Where: § 6, and § Architecture Insights ("the exemption is real and probably correct")
- Why: the trace covered only the *consumer* side — whether a graph head accepts a packed
  argument. `SPEC.md` §9 names a separate question the sweep asks
  (`tools/nd_fastpath_sweep.py:273`, `--survival`): which *producers* return a plain
  `List` that would have packed. `AdjacencyMatrix` builds its matrix with
  `expr_new_integer` element by element (`src/graph/adjmat.c:54`);
  `ConnectedComponents`/`VertexDegree`/`GraphDistance` share that shape. Worse,
  `src/ndarray.c:540-547` records the packed gather as "the operation every sparse-matrix
  and graph kernel is built out of" with a measured PageRank of 14.1 s against 486 ms —
  a live in-tree signal that directly undercuts the blanket-exemption conclusion a plan
  would have inherited.
- Resolved 2026-08-31: verified at source. § 6 now scopes the trace to the consumer side
  explicitly, names the unexamined `--survival` producer question and the four heads with
  that shape, quotes the `ndarray.c:540-547` counter-signal, and replaces "exemption is
  real and probably correct" with a revised conclusion that the blanket claim is **not**
  supported. § Architecture Insights and § Requires Approval updated to match.

**[BLOCKING] "24 heads outside every audit's candidate universe" contradicted this document's own `AdjacencyMatrix` exception two paragraphs later**
- Where: § 6 and § Summary
- Why: `AdjacencyMatrix` is one of the 24, and being an explicit member of
  `numeric_coverage.py`'s `NUMERIC` set (`:252`) is precisely being inside an audit's
  candidate universe. The count is load-bearing — it is the quantitative basis for the
  § Requires Approval "exempt de facto" ask.
- Resolved 2026-08-31: recount verified and corrected to 3 backlogged + 1 named-but-
  unsatisfied + 23 outside every audit = 27, in both § 6 and § Summary.

## Requires Approval

Two scope calls a human should make before a plan is built on this. **Edge weights**: the
absence is subsystem-wide, so a weighted algorithm silently expands into a representation
change affecting `AdjacencyMatrix`, `graph_build_adj` and the `Graph[...]` canonical form.
**The ladder**: this research establishes that graph heads are undocumented and unaudited
on the consumer side, but it does **not** establish that they are exempt — the producer
question was never asked, and `src/ndarray.c:547` is measured evidence pointing the other
way. So the approval needed is not "add an `EXEMPT` line" but the prior question: does a
new algorithm returning machine-numeric vectors owe a packed producer path? Deciding that
by writing an exemption would settle it in the direction the evidence does not support.

---

## Research Question

> Map how graph algorithms are structured in this repo — where they live, what they
> share, how graphs are represented, how they're tested, and what a new algorithm would
> have to conform to. I have something specific in mind but I'm not telling you yet;
> don't optimize for a guess.

## Detailed Findings

### 1. Representation — a plain `Expr` tree, no struct

No C struct and no new `EXPR_*` tag. Stated at `src/graph/graph.h:6-14`:

```
 *     Graph[ List[v1, v2, ...], List[edge1, edge2, ...] ]
 *
 * where each edge is DirectedEdge[u, v] or UndirectedEdge[u, v].
```

- Canonical form produced by `builtin_graph` (`src/graph/construct.c:134-142`), validated
  by `graph_is_valid` (`src/graph/graph_util.c:327-340`).
- **Vertices are arbitrary expressions**, not integers `1..n` (`src/graph/graph.h:14`).
  Only `AdjacencyGraph[m]` constructs integer vertices, and that is local to it
  (`src/graph/adjgraph.c:51-52`).
- **Edge sugar is normalized once, at construction**: `Rule`/`->` →`DirectedEdge`,
  `TwoWayRule`/`<->` → `UndirectedEdge`, in `normalize_edge()`
  (`src/graph/construct.c:39-52`). Algorithms only ever see the two canonical heads.
- Directed/undirected discrimination is `graph_edge_kind()`
  (`src/graph/graph_util.c:36-42`), which returns an **interned pointer** compared by
  identity, not `strcmp`.
- Simple graphs only: self-loops and parallel edges are rejected by `graph_check`
  (`src/graph/graph_util.c:299-325`), parallel-edge detection via an `EKSet` hash of
  encoded endpoint pairs (`src/graph/graph_util.c:155-182`).

**Edge weights do not exist.** The only trace is an aspirational comment,
`src/graph/adjmat.c:9-10`: "Future hook: a WeightedAdjacencyMatrix would fill entries
with edge weights instead of 1 (Locked Decision 2); not implemented in the MVP."
`src/graph/shortestpath.c:3` documents its BFS as "Unweighted" for the same reason.

### 2. The shared substrate — `graph_util.c`

This is the reuse answer, and it is more substantial than the one-shared-helper
impression the algorithm files give. Two data structures, both built per call, both
caller-freed.

`GraphVIdx` (`src/graph/graph_util.c:55-60`) — open-addressed hash, vertex `Expr*` →
integer index, keyed on `expr_hash`/`expr_eq`. **Keys are borrowed**, so the graph must
outlive the index (`src/graph/graph.h:60-67`). API: `graph_vidx_new`/`_get`/`_put`/`_free`
(`graph_util.c:80-130`), plus internal `vidx_build` (`:134-141`) and `vidx_grow` (`:100-116`).

`GraphAdj` (`src/graph/graph.h:112-117`) — the int-indexed adjacency scaffold:

```c
typedef struct GraphAdj {
    int   n;
    const Expr* verts;         /* borrowed */
    int*  outdeg; int** out;   /* successors   */
    int*  indeg;  int** in;    /* predecessors */
} GraphAdj;
```

Built by `graph_build_adj` (`graph_util.c:206-260`) in two passes — count degrees
(`:230-238`), then allocate and fill (`:239-256`); an `UndirectedEdge` contributes
symmetrically to both `out` and `in`. Freed by `graph_adj_free` (`:198-204`).

The one shared *traversal*: `graph_count_components` (`graph_util.c:262-292`), iterative
DFS over `out ∪ in`, with an optional `removed` mask — which is what makes
`VertexConnectivity`'s subset enumeration possible.

Full `graph_util.c` inventory: `head_is_sym` (`:25`), `graph_is_list` (`:32`),
`graph_edge_kind` (`:36`), `graph_vidx_free` (`:62`), `vidx_alloc` (`:70`),
`graph_vidx_new` (`:80`), `vidx_slot` (`:88`), `graph_vidx_get` (`:94`), `vidx_grow`
(`:100`), `graph_vidx_put` (`:118`), `vidx_build` (`:134`), `ekset_init` (`:158`),
`ekset_insert` (`:167`), `graph_vertex_index` (`:184`, linear scan for single lookups),
`graph_check` (`:299`), `graph_adj_free` (`:198`), `graph_build_adj` (`:206`),
`graph_count_components` (`:262`), `graph_is_valid` (`:327`).

### 3. The canonical algorithm skeleton

Five steps, in this order, in every algorithm head:

1. Validate arg count/shape → `return NULL` on mismatch.
2. Convert the `Graph` Expr to a `GraphAdj` via `graph_build_adj` → `NULL` if invalid.
3. Run the algorithm over freshly `calloc`'d scratch arrays.
4. Build a **fresh** result Expr, `expr_copy`-ing any vertex sub-expression that must
   survive.
5. Free every scratch array, then `graph_adj_free`, then return.

Canonical example — `builtin_find_shortest_path`, `src/graph/shortestpath.c:49-71`. Note
step 4: vertices are **copied out** of the borrowed `a->verts`, never transplanted, so
`res` is left wholly intact for the evaluator to free. Zero graph builtin calls
`expr_free(res)` and none uses the NULL-out-before-free pattern from `SPEC.md` §4 —
copying makes it unnecessary.

**Do not copy the file header comment verbatim.** RG-1's touch surface includes "mirror
the header comment", and half of these headers state the ownership rule in a way that
reads backwards. Eight files end the line with a bare "…; frees res" —
`shortestpath.c:11`, `components.c:12`, `adjmat.c:12`, `adjgraph.c:12`, `incmat.c:7`,
`spanningtree.c:8`, `connectivity.c:12`, `generators.c:15` — which reads as *the builtin*
freeing `res`, the opposite of both the code above and `SPEC.md` §4. The other eight say
it correctly ("the evaluator frees res"): `adjlist.c:11`, `degree.c:15`, `edgelist.c:3`,
`directedq.c:3`, `counts.c:3`, `vertexlist.c:3`, `graphplot.c:12`, and `graph.h:21`, which
is explicit and right. The code is uniform; only the boilerplate is inconsistent. A new
file should follow the `graph.h:21` wording.

Per-file algorithm inventory:

| File | Head(s) | Algorithm | Core loop |
|---|---|---|---|
| `shortestpath.c` | `FindShortestPath`, `GraphDistance` | unweighted BFS over `out[]` | `:26-32` |
| `components.c` | `ConnectedComponents`, `WeaklyConnectedComponents` | iterative stack DFS over `out ∪ in` | `:60-69` |
| `components.c` | `StronglyConnectedComponents` | **Tarjan**, iterative (explicit callstack + child cursor) | `:95-116` |
| `spanningtree.c` | `FindSpanningTree` | BFS spanning **forest**, preserves original edge orientation via `original_edge_copy` (`:18-28`) | `:44-60` |
| `connectivity.c` | `VertexConnectivity` | brute-force enumeration of all size-`k` vertex subsets, each tested by re-running `graph_count_components` | `:37-49` |
| `connectivity.c` | `ConnectedGraphQ` | delegates to `graph_count_components` | `:20-28` |
| `degree.c` | `VertexDegree`/`In`/`Out` | O(E) edge scan; all-vertices form uses `GraphVIdx`, **not** `GraphAdj` | `:77-90` |
| `counts.c`, `graphq.c`, `directedq.c` | counts + predicates | O(1) reads / linear scan, no adjacency built | — |

**Duplication worth knowing**: `weak_components` (`components.c:51-73`) reimplements
the same out∪in DFS that `graph_count_components` (`graph_util.c:262-292`) already
provides, rather than calling it. Not a bug; it means "each file reimplements its
traversal" is the working norm and a reviewer should not read a new local BFS as
sloppiness.

### 4. Conventions — uniform to the point of being rules

- **`return NULL` is the only error channel.** `Message(`/`symtab_message` across all of
  `src/graph/`: **0 matches**. A new head emitting a diagnostic would be the first.
- **No options anywhere.** `Method`/`OptionValue`/`OptionsPattern` across `src/graph/`:
  **0 matches**. Every head is option-free.
- **`ATTR_PROTECTED` only.** 27 `ATTR_` hits in the directory, all `ATTR_PROTECTED`
  (`graph.c:21,30,36,…,174`). No `ATTR_LISTABLE`, `ATTR_NUMERICFUNCTION`, `ATTR_HOLDALL`.
- **Degenerate sizes are guarded at the allocation**, via the
  `calloc((size_t)(n > 0 ? n : 1), …)` idiom (`components.c:56,58`; `connectivity.c:34`;
  `graph_util.c:225-228,245-246,264-265`), so a 0-vertex graph never hits a 0-byte
  `calloc`.
- **Disconnection is an answer, not an error**: `ConnectedGraphQ` → `False`,
  `FindShortestPath` → `{}`, `GraphDistance` → `Infinity`, `FindSpanningTree` → a forest.

### 5. Registration — four files, seven sites

| Site | What | Where |
|---|---|---|
| 1 | `symtab_add_builtin("H", builtin_h)` | `src/graph/graph.c` |
| 2 | `symtab_get_def("H")->attributes \|= ATTR_PROTECTED` — **this idiom, never `symtab_set_attributes`** | `src/graph/graph.c` |
| 3 | `symtab_set_docstring("H", "...")` | `src/graph/graph.c` |
| 4 | prototype under the right `Phase` banner | `src/graph/graph.h` |
| 5 | `extern const char* SYM_H;` | `src/sym_names.h:917-946` |
| 6 | `const char* SYM_H = NULL;` | `src/sym_names.c:859-890` |
| 7 | `SYM_H = intern_symbol("H");` | `src/sym_names.c:1745-1775` |

Sites 5-7 are the ones RG-1's plan does not mention, because `RandomGraph` already had
its symbol. A new head does not, and omitting site 7 leaves a `NULL` pointer that
identity-compares equal to nothing — a silent failure, since `graph_edge_kind` and
friends compare interned pointers.

`graph.h` Phase banners: Phase 2 query/representation (`:83`), Phase 3 matrix views
(`:94`), Phase 4 generators (`:99`), Phase 5 shared scaffolding (`:105`), Phase 5
search & computation (`:127`), Phase 6 visualization (`:137`). A new algorithm belongs
under Phase 5 search & computation.

`graph_init()` is already wired at `src/core.c:944-945`, between `graphics_init()` and
`fourier_init()`. No new init call is needed.

**No graph functionality lives in `.m`**: grepping `src/internal/` for graph names
returns nothing. The subsystem is C-only.

`DirectedEdge`/`UndirectedEdge` are **not registered builtins** — they have interned
`SYM_*` pointers but no `symtab_add_builtin` or docstring. They are structural tags.

### 6. The numeric-surface ladder — consumer side clear and undocumented, producer side unexamined

Traced fully, per the scope answer. Every result here is source-derived; the binary was
stale so nothing was probed.

- **`src/pack.c`** — `AWARE` list `:489-787`, `NOT_AWARE` `:818`, `INT64_OK` `:878`.
  Grep for `Graph|Vertex|Edge`: **no matches in the file at all.**
- **`src/ndkernels.c`, `src/ndinteger.c`, `src/ndreduce.c`, `src/ndstruct.c`** — no
  matches in any of the four. No graph head has an ND kernel.
- **`src/compile/`** (incl. `compile_infer.c`, every `compile_emit_*.c`) — no matches.
  No `Compile[]` lowering, scalar or array.

The audit tools, and the distinction that matters:

- `tools/nd_fastpath_sweep.py` `OFF_BUFFER` (`:608-691`) contains **three** graph heads:
  `PathGraph` (`:690`, comment `:689` "Graph constructor over a coordinate/edge vector"),
  and `GraphQ` + `DirectedGraphQ` (`:666`, under the "Type / identity predicates"
  block). `OFF_BUFFER` is documented at `:50-52` as a **measured backlog ratchet** — a
  known-unfixed gap, explicitly *not* an exemption.
- `tools/check_packed_aware.py` `EXEMPT` (`:47-197`) — **no graph head**. Compare the
  `GEO-1` geometry entries (`:79-104`, `:111-116`), which carry named, dated written
  justifications. That is what a documented exemption looks like, and no graph head has one.
- `tools/check_array_exactness.py` `EXEMPT` (`:154`+), `tools/nd_surface_audit.py`,
  `tools/compile_coverage.py` — no graph head in any of them.

**Why they are silent, mechanically** — this is the part that makes the absence a
finding rather than a non-event. `check_packed_aware.py`'s `heads_with_fast_paths()`
(`:261-304`) walks every `.c` under `src/` including `src/graph/` (`all_sources()`,
`:205-211`, has no graph exclusion) and attributes a head only when it finds a
`DISPATCH_MARKERS` token (`is_ndarray(`, `EXPR_NDARRAY`, …, `:246-255`). None of those
tokens appears anywhere in `src/graph/`, so the tool never *nominates* a graph head.
`compile_coverage.py` (`:19-23`) draws its probe pool from the kernel registry and
`AWARE` — both empty here — so graph heads are never probed.
`tools/numeric_coverage.py`'s `category()` (`:266-275`) classifies **25 of the 27**
`symbolic`, because the `io` keyword test matches `"Graphics"`, not `"Graph"` (`:271-273`).
The two exceptions: `AdjacencyMatrix` returns `numeric` at the first branch (`:267-268`),
and `GraphPlot` returns `io` because `"Plot"` is itself in the keyword tuple (`:272`).

So on the consumer side: 3 heads noticed and parked in a backlog (`PathGraph`, `GraphQ`,
`DirectedGraphQ`), 1 head named by an audit but satisfied by none
(`AdjacencyMatrix`, below), 23 heads outside every audit's candidate universe, and 0
heads with a written justification anywhere.

**`AdjacencyMatrix` is named but unsatisfied.** `tools/numeric_coverage.py:252` puts it
in the `NUMERIC` set — the only graph head so classified — yet it has no `AWARE` entry,
no kernel and no lowering. It builds its matrix element-by-element with `expr_new_integer`
(`src/graph/adjmat.c:54`), so it is a producer that hands back a plain `List` which would
have packed. That is the one place the subsystem's de-facto exemption and the tooling's
own classification already disagree.

**The trace above is the CONSUMER side only, and that is a real limit on this finding.**
`SPEC.md` §9 describes a second, separate question the sweep asks — its `--survival` half
(`tools/nd_fastpath_sweep.py:273`): which *producers* hand back a plain `List` that would
have packed, since packing is a chain and a producer that drops it makes its **consumers**
slow rather than itself. That question was never asked here. `AdjacencyMatrix`,
`ConnectedComponents`, `VertexDegree` and `GraphDistance` all have the producer shape —
they return lists of machine integers built one `expr_new_integer` at a time.

**And there is a measured in-tree counter-signal.** `src/ndarray.c:540-547`, on the
packed-position arm of `Part`:

> "`x[[idx]]` — the operation every sparse-matrix and **graph** kernel is built out of —
> degraded the whole Part and materialised BOTH arrays: a 1.6e6-index gather out of a
> 100000-element vector cost 389 ms against 15.7 ms once the index buffer is read
> directly, and a **PageRank built on it 14.1 s against 486 ms**."

So the substrate itself already records graph work as a packed-buffer *consumer*, with a
29× measured difference. Note `PageRank` is not a registered head (§ 8) — this is the
substrate anticipating one. **Revised conclusion**: the consumer side is genuinely clear
and undocumented; the producer side is unexamined; and a blanket "a `Graph` object is the
non-machine-object case, therefore graphs are exempt" claim is *not* supported — a graph
algorithm returning machine-numeric vectors plausibly has real ladder obligations, and
`ndarray.c:547` is the evidence a plan should start from rather than the exemption.

### 7. Tests

- **One file**: `tests/test_graph.c`, 320 lines, 15 test functions, `main()` at `:297-319`.
- **Helper**: `assert_eval_eq(input, expected, is_fullform)` — shared, defined
  `tests/test_utils.h:18-30`, not graph-specific.
- **CMake**: `tests/CMakeLists.txt:972-975` — target `graph_tests`, links
  `$<TARGET_OBJECTS:mathilda_common>` (the whole interpreter, so tests run the real
  evaluator, not mocks). Run: `make graph_tests && ./graph_tests`, or
  `ctest -R graph_tests`.
- **Coverage**: all 27 registered heads are referenced by name. No zero-coverage head.
- **Strong convention**: directed/undirected duplication — most query and algorithm
  tests appear twice, once per edge kind (`test_query_builtins` `:105` vs
  `test_query_undirected` `:139`).
- **Thin or absent**: empty graph (`Graph[{},{}]`) is **tested nowhere**; single-vertex
  appears only in the printing test (`:49`); `IncidenceMatrix` has one assertion and no
  directed `-1/+1` orientation case despite its docstring promising it
  (`graph.c:94-95`); `VertexInDegree`/`VertexOutDegree` are tested only in the
  single-vertex-argument form; `VertexConnectivity` has no articulation-point case.

**Finding — under `NDEBUG` this suite exits 0 on failure.** `assert_eval_eq` uses libc
`assert()`, a no-op under `-DNDEBUG` (CMake `Release`). `tests/test_utils.h:54-57` warns
about exactly this: the `ASSERT`/`ASSERT_MSG` macros (`:74-99`) `exit(1)` unconditionally,
`assert_eval_eq` does not. Only `test_inputform_roundtrip` (`:78-102`) and
`test_random_graph` (`:213-228`) use the surviving macro.

Precisely what is and isn't lost, because it changes the fix: `assert_eval_eq` writes
`FAIL: <input> / Expected / Actual` to **stderr before** the `assert()`
(`tests/test_utils.h:24-28`). So under `NDEBUG` the mismatch is still printed and a human
reading the log sees every failure — but the process keeps running and exits 0, so
**CTest reports green on a failing suite**. The suite is not blind; its exit code is. The
fix is wrapping the call in `ASSERT`, not rewriting the assertions.

**No memory test for graph heads.** No valgrind/memcheck CTest configuration references
`graph_tests`; the only leak script in the tree is
`tests/scripts/geometry_leakcheck.sh`, for a different subsystem.

**Other surfaces**: `benchmarks/29-graph-ops/graph_ops.m` is a *performance* benchmark
(20000-vertex sparse graph) whose `check[...]` calls are sanity values, not assertions
(`:29,32-33,36,42-43,47-48,52-53,57-58`); results in `benchmarks/REPORT.md:180-185`. No
`.m` correctness script, no `docs/experiments/` folder touching graphs.

### 8. What is absent — the gap inventory

27 heads registered: `Graph`, `GraphQ`, `VertexList`, `EdgeList`, `VertexCount`,
`EdgeCount`, `AdjacencyList`, `VertexDegree`, `VertexInDegree`, `VertexOutDegree`,
`DirectedGraphQ`, `AdjacencyMatrix`, `IncidenceMatrix`, `AdjacencyGraph`,
`CompleteGraph`, `CycleGraph`, `PathGraph`, `RandomGraph`, `FindShortestPath`,
`GraphDistance`, `ConnectedComponents`, `WeaklyConnectedComponents`,
`StronglyConnectedComponents`, `FindSpanningTree`, `ConnectedGraphQ`,
`VertexConnectivity`, `GraphPlot`.

Grepped for (source only — nothing was probed at runtime) and **absent as registered
heads**, 25 of 25: `FindVertexColoring`,
`VertexColoring`, `ChromaticNumber`, `FindEdgeColoring`, `FindClique`, `CliqueNumber`,
`FindIndependentVertexSet`, `FindMaximumFlow`, `FindMinimumCut`,
`FindHamiltonianCycle`, `FindEulerianCycle`, `TopologicalSort`, `FindCycle`,
`FindShortestTour`, `PageRank`, `BetweennessCentrality`, `ClosenessCentrality`,
`FindVertexCut`, `FindKClique`, `IsomorphicGraphQ`, `FindGraphIsomorphism`,
`MaximalBipartiteMatching`, `FindPostmanTour`, `TreeGraphQ`, `PlanarGraphQ`.

One qualifier on that list: `PageRank` **does** appear in the tree as text — the substrate
comment at `src/ndarray.c:547` and `docs/spec/changelog/2026-07-27.md` — but as prose, not
as a registered head. "Absent" above means absent from the symbol table, which is what
matters for a new implementation; it does not mean the name appears nowhere.

So the subsystem today is **representation + traversal + connectivity**. A
search/optimization head would be the first in the subsystem to need any of: an option
(`Method ->`), a returned assignment/certificate rather than a vertex or edge list, an
exponential worst case with a heuristic fallback, or a weight. The only precedent for
"expensive combinatorial search" is `VertexConnectivity`'s brute-force subset
enumeration (`connectivity.c:32-52`), which is honest but not a scalable model.

### 9. Docs obligations

- `docs/spec/builtins/graphs.md` — per-head `##` section: prose description, then a
  fenced block of `Input (* Output *)` lines. Canonical entry: `GraphQ` at `:48-58`. The
  file mixes granularity: larger groups (Query/representation `:60-84`, Generators
  `:106-123`) put several heads under one heading with bullets and a shared example
  block. Subsystem overview at `:1-21`.
- `docs/spec/changelog/2026-08-31.md` — today is **Monday 2026-08-31**, so this is the
  current ISO-week file. It **already exists** (heading: "Changelog: week of 2026-08-31
  (Mon) – 2026-09-06 (Sun)"), so a change appends a `##` section rather than creating it.
- `Mathilda_spec.md:58` links the graphs page in its navigational table; per `CLAUDE.md`
  it needs updating only if a new top-level category appears — a new algorithm does not
  qualify.

### 10. RG-1 as the conformance template

Per the scope decision, RG-1 is the pattern to follow, not just history. Its artifact set
lives on branch `rg-1-random-graph-count` under `thoughts/shared/tickets/RG-1/`:
`research.md`, `research-summary.md`, `plan.md`, `plan-summary.md`, `adversarial.md`,
`validation.md`. What it establishes procedurally:

- **A seven-file touch surface** (its `Components & Files Affected` table): the
  implementation `.c`, its header comment, the `graph.h` signature comment, `graph.c`
  registration/docstring, `tests/test_graph.c`, `docs/spec/builtins/graphs.md`, and the
  weekly changelog.
- **An Acceptance Criteria table with executable rows** — 19 of them, each
  `Given/When/Then` plus a literal input expression and expected output, so each maps to
  an `assert_eval_eq` row.
- **An `## Entry Points` section** naming the signature, which the plan contract then
  cross-checks: every argument named there must appear in at least one AC row.
- **A validation document that reports honestly** — RG-1's marks AC-17 `NOT VERIFIED`
  and the full suite `NOT RUN` rather than rounding up (`validation.md`, Automated
  Verification Results).
- **`## Plan Review` findings transcribed and resolved in place**, with two BLOCKING
  findings moved to `### Resolved` with a note on how each was addressed.

## Code References

- `src/graph/graph.h:6-14` — canonical `Graph[List, List]` form, "vertices are arbitrary expressions"
- `src/graph/graph.h:112-117` — the `GraphAdj` struct
- `src/graph/graph.h:60-67` — `GraphVIdx` borrowed-key warning
- `src/graph/graph.h:83,94,99,105,127,137` — the six Phase banners
- `src/graph/graph_util.c:36-42` — `graph_edge_kind`, interned-pointer comparison
- `src/graph/graph_util.c:206-260` — `graph_build_adj`, the two-pass adjacency build
- `src/graph/graph_util.c:262-292` — `graph_count_components`, the one shared traversal
- `src/graph/graph_util.c:299-325` — `graph_check`: self-loop, parallel-edge, endpoint validation
- `src/graph/graph_util.c:327-340` — `graph_is_valid`
- `src/graph/shortestpath.c:49-71` — **the canonical skeleton**
- `src/graph/shortestpath.c:21-34` — file-local BFS
- `src/graph/components.c:86-119` — iterative Tarjan
- `src/graph/components.c:51-73` — `weak_components`, the duplicated traversal
- `src/graph/connectivity.c:32-52` — brute-force subset enumeration, the only combinatorial-search precedent
- `src/graph/construct.c:39-52` — `normalize_edge`, edge-sugar canonicalization
- `src/graph/adjmat.c:9-10` — the only mention of edge weights anywhere
- `src/graph/graph.c:20-26` — the three-line registration idiom
- `src/core.c:944-945` — `graph_init()` call site
- `src/sym_names.h:917-946` / `src/sym_names.c:859-890` / `src/sym_names.c:1745-1775` — the three `SYM_*` sites
- `src/pack.c:489-787` — the `AWARE` list; no graph head appears
- `tools/nd_fastpath_sweep.py:666,690` — `GraphQ`, `DirectedGraphQ`, `PathGraph` in `OFF_BUFFER`
- `tools/nd_fastpath_sweep.py:50-52` — `OFF_BUFFER` documented as a backlog, not an exemption
- `tools/check_packed_aware.py:47-197` — `EXEMPT`; no graph head. `:79-104` is the `GEO-1` precedent
- `tools/numeric_coverage.py:252` — `AdjacencyMatrix` classified `NUMERIC`
- `tools/numeric_coverage.py:266-275` — `category()`, why graph heads fall through to `symbolic`
- `tests/test_graph.c:297-319` — the 15-test `main()`
- `tests/test_utils.h:18-30` — `assert_eval_eq`
- `tests/test_utils.h:54-57` — the `NDEBUG` warning
- `tests/CMakeLists.txt:972-975` — the `graph_tests` target
- `benchmarks/29-graph-ops/graph_ops.m` — the performance benchmark

## Architecture Insights

- **The subsystem mirrors `src/linalg/`**: one builtin per TU, prototypes grouped by
  Phase banner, registration in a hub file. A new algorithm is a **new `.c` file** plus
  edits to four existing files — not a change to any existing algorithm.
- **Uniformity is the design.** Zero `Message()`, zero options, one attribute, one error
  channel. This makes the conformance surface small, and makes any deviation (a first
  option, a first diagnostic) a visible architectural decision rather than a detail.
- **Working memory is rebuilt per call, never cached.** `graph_build_adj` runs fresh on
  every invocation; there is no memoization on the `Graph` object. For an algorithm
  called in a loop this is the dominant cost, and there is no existing mechanism to
  amortize it.
- **The copy-out ownership pattern avoids `SPEC.md` §4's trickiest rule.** Because
  results are built from `expr_copy` of borrowed vertices rather than transplanted
  sub-expressions, no graph builtin needs the NULL-out-before-free dance. Worth
  preserving in anything new.
- **Traversal duplication is normalized.** One shared component counter, but each file
  writes its own BFS/DFS. A new local traversal is idiomatic here, not technical debt.
- **The ladder's silence here is structural, not a verdict.** The audits cannot see this
  subsystem because it contains no NDArray dispatch markers — so "green" here means
  "never looked", which is a different claim from "checked and exempt". And the silence
  is only about the *consumer* side: `src/ndarray.c:547` records graph kernels as
  packed-buffer consumers with a 29× measured difference, so a head returning numeric
  vectors should be assumed to have obligations until the producer question is actually
  asked.

## Historical Context (from thoughts/)

`thoughts/` on this branch holds four documents, none about graphs:
`research/2026-08-18-unitbox-builtin.md`,
`research/2026-08-21-DEMO-1-nminimize-nmaximize-benchmarking.md`,
`plans/2026-08-21-DEMO-1-nminimize-nmaximize-benchmarks.md`, and the `DEMO-2` ticket
folder (an `NMinimize` feasibility fix, blocked at its plan gate on a missing
`## Plan Review` section).

The graph-relevant history is on **branch `rg-1-random-graph-count`**, not here:
`thoughts/shared/tickets/RG-1/` (six artifacts, see Finding 10). RG-1's own research
records that `src/graph/generators.c` had a single commit at that point — `56035303`,
the graph-subsystem bulk add — so the subsystem is essentially greenfield with one
feature landed on top.

RG-1 also recorded two measured defects that remain open and are worth knowing before
touching the generators: the candidate-edge tree leaks ~364 KB per `RandomGraph` call at
n=50, and `int_vertices`/`edges` C arrays leak in all four generators because
`expr_new_function` memcpys rather than adopting (`src/expr.c:257`). Both were
deliberately scoped out of RG-1 rather than fixed.

No subsystem doc exists at `thoughts/shared/subsystems/graph.md`. Given the Phase
convention, the registration idiom and the ladder-exemption question documented here,
writing one is a reasonable follow-on — RG-1's research made the same recommendation and
it has not been acted on.

## Related Research

- `thoughts/shared/tickets/RG-1/research.md` (branch `rg-1-random-graph-count`) — the
  only prior graph research; the generator family, the trailing-count convention, and the
  two open leaks.
- Nothing in `thoughts/shared/research/` on this branch touches the graph subsystem.
