---
ticket: RG-2
created: 2026-08-31T09:07:35-0400
researcher: Michael Sollami
topic: "How are graph algorithms structured in this repo — location, shared substrate, representation, tests, and the conformance surface for a new algorithm?"
type: research
lifecycle: active
full_research: thoughts/shared/tickets/RG-2/research.md
---

# Research Summary: graph algorithm structure in `src/graph/`

**Full research (appendix)**: `thoughts/shared/tickets/RG-2/research.md`

## Recommendation

A new graph algorithm is a **new `.c` file in `src/graph/` plus edits to four existing
files**, and the conformance surface is small and rigid enough to state exactly (below).
The subsystem is representation + traversal only — every classic combinatorial algorithm
is absent — so a search/optimization head would be the first to need an option, a
returned assignment, or a weight, and none of those has a precedent to copy.

## The conformance surface, concretely

**Seven registration sites across four files** (RG-1's plan covers four of them; sites
5–7 it omits because `RandomGraph` already existed):

1–3. `src/graph/graph.c` — `symtab_add_builtin`, then
`symtab_get_def("H")->attributes |= ATTR_PROTECTED` (**this idiom, never
`symtab_set_attributes`**), then `symtab_set_docstring`.
4. `src/graph/graph.h` — prototype under the Phase 5 "search & computation" banner (`:127`).
5–7. `src/sym_names.h:917-946` (`extern`), `src/sym_names.c:859-890` (`= NULL`),
`src/sym_names.c:1745-1775` (`intern_symbol`). **Omitting site 7 leaves a `NULL` that
identity-compares equal to nothing** — silent, because the subsystem compares interned
pointers, not strings.

`graph_init()` is already wired (`src/core.c:944-945`); no new init call.

**The algorithm skeleton** — validate → `graph_build_adj` → run over `calloc`'d scratch →
build a fresh result with `expr_copy` → free scratch + `graph_adj_free`. Canonical
example: `src/graph/shortestpath.c:49-71`.

**Non-negotiable conventions**, all verified as zero-exception across the directory:
`return NULL` is the only error channel (0 `Message()` calls), no head takes any option
(0 `Method`/`OptionValue` hits), `ATTR_PROTECTED` is the only attribute (27/27).

**Plus** `tests/test_graph.c` (target `graph_tests`, `tests/CMakeLists.txt:972-975`),
`docs/spec/builtins/graphs.md`, and `docs/spec/changelog/2026-08-31.md` — which already
exists, so append.

## Options Considered

Not a build-vs-build choice — these are the three genuine forks a new algorithm faces:

1. **Reuse `GraphAdj` vs. write a bespoke structure** — `graph_build_adj`
   (`graph_util.c:206-260`) gives int-indexed `out`/`in` adjacency, rebuilt per call.
   Tradeoff: free and idiomatic, but no caching exists, so a head called in a loop pays
   the rebuild every time and there is no mechanism to amortize it.
2. **Reuse `graph_count_components` vs. a local traversal** — one shared traversal
   exists; every file otherwise writes its own BFS/DFS, and `weak_components`
   (`components.c:51-73`) duplicates the shared one outright. Tradeoff: a new local
   traversal is idiomatic here, not debt — a reviewer should not read it as sloppiness.
3. **Declare a ladder exemption vs. stay silent** — see Decision Criteria.

## Decision Criteria

- **The ladder does not reach this subsystem on the CONSUMER side, and not by decision.**
  No graph head appears in `src/pack.c`'s `AWARE`, any `nd*` kernel file, or
  `src/compile/`. Three (`PathGraph`, `GraphQ`, `DirectedGraphQ`) sit in
  `nd_fastpath_sweep.py`'s `OFF_BUFFER`, documented at `:50-52` as a **backlog ratchet,
  not an exemption**. `AdjacencyMatrix` is named by one audit and satisfied by none. The
  other 23 are outside every audit's candidate pool because `check_packed_aware.py` only
  nominates heads containing an NDArray dispatch marker, and `src/graph/` has none. So
  "green" here means "never looked" — a different claim from "checked and exempt". No
  graph head has a written justification, unlike the `GEO-1` precedent
  (`check_packed_aware.py:79-104`).
- **But the PRODUCER side is unexamined, and the evidence points against exemption.**
  `SPEC.md` §9's separate `--survival` question — which producers return a plain `List`
  that would have packed — was never asked. `AdjacencyMatrix`, `ConnectedComponents`,
  `VertexDegree` and `GraphDistance` all have that shape, building results one
  `expr_new_integer` at a time (`src/graph/adjmat.c:54`). And `src/ndarray.c:540-547`
  records the packed gather as "the operation every sparse-matrix and **graph** kernel is
  built out of", with a measured **PageRank of 14.1 s against 486 ms**. So do not plan on
  the assumption that graphs are exempt — that conclusion is not supported.
- **The gap inventory sets the bar.** All 25 classic algorithms probed are absent —
  colouring, cliques, flow, matching, Hamiltonian/Eulerian, topological sort, centrality,
  isomorphism, planarity. The only precedent for expensive combinatorial search is
  `VertexConnectivity`'s brute-force subset enumeration (`connectivity.c:32-52`): honest,
  but not a scalable model.
- **Edge weights do not exist.** One aspirational comment (`adjmat.c:9-10`), nothing
  else. A weighted algorithm expands into a subsystem-wide representation change.
- **RG-1 is the procedural template** (your call): seven-file touch table, an Acceptance
  Criteria table with executable rows, an `## Entry Points` section cross-checked against
  those rows, and a validation doc that reports `NOT RUN`/`NOT VERIFIED` honestly rather
  than rounding up.

## Open Questions

### Unresolved

_None._

### Resolved

- [x] Does a numeric-returning graph head owe a packed producer path? — **By return
  type.** A `Graph` object or an arbitrary-expression vertex list carries no obligation; a
  dense numeric vector or matrix does. So the return type must be settled in the plan's
  Decisions before the ladder obligation can be stated. _(picked "Not applicable —
  depends on return type")_
- [x] `AdjacencyMatrix` NUMERIC but unpacked? — **Real gap, record it.** `NUMERIC` is the
  right classification; the missing packed path is a genuine unrecorded gap, not a
  tooling error to paper over. _(picked "Real gap — record it")_
- [x] Edge weights? — **Separate prior change.** Own ticket, touching `AdjacencyMatrix`,
  `graph_build_adj` and the canonical `Graph[...]` form. Blocks weighted work; unweighted
  work unaffected. _(picked "Separate prior change")_
- [x] The eight backwards-reading "frees res" headers? — **Leave them**, separate change;
  the research records which eight and points a new file at `graph.h:21`. _(picked
  "Leave them — separate change")_
- [x] Prior/abandoned graph work not obvious from git? — RG-1 is to be treated as the
  conformance template, not just history. _(chose "Treat RG-1 as the template to conform to")_
- [x] How deep into the numeric-surface ladder? — Fully; done. _(chose "Trace it fully for graph heads")_
- [x] Anything out of scope inside `src/graph/`? — `graphplot.c`/rendering, excluded. _(chose "graphplot.c / rendering")_
- [x] Build and probe the binary? — Only if current; it was **stale** (`src/facint.c`,
  `src/info.c`, `src/geometry.c` newer than `./Mathilda`), so nothing was probed and all
  findings are source-derived. _(chose "Build only if it's already current")_

## Requires Approval

Two scope calls before a plan is built on this. **Edge weights** — the absence is
subsystem-wide, so a weighted algorithm silently becomes a representation change.
**The ladder** — this research establishes graph heads are undocumented and unaudited on
the consumer side, but it does **not** establish that they are exempt. The approval needed
is the prior question — does a new numeric-returning head owe a packed producer path? —
not "add an `EXEMPT` line", which would settle it in the direction the evidence does not
support.

## Findings recorded, not acted on

- **Under `NDEBUG` the graph suite exits 0 on failure.** `assert_eval_eq` uses libc
  `assert()`, a no-op under `-DNDEBUG` (CMake Release); `tests/test_utils.h:54-57` warns
  about it. The failure text *is* still printed to stderr (`:24-28`) — so a human reading
  the log sees it — but the process exits 0 and **CTest reports green on a failing
  suite**. Only `test_inputform_roundtrip` and `test_random_graph` use the surviving
  `ASSERT` macro. Fix is wrapping the call, not rewriting assertions.
- **Half the file headers state the ownership rule backwards.** Eight files end with a
  bare "…; frees res" (`shortestpath.c:11`, `components.c:12`, `adjmat.c:12`,
  `adjgraph.c:12`, `incmat.c:7`, `spanningtree.c:8`, `connectivity.c:12`,
  `generators.c:15`), which reads as the builtin freeing `res` — the opposite of the code
  and of `SPEC.md` §4. The other eight say it correctly. Since RG-1's touch surface
  includes mirroring the header comment, a new file should copy `graph.h:21`, not its
  neighbour.
- **Empty graphs are tested nowhere.** No assertion in `tests/test_graph.c` constructs
  `Graph[{},{}]` or a zero-vertex graph, though the code guards for it throughout via
  the `calloc(n > 0 ? n : 1)` idiom.
