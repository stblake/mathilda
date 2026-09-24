# Hypergraphs

A native hypergraph subsystem, implemented in `src/graph/hyp_*.c` (header
`src/graph/graph_hyper.h`, registered by `graph_hyper_init()` at the end of
`graph_init()`).

Mathematica 15 has **no** built-in hypergraph symbols (`Names["*Hypergraph*"]`
is empty). Hypergraphs exist in the Wolfram ecosystem only as add-ons: the
WolframInstitute/Hypergraph paclet (a `Hypergraph` object with `VertexDegree`,
`VertexList`, `EdgeList`, `HypergraphIncidenceMatrix`, …), Wolfram Function
Repository functions (`ResourceFunction["HypergraphToGraph"]`,
`["ConnectedHypergraphQ"]`, `["RandomHypergraph"]`, `["TransversalHypergraph"]`,
`["HypergraphAdjacencyMatrix"]`, …) and the Wolfram Physics Project tooling
(`["WolframModel"]`, `["HypergraphPlot"]`). All of those run as top-level WL.
Mathilda's heads are named after them where a Wolfram name exists
(`Hypergraph`, `HypergraphQ`, `HypergraphToGraph`, `ConnectedHypergraphQ`,
`RandomHypergraph`, `TransversalHypergraph`), accept their argument forms
(including a bare List of hyperedges where the FR function takes one), and are
a strict superset in capability — see the per-head notes for where the FR
functions are narrower (they leave several degenerate inputs unevaluated).

## Representation

Hypergraphs are ordinary `Expr` trees, like `Graph` — no new `EXPR_*` tag:

```
Hypergraph[ List[v1, ..., vn], List[e1, ..., em] ]
```

Each hyperedge `e_j` is a `List` of vertices. Vertices are arbitrary
expressions and pairwise distinct. Hyperedges may **repeat** (a
multi-hypergraph), **overlap** or **nest** arbitrarily, be **empty**, and may
**repeat a vertex** (`{1, 1, 2}`, as Wolfram-model states do).

**Ordered vs. unordered — the documented choice.** A hyperedge `List` is kept
exactly as written, so its order is available to the order-sensitive heads:
`HypergraphToGraph` (the FR convention: ordered hyperedge `{a, b, c}` gives
`a->b, a->c, b->c`), `EdgeList`/`InputForm` round-trip, `HyperedgeSizes` /
`HypergraphRank` / `HypergraphCorank` / `UniformHypergraphQ` (arity = `Length`,
counting repeats), and `HypergraphEdgeDelete` (`SameQ` on the List). Every
**set-theoretic** head — degrees, incidence, duals, expansions, line graphs,
components, distances, restrictions, transversals — reads a hyperedge as the
**set of its distinct vertices**. So a Wolfram-model state such as
`{{1,2,3},{3,4}}` is a valid hyperedge list as it stands, and one object answers
both kinds of question. (There is no separate directed-hyperedge head: an
ordered List *is* the Wolfram-model directed hyperedge.)

Hyperedge indices, in arguments and results, are 1-based positions in
`EdgeList[h]`.

### Performance model: the validated-hypergraph memo

`src/graph/hyp_util.c` re-implements `graph_util.c`'s memo for hypergraphs (the
Graph memo is Graph-specific and was not changed): the last 4 valid hypergraph
**nodes**, keyed by pointer, each slot holding a reference (so the node stays
alive, hence unrecyclable and immutable — mutators unshare it first). A slot
keeps the vertex hash index, the hyperedges as a vertex-index CSR (raw, and
de-duplicated when some hyperedge repeats a vertex), and — built lazily on
first need — the vertex→hyperedge incidence CSR. `Hypergraph[...]` seeds the
memo during construction. So validation is `O(1)` after the first call on an
object, and every algorithm runs on integer arrays with no hashing. Invalid
objects are never memoized. When every vertex is a machine integer in a range
at most `4n + 1024` wide (e.g. `Range[n]`, Wolfram-model states), hyperedge
elements resolve by direct addressing instead of through the hash index, which
roughly halves construction time.

Algorithms are linear in the total incidence `Σ|e|`, except the s-overlap family
(`HypergraphLineGraph`, `HyperedgeConnectedComponents[h, s ≥ 2]`, s-walk
distances with `s ≥ 2`), which is `O(Σ_v deg(v)²)` — the number of hyperedge
pairs meeting at a vertex, which bounds the 1-line graph's size anyway.
`HyperedgeDistance[h, i]` with `s = 1` is linear: each vertex's incidence list is
expanded at most once per BFS.

## Hypergraph
- `Hypergraph[{e1, e2, ...}]`: vertices derived in first-appearance order.
- `Hypergraph[{v1, ...}, {e1, ...}]`: explicit vertices (duplicates dropped,
  first occurrence kept; isolated vertices allowed). Every element of every
  hyperedge must be one of them.
- `Hypergraph[g]`: converts a `Graph`, each edge `u<->v`/`u->v` becoming `{u, v}`.
- Anything malformed (a non-List hyperedge, an unknown vertex, other arities)
  is left unevaluated. A valid hypergraph is its own fixed point.
- Standard output is the terse summary `Hypergraph[<n vertices, m hyperedges>]`;
  `InputForm`/`FullForm` print the literal, which round-trips through the parser.

```
In:  h = Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]
Out: Hypergraph[<7 vertices, 4 hyperedges>]
In:  InputForm[h]
Out: Hypergraph[{1, 2, 3, 4, 5, 6, 7}, {{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
```

(`h` below always refers to this hypergraph.)

## HypergraphQ
`HypergraphQ[h]` gives `True` for a valid Hypergraph, `False` otherwise
(including for a `Graph`; `GraphQ` of a hypergraph is `False`).

## Shared accessors
The Graph accessors dispatch to the hypergraph code when handed a Hypergraph;
Graph behaviour is unchanged.
- `VertexList[h]`, `EdgeList[h]` (the hyperedge Lists as written),
  `VertexCount[h]`, `EdgeCount[h]`.
- `VertexDegree[h]` / `VertexDegree[h, v]`: the number of hyperedges containing
  the vertex; repeated hyperedges count separately, a vertex repeated inside one
  hyperedge counts once. Packed above the packing threshold.
  `VertexInDegree`/`VertexOutDegree` are Graph-only (unevaluated).
- `IncidenceMatrix[h]`: the `n × m` matrix whose `(i, j)` entry is the number of
  times vertex `i` occurs in hyperedge `j` — 0/1 for set-like hyperedges, the
  multiplicity otherwise (as `IncidenceMatrix` gives 2 for a Graph self-loop).
  Dense and packed (Mathilda has no `SparseArray`).

```
In:  VertexDegree[h]
Out: {1, 1, 2, 2, 1, 1, 1}
In:  IncidenceMatrix[Hypergraph[{{1,2},{2,3,3}}]]
Out: {{1, 0}, {1, 1}, {0, 2}}
```

## HyperedgeSizes / HypergraphRank / HypergraphCorank / UniformHypergraphQ
- `HyperedgeSizes[h]`: the arity (`Length`) of each hyperedge, packed.
- `HypergraphRank[h]`, `HypergraphCorank[h]`: largest / smallest arity (0 when
  there are no hyperedges).
- `UniformHypergraphQ[h]`: all arities equal (`True` with no hyperedges).
  `UniformHypergraphQ[h, k]`: `h` is k-uniform. `False` for a non-hypergraph.

```
In:  {HyperedgeSizes[h], HypergraphRank[h], HypergraphCorank[h], UniformHypergraphQ[h]}
Out: {{3, 2, 3, 1}, 3, 1, False}
```

## HypergraphDual
`HypergraphDual[h]`: vertices `1..m` (one per hyperedge), and for each vertex of
`h` in VertexList order the hyperedge of the (ascending) indices of the
hyperedges containing it. The dual of the dual recovers the incidence structure
on `1..n`.

```
In:  InputForm[HypergraphDual[h]]
Out: Hypergraph[{1, 2, 3, 4}, {{1}, {1}, {1, 2}, {2, 3}, {3}, {3}, {4}}]
```

## HypergraphCliqueExpansion / HypergraphStarExpansion / HypergraphToGraph
All three return a `Graph` (simple, so it is validated and memoized as usual).
- `HypergraphCliqueExpansion[h]`: the 2-section — `u <-> v` whenever `u` and `v`
  share a hyperedge; edges in order of first co-occurrence.
- `HypergraphStarExpansion[h]`: the incidence (star) bipartite graph on
  `VertexList[h]` plus nodes `Hyperedge[1], ..., Hyperedge[m]`, with
  `v <-> Hyperedge[j]` for `v ∈ e_j`. Unevaluated if some vertex of `h` is itself
  such a `Hyperedge[j]` node.
- `HypergraphToGraph[h]` (FR name and semantics): `h` read as an ordered
  hypergraph; hyperedge `{v1, ..., vk}` contributes `v_a -> v_b` for every
  `a < b`. Accepts a bare List of hyperedges as the FR function does.
  *Deviation:* Mathilda graphs are simple, so self-loops (from a repeated vertex)
  are dropped and parallel copies merged, and every vertex is kept (the FR
  function returns a multigraph and drops vertices that occur only in unary
  hyperedges).

```
In:  InputForm[HypergraphCliqueExpansion[h]]
Out: Graph[{1, 2, 3, 4, 5, 6, 7}, {1 <-> 2, 1 <-> 3, 2 <-> 3, 3 <-> 4, 4 <-> 5, 4 <-> 6, 5 <-> 6}]
In:  InputForm[HypergraphToGraph[{{1,2,3},{3,4}}]]
Out: Graph[{1, 2, 3, 4}, {1 -> 2, 1 -> 3, 2 -> 3, 3 -> 4}]
```

## HypergraphLineGraph
`HypergraphLineGraph[h]`: Graph on `1..m`, `i <-> j` when hyperedges `i` and `j`
intersect. `HypergraphLineGraph[h, s]`: the **s-line graph** (Aksoy, Joslyn,
Ortiz Marrero, Praggastis, Purvine, "Hypernetwork science via high-order
hypergraph walks", 2020): `i <-> j` when they share at least `s` vertices.
Edges ordered by `i`, then `j`.

```
In:  EdgeList[HypergraphLineGraph[Hypergraph[{{1,2,3},{2,3,4},{3,4,5},{1,5}}], 2]]
Out: {1 <-> 2, 2 <-> 3}
```

## HypergraphConnectedComponents / ConnectedHypergraphQ
- `HypergraphConnectedComponents[h]`: vertex components (a chain of hyperedges
  joins two vertices), as vertex Lists; components ordered by their first
  vertex (Mathilda's `ConnectedComponents` convention), vertices in VertexList
  order. Isolated vertices are singleton components. Union–find, near-linear.
- `ConnectedHypergraphQ[h]` (FR name): `True` iff `h` has at least one vertex and
  one component. Accepts a bare List (the FR function leaves `{}` unevaluated;
  this gives `False`).

```
In:  HypergraphConnectedComponents[h]
Out: {{1, 2, 3, 4, 5, 6}, {7}}
```

## HyperedgeConnectedComponents
`HyperedgeConnectedComponents[h]` / `[h, s]`: the **s-connected components**
(Aksoy et al.) as Lists of 1-based hyperedge indices — hyperedges with at least
`s` distinct vertices, joined when they share at least `s`, closed
transitively. Hyperedges with fewer than `s` vertices lie on no s-walk and are
omitted (so empty hyperedges never appear). `s = 1` is union–find over the
incidence lists; `s ≥ 2` counts overlaps per hyperedge.

```
In:  HyperedgeConnectedComponents[Hypergraph[{{1,2,3},{2,3,4},{3,4,5},{1,5}}], 2]
Out: {{1, 2, 3}, {4}}
```

## HyperedgeDistance / HypergraphDistance
- `HyperedgeDistance[h, i, j]` / `[h, i, j, s]`: the length of a shortest s-walk
  (consecutive hyperedges sharing at least `s` vertices) from hyperedge `i` to
  hyperedge `j`; `Infinity` if none; `0` for `i == j`.
  `HyperedgeDistance[h, i]` / `[h, i, All, s]`: distances to every hyperedge
  (packed when all finite). BFS directly on the incidence structure — the line
  graph is never materialised.
- `HypergraphDistance[h, u, v]` / `[h, u]`: least number of hyperedges in a chain
  of vertices from `u` to `v` (the distance in the clique expansion), `Infinity`
  if none, in VertexList order for the one-source form.

```
In:  {HyperedgeDistance[h, 1], HypergraphDistance[h, 1]}
Out: {{0, 1, 2, Infinity}, {0, 1, 1, 2, 3, 3, Infinity}}
```

## Edits and sub-hypergraphs
In the Edge heads, a List whose every element is a List is a list of hyperedges;
otherwise it is one hyperedge. In the Vertex heads any List is a list of
vertices (wrap a List-valued vertex as `{{...}}`).
- `HypergraphVertexAdd[h, v | {v1, ...}]`: appends vertices not already present.
- `HypergraphVertexDelete[h, v | {v1, ...}]`: removes the vertices **and every
  hyperedge containing one of them**, as `VertexDelete` does for a Graph.
  Unevaluated if a named vertex is absent.
- `HypergraphEdgeAdd[h, e | {e1, ...}]`: appends hyperedges (repeats allowed),
  adding new vertices in first-appearance order.
- `HypergraphEdgeDelete[h, e | {e1, ...}]`: removes every hyperedge `SameQ` to a
  named one (as written: `{2, 1}` does not delete `{1, 2}`). Unevaluated if a
  named hyperedge does not occur.
- `Subhypergraph[h, {v1, ...}]`: `Subgraph`'s rule — the named vertices that are
  in `h` (VertexList order) and the hyperedges lying entirely among them.
- `HypergraphRestriction[h, {v1, ...}]`: Berge's induced sub-hypergraph — every
  hyperedge intersected with the vertex set (order and repeats of kept vertices
  preserved), hyperedges that miss it dropped.

```
In:  InputForm[HypergraphVertexDelete[h, 3]]
Out: Hypergraph[{1, 2, 4, 5, 6, 7}, {{4, 5, 6}, {7}}]
In:  InputForm[HypergraphRestriction[h, {1, 2, 3, 4}]]
Out: Hypergraph[{1, 2, 3, 4}, {{1, 2, 3}, {3, 4}, {4}}]
```

## RandomHypergraph
- `RandomHypergraph[{n, m}, k]`: a k-uniform hypergraph on `1..n` with `m`
  independent, uniformly random k-subsets (listed increasing) as hyperedges;
  a hyperedge may repeat. `0 ≤ k ≤ n`, else unevaluated. Floyd's sampling, O(k)
  per hyperedge. `RandomHypergraph[{n, m}, k, c]` gives `c` of them.
- `RandomHypergraph[{n, {e, a}}]`, `RandomHypergraph[{n, {{e1, a1}, ...}}]`: the FR
  form — `e_i` hyperedges of arity `a_i`, vertices drawn uniformly **with
  replacement** from `1..n` (so a vertex can repeat inside a hyperedge, as in
  the FR output `{8, 1, 1}`); vertex list = the vertices that occur.
  *Deviation:* returns a `Hypergraph` where the FR function returns the bare
  edge List (`EdgeList` recovers it).
- Both draw from the user-visible stream, so `SeedRandom` reproduces them.

```
In:  SeedRandom[1]; EdgeList[RandomHypergraph[{10, 4}, 3]]
Out: {{2, 7, 9}, {2, 6, 10}, {1, 5, 8}, {2, 4, 9}}
```

## TransversalHypergraph / FindMinimumTransversal
A transversal (hitting set) meets every hyperedge; it is minimal when no proper
subset is one.
- `TransversalHypergraph[h]` (FR name): every **minimal** transversal —
  `Hypergraph[VertexList[h], Tr]` for a Hypergraph, the List `Tr` for a bare List
  of hyperedges (non-List elements, the FR's "isolated vertices", are
  ignored). No hyperedges → `{{}}`; an empty hyperedge → `{}`. The FR function
  leaves `{}`, `{{}}`, `{{1},{1,2}}` and `{{1,2},{1,2}}` unevaluated; these are
  answered. Each transversal lists its vertices in VertexList order;
  transversals are sorted by size, then lexicographically by position.
  Algorithm: **MMCS** (Murakami & Uno, 2014), run iteratively, with `O(deg v)`
  critical-hyperedge bookkeeping (per hyperedge, `|F ∩ S|` and the sum of `S`'s
  members in `F`, which *is* the unique member when the count is 1).
- `FindMinimumTransversal[h]`: one minimum-cardinality transversal (NP-hard).
  Exact branch and bound: branch on the uncovered hyperedge with fewest open
  vertices (include, then exclude for later siblings); prune with the best of a
  degree bound, a low-degree-first disjoint-hyperedge packing, and a
  fractional-packing (LP-dual) bound with one dual-ascent pass; seeded by the
  greedy max-coverage bound (lazy heap). Unevaluated if a hyperedge is empty.

**Budgets (correct-or-unevaluated).** Both searches count nodes and return
unevaluated — never a partial or unproven answer — past 5·10⁷ MMCS nodes or
2·10⁷ output entries (TransversalHypergraph), or 2·10⁷ branch-and-bound nodes
(FindMinimumTransversal, which also refuses a greedy bound above 20000 to cap
C-stack depth). Node counts, not wall clock, keep answers machine-independent;
both poll `tc_check_deadline()` every 4096 nodes, so `TimeConstrained` interrupts
them.

```
In:  TransversalHypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]
Out: {{1, 4, 7}, {2, 4, 7}, {3, 4, 7}, {3, 5, 7}, {3, 6, 7}}
In:  FindMinimumTransversal[{{1,2,3},{3,4},{4,5,6},{7}}]
Out: {3, 4, 7}
```

## Benchmarks (experiment 96)
`benchmarks/96-hypergraphs/` (Apple silicon, times in ms, best of 3; cold = a
fresh hypergraph object before every timed call, construction untimed).
Mathematica column: the WolframInstitute/Hypergraph paclet, the Function
Repository ResourceFunctions, or plain WL where neither has the function —
each definition in `mma_baseline.m` is tagged with its source (its
single-source distances use `BreadthFirstScan`, since the kernel's
`GraphDistance[g, s]` is quadratic here). Python column: xgi 0.10, networkx
BFS where xgi's own shortest paths are quadratic, scipy HiGHS for the MILP.
Every case's value check agrees across all three systems.

| Case | Mathilda | Mathematica | xgi | Mma / Mathilda | xgi / Mathilda |
|---|---:|---:|---:|---:|---:|
| Hypergraph construction, 10^5 hyperedges | 19.9 | 2912 | 301 | 146× | 15× |
| VertexDegree, 10^5 hyperedges | 0.005 | 577 | 1.8 | 115373× | 354× |
| HyperedgeSizes, 10^5 hyperedges | 0.024 | 5.6 | 9.6 | 234× | 401× |
| HypergraphDual, 10^5 hyperedges | 8.1 | 236 | 361 | 29× | 44× |
| HypergraphCliqueExpansion, 10^5 hyperedges | 32.3 | 278 | 729 | 8.60× | 23× |
| HypergraphStarExpansion, 10^5 hyperedges | 38.5 | 313 | 466 | 8.13× | 12× |
| HypergraphToGraph (FR), 10^5 hyperedges | 59.3 | 233 | 376 | 3.93× | 6.34× |
| HypergraphLineGraph, 10^5 hyperedges | 47.0 | 819 | 1117 | 17× | 24× |
| HypergraphConnectedComponents, 10^5 hyperedges | 2.3 | 160 | 282 | 70× | 123× |
| ConnectedHypergraphQ (FR), 10^5 hyperedges | 20.1 | 261 | 546 | 13× | 27× |
| HyperedgeConnectedComponents s=2, 10^5 hyperedges | 18.4 | 1306 | 2466 | 71× | 134× |
| HyperedgeDistance from e1, 10^5 hyperedges | 2.8 | 7255 | 1048 | 2587× | 374× |
| HypergraphDistance from v1, 10^5 hyperedges | 2.9 | 5791 | 1057 | 2023× | 369× |
| IncidenceMatrix, 1000 x 10^4 | 0.365 | 37.1 | 5.6 | 102× | 15× |
| RandomHypergraph (FR form), 10^5 triples | 10.6 | 1.6 | 438 | 0.15× | 41× |
| TransversalHypergraph (FR), 16 vertices | 4.2 | 5219 | n/a | 1235× | n/a |
| FindMinimumTransversal, 80 vertices x 200 | 245 | 874 | 784 | 3.56× | 3.20× |
| VertexDegree, 10^5 hyperedges (cold) | 0.432 | 559 | 2.4 | 1293× | 5.48× |
| HyperedgeSizes, 10^5 hyperedges (cold) | 0.025 | 5.6 | 10.4 | 224× | 415× |
| HypergraphDual, 10^5 hyperedges (cold) | 10.1 | 221 | 364 | 22× | 36× |
| HypergraphCliqueExpansion, 10^5 hyperedges (cold) | 33.8 | 273 | 851 | 8.09× | 25× |
| HypergraphStarExpansion, 10^5 hyperedges (cold) | 38.3 | 427 | 577 | 11× | 15× |
| HypergraphLineGraph, 10^5 hyperedges (cold) | 48.5 | 873 | 1160 | 18× | 24× |
| HypergraphConnectedComponents, 10^5 hyperedges (cold) | 2.3 | 177 | 238 | 79× | 106× |
| HyperedgeConnectedComponents s=2, 10^5 hyperedges (cold) | 18.6 | 1357 | 2862 | 73× | 154× |
| HyperedgeDistance from e1, 10^5 hyperedges (cold) | 3.3 | 7388 | 1241 | 2257× | 379× |
| HypergraphDistance from v1, 10^5 hyperedges (cold) | 3.5 | 6068 | 1118 | 1737× | 320× |
| IncidenceMatrix, 1000 x 10^4 (cold) | 0.498 | 40.3 | 5.7 | 81× | 11× |

The one loss is `RandomHypergraph` in the FR form: the ResourceFunction returns
a single packed `RandomInteger` array, while Mathilda returns a validated
`Hypergraph` whose hyperedges are 10⁵ boxed `List` nodes (Mathilda's
representation invariant keeps packed buffers out of expression trees); the
allocation of those nodes alone exceeds 1.6 ms. For scale, Mathilda's own
`RandomInteger[{1, 20000}, {100000, 3}]` takes 16 ms.
`FindMinimumTransversal`'s baselines are MILP solvers (Mathematica's
`LinearOptimization`, scipy/HiGHS); the combinatorial branch and bound is
3.2–3.6× faster at this size, and comparable to the MILPs at 100 vertices.
