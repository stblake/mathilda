# Hypergraphs

A native hypergraph subsystem, implemented in `src/graph/hyp_*.c` (header
`src/graph/graph_hyper.h`, registered by `graph_hyper_init()` at the end of
`graph_init()`).

### Relation to the Wolfram ecosystem

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

### Representation

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

### Shared accessors

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

Examples of these accessors are on the `Hypergraph` page.

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

### Benchmarks (experiment 96)

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

In the examples on the individual pages, `h` denotes
`Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]` (7 vertices, 4 hyperedges).

## Hypergraph

- `Hypergraph[{e1, e2, ...}]`: the hypergraph with hyperedges `e_i`; vertices
  derived in first-appearance order.
- `Hypergraph[{v1, ...}, {e1, ...}]`: explicit vertices (duplicates dropped,
  first occurrence kept; isolated vertices allowed). Every element of every
  hyperedge must be one of them.
- `Hypergraph[g]`: converts a `Graph`, each edge `u<->v`/`u->v` becoming `{u, v}`.

**Features**:
- Hypergraphs are ordinary `Expr` trees of the canonical form
  `Hypergraph[{v1, ..., vn}, {e1, ..., em}]` — no new `EXPR_*` tag. Each
  hyperedge is a `List` of vertices; vertices are arbitrary, pairwise distinct
  expressions. Hyperedges may repeat (a multi-hypergraph), overlap or nest, be
  empty, and may repeat a vertex (`{1, 1, 2}`, as Wolfram-model states do).
- A hyperedge List is kept exactly as written: order-sensitive heads
  (`HypergraphToGraph`, `EdgeList`/`InputForm`, `HyperedgeSizes` and the arity
  heads, `HypergraphEdgeDelete`) see the order and repeats; every set-theoretic
  head reads a hyperedge as the set of its distinct vertices (see the category
  preamble).
- Anything malformed (a non-List hyperedge, an unknown vertex, other arities)
  is left unevaluated. A valid hypergraph is its own fixed point.
- Standard output is the terse summary `Hypergraph[<n vertices, m hyperedges>]`;
  `InputForm`/`FullForm` print the literal, which round-trips through the parser.
- Construction seeds the validated-hypergraph memo (vertex hash index plus a
  vertex-index CSR of the hyperedges), so validation by later heads is `O(1)`.
  When every vertex is a machine integer in a range at most `4n + 1024` wide,
  elements resolve by direct addressing, roughly halving construction time.
- Shared Graph accessors accept a Hypergraph: `VertexList`, `EdgeList` (the
  hyperedge Lists as written), `VertexCount`, `EdgeCount`; `VertexDegree[h]` /
  `VertexDegree[h, v]` (hyperedges containing the vertex — repeated hyperedges
  count separately, a vertex repeated inside one hyperedge counts once; packed
  above the packing threshold); `IncidenceMatrix[h]` (the dense, packed
  `n × m` matrix of occurrence counts — 0/1 for set-like hyperedges, the
  multiplicity otherwise, as `IncidenceMatrix` gives 2 for a Graph self-loop;
  Mathilda has no `SparseArray`). `VertexInDegree`/`VertexOutDegree` are
  Graph-only and stay unevaluated.
- Mathematica has no built-in `Hypergraph`; the name follows the
  WolframInstitute/Hypergraph paclet.
- Benchmark (experiment 96): construction of 10⁵ hyperedges 19.9 ms (Mathematica
  paclet 2912 ms, xgi 301 ms); `VertexDegree` 0.005 ms warm / 0.432 ms cold;
  `IncidenceMatrix` 1000 × 10⁴ 0.365 ms warm / 0.498 ms cold.

```mathematica
In[1]:= Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]
Out[1]= Hypergraph[<7 vertices, 4 hyperedges>]

In[2]:= InputForm[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]
Out[2]= Hypergraph[{1, 2, 3, 4, 5, 6, 7}, {{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]

In[3]:= InputForm[Hypergraph[{1,2,2,3},{{1,2}}]]
Out[3]= Hypergraph[{1, 2, 3}, {{1, 2}}]

In[4]:= InputForm[Hypergraph[Graph[{1<->2, 2->3}]]]
Out[4]= Hypergraph[{1, 2, 3}, {{1, 2}, {2, 3}}]

In[5]:= VertexDegree[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]
Out[5]= {1, 1, 2, 2, 1, 1, 1}

In[6]:= IncidenceMatrix[Hypergraph[{{1,2},{2,3,3}}]]
Out[6]= {{1, 0}, {1, 1}, {0, 2}}
```

An unknown vertex in a hyperedge leaves the call unevaluated:

```mathematica
In[7]:= InputForm[Hypergraph[{1,2},{{1,5}}]]
Out[7]= Hypergraph[{1, 2}, {{1, 5}}]
```

## HypergraphQ

- `HypergraphQ[expr]`: gives `True` if `expr` is a valid `Hypergraph`, and
  `False` otherwise.

**Features**:
- `False` for a `Graph`, for a bare List of hyperedges, and for a malformed
  (unevaluated) `Hypergraph[...]`. Conversely `GraphQ` of a hypergraph is `False`.
- `O(1)` on a memoized hypergraph object.

```mathematica
In[1]:= HypergraphQ[Hypergraph[{{1,2,3},{3,4}}]]
Out[1]= True

In[2]:= HypergraphQ[Graph[{1<->2}]]
Out[2]= False

In[3]:= GraphQ[Hypergraph[{{1,2,3},{3,4}}]]
Out[3]= False

In[4]:= HypergraphQ[Hypergraph[{1,2},{{1,5}}]]
Out[4]= False
```

## HyperedgeSizes / HypergraphRank / HypergraphCorank / UniformHypergraphQ

- `HyperedgeSizes[h]`: the arity (`Length`) of each hyperedge, as a packed list.
- `HypergraphRank[h]`: the largest hyperedge arity.
- `HypergraphCorank[h]`: the smallest hyperedge arity.
- `UniformHypergraphQ[h]`: `True` if all hyperedges have the same arity.
- `UniformHypergraphQ[h, k]`: `True` if `h` is k-uniform.

**Features**:
- Arity is the hyperedge's `Length` as written, counting a repeated vertex and
  giving 0 for an empty hyperedge.
- With no hyperedges, `HypergraphRank` and `HypergraphCorank` are 0 and
  `UniformHypergraphQ` is `True`.
- `UniformHypergraphQ` gives `False` for a non-hypergraph; the other three are
  left unevaluated on one.
- Benchmark (experiment 96): `HyperedgeSizes` on 10⁵ hyperedges 0.024 ms warm /
  0.025 ms cold (Mathematica 5.6 ms, xgi 9.6 ms).

```mathematica
In[1]:= HyperedgeSizes[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]
Out[1]= {3, 2, 3, 1}

In[2]:= {HypergraphRank[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]], HypergraphCorank[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]}
Out[2]= {3, 1}

In[3]:= UniformHypergraphQ[Hypergraph[{{1,2,3},{2,3,4}}], 3]
Out[3]= True

In[4]:= HyperedgeSizes[Hypergraph[{{1,1,2},{}}]]
Out[4]= {3, 0}

In[5]:= {HypergraphRank[Hypergraph[{1,2},{}]], HypergraphCorank[Hypergraph[{1,2},{}]], UniformHypergraphQ[Hypergraph[{1,2},{}]]}
Out[5]= {0, 0, True}

In[6]:= HyperedgeSizes[{{1, 2}}]
Out[6]= HyperedgeSizes[{{1, 2}}]
```

## HypergraphDual

- `HypergraphDual[h]`: the dual hypergraph, with one vertex per hyperedge of
  `h` and one hyperedge per vertex of `h`.

**Features**:
- Vertices are `1..m`; for each vertex of `h`, in VertexList order, the dual has
  the hyperedge of the (ascending) indices of the hyperedges containing it. An
  isolated vertex gives an empty hyperedge.
- Hyperedges are read as sets. The dual of the dual recovers the incidence
  structure on `1..n`.
- Linear in the total incidence. Benchmark (experiment 96): 10⁵ hyperedges
  8.1 ms warm / 10.1 ms cold (Mathematica 236 ms, xgi 361 ms).
- A bare List of hyperedges is not accepted (unevaluated).

```mathematica
In[1]:= InputForm[HypergraphDual[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]]
Out[1]= Hypergraph[{1, 2, 3, 4}, {{1}, {1}, {1, 2}, {2, 3}, {3}, {3}, {4}}]

In[2]:= InputForm[HypergraphDual[HypergraphDual[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]]]
Out[2]= Hypergraph[{1, 2, 3, 4, 5, 6, 7}, {{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]

In[3]:= InputForm[HypergraphDual[Hypergraph[{a,b,c},{{a,b},{a,b}}]]]
Out[3]= Hypergraph[{1, 2}, {{1, 2}, {1, 2}, {}}]

In[4]:= HypergraphDual[{{1, 2}}]
Out[4]= HypergraphDual[{{1, 2}}]
```

## HypergraphCliqueExpansion

- `HypergraphCliqueExpansion[h]`: the 2-section of `h` — the `Graph` with
  `u <-> v` whenever `u` and `v` share a hyperedge.

**Features**:
- Returns a simple `Graph` (validated and memoized as usual) on `VertexList[h]`;
  edges in order of first co-occurrence. A vertex repeated inside a hyperedge
  gives no self-loop.
- Also equals the graph whose distance `HypergraphDistance` measures.
- Benchmark (experiment 96): 10⁵ hyperedges 32.3 ms warm / 33.8 ms cold
  (Mathematica 278 ms, xgi 729 ms).

```mathematica
In[1]:= InputForm[HypergraphCliqueExpansion[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]]
Out[1]= Graph[{1, 2, 3, 4, 5, 6, 7}, {1 <-> 2, 1 <-> 3, 2 <-> 3, 3 <-> 4, 4 <-> 5, 4 <-> 6, 5 <-> 6}]

In[2]:= InputForm[HypergraphCliqueExpansion[Hypergraph[{{1,1,2},{2,3}}]]]
Out[2]= Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3}]

In[3]:= HypergraphCliqueExpansion[{{1, 2}}]
Out[3]= Graph[<2 vertices, 1 edge>]
```

## HypergraphStarExpansion

- `HypergraphStarExpansion[h]`: the incidence (star) bipartite `Graph` of `h`.

**Features**:
- Vertices are `VertexList[h]` followed by nodes `Hyperedge[1], ...,
  Hyperedge[m]`, with an edge `v <-> Hyperedge[j]` for each `v ∈ e_j`.
- Returns a simple `Graph`, validated and memoized as usual.
- Unevaluated if some vertex of `h` is itself such a `Hyperedge[j]` node.
- Benchmark (experiment 96): 10⁵ hyperedges 38.5 ms warm / 38.3 ms cold
  (Mathematica 313 ms, xgi 466 ms).

```mathematica
In[1]:= InputForm[HypergraphStarExpansion[Hypergraph[{{1,2},{2,3}}]]]
Out[1]= Graph[{1, 2, 3, Hyperedge[1], Hyperedge[2]}, {1 <-> Hyperedge[1], 2 <-> Hyperedge[1], 2 <-> Hyperedge[2], 3 <-> Hyperedge[2]}]

In[2]:= EdgeCount[HypergraphStarExpansion[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]]
Out[2]= 9

In[3]:= InputForm[HypergraphStarExpansion[Hypergraph[{Hyperedge[1], 2},{{Hyperedge[1],2}}]]]
Out[3]= HypergraphStarExpansion[Hypergraph[{Hyperedge[1], 2}, {{Hyperedge[1], 2}}]]
```

## HypergraphToGraph

- `HypergraphToGraph[h]`: the directed `Graph` obtained by reading `h` as an
  ordered hypergraph.
- `HypergraphToGraph[{e1, e2, ...}]`: the same for a bare List of hyperedges.

**Features**:
- Wolfram Function Repository name and semantics: hyperedge `{v1, ..., vk}`
  contributes `v_a -> v_b` for every `a < b`.
- Accepts a bare List of hyperedges, as the FR function does.
- *Deviation:* Mathilda graphs are simple, so self-loops (from a repeated
  vertex) are dropped and parallel copies merged, and every vertex is kept (the
  FR function returns a multigraph and drops vertices that occur only in unary
  hyperedges).
- Benchmark (experiment 96): 10⁵ hyperedges 59.3 ms (FR function 233 ms, xgi
  376 ms).

```mathematica
In[1]:= InputForm[HypergraphToGraph[{{1,2,3},{3,4}}]]
Out[1]= Graph[{1, 2, 3, 4}, {1 -> 2, 1 -> 3, 2 -> 3, 3 -> 4}]

In[2]:= InputForm[HypergraphToGraph[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]]
Out[2]= Graph[{1, 2, 3, 4, 5, 6, 7}, {1 -> 2, 1 -> 3, 2 -> 3, 3 -> 4, 4 -> 5, 4 -> 6, 5 -> 6}]

In[3]:= InputForm[HypergraphToGraph[{{1,1,2},{1,2},{5}}]]
Out[3]= Graph[{1, 2, 5}, {1 -> 2}]

In[4]:= HypergraphToGraph[5]
Out[4]= HypergraphToGraph[5]
```

## HypergraphLineGraph

- `HypergraphLineGraph[h]`: the `Graph` on `1..m` with `i <-> j` when hyperedges
  `i` and `j` intersect.
- `HypergraphLineGraph[h, s]`: the **s-line graph** — `i <-> j` when hyperedges
  `i` and `j` share at least `s` vertices.

**Features**:
- s-line graphs follow Aksoy, Joslyn, Ortiz Marrero, Praggastis, Purvine,
  "Hypernetwork science via high-order hypergraph walks" (2020). Hyperedges
  are read as sets.
- Edges ordered by `i`, then `j`. Every hyperedge index is a vertex, even when
  isolated.
- `O(Σ_v deg(v)²)` — the number of hyperedge pairs meeting at a vertex, which
  bounds the 1-line graph's size anyway.
- `s` must be a positive integer; otherwise unevaluated.
- Benchmark (experiment 96): 10⁵ hyperedges 47.0 ms warm / 48.5 ms cold
  (Mathematica 819 ms, xgi 1117 ms).

```mathematica
In[1]:= EdgeList[HypergraphLineGraph[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]]
Out[1]= {1 <-> 2, 2 <-> 3}

In[2]:= EdgeList[HypergraphLineGraph[Hypergraph[{{1,2,3},{2,3,4},{3,4,5},{1,5}}], 2]]
Out[2]= {1 <-> 2, 2 <-> 3}

In[3]:= InputForm[HypergraphLineGraph[Hypergraph[{{1,2,3},{2,3,4},{3,4,5},{1,5}}], 3]]
Out[3]= Graph[{1, 2, 3, 4}, {}]

In[4]:= HypergraphLineGraph[Hypergraph[{{1,2,3},{2,3,4},{3,4,5},{1,5}}], 0]
Out[4]= HypergraphLineGraph[Hypergraph[<5 vertices, 4 hyperedges>], 0]
```

## HypergraphConnectedComponents / ConnectedHypergraphQ

- `HypergraphConnectedComponents[h]`: the vertex components of `h` (two vertices
  are joined by a chain of hyperedges), as vertex Lists.
- `ConnectedHypergraphQ[h]`: `True` iff `h` has at least one vertex and exactly
  one component.

**Features**:
- Components are ordered by their first vertex (Mathilda's
  `ConnectedComponents` convention), vertices within a component in VertexList
  order. Isolated vertices are singleton components.
- Union–find, near-linear in the total incidence.
- `ConnectedHypergraphQ` is the Wolfram Function Repository name. It accepts a
  bare List of hyperedges; the FR function leaves `{}` unevaluated, Mathilda
  gives `False`. It gives `False` for any non-hypergraph.
- Benchmark (experiment 96): `HypergraphConnectedComponents` on 10⁵ hyperedges
  2.3 ms warm and cold (Mathematica 160 ms, xgi 282 ms);
  `ConnectedHypergraphQ` 20.1 ms (FR function 261 ms, xgi 546 ms).

```mathematica
In[1]:= HypergraphConnectedComponents[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]
Out[1]= {{1, 2, 3, 4, 5, 6}, {7}}

In[2]:= HypergraphConnectedComponents[Hypergraph[{c,b,a},{{a,c}}]]
Out[2]= {{c, a}, {b}}

In[3]:= ConnectedHypergraphQ[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]
Out[3]= False

In[4]:= ConnectedHypergraphQ[{{1,2},{2,3}}]
Out[4]= True

In[5]:= ConnectedHypergraphQ[{}]
Out[5]= False
```

## HyperedgeConnectedComponents

- `HyperedgeConnectedComponents[h]`: the connected components of the
  hyperedges of `h`, as Lists of 1-based hyperedge indices.
- `HyperedgeConnectedComponents[h, s]`: the **s-connected components** —
  hyperedges with at least `s` distinct vertices, joined when they share at
  least `s` vertices, closed transitively.

**Features**:
- s-connectivity follows Aksoy et al. (2020). The default is `s = 1`.
- Hyperedges with fewer than `s` vertices lie on no s-walk and are omitted (so
  empty hyperedges never appear).
- `s = 1` is union–find over the incidence lists; `s ≥ 2` counts overlaps per
  hyperedge, `O(Σ_v deg(v)²)`.
- `s` must be a positive integer; otherwise unevaluated.
- Benchmark (experiment 96): `s = 2` on 10⁵ hyperedges 18.4 ms warm / 18.6 ms
  cold (Mathematica 1306 ms, xgi 2466 ms).

```mathematica
In[1]:= HyperedgeConnectedComponents[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]
Out[1]= {{1, 2, 3}, {4}}

In[2]:= HyperedgeConnectedComponents[Hypergraph[{{1,2,3},{2,3,4},{3,4,5},{1,5}}], 2]
Out[2]= {{1, 2, 3}, {4}}

In[3]:= HyperedgeConnectedComponents[Hypergraph[{{1,2,3},{2,3,4},{3,4,5},{1,5}}], 3]
Out[3]= {{1}, {2}, {3}}

In[4]:= HyperedgeConnectedComponents[Hypergraph[{1,2},{{1,2},{},{1}}]]
Out[4]= {{1, 3}}

In[5]:= HyperedgeConnectedComponents[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], 0]
Out[5]= HyperedgeConnectedComponents[Hypergraph[<7 vertices, 4 hyperedges>], 0]
```

## HyperedgeDistance / HypergraphDistance

- `HyperedgeDistance[h, i, j]`: the length of a shortest walk from hyperedge `i`
  to hyperedge `j`, consecutive hyperedges intersecting.
- `HyperedgeDistance[h, i, j, s]`: the length of a shortest s-walk
  (consecutive hyperedges sharing at least `s` vertices).
- `HyperedgeDistance[h, i]` / `HyperedgeDistance[h, i, All, s]`: distances from
  hyperedge `i` to every hyperedge.
- `HypergraphDistance[h, u, v]`: the least number of hyperedges in a chain of
  vertices from `u` to `v`.
- `HypergraphDistance[h, u]`: distances from `u` to every vertex, in
  VertexList order.

**Features**:
- `Infinity` when no walk exists; `0` for `i == j` (resp. `u == v`).
- `HyperedgeDistance` runs a BFS directly on the incidence structure — the line
  graph is never materialised. With `s = 1` it is linear: each vertex's
  incidence list is expanded at most once per BFS; with `s ≥ 2` it is
  `O(Σ_v deg(v)²)`. The all-targets form is packed when all distances are
  finite.
- `HypergraphDistance` is the distance in the clique expansion
  (`HypergraphCliqueExpansion`).
- An out-of-range hyperedge index or an unknown vertex leaves the call
  unevaluated.
- Benchmark (experiment 96), single source on 10⁵ hyperedges:
  `HyperedgeDistance` 2.8 ms warm / 3.3 ms cold (Mathematica 7255 ms, xgi
  1048 ms); `HypergraphDistance` 2.9 ms warm / 3.5 ms cold (Mathematica
  5791 ms, xgi 1057 ms).

```mathematica
In[1]:= HyperedgeDistance[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], 1, 3]
Out[1]= 2

In[2]:= HyperedgeDistance[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], 1]
Out[2]= {0, 1, 2, Infinity}

In[3]:= HyperedgeDistance[Hypergraph[{{1,2,3},{2,3,4},{3,4,5},{1,5}}], 1, All, 2]
Out[3]= {0, 1, 2, Infinity}

In[4]:= HypergraphDistance[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], 1, 6]
Out[4]= 3

In[5]:= HypergraphDistance[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], 1]
Out[5]= {0, 1, 1, 2, 3, 3, Infinity}

In[6]:= HyperedgeDistance[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], 1, 9]
Out[6]= HyperedgeDistance[Hypergraph[<7 vertices, 4 hyperedges>], 1, 9]
```

## HypergraphVertexAdd / HypergraphVertexDelete

- `HypergraphVertexAdd[h, v]` / `HypergraphVertexAdd[h, {v1, ...}]`: appends the
  vertices not already present.
- `HypergraphVertexDelete[h, v]` / `HypergraphVertexDelete[h, {v1, ...}]`:
  removes the vertices **and every hyperedge containing one of them**, as
  `VertexDelete` does for a Graph.

**Features**:
- In the Vertex heads any List is a list of vertices; wrap a List-valued vertex
  as `{{...}}`.
- `HypergraphVertexDelete` is unevaluated if a named vertex is absent.
- Mutators unshare a memoized hypergraph node before editing, so the original
  object is never changed.

```mathematica
In[1]:= InputForm[HypergraphVertexAdd[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], {1, 8, 9}]]
Out[1]= Hypergraph[{1, 2, 3, 4, 5, 6, 7, 8, 9}, {{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]

In[2]:= InputForm[HypergraphVertexAdd[Hypergraph[{{1,2}}], {{a,b}}]]
Out[2]= Hypergraph[{1, 2, {a, b}}, {{1, 2}}]

In[3]:= InputForm[HypergraphVertexDelete[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], 3]]
Out[3]= Hypergraph[{1, 2, 4, 5, 6, 7}, {{4, 5, 6}, {7}}]

In[4]:= InputForm[HypergraphVertexDelete[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], {1, 7}]]
Out[4]= Hypergraph[{2, 3, 4, 5, 6}, {{3, 4}, {4, 5, 6}}]

In[5]:= HypergraphVertexDelete[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], 99]
Out[5]= HypergraphVertexDelete[Hypergraph[<7 vertices, 4 hyperedges>], 99]
```

## HypergraphEdgeAdd / HypergraphEdgeDelete

- `HypergraphEdgeAdd[h, e]` / `HypergraphEdgeAdd[h, {e1, ...}]`: appends
  hyperedges, adding new vertices in first-appearance order.
- `HypergraphEdgeDelete[h, e]` / `HypergraphEdgeDelete[h, {e1, ...}]`: removes
  every hyperedge `SameQ` to a named one.

**Features**:
- In the Edge heads, a List whose every element is a List is a list of
  hyperedges; otherwise it is one hyperedge.
- `HypergraphEdgeAdd` allows repeats (a multi-hypergraph).
- `HypergraphEdgeDelete` compares hyperedges as written (`{2, 1}` does not
  delete `{1, 2}`) and removes all copies of a repeated hyperedge. It is
  unevaluated if a named hyperedge does not occur.

```mathematica
In[1]:= InputForm[HypergraphEdgeAdd[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], {7, 8}]]
Out[1]= Hypergraph[{1, 2, 3, 4, 5, 6, 7, 8}, {{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}, {7, 8}}]

In[2]:= InputForm[HypergraphEdgeAdd[Hypergraph[{{1,2}}], {1,2}]]
Out[2]= Hypergraph[{1, 2}, {{1, 2}, {1, 2}}]

In[3]:= InputForm[HypergraphEdgeDelete[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], {{3,4},{7}}]]
Out[3]= Hypergraph[{1, 2, 3, 4, 5, 6, 7}, {{1, 2, 3}, {4, 5, 6}}]

In[4]:= InputForm[HypergraphEdgeDelete[Hypergraph[{{1,2},{2,3},{1,2}}], {1,2}]]
Out[4]= Hypergraph[{1, 2, 3}, {{2, 3}}]

In[5]:= HypergraphEdgeDelete[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], {4, 3}]
Out[5]= HypergraphEdgeDelete[Hypergraph[<7 vertices, 4 hyperedges>], {4, 3}]
```

## Subhypergraph / HypergraphRestriction

- `Subhypergraph[h, {v1, ...}]`: the sub-hypergraph on the named vertices, keeping
  only the hyperedges lying entirely among them.
- `HypergraphRestriction[h, {v1, ...}]`: Berge's induced sub-hypergraph — every
  hyperedge intersected with the vertex set.

**Features**:
- `Subhypergraph` follows `Subgraph`'s rule: the named vertices that are in `h`
  (in VertexList order; absent names are ignored) and the hyperedges lying
  entirely among them.
- `HypergraphRestriction` preserves the order and repeats of the kept vertices
  within each hyperedge, and drops hyperedges that miss the vertex set.
- A bare List of hyperedges is not accepted (unevaluated).

```mathematica
In[1]:= InputForm[Subhypergraph[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], {1, 2, 3, 4}]]
Out[1]= Hypergraph[{1, 2, 3, 4}, {{1, 2, 3}, {3, 4}}]

In[2]:= InputForm[Subhypergraph[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], {4, 3, 99}]]
Out[2]= Hypergraph[{3, 4}, {{3, 4}}]

In[3]:= InputForm[HypergraphRestriction[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], {1, 2, 3, 4}]]
Out[3]= Hypergraph[{1, 2, 3, 4}, {{1, 2, 3}, {3, 4}, {4}}]

In[4]:= InputForm[HypergraphRestriction[Hypergraph[{{3,1,3,2},{5}}], {3, 2}]]
Out[4]= Hypergraph[{3, 2}, {{3, 3, 2}}]

In[5]:= Subhypergraph[{{1,2}}, {1}]
Out[5]= Subhypergraph[{{1, 2}}, {1}]
```

## RandomHypergraph

- `RandomHypergraph[{n, m}, k]`: a random k-uniform hypergraph on `1..n` with `m`
  hyperedges.
- `RandomHypergraph[{n, m}, k, c]`: a list of `c` such hypergraphs.
- `RandomHypergraph[{n, {e, a}}]`: the Wolfram Function Repository form — `e`
  hyperedges of arity `a`.
- `RandomHypergraph[{n, {{e1, a1}, ...}}]`: `e_i` hyperedges of arity `a_i`.

**Features**:
- k-uniform form: `m` independent, uniformly random k-subsets of `1..n`
  (listed increasing); a hyperedge may repeat. The vertex list is all of
  `1..n`. Requires `0 ≤ k ≤ n`, else unevaluated. Floyd's sampling, `O(k)` per
  hyperedge.
- FR form: vertices drawn uniformly **with replacement** from `1..n` (so a
  vertex can repeat inside a hyperedge, as in the FR output `{8, 1, 1}`); the
  vertex list is the vertices that occur, in first-appearance order.
  *Deviation:* returns a `Hypergraph` where the FR function returns the bare
  edge List (`EdgeList` recovers it).
- Both forms draw from the user-visible random stream, so `SeedRandom`
  reproduces them.
- Benchmark (experiment 96), FR form with 10⁵ triples: 10.6 ms (FR function
  1.6 ms, xgi 438 ms) — the one loss in the experiment. The ResourceFunction
  returns a single packed `RandomInteger` array, while Mathilda returns a
  validated `Hypergraph` whose hyperedges are 10⁵ boxed `List` nodes
  (Mathilda's representation invariant keeps packed buffers out of expression
  trees); allocating those nodes alone exceeds 1.6 ms. For scale, Mathilda's
  own `RandomInteger[{1, 20000}, {100000, 3}]` takes 16 ms.

```mathematica
In[1]:= SeedRandom[1]; EdgeList[RandomHypergraph[{10, 4}, 3]]
Out[1]= {{2, 7, 9}, {2, 6, 10}, {1, 5, 8}, {2, 4, 9}}

In[2]:= SeedRandom[1]; RandomHypergraph[{10, 4}, 3]
Out[2]= Hypergraph[<10 vertices, 4 hyperedges>]

In[3]:= SeedRandom[3]; Map[EdgeList, RandomHypergraph[{5, 2}, 2, 3]]
Out[3]= {{{1, 4}, {4, 5}}, {{1, 3}, {1, 4}}, {{1, 3}, {1, 3}}}

In[4]:= SeedRandom[1]; InputForm[RandomHypergraph[{10, {3, 3}}]]
Out[4]= Hypergraph[{9, 8, 2, 6, 10, 1}, {{9, 8, 2}, {8, 2, 6}, {10, 6, 1}}]

In[5]:= SeedRandom[2]; EdgeList[RandomHypergraph[{6, {{2, 3}, {1, 1}}}]]
Out[5]= {{5, 4, 4}, {2, 3, 5}, {3}}

In[6]:= RandomHypergraph[{3, 2}, 4]
Out[6]= RandomHypergraph[{3, 2}, 4]
```

## TransversalHypergraph / FindMinimumTransversal

- `TransversalHypergraph[h]`: every **minimal** transversal of `h`, as
  `Hypergraph[VertexList[h], Tr]`.
- `TransversalHypergraph[{e1, e2, ...}]`: the List `Tr` of minimal transversals
  of a bare List of hyperedges.
- `FindMinimumTransversal[h]` / `FindMinimumTransversal[{e1, ...}]`: one
  minimum-cardinality transversal.

**Features**:
- A transversal (hitting set) meets every hyperedge; it is minimal when no
  proper subset is one.
- `TransversalHypergraph` is the Wolfram Function Repository name. In the
  bare-List form, non-List elements (the FR's "isolated vertices") are
  ignored. No hyperedges → `{{}}`; an empty hyperedge → `{}`. The FR function
  leaves `{}`, `{{}}`, `{{1},{1,2}}` and `{{1,2},{1,2}}` unevaluated; these are
  answered.
- Each transversal lists its vertices in VertexList order; transversals are
  sorted by size, then lexicographically by position.
- `TransversalHypergraph` algorithm: **MMCS** (Murakami & Uno, 2014), run
  iteratively, with `O(deg v)` critical-hyperedge bookkeeping (per hyperedge,
  `|F ∩ S|` and the sum of `S`'s members in `F`, which *is* the unique member
  when the count is 1).
- `FindMinimumTransversal` (NP-hard): exact branch and bound — branch on the
  uncovered hyperedge with fewest open vertices (include, then exclude for
  later siblings); prune with the best of a degree bound, a low-degree-first
  disjoint-hyperedge packing, and a fractional-packing (LP-dual) bound with one
  dual-ascent pass; seeded by the greedy max-coverage bound (lazy heap).
  Unevaluated if a hyperedge is empty.
- **Budgets (correct-or-unevaluated).** Both searches count nodes and return
  unevaluated — never a partial or unproven answer — past 5·10⁷ MMCS nodes or
  2·10⁷ output entries (`TransversalHypergraph`), or 2·10⁷ branch-and-bound
  nodes (`FindMinimumTransversal`, which also refuses a greedy bound above
  20000 to cap C-stack depth). Node counts, not wall clock, keep answers
  machine-independent; both poll `tc_check_deadline()` every 4096 nodes, so
  `TimeConstrained` interrupts them.
- Benchmark (experiment 96): `TransversalHypergraph` on 16 vertices 4.2 ms (FR
  function 5219 ms); `FindMinimumTransversal` on 80 vertices × 200 hyperedges
  245 ms against MILP baselines (Mathematica `LinearOptimization` 874 ms,
  scipy/HiGHS 784 ms) — 3.2–3.6× faster at this size, and comparable to the
  MILPs at 100 vertices.

```mathematica
In[1]:= TransversalHypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]
Out[1]= {{1, 4, 7}, {2, 4, 7}, {3, 4, 7}, {3, 5, 7}, {3, 6, 7}}

In[2]:= InputForm[TransversalHypergraph[Hypergraph[{{a,b},{b,c}}]]]
Out[2]= Hypergraph[{a, b, c}, {{b}, {a, c}}]

In[3]:= {TransversalHypergraph[{}], TransversalHypergraph[{{}}], TransversalHypergraph[{{1},{1,2}}], TransversalHypergraph[{{1,2},{1,2}}]}
Out[3]= {{{}}, {}, {{1}}, {{1}, {2}}}

In[4]:= FindMinimumTransversal[{{1,2,3},{3,4},{4,5,6},{7}}]
Out[4]= {3, 4, 7}

In[5]:= FindMinimumTransversal[Hypergraph[{{a,b},{b,c},{c,d}}]]
Out[5]= {b, c}

In[6]:= FindMinimumTransversal[{{1,2},{}}]
Out[6]= FindMinimumTransversal[{{1, 2}, {}}]
```
