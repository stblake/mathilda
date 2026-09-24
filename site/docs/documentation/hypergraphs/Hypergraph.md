# Hypergraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Hypergraph[{e1, e2, ...}] represents a hypergraph whose hyperedges e_i are Lists of vertices; the vertices are derived in first-appearance order. Hypergraph[{v1, ...}, {e1, ...}] gives the vertices explicitly. Hypergraph[g] converts a Graph. Hyperedges may repeat, overlap, nest, be empty, or repeat a vertex; their order is kept for ordered (Wolfram-model) use, while set-based heads read each as the set of its vertices.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

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

### Scope (1)

```mathematica
In[7]:= InputForm[Hypergraph[{1,2},{{1,5}}]]
Out[7]= Hypergraph[{1, 2}, {{1, 5}}]
```

## Options & behaviour

An unknown vertex in a hyperedge leaves the call unevaluated:

## Implementation notes

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

**Attributes:** `Protected`.

## References

**See also:** [Graph](../../graphs/Graph/), [List](../../other-advanced/List/), [HypergraphToGraph](../../hypergraphs/HypergraphToGraph/), [EdgeList](../../graphs/EdgeList/), [InputForm](../../expression-information/InputForm/), [HyperedgeSizes](../../hypergraphs/HyperedgeSizes/), [HypergraphEdgeDelete](../../hypergraphs/HypergraphEdgeDelete/), [FullForm](../../expression-information/FullForm/)

- Source: [`src/graph/hyp_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_init.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)
