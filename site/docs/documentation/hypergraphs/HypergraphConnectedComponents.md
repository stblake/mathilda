# HypergraphConnectedComponents

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphConnectedComponents[h] gives the connected components of h as Lists of vertices; two vertices are connected when a chain of hyperedges joins them. Components are ordered by their first vertex.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

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

### Applications (3)

```mathematica
In[6]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
Out[6]= Hypergraph[<7 vertices, 4 hyperedges>]
```

Two components: {1..6} and the isolated 7

```mathematica
In[7]:= Length[HypergraphConnectedComponents[h]]
Out[7]= 2
```

Ordered by first vertex

```mathematica
In[8]:= HypergraphConnectedComponents[Hypergraph[{a, b, c, d}, {{a, b}, {c, d}}]]
Out[8]= {{a, b}, {c, d}}
```

## Implementation notes

**Algorithm.** `builtin_hypergraph_connected_components` gives the vertex
components — two vertices are connected when a chain of hyperedges joins them.
`vertex_uf` runs union–find over the vertices: within each hyperedge's distinct
set it unions the first member with every other. `groups_to_list` then collects
the classes, ordered by their smallest member (first vertex in `VertexList`
order), vertices within a class ascending. An isolated vertex is a singleton
component. Accepts a `Hypergraph` or a bare List of hyperedges (`hyp_arg`).

**Data structures.** The distinct-vertex CSR `soff/sv`; a union–find parent array
with path-halving and union-by-smaller-root; `groups_to_list` scratch (`id`,
`size`, `cid`, per-class member buffers). The result is a `List` of vertex Lists.

**Complexity / limits.** Near-linear, `O(Σ|e| · α(n))`. Note the ordering differs
from Graph `ConnectedComponents` (which is largest-first / strong components); the
hypergraph form is deterministically ordered by first vertex.

- Components are ordered by their first vertex, vertices within a component in
  VertexList order. (Graph `ConnectedComponents` instead follows Mathematica:
  largest first when undirected, strong components when directed.) Isolated
  vertices are singleton components.
- Union–find, near-linear in the total incidence.
- `ConnectedHypergraphQ` is the Wolfram Function Repository name. It accepts a
  bare List of hyperedges; the FR function leaves `{}` unevaluated, Mathilda
  gives `False`. It gives `False` for any non-hypergraph.
- Benchmark (experiment 96): `HypergraphConnectedComponents` on 10⁵ hyperedges
  2.3 ms warm and cold (Mathematica 160 ms, xgi 282 ms);
  `ConnectedHypergraphQ` 20.1 ms (FR function 261 ms, xgi 546 ms).

**Attributes:** `Protected`.

## References

**See also:** [ConnectedHypergraphQ](../../hypergraphs/ConnectedHypergraphQ/), [ConnectedComponents](../../graphs/ConnectedComponents/)

- Source: [`src/graph/hyp_ops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_ops.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

Two vertices are in the same component when a chain of hyperedges joins them;
`HypergraphConnectedComponents` returns these classes as vertex Lists. Isolated
vertices are singleton components.

The ordering is deterministic and differs from Graph `ConnectedComponents`
(largest-first when undirected): hypergraph components are ordered by their first
vertex in `VertexList` order, with vertices inside a component ascending. The
engine is union–find, near-linear in the total incidence. `ConnectedHypergraphQ`
is the companion predicate — one non-empty component.
