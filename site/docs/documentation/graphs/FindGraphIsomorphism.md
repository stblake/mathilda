# FindGraphIsomorphism

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindGraphIsomorphism[g1, g2] gives {assoc}, an isomorphism from g1 to g2 as an association of vertices, or {} if the graphs are not isomorphic. FindGraphIsomorphism[g1, g2, n] / [g1, g2, All] gives up to n / all isomorphisms.`**

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= FindGraphIsomorphism[CycleGraph[4], Graph[{a,b,c,d},{UndirectedEdge[a,c],UndirectedEdge[c,b],UndirectedEdge[b,d],UndirectedEdge[d,a]}]]
Out[1]= {<|1 -> a, 2 -> c, 3 -> b, 4 -> d|>}

In[2]:= FindGraphIsomorphism[CycleGraph[4], CycleGraph[4], 2]
Out[2]= {<|1 -> 1, 2 -> 2, 3 -> 3, 4 -> 4|>, <|1 -> 1, 2 -> 4, 3 -> 3, 4 -> 2|>}

In[3]:= Length[FindGraphIsomorphism[CycleGraph[4], CycleGraph[4], All]]
Out[3]= 8

In[4]:= FindGraphIsomorphism[CycleGraph[4], StarGraph[4]]
Out[4]= {}

In[5]:= FindGraphIsomorphism[Graph[{},{}], Graph[{},{}]]
Out[5]= {<||>}
```

### Applications (5)

One map, here the identity

```mathematica
In[6]:= FindGraphIsomorphism[CycleGraph[4], CycleGraph[4]]
Out[6]= {<|1 -> 1, 2 -> 2, 3 -> 3, 4 -> 4|>}
```

Relabels the vertices of the first graph

```mathematica
In[7]:= FindGraphIsomorphism[PathGraph[{1, 2, 3}], PathGraph[{b, a, c}]]
Out[7]= {<|1 -> b, 2 -> a, 3 -> c|>}
```

At most two maps

```mathematica
In[8]:= FindGraphIsomorphism[PathGraph[{1, 2, 3}], PathGraph[{a, b, c}], 2]
Out[8]= {<|1 -> a, 2 -> b, 3 -> c|>, <|1 -> c, 2 -> b, 3 -> a|>}
```

The dihedral group of the square has eight elements

```mathematica
In[9]:= Length[FindGraphIsomorphism[CycleGraph[4], CycleGraph[4], All]]
Out[9]= 8
```

Non-isomorphic graphs give the empty list

```mathematica
In[10]:= FindGraphIsomorphism[CycleGraph[4], PathGraph[{1, 2, 3, 4}]]
Out[10]= {}
```

## Implementation notes

**Algorithm.** `builtin_find_graph_isomorphism` accepts `FindGraphIsomorphism[g, h]` (one map) or `[g, h, n]` / `[g, h, All]` (up to `n` maps). Each graph is reduced by `gi_build` to a vertex-coloured structure with up to three relations (undirected, directed out, directed in): self-loops are folded into the vertex colour, and an edge class of multiplicity `k > 1` is subdivided by a new vertex coloured `(k, kind)`, so mixed graphs, loops and multigraphs are all handled; weights are ignored. `gi_compatible` rejects on vertex/edge counts, degree histograms and colour multisets and returns `{}`. Otherwise the engine `galg_iso_find` (one map) or `galg_iso_enumerate` (several) in `galg_iso.c` runs individualization-refinement: both graphs are refined to the coarsest equitable partition (1-dimensional Weisfeiler-Leman) in lockstep, a path is followed in `g` and a trace-identical path searched in `h`, aborting on the first deviating event and verifying the leaf edge by edge. If the direct search passes about 4n nodes it falls back to comparing canonical forms. The result is a list of `Association`s `vertex of g -> vertex of h`.

**Data structures.** `Graph[List, List]` expression trees are flattened to CSR relations (`GalgIsoGraph`); the engine keeps an ordered partition with exact undo logs, so memory is linear in the graph plus the current path. Enumerated maps are collected in an `int` buffer capped at 2^27 cells.

**Complexity / limits.** Graph isomorphism has no known polynomial algorithm; one refinement is `O((n + m) log n)` and the tree is small for most graphs. The search is budgeted at 5e7 refinement nodes and polls the `TimeConstrained` deadline; on exhaustion, or if enumeration exceeds its cap, the head stays unevaluated rather than guessing. `All` on a highly symmetric graph returns the whole automorphism coset, which can be huge.

- `Protected`; unevaluated on non-graph arguments.
- The Association runs over `VertexList[g1]` in order.
- The order of the list of isomorphisms is the engine's search order;
  Mathematica's differs.
- Two empty graphs give `{<||>}` (Mathematica gives `{}`, which contradicts its
  own `IsomorphicGraphQ` answer `True`).
- Uses the `IsomorphicGraphQ` individualization-refinement engine; every map
  returned is verified edge by edge. Directed and mixed graphs are supported;
  edge weights are ignored.

**Attributes:** `Protected`.

## References

**See also:** [IsomorphicGraphQ](../../graphs/IsomorphicGraphQ/)

- B. D. McKay and A. Piperno, *Practical graph isomorphism, II*, J. Symbolic Comput. **60** (2014) 94-112.
- Source: [`src/graph/galg_isoheads.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_isoheads.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

The result is a list of associations sending each vertex of the first graph to a vertex of the second; `{}` means the graphs are not isomorphic. A third argument `n` or `All` asks for up to that many maps, and `All` returns the whole coset of the automorphism group, which can be very large for symmetric graphs.

Mixed graphs, self-loops and multigraphs are supported; edge weights and other properties are ignored. If the search budget or a `TimeConstrained` limit runs out the call stays unevaluated rather than guessing.
