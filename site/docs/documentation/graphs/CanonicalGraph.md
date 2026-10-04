# CanonicalGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CanonicalGraph[g] gives a canonical form of g on vertices 1..n: two graphs are isomorphic iff their canonical graphs are identical.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= EdgeList[CanonicalGraph[Graph[{UndirectedEdge[3,1],UndirectedEdge[1,2]}]]]
Out[1]= {1 <-> 3, 2 <-> 3}

In[2]:= CanonicalGraph[PathGraph[{1,2,3}]] === CanonicalGraph[Graph[{UndirectedEdge[2,1],UndirectedEdge[3,1]}]]
Out[2]= True

In[3]:= CanonicalGraph[CycleGraph[4]] === CanonicalGraph[StarGraph[4]]
Out[3]= False

In[4]:= CanonicalGraph[x]
Out[4]= CanonicalGraph[x]
```

### Options (1)

```mathematica
In[5]:= CanonicalGraph[Graph[{1,2,3},{UndirectedEdge[1,2]}, EdgeWeight->{7}]] === CanonicalGraph[Graph[{1,2,3},{UndirectedEdge[1,3]}]]
Out[5]= True
```

### Applications (3)

Relabelled onto 1..n, edges sorted

```mathematica
In[6]:= EdgeList[CanonicalGraph[CompleteGraph[3]]]
Out[6]= {1 <-> 2, 1 <-> 3, 2 <-> 3}
```

A triangle on {2,3,4} canonicalises to one on {1,2,3}

```mathematica
In[7]:= EdgeList[CanonicalGraph[Graph[{2 <-> 3, 3 <-> 4, 4 <-> 2}]]]
Out[7]= {1 <-> 2, 1 <-> 3, 2 <-> 3}
```

Isomorphic graphs share a canonical form

```mathematica
In[8]:= CanonicalGraph[Graph[{1 <-> 2, 2 <-> 3}]] === CanonicalGraph[Graph[{5 <-> 9, 9 <-> 7}]]
Out[8]= True
```

## Implementation notes

**Algorithm.** `builtin_canonical_graph` returns a relabelling of `g` onto the
vertices `1..n` such that two graphs are isomorphic iff their canonical graphs
are identical. It first reduces `g` to a vertex-coloured simple structure
(`gi_build`): self-loops fold into the vertex colour, and an edge class of
multiplicity `k > 1` is subdivided through a new colour-`(k, kind)` vertex, so
mixed graphs, loops and multigraphs are all handled. `gi_canonical_rank` then
runs `galg_iso_canon` — the individualization-refinement engine of
`galg_iso.c`, the nauty / bliss / Traces family: a 1-dimensional Weisfeiler-Leman
colour refinement (Hopcroft / Paige-Tarjan splitter queue) to the coarsest
equitable partition, a search tree that individualizes one vertex of the target
cell per level, and a leaf key `(trace₁, …, relabelled graph)` whose best leaf
defines the canonical labelling. The ranks of the original vertices are read off
the canonical leaf, every edge is rewritten on `1..n` (undirected endpoints
sorted), the edge list is sorted, and a fresh `Graph[{1,…,n}, {…}]` is built.

**Data structures.** `GalgIsoGraph` carries up to three CSR relations
(undirected, directed-out, directed-in) plus a ranked `vcol` colour array; the
engine is pure C with no `Expr` dependency. The ordered partition is kept
nauty-style in a `lab[]` position array (a cell named by its start), refinement
is driven by a FIFO of splitter cells, and each refinement emits a 64-bit
**trace** compared on the fly so a losing branch is abandoned early. Edges are
sorted with `qsort` over `(a, b, directed)` triples.

**Complexity / limits.** One refinement is `O((n+m) log n)`; the search tree is
small for almost all graphs but can be super-polynomial on adversarial
constructions (CFI-type), so the whole search is budgeted at `GI_BUDGET = 5·10⁷`
nodes and polls the `TimeConstrained` deadline — on exhaustion the head stays
unevaluated rather than guessing. Edge weights and other properties are ignored,
as in Mathematica.

- `Protected`; unevaluated on a non-graph.
- `CanonicalGraph[g] === CanonicalGraph[h]` iff `g` and `h` are isomorphic.
- The particular representative differs from Mathematica's (both are
  arbitrary).
- Properties such as `EdgeWeight` are dropped, as in Mathematica.
- Canonical labeling keeps the best leaf by (trace, relabeled graph) with
  automorphism pruning (leaf and internal-node automorphisms, jump-back,
  first-path orbits) on the `IsomorphicGraphQ` engine. Directed and mixed graphs
  are supported.

**Attributes:** `Protected`.

## References

**See also:** [EdgeWeight](../../graphs/EdgeWeight/), [IsomorphicGraphQ](../../graphs/IsomorphicGraphQ/)

- B. D. McKay and A. Piperno, *Practical graph isomorphism, II*, J. Symbolic Comput. **60** (2014) 94-112.
- B. Weisfeiler and A. Leman, *The reduction of a graph to canonical form and the algebra which appears therein*, Nauchno-Techn. Inform. Ser. 2, **9** (1968) 12-16.
- Source: [`src/graph/galg_isoheads.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_isoheads.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

`CanonicalGraph[g]` returns a relabelling of `g` onto the vertices `1..n` that is
an isomorphism invariant: two graphs are isomorphic exactly when their canonical
graphs are identical. The last example is the defining property — two differently
labelled paths of length two collapse to the same canonical graph, so `===`
(verbatim equality) is `True`.

The canonical labelling comes from an individualization-refinement engine (the
nauty / Traces family) run on a coloured reduction of `g`, so directed graphs,
self-loops and multigraphs are all supported; edge weights and other properties
are ignored. Since the result is an opaque `Graph` object, the examples read it
back through `EdgeList`. The search is budgeted and leaves the call unevaluated
on exhaustion.
