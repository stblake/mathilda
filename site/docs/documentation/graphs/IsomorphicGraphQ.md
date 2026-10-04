# IsomorphicGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`IsomorphicGraphQ[g1, g2, ...] gives True if all the graphs are isomorphic. Colour refinement with individualization-refinement search.`**

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= IsomorphicGraphQ[CycleGraph[5], Graph[{UndirectedEdge[1,3],UndirectedEdge[3,5],UndirectedEdge[5,2],UndirectedEdge[2,4],UndirectedEdge[4,1]}]]
Out[1]= True

In[2]:= IsomorphicGraphQ[CycleGraph[6], PathGraph[Range[6]]]
Out[2]= False

In[3]:= IsomorphicGraphQ[CycleGraph[4], CycleGraph[4], CycleGraph[4]]
Out[3]= True

In[4]:= IsomorphicGraphQ[Graph[{1->2,2->3}], Graph[{1->2,1->3}]]
Out[4]= False

In[5]:= IsomorphicGraphQ[CycleGraph[4], x]
Out[5]= False

In[6]:= IsomorphicGraphQ[CycleGraph[4]]
Out[6]= IsomorphicGraphQ[Graph[<4 vertices, 4 edges>]]
```

### Applications (3)

A relabelled 4-cycle

```mathematica
In[7]:= IsomorphicGraphQ[CycleGraph[4], Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}]]
Out[7]= True
```

Different degree sequences

```mathematica
In[8]:= IsomorphicGraphQ[PathGraph[4], CycleGraph[4]]
Out[8]= False
```

The triangle is both

```mathematica
In[9]:= IsomorphicGraphQ[CompleteGraph[3], CycleGraph[3]]
Out[9]= True
```

## Implementation notes

**Algorithm.** `builtin_isomorphic_graph_q` decides whether all its arguments are pairwise
isomorphic by an individualization-refinement canonical-labelling engine, the family of nauty /
bliss / Traces. The core is colour refinement — one-dimensional Weisfeiler–Leman implemented in
the Hopcroft / Paige–Tarjan style (a FIFO of splitter cells, skipping the largest fragment for
an `O((n+m) log n)` pass). A running 64-bit "trace" hash is emitted per split and compared
on-the-fly against a reference, so a branch aborts the moment it can no longer match. Over that
sits an individualization-refinement search tree with exact per-level undo, automorphism
detection and orbit pruning. The two graphs are compared first by a lockstep `G -> H` trace
search; if that exceeds a node limit it falls back to comparing the two canonical forms. Cheap
rejections (vertex/edge counts, colour-class sizes, degree histograms) run first.

**Data structures.** `GalgIsoGraph` holds up to three CSR relations — undirected, out and in —
plus optional vertex colours, so directed graphs, loops and multigraphs are handled via a
coloured subdivision. The ordered partition carries label/inverse arrays, cell start/size
arrays, and a min segment tree over positions for target-cell selection; automorphisms are kept
as sparse generators with union-find orbits.

**Complexity / limits.** One refinement is `O((n+m) log n)`; the search tree is small for random
or structured graphs, with a worst-case refinement-node budget `GI_BUDGET = 5·10^7` (polled
against `TimeConstrained`). Variadic, at least two arguments; `True` iff all are mutually
isomorphic. Any non-graph argument gives `False`; only budget exhaustion or allocation failure
returns unevaluated.

- `Protected`. `False` if any argument is not a graph; one argument stays
  unevaluated, as in Mathematica.
- Directed and mixed graphs are supported (Mathematica leaves mixed graphs
  unevaluated). Edge weights are ignored, as in Mathematica.
- The reduction to the engine also folds self-loops into vertex colours and
  turns an edge of multiplicity k > 1 into a coloured subdivision vertex, so
  loops and multigraphs will work unchanged once `Graph` accepts them (today's
  validator rejects both). This applies to the whole isomorphism family
  (`FindGraphIsomorphism`, `CanonicalGraph`, `GraphAutomorphismGroup`).
- Engine (`galg_iso.c`): individualization-refinement in the nauty / bliss /
  Traces family. Hopcroft-style equitable refinement (U, out and in relations
  counted separately, largest fragment skipped, every split decision a function
  of cell positions, sizes and counts only) emits a 64-bit trace that is
  compared event by event, so a branch that cannot match dies at its first
  deviating split. The search tree is iterative with exact undo (no recursion at
  depth n). `g -> h` searches run in lockstep against `h`'s trace.
- Exact (see the exactness policy under `FindVertexCover`): complete answer or
  unevaluated on budget exhaustion; polls `TimeConstrained`.
- Performance: on an exactly 3-regular random graph with 10^4 vertices (where
  colour refinement alone learns nothing) `IsomorphicGraphQ` takes about 4 ms
  against about 0.8 s in Mathematica 15; at 10^5 vertices about 0.2 s against
  about 60 s.

**Attributes:** `Protected`.

## References

**See also:** [Graph](../../graphs/Graph/), [FindGraphIsomorphism](../../graphs/FindGraphIsomorphism/), [CanonicalGraph](../../graphs/CanonicalGraph/), [GraphAutomorphismGroup](../../graphs/GraphAutomorphismGroup/), [FindVertexCover](../../graphs/FindVertexCover/), [TimeConstrained](../../time-and-date/TimeConstrained/)

- B. D. McKay and A. Piperno, *Practical graph isomorphism, II*, J. Symbolic Comput. **60** (2014) 94-112.
- R. Paige and R. E. Tarjan, *Three partition refinement algorithms*, SIAM J. Comput. **16** (1987) 973-989.
- Source: [`src/graph/galg_iso.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_iso.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

`IsomorphicGraphQ[g1, g2, ...]` returns `True` when the graphs are all isomorphic — related by a
relabelling of vertices that preserves adjacency. It is variadic: all the arguments must be
mutually isomorphic.

Directed graphs are handled (in- and out-adjacency are matched separately), as are disconnected
graphs. A non-graph argument gives `False`.
