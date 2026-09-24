# CanonicalGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CanonicalGraph[g] gives a canonical form of g on vertices 1..n: two graphs are isomorphic iff their canonical graphs are identical.`**

## Examples (5)

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

## Implementation notes

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

- Source: [`src/graph/galg_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
