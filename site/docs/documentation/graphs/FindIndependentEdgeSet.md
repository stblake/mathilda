# FindIndependentEdgeSet

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindIndependentEdgeSet[g] gives a maximum independent edge set (maximum matching) of g, ignoring edge direction. Hopcroft-Karp for bipartite graphs, Edmonds' blossom algorithm otherwise.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= FindIndependentEdgeSet[PetersenGraph[]]
Out[1]= {1 <-> 3, 2 <-> 4, 5 <-> 10, 6 <-> 7, 8 <-> 9}

In[2]:= FindIndependentEdgeSet[StarGraph[4]]
Out[2]= {1 <-> 2}

In[3]:= FindIndependentEdgeSet[Graph[{1->2,3->2,3->4}]]
Out[3]= {1 -> 2, 3 -> 4}

In[4]:= FindIndependentEdgeSet[x]
Out[4]= FindIndependentEdgeSet[x]
```

## Implementation notes

- `Protected`; unevaluated on a non-graph.
- Edge direction is ignored; edges are returned in `EdgeList` order.
- Algorithm: Karp-Sipser greedy start, then Hopcroft-Karp when the graph is
  bipartite and Edmonds' blossom algorithm otherwise (union-find blossom bases;
  failed searches retire their Hungarian trees). Matching on a triangulated grid:
  about 2.8 ms against 190 ms in Mathematica.
- Weighted graphs: optimizes cardinality, as the Wolfram documentation states.
  (The differential test found Mathematica returning non-maximum matchings on
  weighted graphs.)

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/galg_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
