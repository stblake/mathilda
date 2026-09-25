# HypergraphVertexDelete

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphVertexDelete[h, v] or [h, {v1, ...}] deletes the vertices and every hyperedge containing one of them, as VertexDelete does for a Graph.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

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

## Implementation notes

- In the Vertex heads any List is a list of vertices; wrap a List-valued vertex
  as `{{...}}`.
- `HypergraphVertexDelete` is unevaluated if a named vertex is absent.
- Mutators unshare a memoized hypergraph node before editing, so the original
  object is never changed.

**Attributes:** `Protected`.

## References

**See also:** [HypergraphVertexAdd](../../hypergraphs/HypergraphVertexAdd/), [VertexDelete](../../graphs/VertexDelete/)

- Source: [`src/graph/hyp_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_init.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)
