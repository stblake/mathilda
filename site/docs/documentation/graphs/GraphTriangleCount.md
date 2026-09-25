# GraphTriangleCount

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphTriangleCount[g] gives the number of triangles of an undirected graph g, or of directed 3-cycles of a directed graph. O(m^1.5) degree-ordered triangle listing.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= GraphTriangleCount[Graph[{1<->2,2<->3,3<->1,3<->4}]]
Out[1]= 1

In[2]:= GraphTriangleCount[CompleteGraph[5]]
Out[2]= 10

In[3]:= GraphTriangleCount[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[3]= 1

In[4]:= GraphTriangleCount[Graph[{1->2,2->3,1->3}]]
Out[4]= 0

In[5]:= GraphTriangleCount[Graph[{1<->2, 2->3}]]
Out[5]= GraphTriangleCount[Graph[<3 vertices, 2 edges>]]
```

## Implementation notes

- Exact, weights ignored, mixed graphs unevaluated (as in Wolfram).
- Undirected: the number of 3-cliques. Directed (*reverse-engineered*):
  triangles are directed 3-cycles, so a transitive triple `1->2, 2->3, 1->3`
  does not count.
- **Algorithm:** triangle listing with an acyclic orientation (degree order,
  or index order when `maxdeg² <= 4m`), O(m^1.5) worst case; directed 3-cycles
  via two-bit arc flags per edge.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
