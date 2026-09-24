# DegreeCentrality

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DegreeCentrality[g] gives the list of vertex degrees of g; DegreeCentrality[g, "In"] and [g, "Out"] count incoming and outgoing edges. In a mixed graph an undirected edge counts as one edge in each direction.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= DegreeCentrality[StarGraph[5]]
Out[1]= {4, 1, 1, 1, 1}

In[2]:= DegreeCentrality[Graph[{1->2, 2->3, 3->1, 3->4}], "In"]
Out[2]= {1, 1, 1, 1}

In[3]:= DegreeCentrality[Graph[{1->2, 2->3, 3->1, 3->4}], "Out"]
Out[3]= {1, 1, 2, 0}

In[4]:= DegreeCentrality[Graph[{1<->2, 2->3}]]
Out[4]= {2, 3, 1}

In[5]:= DegreeCentrality[Graph[{1->2, 2->3, 3->1, 3->4}], "Foo"]
Out[5]= DegreeCentrality[Graph[<4 vertices, 4 edges>], "Foo"]
```

## Implementation notes

- Exact degrees, in `VertexList` order; weights ignored.
- On a mixed graph an undirected edge counts once in each direction (so it
  adds 2 to the default total), as in Wolfram.
- An unknown direction string leaves the call unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [VertexList](../../graphs/VertexList/)

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
