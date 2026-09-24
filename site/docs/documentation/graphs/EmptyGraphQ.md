# EmptyGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EmptyGraphQ[g] gives True if g is a graph with no edges (it may have vertices), and False otherwise.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= EmptyGraphQ[Graph[{1,2,3},{}]]
Out[1]= True

In[2]:= EmptyGraphQ[Graph[{},{}]]
Out[2]= True

In[3]:= EmptyGraphQ[PathGraph[2]]
Out[3]= False

In[4]:= EmptyGraphQ[5]
Out[4]= False
```

## Implementation notes

- `Protected`. Any number of vertices, including none. `False` for a non-graph
  (see `UndirectedGraphQ`).

**Attributes:** `Protected`.

## References

**See also:** [UndirectedGraphQ](../../graphs/UndirectedGraphQ/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
