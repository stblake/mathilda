# EmptyGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EmptyGraphQ[g] gives True if g is a graph with no edges (it may have vertices), and False otherwise.`**

## Examples (6)

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

### Applications (2)

Isolated vertices, no edges

```mathematica
In[5]:= EmptyGraphQ[Graph[{1, 2, 3}, {}]]
Out[5]= True
```

Has edges

```mathematica
In[6]:= EmptyGraphQ[Graph[{1 -> 2, 2 -> 3, 3 -> 1}]]
Out[6]= False
```

## Implementation notes

**Algorithm.** `builtin_empty_graph_q` is a direct structural predicate: it returns `True` iff
its argument is a valid graph whose edge list is empty. The vertex count is irrelevant, so a
graph of isolated vertices counts as empty. No adjacency is built; it reads the length of the
edge `List` (the graph's second argument) after validation.

**Data structures.** None beyond the validated-graph node — the edge-list length comes straight
from the expression, and validation itself is served from the per-node memo.

**Complexity / limits.** `O(1)` after the graph's one-time validation. Exactly one argument. A
non-graph argument returns `False` (it is a predicate), never unevaluated.

- `Protected`. Any number of vertices, including none. `False` for a non-graph
  (see `UndirectedGraphQ`).

**Attributes:** `Protected`.

## References

**See also:** [UndirectedGraphQ](../../graphs/UndirectedGraphQ/)

- Source: [`src/graph/graphprops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graphprops.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)

## Notes & additional examples

### Notes

A graph is empty when it has no edges, whatever its vertices — a set of isolated vertices is an
empty graph. The predicate is about edges only; `VertexCount` can be any value.

Like the other `...Q` predicates it returns `False` for a non-graph argument rather than staying
unevaluated.
