# VertexReplace

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexReplace[g, {v1 -> w1, ...}] replaces vertices of g according to the rules (applied to each vertex as by Replace, so patterns work). Vertices mapped together merge; a merge that would create a self-loop or a parallel edge leaves the call unevaluated. Edge weights are kept.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= EdgeList[VertexReplace[CycleGraph[3], {1->a, 2->b}]]
Out[1]= {a <-> b, b <-> 3, 3 <-> a}

In[2]:= VertexList[VertexReplace[PathGraph[Range[4]], x_?EvenQ :> x^2]]
Out[2]= {1, 4, 3, 16}

In[3]:= EdgeList[VertexReplace[Graph[{1<->2, 3<->4}], {3->1}]]
Out[3]= {1 <-> 2, 1 <-> 4}

In[4]:= VertexReplace[PathGraph[Range[3]], 2->1]
Out[4]= VertexReplace[Graph[<3 vertices, 2 edges>], 2 -> 1]
```

## Implementation notes

- `Protected`. A non-graph first argument is left unevaluated.
- Patterns and `RuleDelayed` work, as in `Replace`.
- Vertices mapped to the same name merge.
- Deviation: Mathilda graphs are simple, so merging two adjacent vertices (which
  would create a self-loop) is left unevaluated; Mathematica returns a
  multigraph.

**Attributes:** `Protected`.

## References

**See also:** [Replace](../../assignment-and-rules/Replace/), [RuleDelayed](../../assignment-and-rules/RuleDelayed/)

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
