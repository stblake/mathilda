# VertexReplace

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexReplace[g, {v1 -> w1, ...}] replaces vertices of g according to the rules (applied to each vertex as by Replace, so patterns work). Vertices mapped together merge; a merge that would create a self-loop or a parallel edge leaves the call unevaluated. Edge weights are kept.`**

## Examples (9)

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

### Applications (5)

Rename one vertex

```mathematica
In[5]:= VertexList[VertexReplace[PathGraph[{1, 2, 3}], 2 -> x]]
Out[5]= {1, x, 3}
```

Edges follow the renamed vertices

```mathematica
In[6]:= EdgeList[VertexReplace[CycleGraph[4], {1 -> a, 2 -> b}]]
Out[6]= {a <-> b, b <-> 3, 3 <-> 4, 4 <-> a}
```

Patterns and delayed rules are allowed

```mathematica
In[7]:= VertexList[VertexReplace[PathGraph[{1, 2, 3}], n_Integer :> n^2]]
Out[7]= {1, 4, 9}
```

Directions are preserved

```mathematica
In[8]:= EdgeList[VertexReplace[Graph[{1 -> 2, 2 -> 3}], {1 -> a, 3 -> c}]]
Out[8]= {a -> 2, 2 -> c}
```

A swap permutes the labels

```mathematica
In[9]:= EdgeList[VertexReplace[Graph[{1 -> 2, 3 -> 2}], {1 -> 3, 3 -> 1}]]
Out[9]= {3 -> 2, 1 -> 2}
```

## Implementation notes

**Algorithm.** `builtin_vertex_replace` takes `VertexReplace[g, rule]` or `VertexReplace[g, {rules}]`, where every rule must be a `Rule` or `RuleDelayed`. It evaluates `Replace[VertexList[g], rules, {1}]` once, so first-match-wins, patterns and delayed right-hand sides behave exactly as in `Replace`. The images become the new vertex list, de-duplicated through a scratch hash with a position map from old to new index. Each edge is rebuilt with its original head and orientation from the mapped endpoints (untouched edges are shared), and weights are copied.

**Data structures.** `Graph[List, List]` expression tree; the old-to-new vertex map is an `int` array, and the rebuilt endpoint arrays go to `gops_graph_new`.

**Complexity / limits.** `O(V + E)` plus one `Replace` evaluation over the vertex list. If the replacement merges two vertices so that a self-loop or parallel edge would result (for example renaming 2 to 1 on a path 1-2-3), the graph is not simple and the call stays unevaluated. Non-rule second arguments are also left alone.

- `Protected`. A non-graph first argument is left unevaluated.
- Patterns and `RuleDelayed` work, as in `Replace`.
- Vertices mapped to the same name merge.
- Deviation: Mathilda graphs are simple, so merging two adjacent vertices (which
  would create a self-loop) is left unevaluated; Mathematica returns a
  multigraph.

**Attributes:** `Protected`.

## References

**See also:** [Replace](../../assignment-and-rules/Replace/), [RuleDelayed](../../assignment-and-rules/RuleDelayed/)

- Source: [`src/graph/gops_edit.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_edit.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

The rules are applied to the vertex list with the semantics of `Replace` at level 1, so the first matching rule wins and vertices with no match keep their name.

If the renaming sends two vertices to one name and that would create a self-loop or parallel edge, the result would not be a simple graph and the call is left unevaluated.
