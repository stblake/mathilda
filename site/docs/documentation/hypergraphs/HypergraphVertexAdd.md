# HypergraphVertexAdd

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphVertexAdd[h, v] adds the vertex v to h; HypergraphVertexAdd[h, {v1, ...}] adds several. Existing vertices are left alone.`**

## Examples (8)

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

### Applications (3)

A single new vertex

```mathematica
In[6]:= InputForm[HypergraphVertexAdd[Hypergraph[{{1, 2}}], 3]]
Out[6]= Hypergraph[{1, 2, 3}, {{1, 2}}]
```

Several at once

```mathematica
In[7]:= InputForm[HypergraphVertexAdd[Hypergraph[{{1, 2}}], {3, 4}]]
Out[7]= Hypergraph[{1, 2, 3, 4}, {{1, 2}}]
```

An existing vertex is left alone

```mathematica
In[8]:= InputForm[HypergraphVertexAdd[Hypergraph[{{1, 2}}], 2]]
Out[8]= Hypergraph[{1, 2}, {{1, 2}}]
```

## Implementation notes

**Algorithm.** `builtin_hypergraph_vertex_add` appends vertices to `h`, leaving
the hyperedges untouched. `edit_items` interprets the second argument: for the
Vertex heads any `List` is a list of vertices (so a List-valued vertex must be
wrapped as `{{...}}`), otherwise it is a single vertex. Each candidate not already
in the memo's vertex index — and not a duplicate among the candidates, screened by
a throwaway `GraphVIdx extra` — is appended in order after the existing vertices.
A vertex already present is silently skipped. The result is a fresh
`Hypergraph[{...}, edges]`.

**Data structures.** The memo's vertex `GraphVIdx` (membership test) and a
temporary `GraphVIdx extra` (de-duplicating the new candidates); the vertex and
(shared) edge Lists are rebuilt with `mk_hyp`/`mk_list`. Because the memo holds an
immutable reference, the original object is never mutated.

**Complexity / limits.** `O(n + k)` for `k` candidates. Only the vertex set grows;
to add hyperedges use `HypergraphEdgeAdd`.

- In the Vertex heads any List is a list of vertices; wrap a List-valued vertex
  as `{{...}}`.
- `HypergraphVertexDelete` is unevaluated if a named vertex is absent.
- Mutators unshare a memoized hypergraph node before editing, so the original
  object is never changed.

**Attributes:** `Protected`.

## References

**See also:** [HypergraphVertexDelete](../../hypergraphs/HypergraphVertexDelete/), [VertexDelete](../../graphs/VertexDelete/)

- Source: [`src/graph/hyp_ops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_ops.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`HypergraphVertexAdd[h, v]` appends the vertex `v`, and `HypergraphVertexAdd[h, {v1, ...}]`
appends several; the hyperedges are untouched, so this introduces isolated
vertices. A vertex already present is silently skipped, and the new vertices keep
their given order after the existing ones.

In the Vertex heads any `List` second argument is read as a *list of vertices*, so
to add a single List-valued vertex wrap it as `{{...}}`. Because a memoized
hypergraph is held immutably, the edit builds a fresh object and never changes the
original. To add hyperedges instead, use `HypergraphEdgeAdd`.
