# HypergraphDual

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphDual[h] gives the dual hypergraph: vertices 1..m (one per hyperedge of h) and, for each vertex v of h in VertexList order, the hyperedge of indices of the hyperedges containing v.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= InputForm[HypergraphDual[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]]
Out[1]= Hypergraph[{1, 2, 3, 4}, {{1}, {1}, {1, 2}, {2, 3}, {3}, {3}, {4}}]

In[2]:= InputForm[HypergraphDual[HypergraphDual[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]]]
Out[2]= Hypergraph[{1, 2, 3, 4, 5, 6, 7}, {{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]

In[3]:= InputForm[HypergraphDual[Hypergraph[{a,b,c},{{a,b},{a,b}}]]]
Out[3]= Hypergraph[{1, 2}, {{1, 2}, {1, 2}, {}}]

In[4]:= HypergraphDual[{{1, 2}}]
Out[4]= HypergraphDual[{{1, 2}}]
```

### Applications (3)

```mathematica
In[5]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
Out[5]= Hypergraph[<7 vertices, 4 hyperedges>]
```

Vertices 1..m, one hyperedge per vertex of h

```mathematica
In[6]:= InputForm[HypergraphDual[h]]
Out[6]= Hypergraph[{1, 2, 3, 4}, {{1}, {1}, {1, 2}, {2, 3}, {3}, {3}, {4}}]
```

The dual has m vertices and n hyperedges

```mathematica
In[7]:= {VertexCount[HypergraphDual[h]], EdgeCount[HypergraphDual[h]]}
Out[7]= {4, 7}
```

## Implementation notes

**Algorithm.** `builtin_hypergraph_dual` builds the dual: vertices `1..m`, one per
hyperedge of `h`, and one hyperedge per vertex of `h`. It requests the incidence
CSR (`hyp_view(..., want_incidence = 1)`) so that, for vertex `i` in `VertexList`
order, the dual hyperedge is the list of hyperedge indices `ve[voff[i]..voff[i+1])`
— already ascending and de-duplicated, since the incidence is built from the
distinct-vertex sets. An isolated vertex gives an empty hyperedge. The result is
`Hypergraph[{1, ..., m}, {...}]`, so the dual of the dual recovers the incidence
structure on `1..n`.

**Data structures.** The incidence CSR `voff/ve` from the memo. The `1..m` vertex
`Integer`s are built once and shared (via `expr_copy`) into both the dual's vertex
List and its hyperedges; the dual is assembled with `mk_hyp`/`mk_list`.

**Complexity / limits.** `O(Σ|e|)` — linear in the total incidence, which equals
the total size of the dual's hyperedges. A bare List of hyperedges is not accepted
(the head needs the memoized incidence), so it is left unevaluated.

- Vertices are `1..m`; for each vertex of `h`, in VertexList order, the dual has
  the hyperedge of the (ascending) indices of the hyperedges containing it. An
  isolated vertex gives an empty hyperedge.
- Hyperedges are read as sets. The dual of the dual recovers the incidence
  structure on `1..n`.
- Linear in the total incidence. Benchmark (experiment 96): 10⁵ hyperedges
  8.1 ms warm / 10.1 ms cold (Mathematica 236 ms, xgi 361 ms).
- A bare List of hyperedges is not accepted (unevaluated).

**Attributes:** `Protected`.

## References

- C. Berge, *Hypergraphs: Combinatorics of Finite Sets*, North-Holland Mathematical Library 45 (Elsevier, 1989), ch. 1 (dual hypergraph).
- Source: [`src/graph/hyp_ops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_ops.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

The dual exchanges the roles of vertices and hyperedges. Its vertices are
`1, ..., m` (one per hyperedge of `h`), and for each vertex `v` of `h`, in
`VertexList` order, it has the hyperedge listing the ascending indices of the
hyperedges of `h` that contain `v`. An isolated vertex of `h` becomes an empty
hyperedge of the dual. Hence the dual of a hypergraph with `n` vertices and `m`
hyperedges has `m` vertices and `n` hyperedges.

Hyperedges are read as sets, so `HypergraphDual[HypergraphDual[h]]` recovers the
incidence structure of `h` on `1..n`. The build is linear in the total incidence.
A bare List of hyperedges is not accepted — the head needs the memoized incidence
of a real `Hypergraph` object.
