# CirculantGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CirculantGraph[n, j] gives the circulant graph on n vertices with i joined to i+j and i-j (mod n); CirculantGraph[n, {j1, j2, ...}] uses every jump ji.`**

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= CirculantGraph[5, 1]
Out[1]= Graph[<5 vertices, 5 edges>]

In[2]:= EdgeList[CirculantGraph[6, 2]]
Out[2]= {1 <-> 3, 1 <-> 5, 2 <-> 4, 2 <-> 6, 3 <-> 5, 4 <-> 6}

In[3]:= EdgeList[CirculantGraph[6, {1, 3}]]
Out[3]= {1 <-> 2, 1 <-> 4, 1 <-> 6, 2 <-> 3, 2 <-> 5, 3 <-> 4, 3 <-> 6, 4 <-> 5, 5 <-> 6}

In[4]:= VertexDegree[CirculantGraph[8, {1, 2}]]
Out[4]= {4, 4, 4, 4, 4, 4, 4, 4}
```

### Applications (6)

A single offset of 1 gives the 6-cycle

```mathematica
In[5]:= EdgeList[CirculantGraph[6, 1]]
Out[5]= {1 <-> 2, 1 <-> 6, 2 <-> 3, 3 <-> 4, 4 <-> 5, 5 <-> 6}
```

The offset n/2 joins antipodes, each edge counted once

```mathematica
In[6]:= EdgeList[CirculantGraph[6, 3]]
Out[6]= {1 <-> 4, 2 <-> 5, 3 <-> 6}
```

A negative offset is reduced mod n, so this is the 5-cycle again

```mathematica
In[7]:= EdgeList[CirculantGraph[5, -1]]
Out[7]= {1 <-> 2, 1 <-> 5, 2 <-> 3, 3 <-> 4, 4 <-> 5}
```

A list of offsets gives the square of the 8-cycle

```mathematica
In[8]:= EdgeList[CirculantGraph[8, {1, 2}]]
Out[8]= {1 <-> 2, 1 <-> 3, 1 <-> 7, 1 <-> 8, 2 <-> 3, 2 <-> 4, 2 <-> 8, 3 <-> 4, 3 <-> 5, 4 <-> 5, 4 <-> 6, 5 <-> 6, 5 <-> 7, 6 <-> 7, 6 <-> 8, 7 <-> 8}
```

Each distinct offset contributes two to every degree

```mathematica
In[9]:= VertexDegree[CirculantGraph[10, {1, 3}]]
Out[9]= {4, 4, 4, 4, 4, 4, 4, 4, 4, 4}
```

The count is n times the number of offsets when none is n/2

```mathematica
In[10]:= EdgeCount[CirculantGraph[12, {1, 2, 3}]]
Out[10]= 36
```

## Implementation notes

**Algorithm.** `builtin_circulant_graph` takes exactly two arguments: a positive machine-integer vertex count `n` and either one integer offset `j` or a `List` of integer offsets. Each offset is reduced into `0..n-1` (`j mod n`, negatives wrapped), and for every vertex `i` and every non-zero offset `s` the edge `i ~ (i + s) mod n` is emitted. An offset of `0` (or a multiple of `n`) contributes nothing, and an offset `s` and its mirror `n - s` give the same edges.

**Data structures.** Edges are accumulated as packed 64-bit keys `(min << 32) | max` in the shared `Pairs` buffer (`pairs_add`). `pairs_graph` sorts them lexicographically, drops duplicates (for example `s = n/2` met from both ends), and builds the canonical `Graph[Range[n], {UndirectedEdge[i, j], ...}]` expression tree. All edges share one `UndirectedEdge` head node and reuse the vertex integer nodes by reference.

**Complexity / limits.** `O(n * |offsets|)` edge generation plus an `O(E log E)` sort. Calls with a non-integer `n` or offset, a non-positive `n`, `n` above 10^8 vertices, or more than 5 x 10^7 edges are left unevaluated. Options are not supported; the result is always undirected.

- Undirected on `1..n`; the edge list is the sorted list of pairs `{i, j}`,
  `i < j`, identical to Wolfram's `EdgeList`.
- Options (`DirectedEdges`, layout options) are not supported.
- Resource limit: more than 10^8 vertices or 5×10^7 edges is left
  unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/gmet_generators.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_generators.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

`CirculantGraph[n, j]` joins `i` to `i ± j (mod n)` on vertices `1..n`; the offset may be one integer or a list. Offsets are reduced modulo `n`, so `j` and `n - j` give the same graph, and an offset that is a multiple of `n` adds nothing.

The result is an undirected graph whose edge list is sorted lexicographically. Non-integer arguments, or a non-positive `n`, leave the call unevaluated.
