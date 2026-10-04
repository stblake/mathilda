# HararyGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HararyGraph[k, n] gives the Harary graph: a k-connected graph on n vertices with the minimum number of edges (k >= 2, n > k).`**

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= HararyGraph[2, 5]
Out[1]= Graph[<5 vertices, 5 edges>]

In[2]:= EdgeList[HararyGraph[3, 6]]
Out[2]= {1 <-> 2, 1 <-> 4, 1 <-> 6, 2 <-> 3, 2 <-> 5, 3 <-> 4, 3 <-> 6, 4 <-> 5, 5 <-> 6}

In[3]:= EdgeList[HararyGraph[3, 5]]
Out[3]= {1 <-> 2, 1 <-> 3, 1 <-> 4, 1 <-> 5, 2 <-> 3, 2 <-> 5, 3 <-> 4, 4 <-> 5}

In[4]:= HararyGraph[1, 5]
Out[4]= HararyGraph[1, 5]

In[5]:= HararyGraph[4, 4]
Out[5]= HararyGraph[4, 4]
```

### Applications (4)

Ceil(k n / 2) = 12, the minimum for 3-connectivity

```mathematica
In[6]:= EdgeCount[HararyGraph[3, 8]]
Out[6]= 12
```

Connectivity is exactly k

```mathematica
In[7]:= VertexConnectivity[HararyGraph[4, 7]]
Out[7]= 4
```

N vertices, labelled 1..n

```mathematica
In[8]:= VertexCount[HararyGraph[4, 7]]
Out[8]= 7
```

Nearly regular of degree k

```mathematica
In[9]:= VertexDegree[HararyGraph[3, 8]]
Out[9]= {3, 3, 3, 3, 3, 3, 3, 3}
```

## Implementation notes

**Algorithm.** `builtin_harary_graph` constructs `HararyGraph[k, n]`, the
`k`-connected graph on `n` vertices with the fewest possible edges (Harary 1962).
With `r = ⌊k/2⌋` the core is the circulant `C_n(1, 2, …, r)`: each vertex joined to
its `r` nearest neighbours on each side of the cycle `0, 1, …, n−1`
(`circulant_add`). When `k` is odd two cases add the last ring of edges:
`n` even adds the `n/2` "diameter" chords `i — i+n/2`; `n` odd adds a near-perfect
matching across the half-cycle (`0 — h`, `0 — h+1`, and `i — i+h+1` for
`i = 1 … h−1`, with `h = (n−1)/2`). Vertices are emitted to `pairs_graph`, which
renumbers them to `1 … n` and deduplicates parallel edges, yielding an undirected
graph that is exactly `k`-vertex-connected with `⌈kn/2⌉` edges.

**Data structures.** A `Pairs` edge accumulator (a growable array of endpoint
pairs) filled by `circulant_add`/`pairs_add`, then turned into the `Graph`
expression by `pairs_graph`, which sorts and dedups so the doubled edges of small
cases collapse.

**Complexity / limits.** `O(kn)` edges are emitted and the build is linear in that
size. The arguments must satisfy `k ≥ 0` and `n > k` (a `k`-connected graph needs
at least `k+1` vertices); `n` is capped at `GEN_MAX_VERTICES`. Anything outside
those bounds leaves the call unevaluated. Note the argument order is
`HararyGraph[k, n]` — connectivity first, vertex count second.

- Requires `k >= 2` and `n > k`; otherwise unevaluated.
- Undirected on `1..n`; the edge list is the sorted list of pairs `{i, j}`,
  `i < j`, identical to Wolfram's `EdgeList`.
- Options (`DirectedEdges`, layout options) are not supported.
- Resource limit: more than 10^8 vertices or 5×10^7 edges is left
  unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- F. Harary, *The maximum connectivity of a graph*, Proc. Natl. Acad. Sci. USA **48** (1962) 1142-1146.
- Source: [`src/graph/gmet_generators.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_generators.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

`HararyGraph[k, n]` is the extremal graph of Harary's theorem: the `k`-connected
graph on `n` vertices with as few edges as possible, namely `⌈kn/2⌉`. It is built
as the circulant on the `⌊k/2⌋` nearest neighbours of a cycle, plus one extra ring
of chords when `k` is odd. The argument order is connectivity first, so
`HararyGraph[4, 7]` is 4-connected on 7 vertices.

The result is an opaque `Graph` object; query a property (`EdgeCount`,
`VertexConnectivity`, `VertexDegree`) to inspect it. The construction requires
`n > k`, since a `k`-connected graph needs at least `k + 1` vertices; smaller `n`
is left unevaluated.
