# RandomGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RandomGraph[{n, m}] gives a random undirected graph with n vertices and m edges. RandomGraph[{n, m}, k] gives a list of k such graphs; memory use grows with k.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= SeedRandom[42]; EdgeList[RandomGraph[{5, 4}]]
Out[1]= {3 <-> 5, 1 <-> 5, 2 <-> 5, 3 <-> 4}

In[2]:= Length[RandomGraph[{6, 5}, 3]]
Out[2]= 3

In[3]:= RandomGraph[{3, 4}]
Out[3]= RandomGraph[{3, 4}]

In[4]:= RandomGraph[{6, 5}, 0]
Out[4]= {}

In[5]:= Head[RandomGraph[{6, 5}, 1]]
Out[5]= List

In[6]:= RandomGraph[{6, 5}, -1]
Out[6]= RandomGraph[{6, 5}, -1]
```

## Implementation notes

- `Protected`. Uses the seeded system RNG, so `SeedRandom` makes it
  reproducible. Vertices are `1..n`, edges undirected (see `CycleGraph`).
- Returns unevaluated if `m` exceeds `n(n-1)/2`. For `n <= 1` the only possible
  graph is the edgeless one, which is what is returned.
- `k = 0` gives `{}`; `k = 1` gives a one-element list, **not** a bare `Graph`.
  A negative, non-integer, or symbolic `k` leaves the expression unevaluated,
  silently — the convention the count-taking `Random*` heads share.
- Algorithm: the `n(n-1)/2` candidate edges are never materialised. Each graph
  draws exactly what `RandomSample` over the row-major candidate list would (so
  seeded output is that of `RandomSample`), in `O(m)` time and memory per graph
  — `RandomGraph[{200000, 300000}]` takes about 80 ms. (Until v0.185 each
  element built the full `O(n^2)` candidate list, and leaked it.)

**Attributes:** `Protected`.

## References

**See also:** [SeedRandom](../../random-number-generation/SeedRandom/), [CycleGraph](../../graphs/CycleGraph/), [Graph](../../graphs/Graph/), [RandomSample](../../random-number-generation/RandomSample/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_slow.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_slow.c)
