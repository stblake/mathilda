# RandomGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RandomGraph[{n, m}] gives a random undirected graph with n vertices and m edges. RandomGraph[{n, m}, k] gives a list of k such graphs; memory use grows with k.`**

## Examples (9)

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

### Applications (3)

```mathematica
In[7]:= r = RandomGraph[{6, 8}];
```

G(n, m): exactly n vertices and m edges

```mathematica
In[8]:= {VertexCount[r], EdgeCount[r]}
Out[8]= {6, 8}
```

A list of k independent graphs

```mathematica
In[9]:= Length[RandomGraph[{4, 3}, 5]]
Out[9]= 5
```

## Implementation notes

**Algorithm.** `builtin_random_graph` samples from the Erdős–Rényi `G(n, m)` model: a simple
undirected graph on vertices `1..n` with exactly `m` edges, every `m`-edge graph equally likely.
The `m` edges are drawn without replacement from the `n(n-1)/2` candidate pairs, but the
candidate list is never materialised: `random_sample_indices(maxe, m)` produces exactly the
draws `RandomSample` would (so it honours `SeedRandom`), and each sampled index is decoded to a
vertex pair by binary search over the triangular row offsets. The assembled `Graph[...]` is
re-validated by the evaluator's `builtin_graph`.

**Data structures.** Only the `m` sampled indices and the resulting edge `List`; the decode is
`O(m log n)` time and `O(m)` memory — no `O(n^2)` candidate array.

**Complexity / limits.** `O(m log n)`. Guards: `n` or `m` negative returns unevaluated; `n`
above `2^31-1` returns unevaluated (overflow); `m` greater than `n(n-1)/2` returns unevaluated.
`RandomGraph[{n, m}, k]` returns a list of `k` independent graphs (`k = 0` gives `{}`); on any
sub-failure the whole list is discarded. Always undirected.

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

- P. Erdős and A. Rényi, *On random graphs I*, Publ. Math. Debrecen **6** (1959) 290-297.
- Source: [`src/graph/generators.c`](https://github.com/stblake/mathilda/blob/main/src/graph/generators.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_slow.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_slow.c)
- Tests: [`tests/test_graphplot.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphplot.c)

## Notes & additional examples

### Notes

`RandomGraph[{n, m}]` draws a uniformly random simple undirected graph on `n` vertices with
exactly `m` edges (the Erdős–Rényi `G(n, m)` model). The vertex and edge counts are therefore
fixed by the arguments; the structure is what varies, and it is reproducible under
`SeedRandom`.

`RandomGraph[{n, m}, k]` returns a list of `k` independent such graphs. `m` may not exceed
`n(n-1)/2`, the number of possible edges.
