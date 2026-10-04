# EdgeBetweennessCentrality

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeBetweennessCentrality[g] gives for each edge, in EdgeList order, the sum over ordered pairs of vertices s, t of the fraction of shortest s-t paths using that edge. Uses EdgeWeight as lengths.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= EdgeBetweennessCentrality[PathGraph[{1,2,3,4}]]
Out[1]= {6.0, 8.0, 6.0}

In[2]:= EdgeBetweennessCentrality[PathGraph[2]]
Out[2]= {2.0}

In[3]:= EdgeBetweennessCentrality[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[3]= {4.0, 5.0, 3.0, 3.0}
```

### Options (1)

```mathematica
In[4]:= EdgeBetweennessCentrality[Graph[{1,2,3,4},{1<->2,2<->3,3<->4,1<->4}, EdgeWeight->{1,1,1,5}]]
Out[4]= {6.0, 8.0, 6.0, 0.0}
```

### Applications (2)

The central edge is on the most shortest paths

```mathematica
In[5]:= EdgeBetweennessCentrality[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4}]]
Out[5]= {6.0, 8.0, 6.0}
```

In EdgeList order

```mathematica
In[6]:= EdgeBetweennessCentrality[Graph[{1 -> 2, 2 -> 3, 1 -> 3, 3 -> 4}]]
Out[6]= {1.0, 2.0, 2.0, 3.0}
```

## Implementation notes

**Algorithm.** `builtin_edge_betweenness_centrality` is the edge variant of **Brandes'
algorithm**: the same reverse-order dependency accumulation, but each unit of dependency is
charged to `acc[eid]` — the id of the arc carrying it — instead of to a vertex. Unlike the
vertex form, it *does* use `EdgeWeight` as edge lengths, so when the graph carries usable
weights it takes the Dijkstra variant (`gmet_csr_build(..., ew, want_eid=1, ...)`), and with
equal path lengths within a relative tolerance `BRANDES_TIE = 1e-12` it splits the flow across
the tied shortest paths. It sums over ordered pairs on undirected graphs too (no ×0.5),
matching the convention of the reference implementation.

**Data structures.** The per-source `GmetCSR` carries an `eid[]` array mapping each arc back to
its `EdgeList` position; the accumulator has one slot per edge rather than per vertex. The rest
is shared with the vertex form: `dist[]`/`sigma[]`/`delta[]` and the settled-vertex stack, run
per source across a thread team and merged. The result is a packed machine-real vector in
`EdgeList` order.

**Complexity / limits.** `O(nm)` unweighted, `O(nm + n^2 log n)` weighted, `n = |V|`,
`m = |E|`. No size cap; non-graph input returns unevaluated.

- *(w)* weight-aware; machine reals.
- Summed over **ordered** pairs even on undirected graphs (the edge of a `K2`
  scores 2).
- Weighted: ties within a relative `1e-12` share paths.
- **Algorithm:** Brandes, O(nm) unweighted, O(nm + n² log n) weighted, sources
  spread over threads with per-thread accumulators.

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- U. Brandes, *A faster algorithm for betweenness centrality*, J. Math. Sociology **25** (2001) 163-177.
- Source: [`src/graph/gmet_centrality.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_centrality.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

The result is a list of values in `EdgeList` order, one per edge: the sum over ordered pairs
of vertices `s`, `t` of the fraction of shortest `s`-`t` paths that use that edge.

Unlike the vertex `BetweennessCentrality`, this head uses `EdgeWeight` as edge lengths when the
graph carries usable weights, and it always sums over ordered pairs (there is no undirected ×1/2
factor). Tied shortest paths split the flow evenly.
