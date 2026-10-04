# BetweennessCentrality

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BetweennessCentrality[g] gives for each vertex v the sum over pairs of other vertices s, t of the fraction of shortest s-t paths passing through v (unordered pairs for undirected graphs, ordered pairs for directed ones). Brandes' algorithm; edge weights are ignored, as in Wolfram. Mixed graphs are not supported.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= BetweennessCentrality[StarGraph[5]]
Out[1]= {6.0, 0.0, 0.0, 0.0, 0.0}

In[2]:= BetweennessCentrality[PathGraph[4]]
Out[2]= {0.0, 2.0, 2.0, 0.0}

In[3]:= BetweennessCentrality[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[3]= {1.0, 2.0, 3.0, 0.0}

In[4]:= BetweennessCentrality[Graph[{1<->2, 2->3}]]
Out[4]= BetweennessCentrality[Graph[<3 vertices, 2 edges>]]
```

### Applications (2)

The interior of a path carries every cross-path

```mathematica
In[5]:= BetweennessCentrality[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 5}]]
Out[5]= {0.0, 3.0, 4.0, 3.0, 0.0}
```

Directed: ordered pairs, no 1/2 factor

```mathematica
In[6]:= BetweennessCentrality[Graph[{1 -> 2, 2 -> 3, 1 -> 3, 3 -> 4}]]
Out[6]= {0.0, 0.0, 2.0, 0.0}
```

## Implementation notes

**Algorithm.** `builtin_betweenness_centrality` computes the (unnormalised) Freeman
betweenness of every vertex with **Brandes' algorithm**: one shortest-path exploration per
source accumulates each vertex's dependency `delta` by walking the settled vertices back in
reverse order of discovery, so the pair-counting is done without materialising predecessor
lists. The vertex form deliberately ignores `EdgeWeight` — it builds its CSR with
`NULL` weights and always takes the unweighted BFS variant (`brandes_bfs_work`); the sibling
edge form accumulates onto edge ids instead and does use the Dijkstra variant
(`brandes_dijkstra_work`, binary heap `GmetHeap`) when the graph carries usable weights. A
mixed directed/undirected graph is left unevaluated, because no standard betweenness
definition fits it.

**Data structures.** Its own compressed-row adjacency `GmetCSR` (`off[]`/`adj[]`, out-arcs),
not the shared `GraphAdj`. Per source it fills `dist[]` (BFS layers), `sigma[]` (shortest-path
counts, as `double`), `delta[]` (dependency), and an `order[]` stack of settled vertices.
Sources are split across a thread team (up to 32) with per-thread accumulators merged at the
end. An undirected graph counts each unordered pair once (a final ×0.5); a directed graph
counts ordered pairs. Results come back as a packed machine-real vector.

**Complexity / limits.** `O(nm)` for the unweighted (BFS) path and `O(nm + n^2 log n)` for the
weighted Dijkstra path, `n = |V|`, `m = |E|`. No size cap. Non-graph or mixed input returns
unevaluated.

- Unnormalized; summed over unordered pairs on undirected graphs, ordered
  pairs on directed ones; weights ignored. Machine reals.
- **Deviation:** mixed graphs are left unevaluated (Wolfram's values there
  match no standard definition, e.g. `{1<->2, 2->3}` gives vertex 2 a
  betweenness of 0).
- **Algorithm:** Brandes, O(nm), with the sources spread over threads with
  per-thread accumulators.

**Attributes:** `Protected`.

## References

- U. Brandes, *A faster algorithm for betweenness centrality*, J. Math. Sociology **25** (2001) 163-177.
- Source: [`src/graph/gmet_centrality.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_centrality.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

The value is unnormalised Freeman betweenness: for each vertex `v`, the sum over pairs of
other vertices `s`, `t` of the fraction of shortest `s`-`t` paths that pass through `v`.
Undirected graphs count each pair once; directed graphs count ordered pairs.

The vertex form ignores `EdgeWeight` and always measures shortest paths by edge count. If you
need weighted shortest paths for the pair fractions, that distinction lives in
`EdgeBetweennessCentrality`, which does honour `EdgeWeight`. A mixed graph (both directed and
undirected edges) is left unevaluated.
