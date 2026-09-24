# BetweennessCentrality

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BetweennessCentrality[g] gives for each vertex v the sum over pairs of other vertices s, t of the fraction of shortest s-t paths passing through v (unordered pairs for undirected graphs, ordered pairs for directed ones). Brandes' algorithm; edge weights are ignored, as in Wolfram. Mixed graphs are not supported.`**

## Examples (4)

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

## Implementation notes

- Unnormalized; summed over unordered pairs on undirected graphs, ordered
  pairs on directed ones; weights ignored. Machine reals.
- **Deviation:** mixed graphs are left unevaluated (Wolfram's values there
  match no standard definition, e.g. `{1<->2, 2->3}` gives vertex 2 a
  betweenness of 0).
- **Algorithm:** Brandes, O(nm), with the sources spread over threads with
  per-thread accumulators.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
