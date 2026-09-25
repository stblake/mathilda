# EdgeBetweennessCentrality

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeBetweennessCentrality[g] gives for each edge, in EdgeList order, the sum over ordered pairs of vertices s, t of the fraction of shortest s-t paths using that edge. Uses EdgeWeight as lengths.`**

## Examples (4)

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

## Implementation notes

- *(w)* weight-aware; machine reals.
- Summed over **ordered** pairs even on undirected graphs (the edge of a `K2`
  scores 2).
- Weighted: ties within a relative `1e-12` share paths.
- **Algorithm:** Brandes, O(nm) unweighted, O(nm + n² log n) weighted, sources
  spread over threads with per-thread accumulators.

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
