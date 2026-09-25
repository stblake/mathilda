# KatzCentrality

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KatzCentrality[g, a] gives the Katz centrality x = a A^T x + 1 of each vertex; KatzCentrality[g, a, b] uses b (a number or a list) in place of 1. Edge weights are ignored.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= KatzCentrality[PathGraph[3], 0.1]
Out[1]= {1.12245, 1.22449, 1.12245}

In[2]:= KatzCentrality[CompleteGraph[4], 1/2]
Out[2]= {-2.0, -2.0, -2.0, -2.0}

In[3]:= KatzCentrality[PathGraph[3], 0.1, 2]
Out[3]= {2.2449, 2.44898, 2.2449}

In[4]:= KatzCentrality[PathGraph[3], 0.1, {1, 2, 3}]
Out[4]= {1.2449, 2.44898, 3.2449}

In[5]:= KatzCentrality[PathGraph[3], 0]
Out[5]= {1, 1, 1}

In[6]:= KatzCentrality[Graph[{1,2,3},{}], 0.3]
Out[6]= {1, 1, 1}

In[7]:= KatzCentrality[CompleteGraph[4], 1/3]
Out[7]= KatzCentrality[Graph[<4 vertices, 6 edges>], 1/3]
```

## Implementation notes

- Solves `x = a A^T x + b`, with `b` = 1 by default, a number, or a list.
  Weights ignored.
- Answered even beyond the convergence radius, like Wolfram (`K4` with
  `a = 1/2` gives `-2`); uses a dense LAPACK solve, `n <= 4000`. A system
  singular up to rounding (LU pivot ratio below `1e-12`, e.g. `K4` with
  `a = 1/3`) is left unevaluated. Mathematica 15 returns rounding noise there
  (about `1.8*10^16` per vertex for `K4`, `a = 1/3`).
- No edges, or an exact `a = 0`, return `b` exactly.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
