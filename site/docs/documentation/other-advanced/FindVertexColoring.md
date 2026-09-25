# FindVertexColoring

!!! warning "Status: Partial"
    implemented with documented limitations or caveats; some argument forms fall through to symbolic/unevaluated output.

## Description

**`FindVertexColoring[g] gives a MINIMAL vertex colouring of g as a list of integers in VertexList order: the number of distinct colours equals the chromatic number, and no edge joins two vertices of equal colour. Minimality is proven by exact search (Wolfram's "BacktrackingDS" method) seeded by DSATUR upper and clique lower bounds. Exact colouring is NP-hard, so there are two guards, and exceeding EITHER returns the expression UNEVALUATED -- never a valid-but-larger colouring. (1) Graphs of more than 128 vertices are refused outright. (2) The search gives up after 8 million nodes, on the order of 100 seconds; a node count rather than a clock, so the answer does not depend on how fast the host is. To bound how long a call may take, wrap it in TimeConstrained -- that is the intended lever, and the search polls for the deadline so it is honoured; the 8-million-node ceiling is a last-resort backstop for an unattended run, not the responsiveness mechanism. Cost is driven by DENSITY, not by vertex count: a sparse 128-vertex graph answers instantly, while a dense one may exhaust the budget and refuse. FindVertexColoring[g, {c1, ...}] and FindVertexColoring[g, l] are not implemented.`**

## Examples

_No verified examples yet for this function._

## Algorithm

vertexcoloring.c - FindVertexColoring[g]: a MINIMAL vertex colouring.

Wolfram's FindVertexColoring returns a colouring whose number of distinct colours equals the chromatic number. That word "minimal" is the whole difficulty: computing it is NP-hard. A greedy or DSATUR-only implementation returns a valid-but-frequently-larger colouring, and it fails SILENTLY -- a plausible list of integers that quietly contradicts the documented meaning. So the search here is exact:

```text
  ub = fvc_dsatur_bound()  -- a good upper bound, cheap, and a real colouring
  lb = fvc_clique_bound()  -- a greedy-clique lower bound
  if lb == ub              -- ub is proven optimal; answer with no search
  else                     -- DSATUR branch-and-bound (fvc_bb), improving the
                              incumbent until it meets lb or the tree is
                              exhausted
```

The lower bound is not an optimisation. Without it CompleteGraph[128] -- under the vertex cap, so accepted -- would search instead of answering from the bounds, which for a complete graph is hopeless. With it lb == ub and the answer is immediate, at zero search nodes.

Neither bound makes the search cheap in general, though: see FVC_MAX_STEPS. Exactness is guaranteed by REFUSING (unevaluated) whenever it cannot be proven, never by returning a merely-valid colouring.

This is Wolfram's own "BacktrackingDS" method, so shipping only it is a documented subset rather than a divergence. No Method option is offered.

Adjacency is the UNDIRECTED neighbourhood: an edge constrains its endpoints whichever way it points. GraphAdj stores successors in out[] and predecessors in in[], with an UndirectedEdge contributing to both, so the neighbourhood of v is out[v] together with in[v] -- the same walk graph_count_components does. It is walked IN PLACE; no union structure is materialised, so there is nothing beyond the GraphAdj itself for a caller to free.

Memory (SPEC section 4): returns freshly-allocated results; the evaluator frees res.

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_slow.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_slow.c)
