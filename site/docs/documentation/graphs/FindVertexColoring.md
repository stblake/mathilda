# FindVertexColoring

!!! warning "Status: Partial"
    implemented with documented limitations or caveats; some argument forms fall through to symbolic/unevaluated output.

## Description

**`FindVertexColoring[g] gives a MINIMAL vertex colouring of g as a list of integers in VertexList order: the number of distinct colours equals the chromatic number, and no edge joins two vertices of equal colour. Minimality is proven by exact search (Wolfram's "BacktrackingDS" method) seeded by DSATUR upper and clique lower bounds. Exact colouring is NP-hard, so there are two guards, and exceeding EITHER returns the expression UNEVALUATED -- never a valid-but-larger colouring. (1) Graphs of more than 128 vertices are refused outright. (2) The search gives up after 8 million nodes, on the order of 100 seconds; a node count rather than a clock, so the answer does not depend on how fast the host is. To bound how long a call may take, wrap it in TimeConstrained -- that is the intended lever, and the search polls for the deadline so it is honoured; the 8-million-node ceiling is a last-resort backstop for an unattended run, not the responsiveness mechanism. Cost is driven by DENSITY, not by vertex count: a sparse 128-vertex graph answers instantly, while a dense one may exhaust the budget and refuse. FindVertexColoring[g, {c1, ...}] and FindVertexColoring[g, l] are not implemented.`**

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= FindVertexColoring[CycleGraph[5]]
Out[1]= {1, 2, 1, 2, 3}

In[2]:= FindVertexColoring[CycleGraph[4]]
Out[2]= {1, 2, 1, 2}

In[3]:= FindVertexColoring[CompleteGraph[4]]
Out[3]= {1, 2, 3, 4}

In[4]:= FindVertexColoring[Graph[{1,2,3},{}]]
Out[4]= {1, 1, 1}

In[5]:= FindVertexColoring[Graph[{},{}]]
Out[5]= {}

In[6]:= FindVertexColoring[5]
Out[6]= FindVertexColoring[5]
```

### Applications (5)

A clique needs a distinct colour per vertex

```mathematica
In[7]:= FindVertexColoring[CompleteGraph[4]]
Out[7]= {1, 2, 3, 4}
```

An even cycle is 2-chromatic

```mathematica
In[8]:= FindVertexColoring[CycleGraph[4]]
Out[8]= {1, 2, 1, 2}
```

An odd cycle needs three colours

```mathematica
In[9]:= FindVertexColoring[CycleGraph[5]]
Out[9]= {1, 2, 1, 2, 3}
```

A path is bipartite

```mathematica
In[10]:= FindVertexColoring[PathGraph[{1, 2, 3, 4}]]
Out[10]= {2, 1, 2, 1}
```

The chromatic number of the Petersen graph is 3

```mathematica
In[11]:= Max[FindVertexColoring[PetersenGraph[]]]
Out[11]= 3
```

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

**Algorithm.** `builtin_find_vertex_coloring` returns a colour assignment whose
number of distinct colours equals the **chromatic number** — Wolfram's documented
meaning, and the reason the search is exact rather than greedy: a merely-valid
colouring with too many colours would be a plausible list of integers that
silently contradicts that contract. It brackets the chromatic number between a
cheap upper bound `ub` from a **DSATUR** pass (`fvc_dsatur_bound`, which also
exhibits a real colouring) and a lower bound `lb` from a multi-start greedy clique
(`fvc_clique_bound`). When `lb == ub` the DSATUR colouring is proven optimal and
returned with zero search — this is what makes `CompleteGraph[128]` immediate
rather than a hang. Otherwise it runs a **DSATUR branch-and-bound** (`fvc_bb`) that
picks the next vertex dynamically as the uncoloured one of maximum saturation and
prunes on three grounds: bound (a partial colouring already at `best_k` colours
cannot win), symmetry breaking (a vertex tries only colours `1..used+1`), and
optimality (stop once `best_k == lb`). Adjacency is the **undirected**
neighbourhood — an edge constrains its endpoints whichever way it points. The empty
graph colours to `{}`. Exactness is preserved by **refusing** (returning `NULL`,
leaving the call unevaluated) whenever minimality cannot be proven: above
`FVC_MAX_VERTICES = 128`, when the node budget `FVC_MAX_STEPS = 8,000,000` is
spent, or on allocation failure — never by returning the incumbent. A second
argument also declines (the `FindVertexColoring[g, {c1, ...}]` forms are a later
layer). This is Wolfram's own `"BacktrackingDS"` method, so shipping only it is a
documented subset.

**Data structures.** `graph_build_adj` yields a `GraphAdj` with successor `out[]`
and predecessor `in[]` lists (an undirected edge appears in both), walked in place
so there is nothing beyond the `GraphAdj` to free. The working colouring `col[]`
is 1-based (0 = uncoloured); the branch-and-bound threads an `FvcBB` struct
(working colouring, best complete colouring, incumbent `best_k`, lower bound `lb`)
by pointer, with `seen[]`/`forbid[]` as stack arrays sized to the vertex cap so the
hot path never allocates. The result is a `List` of positive-integer colour labels,
one per vertex in `VertexList` order.

**Complexity / limits.** Exact graph colouring is NP-hard; the vertex cap bounds
size and the **node count** bounds cost *deterministically* — a wall-clock cutoff
would make the answer machine-dependent (a fast host proves minimality, a slow one
refuses the same graph), whereas a fixed node budget gives every machine the same
answer. `fvc_bb` polls `tc_check_deadline()` every 4096 nodes so an interactive
`TimeConstrained[...]` is honoured even where `SIGPROF` is unreliable; the node
budget is only the backstop for an unattended run.

- `Protected`. The number of distinct colours equals the chromatic number, and
  no edge joins two vertices of equal colour. Edge direction is ignored.
- Minimality is proven by exact search (Wolfram's `"BacktrackingDS"` method)
  seeded by DSATUR upper and clique lower bounds.
- Exact colouring is NP-hard, so there are two guards, and exceeding either
  returns the expression unevaluated — never a valid-but-larger colouring:
  graphs of more than 128 vertices are refused outright, and the search gives
  up after 8 million nodes (on the order of 100 seconds). The limit is a node
  count rather than a clock, so the answer does not depend on host speed.
- To bound how long a call may take, wrap it in `TimeConstrained` — that is the
  intended lever, and the search polls for the deadline. The node ceiling is a
  last-resort backstop for unattended runs.
- Cost is driven by density, not vertex count: a sparse 128-vertex graph
  answers instantly, while a dense one may exhaust the budget and refuse.
- `FindVertexColoring[g, {c1, ...}]` and `FindVertexColoring[g, l]`
  (Mathematica's colour-list forms) are not implemented and stay unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [VertexList](../../graphs/VertexList/), [TimeConstrained](../../time-and-date/TimeConstrained/)

- D. Brélaz, *New methods to color the vertices of a graph*, Communications of the ACM **22** (1979) 251-256.
- Source: [`src/graph/vertexcoloring.c`](https://github.com/stblake/mathilda/blob/main/src/graph/vertexcoloring.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_slow.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_slow.c)

## Notes & additional examples

### Notes

`FindVertexColoring[g]` returns a list of positive-integer colour labels, one per
vertex in `VertexList` order, such that adjacent vertices differ — and, crucially,
using as few distinct colours as possible. The number of distinct colours is the
chromatic number of `g`, so `Max` of the result reads off that number. An edge
constrains its endpoints in either direction, so direction is ignored.

Because minimal colouring is NP-hard, the search is exact and may **refuse** rather
than return a merely-valid colouring: the call is left unevaluated for a graph of
more than 128 vertices, or when an internal node budget is exhausted before
minimality can be proven. Bound an interactive call with `TimeConstrained` if
needed. Only the one-argument form is supported.
