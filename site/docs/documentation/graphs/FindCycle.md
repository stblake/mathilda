# FindCycle

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindCycle[g] finds a cycle of g, as {list of edges}, or {}. FindCycle[g, k] finds a cycle of length at most k, FindCycle[g, {k}] of length exactly k, FindCycle[g, {kmin, kmax}] of length in that range; FindCycle[g, kspec, n] finds at most n (n an integer or All). Each cycle is reported once. Length-bounded searches are exhaustive and give up (unevaluated) after 5*10^7 steps; use TimeConstrained to bound them.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= FindCycle[Graph[{1->2,2->3,3->4,4->2,3->1}]]
Out[1]= {{1 -> 2, 2 -> 3, 3 -> 1}}

In[2]:= FindCycle[PathGraph[Range[4]]]
Out[2]= {}

In[3]:= FindCycle[CompleteGraph[4], 3, All]
Out[3]= {{1 <-> 2, 2 <-> 3, 3 <-> 1}, {1 <-> 2, 2 <-> 4, 4 <-> 1}, {1 <-> 3, 3 <-> 4, 4 <-> 1}, {2 <-> 3, 3 <-> 4, 4 <-> 2}}

In[4]:= FindCycle[CompleteGraph[4], {3,4}, 2]
Out[4]= {{1 <-> 2, 2 <-> 3, 3 <-> 1}, {1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}}

In[5]:= FindCycle[{Graph[{1<->2,2<->3,3<->1,3<->4,4<->5,5<->3}], 4}]
Out[5]= {{4 <-> 3, 3 <-> 5, 5 <-> 4}}

In[6]:= FindCycle[Graph[{1->2,2<->3}]]
Out[6]= FindCycle[Graph[<3 vertices, 2 edges>]]
```

## Implementation notes

- `Protected`. A non-graph argument is left unevaluated.
- Plain form: linear time, by Mathematica's own search order (a stack DFS that
  scans a vertex's edges when visiting it), so the reported cycle matches
  Mathematica's.
- `FindCycle[{g, v}]` (plain form) is found by BFS from `v`, in linear time.
- Each cycle is reported once. Cycle length counts edges; undirected cycles have
  length `>= 3`, directed `>= 2`.
- The length-bounded / enumerating forms backtrack from each vertex as the
  cycle's lowest vertex — exponential in the worst case (a cycle of exact length
  `n` is a Hamiltonian cycle) — poll `TimeConstrained`, and give up (unevaluated)
  after `5*10^7` steps. Their choice and order of cycles is Mathilda's own
  (Mathematica's differs).
- Mixed graphs, and weighted graphs with a length spec, are left unevaluated.
- Reuses the per-graph incidence-list cache shared with `FindPath` and
  `FindEulerianCycle`.

**Attributes:** `Protected`.

## References

**See also:** [TimeConstrained](../../time-and-date/TimeConstrained/), [FindPath](../../graphs/FindPath/), [FindEulerianCycle](../../graphs/FindEulerianCycle/)

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
