# FindCycle

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindCycle[g] finds a cycle of g, as {list of edges}, or {}. FindCycle[g, k] finds a cycle of length at most k, FindCycle[g, {k}] of length exactly k, FindCycle[g, {kmin, kmax}] of length in that range; FindCycle[g, kspec, n] finds at most n (n an integer or All). Each cycle is reported once. Length-bounded searches are exhaustive and give up (unevaluated) after 5*10^7 steps; use TimeConstrained to bound them.`**

## Examples (13)

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

### Applications (7)

The only cycle of a pentagon

```mathematica
In[7]:= FindCycle[CycleGraph[5]]
Out[7]= {{1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 5, 5 <-> 1}}
```

A tree has no cycle, so the result is empty

```mathematica
In[8]:= FindCycle[PathGraph[{1, 2, 3}]]
Out[8]= {}
```

A cycle of exactly four edges

```mathematica
In[9]:= FindCycle[CompleteGraph[4], {4}]
Out[9]= {{1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}}
```

Every triangle, each reported once

```mathematica
In[10]:= FindCycle[CompleteGraph[4], {3}, All]
Out[10]= {{1 <-> 2, 2 <-> 3, 3 <-> 1}, {1 <-> 2, 2 <-> 4, 4 <-> 1}, {1 <-> 3, 3 <-> 4, 4 <-> 1}, {2 <-> 3, 3 <-> 4, 4 <-> 2}}
```

A cycle through a given vertex

```mathematica
In[11]:= FindCycle[{CycleGraph[5], 2}]
Out[11]= {{2 <-> 1, 1 <-> 5, 5 <-> 4, 4 <-> 3, 3 <-> 2}}
```

Directed cycle

```mathematica
In[12]:= FindCycle[Graph[{1 -> 2, 2 -> 3, 3 -> 1}]]
Out[12]= {{1 -> 2, 2 -> 3, 3 -> 1}}
```

Up to two cycles of any length

```mathematica
In[13]:= FindCycle[CompleteGraph[4], Infinity, 2]
Out[13]= {{1 <-> 2, 2 <-> 3, 3 <-> 1}, {1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}}
```

## Implementation notes

**Algorithm.** `builtin_find_cycle` takes `FindCycle[g]`, `FindCycle[g, kspec]`, `FindCycle[g, kspec, n]` and `FindCycle[{g, v}, ...]`, with `kspec` being `k`, `Infinity`, `{k}` or `{kmin, kmax}` (cycle length counts edges) and `n` a positive integer or `All`. The plain one-cycle form is a stack-based DFS (`dfs_any_cycle`) that reports the first back edge to a vertex on the current root path, written from that ancestor, following Mathematica's search order. With `{g, v}` it runs a BFS from `v` (`bfs_cycle_through`). Length-bounded and multi-cycle forms backtrack from each vertex as the lowest-positioned vertex of the cycle, so each cycle is reported once rather than once per rotation or direction. Undirected cycles need length 3 or more, directed ones length 2 or more.

**Data structures.** The search runs over the per-graph edge-index views (`eu`/`ev`/`directed`) and a CSR incidence structure, with `int` arrays for the DFS stack, visited marks and current path; results accumulate in a growable list of cycle expressions, each an edge list.

**Complexity / limits.** `O(V + E)` for the single-cycle and through-a-vertex forms. The length-bounded and enumerating forms are exponential in the worst case (a cycle of exact length k is NP-hard), so they poll `TimeConstrained` and give up, leaving the call unevaluated, after a fixed step budget. Mixed directed/undirected graphs, and weighted graphs with a length spec, are left unevaluated. No cycle gives `{}`. For the bounded forms the choice and order of cycles is Mathilda's own deterministic one.

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

- Source: [`src/graph/gops_cycles.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_cycles.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

A cycle is returned as a list of edges, and the result is a list of cycles: `{}` when there is none. A length specification counts edges and may be `k` (at most), `{k}` (exactly), `{kmin, kmax}` or `Infinity`; a third argument asks for up to `n` cycles or `All`.

The plain form is linear time. Exact-length and enumerating forms are exponential in the worst case, so they give up under a step budget and then stay unevaluated. For those forms, which cycles are reported first is Mathilda's own deterministic choice.
