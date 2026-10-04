# PlanarGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PlanarGraphQ[g] gives True if g can be drawn in the plane without edge crossings. Linear-time left-right planarity test; gives False for non-graphs.`**

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= PlanarGraphQ[CompleteGraph[4]]
Out[1]= True

In[2]:= PlanarGraphQ[CompleteGraph[5]]
Out[2]= False

In[3]:= PlanarGraphQ[CompleteGraph[{3,3}]]
Out[3]= False

In[4]:= PlanarGraphQ[PetersenGraph[]]
Out[4]= False

In[5]:= PlanarGraphQ[Graph[{1->2,2->3,3->1}]]
Out[5]= True

In[6]:= PlanarGraphQ[x]
Out[6]= False
```

### Applications (7)

K4 draws in the plane without crossings

```mathematica
In[7]:= PlanarGraphQ[CompleteGraph[4]]
Out[7]= True
```

K5 is the first forbidden minor

```mathematica
In[8]:= PlanarGraphQ[CompleteGraph[5]]
Out[8]= False
```

So is the utility graph k3,3

```mathematica
In[9]:= PlanarGraphQ[CompleteGraph[{3, 3}]]
Out[9]= False
```

Grids are planar

```mathematica
In[10]:= PlanarGraphQ[GridGraph[{3, 3}]]
Out[10]= True
```

Ten vertices, fifteen edges, but not planar

```mathematica
In[11]:= PlanarGraphQ[PetersenGraph[]]
Out[11]= False
```

Direction is ignored

```mathematica
In[12]:= PlanarGraphQ[Graph[{1, 2, 3}, {1 -> 2, 2 -> 3}]]
Out[12]= True
```

A non-graph gives false

```mathematica
In[13]:= PlanarGraphQ[x]
Out[13]= False
```

## Implementation notes

**Algorithm.** `builtin_planar_graph_q` gives `False` for a non-graph and otherwise runs `galg_planar_test` on the underlying simple undirected graph (direction and anti-parallel pairs do not matter). It is the left-right planarity test of de Fraysseix and Rosenstiehl in Brandes's formulation. Trivial cases first: at most 4 vertices or fewer than 9 edges is planar (K5 and K3,3 need 9), and `m > 3n - 6` is rejected by the Euler bound. Then per connected component two DFS passes: *orientation* orients edges away from the root and computes `lowpt`, `lowpt2` and a nesting depth `2*lowpt + [lowpt2 < height]`, ordering each vertex's outgoing edges by it; *testing* walks the tree in nesting order keeping a stack of conflict pairs (L, R) of return-edge intervals, merging constraints after each child (`add_constraints`) and trimming returns when retreating (`remove_back_edges`). A contradiction means non-planar. The embedding phase is omitted.

**Data structures.** CSR adjacency from the `Graph[List, List]` tree; one malloc'd block of ints (heights, parent edges, lowpoints, `ref` chains, a flat conflict-pair stack of four edge ids per pair). Both DFS passes are iterative with an explicit stack and per-vertex cursors, so a path of 10^6 vertices needs no call stack; phase 2 runs on preorder numbers for cache locality.

**Complexity / limits.** `O(n + m)` time and memory. Returns a definite `True`/`False`; only an allocation failure leaves the head unevaluated. The test decides planarity but constructs no embedding or Kuratowski subgraph.

- `Protected`; gives `False` for a non-graph. Edge directions are ignored.
- Algorithm (`galg_planar.c`): the linear-time left-right planarity test of de
  Fraysseix and Rosenstiehl in Brandes' formulation (the same algorithm as
  networkx's `check_planarity`), without the embedding phase. Two iterative DFS
  passes per connected component (orientation with lowpoints and nesting
  depths, then conflict-pair testing), so path-like graphs with 10^6 vertices
  need no call stack. The Euler bound rejects `m > 3n - 6` up front.
- O(n + m) time and memory.

**Attributes:** `Protected`.

## References

- U. Brandes, *The Left-Right Planarity Test*, manuscript (2009), algorithms 2-5.
- Source: [`src/graph/galg_planar.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_planar.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

The test is the linear-time left-right criterion, so it is cheap even on very large sparse graphs. It decides planarity only; it does not return an embedding or a Kuratowski subgraph.

Edge direction, anti-parallel pairs and loops are ignored: the test runs on the underlying simple undirected graph. Graphs with more than `3n - 6` edges are rejected immediately by the Euler bound.
