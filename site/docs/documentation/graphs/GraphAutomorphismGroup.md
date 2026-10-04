# GraphAutomorphismGroup

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphAutomorphismGroup[g] gives the automorphism group of g as PermutationGroup[{Cycles[...], ...}] acting on vertex positions.`**

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= GraphAutomorphismGroup[CycleGraph[4]]
Out[1]= PermutationGroup[{Cycles[{{2, 4}}], Cycles[{{1, 3}}], Cycles[{{1, 2}, {3, 4}}]}]

In[2]:= GraphAutomorphismGroup[PathGraph[{1,2,3}]]
Out[2]= PermutationGroup[{Cycles[{{1, 3}}]}]

In[3]:= GraphAutomorphismGroup[Graph[{1->2,2->3}]]
Out[3]= PermutationGroup[{}]

In[4]:= GraphAutomorphismGroup[x]
Out[4]= GraphAutomorphismGroup[x]
```

### Applications (5)

Only the reversal

```mathematica
In[5]:= GraphAutomorphismGroup[PathGraph[{1, 2, 3, 4}]]
Out[5]= PermutationGroup[{Cycles[{{1, 4}, {2, 3}}]}]
```

Generators of the dihedral group of the square

```mathematica
In[6]:= GraphAutomorphismGroup[CycleGraph[4]]
Out[6]= PermutationGroup[{Cycles[{{2, 4}}], Cycles[{{1, 3}}], Cycles[{{1, 2}, {3, 4}}]}]
```

The leaves permute freely

```mathematica
In[7]:= GraphAutomorphismGroup[StarGraph[4]]
Out[7]= PermutationGroup[{Cycles[{{2, 3}}], Cycles[{{3, 4}}]}]
```

A handful of generators describes all 120 symmetries

```mathematica
In[8]:= Length[First[GraphAutomorphismGroup[PetersenGraph[]]]]
Out[8]= 4
```

Transpositions generate the symmetric group

```mathematica
In[9]:= GraphAutomorphismGroup[CompleteGraph[4]]
Out[9]= PermutationGroup[{Cycles[{{2, 3}}], Cycles[{{3, 4}}], Cycles[{{1, 4}}]}]
```

## Implementation notes

**Algorithm.** `builtin_graph_automorphism_group` reduces the graph with `gi_build` (loops into colours, parallel edges by subdivision vertices, as in `FindGraphIsomorphism`) and calls `galg_iso_automorphisms`, which walks the individualization-refinement search tree of `galg_iso.c` with nauty-style pruning: trace pruning, automorphisms detected from leaves equal to the first or best leaf and from internal nodes whose trace matches the first path, and orbit pruning on the first path. The automorphisms found generate the full group. Each generator is turned into `Cycles[{...}]` (1-based points, each cycle starting at its least point, fixed points dropped, the identity skipped) and the result is `PermutationGroup[{gens}]`.

**Data structures.** CSR relations (`GalgIsoGraph`) built from the `Graph[List, List]` tree, a flat `int` array of generators (one permutation of the reduced vertex set per row), and a `seen` byte mask used while extracting cycles. Only generators acting on the original vertices are exposed, since subdivision vertices are determined by their ends.

**Complexity / limits.** Computing a generating set is not known to be polynomial in general, but the pruned search is fast on most graphs, and only a generating set is returned (never the elements), so the output stays compact even when the group is enormous. The search is budgeted at 5e7 nodes and respects `TimeConstrained`; on exhaustion the head stays unevaluated. A graph with no symmetry yields `PermutationGroup[{}]`. The set of generators is not unique.

- `Protected`; unevaluated on a non-graph.
- The group is given by a generating set: the group is Mathematica's, but the
  generators may differ.
- Mathilda has no permutation-group functions, so the result is an inert
  expression.
- The automorphisms found during canonical labeling (leaf and internal-node
  automorphisms, with jump-back and first-path orbit pruning) generate the whole
  group. Directed and mixed graphs are supported; edge weights are ignored.

**Attributes:** `Protected`.

## References

- B. D. McKay and A. Piperno, *Practical graph isomorphism, II*, J. Symbolic Comput. **60** (2014) 94-112.
- Source: [`src/graph/galg_isoheads.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_isoheads.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

The group is returned as `PermutationGroup[{Cycles[...], ...}]` acting on the vertices numbered by their position (1-based) in the vertex list. Only a generating set is produced, and it is not unique, so do not compare groups by their generators.

A graph with no symmetry gives `PermutationGroup[{}]`. The engine is the individualization-refinement search also used by `FindGraphIsomorphism`, so large, sparse graphs are fast.
