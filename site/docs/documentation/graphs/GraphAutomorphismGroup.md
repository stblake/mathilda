# GraphAutomorphismGroup

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphAutomorphismGroup[g] gives the automorphism group of g as PermutationGroup[{Cycles[...], ...}] acting on vertex positions.`**

## Examples (4)

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

## Implementation notes

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

- Source: [`src/graph/galg_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
