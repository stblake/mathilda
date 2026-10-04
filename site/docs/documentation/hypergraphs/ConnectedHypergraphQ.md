# ConnectedHypergraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ConnectedHypergraphQ[h] gives True if h has at least one vertex and is connected. h may be a plain List of hyperedges.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= HypergraphConnectedComponents[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]
Out[1]= {{1, 2, 3, 4, 5, 6}, {7}}

In[2]:= HypergraphConnectedComponents[Hypergraph[{c,b,a},{{a,c}}]]
Out[2]= {{c, a}, {b}}

In[3]:= ConnectedHypergraphQ[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]
Out[3]= False

In[4]:= ConnectedHypergraphQ[{{1,2},{2,3}}]
Out[4]= True

In[5]:= ConnectedHypergraphQ[{}]
Out[5]= False
```

### Applications (3)

A single component

```mathematica
In[6]:= ConnectedHypergraphQ[Hypergraph[{{1, 2}, {2, 3}}]]
Out[6]= True
```

Vertex 7 is isolated

```mathematica
In[7]:= ConnectedHypergraphQ[Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]]
Out[7]= False
```

Accepts a bare list of hyperedges

```mathematica
In[8]:= ConnectedHypergraphQ[{{1, 2}, {2, 3}, {3, 4}}]
Out[8]= True
```

## Implementation notes

**Algorithm.** `builtin_connected_hypergraph_q` returns `True` iff `h` has at
least one vertex and exactly one connected component. It runs the same
`vertex_uf` union–find as `HypergraphConnectedComponents` (unioning the members of
each hyperedge's distinct set) and counts the roots `p[i] == i`; the answer is
`roots == 1`. It is the Wolfram Function Repository name and accepts a bare List
of hyperedges (`hyp_arg`); a non-hypergraph argument gives `False` (a `*Q`
predicate), and the empty List `{}` gives `False` (the FR function leaves it
unevaluated).

**Data structures.** The distinct-vertex CSR `soff/sv`; a union–find parent array
(path-halving). The result is a `True`/`False` symbol.

**Complexity / limits.** Near-linear, `O(Σ|e| · α(n))` for the union step plus
`O(n)` to count roots.

- Components are ordered by their first vertex, vertices within a component in
  VertexList order. (Graph `ConnectedComponents` instead follows Mathematica:
  largest first when undirected, strong components when directed.) Isolated
  vertices are singleton components.
- Union–find, near-linear in the total incidence.
- `ConnectedHypergraphQ` is the Wolfram Function Repository name. It accepts a
  bare List of hyperedges; the FR function leaves `{}` unevaluated, Mathilda
  gives `False`. It gives `False` for any non-hypergraph.
- Benchmark (experiment 96): `HypergraphConnectedComponents` on 10⁵ hyperedges
  2.3 ms warm and cold (Mathematica 160 ms, xgi 282 ms);
  `ConnectedHypergraphQ` 20.1 ms (FR function 261 ms, xgi 546 ms).

**Attributes:** `Protected`.

## References

**See also:** [HypergraphConnectedComponents](../../hypergraphs/HypergraphConnectedComponents/), [ConnectedComponents](../../graphs/ConnectedComponents/)

- Source: [`src/graph/hyp_ops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_ops.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`ConnectedHypergraphQ[h]` is `True` when `h` has at least one vertex and exactly
one connected component — the Boolean companion of
`HypergraphConnectedComponents`. It is the Wolfram Function Repository name and,
like that function, accepts a bare List of hyperedges as well as a `Hypergraph`
object.

It gives `False` for any non-hypergraph and for the empty List `{}` (which the FR
function leaves unevaluated). Internally it reuses the same vertex union–find and
simply checks for a single root, so it is near-linear in the total incidence.
