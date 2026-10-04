# HypergraphRank

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphRank[h] gives the largest hyperedge arity of h (0 if h has no hyperedges).`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= HyperedgeSizes[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]
Out[1]= {3, 2, 3, 1}

In[2]:= {HypergraphRank[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]], HypergraphCorank[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]}
Out[2]= {3, 1}

In[3]:= UniformHypergraphQ[Hypergraph[{{1,2,3},{2,3,4}}], 3]
Out[3]= True

In[4]:= HyperedgeSizes[Hypergraph[{{1,1,2},{}}]]
Out[4]= {3, 0}

In[5]:= {HypergraphRank[Hypergraph[{1,2},{}]], HypergraphCorank[Hypergraph[{1,2},{}]], UniformHypergraphQ[Hypergraph[{1,2},{}]]}
Out[5]= {0, 0, True}

In[6]:= HyperedgeSizes[{{1, 2}}]
Out[6]= HyperedgeSizes[{{1, 2}}]
```

### Applications (2)

The largest arity

```mathematica
In[7]:= HypergraphRank[Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]]
Out[7]= 3
```

No hyperedges -> 0

```mathematica
In[8]:= HypergraphRank[Hypergraph[{1, 2}, {}]]
Out[8]= 0
```

## Implementation notes

**Algorithm.** `builtin_hypergraph_rank` is `rank_impl(res, want_max = 1)`: over
the raw hyperedge CSR `eoff` it tracks the largest arity `eoff[j+1] - eoff[j]`,
counting a repeated vertex (the `Length`, not the distinct-vertex count). With no
hyperedges it is `0`. `HypergraphCorank` is the same routine with
`want_max = 0`.

**Data structures.** The borrowed `HypView` from the validated-hypergraph memo;
the raw hyperedge offsets `eoff` are all that is read. The result is a single
`Integer`.

**Complexity / limits.** `O(m)` — one pass over the hyperedge offsets. Left
unevaluated on a non-hypergraph.

- Arity is the hyperedge's `Length` as written, counting a repeated vertex and
  giving 0 for an empty hyperedge.
- With no hyperedges, `HypergraphRank` and `HypergraphCorank` are 0 and
  `UniformHypergraphQ` is `True`.
- `UniformHypergraphQ` gives `False` for a non-hypergraph; the other three are
  left unevaluated on one.
- Benchmark (experiment 96): `HyperedgeSizes` on 10⁵ hyperedges 0.024 ms warm /
  0.025 ms cold (Mathematica 5.6 ms, xgi 9.6 ms).

**Attributes:** `Protected`.

## References

**See also:** [HyperedgeSizes](../../hypergraphs/HyperedgeSizes/), [HypergraphCorank](../../hypergraphs/HypergraphCorank/), [UniformHypergraphQ](../../hypergraphs/UniformHypergraphQ/), [Length](../../structural-manipulation/Length/)

- Source: [`src/graph/hyp_ops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_ops.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`HypergraphRank[h]` is the largest hyperedge arity — `Max[HyperedgeSizes[h]]` — and
`HypergraphCorank[h]` the smallest. Arity is the `Length` as written, so a repeated
vertex contributes to the rank. A hypergraph with no hyperedges has rank `0`.

A hypergraph is `k`-uniform exactly when its rank equals its corank equals `k`;
`UniformHypergraphQ` tests that directly. Like the other arity heads this reads
only the hyperedge offsets, so it is `O(m)` and independent of how large the
hyperedges are.
