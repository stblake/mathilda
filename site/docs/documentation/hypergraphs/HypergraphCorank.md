# HypergraphCorank

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphCorank[h] gives the smallest hyperedge arity of h (0 if h has no hyperedges).`**

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

The smallest arity

```mathematica
In[7]:= HypergraphCorank[Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]]
Out[7]= 1
```

A 3-uniform hypergraph: corank = rank

```mathematica
In[8]:= HypergraphCorank[Hypergraph[{{1, 2, 3}, {2, 3, 4}}]]
Out[8]= 3
```

## Implementation notes

**Algorithm.** `builtin_hypergraph_corank` is `rank_impl(res, want_max = 0)`, the
mirror of `HypergraphRank`: over the raw hyperedge CSR `eoff` it keeps the
smallest arity `eoff[j+1] - eoff[j]` (the `Length`, counting a repeated vertex).
With no hyperedges it is `0`.

**Data structures.** The borrowed `HypView` from the validated-hypergraph memo;
only the raw hyperedge offsets `eoff` are read. The result is a single `Integer`.

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

**See also:** [HyperedgeSizes](../../hypergraphs/HyperedgeSizes/), [HypergraphRank](../../hypergraphs/HypergraphRank/), [UniformHypergraphQ](../../hypergraphs/UniformHypergraphQ/), [Length](../../structural-manipulation/Length/)

- Source: [`src/graph/hyp_ops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_ops.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`HypergraphCorank[h]` is the smallest hyperedge arity — `Min[HyperedgeSizes[h]]` —
the companion of `HypergraphRank` (the largest). Arity is the `Length` as written,
counting a repeated vertex; with no hyperedges the corank is `0`.

When rank and corank coincide the hypergraph is uniform, and their common value is
the `k` of `k`-uniformity. The computation reads only the hyperedge offsets, so it
is `O(m)`.
