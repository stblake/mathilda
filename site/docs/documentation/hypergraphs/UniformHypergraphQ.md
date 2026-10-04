# UniformHypergraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`UniformHypergraphQ[h] gives True if every hyperedge of h has the same arity. UniformHypergraphQ[h, k] gives True if h is k-uniform.`**

## Examples (9)

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

### Applications (3)

All hyperedges have arity 3

```mathematica
In[7]:= UniformHypergraphQ[Hypergraph[{{1, 2, 3}, {2, 3, 4}}]]
Out[7]= True
```

Mixed arities

```mathematica
In[8]:= UniformHypergraphQ[Hypergraph[{{1, 2, 3}, {3, 4}}]]
Out[8]= False
```

3-uniform, so not 2-uniform

```mathematica
In[9]:= UniformHypergraphQ[Hypergraph[{{1, 2, 3}, {2, 3, 4}}], 2]
Out[9]= False
```

## Implementation notes

**Algorithm.** `builtin_uniform_hypergraph_q` scans the raw hyperedge arities
`eoff[j+1] - eoff[j]` and returns `True` iff they are all equal. With a second
argument `k` (a non-negative integer), the common arity must equal `k`
(`k`-uniform). The first arity seen becomes the reference; a hypergraph with no
hyperedges is `True`. On a non-hypergraph the one-argument form returns `False`
(it is a `*Q` predicate); the `k` form with a malformed `k` is left unevaluated.

**Data structures.** The borrowed `HypView` from the validated-hypergraph memo;
only the raw hyperedge offsets `eoff` are read. The result is a `True`/`False`
symbol.

**Complexity / limits.** `O(m)`, short-circuiting on the first mismatch. Arity is
the `Length` as written, so a repeated vertex counts toward uniformity.

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

**See also:** [HyperedgeSizes](../../hypergraphs/HyperedgeSizes/), [HypergraphRank](../../hypergraphs/HypergraphRank/), [HypergraphCorank](../../hypergraphs/HypergraphCorank/), [Length](../../structural-manipulation/Length/)

- Source: [`src/graph/hyp_ops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_ops.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`UniformHypergraphQ[h]` is `True` when every hyperedge has the same arity;
`UniformHypergraphQ[h, k]` additionally pins that common arity to `k`. Arity is the
`Length` as written, so repeats count toward uniformity, and a hypergraph with no
hyperedges is vacuously uniform.

Equivalently, `h` is uniform when `HypergraphRank[h] == HypergraphCorank[h]`. As a
`*Q` predicate the one-argument form returns `False` on a non-hypergraph rather
than staying unevaluated; the `k` form with a malformed `k` is left unevaluated.
