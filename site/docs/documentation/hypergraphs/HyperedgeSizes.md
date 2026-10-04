# HyperedgeSizes

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HyperedgeSizes[h] gives the arity (Length) of each hyperedge of h, in EdgeList order.`**

## Examples (10)

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

### Applications (4)

```mathematica
In[7]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
Out[7]= Hypergraph[<7 vertices, 4 hyperedges>]
```

The arity of each hyperedge, in EdgeList order

```mathematica
In[8]:= HyperedgeSizes[h]
Out[8]= {3, 2, 3, 1}
```

Their total is the total incidence of h

```mathematica
In[9]:= Total[HyperedgeSizes[h]]
Out[9]= 9
```

Arity counts repeats; an empty hyperedge is 0

```mathematica
In[10]:= HyperedgeSizes[Hypergraph[{{1, 1, 2}, {}}]]
Out[10]= {3, 0}
```

## Implementation notes

**Algorithm.** `builtin_hyperedge_sizes` takes a `HypView` of the hypergraph
(incidence not needed) and reads the **raw** hyperedge CSR `eoff`: the arity of
hyperedge `j` is `eoff[j+1] - eoff[j]`, its `Length` exactly as written, so a
repeated vertex counts and an empty hyperedge is `0`. The arities are returned in
`EdgeList` order.

**Data structures.** The borrowed integer view from the validated-hypergraph memo
(no hashing, no copy of the tree). The result is assembled through `hyp_int_list`,
which hands back a packed `int64` buffer above the packing threshold and a plain
`List` of integers otherwise.

**Complexity / limits.** `O(m)` for `m` hyperedges — the offset differences are a
single pass, independent of the total incidence. Left unevaluated on a
non-hypergraph.

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

**See also:** [HypergraphRank](../../hypergraphs/HypergraphRank/), [HypergraphCorank](../../hypergraphs/HypergraphCorank/), [UniformHypergraphQ](../../hypergraphs/UniformHypergraphQ/), [Length](../../structural-manipulation/Length/)

- Source: [`src/graph/hyp_ops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_ops.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`HyperedgeSizes` reports each hyperedge's `Length` **as written** — a repeated
vertex counts, and an empty hyperedge is `0` — because arity is an ordered,
multiset notion, not the distinct-vertex count. The sum of the sizes is the total
incidence `Σ|e|`, which is also the edge count of the star expansion.

The result is a packed `int64` list above the packing threshold, so it is cheap to
feed into `Max`, `Min`, `Total`, `Tally` and the rest of the numeric pipeline —
indeed `HypergraphRank` and `HypergraphCorank` are just its maximum and minimum.
