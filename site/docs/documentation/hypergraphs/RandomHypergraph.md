# RandomHypergraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RandomHypergraph[{n, m}, k] gives a k-uniform hypergraph on 1..n with m independent uniformly random k-subsets as hyperedges; RandomHypergraph[{n, m}, k, c] gives c of them. RandomHypergraph[{n, {e, a}}] and [{n, {{e1, a1}, ...}}] draw e_i hyperedges of arity a_i with vertices chosen with replacement from 1..n. Honors SeedRandom.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= SeedRandom[1]; EdgeList[RandomHypergraph[{10, 4}, 3]]
Out[1]= {{2, 7, 9}, {2, 6, 10}, {1, 5, 8}, {2, 4, 9}}

In[2]:= SeedRandom[1]; RandomHypergraph[{10, 4}, 3]
Out[2]= Hypergraph[<10 vertices, 4 hyperedges>]

In[3]:= SeedRandom[3]; Map[EdgeList, RandomHypergraph[{5, 2}, 2, 3]]
Out[3]= {{{1, 4}, {4, 5}}, {{1, 3}, {1, 4}}, {{1, 3}, {1, 3}}}

In[4]:= SeedRandom[1]; InputForm[RandomHypergraph[{10, {3, 3}}]]
Out[4]= Hypergraph[{9, 8, 2, 6, 10, 1}, {{9, 8, 2}, {8, 2, 6}, {10, 6, 1}}]

In[5]:= SeedRandom[2]; EdgeList[RandomHypergraph[{6, {{2, 3}, {1, 1}}}]]
Out[5]= {{5, 4, 4}, {2, 3, 5}, {3}}

In[6]:= RandomHypergraph[{3, 2}, 4]
Out[6]= RandomHypergraph[{3, 2}, 4]
```

## Implementation notes

- k-uniform form: `m` independent, uniformly random k-subsets of `1..n`
  (listed increasing); a hyperedge may repeat. The vertex list is all of
  `1..n`. Requires `0 ≤ k ≤ n`, else unevaluated. Floyd's sampling, `O(k)` per
  hyperedge.
- FR form: vertices drawn uniformly **with replacement** from `1..n` (so a
  vertex can repeat inside a hyperedge, as in the FR output `{8, 1, 1}`); the
  vertex list is the vertices that occur, in first-appearance order.
  *Deviation:* returns a `Hypergraph` where the FR function returns the bare
  edge List (`EdgeList` recovers it).
- Both forms draw from the user-visible random stream, so `SeedRandom`
  reproduces them.
- Benchmark (experiment 96), FR form with 10⁵ triples: 10.6 ms (FR function
  1.6 ms, xgi 438 ms) — the one loss in the experiment. The ResourceFunction
  returns a single packed `RandomInteger` array, while Mathilda returns a
  validated `Hypergraph` whose hyperedges are 10⁵ boxed `List` nodes
  (Mathilda's representation invariant keeps packed buffers out of expression
  trees); allocating those nodes alone exceeds 1.6 ms. For scale, Mathilda's
  own `RandomInteger[{1, 20000}, {100000, 3}]` takes 16 ms.

**Attributes:** `Protected`.

## References

**See also:** [Hypergraph](../../hypergraphs/Hypergraph/), [EdgeList](../../graphs/EdgeList/), [SeedRandom](../../random-number-generation/SeedRandom/), [RandomInteger](../../random-number-generation/RandomInteger/), [List](../../other-advanced/List/)

- Source: [`src/graph/hyp_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_init.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)
