# RandomHypergraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RandomHypergraph[{n, m}, k] gives a k-uniform hypergraph on 1..n with m independent uniformly random k-subsets as hyperedges; RandomHypergraph[{n, m}, k, c] gives c of them. RandomHypergraph[{n, {e, a}}] and [{n, {{e1, a1}, ...}}] draw e_i hyperedges of arity a_i with vertices chosen with replacement from 1..n. Honors SeedRandom.`**

## Examples (9)

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

### Applications (3)

3 random 2-subsets of 1..6

```mathematica
In[7]:= SeedRandom[42]; EdgeList[RandomHypergraph[{6, 3}, 2]]
Out[7]= {{2, 5}, {5, 6}, {4, 6}}
```

Every hyperedge has arity 2, uniform by construction

```mathematica
In[8]:= SeedRandom[42]; UniformHypergraphQ[RandomHypergraph[{6, 3}, 2], 2]
Out[8]= True
```

The FR form: 4 hyperedges of arity 3

```mathematica
In[9]:= SeedRandom[7]; HyperedgeSizes[RandomHypergraph[{10, {4, 3}}]]
Out[9]= {3, 3, 3, 3}
```

## Implementation notes

**Algorithm.** `builtin_random_hypergraph` has two families.
`RandomHypergraph[{n, m}, k]` (and the `..., c` repeat form) builds a `k`-uniform
hypergraph on `1..n`: each of the `m` hyperedges is an independent uniformly
random `k`-subset drawn by **Floyd's algorithm** — for `t = n-k .. n-1`, pick
`r` in `[0, t]` and take `r` unless it is already in the subset, in which case take
`t` — then sorted (insertion sort for `k <= 16`, else `qsort`). The
Function-Repository form `RandomHypergraph[{n, {e, a}}]` (or a list of `{e_i, a_i}`
pairs) draws vertices uniformly **with replacement**, so a vertex may repeat inside
a hyperedge; the vertex set is the vertices that actually occur, in first-appearance
order (the constructor derives them).

**Data structures.** A reusable `stamp` membership array (compared to the hyperedge
index `j`, so no clearing between draws) and a `buf` of chosen indices; the shared
`List` head and the `1..n` integer nodes are built once and shared by `expr_copy`.
All randomness is `random_uniform_01`, the user-visible stream, so both forms
reproduce under `SeedRandom`.

**Complexity / limits.** `O(k)` per `k`-uniform hyperedge (no rejection). Requires
`0 <= k <= n`, else unevaluated. Counts are capped at `INT32_MAX`. Where the FR
function returns the bare edge List, this returns a validated `Hypergraph`
(`EdgeList` recovers the List).

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

- J. Bentley and B. Floyd, *Programming pearls: a sample of brilliance*, Communications of the ACM **30**(9) (1987) 754-757 (Floyd's sampling).
- Source: [`src/graph/hyp_random.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_random.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_graphplot.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphplot.c)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`RandomHypergraph[{n, m}, k]` draws `m` independent, uniformly random `k`-subsets
of `1..n` as hyperedges (so a hyperedge may repeat); the vertex list is all of
`1..n`, and the result is `k`-uniform by construction. `RandomHypergraph[{n, m}, k, c]`
gives `c` such hypergraphs.

`RandomHypergraph[{n, {e, a}}]` is the Function Repository form: `e` hyperedges of
arity `a` whose vertices are drawn **with replacement**, so a vertex may repeat
inside a hyperedge; `{n, {{e1, a1}, ...}}` mixes several arities. All forms draw
from the user-visible random stream, so `SeedRandom` reproduces them exactly — as
above, where the same seed makes the first two examples agree. Where the FR
function returns the bare edge List, this returns a validated `Hypergraph`;
`EdgeList` recovers the List.
