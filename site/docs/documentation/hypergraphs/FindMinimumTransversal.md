# FindMinimumTransversal

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindMinimumTransversal[h] gives a smallest set of vertices meeting every hyperedge of h (a minimum hitting set). Left unevaluated if some hyperedge is empty or the exact search exceeds its budget.`**

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= TransversalHypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]
Out[1]= {{1, 4, 7}, {2, 4, 7}, {3, 4, 7}, {3, 5, 7}, {3, 6, 7}}

In[2]:= InputForm[TransversalHypergraph[Hypergraph[{{a,b},{b,c}}]]]
Out[2]= Hypergraph[{a, b, c}, {{b}, {a, c}}]

In[3]:= {TransversalHypergraph[{}], TransversalHypergraph[{{}}], TransversalHypergraph[{{1},{1,2}}], TransversalHypergraph[{{1,2},{1,2}}]}
Out[3]= {{{}}, {}, {{1}}, {{1}, {2}}}

In[4]:= FindMinimumTransversal[{{1,2,3},{3,4},{4,5,6},{7}}]
Out[4]= {3, 4, 7}

In[5]:= FindMinimumTransversal[Hypergraph[{{a,b},{b,c},{c,d}}]]
Out[5]= {b, c}

In[6]:= FindMinimumTransversal[{{1,2},{}}]
Out[6]= FindMinimumTransversal[{{1, 2}, {}}]
```

### Applications (3)

One vertex, 2, hits both hyperedges

```mathematica
In[7]:= FindMinimumTransversal[{{1, 2}, {2, 3}}]
Out[7]= {2}
```

A smallest hitting set has size 2

```mathematica
In[8]:= FindMinimumTransversal[{{1, 2}, {2, 3}, {3, 4}}]
Out[8]= {2, 3}
```

The minimum transversal number

```mathematica
In[9]:= Length[FindMinimumTransversal[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]]
Out[9]= 3
```

## Implementation notes

**Algorithm.** `builtin_find_minimum_transversal` returns one minimum-cardinality
transversal (a minimum hitting set, NP-hard) by exact **branch and bound**. It
branches on the uncovered hyperedge with the fewest open vertices — include a
vertex, then exclude it for later siblings. The search is seeded with a greedy
**max-coverage upper bound** (a lazy max-heap on coverage gains, stale entries
re-pushed, so the greedy pass is `O(Σ|e| log n)`) and pruned by
`mn_lower_bound`, the larger of three bounds: a **degree bound** (fewest
high-degree open vertices to cover all uncovered hyperedges), a **greedy disjoint
packing** (pairwise-disjoint uncovered hyperedges, taken in increasing total open
degree, one vertex each), and a **fractional-packing / LP-dual** bound
(`y_F = min_{v∈F} 1/d(v)` dual-feasible, raised by one low-degree-first
dual-ascent pass, ceiling-rounded). When a packing already matches the greedy
bound the answer is proved optimal with no search at all.

**Data structures.** The distinct-vertex CSR `soff/sv` and incidence CSR
`voff/ve`; `MinState` with the covered-count array, swap-removed uncovered set,
`excl` mask, current/best solution stacks, and the lower-bound scratch
(`odeg`/`touched`/`dcount` for the degree bound, `ekey`/`eord` ordering, `mark`
stamping for the packing, `ey`/`slack` for the LP dual). The result is a sorted
`List` of vertices.

**Complexity / limits.** Exponential worst case, kept tractable by the bounds.
**Correct-or-unevaluated:** left unevaluated if any hyperedge is empty
(unhittable), past `HYP_MIN_MAX_NODES = 2·10⁷` branch-and-bound nodes, or if the
greedy upper bound exceeds 20000 (a C-stack-depth guard); it polls
`tc_check_deadline()` every 4096 nodes for `TimeConstrained`.

- A transversal (hitting set) meets every hyperedge; it is minimal when no
  proper subset is one.
- `TransversalHypergraph` is the Wolfram Function Repository name. In the
  bare-List form, non-List elements (the FR's "isolated vertices") are
  ignored. No hyperedges → `{{}}`; an empty hyperedge → `{}`. The FR function
  leaves `{}`, `{{}}`, `{{1},{1,2}}` and `{{1,2},{1,2}}` unevaluated; these are
  answered.
- Each transversal lists its vertices in VertexList order; transversals are
  sorted by size, then lexicographically by position.
- `TransversalHypergraph` algorithm: **MMCS** (Murakami & Uno, 2014), run
  iteratively, with `O(deg v)` critical-hyperedge bookkeeping (per hyperedge,
  `|F ∩ S|` and the sum of `S`'s members in `F`, which *is* the unique member
  when the count is 1).
- `FindMinimumTransversal` (NP-hard): exact branch and bound — branch on the
  uncovered hyperedge with fewest open vertices (include, then exclude for
  later siblings); prune with the best of a degree bound, a low-degree-first
  disjoint-hyperedge packing, and a fractional-packing (LP-dual) bound with one
  dual-ascent pass; seeded by the greedy max-coverage bound (lazy heap).
  Unevaluated if a hyperedge is empty.
- **Budgets (correct-or-unevaluated).** Both searches count nodes and return
  unevaluated — never a partial or unproven answer — past 5·10⁷ MMCS nodes or
  2·10⁷ output entries (`TransversalHypergraph`), or 2·10⁷ branch-and-bound
  nodes (`FindMinimumTransversal`, which also refuses a greedy bound above
  20000 to cap C-stack depth). Node counts, not wall clock, keep answers
  machine-independent; both poll `tc_check_deadline()` every 4096 nodes, so
  `TimeConstrained` interrupts them.
- Benchmark (experiment 96): `TransversalHypergraph` on 16 vertices 4.2 ms (FR
  function 5219 ms); `FindMinimumTransversal` on 80 vertices × 200 hyperedges
  245 ms against MILP baselines (Mathematica `LinearOptimization` 874 ms,
  scipy/HiGHS 784 ms) — 3.2–3.6× faster at this size, and comparable to the
  MILPs at 100 vertices.

**Attributes:** `Protected`.

## References

**See also:** [TransversalHypergraph](../../hypergraphs/TransversalHypergraph/), [Tr](../../linear-algebra/Tr/), [TimeConstrained](../../time-and-date/TimeConstrained/)

- C. Berge, *Hypergraphs: Combinatorics of Finite Sets*, North-Holland Mathematical Library 45 (Elsevier, 1989), ch. 2 (transversals / minimum hitting set).
- Source: [`src/graph/hyp_transversal.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_transversal.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`FindMinimumTransversal[h]` returns a single smallest vertex set meeting every
hyperedge — a minimum hitting set, so its `Length` is the transversal number of
`h`. Where `TransversalHypergraph` enumerates all the *minimal* transversals,
this returns just one of *minimum* cardinality.

The minimum hitting set is NP-hard, so the engine is exact branch and bound: a
greedy max-coverage upper bound seeds the search, which is then pruned by a degree
bound, a disjoint-hyperedge packing bound, and a fractional LP-dual bound — and
when a packing already matches the greedy bound the answer is proved optimal with
no search. Like `TransversalHypergraph` it is correct-or-unevaluated: it is left
unevaluated if some hyperedge is empty (nothing can hit it) or if the exact search
exceeds its node budget, and `TimeConstrained` can interrupt it. Both a
`Hypergraph` and a bare List of hyperedges are accepted.
