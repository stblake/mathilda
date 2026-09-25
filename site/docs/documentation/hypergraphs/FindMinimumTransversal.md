# FindMinimumTransversal

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindMinimumTransversal[h] gives a smallest set of vertices meeting every hyperedge of h (a minimum hitting set). Left unevaluated if some hyperedge is empty or the exact search exceeds its budget.`**

## Examples (6)

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

## Implementation notes

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

- Source: [`src/graph/hyp_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_init.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)
