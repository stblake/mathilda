# TransversalHypergraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`TransversalHypergraph[h] gives the hypergraph of all minimal transversals (minimal hitting sets) of h, on the vertices of h. For a plain List of hyperedges the List of minimal transversals is returned. Left unevaluated if the search exceeds its budget.`**

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

{2} hits both; so does the minimal {1, 3}

```mathematica
In[7]:= TransversalHypergraph[{{1, 2}, {2, 3}}]
Out[7]= {{2}, {1, 3}}
```

Returned as a Hypergraph on h's vertices

```mathematica
In[8]:= InputForm[TransversalHypergraph[Hypergraph[{{a, b}, {b, c}}]]]
Out[8]= Hypergraph[{a, b, c}, {{b}, {a, c}}]
```

Disjoint hyperedges: pick one vertex from each

```mathematica
In[9]:= TransversalHypergraph[{{1, 2}, {3, 4}}]
Out[9]= {{1, 3}, {1, 4}, {2, 3}, {2, 4}}
```

## Implementation notes

**Algorithm.** `builtin_transversal_hypergraph` enumerates every **minimal**
transversal (minimal hitting set) of `h` by **MMCS** (Murakami & Uno, 2014), run
iteratively with explicit frames so a huge transversal cannot overflow the C
stack. The search is depth-first over partial solutions `S` that are kept minimal
at every node — each `s ∈ S` must own a *critical* hyperedge, one met by `S` only
in `s`. Criticality is maintained in `O(deg v)` per add/remove by keeping, per
hyperedge, `|F ∩ S|` and the **sum** of `S`'s members in `F` (so when the count is
1 the sum *is* the unique member); a node with no non-critical member (`zc == 0`)
and no uncovered hyperedge emits `S`. At each node the uncovered hyperedge with the
fewest remaining candidates is branched on (`tr_choose`). Results are sorted by
size, then lexicographically by `VertexList` position. Degenerate inputs the FR
function declines are answered: no hyperedges → `{{}}`, an empty hyperedge → `{}`.

**Data structures.** The distinct-vertex CSR `soff/sv` and incidence CSR
`voff/ve`; the `TrState` bookkeeping (`cnt`, `sum`, `crit`, the swap-removed
`uncov`/`upos` uncovered set, the `incand` candidate mask, the solution stack `S`);
an `IVec` output pool with `(off, len)` records sorted by `cmp_rec`. A `Hypergraph`
input yields `Hypergraph[VertexList[h], Tr]`, a bare List yields the List `Tr`.

**Complexity / limits.** Output-sensitive — `Tr(h)` can be exponentially large.
**Correct-or-unevaluated:** the search counts nodes and returns unevaluated (never
a partial answer) past `HYP_TR_MAX_NODES = 5·10⁷` nodes or `HYP_TR_MAX_OUT = 2·10⁷`
output entries, and polls `tc_check_deadline()` every 4096 nodes so
`TimeConstrained` interrupts it. Node counts (not wall-clock) keep the verdict
machine-independent.

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

**See also:** [FindMinimumTransversal](../../hypergraphs/FindMinimumTransversal/), [Tr](../../linear-algebra/Tr/), [TimeConstrained](../../time-and-date/TimeConstrained/)

- K. Murakami and T. Uno, *Efficient algorithms for dualizing large-scale hypergraphs*, Discrete Applied Mathematics **170** (2014) 83-94 (the MMCS algorithm).
- C. Berge, *Hypergraphs: Combinatorics of Finite Sets*, North-Holland Mathematical Library 45 (Elsevier, 1989), ch. 2 (transversals).
- Source: [`src/graph/hyp_transversal.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_transversal.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

A transversal (hitting set) is a vertex set meeting every hyperedge; it is minimal
when no proper subset still does. `TransversalHypergraph[h]` returns all minimal
transversals — the transversal hypergraph `Tr(h)` of Berge — as
`Hypergraph[VertexList[h], Tr]`, or, for a bare List of hyperedges, as the List
`Tr`. Each transversal lists its vertices in `VertexList` order, and the
transversals are sorted by size then lexicographically.

The algorithm is MMCS (Murakami & Uno, 2014), the practical state of the art, run
iteratively so a large output cannot overflow the C stack. Because `Tr(h)` can be
exponentially large, the search is correct-or-unevaluated: it counts nodes and
returns unevaluated — never a partial answer — past its budget, and
`TimeConstrained` can interrupt it. Degenerate inputs the FR function declines are
answered here: no hyperedges gives `{{}}`, an empty hyperedge gives `{}`.
