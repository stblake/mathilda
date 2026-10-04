---
references:
  - "K. Murakami and T. Uno, *Efficient algorithms for dualizing large-scale hypergraphs*, Discrete Applied Mathematics **170** (2014) 83-94 (the MMCS algorithm)."
  - "C. Berge, *Hypergraphs: Combinatorics of Finite Sets*, North-Holland Mathematical Library 45 (Elsevier, 1989), ch. 2 (transversals)."
source: src/graph/hyp_transversal.c
---
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
