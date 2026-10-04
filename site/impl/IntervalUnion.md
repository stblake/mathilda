---
source: src/interval.c
---
**Algorithm.** `IntervalUnion[i1, i2, ...]` gives the interval representing the set union of
its argument intervals. `builtin_intervalunion` (`src/interval.c`) collects the ranges of
every argument interval and hands them to the same canonicaliser `Interval[...]` uses:
endpoint pairs are ordered, ranges are sorted by lower endpoint, and overlapping or touching
ranges are merged. The result is one `Interval[...]` with the minimal list of disjoint
ranges covering the union — a single range when the inputs overlap, several when they do
not.

**Data structures.** Ordinary `Expr` trees. `IntervalUnion` is registered
`Flat | Orderless | Protected`: `Flat` lets nested `IntervalUnion[...]` calls splice
together and `Orderless` lets the evaluator present the arguments in canonical order, both
of which match the algebra of set union (associative and commutative) and let the
canonicaliser see every range in one pass. Endpoint ordering uses the same exact
(200-bit MPFR) comparison as the rest of the interval module.

**Complexity / limits.** `O(k log k)` in the total number `k` of ranges, dominated by the
sort. An empty call `IntervalUnion[]` is the empty interval `Interval[]`. The companions are
`IntervalIntersection` (also `Flat | Orderless`, giving `Interval[]` when the inputs are
disjoint) and `IntervalMemberQ`.
