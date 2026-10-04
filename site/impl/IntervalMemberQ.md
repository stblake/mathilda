---
source: src/interval.c
---
**Algorithm.** `builtin_intervalmemberq` takes `IntervalMemberQ[interval, x]`. The
first argument must be an `Interval[...]` (else it declines, `NULL`); the second is
either a scalar or another `Interval`. Both forms reduce to endpoint comparisons via
`interval_endpoint_cmp`, which sets an `undecided` flag when a comparison cannot be
decided (symbolic or incomparable endpoints):

- **Scalar membership** — `x` is a member iff for some pair `[lo_i, hi_i]` of the
  interval both `lo_i <= x` and `x <= hi_i`. If no pair contains it and every
  comparison was decidable, return `False`; if any comparison was undecidable and none
  matched, decline (`NULL`).
- **Subset test** (`x` is itself an `Interval`) — every pair of `x` must sit inside
  some pair of the interval; a piece that escapes returns `False`, an undecidable
  containment declines, and all-contained returns `True`.

**Data structures.** The multi-pair `Interval[{a1,b1},{a2,b2},...]` is read through the
`interval_pair_count` / `interval_pair_lo` / `interval_pair_hi` accessors shared across
the interval module; endpoints are kept as exact `Expr` when exact.

**Complexity / limits.** `O(n)` for scalar membership over `n` interval pieces, `O(m·n)`
for the subset test of an `m`-piece interval. Attributes: `Protected`. Comparisons are
exact where the endpoints are exact (so `IntervalMemberQ[Interval[{0, 3}], Pi]` is
`False`, Pi being about 3.14159); an endpoint whose order cannot be decided leaves the
call unevaluated rather than guessing.
