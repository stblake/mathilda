# IntervalMemberQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`IntervalMemberQ[interval, x] gives True if x lies within interval, and False otherwise. IntervalMemberQ[interval, other] tests whether the interval other is wholly contained in interval.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

2 lies in [1, 3]

```mathematica
In[1]:= IntervalMemberQ[Interval[{1, 3}], 2]
Out[1]= True
```

```mathematica
In[2]:= IntervalMemberQ[Interval[{1, 3}], 5]
Out[2]= False
```

Subset test: [2,3] inside [0,10]

```mathematica
In[3]:= IntervalMemberQ[Interval[{0, 10}], Interval[{2, 3}]]
Out[3]= True
```

Pi ~ 3.14159 is just outside

```mathematica
In[4]:= IntervalMemberQ[Interval[{0, 3}], Pi]
Out[4]= False
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

- Source: [`src/interval.c`](https://github.com/stblake/mathilda/blob/main/src/interval.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_interval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_interval.c)

## Notes & additional examples

### Notes

`IntervalMemberQ[interval, x]` tests membership. With a scalar `x` it asks whether `x`
lies in one of the interval's `[lo, hi]` pieces; with a second `Interval` it asks
whether that interval is wholly contained in the first.

Comparisons are exact when the endpoints are exact, so `IntervalMemberQ[Interval[{0,
3}], Pi]` is `False` because `Pi` exceeds `3`. If an endpoint ordering cannot be
decided (a symbolic bound), the call is left unevaluated rather than guessing.
`IntervalMemberQ` is `Protected`.
