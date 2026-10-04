# IntervalIntersection

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`IntervalIntersection[i1, i2, ...] gives the interval representing the intersection of the intervals ij (Interval[] if they are disjoint).`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

The common range

```mathematica
In[1]:= IntervalIntersection[Interval[{1, 5}], Interval[{3, 8}]]
Out[1]= Interval[{3, 5}]
```

Folds over any number of arguments

```mathematica
In[2]:= IntervalIntersection[Interval[{1, 10}], Interval[{2, 8}], Interval[{3, 5}]]
Out[2]= Interval[{3, 5}]
```

Disjoint inputs give the empty interval

```mathematica
In[3]:= IntervalIntersection[Interval[{1, 2}], Interval[{5, 6}]]
Out[3]= Interval[]
```

A multi-segment interval keeps each overlapping piece

```mathematica
In[4]:= IntervalIntersection[Interval[{0, 4}], Interval[{1, 2}, {3, 5}]]
Out[4]= Interval[{1, 2}, {3, 4}]
```

## Implementation notes

**Algorithm.** `builtin_intervalintersection` intersects a sequence of intervals (a
bare numeric scalar is treated as a degenerate point interval). It folds left: the
accumulator starts as the pair list of the first argument (`iv_collect`, which
pushes each `[lo, hi]` pair, or `[x, x]` for a scalar/infinity, and declines on a
non-interval). For each later argument it intersects every pair of the new interval
against the whole accumulator with `iv_intersect_pair`, replacing the accumulator
with the collected overlaps. An empty accumulator yields `Interval[]` (the empty
interval — the arguments are disjoint); otherwise the surviving `[lo, hi]` pairs are
handed to `iv_canonicalize_pairs`, which sorts and merges them into the normal form.

The head is `Orderless` and `Flat`, so nested intersections flatten and argument
order does not matter before the fold runs.

**Data structures.** An `IvBuild` growable array of parallel `Expr*` lo/hi
endpoints (`ivb_init`/`ivb_push`/`ivb_free`). Endpoints are kept as exact `Expr`
leaves when exact; comparisons go through `interval_endpoint_cmp`, which can report
"undecided" for symbolic endpoints.

**Complexity / limits.** `O(p·q)` endpoint comparisons per fold step for interval
multiplicities `p`, `q`; linear in the number of arguments overall. A non-interval,
non-numeric argument makes the whole call decline (`NULL`) and stay symbolic.
Disjoint inputs correctly produce the empty `Interval[]`.

**Attributes:** `Flat`, `Orderless`, `Protected`.

## References

- Source: [`src/interval.c`](https://github.com/stblake/mathilda/blob/main/src/interval.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_interval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_interval.c)

## Notes & additional examples

### Notes

`IntervalIntersection[i1, i2, ...]` gives the interval of points common to all its
arguments. It folds left over the arguments, intersecting each segment against the
accumulator and canonicalising the surviving `[lo, hi]` pairs into sorted,
non-overlapping normal form; disjoint inputs collapse to the empty `Interval[]`. A
multi-segment `Interval[{a1, b1}, {a2, b2}, ...]` is a union of ranges, and each
surviving overlap is kept. The head is `Orderless` and `Flat`. A bare number counts
as a degenerate point interval; a non-interval, non-numeric argument leaves the call
unevaluated.
