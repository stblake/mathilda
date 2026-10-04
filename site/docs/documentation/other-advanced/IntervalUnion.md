# IntervalUnion

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`IntervalUnion[i1, i2, ...] gives the interval representing the union of the intervals ij.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

Overlapping ranges merge

```mathematica
In[1]:= IntervalUnion[Interval[{1, 3}], Interval[{2, 5}]]
Out[1]= Interval[{1, 5}]
```

Disjoint ranges stay separate

```mathematica
In[2]:= IntervalUnion[Interval[{1, 2}], Interval[{5, 6}]]
Out[2]= Interval[{1, 2}, {5, 6}]
```

Flat + Orderless -- the algebra of set union

```mathematica
In[3]:= Attributes[IntervalUnion]
Out[3]= {Flat, Orderless, Protected}
```

## Implementation notes

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

**Attributes:** `Flat`, `Orderless`, `Protected`.

## References

- Source: [`src/interval.c`](https://github.com/stblake/mathilda/blob/main/src/interval.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_interval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_interval.c)

## Notes & additional examples

### Notes

`IntervalUnion[i1, i2, ...]` gives the interval representing the set union of its arguments.
Overlapping or touching ranges merge into one (`{1,3} ∪ {2,5}` becomes `Interval[{1, 5}]`),
while disjoint ranges are kept as a multi-range interval (`Interval[{1, 2}, {5, 6}]`).

It is `Flat` and `Orderless`, matching set union's associativity and commutativity: nested
unions splice together and argument order does not matter, so the canonicaliser sees every
range at once. An empty `IntervalUnion[]` is the empty interval `Interval[]`. Its companions
are `IntervalIntersection` (which gives `Interval[]` for disjoint inputs) and
`IntervalMemberQ`.
