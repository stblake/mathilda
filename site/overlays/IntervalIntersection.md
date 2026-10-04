### Worked examples

```mathematica
In[1]:= IntervalIntersection[Interval[{1, 5}], Interval[{3, 8}]]  (* the common range *)
```

```mathematica
In[1]:= IntervalIntersection[Interval[{1, 10}], Interval[{2, 8}], Interval[{3, 5}]]  (* folds over any number of arguments *)
```

```mathematica
In[1]:= IntervalIntersection[Interval[{1, 2}], Interval[{5, 6}]]  (* disjoint inputs give the empty interval *)
```

```mathematica
In[1]:= IntervalIntersection[Interval[{0, 4}], Interval[{1, 2}, {3, 5}]]  (* a multi-segment interval keeps each overlapping piece *)
```

### Notes

`IntervalIntersection[i1, i2, ...]` gives the interval of points common to all its
arguments. It folds left over the arguments, intersecting each segment against the
accumulator and canonicalising the surviving `[lo, hi]` pairs into sorted,
non-overlapping normal form; disjoint inputs collapse to the empty `Interval[]`. A
multi-segment `Interval[{a1, b1}, {a2, b2}, ...]` is a union of ranges, and each
surviving overlap is kept. The head is `Orderless` and `Flat`. A bare number counts
as a degenerate point interval; a non-interval, non-numeric argument leaves the call
unevaluated.
