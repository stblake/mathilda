### Worked examples

```mathematica
In[1]:= IntervalUnion[Interval[{1, 3}], Interval[{2, 5}]]  (* overlapping ranges merge *)
```

```mathematica
In[1]:= IntervalUnion[Interval[{1, 2}], Interval[{5, 6}]]  (* disjoint ranges stay separate *)
```

```mathematica
In[1]:= Attributes[IntervalUnion]  (* Flat + Orderless -- the algebra of set union *)
```

### Notes

`IntervalUnion[i1, i2, ...]` gives the interval representing the set union of its arguments.
Overlapping or touching ranges merge into one (`{1,3} ∪ {2,5}` becomes `Interval[{1, 5}]`),
while disjoint ranges are kept as a multi-range interval (`Interval[{1, 2}, {5, 6}]`).

It is `Flat` and `Orderless`, matching set union's associativity and commutativity: nested
unions splice together and argument order does not matter, so the canonicaliser sees every
range at once. An empty `IntervalUnion[]` is the empty interval `Interval[]`. Its companions
are `IntervalIntersection` (which gives `Interval[]` for disjoint inputs) and
`IntervalMemberQ`.
