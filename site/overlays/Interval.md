### Worked examples

```mathematica
In[1]:= Interval[{1, 2}] + Interval[{3, 4}]  (* arithmetic threads to a rigorous enclosure *)
```

```mathematica
In[1]:= Interval[{2, 3}]*Interval[{-1, 1}]  (* the product spans the sign change *)
```

```mathematica
In[1]:= Sqrt[Interval[{1, 4}]]  (* monotone functions thread; exact endpoints stay exact *)
```

```mathematica
In[1]:= Sin[Interval[{0, Pi}]]  (* the crest at Pi/2 pins the upper endpoint to exact 1 *)
```

```mathematica
In[1]:= IntervalMemberQ[Interval[{1, 5}], 3]  (* membership test *)
```

### Notes

`Interval[{min, max}]` represents the inclusive range of reals between `min` and `max`, and
`Interval[{a1,b1}, {a2,b2}, ...]` the union of several ranges; the input is canonicalised
(pairs ordered, ranges sorted and overlapping ones merged). It is a *set*, not a number — it
carries `Protected` only, not `Listable` or `NumericFunction`.

Arithmetic and elementary functions thread through to produce **rigorous enclosures**:
endpoints are computed by the ordinary evaluator, so exact endpoints stay exact (that is why
`Sqrt[Interval[{1,4}]]` is `Interval[{1,2}]` exactly), and only inexact endpoints are nudged
one ULP outward so containment is never violated. Non-monotone functions pin an interior
extremum — `Sin` over `{0, Pi}` reaches its crest at `Pi/2`, giving an upper endpoint of
exactly `1`. Use `IntervalUnion`, `IntervalIntersection` and `IntervalMemberQ` for set
operations.
