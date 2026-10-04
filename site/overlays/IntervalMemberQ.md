### Worked examples

```mathematica
In[1]:= IntervalMemberQ[Interval[{1, 3}], 2]  (* 2 lies in [1, 3] *)
```

```mathematica
In[1]:= IntervalMemberQ[Interval[{1, 3}], 5]
```

```mathematica
In[1]:= IntervalMemberQ[Interval[{0, 10}], Interval[{2, 3}]]  (* subset test: [2,3] inside [0,10] *)
```

```mathematica
In[1]:= IntervalMemberQ[Interval[{0, 3}], Pi]  (* Pi ~ 3.14159 is just outside *)
```

### Notes

`IntervalMemberQ[interval, x]` tests membership. With a scalar `x` it asks whether `x`
lies in one of the interval's `[lo, hi]` pieces; with a second `Interval` it asks
whether that interval is wholly contained in the first.

Comparisons are exact when the endpoints are exact, so `IntervalMemberQ[Interval[{0,
3}], Pi]` is `False` because `Pi` exceeds `3`. If an endpoint ordering cannot be
decided (a symbolic bound), the call is left unevaluated rather than guessing.
`IntervalMemberQ` is `Protected`.
