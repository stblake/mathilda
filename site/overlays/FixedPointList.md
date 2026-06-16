### Worked examples

```mathematica
In[1]:= FixedPointList[Floor[#/2] &, 100]
Out[1]= {100, 50, 25, 12, 6, 3, 1, 0, 0}
```

```mathematica
In[1]:= FixedPointList[(# + 2/#)/2 &, 1.0]
Out[1]= {1.0, 1.5, 1.41667, 1.41422, 1.41421, 1.41421, 1.41421}
```

```mathematica
In[1]:= Length[FixedPointList[If[EvenQ[#], #/2, 3 # + 1] &, 27, SameTest -> (#2 == 1 &)]]
Out[1]= 112
```

### Notes

`FixedPointList[f, expr]` records the whole orbit `{expr, f[expr], ...}`
until two consecutive entries are `SameQ`, so the last two elements always
coincide. The second example shows Newton/Heron iteration converging
quadratically to `Sqrt[2]` (note the doubling of correct digits each step).
The third counts the Collatz `3n+1` trajectory length for the seed 27 by
supplying a custom `SameTest` that halts when the value reaches 1 — the
famous 111-step descent (112 entries including the seed).
