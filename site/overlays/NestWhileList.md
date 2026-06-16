### Worked examples

```mathematica
In[1]:= NestWhileList[#/2 &, 256, EvenQ]
Out[1]= {256, 128, 64, 32, 16, 8, 4, 2, 1}
```

The full Collatz (3n+1) trajectory starting from 27 — it climbs as high as 9232
before crashing to 1, taking 112 entries to get there:

```mathematica
In[1]:= Length[NestWhileList[If[EvenQ[#], #/2, 3 # + 1] &, 27, # > 1 &]]
Out[1]= 112
```

With `UnsameQ` and a history of `2`, `NestWhileList` is `FixedPointList`. Here it
records every Newton convergent of `Sqrt[2]` up to the machine fixed point:

```mathematica
In[1]:= NestWhileList[(# + 2/#)/2 &, 1.0, UnsameQ, 2]
Out[1]= {1.0, 1.5, 1.41667, 1.41422, 1.41421, 1.41421, 1.41421}
```

### Notes

`NestWhileList[f, expr, test]` is the value-collecting companion of `NestWhile`:
it returns `{expr, f[expr], f[f[expr]], ...}`, continuing while `test` applied to
the most recent result(s) yields `True`. The optional fourth argument sets how
many recent results `test` receives (an integer, `All`, or `{mmin, mmax}`); a
fifth caps applications; a sixth appends `n` extra iterates or drops the last `n`.
With `UnsameQ` and a history of `2` it is equivalent to `FixedPointList`.
