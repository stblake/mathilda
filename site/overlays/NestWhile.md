### Worked examples

```mathematica
In[1]:= NestWhile[#/2 &, 256, EvenQ]
Out[1]= 1
```

Count how many `+1` steps are taken before a predicate fails — here it lands
exactly on the boundary value:

```mathematica
In[1]:= NestWhile[# + 1 &, 1, # < 100 &]
Out[1]= 100
```

With test `UnsameQ` and history length `2`, `NestWhile` becomes `FixedPoint`.
Newton's iteration for `Sqrt[2]` converges to the machine fixed point:

```mathematica
In[1]:= NestWhile[(# + 2/#)/2 &, 1.0, UnsameQ, 2]
Out[1]= 1.41421
```

### Notes

`NestWhile[f, expr, test]` repeatedly applies `f`, stopping as soon as `test`
applied to the most recent result no longer yields `True`. The optional fourth
argument controls how many recent results are handed to `test` (an integer,
`All`, or `{mmin, mmax}`); a fifth caps the number of applications, and a sixth
shifts the returned result forward (`n`) or back (`-n`). With `UnsameQ` and a
history of `2`, `NestWhile` is equivalent to `FixedPoint`.
