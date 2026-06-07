### Worked examples

```mathematica
In[1]:= ToRadicals[Root[#^2 - 2 &, 1]]
Out[1]= -Sqrt[2]

In[2]:= ToRadicals[Root[1 + #1 + #1^3 &, 1]]
Out[2]= -1/3 ((1/2 (27 + 3 Sqrt[93]))^(1/3) - 3/(1/2 (27 + 3 Sqrt[93]))^(1/3))
```

### Notes

`ToRadicals[expr]` rewrites `Root` objects in `expr` using radicals. It always
succeeds when the underlying polynomial has degree at most four (and for
binomial `Root[a #^n + b &, k]` of any degree); degree-five-and-higher Root
objects are returned unchanged. It threads automatically over lists, equations,
and the results of `Solve`.
