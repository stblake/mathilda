### Worked examples

```mathematica
In[1]:= ToRadicals[Root[#^2 - 2 &, 1]]
Out[1]= -Sqrt[2]
```

The smaller root of `x^2 + x - 1` is the negative reciprocal of the golden ratio,
recovered exactly in radicals:

```mathematica
In[1]:= ToRadicals[Root[#^2 + # - 1 &, 1]]
Out[1]= 1/2 (-1 - Sqrt[5])
```

Cardano's formula appears automatically for the real root of a depressed cubic:

```mathematica
In[1]:= ToRadicals[Root[1 + #1 + #1^3 &, 1]]
Out[1]= -1/3 ((1/2 (27 + 3 Sqrt[93]))^(1/3) - 3/(1/2 (27 + 3 Sqrt[93]))^(1/3))
```

It threads through the implicit `Root` objects produced by `Solve`, here giving
the three cube roots of two:

```mathematica
In[1]:= ToRadicals[Solve[x^3 - 2 == 0, x]]
Out[1]= {{x -> 2^(1/3)}, {x -> -(-1)^(1/3) 2^(1/3)}, {x -> (-1)^(2/3) 2^(1/3)}}
```

### Notes

`ToRadicals[expr]` rewrites `Root` objects in `expr` using radicals. It always
succeeds when the underlying polynomial has degree at most four (and for
binomial `Root[a #^n + b &, k]` of any degree); degree-five-and-higher Root
objects are returned unchanged. It threads automatically over lists, equations,
and the results of `Solve`.
