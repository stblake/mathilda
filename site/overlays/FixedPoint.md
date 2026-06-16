### Worked examples

```mathematica
In[1]:= FixedPoint[Floor[#/2] &, 100]
Out[1]= 0
```

```mathematica
In[1]:= FixedPoint[Function[x, (x + 2/x)/2], 1.0]
Out[1]= 1.41421
```

```mathematica
In[1]:= FixedPoint[1 + 1/# &, 1.0]
Out[1]= 1.61803
```

### Notes

`FixedPoint[f, expr]` applies `f` repeatedly until the result stops changing. The second example is Newton's iteration converging to `Sqrt[2]`; the third iterates `x -> 1 + 1/x`, whose attracting fixed point is the golden ratio `(1 + Sqrt[5])/2`. Convergence relies on exact equality of successive results, so use machine reals (or a `SameTest`) for numeric iterations.
