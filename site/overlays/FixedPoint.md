### Worked examples

```mathematica
In[1]:= FixedPoint[Floor[#/2]&, 100]
Out[1]= 0

In[2]:= FixedPoint[Function[x, (x + 2/x)/2], 1.0]
Out[2]= 1.41421
```

### Notes

`FixedPoint[f, expr]` applies `f` repeatedly until the result stops changing. The second example is Newton's iteration converging to `Sqrt[2]`; convergence relies on exact equality of successive results, so use machine reals (or a `SameTest`) for numeric iterations.
