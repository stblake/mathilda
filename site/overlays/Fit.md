### Worked examples

```mathematica
In[1]:= Fit[{1, 4, 9, 16}, {1, x, x^2}, x]
Out[1]= 0.0 + 0.0 x + 1.0 x^2

In[2]:= Fit[{{1, 1}, {2, 4}, {3, 9}}, {1, x, x^2}, x]
Out[2]= 0.0 + 0.0 x + 1.0 x^2

In[3]:= Fit[{1, 2, 1.3, 3.75, 2.25}, {1, x}, x]
Out[3]= 0.785 + 0.425 x
```

### Notes

`Fit[data, {f1, ..., fn}, x]` returns the least-squares linear combination
`a1 f1 + ... + an fn` of the basis functions. Plain `{v1, v2, ...}` data is
taken at abscissae `1, 2, ...`, while `{{x, v}, ...}` pairs supply explicit
abscissae; the perfect square data above recovers `x^2` exactly.
