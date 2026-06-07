### Worked examples

```mathematica
In[1]:= LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {1, 2, 2}]
Out[1]= {2/3, 1/2}

In[2]:= LeastSquares[{{1, 0}, {0, 1}, {1, 1}}, {1, 1, 3}]
Out[2]= {4/3, 4/3}
```

### Notes

`LeastSquares[m, b]` returns the `x` minimising `Norm[m . x - b]` for the
overdetermined system `m . x == b`. With exact (rational) input it gives an
exact rational answer via the pseudoinverse; pass `Method ->` or
`Tolerance ->` to control the solver and singular-value truncation.
