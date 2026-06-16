### Worked examples

```mathematica
In[1]:= ExponentialMovingAverage[{1, 2, 3, 4, 5}, 1/2]
Out[1]= {1, 3/2, 9/4, 25/8, 65/16}
```

```mathematica
In[1]:= ExponentialMovingAverage[{1, 2, 3, 4, 5}, 0.2]
Out[1]= {1.0, 1.2, 1.56, 2.048, 2.6384}
```

```mathematica
In[1]:= N[ExponentialMovingAverage[{10, 20, 30, 40}, 1/3], 20]
Out[1]= {10.0, 13.3333333333333333334, 18.888888888888888889, 25.9259259259259259259}
```

```mathematica
In[1]:= ExponentialMovingAverage[{a, b, c}, alpha]
Out[1]= {a, a + alpha (-a + b), a + alpha (-a + b) + alpha (-a - alpha (-a + b) + c)}
```

### Notes

`ExponentialMovingAverage[list, alpha]` applies the recurrence
`y_1 = x_1`, `y_{i+1} = y_i + alpha (x_{i+1} - y_i)`. Exact rationals stay
exact, inexact data evaluates at the requested precision, and a symbolic
smoothing constant `alpha` produces the unrolled closed form term by term.
