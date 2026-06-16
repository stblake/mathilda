### Worked examples

```mathematica
In[1]:= Variance[{1, 2, 3, 4, 5}]
Out[1]= 5/2
```

`Variance` is the *unbiased* (sample) estimator, dividing the sum of squared deviations by `n - 1`. On exact rational data the answer stays exact:

```mathematica
In[1]:= Variance[{2, 4, 4, 4, 5, 5, 7, 9}]
Out[1]= 32/7
```

Feeding the same integers as 40-digit reals propagates the precision all the way through the deviation sum:

```mathematica
In[1]:= Variance[N[{1, 1, 2, 3, 5, 8, 13}, 40]]
Out[1]= 19.571428571428571428571428571428571428568
```

### Notes

`Variance[data]` gives the unbiased variance estimate (Bessel-corrected, `1/(n-1)` normalization) of the elements in `data`. Exact inputs yield exact rational results; arbitrary-precision inputs carry their precision through the computation.
