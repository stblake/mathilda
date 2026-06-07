### Worked examples

```mathematica
In[1]:= RepeatedTiming[Sum[i, {i, 100}]]
Out[1]= {4.53806e-05, 5050}
```

### Notes

`RepeatedTiming` evaluates the expression many times and returns `{averageSeconds, result}`, giving a steadier estimate than `Timing` for fast operations. The timing value varies between runs; only the result element is reproducible.
