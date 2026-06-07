### Worked examples

```mathematica
In[1]:= SetAccuracy[1.5, 30]
Out[1]= 1.5

In[2]:= Accuracy[SetAccuracy[1.5, 30]]
Out[2]= 30.2279
```

### Notes

`SetAccuracy[x, n]` returns a value equal to `x` with `n` digits of accuracy (digits past the decimal point); use `Accuracy` to confirm, since the printed form often looks unchanged. It is the absolute-magnitude counterpart to `SetPrecision`.
