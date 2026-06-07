### Worked examples

```mathematica
In[1]:= Timing[Sum[i, {i, 1000}]]
Out[1]= {0.000244, 500500}
```

### Notes

`Timing` returns `{seconds, result}`. The first element is the CPU time spent and varies between runs and machines; only the second element (the computed result) is reproducible.
