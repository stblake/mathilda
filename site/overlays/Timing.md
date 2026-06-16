### Worked examples

```mathematica
In[1]:= Timing[Sum[i, {i, 1000}]]
Out[1]= {0.000244, 500500}
```

The timing is non-deterministic, so extract the reproducible result with `Part`:

```mathematica
In[1]:= Timing[Sum[i, {i, 1, 1000000}]][[2]]
Out[1]= 500000500000
```

```mathematica
In[1]:= Timing[D[Tan[x]^10, x]][[2]]
Out[1]= 10 Sec[x]^2 Tan[x]^9
```

### Notes

`Timing` returns `{seconds, result}`. The first element is the CPU time spent and varies between runs and machines; only the second element (the computed result) is reproducible.
