### Worked examples

```mathematica
In[1]:= SetPrecision[1.5, 30]
Out[1]= 1.5

In[2]:= Precision[SetPrecision[1.5, 30]]
Out[2]= 30.103
```

### Notes

`SetPrecision[x, n]` returns a value equal to `x` but carrying `n` digits of precision; the printed form may look unchanged while the internal precision is raised (confirm with `Precision`). Padding extra digits onto a machine-precision number introduces meaningless trailing bits, so only widen precision when the original value genuinely has them.
