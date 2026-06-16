### Worked examples

```mathematica
In[1]:= SetPrecision[1.5, 30]
Out[1]= 1.5

In[2]:= Precision[SetPrecision[1.5, 30]]
Out[2]= 30.103
```

Applied to an exact constant, `SetPrecision` produces a high-precision numeric
value — here 50 correct digits of `Pi`:

```mathematica
In[1]:= SetPrecision[Pi, 50]
Out[1]= 3.1415926535897932384626433832795028841971693993751
```

It also re-rounds an exact rational to a chosen precision; note the trailing
digit shows where the finite expansion is cut off:

```mathematica
In[1]:= SetPrecision[1/3, 40]
Out[1]= 0.33333333333333333333333333333333333333332
```

### Notes

`SetPrecision[x, n]` returns a value equal to `x` but carrying `n` digits of precision; the printed form may look unchanged while the internal precision is raised (confirm with `Precision`). Padding extra digits onto a machine-precision number introduces meaningless trailing bits, so only widen precision when the original value genuinely has them.
