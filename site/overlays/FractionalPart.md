---
status: Stable
---
### Worked examples

```mathematica
In[1]:= FractionalPart[2.7]
Out[1]= 0.7
```

```mathematica
In[1]:= FractionalPart[-2.7]
Out[1]= -0.7
```

```mathematica
In[1]:= FractionalPart[7/2]
Out[1]= 1/2
```

```mathematica
In[1]:= IntegerPart[-2.7] + FractionalPart[-2.7]
Out[1]= -2.7
```

### Notes

`FractionalPart[x]` is `x - IntegerPart[x]`, so it carries the **sign of `x`**:
`FractionalPart[-2.7] = -0.7`, not `0.3`. It preserves the input's precision and
keeps exact inputs exact (`FractionalPart[7/2] = 1/2`). `FractionalPart` is
`Listable`, and reconstructs the number with `IntegerPart`.
