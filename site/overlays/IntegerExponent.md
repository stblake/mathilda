### Worked examples

```mathematica
In[1]:= IntegerExponent[1000]
Out[1]= 3
```

```mathematica
In[1]:= IntegerExponent[1000, 2]
Out[1]= 3
```

```mathematica
In[1]:= IntegerExponent[100!, 5]
Out[1]= 24
```

```mathematica
In[1]:= IntegerExponent[20!, 2]
Out[1]= 18
```

### Notes

`IntegerExponent[n, b]` gives the highest power of `b` dividing `n`; with no
base it counts trailing decimal zeros. The factorial examples reproduce
Legendre's formula: `100!` ends in `IntegerExponent[100!, 5] = 24` zeros (the
number of factors of five), and `20!` contains `2^18`. The sign of `n` is
ignored and `IntegerExponent[0, b]` is `Infinity`.
