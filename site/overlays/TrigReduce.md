### Worked examples

```mathematica
In[1]:= TrigReduce[Sin[x] Cos[x]]
Out[1]= 1/2 Sin[2 x]
```

A square is linearised through the power-reduction formula:

```mathematica
In[1]:= TrigReduce[Sin[x]^2]
Out[1]= 1/2 (1 - Cos[2 x])
```

Higher powers spread across several harmonics:

```mathematica
In[1]:= TrigReduce[Cos[x]^3]
Out[1]= 1/4 (3 Cos[x] + Cos[3 x])
```

A product of squares reduces to a single fourth harmonic — exactly the
integrand identity behind `Integrate[Sin[x]^2 Cos[x]^2, x]`:

```mathematica
In[1]:= TrigReduce[Sin[x]^2 Cos[x]^2]
Out[1]= 1/8 (1 - Cos[4 x])
```

The product-to-sum identity for two sines appears directly:

```mathematica
In[1]:= TrigReduce[2 Sin[x] Sin[y]]
Out[1]= -Cos[x + y] + Cos[x - y]
```
