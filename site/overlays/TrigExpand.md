### Worked examples

```mathematica
In[1]:= TrigExpand[Sin[a + b]]
Out[1]= Cos[a] Sin[b] + Sin[a] Cos[b]
```

Integer multiples in the argument are expanded via the multiple-angle formulas:

```mathematica
In[1]:= TrigExpand[Sin[3 x]]
Out[1]= -Sin[x]^3 + 3 Cos[x]^2 Sin[x]
```

The double-angle identity for the cosine drops out automatically:

```mathematica
In[1]:= TrigExpand[Cos[2 x]]
Out[1]= Cos[x]^2 - Sin[x]^2
```

`TrigExpand` works on hyperbolic functions as well, giving the addition formula
for `sinh`:

```mathematica
In[1]:= TrigExpand[Sinh[x + y]]
Out[1]= Cosh[x] Sinh[y] + Sinh[x] Cosh[y]
```
