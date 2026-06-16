### Worked examples

```mathematica
In[1]:= Cosh[0]
Out[1]= 1

In[2]:= Cosh[-x]
Out[2]= Cosh[x]
```

```mathematica
In[1]:= Cosh[Pi I]
Out[1]= -1
```

```mathematica
In[1]:= Cosh[ArcSinh[x]]
Out[1]= Sqrt[1 + x^2]
```

```mathematica
In[1]:= TrigExpand[Cosh[a + b]]
Out[1]= Cosh[a] Cosh[b] + Sinh[a] Sinh[b]
```

```mathematica
In[1]:= Series[Cosh[x], {x, 0, 6}]
Out[1]= 1 + 1/2 x^2 + 1/24 x^4 + 1/720 x^6 + O[x]^7

In[2]:= N[Cosh[1], 40]
Out[2]= 1.5430806348152437784779056207570616826015
```

### Notes

`Cosh[z]` is the hyperbolic cosine, `(Exp[z] + Exp[-z])/2`. It is even, so the sign of the argument is dropped. `Cosh` is Listable.
