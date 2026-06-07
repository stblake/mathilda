### Worked examples

```mathematica
In[1]:= Cosh[0]
Out[1]= 1

In[2]:= N[Cosh[1]]
Out[2]= 1.54308

In[3]:= Cosh[-x]
Out[3]= Cosh[x]
```

### Notes

`Cosh[z]` is the hyperbolic cosine, `(Exp[z] + Exp[-z])/2`. It is even, so the sign of the argument is dropped. `Cosh` is Listable.
