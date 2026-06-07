### Worked examples

```mathematica
In[1]:= Sinh[0]
Out[1]= 0

In[2]:= N[Sinh[1]]
Out[2]= 1.1752

In[3]:= Sinh[-x]
Out[3]= -Sinh[x]

In[4]:= Sinh[ArcSinh[x]]
Out[4]= x
```

### Notes

`Sinh[z]` is the hyperbolic sine, `(Exp[z] - Exp[-z])/2`. It is odd, so negative arguments pull the sign out front. `Sinh` is Listable.
