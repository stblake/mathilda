### Worked examples

```mathematica
In[1]:= Exp[0]
Out[1]= 1

In[2]:= Exp[I Pi]
Out[2]= -1

In[3]:= Exp[Log[x]]
Out[3]= x
```

### Notes

`Exp[z]` is `E^z`; it inverts `Log` and evaluates Euler's identity exactly. Exp is Listable, and numeric inputs route to libm / MPFR.
