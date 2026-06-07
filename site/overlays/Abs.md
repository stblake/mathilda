### Worked examples

```mathematica
In[1]:= Abs[-5]
Out[1]= 5

In[2]:= Abs[3 + 4 I]
Out[2]= 5

In[3]:= Abs[-3/4]
Out[3]= 3/4

In[4]:= Abs[{-1, 2, -3}]
Out[4]= {1, 2, 3}
```

### Notes

For complex `z`, `Abs[z]` returns the modulus `Sqrt[Re[z]^2 + Im[z]^2]`. Abs is Listable, so it threads element-wise over lists.
