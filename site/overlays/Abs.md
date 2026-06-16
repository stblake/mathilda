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

```mathematica
In[1]:= Abs[(1 + I)^10]
Out[1]= 32

In[2]:= Abs[Sqrt[2] + Sqrt[3] I]
Out[2]= Sqrt[5]

In[3]:= N[Abs[Gamma[1/3 + 2 I]], 30]
Out[3]= 0.09665959425732664141022797859867
```

### Notes

For complex `z`, `Abs[z]` returns the modulus `Sqrt[Re[z]^2 + Im[z]^2]`. Abs is Listable, so it threads element-wise over lists. Exact arguments give exact results: `Abs[(1 + I)^10]` is `32` (since `|1 + I| = Sqrt[2]` and `(Sqrt[2])^10 = 32`), and `Abs[Sqrt[2] + Sqrt[3] I]` collapses to `Sqrt[5]`. The modulus of a complex value of a special function such as `Gamma` is evaluated to the requested precision under MPFR.
