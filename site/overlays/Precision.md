### Worked examples

```mathematica
In[1]:= Precision[N[Pi, 30]]
Out[1]= 30.103

In[2]:= Precision[1]
Out[2]= Infinity

In[3]:= Precision[1.5]
Out[3]= MachinePrecision
```

### Notes

`Precision` gives the number of significant decimal digits in a number. Exact quantities (integers, rationals, symbols) have `Infinity` precision, machine-precision reals report `MachinePrecision`, and arbitrary-precision reals report their actual digit count.
