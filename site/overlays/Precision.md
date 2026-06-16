### Worked examples

```mathematica
In[1]:= Precision[N[Pi, 30]]
Out[1]= 30.103

In[2]:= Precision[1]
Out[2]= Infinity

In[3]:= Precision[1.5]
Out[3]= MachinePrecision
```

Arithmetic is precision-contagious: a sum is no more precise than its least precise operand, so adding a 30-digit number to a 50-digit number yields about 30 digits:

```mathematica
In[1]:= Precision[N[Pi, 50] + N[E, 30]]
Out[1]= 30.103
```

Squaring a 100-digit square root *gains* a fraction of a digit, reflecting the conditioning of the operation:

```mathematica
In[1]:= Precision[N[Sqrt[2], 100]^2]
Out[1]= 100.243
```

### Notes

`Precision` gives the number of significant decimal digits in a number. Exact quantities (integers, rationals, symbols) have `Infinity` precision, machine-precision reals report `MachinePrecision`, and arbitrary-precision reals report their actual digit count.
