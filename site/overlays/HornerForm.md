### Worked examples

```mathematica
In[1]:= HornerForm[1 + x + x^2 + x^3]
Out[1]= 1 + x (1 + x (1 + x))

In[2]:= HornerForm[a x^3 + b x^2 + c x + d, x]
Out[2]= d + x (c + x (b + a x))
```

### Notes

The nested (Horner) form evaluates a degree-n polynomial in n multiplications and n additions instead of the naive 2n. Pass an explicit variable as the second argument when the coefficients are themselves symbolic.
