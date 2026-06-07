### Worked examples

```mathematica
In[1]:= Complex[3, 4]
Out[1]= 3 + 4*I

In[2]:= (3 + 4 I) + (1 - 2 I)
Out[2]= 4 + 2*I

In[3]:= Complex[5, 0]
Out[3]= 5
```

### Notes

`Complex[re, im]` is the canonical head for `re + im I`; a zero imaginary part collapses back to the underlying real number.
