### Worked examples

```mathematica
In[1]:= Complex[3, 4]
Out[1]= 3 + 4*I

In[2]:= Complex[5, 0]
Out[2]= 5
```

```mathematica
In[1]:= (1 + I)^8
Out[1]= 16
```

```mathematica
In[1]:= (2 + 3 I)/(1 - I)
Out[1]= -1/2 + 5/2*I
```

```mathematica
In[1]:= N[(3 + 4 I)^(1/3), 40]
Out[1]= 1.6289371459221758752146093717175049715341 + 0.52017450230454583954569417015944746788805*I
```

### Notes

`Complex[re, im]` is the canonical head for `re + im I`; a zero imaginary part collapses back to the underlying real number.
