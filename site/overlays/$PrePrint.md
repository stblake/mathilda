### Worked examples

```mathematica
In[1]:= $PrePrint = Framed
Out[1]= Framed[Framed]

In[2]:= 3 + 4
Out[2]= Framed[7]
```

### Notes

`$PrePrint`, if set, is applied to every expression just before it is printed.
The displayed form reflects the `$PrePrint` value (here wrapped in `Framed`),
but `Out[n]` is assigned the unmodified result, so later references see the
plain value. Unset by default.
