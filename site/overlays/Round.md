### Worked examples

```mathematica
In[1]:= Round[7/2]
Out[1]= 4

In[2]:= Round[5/2]
Out[2]= 2

In[3]:= Round[17, 5]
Out[3]= 15

In[4]:= Round[{1.4, 2.5, 3.6}]
Out[4]= {1, 2, 4}
```

### Notes

`Round` breaks ties to the nearest even integer (banker's rounding), so `Round[5/2]` and `Round[2.5]` both give 2. `Round[x, a]` rounds to the nearest multiple of `a`; Round is Listable.
