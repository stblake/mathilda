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

```mathematica
In[1]:= Round[2.5] + Round[3.5] + Round[4.5]
Out[1]= 10
```

```mathematica
In[1]:= Round[GoldenRatio, 1/100]
Out[1]= 81/50
```

```mathematica
In[1]:= Round[N[Pi, 40], 1/10^20]
Out[1]= 157079632679489661923/50000000000000000000
```

### Notes

`Round` breaks ties to the nearest even integer (banker's rounding), so `Round[5/2]` and `Round[2.5]` both give 2. `Round[x, a]` rounds to the nearest multiple of `a`; Round is Listable.
