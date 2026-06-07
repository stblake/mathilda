### Worked examples

```mathematica
In[1]:= Arg[1]
Out[1]= 0

In[2]:= Arg[-1]
Out[2]= Pi

In[3]:= Arg[I]
Out[3]= 1/2 Pi

In[4]:= Arg[1 + I]
Out[4]= 1/4 Pi
```

### Notes

`Arg[z]` gives the phase angle in the range `(-Pi, Pi]`: 0 for positive reals, `Pi` for negative reals.
