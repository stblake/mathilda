### Worked examples

```mathematica
In[1]:= Arg[1]
Out[1]= 0

In[2]:= Arg[-1]
Out[2]= Pi

In[3]:= Arg[1 + I]
Out[3]= 1/4 Pi
```

Powers wrap the phase into the principal branch:

```mathematica
In[1]:= Arg[(1 + I)^10]
Out[1]= 1/2 Pi

In[2]:= Arg[-2 + 2 I]
Out[2]= 3/4 Pi
```

Stepping the argument of `(1 + I)^k` traces the unit circle and folds back at `Pi`:

```mathematica
In[1]:= Table[Arg[(1 + I)^k], {k, 0, 8}]
Out[1]= {0, 1/4 Pi, 1/2 Pi, 3/4 Pi, Pi, -3/4 Pi, -1/2 Pi, -1/4 Pi, 0}
```

Generic complex numbers evaluate to high precision:

```mathematica
In[1]:= N[Arg[2 + 3 I], 40]
Out[1]= 0.98279372324732906798571061101466601449686
```

### Notes

`Arg[z]` gives the phase angle in the range `(-Pi, Pi]`: 0 for positive reals, `Pi` for negative reals.
