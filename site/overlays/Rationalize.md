### Worked examples

```mathematica
In[1]:= Rationalize[0.5]
Out[1]= 1/2

In[2]:= Rationalize[N[Pi], 10^-4]
Out[2]= 333/106

In[3]:= Rationalize[1.2 + 6.7 x]
Out[3]= 6/5 + 67/10 x
```

Tightening the tolerance recovers the continued-fraction convergents. With the
golden ratio it returns a ratio of consecutive Fibonacci numbers, and with `π`
the celebrated convergent `355/113` and beyond:

```mathematica
In[1]:= Rationalize[N[GoldenRatio], 10^-6]
Out[1]= 1597/987

In[2]:= Rationalize[N[Pi, 40], 10^-10]
Out[2]= 312689/99532
```

`Rationalize[x, 0]` forces conversion of an inexact number using a tolerance
derived from its precision:

```mathematica
In[1]:= Rationalize[N[E], 0]
Out[1]= 325368125/119696244
```

### Notes

`Rationalize[x]` finds a nearby rational with small denominator (within `c/q^2`, `c = 10^-4`); the two-argument form `Rationalize[x, dx]` returns the smallest-denominator rational within `dx` of `x`. It threads over compound expressions.
