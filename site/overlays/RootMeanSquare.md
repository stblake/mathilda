### Worked examples

```mathematica
In[1]:= RootMeanSquare[{3, 4}]
Out[1]= 5/Sqrt[2]
```

```mathematica
In[1]:= RootMeanSquare[{a, b}]
Out[1]= Sqrt[1/2 (a^2 + b^2)]
```

```mathematica
In[1]:= RootMeanSquare[Range[10]]
Out[1]= Sqrt[77/2]

In[2]:= N[RootMeanSquare[{1, 2, 3, 4, 5}], 30]
Out[2]= 3.316624790355399849114932736672
```

```mathematica
In[1]:= N[RootMeanSquare[Table[Sin[n], {n, 1, 1000}]], 20]
Out[1]= 0.707242937053949660224
```

### Notes

`RootMeanSquare[list]` returns `Sqrt[Mean[list^2]]` and stays **exact** on exact
input — a list of two symbols yields the closed form `Sqrt[(a^2 + b^2)/2]`. The
final example is a numerical demonstration of an analytic fact: the RMS of
`Sin` sampled over many integer arguments approaches `1/Sqrt[2] ≈ 0.70711`,
the continuous RMS of a sinusoid.
