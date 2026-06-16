### Worked examples

```mathematica
In[1]:= N[GoldenAngle]
Out[1]= 2.39996
```

```mathematica
In[1]:= N[GoldenAngle, 40]
Out[1]= 2.399963229728653322231555506633613853125
```

```mathematica
In[1]:= N[GoldenAngle/Degree, 30]
Out[1]= 137.5077640500378546463487396284
```

```mathematica
In[1]:= N[GoldenAngle - (3 - Sqrt[5]) Pi, 40]
Out[1]= 0.0
```

### Notes

`GoldenAngle` is `(3 - Sqrt[5]) Pi = 2 Pi / GoldenRatio^2`, the angle that
divides a full turn in the golden ratio — about `137.5` degrees, the divergence
angle that governs optimal phyllotactic spiral packing in plants. Dividing by
`Degree` shows the familiar `137.5...`, and subtracting the closed form
`(3 - Sqrt[5]) Pi` numerically returns `0.0`. It is a protected `Constant`
(so `D[GoldenAngle, x]` is `0`) that `N` evaluates to any precision.
