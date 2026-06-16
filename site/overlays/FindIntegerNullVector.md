### Worked examples

```mathematica
In[1]:= FindIntegerNullVector[{N[Zeta[2], 40], N[Pi^2, 40]}]
Out[1]= {-6, 1}
```

```mathematica
In[1]:= FindIntegerNullVector[{N[GoldenRatio, 40]^2, N[GoldenRatio, 40], 1}]
Out[1]= {-1, 1, 1}
```

```mathematica
In[1]:= FindIntegerNullVector[{N[Log[2], 40], N[Log[3], 40], N[Log[6], 40]}]
Out[1]= {-1, -1, 1}
```

```mathematica
In[1]:= FindIntegerNullVector[{N[Cos[Pi/7], 40]^3, N[Cos[Pi/7], 40]^2, N[Cos[Pi/7], 40], 1}]
Out[1]= {8, -4, -4, 1}
```

### Notes

`FindIntegerNullVector` is integer-relation detection (PSLQ-style): given numerical values it recovers an exact integer combination summing to zero. Feeding `{ζ(2), π²}` recovers `−6 ζ(2) + π² = 0`, i.e. `ζ(2) = π²/6`. The powers of the golden ratio return `{-1, 1, 1}`, the minimal polynomial `−φ² + φ + 1 = 0`. The logarithms recover `log 6 = log 2 + log 3`. Most strikingly, the powers of `cos(π/7)` recover its minimal polynomial `8 x³ − 4 x² − 4 x + 1 = 0`, reconstructing exact algebraic structure from 40-digit floating-point samples. Working precision must comfortably exceed the size of the relation sought.
