### Worked examples

```mathematica
In[1]:= NLimit[Sin[x]/x, x -> 0]
Out[1]= 1.0
```

```mathematica
In[1]:= NLimit[(1 + 1/n)^n, n -> Infinity]
Out[1]= 2.71828
```

```mathematica
In[1]:= NLimit[n (2^(1/n) - 1), n -> Infinity]
Out[1]= 0.693147
```

```mathematica
In[1]:= NLimit[Zeta[x] - 1/(x - 1), x -> 1]
Out[1]= 0.577216
```

```mathematica
In[1]:= NLimit[Sin[x]/x, x -> 0, Method -> "Levin"]
Out[1]= 1.0
```

### Notes

`NLimit[expr, z -> z0]` builds a geometric sequence of sample points approaching
`z0` and recovers the limit by sequence acceleration. The first three cases give
`1`, the constant `E = 2.71828...`, and `Log[2] = 0.693147...`. The fourth is the
classic Laurent-expansion limit of the Riemann zeta function at its pole: the
constant term is the Euler–Mascheroni constant `EulerGamma = 0.577216...`. `z0`
may be finite, complex, or an infinite point such as `Infinity` or `I Infinity`.
`Method -> Automatic` (default) keeps the most self-consistent of Richardson
extrapolation (`EulerSum`), Wynn's epsilon (`SequenceLimit`) and Levin's
u-transform; `Method -> "Levin"` (`"LevinU"`/`"LevinT"`/`"LevinV"`) forces
Levin's transformation. `Chop` small spurious residuals.
