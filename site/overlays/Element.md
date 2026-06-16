### Worked examples

```mathematica
In[1]:= Element[5, Integers]
Out[1]= True
```

```mathematica
In[1]:= Element[5/2, Integers]
Out[1]= False
```

```mathematica
In[1]:= Element[7, Primes]
Out[1]= True

In[2]:= Element[1 + I, Algebraics]
Out[2]= True
```

```mathematica
In[1]:= Element[{2, 3, 5, 7}, Primes]
Out[1]= True
```

```mathematica
In[1]:= Assuming[Element[x, Integers], Element[x, Reals]]
Out[1]= True
```

### Notes

`Element[x, dom]` returns `True` if `x` is provably in the domain, `False` if
provably not, and stays unevaluated otherwise. Supported domains are `Integers`,
`Rationals`, `Reals`, `Algebraics`, `Complexes`, `Booleans`, `Primes`, and
`Composites`. Numeric and structural literals decide directly — including
membership in `Primes` and in `Algebraics` (e.g. the Gaussian integer `1 + I`).
A list or `Alternatives` of variables is shorthand for the conjunction over the
components. For symbolic queries `Element` consults `$Assumptions`, so under
`Assuming[Element[x, Integers], ...]` it climbs the domain lattice and reports
that an integer is also real.
