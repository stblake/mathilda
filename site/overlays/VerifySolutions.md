### Worked examples

```mathematica
In[1]:= VerifySolutions
Out[1]= VerifySolutions

In[2]:= Solve[x^2 - 3 x + 2 == 0, x, VerifySolutions -> True]
Out[2]= {{x -> 1}, {x -> 2}}
```

Verification survives radical-introducing transformations: solving an equation with a square root can generate a candidate that fails back-substitution, and with `VerifySolutions -> True` the spurious branch is discarded — here leaving the empty solution set for an equation with no real root:

```mathematica
In[1]:= Solve[Sqrt[x] == -2, x, VerifySolutions -> True]
Out[1]= {}
```

Confirmation extends to complex roots; all four fourth-roots of unity pass back-substitution and are kept:

```mathematica
In[1]:= Solve[x^4 == 1, x, VerifySolutions -> True]
Out[1]= {{x -> -1}, {x -> 1}, {x -> -I}, {x -> I}}
```

### Notes

`VerifySolutions` is a `Solve` option (default `Automatic`), not a function;
the bare symbol evaluates to itself. It controls whether each returned solution
is checked by back-substitution before being reported, discarding spurious
roots introduced during the solving process.
