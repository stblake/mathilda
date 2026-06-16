### Worked examples

```mathematica
In[1]:= N[Glaisher]
Out[1]= 1.28243
```

```mathematica
In[1]:= N[Glaisher, 40]
Out[1]= 1.2824271291006226368753425688697917277676
```

```mathematica
In[1]:= NumericQ[Glaisher]
Out[1]= True
```

```mathematica
In[1]:= D[Glaisher, x]
Out[1]= 0
```

### Notes

`Glaisher` is the Glaisher-Kinkelin constant `A`, defined by
`Log[A] == 1/12 - Zeta'[-1]` and appearing in the asymptotics of the
hyperfactorial and in many `Zeta`-derivative identities. It is held symbolic
(attributes `Constant` and `Protected`, so `D[Glaisher, x]` is `0` and
`NumericQ` is `True`) until `N` forces a value; `N[Glaisher, 40]` returns it
to 40 digits via its MPFR series.
