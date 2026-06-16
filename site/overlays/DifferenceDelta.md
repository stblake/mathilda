### Worked examples

```mathematica
In[1]:= DifferenceDelta[n^2, n]
Out[1]= 1 + 2 n
```

```mathematica
In[1]:= DifferenceDelta[f[n], n]
Out[1]= -f[n] + f[1 + n]
```

```mathematica
In[1]:= DifferenceDelta[n^3, n]
Out[1]= 1 + 3 n + 3 n^2
```

```mathematica
In[1]:= DifferenceDelta[Binomial[n, k], n]
Out[1]= -Binomial[n, k] + Binomial[1 + n, k]
```

### Notes

`DifferenceDelta[f, n]` is the forward difference operator `Δ f = (f /. n -> n+1) - f`,
the discrete analogue of the derivative `D`. On `n^2` it returns `2 n + 1`, the
discrete counterpart of `2 n`; applied to `n^3` it gives `3 n^2 + 3 n + 1`. For an
unknown function head it expands literally to `f[n+1] - f[n]`. The fourth example
is Pascal's rule in disguise: `Binomial[n+1, k] - Binomial[n, k] = Binomial[n, k-1]`.
`DifferenceDelta` is the left inverse of indefinite `Sum`, mirroring the way `D`
inverts the indefinite integral.
