### Worked examples

```mathematica
In[1]:= ExpIntegralEi[0]
Out[1]= -Infinity
```

```mathematica
In[1]:= D[ExpIntegralEi[z], z]
Out[1]= E^z/z
```

```mathematica
In[1]:= N[ExpIntegralEi[1], 40]
Out[1]= 1.8951178163559367554665209343316342690171
```

```mathematica
In[1]:= N[ExpIntegralEi[I], 30]
Out[1]= 0.3374039229009681346626462038893 + 2.516879397162079634172675005462*I
```

### Notes

`ExpIntegralEi[z]` is the exponential integral `Ei(z)`, with a branch cut on
`(-Infinity, 0)` and derivative `E^z/z`. On the imaginary axis it ties to the
cosine/sine integrals via `Ei(I) = Ci(1) + I (Pi/2 + Si(1))`. Real and complex
arguments evaluate at machine or arbitrary (MPFR) precision. Listable.
