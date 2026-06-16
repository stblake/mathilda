### Worked examples

```mathematica
In[1]:= Beta[3, 5]
Out[1]= 1/105
```

The central value is exactly `Pi`, and rational orders fold into Gamma quotients:

```mathematica
In[1]:= Beta[1/2, 1/2]
Out[1]= Pi

In[2]:= Beta[1/3, 1/3]
Out[2]= Gamma[1/3]^2/Gamma[2/3]
```

Positive-integer orders give the reciprocal binomial relation `1/B(7, 3) = 9 C(8, 2)`:

```mathematica
In[1]:= Beta[7, 3]
Out[1]= 1/252

In[2]:= 1/Beta[7, 3] - 9 Binomial[8, 2]
Out[2]= 0
```

The three-argument incomplete beta, and arbitrary-precision numerics:

```mathematica
In[1]:= Beta[5, 2, 3]
Out[1]= 1025/12

In[2]:= N[Beta[2.5, 3.5], 30]
Out[2]= 0.036815538909255388078101134397
```

### Notes

`Beta[a, b] = Gamma[a] Gamma[b]/Gamma[a+b]` is the Euler beta function. `Beta[z, a, b]` is the incomplete beta integral, and `Beta[z0, z1, a, b]` the generalized incomplete form. Exact for rational arguments via Pochhammer; non-positive integer poles give `ComplexInfinity`; the incomplete form reduces through `Hypergeometric2F1`. Real and complex inputs evaluate at machine or MPFR precision. Listable.
