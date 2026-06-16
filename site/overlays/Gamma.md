### Worked examples

```mathematica
In[1]:= Gamma[5]
Out[1]= 24
```

```mathematica
In[1]:= Gamma[7/2]
Out[1]= 15/8 Sqrt[Pi]
```

```mathematica
In[1]:= Gamma[-1/2]
Out[1]= -2 Sqrt[Pi]
```

```mathematica
In[1]:= N[Gamma[1/3], 40]
Out[1]= 2.6789385347077476336556929409746776441289
```

```mathematica
In[1]:= N[Gamma[3 + 4 I], 20]
Out[1]= 0.00522553847136921419473 - 0.172547079294300187719*I
```

### Notes

`Gamma` is the Euler gamma function, the analytic continuation of the
factorial: `Gamma[n] = (n-1)!`, so `Gamma[5] = 24`. Half-integer arguments
collapse to exact rational multiples of `Sqrt[Pi]` — `Gamma[7/2] = 15/8 Sqrt[Pi]`
— and this continues through the poles at non-positive integers into negative
half-integers, where `Gamma[-1/2] = -2 Sqrt[Pi]`. For arguments with no
closed form it evaluates to arbitrary precision via MPFR (`Gamma[1/3]` to 40
digits) and across the complex plane (`Gamma[3 + 4 I]`). The two- and
three-argument forms give the upper incomplete and generalized incomplete gamma
functions.
