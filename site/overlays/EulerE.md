### Worked examples

```mathematica
In[1]:= EulerE[4]
Out[1]= 5
```

```mathematica
In[1]:= Table[EulerE[2 n], {n, 0, 6}]
Out[1]= {1, -1, 5, -61, 1385, -50521, 2702765}
```

```mathematica
In[1]:= EulerE[6, x]
Out[1]= -3 x + 5 x^3 - 3 x^5 + x^6
```

```mathematica
In[1]:= Sum[EulerE[2 k] / (2 k)! Pi^(2 k), {k, 0, 5}]
Out[1]= 1 - 1/2 Pi^2 + 5/24 Pi^4 - 61/720 Pi^6 + 277/8064 Pi^8 - 50521/3628800 Pi^10
```

```mathematica
In[1]:= N[EulerE[5, 1/3], 40]
Out[1]= -0.24897119341563786008230452674897119341565
```

### Notes

`EulerE[n]` is the integer Euler number `E_n` (odd `n` vanish, `E_0 = 1`),
and `EulerE[n, x]` is the degree-`n` Euler polynomial with exact rational
coefficients. The truncated `secant`-style series above is the partial sum
of `sec(Pi/2)`'s generating expansion; the polynomial form stays symbolic in
`x` or, given an inexact argument, evaluates to arbitrary precision.
