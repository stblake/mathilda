### Worked examples

```mathematica
In[1]:= HarmonicNumber[4]  (* the exact fourth harmonic number 25/12 *)
```

```mathematica
In[1]:= HarmonicNumber[10, 2]  (* generalized order-2 harmonic number, an exact rational *)
```

```mathematica
In[1]:= HarmonicNumber[Infinity, 2]  (* the limit is Zeta[2] = Pi^2/6 *)
```

```mathematica
In[1]:= N[HarmonicNumber[10], 20]  (* numericalised via EulerGamma + PolyGamma[0, n+1] *)
```

```mathematica
In[1]:= HarmonicNumber[n, -2]  (* non-positive order gives the Faulhaber polynomial in n *)
```

```mathematica
In[1]:= HarmonicNumber[Range[4]]  (* threads over the list *)
```

### Notes

`HarmonicNumber[n] = Sum_{i=1}^n 1/i` and `HarmonicNumber[n, r] = Sum_{i=1}^n 1/i^r`.
An exact non-negative integer `n` expands to the explicit finite sum (an exact
rational for integer `r`); `n -> Infinity` gives `Zeta[r]`; a non-positive integer
order `r = -m` gives a Faulhaber polynomial in `n` built from `BernoulliB`.

Inexact or numericizable arguments reduce through the analytic identity
`H_n^{(r)} = Zeta[r] - Zeta[r, n+1]` (and, for `r == 1`, `EulerGamma + PolyGamma[0, n+1]`),
so arbitrary precision and complex arguments pass straight through `Zeta` and
`PolyGamma`. `HarmonicNumber` carries a real `NDArray` kernel and lowers under
`Compile[]` at both scalar and rank-1 array shapes.
