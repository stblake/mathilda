### Worked examples

```mathematica
In[1]:= Chop[NSeries[Cos[x], {x, 0, 6}]]
Out[1]= 1.0 - 0.5 x^2 + 0.0416667 x^4 - 0.00138889 x^6 + O[x]^7
```

```mathematica
In[1]:= Chop[NSeries[Sin[x + 1/x], {x, 0, 3}]]
Out[1]= -0.128943/x^3 + 0.576725/x + 0.576725 x - 0.128943 x^3 + O[x]^4
```

### Notes

`NSeries[f, {x, x0, n}]` recovers Taylor or Laurent coefficients by sampling `f`
on a circle in the complex plane around `x0` and taking a discrete Fourier
transform (Cauchy's integral formula). The `Cos` example reproduces the familiar
machine-precision coefficients `1, -1/2, 1/24, -1/720`. The second case is the
function's headline capability: `Sin[x + 1/x]` has an essential singularity at
the origin, so the ordinary symbolic `Series` cannot expand it, yet `NSeries`
returns its full Laurent expansion. The coefficients are Bessel values
(`0.576725 = BesselJ[1, 2]`, the `±1` terms, and `-0.128943` for the `±3`
terms). `Chop` clears the small spurious residuals from the numerical transform.
