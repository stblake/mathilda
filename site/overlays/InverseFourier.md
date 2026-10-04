### Worked examples

```mathematica
In[1]:= InverseFourier[{1, 1, 1, 1}]  (* a constant spectrum inverts to a single spike *)
```

```mathematica
In[1]:= InverseFourier[{2, 0, 0, 0}]  (* only the zero-frequency term, so a constant signal *)
```

```mathematica
In[1]:= InverseFourier[Fourier[{1, 2, 3, 4}]]  (* InverseFourier is the exact inverse of Fourier *)
```

```mathematica
In[1]:= Fourier[InverseFourier[{1., 2., 3., 4.}]]  (* and the composition the other way is identity too *)
```

### Notes

`InverseFourier` is `Fourier` with the opposite exponent sign and the inverse normalisation:
with `FourierParameters -> {a, b}` (default `{0, 1}`) it computes
`u[r] = n^{-(1+a)/2} Sum_s v[s] Exp[-2 Pi I b (r-1)(s-1)/n]`. Under the default convention
both transforms are unitary, so `InverseFourier[Fourier[list]] === list` (and the reverse
composition too).

Everything else matches `Fourier`: the same three regimes (machine FFT, arbitrary-precision
MPFR FFT, exact symbolic), the same `FourierParameters` conventions, the same multidimensional
and `Fourier[list, {p1, …}]` position forms, the same `NDArray` buffer fast path, and the same
`Compile[]` lowering (always a complex array). A machine result real to within roundoff
collapses to a real list.
