### Worked examples

```mathematica
In[1]:= FourierDCT[{1, 2, 3, 4}]  (* the default is the type-II cosine transform *)
```

```mathematica
In[1]:= FourierDCT[{1, 2, 3, 4}, 1]  (* type I; a type is one of 1..4 or "I".."IV" *)
```

```mathematica
In[1]:= FourierDCT[FourierDCT[{1., 2., 3., 4.}, 2], 3]  (* type III inverts type II *)
```

### Notes

`FourierDCT[list]` is the real orthonormal (unitary) discrete cosine transform, type II by
default. `FourierDCT[list, m]` selects a type: the integer `1`, `2`, `3`, `4` or the string
`"I"`, `"II"`, `"III"`, `"IV"`.

The types pair as inverses: **I and IV are their own inverses**, while **II and III invert
each other**. So `FourierDCT[FourierDCT[x, 2], 3]` and `FourierDCT[FourierDCT[x, 1], 1]` both
return `x` (up to roundoff). Because the transform is orthonormal there is no separate inverse
normalisation to track.

Real input gives a real result; complex input is transformed on its real and imaginary parts.
Exact input is numericalised with `N` first; a genuinely symbolic list stays unevaluated
(unlike `Fourier`, which has a symbolic regime). The transform is separable, so a rectangular
nested array is transformed per axis. Machine precision uses FFTW's real-to-real plans when
available (a direct matrix product otherwise), and MPFR input takes a direct arbitrary-
precision matrix path. A visible `NDArray[...]` is transformed on its buffer; `FourierDCT`
also lowers inside `Compile[]` (real → real).
