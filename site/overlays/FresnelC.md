### Worked examples

```mathematica
In[1]:= FresnelC[0]  (* odd and entire, so the value at zero is zero *)
```

```mathematica
In[1]:= FresnelC[Infinity]  (* the limit value is 1/2 *)
```

```mathematica
In[1]:= N[FresnelC[1], 20]  (* arbitrary precision via the MPFR series *)
```

```mathematica
In[1]:= D[FresnelC[x], x]  (* the integrand Cos[Pi x^2/2] *)
```

```mathematica
In[1]:= FresnelC[{0.5, 1.0, 1.5}]  (* threads element-wise over the packed real list *)
```

### Notes

Mathilda uses the Pi/2-normalized (Wolfram) convention
`FresnelC[z] = Int_0^z Cos[Pi t^2/2] dt`. The function is entire and odd, with
`FresnelC[±Infinity] = ±1/2` and `FresnelC[±I Infinity] = ±I/2`.

`FresnelC` and `FresnelS` share one numeric kernel: the pair `(C, S)` is computed
together and each builtin returns its component. The real path uses a convergent
Maclaurin series for small/moderate arguments and an asymptotic expansion (DLMF
7.12) for large ones; complex arguments always use the convergent paired `A/B`
series, correct across the whole plane. `FresnelC` carries a real `NDArray`
kernel and lowers under `Compile[]` at both scalar and rank-1 array shapes.
