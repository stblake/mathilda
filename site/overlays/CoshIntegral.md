### Worked examples

```mathematica
In[1]:= CoshIntegral[0]  (* logarithmic singularity at the origin *)
```

```mathematica
In[1]:= CoshIntegral[Infinity]  (* grows without bound on the positive real axis *)
```

```mathematica
In[1]:= N[CoshIntegral[1], 20]  (* arbitrary-precision value through the MPFR kernel *)
```

```mathematica
In[1]:= D[CoshIntegral[x], x]  (* the derivative is Cosh[x]/x *)
```

```mathematica
In[1]:= CoshIntegral[{1.0, 2.0, 3.0}]  (* threads element-wise over the packed real list *)
```

### Notes

`CoshIntegral[z] = Chi(z) = EulerGamma + Log[z] + Int_0^z (Cosh[t] - 1)/t dt`
is the hyperbolic sibling of `CosIntegral` (`Chi(z) = Ci(i z) - i Pi/2`). It has
a logarithmic singularity at `0` and a branch cut along the negative real axis;
a negative real argument returns the from-above branch value `Chi(|x|) + i Pi`.

The numeric kernel sums a convergent Maclaurin series for moderate `|z|` and an
asymptotic expansion for large `|z|` (where `Chi` grows like `e^{|z|}`;
machine-real results that overflow a C `double` are kept as extended-exponent
reals). `CoshIntegral` carries a real `NDArray` kernel and lowers under
`Compile[]` at both scalar and rank-1 array shapes.
