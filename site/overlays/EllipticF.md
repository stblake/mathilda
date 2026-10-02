### Worked examples

```mathematica
In[1]:= EllipticF[Pi/2, 1/2]  (* a full quarter period is the complete integral *)
Out[1]= (8 Pi^(3/2))/Gamma[-1/4]^2
```

```mathematica
In[1]:= EllipticF[phi, 0]  (* with m = 0 the integrand is 1 *)
Out[1]= phi
```

```mathematica
In[1]:= EllipticF[-2 x, m]  (* odd in the amplitude, by the same superficial-negativity test the trig heads use *)
Out[1]= -EllipticF[2 x, m]
```

```mathematica
In[1]:= D[EllipticF[phi, m], phi]  (* the amplitude derivative is the integrand itself *)
Out[1]= 1/Sqrt[1 - m Sin[phi]^2]
```

```mathematica
In[1]:= D[EllipticF[phi, m], m]  (* the parameter derivative, checked against a 30-digit central difference before it was allowed in *)
Out[1]= -1/2 EllipticF[phi, m]/m - 1/2 EllipticE[phi, m]/(m (-1 + m)) + 1/4 Sin[2 phi]/((-1 + m) Sqrt[1 - m Sin[phi]^2])
```

```mathematica
In[1]:= N[EllipticF[1/5, 1/2], 25]
Out[1]= 0.2006673105648029970320719

In[1]:= NIntegrate[1/Sqrt[1 - 1/2 Sin[t]^2], {t, 0, 1/5}]
Out[1]= 0.200667
```

```mathematica
In[1]:= N[EllipticF[1/5 + 3 Pi, 1/2] - (EllipticF[1/5, 1/2] + 6 EllipticK[1/2]), 20]  (* the quasi-period F(phi + k Pi) = F(phi) + 2k K, residual at 20 digits *)
Out[1]= 1.17549435082228750797e-38
```

```mathematica
In[1]:= N[EllipticF[ArcSin[3/2], 1/2], 20]  (* |z| > 1 makes ArcSin[z] complex -- the normal case for an elliptic pencil *)
Out[1]= 1.3820851603910208068 - 1.85407467730137191843*I
```

```mathematica
In[1]:= N[EllipticF[Pi/2 - 10^-8, 99/100], 20]  (* the amplitude that used to cost eight digits *)
Out[1]= 3.69563726298987467781
```

```mathematica
In[1]:= EllipticF[{0.3, 0.6, 0.9}, 0.5]  (* Listable, and the list rides the packed buffer *)
Out[1]= {0.302255, 0.618108, 0.960966}
```

### Notes

The second argument is the **parameter** `m = k²`, not the modulus `k`.

The amplitude may be of any size and may be complex, and both cases matter in practice.
Size is handled by the quasi-period `F(φ + kπ | m) = F(φ | m) + 2k K(m)`: the machine
kernel reduces `φ` into `|r| ≤ π/2`, evaluates there, and adds back `2k K(m)` — so a
shifted amplitude declines exactly when the complete integral does. Complex amplitudes
arrive unbidden, since `EllipticF[ArcSin[z], m]` is complex as soon as `|z| > 1`, which is
the normal spelling of an elliptic pencil; those go to Arb, which is defined on the whole
plane.

The `π/2 − 10⁻⁸` example is a regression guard. `R_F`'s first argument is `Cos[φ]²`, and
computing it as `1 − Sin[φ]²` gave ~100% relative error there: the result was wrong by
2.7e-08 relative, an eight-digit loss, against 6.7e-16 at a generic amplitude. Squaring
`cos` instead brings it to 1.7e-15 worst case.
