### Worked examples

```mathematica
In[1]:= EllipticPi[0, 1/2]  (* with n = 0 the third factor vanishes, leaving the first kind *)
Out[1]= (8 Pi^(3/2))/Gamma[-1/4]^2
```

```mathematica
In[1]:= EllipticPi[0, phi, m]  (* the same reduction, incomplete *)
Out[1]= EllipticF[phi, m]
```

```mathematica
In[1]:= EllipticPi[n, 0]  (* at m = 0 there is a closed form in n alone *)
Out[1]= (1/2 Pi)/Sqrt[1 - n]
```

```mathematica
In[1]:= EllipticPi[n, phi, 0]  (* and its incomplete counterpart *)
Out[1]= ArcTanh[Sqrt[-1 + n] Tan[phi]]/Sqrt[-1 + n]
```

```mathematica
In[1]:= EllipticPi[1, m]  (* a pole of the COMPLETE form only: the integrand carries 1/Cos[t]^2 at the upper limit *)
Out[1]= ComplexInfinity
```

```mathematica
In[1]:= N[EllipticPi[1/2, 1/4], 30]  (* complete, via R_F + (n/3) R_J *)
Out[1]= 2.41367150420119464066692352054
```

```mathematica
In[1]:= N[EllipticPi[1/2, Pi/3, 1/4], 25]  (* incomplete, amplitude Pi/3 *)
Out[1]= 1.3101681612463965511336307
```

```mathematica
In[1]:= N[EllipticPi[3/2, 1/2], 20]  (* for n > 1 the path crosses the pole at Sin[t]^2 == 1/n: the value is complex, not a real principal value *)
Out[1]= -0.456720313452909897007 - 2.7206990463513267759*I
```

```mathematica
In[1]:= D[EllipticPi[n, phi, m], phi]  (* the amplitude derivative is the integrand *)
Out[1]= 1/(Sqrt[1 - m Sin[phi]^2] (1 - n Sin[phi]^2))
```

```mathematica
In[1]:= D[EllipticPi[n, m], m]  (* the complete form's parameter derivative is closed form; the incomplete one's stays inert *)
Out[1]= (1/2 (EllipticE[m]/(-1 + m) + EllipticPi[n, m]))/(-m + n)
```

```mathematica
In[1]:= EllipticPi[0.5, NDArray[{0.1, 0.2, 0.3}], 0.25]  (* a visible NDArray in the amplitude slot *)
Out[1]= {0.100209, 0.201675, 0.305686}
```

```mathematica
In[1]:= CompileDiagnostics[{{x, _Real}}, EllipticPi[0.5, x, 0.25]]  (* the three-argument shape lowers too, through the n-ary kernel opcode *)
Out[1]= {"Compiled" -> True, "ResultType" -> "Real", "Instructions" -> 7, "CommonSubexpressions" -> 0, "InstructionsUnoptimized" -> 7}
```

### Notes

`EllipticPi` is arity-overloaded: **two** arguments are the complete integral `Π(n|m)`,
**three** the incomplete `Π(n; φ|m)`. The parameter argument is `m = k²`, not the modulus.
A wrong count emits `EllipticPi::argt`.

`n > 1` is the case worth understanding. The integration path crosses the pole at
`sin²t = 1/n`, and the value there is **genuinely complex** — `Π(3/2 | 1/2)` is
`−0.456720313453 − 2.72069904635 i`, with mpmath and Arb agreeing — not the real Cauchy
principal value it was long documented to be. So the `double` kernel declines there, as
`K`, `E` and `F` decline outside their own real domains, and Arb answers.

The pole at `n = 1` belongs to the complete form alone, because that is where the
`1/cos²t` in the integrand meets the upper limit `π/2`; it diverges as `1/√ε` (21.5, 221,
2221, 22214 at `ε = 10⁻², 10⁻⁴, 10⁻⁶, 10⁻⁸`). The incomplete form is finite there:
`Π[1, 1, 1/2]` is `1.73199154202`.

Both arities have machine kernels on Carlson's `R_J` (and `R_C`, which `R_J` needs): the
complete form at 41 ns/element and 2 ulp, the incomplete at 1 ulp with the quasi-period
`Π(n; φ+kπ|m) = Π(n; φ|m) + 2k Π(n|m)`. Because the element-wise NDArray layer tops out at
arity two, only the complete form rides that buffer — but `Compile[]` lowers all three
shapes, and a visible `NDArray` in any argument position of the three-argument form is
delisted and re-evaluated rather than left standing.
