---
status: Partial
references:
  - "Bronstein, \"Symbolic Integration I: Transcendental Functions\", 2nd ed. (Springer, 2005)."
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (Kluwer, 1992), ch. 11–12."
---
### Worked examples

```mathematica
In[1]:= Integrate[1/(x^3 + 1), x]  (* rational: Hermite reduction + Lazard–Rioboo–Trager, real Log/ArcTan form *)
```

```mathematica
In[1]:= Integrate[Log[x]^2/(1 + x), x]  (* the PolyLogLadder stage: R(x) Log[linear]^k to PolyLog of rising order *)
```

```mathematica
In[1]:= Integrate[x/Sin[x], x]  (* polynomial times trig rational: TrigToExp, then the exponential-tower polylog ladder *)
```

```mathematica
In[1]:= Integrate[x^n Exp[-x], x]  (* symbolic exponent: the GammaPower recogniser gives the incomplete Gamma at once *)
```

```mathematica
In[1]:= Integrate[Log[x] Sin[x], x]  (* one integration by parts closes to CosIntegral *)
```

```mathematica
In[1]:= Integrate[(1 + x^2)/Sqrt[1 - x^4], x]  (* non-elementary: the last stage, ParallelMixedSpecial, answers with an elliptic integral *)
```

```mathematica
In[1]:= Integrate[Cos[x]/(1 + x^2), {x, -Infinity, Infinity}]  (* residue theorem, Fourier (Jordan) family: upper half-plane pole at I *)
```

```mathematica
In[1]:= Integrate[1/(2 + Cos[x]), {x, 0, 2 Pi}]  (* residue theorem on the unit circle for a periodic integrand *)
```

```mathematica
In[1]:= Integrate[x^(s - 1)/(1 + x), {x, 0, Infinity}, Assumptions -> 0 < s < 1]  (* symbolic exponent: convergence for 0 < s < 1 is proved from the Assumptions *)

In[2]:= Integrate[x^(s - 1)/(1 + x), {x, 0, Infinity}, Assumptions -> 0 < s < 1, Method -> "Mellin"]  (* the Mellin table gives the same value as Beta[s, 1 - s] = Pi Csc[Pi s] *)
```

```mathematica
In[1]:= Integrate[Exp[-a x] BesselJ[0, b x], {x, 0, Infinity}, Assumptions -> a > 0 && b > 0]  (* Laplace transform of a Bessel function *)
```

```mathematica
In[1]:= Integrate[Exp[-z Cosh[t]], {t, 0, Infinity}, Assumptions -> z > 0]  (* recognised as the integral representation of BesselK *)
```

```mathematica
In[1]:= Integrate[Log[x]/(1 + x^2), {x, 0, 1}, Method -> "SeriesExpansion"]  (* expand the geometric series, integrate term by term, resum to Catalan's constant *)
```

```mathematica
In[1]:= Integrate[(x^a - x^b)/Log[x], {x, 0, 1}, Assumptions -> a > 0 && b > 0]  (* Frullani-type logarithm, by differentiating under the integral sign in the exponent *)
```

```mathematica
In[1]:= Integrate[1/(1 + Tan[x]^Sqrt[2]), {x, 0, Pi/2}]  (* King reflection: f(x) + f(Pi/2 - x) = 1, with no elementary antiderivative needed *)

In[2]:= Integrate[x Sin[x]/(1 + Cos[x]^2), {x, 0, Pi}]  (* reflection x -> Pi - x removes the factor x *)
```

```mathematica
In[1]:= Integrate[Max[x, x^2], {x, 0, 2}]  (* piecewise: split at the breakpoint x = 1 and integrate each branch *)
```

```mathematica
In[1]:= Integrate[1/(x^3 - 1), {x, 0, 2}, PrincipalValue -> True]  (* Cauchy principal value across the pole at x = 1 *)
```

```mathematica
In[1]:= Integrate[Log[1 + a^2 x^2]/(1 + x^2), {x, 0, Infinity}, Assumptions -> a > 0]  (* differentiation under the integral sign in a, then integration back *)
```

```mathematica
In[1]:= Integrate[BesselJ[0, x] E^(-x^2), {x, 0, Infinity}, Method -> "MeijerConvolution"]  (* two kernels, each a Meijer G of x^2, combined by the convolution theorem *)
```

```mathematica
In[1]:= Integrate[Exp[-x] Sin[x]/x, {x, 0, Infinity}]  (* a known gap: the value is Pi/4 but no stage closes it, so the integral is returned unevaluated *)
```

### Notes

- **Forms.** `Integrate[f, x]` (indefinite, no constant of integration), `Integrate[f, {x, a, b}]` (definite, including improper bounds), iterated `Integrate[f, {x, a, b}, {y, c, d}]` (innermost spec last), and complex line integrals `Integrate[f, {z, z0, z1, …}]` along a polyline. A `List` integrand threads; `Integrate` itself is not `Listable`, so a range spec is never threaded over. Attributes: `Protected`.
- **Options.** `Method -> Automatic | "name"` (a named method runs that one stage with no fallback), `Assumptions -> cond` (used by the definite stages to prove convergence strips and parameter signs; the result is then `Refine`d under it), `PrincipalValue -> True` (Cauchy principal value at interior poles), and `TimeConstraint -> 3` (seconds for the automatic indefinite cascade; `Infinity` removes the bound; a pinned `Method` is unbounded unless you pass one).
- **Method names.** Indefinite: `"BronsteinRational"`, `"Undefined"`, `"LinearRadicals"`, `"QuadraticRadicals"`, `"LinearRatioRadicals"`, `"ChebychevAlgebraic"`, `"Weierstrass"`, `"DerivativeDivides"` (or `{"DerivativeDivides", "Substitution" -> u}`), `"RischTranscendental"`, `"CRCTable"`, `"PolyLogLadder"`, `"ParallelMixedTower"`, `"GeneralisedGoursat"`, `"ParallelMixedSpecial"`. Definite: `"NewtonLeibniz"`, `"Residue"`, `"Symmetry"`, `"Piecewise"`, `"Beta"`, `"TrigPower"`, `"SeriesExpansion"`, `"FeynmanParameter"`, `"King"`, `"Mellin"` (= `"RamanujanMasterTheorem"`), `"IntegralRepresentation"`, `"MeijerConvolution"`, `"SinPowerMonomial"`, `"OscillatoryPower"`, `"RationalLog"`, `"DiffUnderInt"`, `"LineIntegral"`.
- **Answer quality.** The cascade prefers an elementary answer; special functions (`Erf`, `ExpIntegralEi`, `PolyLog`, `FresnelS`, the elliptic family, …) come only from stages that need them. Antiderivatives are not post-simplified, so the form can differ from Mathematica's by a constant or by an equivalent rearrangement. `Integrate::nonelem` reports a proof of non-elementarity only when no later stage answers.
- **Soundness.** Every definite stage either is correct by construction behind proved gates (residue, Beta, Mellin strips) or is cross-checked numerically; a value still containing the integration variable is discarded, a pole on the path gives `Integrate::idiv` and stays unevaluated, and an undecidable case is left unevaluated rather than guessed.
- **Inexact input.** Machine or arbitrary-precision coefficients are rationalised, integrated exactly, and numericalised back at the input precision. For a purely numeric definite integral, `NIntegrate` is usually the better tool.
- See also: `NIntegrate`, `D`, `Sum`, `Limit`, `Residue`, `DSolve`.
