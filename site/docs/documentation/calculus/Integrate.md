# Integrate

!!! warning "Status: Partial"
    implemented with documented limitations or caveats; some argument forms fall through to symbolic/unevaluated output.

## Description

**`Integrate[f, x] gives the indefinite integral of f with respect to x.`**

**`Integrate[f, {x, xmin, xmax}] gives the definite integral by the`**

**`Integrate[f, {x, xmin, xmax}, {y, ymin, ymax}, ...] gives the iterated`**

**`Integrate[f, x, Method -> "<name>"] dispatches directly to a single`**

<details>
<summary>Notes</summary>

fundamental theorem of calculus (Method -\> "NewtonLeibniz"). multiple integral (innermost/last spec integrated first; inner bounds may depend on outer variables).  See also Integrate\`SingularPoints. subroutine, bypassing the default cascade.  Accepted method names: "Automatic"          — try BronsteinRational, then RischNorman, then CRCTable (default) "BronsteinRational"  — Integrate\`BronsteinRational (polynomial / rational) "DerivativeDivides"  — Integrate\`DerivativeDivides (substitution u(x); direct + Eliminate/Solve) "LinearRadicals"     — Integrate\`LinearRadicals (rationalise radicals of a x + b) "QuadraticRadicals"  — Integrate\`QuadraticRadicals (Euler substitution for Sqrt\[a x^2 + b x + c\]) "LinearRatioRadicals" — Integrate\`LinearRatioRadicals (rationalise radicals of (a x + b)/(c x + d)) "ChebychevAlgebraic" — Integrate\`ChebychevAlgebraic (binomial x^p (a x^r + b)^q via Chebychev's theorem) "GoursatAlgebraic"   — Integrate\`GoursatAlgebraic (pseudo-elliptic F/R^p, p in {1/2,1/3,2/3,1/4,3/4}, via Mobius eigendescent) "Weierstrass"        — Integrate\`Weierstrass (continuous tan(x/2) / tanh(x/2) substitution) "RischNorman"        — Integrate\`RischNorman (Bronstein pmint heuristic) "RischTranscendental"       — Integrate\`RischTranscendental (recursive transcendental Risch; correct by construction) "CRCTable"           — Integrate\`CRCTable (lazy-loaded CRC integral table) "Undefined"          — Integrate\`Undefined (unknown functions u\[x\], u'\[x\]; Roach §1.7) "NewtonLeibniz"       — real definite integrals via F(b)-F(a) (implicit for the {x,a,b} form) "LineIntegral"        — complex contour integrals (implicit for the {x,z0,...,zn} form) "Residue"             — improper/periodic real definite integrals by the residue theorem (rational/Fourier on (-Inf,Inf), rational-in-Sin/Cos over a period, principal values, even half-lines); tried before NewtonLeibniz under Automatic "DiffUnderInt"         — parameter-dependent definite integrals by differentiation under the ("DifferentiationUnderIntegral") integral sign (Feynman's trick): Integrate\`DiffUnderInt; Laplace/Fourier, sinc, and even-rational half-line families; tried after Residue and NewtonLeibniz in the definite cascade "RamanujanMasterTheorem" — half-line Int\_0^Inf x^(s-1) f(x) dx by the Mellin transform / ("Mellin")              Ramanujan Master Theorem: Integrate\`RamanujanMasterTheorem; exp/Gaussian/algebraic/Cos/Sin/ArcTan/Log/BesselJ/pFq/PolyLog kernels (monomial x^k substitution; Erf, incomplete Gamma, BesselJ^2 reduced to pFq); also the exp-geometric kernel 1/(E^(c x)+g) (Bose-Einstein / Fermi-Dirac -\> Gamma\*PolyLog), a Frullani pre-pass (f(a x)-f(b x))/x -\> (f(0)-f(Inf)) Log\[b/a\], and a Log\[x\]^k weight; strip-gated, yielding a ConditionalExpression when Assumptions do not prove convergence; after NewtonLeibniz under Automatic Method -\> {"DerivativeDivides", "Substitution" -\> u} pins the kernel u(x), trialing only that substitution. Named methods are strict: failure returns unevaluated, with no fallback. The CRCTable rules are loaded from disk on first use only. An applied 1-D InterpolatingFunction integrates to its antiderivative InterpolatingFunction (mirroring D).

</details>

## Examples (70)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= Integrate[3 + 5 x + 2 x^2, x]
Out[1]= 3 x + 5/2 x^2 + 2/3 x^3

In[2]:= Integrate[2 x/(x^2 + 1), x]
Out[2]= Log[1 + x^2]

In[3]:= Integrate[1/(x - a)^2, x]
Out[3]= 1/(a - x)

In[4]:= Integrate[(2x+3)/(x^2+3x+5)^2, x]
Out[4]= -1/(5 + 3 x + x^2)
```

Phase 2 LRT closes this

```mathematica
In[5]:= Integrate[1/((x-1)(x-2)(x-3)), x]
Out[5]= -Log[-2 + x] + 1/2 Log[3 - 4 x + x^2]
```

Phase 4 LogToReal

```mathematica
In[6]:= Integrate[1/(x^2 + 1), x]
Out[6]= ArcTan[x]
```

Two quadratic factors

```mathematica
In[7]:= Integrate[1/(x^4 + x^2 + 1), x] 1/6 Sqrt[3] ArcTan[(1 + 2 x)/Sqrt[3]] + 1/4 Log[1 + x + x^2] - 1/4 Log[1 - x + x^2]
Out[7]= 1/4 Log[1 + x + x^2] - 1/4 Log[1 - x + x^2] + 1/2 (ArcTan[(1 + 2 x)/Sqrt[3]] (1/4 Log[1 + x + x^2] + 1/2 ArcTan[(-1 + 2 x)/Sqrt[3]]/Sqrt[3] + 1/2 ArcTan[(1 + 2 x)/Sqrt[3]]/Sqrt[3] - 1/4 Log[1 - x + x^2]))/Sqrt[3]
```

### Scope (38)

```mathematica
In[8]:= Integrate[x^2, {x, 0, 1}]
Out[8]= 1/3

In[9]:= Integrate[1/(1 + x^2), {x, 0, Infinity}]
Out[9]= 1/2 Pi
```

Improper, convergent

```mathematica
In[10]:= Integrate[1/Sqrt[x], {x, 0, 1}]
Out[10]= 2
```

Continuous branch form

```mathematica
In[11]:= Integrate[1/(2 + Cos[x]), {x, 0, 2 Pi}]
Out[11]= (2 Pi)/Sqrt[3]
```

Iterated

```mathematica
In[12]:= Integrate[x y, {x, 0, 1}, {y, 0, 1}]
Out[12]= 1/4
```

```mathematica
In[13]:= Integrate`SingularPoints[1/((x - 1)(x - 2)), {x, 0, 3}]
Out[13]= {1, 2}

In[14]:= Integrate[1/x, {x, 1 - I, 2 + 3 I}]
Out[14]= Log[2 + 3*I] - Log[1 - I]

In[15]:= Integrate[z^2, {z, 0, 1 + I}]
Out[15]= -2/3 + 2/3*I
```

Branch-point endpoint

```mathematica
In[16]:= Integrate[1/Sqrt[z], {z, 0, 1 + I}]
Out[16]= 2 Sqrt[1 + I]
```

```mathematica
In[17]:= Integrate`PathSingularPoints[1/z, {z, -1 - I, 1 + I}]
Out[17]= {0}
```

CCW loop about 0

```mathematica
In[18]:= Chop[N[Integrate[1/z, {z, 1, I, -1, -I, 1}]]]
Out[18]= 0.0 + 6.28319*I
```

```mathematica
In[19]:= Integrate[1/(1 + x^4), {x, -Infinity, Infinity}]
Out[19]= Pi/Sqrt[2]
```

Order-2 pole

```mathematica
In[20]:= Integrate[1/(1 + x^2)^2, {x, -Infinity, Infinity}]
Out[20]= 1/2 Pi
```

```mathematica
In[21]:= Integrate[Cos[x]/(1 + x^2), {x, -Infinity, Infinity}]
Out[21]= Pi/E

In[22]:= Integrate[1/(2 + Cos[x]), {x, 0, 2 Pi}]
Out[22]= (2 Pi)/Sqrt[3]
```

Principal value

```mathematica
In[23]:= Integrate[Sin[x]/x, {x, -Infinity, Infinity}]
Out[23]= Pi
```

Even half-line

```mathematica
In[24]:= Integrate[1/(1 + x^4), {x, 0, Infinity}]
Out[24]= (1/2 Pi)/Sqrt[2]
```

```mathematica
In[25]:= Integrate[x f'[x] + f[x], x]
Out[25]= x f[x]

In[26]:= Integrate[f'[x] g[x] + f[x] g'[x], x]
Out[26]= f[x] g[x]

In[27]:= Integrate[f'[x]/f[x], x]
Out[27]= Log[f[x]]

In[28]:= Integrate[(f'[x] g'[x] - f[x] g''[x])/(f[x]^2 + g'[x]^2), x]
Out[28]= -ArcTan[Derivative[1][g][x]/f[x]]
```

Composite argument

```mathematica
In[29]:= Integrate[2 x f'[x^2], x]
Out[29]= f[x^2]
```

```mathematica
In[30]:= Integrate[(f[x] - x f[x] + f[x] Log[x f[x]] + x f'[x])/f[x], x]
Out[30]= -1/2 x^2 + x Log[x f[x]]

In[31]:= Integrate[Log[f[x]] f'[x]/f[x], x]
Out[31]= 1/2 Log[f[x]]^2
```

Direct, correct branch

```mathematica
In[32]:= Integrate[Sin[x] Sqrt[1 - Cos[x]], x]
Out[32]= 2/3 (1 - Cos[x])^(3/2)
```

```mathematica
In[33]:= Integrate[1/(x Log[x]), x]
Out[33]= Log[Log[x]]

In[34]:= Integrate`DerivativeDivides[2 x Exp[x^2], x]
Out[34]= E^x^2
```

U = Sqrt[Tan[x]]

```mathematica
In[35]:= Integrate[Sqrt[Tan[x]], x] - Log[1 + Tan[x] + Sqrt[2] Sqrt[Tan[x]]]/(2 Sqrt[2]) + Log[1 + Tan[x] - Sqrt[2] Sqrt[Tan[x]]]/(2 Sqrt[2])
Out[35]= ArcTan[-1 + Sqrt[2] Sqrt[Tan[x]]]/Sqrt[2] + ArcTan[1 + Sqrt[2] Sqrt[Tan[x]]]/Sqrt[2] - Log[1 + Tan[x] + Sqrt[2] Sqrt[Tan[x]]]/Sqrt[2] + Log[1 + Tan[x] - Sqrt[2] Sqrt[Tan[x]]]/Sqrt[2]
```

U = Sqrt[Cot[x]]

```mathematica
In[36]:= Integrate[Sqrt[Cot[x]], x] - Log[1 + Cot[x] - Sqrt[2] Sqrt[Cot[x]]]/(2 Sqrt[2]) + Log[1 + Cot[x] + Sqrt[2] Sqrt[Cot[x]]]/(2 Sqrt[2])
Out[36]= -ArcTan[-1 + Sqrt[2] Sqrt[Cot[x]]]/Sqrt[2] - ArcTan[1 + Sqrt[2] Sqrt[Cot[x]]]/Sqrt[2] - Log[1 + Cot[x] - Sqrt[2] Sqrt[Cot[x]]]/Sqrt[2] + Log[1 + Cot[x] + Sqrt[2] Sqrt[Cot[x]]]/Sqrt[2]
```

```mathematica
In[37]:= Integrate[1/Sqrt[x + 1], x]
Out[37]= 2 Sqrt[1 + x]

In[38]:= Integrate[Sqrt[x]/(1 + Sqrt[x]), x]
Out[38]= -2 Sqrt[x] + x + 2 Log[1 + Sqrt[x]]

In[39]:= Integrate[1/Sqrt[x^2 + 1], x]
Out[39]= ArcSinh[x]

In[40]:= Integrate[1/Sqrt[1 - x^2], x]
Out[40]= ArcSin[x]
```

Symbolic leading coeff

```mathematica
In[41]:= Integrate[1/Sqrt[a x^2 + 1], x]
Out[41]= -Log[-Sqrt[a] x + Sqrt[1 + a x^2]]/Sqrt[a]
```

```mathematica
In[42]:= Integrate[1/Sqrt[(x + 1)/(x - 1)], x] 2 ArcTanh[Sqrt[(1 + x)/(-1 + x)]]
Out[42]= 2 ArcTanh[Sqrt[(1 + x)/(-1 + x)]] (2 Sqrt[(1 + x)/(-1 + x)]/(-1 + (1 + x)/(-1 + x)) - 2 ArcTanh[Sqrt[(1 + x)/(-1 + x)]])
```

Paper eq. (10)

```mathematica
In[43]:= Integrate[3/(5 - 4 Cos[x]), x]
Out[43]= 2 Pi Floor[(1/2 (-Pi + x))/Pi] + 2 ArcTan[3 Tan[1/2 x]]
```

Continuous form

```mathematica
In[44]:= Integrate[1/(2 + Cos[x]), x]
Out[44]= 2 (Pi Floor[(1/2 (-Pi + x))/Pi])/Sqrt[3] + 2 ArcTan[Tan[1/2 x]/Sqrt[3]]/Sqrt[3]
```

Hyperbolic: no Floor

```mathematica
In[45]:= Integrate[1/(2 + Cosh[x]), x]
Out[45]= (2 ArcTanh[Tanh[1/2 x]/Sqrt[3]])/Sqrt[3]
```

### Options (7)

Strict, no fallback

```mathematica
In[46]:= Integrate[Sin[x], x, Method -> "RischNorman"]
Out[46]= -Cos[x]
```

```mathematica
In[47]:= Integrate[x^3, x, Method -> "BronsteinRational"]
Out[47]= 1/4 x^4

In[48]:= Integrate[Sin[x]/Cos[x]^2, x, Method -> "DerivativeDivides"]
Out[48]= Sec[x]

In[49]:= Integrate[1/(1 + x^(1/3)), x, Method -> "LinearRadicals"]
Out[49]= -3 x^(1/3) + 3/2 x^(2/3) + 3 Log[1 + x^(1/3)]

In[50]:= Integrate[Sqrt[2 x + 3]/x, x, Method -> "LinearRadicals"]
Out[50]= 2 Sqrt[3 + 2 x] - 2 Sqrt[3] ArcCoth[Sqrt[3 + 2 x]/Sqrt[3]]

In[51]:= Integrate[1/Sqrt[x^2 - 1], x, Method -> "QuadraticRadicals"]
Out[51]= Log[2 x + 2 Sqrt[-1 + x^2]]

In[52]:= Integrate[1/Sqrt[(2 x + 1)/(x + 3)], x, Method -> "LinearRatioRadicals"] 5/2 ArcTanh[Sqrt[(1 + 2 x)/(3 + x)]/Sqrt[2]]/Sqrt[2]
Out[52]= (5/2 ArcTanh[Sqrt[(1 + 2 x)/(3 + x)]/Sqrt[2]] (-5 Sqrt[(1 + 2 x)/(3 + x)]/(-4 + 2 (1 + 2 x)/(3 + x)) + 5/2 ArcTanh[Sqrt[(1 + 2 x)/(3 + x)]/Sqrt[2]]/Sqrt[2]))/Sqrt[2]
```

### Worked examples (10)

```mathematica
In[53]:= Integrate[(x^a-1)/Log[x], {x,0,1}]
Out[53]= Log[1 + a]

In[54]:= Integrate[Exp[-a x] Sin[b x]/x, {x,0,Infinity}, Assumptions->a>0]
Out[54]= 1/2 (Pi b)/Sqrt[b^2] - ArcTan[a/b]

In[55]:= Integrate[Sin[a x]^2/x^2, {x,0,Infinity}, Assumptions->a>0]
Out[55]= 1/2 Pi a

In[56]:= Integrate[Log[1+a^2 x^2]/(1+x^2), {x,0,Infinity}, Assumptions->a>0]
Out[56]= Pi Log[1 + a]

In[57]:= Integrate[Exp[-c x](1-Cos[a x])/x^2, {x,0,Infinity}, Assumptions->{a>0,c>0}]
Out[57]= a ArcTan[a/c] - 1/2 c Log[1 + a^2/c^2]

In[58]:= Integrate[Exp[-x^2] Sin[a x]/x, {x,0,Infinity}]
Out[58]= 1/2 Pi Erf[1/2 a]

In[59]:= Integrate[Exp[-x^2], {x,0,Infinity}]
Out[59]= 1/2 Sqrt[Pi]

In[60]:= Integrate[x^(s-1) Exp[-x], {x,0,Infinity}]
Out[60]= ConditionalExpression[Gamma[s], s > 0]

In[61]:= Integrate[x^(s-1) BesselJ[ν,2√x]/x^(ν/2), {x,0,Infinity}]
Out[61]= Integrate[BesselJ[u03bd, 2 u221ax] x^(-1 + s - 1/2 u03bd), {x, 0, Infinity}]

In[62]:= Integrate[x^(s-1) (Γ[a]-Γ[a,x])/x^a, {x,0,Infinity}]
Out[62]= Integrate[x^(-1 - a + s) (u0393[a] - u0393[a, x]), {x, 0, Infinity}]
```

### Applications (8)

```mathematica
In[63]:= Integrate[1/(1 + x^2), x]
Out[63]= ArcTan[x]

In[64]:= Integrate[1/x, x]
Out[64]= Log[x]

In[65]:= Integrate[Cos[x], x]
Out[65]= Sin[x]

In[66]:= Integrate[x^3 + x, x]
Out[66]= 1/2 x^2 + 1/4 x^4

In[67]:= Integrate[1/(x^3 + 1), x]
Out[67]= 1/3 Log[1 + x] + ArcTan[(-1 + 2 x)/Sqrt[3]]/Sqrt[3] - 1/6 Log[1 - x + x^2]

In[68]:= Integrate[(x^2 + 1)/(x^4 + 1), x]
Out[68]= ArcTan[x/Sqrt[2]]/Sqrt[2] + ArcTan[(x + x^3)/Sqrt[2]]/Sqrt[2]

In[69]:= Integrate[x*Exp[x], x]
Out[69]= -E^x + x E^x

In[70]:= Integrate[1/(x*Log[x]), x]
Out[70]= Log[Log[x]]
```

## Options & behaviour

### Differentiation under the integral sign (`Integrate\`DiffUnderInt`)

For a parameter-dependent definite integral `I(p) = Integrate[f(x,p), {x,a,b}]`,
this method (Leibniz rule / "Feynman's trick") differentiates the integrand with
respect to a free parameter `p`, evaluates the resulting simpler definite
integral `J(p) = Integrate[D[f,p], {x,a,b}]`, integrates `J(p)` back over the
parameter, and fixes the constant of integration from an **exact** base value
`I(p0)` (a `p` where `f` vanishes identically, or reduces to a directly-
integrable form). Every case is the first-order ODE `I'(p) = J(p)`
(Boulnois 2023). Verification is symbolic and correct-by-construction
(`Simplify[D[I,p] - J] === 0` plus the exact base) — there is **no** numeric
crosscheck. Assumptions (`a > 0`, …) are honoured and used to clean the closed
forms.

Because the general integrator is slow/hangs on the parameter-dependent inner
integrals Feynman's trick produces, `DiffUnderInt` evaluates the standard
families itself with closed-form formulas: the **Laplace/Fourier half-line**
`∫₀^∞ xⁿ e^{-p x}{1,cos,sin} dx`, the **sinc/Frullani** `∫₀^∞ …/x dx`, the
**even-rational half-line** `∫₀^∞ P(x)/Q(x²) dx`, the **general (non-even)
rational half-line** `∫₀^∞ R(s) ds` (real `ArcTan`/`Log` boundary values — this
is what closes a *decaying* sinc such as `∫₀^∞ e^{-p x} Sin[q x]/x dx = ArcTan[q/p]`
whose Laplace image is non-even), and the **Gaussian moment** family
`∫₀^∞ xⁿ e^{-p x²}{1,cos} dx` in `Sqrt[Pi]`/`e^{-q²/4p}`. The Gaussian
parameter back-integration `∫ c e^{-k p²} dp` is supplied directly as an `Erf`
(the engine does not produce it). Forms outside these families (finite-period
trig, piecewise/`Min`-`Max` results, the Sin-Gaussian Dawson/Erfi moment) are
declined — the integral is returned unevaluated, fast, never a wrong value. The
two rational half-line families gate on the inner integrand being a **rational
function of `x`**: a differentiated exp-geometric/Mellin form (still carrying
`e^x` or `x^{s-1}`) is not, and feeding it to `Apart[·, x]` otherwise drives a
non-terminating rewrite — so such forms are declined up front and left for the
Ramanujan/Mellin method.

.

### Mellin / Ramanujan Master Theorem (`Integrate\`RamanujanMasterTheorem`)

The series/transform-based mechanism for half-line integrals
`∫₀^∞ x^{s-1} f(x) dx` of a *transcendental* `f` (the class residue and FTC do
not close). By Ramanujan's Master Theorem, if `f(x) = Σ (-1)^k φ(k) x^k / k!`
then the Mellin transform is `Γ(s) φ(-s)` on the fundamental strip. The
integrand is Expanded and, term by term, decomposed into `C · x^ρ · f(λx)`; the
kernel `f` is matched against a table of proven base Mellin transforms and the
power prefactor sets `s = ρ + 1`:

| kernel `f` | `∫₀^∞ x^{s-1} f dx` | strip |
|------------|--------------------|-------|
| `Exp[c x]`, `Re c<0` | `Γ(s) (-c)^{-s}` | `0<Re s` |
| `Exp[c x^2]`, `Re c<0` | `½ (-c)^{-s/2} Γ(s/2)` | `0<Re s` |
| `(p + q x^m)^{-a}`, `p,q>0` | `(1/m) p^{s/m-a} q^{-s/m} B(s/m, a-s/m)` | `0<Re s<m Re a` |
| `Cos[λ x]`, `λ>0` | `π / (2 Sin(πs/2) Γ(1-s)) λ^{-s}` | `0<Re s<1` |
| `Sin[λ x]`, `λ>0` | `π / (2 Cos(πs/2) Γ(1-s)) λ^{-s}` | `-1<Re s<1` |
| `Log[1 + λ x]`, `λ>0` | `π / (s Sin(π s)) λ^{-s}` | `-1<Re s<0` |
| `ArcTan[λ x]`, `λ>0` | `-π / (2 s Cos(πs/2)) λ^{-s}` | `-1<Re s<0` |
| `BesselJ[ν, λ x]`, `λ>0` | `2^{s-1} λ^{-s} Γ((ν+s)/2)/Γ((ν-s)/2+1)` | `-Re ν<Re s<3/2` |
| `pFq[{a}; {b}; -λ x]`, `λ>0` | `(∏Γ(b_j)/∏Γ(a_i)) Γ(s) (∏Γ(a_i-s)/∏Γ(b_j-s)) λ^{-s}` | `0<Re s<min Re a_i` |
| `PolyLog[ν, -λ x]`, `λ>0` | `π (-s)^{-ν} λ^{-s} / Sin(π s)` | `-1<Re s<0` |
| `1/(e^{c x}+γ)`, `c>0`, `-1≤γ≤1` | `c^{-s} Γ(s) (-1/γ) PolyLog(s, -γ)` | `0<Re s` (`1<Re s` if `γ=-1`) |

The last row is the **exponential-geometric** kernel of the statistical-mechanics
integrals: expanding `1/(e^{cx}+γ) = (-1/γ) Σ_{j≥1} (-γ)^j e^{-jcx}` and
integrating term by term lands on `PolyLog`. Its two headline specialisations are
`γ=-1` **Bose–Einstein** `∫₀^∞ x^{s-1}/(e^{cx}-1) = c^{-s} Γ(s) ζ(s)` (Planck /
Debye; the denominator zero at `x=0` tightens the strip to `Re s>1`) and `γ=+1`
**Fermi–Dirac** `∫₀^∞ x^{s-1}/(e^{cx}+1) = c^{-s} Γ(s) η(s)` (emitted as
`-Γ(s) PolyLog(s,-1)`, which stays finite at `s=1` where `(1-2^{1-s})ζ(s)` would
be `0·∞`). A **symbolic fugacity** is admitted too — the general Bose integral
`∫₀^∞ x^{s-1}/(z^{-1} e^x - 1) dx = Γ(s) PolyLog(s, z)` closes for a symbolic `z`
whenever the `Assumptions` confine `γ' = -z` to `(-1, 1]`. The built-in
assumption engine only discharges syntactic matches (it proves neither `1/z>0`
nor `-1≤-z≤1` from `0<z<1`), so the `-1<γ'≤1` gate is decided by a small **sound
interval-bound prover** over the parameter box read off the `Assumptions`:
interval arithmetic yields a guaranteed enclosure, so the gate never accepts an
inadmissible fugacity (an unbounded or out-of-range `z` simply declines).

Four operational layers extend the table:

- **Monomial substitution** `g(x^k)` (`k≠1`) via `y = x^k`:
  `∫ x^{s-1} g(x^k) = (1/k) ∫ y^{s/k-1} g(y)`, so `Sin[√x]`, `BesselJ[ν,2√x]`,
  `ArcTan[√x]`, `Cos[x²]` reduce to the linear table at `s/k`.
- **Hypergeometric reduction** (applied before Expand, so a cancellation kernel
  is never split): `Erf[u] → u·₁F₁`, `Γ[a]-Γ[a,x] → x^a/a·₁F₁` (lower incomplete
  gamma), and the product `BesselJ[ν,·]² → ₁F₂` (a Mellin convolution closed via
  the `J²` identity rather than a Barnes integral).
- **Parametric differentiation** for `Log[1+λx]^n (1+λx)^{-w₀}`:
  `M = (-1)^n ∂ⁿ_w[λ^{-s} B(s, w-s)]|_{w=w₀}`, strip `-n<Re s<w₀`.
- **`Log[x]^k` weight** (a bare `Log[x]`, distinct from the `Log[1+λx]` kernel):
  since `∂_s x^{s-1} = x^{s-1} Log x`, a `Log[x]^k` factor is the `k`-th
  `s`-derivative of the base transform `M_R(s)`. The open strip carries unchanged
  (`Log^k` is dominated by `x^{±ε}`), so e.g. `∫₀^∞ Log[x]/(1+x²) dx = 0` and
  `∫₀^∞ x Log[x] e^{-x} dx = Γ'(2) = 1-γ`.
- The `pFq` transform is the master kernel — `1F1`, `2F1`, `3F2`, … close
  uniformly (`Hypergeometric1F1`/`2F1` are stored as `HypergeometricPFQ`).

A **Frullani pre-pass** (run on the whole integrand before Expand, since each
half is individually divergent) recognises `(f(a x)-f(b x))/x` and returns
`(f(0⁺)-f(∞)) Log(b/a)` (`a,b>0`): the scale ratio is read structurally and the
pairing verified by the exact identity `(t₁ /. x→ρx)+t₂ = 0`; the boundary values
are the finite limits of `f`. So `∫₀^∞ (e^{-2x}-e^{-5x})/x dx = Log(5/2)` and
`∫₀^∞ (ArcTan(5x)-ArcTan(2x))/x dx = (π/2) Log(5/2)`.

Each application is **gated on its convergence strip** — checked by `Simplify`
(numerically for a numeric `s`, or against the supplied `Assumptions` for a
symbolic `s`), so every result is correct by construction. **When the
assumptions do not prove the strip, the value is returned as a
`ConditionalExpression[value, strip]`** (matching Wolfram) — it collapses to the
bare value once the strip is proved and to `Undefined` if it is refuted. A
provably-violated strip declines. Verification is symbolic; there is **no**
numeric crosscheck (the trig/PolyLog transforms use reflection-formula forms
regular at `s=0`, so e.g. `∫₀^∞ Sin[x]/x dx = π/2` falls out with no limit). A
sum is integrated term by term (each term must converge on its own). Out of
scope — products of three or more transcendental kernels, finite intervals, and
two-sided reductions — return unevaluated, never a wrong value.

 (Ramanujan's canonical example)

 (Debye)
 (Fermi–Dirac)

 (Frullani)
.

### Definite integration (Newton-Leibniz)

`Integrate[f, {x, xmin, xmax}]` gives the definite integral by the
**fundamental theorem of calculus** (`src/calculus/integrate_newton_leibniz.c`,
exposed explicitly as `Integrate\`NewtonLeibniz[f, {x, xmin, xmax}]` and as
`Method -> "NewtonLeibniz"`).  The mechanism:

1. Antidifferentiate `f` with the ordinary indefinite cascade to get `F`.  If
   no closed antiderivative exists, the definite integral is left unevaluated
   — never assigned a wrong value.
2. Locate the real poles of the **integrand** `f` strictly inside
   `(xmin, xmax)` via `Integrate\`SingularPoints` (roots of
   `Denominator[Together[f]]`).  A continuous `f` has an antiderivative that is
   pole-free where `f` is, so only `f`'s own poles bound the segments.
3. Split `[xmin, xmax]` at those poles and form the telescoping sum
   `Σ (F(p_{i+1}⁻) − F(p_i⁺))`, evaluating each boundary through the `Limit`
   engine: a plain limit at infinite endpoints, a one-sided limit
   (`Direction -> "FromBelow"/"FromAbove"`) at interior poles, and direct
   substitution at ordinary finite endpoints.  Improper integrals thereby
   acquire their correct finite value; a genuinely divergent integral emits
   `Integrate::idiv` and is left unevaluated (matching Mathematica).

`Integrate[f, {x, a, b}, {y, c, d}, …]` gives the iterated multiple integral,
reduced innermost-first (the last spec is the inner integral), so an inner
bound may depend on an outer variable.

`Integrate\`SingularPoints[expr, {x, a, b}]` returns the sorted list of the
interior real poles used for the split — exposed for inspection and reuse.

**Cauchy principal value.** `Integrate[f, {x, a, b}, PrincipalValue -> True]`
computes the Cauchy principal value across an interior pole. It is defined only
for **odd-order** poles (the integrand changes sign across the pole, so the
symmetric one-sided divergences cancel); the value is then `Re[F(b) − F(a)]`
with the real-branch antiderivative — the branch-cut crossing at each pole
contributes a purely imaginary `I Pi` (residue) that `Re` removes. An
**even-order** interior pole has no principal value and emits `Integrate::idiv`.
With no interior pole the option is a no-op. Examples: `∫₋₁¹ dx/x = 0`,
`∫₀³ dx/(x−2) = −Log 2`, `∫₀² x/(x²−1) dx = ½ Log 3`.

**Continuous branch-form antiderivatives.** Many continuous periodic integrands
(e.g. `1/(2 + Cos[x])`) antidifferentiate to a Weierstrass branch form that is
*already continuous* — an `ArcTan[Tan[x/2] …]` term whose jump is exactly
cancelled by an accompanying `Floor` step — so `F(b) − F(a)` by direct
substitution is the correct value with no interior split.  To stay correct when
an antiderivative could instead carry a genuine *uncorrected* branch jump, a
result built from a step head (`Floor`, `Sign`, …) or an inverse-trig node over
a pole-bearing trig head is accepted only when a numeric `NIntegrate`
cross-check agrees; otherwise the integral is left unevaluated rather than
risking a wrong value.

### Complex line / contour integration

When a `{x, a, b}` spec has a **non-real endpoint**, or when a spec lists more
than two points `{x, z0, z1, …, zn}`, `Integrate` evaluates the **contour
integral** of `f` along the straight segments `z0 → z1 → … → zn` in the complex
plane (`src/calculus/integrate_line.c`, exposed explicitly as
`Integrate\`LineIntegral[f, {x, z0, …, zn}]`).  Real-endpoint two-point specs
still go through the real-axis Newton-Leibniz path above.

Each segment `a → b` is parametrised by a **real** parameter,
`γ(t) = a + t (b − a)`, `t ∈ [0, 1]`, which reduces the complex problem to the
real machinery already in place:

1. Antidifferentiate `f` in `x` to get `F` (bail if unknown).
2. **On-path singularities** become real roots `t* ∈ (0, 1)` of
   `Denominator[Together[f(γ(t))]]` — a singularity strictly on the contour
   makes the integral divergent (`Integrate::idiv`, left unevaluated).
   `Integrate\`PathSingularPoints[f, {x, z0, …, zn}]` returns those points.
3. The segment value is the continuous change of `F` along the segment, with
   endpoint values taken as **real one-sided limits in `t`** when substitution
   is singular (so a complex-ray approach is a real one-sided limit the `Limit`
   engine can take).  For rational integrands whose antiderivative is a sum of
   logarithms / inverse-tangents of **affine** arguments, the branch-correct
   value is recovered by combining each into a single principal `Log` of a ratio
   `Log[(u(b))/(u(a))]` — exact because a straight segment subtends an angle
   `< π` at any point off it, so closed-contour residues come out exactly.
4. Every segment value is numerically cross-checked against a complex quadrature
   of `f(γ(t)) γ'(t)`; an uncorrectable branch crossing leaves the integral
   unevaluated rather than returning a wrong branch.

### Contour / residue-theorem definite integration

For the classical families of **improper** and **periodic** real integrals that
complex analysis dispatches by summing residues over the poles enclosed by a
standard contour, `Integrate[f, {x, a, b}]` runs a residue-theorem method
**before** Newton-Leibniz (also reachable as `Method -> "Residue"` or
`Integrate`ContourResidue[f, {x, a, b}]`).  Each answer is
**correct-by-construction**: once a family's structural gates hold, the residue
theorem gives the exact value, so there is *no* numeric quadrature crosscheck.
The only post-hoc gate is a self-consistency check that the residue sum closed to
a scalar (no surviving `x`/`Root`) and — for the real-valued families — that its
imaginary part vanishes (a residual `Im` would betray a mis-classified pole); a
failure returns unevaluated and the Newton-Leibniz path takes over.  Four
recognizers:

- **Rational on `(-∞, ∞)`** — `f = P/Q` with `deg Q ≥ deg P + 2` and **no real
  pole**: value `= 2 π i · Σ Res` over the poles in the upper half-plane.
- **Fourier / Jordan on `(-∞, ∞)`** — `f = R(x) · K` with
  `K ∈ {Cos[a x], Sin[a x], Exp[I a x]}`, `R` rational with `deg`-drop `≥ 1`,
  `a` a nonzero real: with `J = 2 π i · Σ_UHP Res[R Exp[I a x]]` (for `a > 0`),
  `∫ R Cos = Re[J]` and `∫ R Sin = Im[J]`. The bare complex-exponential kernel is
  recognised in both spellings — `Exp[I a x]` and the evaluator-normalised
  `Power[E, I a x]` — and returns `J` directly (`∫ R Exp[I a x] = J`); for
  negative `a` the lower half-plane is closed. Real-exponent `Exp[a x]` is *not*
  a Fourier kernel (it routes to the rectangular family). `Cos`/`Sin` are closed
  via `Re`/`Im` of the single decaying exponential, not split into two.
- **Rational-in-`{Sin, Cos}` over a full period** `(0, 2π)` or `(-π, π)` — via
  `z = Exp[I x]` on the unit circle: value `= 2 π i · Σ Res` over the poles
  inside the unit disk.
- **Removable axis singularity (Fourier)** — a **simple** real-axis pole of `R`
  at which the kernel vanishes (so `f = R·K` is analytic there, e.g. `Sin[x]/x`
  at `0`) contributes a **half residue** `π i · Res`, the indented-contour value,
  giving `∫ Sin[x]/x = π`.  A *genuine* axis pole (kernel nonzero there, e.g.
  `Cos[x]/x`) makes the ordinary integral diverge and returns unevaluated —
  plain `Integrate` does not compute a principal value.

A one-line **half-line** add-on covers even integrands:
`∫₀^∞ f = ½ ∫₋∞^∞ f`.  Three further recognizers handle branch-cut and
symbolic-exponent contours:

- **Keyhole / Mellin on `(0, ∞)`** — a branch power times a rational function,
  `f = x^p R(x)` with `p` non-integer (so `s = p + 1 ∉ ℤ`): value
  `= −π · Σ_k Res[ z^{s-1} R(z), z_k ] · e^{−iπs} / sin(π s)` over the poles
  `z_k` of `R`, with `z^{s-1}` on the branch `arg ∈ (0, 2π)`.  The residue is of
  the **full** integrand `z^{s-1} R(z)` — for a simple pole this equals
  `z_k^{s-1} Res(R, z_k)`, but for an order-≥2 pole it does not (a pure double
  pole has `Res(R)=0`), so poles of any order are handled by shifting `w=z−z_k`
  and expanding the analytic factor `(1+w/z_k)^{s-1}`.  The `e^{−iπs}/sin(π s)`
  prefactor is the exact reduction of the keyhole jump `1/(1 − e^{2πi s})`,
  landing a numeric `s` on an algebraic multiple of `π` (e.g.
  `∫₀^∞ x^{1/3}/(x²+1) = π/√3`, `∫₀^∞ √x/(1+x)² = π/2`).  Requires
  `0 < Re(s) < deg Q − deg P`.
- **Sector on `(0, ∞)`** — `f = x^m/(c + x^n)` with the exponent `n` possibly a
  **symbolic parameter**: the wedge of angle `2π/n` gives
  `(π/n) c^{s/n − 1} csc(π s/n)`, `s = m + 1`.  This is the one family admitting a
  symbolic `n` (the keyhole cannot enumerate `n` poles), powering
  `Integrate[1/(1 + x^n), {x, 0, ∞}, Assumptions -> n > 1] = (π/n) csc(π/n)`.
- **Rectangular / quasi-periodic on `(-∞, ∞)`** — `f = Exp[c x] R(Exp[x])`
  (period `2πi`): reduced to the keyhole core by `w = Exp[x]`
  (`∫_{-∞}^∞ f dx = ∫₀^∞ f(Log w)/w dw = ∫₀^∞ w^{c-1} R(w) dw`), so e.g.
  `Integrate[Exp[a x]/(Exp[x]+1), {x, -∞, ∞}, Assumptions -> 0 < a < 1] = π csc(π a)`.

**Assumptions and symbolic parameters.**  An `Integrate[f, {x, a, b},
Assumptions -> …]` option lets the residue families evaluate integrals whose
parameters are symbolic (`a > 0`, `0 < a < 1`, `n > 1`, …).  The recognizers
classify a parameter-dependent pole or kernel frequency by reading its sign at a
single generic point of the region the assumptions pin (a sign-consistent
instantiation), while the residue arithmetic stays fully symbolic — so the
closed form is still **correct by construction, with no numeric crosscheck**.
Convergence/applicability gates (`n > m + 1`, `0 < s < deg`-drop, `c > 0`) are
verified against the assumption-**guaranteed** interval bounds, not the sample
point, so an under-constrained problem (e.g. only `n > 0` for the sector family)
is refused rather than guessed; a parameter the assumptions leave two-sided
unbounded is likewise refused.  Radical pole locations are `PowerExpand`-cleaned
under all-positive parameters (`Sqrt[-4 a²] → 2 I a`) so a rational answer closes
to `π/a` rather than a `Sqrt[-4 a²]` surface, and real-parameter conjugation
uses `I → −I` (the symbolic `Conjugate` would not reduce).

The residue sums are closed with `RootReduce` (algebraic families) or the
`Conjugate` identities for `Re[J]`/`Im[J]` (the Fourier family); a value that
still contains an unreduced `Root`, or a surviving imaginary part, is treated as
a mis-fire and returned unevaluated.  Negative controls such as
`Integrate[1/(1 + x^3), {x, -Infinity, Infinity}]` and
`Integrate[1/(x^2 - 1), {x, -Infinity, Infinity}, Method -> "Residue"]`
(genuine real-axis poles),
`Integrate[Cos[x]/x, {x, -Infinity, Infinity}, Method -> "Residue"]`
(genuine axis pole, kernel nonzero),
`Integrate[1/Sqrt[1 + x^4], {x, -Infinity, Infinity}]` (branch point, not
rational), and `Integrate[1/(2 + Cos[x]), {x, 0, Pi}]` (not a full period) all
stay unevaluated.  The keyhole/Mellin, sector and rectangular families
described above extend the reach to branch-cut and symbolic-exponent contours
a log-keyhole (`∫₀^∞ Log[x] R(x)`) with symbolic on-circle poles remains out of
scope (it needs assumption-aware `Arg`/`Log` branch reasoning that Mathilda does
not yet have).

The `Integrate`` package also exposes the lower-level helpers
`Integrate`HermiteReduce`, `Integrate`IntegratePolynomial`,
`Integrate`BronsteinRational` (the explicit form),
`Integrate`IntRationalLogPart` (Phase 2's LRT computation),
`Integrate`RischNorman` (Bronstein pmint), `Integrate`LinearRadicals`
(linear-radical substitution), `Integrate`QuadraticRadicals`
(quadratic-radical Euler substitution), `Integrate`LinearRatioRadicals`
(linear-fractional / Möbius radical substitution), `Integrate`Weierstrass`
(continuous `Tan[x/2]` / `Tanh[x/2]` substitution), `Integrate`CRCTable`
(table lookup), `Integrate`Undefined` (unknown-function integration,
Roach §1.7), and the unit-test helpers `Integrate`Helpers`Content`,
`...`Primitive`, `...`Monic`, `...`LeadingCoefficient`,
`...`SquareFree`, `...`ExtractConstants`, `...`ApartList`.  All are
`Protected`; the BronsteinRational helpers additionally have
`ReadProtected`.

### Integrate`CRCTable

`Integrate`CRCTable[f, x]` looks `f` up in the CRC Standard Mathematical
Tables (31st ed., 600+ formulas).  The rules live in
`src/internal/CRCMathTablesIntegrals.m`, internally on the head
`IntegrateTable`; `Integrate`CRCTable` is a thin wrapper around it.
The .m file is `Get`-loaded on the first invocation of the CRCTable
stage rather than at startup, so sessions that never call `Integrate`
pay nothing for the table.

Every recursive rule in the table carries an `IntegerQ[index] && index
> base` (or `< base`) guard sufficient to guarantee termination via
first-principles analysis of the reduction direction.  Without these
guards rules such as Formula 49 (the `1/(x^2 - c^2)^n` reduction)
would diverge on negative or non-integer `n`.  As defence-in-depth,
the C dispatcher caps CRC-rule recursion depth at 256 levels and
emits `Integrate`CRCTable::depth` rather than locking up on any rule
that escapes the audit.

The table covers the inverse-trig and inverse-hyperbolic families
(Formulas 427–464 and their `427h`–`464h` analogs) among many others.
Rules whose denominator carries a *squared* coefficient
(`Sqrt[1 ∓ a^2 x^2]`, `(1 ± a^2 x^2)`) bind that coefficient linearly as
`a_` via `c_ + a_. x_^2` and recover the linear coefficient as `Sqrt[±a]`
on the RHS (with a `Condition` linking it to the numerator coefficient) —
the matcher does not invert a square, so `a_^2` in a pattern will not
bind. Some more elaborate multi-argument `/;`-guarded rules still do not
fire; this is a separate issue tracked under the matcher work.

### Integrate`Undefined

`Integrate`Undefined[f, x]` integrates expressions that are rational in
**undefined functions** `u[x]` and their derivatives, following Kelly
Roach, "Indefinite and Definite Integration" (1992), §1.7 ("Undefined
Functions").  Each undefined function value `u[g]` and its derivative
tower `u'[g], u''[g], …` is treated as a differential-field generator
the integrand is reduced by recognising integration-by-parts /
total-derivative structure in the top generator.  A single inner call to
the rational integrator over a substituted generator symbol subsumes
Roach's polynomial, fraction, and log parts (so `1/u -> Log[u]`,
`1/u^2 -> -1/u`, and `a/(a^2+u'^2) -> ArcTan` all fall out).  Composite
arguments are handled via the chain rule (`g' = D[g, x]`).  A logarithm of
an unknown-function expression, `Log[eta]`, is itself recognised as a
transcendental generator and reduced by parts
(`Integrate[C L + D, x] = G L - Integrate[G L' - D, x]` with
`G = Integrate[C, x]`, `L' = eta'/eta`), with self-referential resolution
for perfect powers (`Integrate[Log[f] f'/f, x] = Log[f]^2/2`).  The stage
is gated to only run when `f` contains an undefined-function derivative or
such a logarithm; genuinely non-elementary integrands (e.g. `f'[x] g'[x]`,
`f'[x]^2`) are left unevaluated, with a cycle guard preventing the by-parts
recursion from looping.  `Protected`, `ReadProtected`.

Known limitations: transcendental generators other than `Log` (e.g.
`ArcTan[eta]`, `Exp[eta]` with `eta` containing an unknown function) are
not yet recognised; nested unknown arguments `f[g[x]]` (with `g` itself
undefined) are deferred.

### Integrate`DerivativeDivides

`Integrate`DerivativeDivides[f, x]` integrates **by substitution** — the
classical "derivative-divides" method (Moses' SIN, Maxima's `diffdiv`).  It
recognises an integrand of the shape `f(x) = c · h(u(x)) · u'(x)`, reduces the
problem to `Integrate[h[u], u]`, and back-substitutes `u -> u(x)`.
Implemented in `src/calculus/integrate_derivdivides.c`.

Candidate kernels `u(x)` are every distinct subexpression of `f` that depends
on `x`, except `x` and `f` themselves (the C analogue of
`Cases[Union[Level[f, {0, Infinity}]], e_ /; !FreeQ[e, x]]`).  Two
complementary strategies are tried per kernel, each followed by a
**verification gate** `PossibleZeroQ[D[result, x] - f] === True` (numeric
sampling; as of 2026-06-09 the former rigorous `Simplify` confirm was dropped —
it cost ~1.1 s of `trigrat` normalization on the radical-trig cases for no gain
over the sampler):

1. **Direct quotient** — `q = Cancel[Together[f / D[u(x), x]]]`, then
   `q /. u(x) -> u`; accepted when free of `x`.  Cheap, emits no diagnostics,
   handles transcendental compositions (`x Exp[x^2]`), and **selects the
   correct radical branch inherently** (no squaring).  Tried first.

2. **Eliminate / Solve** — builds the differential relation
   `Eliminate[{Dt[y] == f Dt[x], u == u(x), Dt[u == u(x)]}, {x, Dt[x]}]`,
   solves for `Dt[y]`, and reduces each branch with `Factor //@ `,
   `PowerExpand` and `Cancel[. / Dt[u]]`.  It closes integrands the direct
   strategy cannot — including **radical substitutions** like
   `u = Sqrt[Tan[x]]` (reducing `Integrate[Sqrt[Tan[x]], x]` to the rational
   integral `2 u^2/(1 + u^4)`), and cases where `Cancel` canonicalises
   `1/Cos[x]` to `Sec[x]`, breaking the syntactic substitution.  Because the
   algebraisation can square radicals and invert functions, `Solve` returns
   several branches; the verification gate is what **selects the branch that
   differentiates back to `f`**.  This strategy runs in the Automatic cascade as
   well as under the explicit method (its Eliminate `::ifun` / `::alg`
   diagnostics are muted while the integrator drives it).  Because it is
   heavyweight (~0.1–1 s per kernel), as of 2026-06-09 it runs **only on the
   outermost integrand** — reduced sub-integrals are finished by the direct
   strategy and the rest of the cascade.

The reduced integral re-enters the full `Integrate`, so substitutions compose.
Three guards keep the recursion finite and cheap: an **integrand memo** that
short-circuits any integrand (canonicalised by renaming the integration variable
to a fixed sentinel) already attempted in the current top-level descent — this
breaks circular substitution chains and collapses overlapping subproblems that
would otherwise fan out exponentially (e.g. `Integrate[x Sin[x^2], x]`); the
**outermost-only** restriction on the Eliminate/Solve search above; and a hard
depth backstop (8) with per-call fresh substitution symbols.  Strict: returns
unevaluated when no substitution closes the integral.  `Protected`,
`ReadProtected`.

Known limitations: kernels must appear **literally** in `f` (so `Tan[x]`,
which Mathilda keeps atomic rather than `Sin[x]/Cos[x]`, exposes no `Cos[x]`
kernel — such integrands are handled by RischNorman instead); the reduced
integral must itself close under the other methods.

### Integrate`LinearRadicals

`Integrate`LinearRadicals[f, x]` integrates a **rational function of `x` and
radicals `(a x + b)^(m/n)` that share one linear argument** by the classical
rationalising substitution `u = (a x + b)^(1/n)`, where `n = LCM` of the radical
denominators.  Implemented in `src/calculus/integrate_linrad.c`.

It first scans `f` for every `Power[base, p/q]` with `q > 1`: each `base` must be
a degree-1 polynomial in `x` and all must be the **same** linear form `a x + b`
(distinct bases — e.g. `Sqrt[x] + Sqrt[x+1]` — are rejected, as are radicals of a
non-linear argument such as `Sqrt[x^2+1]`).  With `x = (u^n - b)/a` and
`dx = (n/a) u^(n-1) du` the integrand becomes

a **rational function of `u`** that `Integrate`BronsteinRational` closes exactly.
The radical substitution reuses `poly_subst_radical_to_gen` (shared with the
algebraic-factoring path); after the substituted integrand is confirmed rational
in `{x, u}`, the reduced integral re-enters the full `Integrate` and the result
is back-substituted `u -> (a x + b)^(1/n)`.  Unlike the heuristic methods, the
result is **not** put through a `Simplify[D[result, x] - f] === 0` gate: the
substitution is an exact bijection that introduces no branch issues, so the
antiderivative is correct by construction once the rational sub-integral closes.
(Skipping the gate also avoids a prohibitively expensive — and on integrands with
symbolic parameters, non-terminating — `PossibleZeroQ`/`Simplify` over the
symbolic-radical residue.)  A depth guard (8) and per-call fresh substitution
symbols keep the recursion finite and collision-free.  Strict: returns
unevaluated when `f` is not of this form or the reduced integral does not close.
`Protected`, `ReadProtected`.

### Integrate`QuadraticRadicals

`Integrate`QuadraticRadicals[f, x]` integrates a **rational function of `x` and
square roots `(a x^2 + b x + c)^(m/2)` that share one quadratic argument** by a
classical **Euler substitution**.  Implemented in
`src/calculus/integrate_quadrad.c`.

It scans `f` for every `Power[base, p/q]` whose `x`-dependent `base` is a
degree-2 polynomial: each must be a **square root** (`q == 2`) and all must be
the **same** radicand `rad = a x^2 + b x + c` (a cube root such as
`(x^2+1)^(1/3)`, a radical of a cubic such as `Sqrt[x^3+1]`, or distinct
radicands are rejected; a purely linear radical belongs to
`Integrate`LinearRadicals`).  Substituting the radicals to a fresh symbol `y`
(`rad^(p/2) -> y^p`, via `poly_subst_radical_to_gen`) must leave a rational
function of `{x, y}`.

To keep the antiderivative **real-valued**, exactly **one** substitution is
chosen by the sign of the leading coefficient `a` (and, when `a < 0`, of the
discriminant `b^2 - 4 a c`) — the routine does *not* try all three Euler forms:

| Condition | Euler substitution | Image `y = Sqrt[rad]` |
|-----------|--------------------|-----------------------|
| `a > 0` | first | `Sqrt[a] x + u` |
| `a < 0` and `b^2 - 4 a c > 0` | third | `(x - alpha) u`, with `alpha` a real root |
| `a` symbolic | first (best-effort `a > 0` branch) | `Sqrt[a] x + u` |

For numeric `a < 0` a real radicand requires real roots, so the third
substitution subsumes the second (the `c > 0` form).  Each branch yields
`x = X(u)`, `dx = X'(u) du` and the radical image above; the rationalised
integrand `f dx /. {y -> ..., x -> X}` re-enters the full `Integrate` and the
result is back-substituted `u -> U(x)`.  As with `Integrate`LinearRadicals`, the
Euler substitution is an exact bijection on the relevant domain, so the result
is **not** put through a `Simplify[D[result, x] - f] === 0` gate — it is correct
by construction once the rational sub-integral closes.  A depth guard (8) and
fresh per-call substitution symbols keep the recursion finite.  Strict: returns
unevaluated when `f` is not of this form or the reduced integral does not close.
`Protected`, `ReadProtected`.

### Integrate`LinearRatioRadicals

`Integrate`LinearRatioRadicals[f, x]` integrates a **rational function of `x`
and radicals `((a x + b)/(c x + d))^(m/n)` that share one linear-fractional
(Möbius) argument** by a rationalising substitution.  Implemented in
`src/calculus/integrate_linratiorad.c`.

It scans `f` for every `Power[base, p/q]` (`|q| > 1`) whose `base` depends on
`x`; all such bases must be **structurally identical** and `n = LCM` of the
radical denominators (distinct bases are rejected).  Unlike
`Integrate`LinearRadicals` / `Integrate`QuadraticRadicals`, the base is **not**
required to be a polynomial — it is the ratio
`Times[a x + b, Power[c x + d, -1]]`, which is exactly what partitions this
method from those two (their scans demand a polynomial base, so a ratio radical
falls through to here).  The shared base is canonicalised with
`Cancel[Together[.]]` and the Möbius coefficients read off its numerator
(`a, b`) and denominator (`c, d`); the denominator must be genuinely linear
(degree 1 — a constant denominator is `Integrate`LinearRadicals`' job) and the
map non-degenerate (`a d - b c != 0`).

With `m = n` the substitution `u = ((a x + b)/(c x + d))^(1/m)` inverts the
Möbius map,

and the radicals rewrite to `u` (`r^(p/q) -> u^(p m/q)`, via
`poly_subst_radical_to_gen`).  The result must be rational in `{x, u}`; the
rationalised integrand re-enters the full `Integrate` and the antiderivative is
back-substituted `u -> ((a x + b)/(c x + d))^(1/m)`.  As with the other radical
stages, the Möbius substitution is an exact bijection on the relevant domain, so
the result is **not** put through a `Simplify[D[result, x] - f] === 0` gate — it
is correct by construction once the rational sub-integral closes.  A depth guard
(8) and fresh per-call substitution symbols keep the recursion finite.  Strict:
returns unevaluated when `f` is not of this form or the reduced integral does
not close.  `Protected`, `ReadProtected`.

### Integrate`Weierstrass

`Integrate`Weierstrass[f, x]` integrates a **rational function of the
trigonometric kernels** `Sin/Cos/Tan/Cot/Sec/Csc[x]` — or the **hyperbolic
kernels** `Sinh/Cosh/Tanh/Coth/Sech/Csch[x]` — by the continuous Weierstrass
substitution of Jeffrey & Rich (*The Evaluation of Trigonometric Integrals
Avoiding Spurious Discontinuities*, ACM TOMS 20(1), 1994).  Added 2026-06-09.

Algorithm: substitute `u = Tan[x/2]` (`u = Tanh[x/2]` for hyperbolic), turning
`f` into a rational function of `u`; integrate that (recursing through
`Integrate`BronsteinRational`); back-substitute; and — for the trigonometric
case — add the secular correction `K Floor[(x - b)/p]` (`b = Pi`, `p = 2 Pi`)
that removes the spurious jump discontinuities the classical substitution
introduces at the poles of `Tan[x/2]` (odd multiples of `Pi`).  The jump `K` is
the difference of the one-sided limits of the `u`-antiderivative at `±Infinity`
if that limit diverges (a *genuine* singularity of the integrand) no correction
is applied.  A `TrigExpand` pre-pass reduces multiple/sum-angle arguments
(`Cos[2 x]`, `Cosh[x] Cosh[2 x]`, ...) to kernels of the bare variable.

### Hyperbolic case needs no `Floor` correction

`Tanh[x/2]` is a smooth,
strictly monotone bijection `R -> (-1, 1)` with no poles, so the substitution
introduces no spurious discontinuity — the back-substituted antiderivative is
already continuous (genuine singularities such as `Cosh[x] = 2` are real poles
and are correctly left in place).

The substitution is an exact identity, the rational sub-integral is closed by a
verified integrator, and `Floor' = 0` almost everywhere, so the result is
**correct by construction** — no differentiate-back gate is applied (and the
`Floor` term would defeat symbolic `D` anyway).  In the Automatic cascade only
genuine rational integrands (a kernel in a denominator) are intercepted, so
polynomial trig such as `Integrate[Sin[x], x]` keeps its cleaner table form; the
explicit `Method -> "Weierstrass"` has no such gate.  Strict: returns unevaluated
when `f` is not a rational function of the trig/hyperbolic kernels of `x` (e.g.
`x` outside a kernel, a kernel of a nonlinear argument, mixed trig + hyperbolic,
or a radical of a kernel).  `Protected`, `ReadProtected`.

### InterpolatingFunction integrands

`Integrate[InterpolatingFunction[...][x], x]` returns the antiderivative as a
fresh applied `InterpolatingFunction[...][x]`, mirroring how
`D[InterpolatingFunction[...][x], x]` differentiates such objects.
Differentiation only bumps the object's derivative-order annotation; the
per-window evaluation kernels evaluate derivatives of order `>= 0` only and
cannot produce an antiderivative, so integration instead builds a genuinely new
interpolant (`src/calculus/integrate_interp.c`):

1. Read the grid x-coordinates from the object's stored table.
2. Sample the original function values `y_i = ifun[x_i]`.
3. Accumulate the antiderivative node values `F_0 = 0`,
   `F_i = F_{i-1} + Integrate_{x_{i-1}}^{x_i} ifun` by 5-point Gauss-Legendre
   quadrature (exact through degree 9 — i.e. the default/Spline/Hermite
   piecewise-cubic interpolants; very high explicit `InterpolationOrder` panels
   incur a small quadrature error).
4. Build a Hermite `InterpolatingFunction` through `{{x_i}, F_i, y_i}` — because
   the antiderivative's exact derivative is the original function, supplying
   `F'(x_i) = y_i` makes `D[Integrate[ifun[x], x], x]` round-trip to `ifun[x]`.

Only the 1-D, direct case (the applied argument is the integration variable
itself) is reduced; `Integrate[ifun[g[x]], x]` is not generally expressible as
an `InterpolatingFunction` and is left to the cascade above. The construction
also handles the derivative-annotated objects produced by `D` (integrating the
sampled derivative recovers the lower-order antiderivative). Computations use
machine doubles.

## Algorithm

integrate.c

```text
`Integrate[f, x]` System` dispatcher.  Cascades through three
```

stages and supports an explicit `Method -> "..."` option:

```text
  1. Integrate`BronsteinRational   — polynomial / rational integrands
  2. Integrate`RischNorman         — parallel-Risch (Bronstein pmint)
  3. Integrate`CRCTable            — CRC integral table (lazy-loaded)
```

Method values: "Automatic" (default, full cascade), "BronsteinRational", "RischNorman", "CRCTable" (strict passthrough, no fallback).

The CRC table is large and most sessions never need it, so its .m file is Get-loaded on first invocation of try_crctable() rather

```text
than at startup.  See LAZY_LOAD_CRC below.
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| NI 50-digit Gaussian | 6.87 s | 2.11 s | 1.1 s |
| integrate Sin[x] Exp[x] | 6.62 s | 0.726 s | 5.47 s |
| integrate Exp[x]/x (non-elementary) | 5.38 s | 0.318 s | 33.5 s |
| integrate Tan[x]^3 | 4.07 s | 0.196 s | 3.62 s |
| NI 2-D ridge | 2.54 s | 43.5 s | 0.685 s |
| Discriminant of deg 20 | 2.51 s | 0.068 s | 0.182 s |

## Implementation notes

**Algorithm.** `Integrate[f, x]` computes an *indefinite, single-variable*
antiderivative only; there is no constant of integration and no definite-bound
support, and results are not guaranteed to be simplified. The dispatcher
`builtin_integrate` (src/calculus/integrate.c) routes through a fixed cascade of
context-qualified sub-integrators, each of which returns the antiderivative on
success or comes back as its own unevaluated head (e.g.
`Integrate\`BronsteinRational[...]`) to signal "fall through". The Automatic
cascade order is:

1. **Undefined-function integrator** (`try_undefined` → src/calculus/integrate_unknown.c),
   Roach 1992 §1.7: integrands rational in unknown `u[x]` and their derivatives
   (e.g. `x f'[x] + f[x] -> x f[x]`). Cheaply gated.
2. **Bronstein rational integration** (`try_rational` → src/calculus/intrat.c),
   gated by `PolynomialQ`/rational-function tests. This is the mathematically
   substantial stage: polynomial part by term-wise power rule, then **Mack's
   linear Hermite reduction** (`HermiteReduce`) to split off the rational
   (algebraic) part of the antiderivative, then the **Lazard–Rioboo–Trager
   log part** (`IntRationalLogPart`) built on a subresultant PRS, with Rioboo's
   recursive `LogToAtan`/`LogToReal` conversions producing real-elementary
   `Log`/`ArcTan` forms (a `RootSum`-based `NaiveLogPart` is the universal
   fallback). A pre-Hermite derivative-recognition fast path catches `c·D'/D^k`.
3. **Radical substitutions**: linear-radical (`try_linrad`), quadratic-radical
   Euler substitution (`try_quadrad`), and linear-ratio Möbius substitution
   (`try_linratiorad`).
4. **Derivative-divides** substitution (`try_derivdivides`): recognises
   `c·h(u(x))·u'(x)` and reduces to `Integrate[h[u], u]`.
5. **Risch–Norman heuristic** (`try_risch` → src/calculus/intrischnorman.c),
   Bronstein's *pmint*: the parallel-Risch ansatz for transcendental
   integrands (exp/log/trig via `tan` rewriting). Builds a candidate monomial
   set, a vector field with `splitFactor`/deflation, a linear system solved by
   `RowReduce`, plus a log-candidate sum; retries over `K = I`.
6. **CRC integral table** (`try_crctable`): a large rule set lazily
   `Get`-loaded from src/internal/CRCMathTablesIntegrals.m on first use only.
   Recursion-capped at `MAX_CRC_DEPTH` (256); a leaked internal
   `IntegrateTable[...]` head anywhere in the result counts as a miss.

An explicit `Method -> "<name>"` option (parsed by `parse_method_option`)
bypasses the cascade and dispatches to one stage strictly (no fallback).
Inexact integrands are force-rationalised by `common_rationalize_input`
(min-precision-aware, ½-ulp fallback), integrated exactly, then numericalised
back via `common_numericalize_result` for inexact-in/inexact-out semantics.
Applied 1-D `InterpolatingFunction` objects integrate to a fresh antiderivative
`InterpolatingFunction` before any rationalisation (`integrate_interp`).

**Data structures.** Everything is a Mathilda `Expr*` tree; the rational stage
leans on Mathilda's own polynomial builtins (`Together`, `Numerator`,
`Denominator`, `SubresultantPolynomialRemainders`, `PolynomialGCD`). Stages
communicate by evaluating context-qualified heads (`call_stage`) and
distinguishing success from passthrough by inspecting the result head
(`result_is_unresolved`, `result_contains_head`).

**Complexity / limits.** Indefinite, single-variable only. The rational stage
is complete for rational integrands (Hermite + LRT); transcendental coverage is
heuristic (pmint may give up) backed by a finite table. No definite integrals,
no multivariate integration, no constant of integration.

- `Protected`, `Listable`.
- Eleven-stage dispatch cascade (`DerivativeDivides`, `LinearRadicals`,
  `QuadraticRadicals` and `LinearRatioRadicals` added 2026-06-06; `Weierstrass`
  added 2026-06-09; `ChebychevAlgebraic` and `GoursatAlgebraic` added 2026-06-29):
  `Integrate[f, x]` (Method -> Automatic, default) tries each subroutine in
  order and returns the first non-`NULL` result:
  1. `Integrate\`Undefined[f, x]` — when `f` contains an undefined-function
     derivative (e.g. `f'[x]`); see below.
  2. `Integrate\`BronsteinRational[f, x]` — when `PolynomialQ[f, x] ||
     rationalQ[f, x]`.
  3. `Integrate\`LinearRadicals[f, x]` — rational functions of `x` and radicals
     `(a x + b)^(m/n)` of one shared linear argument; rationalised by
     `u = (a x + b)^(1/n)`.
  4. `Integrate\`QuadraticRadicals[f, x]` — rational functions of `x` and square
     roots `(a x^2 + b x + c)^(m/2)` of one shared quadratic argument;
     rationalised by a single real-valued Euler substitution.
  5. `Integrate\`LinearRatioRadicals[f, x]` — rational functions of `x` and
     radicals `((a x + b)/(c x + d))^(m/n)` of one shared linear-fractional
     argument; rationalised by `u = ((a x + b)/(c x + d))^(1/n)`.
  6. `Integrate\`ChebychevAlgebraic[f, x]` — Chebychev binomial differentials
     `x^p (a x^r + b)^q` (`p, q, r` rational, `a, b` free of `x`).  Elementary
     iff one of `q`, `(p+1)/r`, `q+(p+1)/r` is an integer (Chebychev's theorem),
     with substitutions `x = u^N` (Type I), `u^s = a x^r + b` (Type II), or
     `u = x^r` then `t^s = (a u + b)/u` (Type III) that rationalise `f`.
     Recognition is a single structural scan, so it runs ahead of
     `DerivativeDivides`'s Eliminate/Solve search.  Non-elementary binomials
     return `NULL` (the cascade falls through to later methods).
  7. `Integrate\`GoursatAlgebraic[f, x]` — pseudo-elliptic integrands
     `F(x) R(x)^q` (`F` rational, `R` a polynomial) with `q` any rational of
     reduced denominator `2`, `3`, or `4` by Goursat's algorithm and its
     cube-/fourth-root generalisations (Blake 2026).  The exponent is split
     `R^q = R^k R^(-p)` with radical order `p in {1/2, 1/3, 2/3, 1/4, 3/4}` and
     the integer `R^k` absorbed into `F`, so positive-power radicals such as
     `(1-x^3)^(1/3)/x` are handled, not only radicals already in a denominator.
     A Mobius automorphism cycling the roots of `R` splits the integrand into
     eigencomponents that descend to genus-0 curves when the elementarity
     criterion holds (`p=1/2`: Klein four-group `V4`, trivial projection
     vanishes; `p=1/3,2/3`: order-3 cycle; `p=1/4,3/4`: order-4 cycle on
     harmonic roots).  For `p=1/2` with `R` a cubic carrying the `t^3-1` higher
     symmetry, when `V4` declines a Section-4 (Goursat 1887) period-3 reduction
     is tried: an order-3 Mobius `S` fixes one ramification point and cycles the
     other three, and the integral is elementary when `F` is a non-trivial
     period-3 character `F(S) = Exp[2 Pi I/3] F` (so `(x-1)/((x+2) Sqrt[x^3-1])`
     integrates).  For `p=1/3` with `R` a cubic, when the order-3 eigendescent is
     obstructed (`H1 != 0`) but `F` has a pole at a non-branch point, a
     constructive third-kind logarithmic-part reduction is tried instead: the
     antiderivative is
     `C Sum_{j=0..2} Exp[4 Pi I j/3] Log[R^(1/3) - Exp[2 Pi I j/3] (lead R)^(1/3) x]`,
     so the parametric family
     `(2-(k+1)x)/((1-(k+1)x) (x(1-x)(1-k x))^(1/3))` integrates over `Q(k)`.
     The rational reductions are integrated recursively and
     back-substituted.  Obstructed (genuinely elliptic) integrands,
     non-harmonic quartics, and the cross-character A4 cases (Section 5, e.g.
     `t/((t^3+8) Sqrt[t^3-1])`) return `NULL`.  A differentiate-back guard rejects the
     rare cases where the eigenspace zero-test misfires on deeply nested radical
     roots; the whole attempt runs under a CPU-time budget so a cyclotomic-root
     `R` with an unlucky cofactor (where algebraic-number `Together`/`Cancel`
     blows up) declines rather than hanging the cascade.  Uses
     `Solve[..., Cubics -> True, Quartics -> True]` (the Ferrari quartic solver,
     added 2026-06-29).  The Boolean global `Integrate\`GoursatDebug` (default
     `False`; added 2026-06-30) traces the descent to stderr when set `True`:
     whether the integrand matches the `F(x) R(x)^(-p)` form (and the recognised
     `F`, `R`, `p`), which involution / eigenspace criterion is tested and
     whether it holds (`V4` trivial projection, the order-3/order-4
     ω-eigencomponents, the period-3 trivial projection per fixed point), and the
     differentiate-back verdict.  The flag is latched once at the outermost call,
     so recursive genus-0 reductions share it and indent by depth.
     Before returning (square-root case), the antiderivative is normalised so its
     one radical is the single generator `Sqrt[R]`: any term left carrying `Sqrt[R]`
     *split* across a numerator/denominator as a product of factor roots
     (`Sqrt[x] Sqrt[(1-x)(1-k^2 x)] ...`) has its `x`-dependent half-power factors
     recombined into one radicand and reduced over `R` (to a rational when it is a
     perfect square, else `rational*Sqrt[R]`), removing the spurious branch point and
     keeping a downstream `Simplify` to a single radical generator.  Constant radicals
     are left intact, and the rewrite is kept only if the differentiate-back guard
     still passes.
     A graded battery of worked examples (every exponent `p`, both numerator and
     denominator radicals, and every involution equation, with the negative
     controls that decline) is collected in
     [`GOURSAT_EXERCISES.md`](../../../GOURSAT_EXERCISES.md) and mirrored as the
     `test_graded` ladder in `tests/test_integrate_goursat.c`.
     - **Fresnel stage** (Automatic cascade, immediately after Goursat): a
       Gaussian-phase trig integrand `K Sin[a x^2 + b x + c]` or `K Cos[...]`
       (`a != 0`, `K` free of `x`) closes to `FresnelS`/`FresnelC` by completing
       the square — the trigonometric sibling of the `K E^(a x^2+b x+c) -> Erf`
       recognizer, deterministic and diff-back verified
       (`src/calculus/integrate_fresnel.c`): `Integrate[Sin[x^2],x] =
       Sqrt[Pi/2] FresnelS[Sqrt[2/Pi] x]`, `Integrate[Cos[x^2],x] =
       Sqrt[Pi/2] FresnelC[Sqrt[2/Pi] x]`, and completing-square cases like
       `Sin[x^2+x+1]`.
  8. `Integrate\`Weierstrass[f, x]` — rational functions of the trig kernels
     `Sin/Cos/Tan/Cot/Sec/Csc[x]` (or hyperbolic `Sinh/Cosh/.../Csch[x]`) with a
     kernel in a denominator; continuous `Tan[x/2]` / `Tanh[x/2]` substitution
     (Jeffrey & Rich 1994).  Runs ahead of `DerivativeDivides`: it is
     domain-specific, deterministic, correct by construction, and yields a real,
     continuous antiderivative rather than a complex-logarithm form.
  9. `Integrate\`DerivativeDivides[f, x]` — substitution `u(x)`; in the
     cascade the quiet, branch-correct **direct quotient** strategy only.
  10. `Integrate\`RischNorman[f, x]` — Bronstein pmint (parallel Risch), all
     integrands.
  11. `Integrate\`RischTranscendental[f, x]` — the **recursive** transcendental
     Risch algorithm; runs after RischNorman and only adds
     closed forms the earlier stages missed.  Correct by construction (no
     differentiation check).  Handles logarithmic polynomials and the
     special-function cases below (Erf, ExpIntegralEi, LogIntegral, PolyLog).
     Its single-extension cases (fractional Rothstein-Trager log-part, Hermite,
     hyperexponential) are decision procedures over the field `C(x)(t)` for the
     kernel `t = Log[u]` or `E^u`, so each first checks — via `rt_is_ratl_in_xt`
     — that the kernelized integrand is genuinely rational in x and t; an
     integrand with a transcendental coefficient of x (`Sin[x]/Log[x]`,
     `Sin[x]/(1+E^x)`, `Gamma[x]/Log[x]`) lies outside that field and is
     **declined**, never mis-certified (this previously produced a wrong `0`).
     Its Risch differential equation is solved by Bronstein's rational one-step
     (SPDE) reduction (polynomial-gcd time, no undetermined-coefficient
     blow-up), closing high-degree `R(x) e^x` forms and — via the
     exponential special-denominator Laurent ansatz and an exact
     tower-variable verification — nested exp/log towers such as
     `∫((e^x−x²+2x)/(x²(x+e^x)²)) e^((x²−1)/x + 1/(x+e^x)) dx = e^(−x+1/(e^x+x)+(x²−1)/x)`
     and `∫(1/(x log(1+e^x)) − …) dx = log(x)/log(1+e^x)`.  Also closes
     polynomial × `Log` with an irreducible-quadratic argument (a new `ArcTan`
     in the answer, e.g. `∫(x⁵−1) log(x²−x+1) dx`), Q-linearly **dependent**
     logarithms (`∫(log(x/(1+x))+log(1+x))/x = log(x)²/2`), mixed **tangent + log**
     towers (`∫12x⁴+5x⁴ log(x¹²cos x)−x⁵ tan x = x⁵ log(x¹²cos x)`), and deeper
     nested exponentials — including a constant-plus-rational inner exponent
     (`∫… = x e^(x e^(1+1/x))`) and a nested exp inside a circular kernel
     (`∫e^(e^x)(1+x e^x) cos(x e^(e^x)) = sin(x e^(e^x))`).  Arithmetic warnings
     from transient internal singular expressions are muted (as in Mathematica).
  12. `Integrate\`CRCTable[f, x]` — CRC integral table lookup (lazy-loaded
     from `src/internal/CRCMathTablesIntegrals.m` on first call).  Integer-power
     reduction rules cover all six circular and six hyperbolic functions
     (`Sin/Cos/Tan/Cot/Sec/Csc` and `Sinh/…/Csch`, argument `a x`), their
     products (`Sinh^m Cosh^n`), canonical-head reciprocal products
     (`Csc^m Sec^n`, `Csch^m Sech^n`), mixed tangent/secant powers
     (`Tan^m Sec^n`, `Cot^m Csc^n` and hyperbolic analogues — the canonical form
     of `Sin^m/Cos^n` quotients), `E^(a x)` times circular or hyperbolic
     powers (the hyperbolic case guards the `a = n b` resonance), and
     polynomial × hyperbolic (`x^n Sinh/Cosh`, `x Sinh^m`, `x/Sinh^n = x Csch^n`).
  If every stage gives up the call bubbles back unevaluated.
- `Method -> "<name>"` option (3rd argument) bypasses the cascade and
  dispatches strictly to a single subroutine, with no fallback:
  - `"Automatic"` — default cascade above.
  - `"BronsteinRational"` — `Integrate\`BronsteinRational[f, x]`.
  - `"DerivativeDivides"` — `Integrate\`DerivativeDivides[f, x]` (direct **and**
    the more thorough Eliminate/Solve branch search).  The list form
    `Method -> {"DerivativeDivides", "Substitution" -> u}` **pins** the kernel
    `u(x)`: instead of collecting and trialing every `x`-dependent subexpression,
    only that one substitution is attempted (still both strategies, still strict
    — no fallback to other kernels if `u` does not close the integral). E.g.
    `Integrate[Sqrt[x]/(1 + Sqrt[x]), x, Method -> {"DerivativeDivides", "Substitution" -> Sqrt[x]}]`.
  - `"LinearRadicals"` — `Integrate\`LinearRadicals[f, x]`.
  - `"QuadraticRadicals"` — `Integrate\`QuadraticRadicals[f, x]`.
  - `"LinearRatioRadicals"` — `Integrate\`LinearRatioRadicals[f, x]`.
  - `"ChebychevAlgebraic"` — `Integrate\`ChebychevAlgebraic[f, x]` (Chebychev
    binomial differential `x^p (a x^r + b)^q`).
  - `"GoursatAlgebraic"` — `Integrate\`GoursatAlgebraic[f, x]` (pseudo-elliptic
    `F/R^p`, `p` in `{1/2, 1/3, 2/3, 1/4, 3/4}`, via Mobius eigendescent).
  - `"Weierstrass"` — `Integrate\`Weierstrass[f, x]` (no denominator gate: applies
    to any rational function of the trig/hyperbolic kernels of `x`, including
    polynomial trig).
  - `"RischNorman"` — `Integrate\`RischNorman[f, x]` (parallel Risch / pmint).
  - `"RischTranscendental"` — `Integrate\`RischTranscendental[f, x]`, the recursive
    transcendental Risch algorithm (`src/calculus/integrate_risch_transcendental.c`).
    A decision procedure over a differential transcendental tower, distinct
    from the parallel-Risch heuristic `"RischNorman"`.  Every case is correct
    by construction — it fires only behind an exact structural certificate, so
    the result is not checked by differentiation.  Cases:
      - rational: delegated to `Integrate\`BronsteinRational`;
      - logarithmic polynomial `P(x, Log[u])`: the recursive primitive-
        polynomial coefficient matching, with a limited-integration oracle that
        folds a would-be new logarithm back into the tower (e.g.
        `Integrate[Log[2 x + 3], x]`, `Integrate[Log[x]/x, x]`);
      - exponential (Laurent) polynomial `sum_i p_i(x) E^(i u)`, `u` polynomial
        in `x`, `i` positive or negative: the powers of `E^u` decouple and each
        `i != 0` term solves the Risch differential equation
        `q_i' + i u' q_i = p_i` by a polynomial ansatz
        (`Integrate[x E^x, x] = (x-1) E^x`, `Integrate[x E^(x^2), x] = E^(x^2)/2`,
        `Integrate[(E^x+E^(-x))/2, x] = Sinh[x]`);
      - Hermite reduction for a repeated pole of `theta = Log[u]` or
        `theta = E^u` (the latter when `D` is coprime to `theta`):
        `Q = H(theta)/Hden(theta) + sum_j c_j Log(g_j)` with
        `Hden = gcd(D, dD/dtheta)`, solved by `SolveAlways` over `theta` and `x`
        (`Integrate[1/(x (1+Log[x])^2), x] = -1/(1+Log[x])`,
        `Integrate[E^x/(1+E^x)^2, x] = -1/(1+E^x)`);
      - a coupled hyperexponential case (a unified ansatz
        `Q = sum_i w_i(x) E^(i u) + H(E^u)/Hden(E^u) + sum_j c_j Log(g_j)` solved
        by `SolveAlways` over `theta` and `x`) that closes mixed
        polynomial-plus-log exponentials such as
        `Integrate[1/(1 + E^x), x] = x - Log[1 + E^x]`, and — with the Hermite
        term `H/Hden` fused in (the `theta`-coprime denominator split into its
        repeated part `Hden = gcd(Dtil, dDtil/dtheta)` and squarefree radical) —
        the repeated / `theta = 0` exponential poles
        `Integrate[1/(1 + E^x)^2, x] = x + 1/(1 + E^x) - Log[1 + E^x]`,
        `Integrate[1/(E^x (1 + E^x)^2), x]`, `Integrate[1/(1 + E^x)^3, x]`;
      - a multi-kernel **sum-of-exponentials** case: an integrand that
        exponentializes to a sum `sum_k p_k(x) E^(W_k)` of NON-commensurate
        exponentials `E^(W_k)` (e.g. the `(1 ± I) x` pair from `E^x Sin[x]`)
        decouples — each term integrates by its own Risch DE
        `q_k' + W_k' q_k = p_k` — closing `Integrate[E^x Sin[x], x]`,
        `Integrate[x E^x Sin[x], x]`, `Integrate[E^(2x) Cos[3x], x]`, ... (via
        the complex exponentials the answer is left in an I-laden `Cosh`/`Sinh`
        form, a `Simplify` opportunity; the diff-back is exactly `0`);
      - nested **logarithmic** and **exponential tower** cases: a rational
        function of a chain of nested logarithms (`Log[x]`, `Log[Log[x]]`, ...) or
        nested exponentials (`E^x`, `E^(E^x)`, ...) is integrated over the tower
        derivation `D = d/dx + sum_i Dt_i d/dt_i` by one unified `SolveAlways`
        ansatz over all tower variables, closing
        `Integrate[1/(x Log[x] Log[Log[x]]), x] = Log[Log[Log[x]]]`,
        `Integrate[Log[Log[x]]/(x Log[x]), x] = Log[Log[x]]^2/2`,
        `Integrate[E^x E^(E^x), x] = E^(E^x)`, ...  The tower cases are
        bounded-ansatz searches and so are diff-back verified (a non-elementary
        integrand declines rather than returning a spurious form); a whole-tower
        rationality gate and a single-kernel nesting gate keep the other cases
        from ever certifying a wrong nested answer;
      - the **genuine one-extension-at-a-time recursion** (Bronstein ch. 5),
        which the flat tower ansatz above cannot express: it builds
        the ordered differential tower (structure-theorem triangularity check) and
        peels one kernel at a time, integrating the polynomial/Laurent part in the
        top kernel *coefficient by coefficient* — each coefficient integral is an
        integration in the LOWER field (the recursion), bottoming out in `C(x)`.
        This closes **mixed exp/log towers** and **rational lower-field
        coefficients** the single-kind flat cases decline, e.g.
        `Integrate[E^x/x + E^x Log[x], x] = E^x Log[x]`,
        `Integrate[Log[1+E^x] + x E^x/(1+E^x), x] = x Log[1+E^x]`,
        `Integrate[1/(x^2 Log[x]) - Log[Log[x]]/x^2, x] = Log[Log[x]]/x`; diff-back
        verified, so non-elementary mixed integrands (`E^x Log[x]`) decline.  A
        **proper rational part** at a logarithmic top level (a genuine `t_n`-pole)
        is integrated by tower **Hermite reduction + Rothstein–Trager**
        (`H(t)/Hden + sum_j c_j Log(g_j)`, `Hden = gcd(den, d den/dt_n)`, constant
        residues), closing `Integrate[1/(x Log[x] (1+Log[Log[x]])^2), x] =
        -1/(1+Log[Log[x]])` and `Integrate[1/(x(1+Log[x])) + E^x, x] = E^x +
        Log[1+Log[x]]` (a non-constant residue declines).  For an **exponential**
        top level the Laurent and log parts couple, so a proper `t_n`-pole is closed
        by a unified coupled-hyperexponential ansatz over the tower derivation, e.g.
        `Integrate[(2 Log[x]/x) E^(Log[x]^2)/(1+E^(Log[x]^2)), x] =
        Log[1+E^(Log[x]^2)]`.  The exponential Laurent step solves a **field Risch
        differential equation** `q_i' + i w_n' q_i = p_i` in the lower field for each
        power; a `q_i` that is rational there — with an arbitrary denominator, by the
        RDE denominator theorem `q_i = h/Denominator[p_i]` with `h` a bounded
        polynomial — is found, closing
        `Integrate[(2/x - 1/(x Log[x]^2)) E^(Log[x]^2), x] = E^(Log[x]^2)/Log[x]`
        (`q = 1/Log[x]`) and
        `Integrate[(2 Log[x]/(x(1+Log[x])) - 1/(x(1+Log[x])^2)) E^(Log[x]^2), x] =
        E^(Log[x]^2)/(1+Log[x])` (`q = 1/(1+Log[x])`, non-monomial denominator).
        A structural pre-pass re-splits an **evaluator-merged exponential monomial**
        `E^(a+b) -> E^a E^b` (the evaluator folds `E^x E^(E^x)` into `E^(x+E^x)`,
        whose exponent carries the foreign kernel `E^x` and breaks tower
        independence), restoring the independent basis `{E^x, E^(E^x)}` and closing
        `Integrate[E^x E^(E^x)/(1+E^(E^x)), x] = Log[1+E^(E^x)]` (the non-elementary
        `E^(E^x)/(1+E^(E^x))` still declines);
      - a trig/hyperbolic front-end (`TrigToExp` -> exponential machinery ->
        `ExpToTrig`) that closes `Sin`, `Cos`, `Sinh`, `Cosh`, `Sin[x]^2`,
        `Sin[x] Cos[x]`, `Tan`, `Tanh`, ...; through the complex substitution
        `Tan`/`Tanh` come out in a correct but I-laden form (e.g.
        `I x - Log[1 + E^(2 I x)] = -Log[Cos[x]]`) that no current simplifier
        reduces to real closed form (a `Simplify` improvement opportunity).
        A structural pre-pass `rt_powers_to_exp` re-exposes a **transcendental
        general power** `b^e` (base `b != e`, non-rational exponent) as a raw
        base-e kernel `E^(e Log b)` — the evaluator collapses `E^(c Log b) -> b^c`,
        so a trig/inverse-trig *of a logarithm* hides its exponential kernels as
        general powers (`Cos[x Log x] -> (x^(I x) + x^(-I x))/2`, kernels
        `E^(±I x Log x)`) and evades every `Exp`/`Power[E, .]`-keyed recognizer.
        Rewriting without going through the evaluator keeps the raw `Power[E, .]`
        spelling the tower collectors/substitution already match, with `Log b` a
        genuine logarithmic sub-kernel; the exponent is split so an
        evaluator-merged polynomial coefficient (`x^2 x^(I x) -> x^(2+I x)`) peels
        back to an algebraic power times the shared primitive kernel `x^(I x)`.
        This closes the elementary mixed log/exp-power towers
        `Integrate[Cos[Log x], x] = (x/2)(Cos[Log x] + Sin[Log x])`,
        `Integrate[Sin[Log x], x] = (x/2)(Sin[Log x] - Cos[Log x])`,
        `Integrate[x^I, x] = (1/2 - I/2) x^(1+I)`,
        `Integrate[x^x (1 + Log x), x] = x^x`, and
        `Integrate[3 x^2 Cos[x Log x] - x^3 (1 + Log x) Sin[x Log x], x] =
        x^3 Cos[x Log x]`; the non-elementary siblings (`Cos[x Log x]`, `x^x`)
        still decline;
      - `K E^(a x^2 + b x + c)` (`a != 0`) → `Erf`/`Erfi` (Gaussian fast-path);
      - `g E^f` with `f = p/q`, `q = s^2` a perfect square → up to two `Erfi`
        error-function terms (Cherry 1989 §3, completing-square), emitted in the
        same solve as the `ei` terms: `Integrate[(1/x+1/x^2)E^(1/x^2),x] =
        -ExpIntegralEi[1/x^2]/2 - Sqrt[Pi] Erfi[1/x]/2`, and (symbolic parameter)
        `Integrate[E^((x^4+a)/x^2),x] = (Sqrt[Pi]/4)(Erfi[(x^2+Sqrt[a])/x]/E^(2Sqrt[a])
        + E^(2Sqrt[a]) Erfi[(x^2-Sqrt[a])/x])`. A non-square `q` gives no erf term;
      - `g E^f` with `g, f` rational in `x` → `ExpIntegralEi` (Cherry 1989, base
        field) — the general rational exponential integral. The `ei` arguments
        `f + alpha_i` are generated by the resultant `Res_x(g1, p + alpha q)`
        (`f = p/q`, `g1` = the `q`-coprime part of `den(g)`) plus the `q`-side
        term `f` itself, and `y` + the constants solved as one linear system over
        the rational matching identity `g = y' + f' y + sum_i c_i f'/(f+alpha_i)`;
        emits `y E^f + sum_i c_i E^(-alpha_i) ExpIntegralEi[f+alpha_i]`. Closes
        e.g. `Integrate[E^(1/x),x] = x E^(1/x) - ExpIntegralEi[1/x]`,
        `Integrate[E^x/(x+1)^2,x] = ExpIntegralEi[x+1]/E - E^x/(x+1)`,
        `Integrate[(x^2+3)E^x/(x^2+3x+2),x] = E^x - 7 ExpIntegralEi[x+2]/E^2 +
        4 ExpIntegralEi[x+1]/E`, and nonlinear exponents
        `Integrate[E^(x^2)/x,x] = ExpIntegralEi[x^2]/2`. A narrow `(M E^(a x + b))/
        (c x + d)` fast-path handles the linear case in O(1). Higher-multiplicity
        poles close (the elementary part `y` carries the pole to multiplicity
        `m-1`), e.g. `Integrate[E^x/(x+1)^3,x] = ExpIntegralEi[x+1]/(2E) -
        E^x(x+2)/(2(x+1)^2)`. **Real algebraic** `ei`-argument constants are
        admitted (Cherry's algebraically-closed constant field): the resultant
        roots may be irrational reals and the coefficients are solved over the
        extension they generate — `Integrate[E^x/(x^2-2),x] =
        (E^Sqrt[2] ExpIntegralEi[x-Sqrt[2]] - E^(-Sqrt[2]) ExpIntegralEi[x+Sqrt[2]])
        /(2 Sqrt[2])` (Cherry p.894), also `E^x/(x^2-x-1)` (golden ratio) etc.
        A single irreducible-quadratic denominator with **complex** conjugate
        roots closes with a conjugate pair of `ExpIntegralEi` over the number field
        the roots generate — `Q(i)` **and** `Q(i sqrt d)`:
        `Integrate[E^x/(x^2+1),x] = (I/2)(ExpIntegralEi[x+I]/E^I - E^I
        ExpIntegralEi[x-I])`, `E^x/(x^2+x+1)` (over `Q(i sqrt3)`), `E^x/(x^2+3)`,
        `E^x/(x^2+2x+5)`, a numerator `(2x+1)E^x/(x^2+x+3)` (over `Q(i sqrt11)`),
        and the polynomial-part case **d12** `Integrate[(x^2+1)E^x/(x^2+x+1),x] =
        E^x + (complex-conjugate ei pair)`. The lone conjugate pair closes over ANY
        `Q(i sqrt d)` — a shifted center (`E^x/(x^2+2x+3)`, `E^x/(x^2+2x+7)`,
        `(3x+1)E^x/(x^2+2x+3)`) that the direct number-field `Solve` cannot crack is
        solved over Q by the `{1, chs}`-basis fallback. Complex poles mixed with a
        `P2`/reciprocal term (`E^(1/x)/(x^2+1)`) and degree-`≥3` constants
        (`E^x/(x^3-2)`) are deferred and decline cleanly. A constant exponent offset
        `E^(c + h(x))` (e.g. `E^(1/x+2) = E^2 E^(1/x)`) is factored out before matching.
        A **multi-term**
        exponential `Sum_i p_i E^(i w)` (rational in a single kernel `E^w`, several
        commensurate Laurent terms) that the single-shape matcher cannot peel is
        integrated term-by-term (Cherry Thm 5.4 case b) —
        `Integrate[(E^x+E^(2x))/(x-1),x] = E ExpIntegralEi[x-1] +
        E^2 ExpIntegralEi[2x-2]`, `E^x(E^x+1)/((x-1)(x-2))`;
      - **error-function integration of transcendental Liouvillian integrands**
        (Knowles 1992/93, the layer above Cherry): an integrand that may itself
        contain `Erf`/`Erfi` is closed in terms of elementary functions plus
        `Erf`/`Erfi` over the `RT_PRIM` Liouvillian-primitive tower. For each
        exponential monomial `e^w` the perfect-square gate (`-w = u^2 -> erf(u)`,
        `+w = u^2 -> erfi(u)`) generates the error-function arguments, and
        `v + Sum_i k_i erf(u_i)` is solved as one linear system over the constants
        (diff-back verified): `Integrate[E^(-x^2-Erf[x]^2),x] = (Pi/4) Erf[Erf[x]]`,
        `Integrate[Erf[x] E^(-x^2-Erf[x]^2),x] = -(Sqrt[Pi]/4) E^(-Erf[x]^2)`,
        `Integrate[2 E^(-x^2) Erf[x] - 3 E^(-1/x^2)/x^2,x] = (Sqrt[Pi]/2) Erf[x]^2 +
        (3 Sqrt[Pi]/2) Erf[1/x]`; `Integrate[x E^(-x^2-Erf[x]^2),x]` correctly
        declines (no erf-elementary antiderivative). The **radical (quasiquadratic)**
        case is also handled (Knowles Part I §6, Ex 8.1): when a top exponential
        hides an algebraic factor `E^(l/2 Log[g]) = g^(l/2)`, it is pulled out first,
        surfacing a half-integer power of a log tower variable; the erf argument is
        then a radical `1/Sqrt[Log g]`, solved in `s = Sqrt[Log g]` — e.g.
        `Integrate[E^(1/2 Log[Log[x]] - 1/Log[x])/(x Log[x]^2),x] =
        -Sqrt[Pi] Erf[1/Sqrt[Log[x]]]` (and its `Erfi` dual). Both **spellings**
        of that integrand reach the engine: the already-reduced
        `Integrate[E^(-1/Log[x])/(x Log[x]^(3/2)),x]` carries the half-integer
        power openly, and the `RischTranscendental` scope gate — which refuses
        every algebraic function of `x`, rightly, since the recursive Risch
        algorithm is not a decision procedure over an algebraic extension —
        admits it as a narrow exception: when *every* algebraic site is a
        fractional power of a transcendental **kernel**, each is rewritten
        `g^(p/q) = g^n E^(r Log g)` (the inverse of the collapse above) and the
        answer is returned only behind a diff-back check, because that rewritten
        tower hides the algebraic relation `(E^(r Log g))^q = g^p`. A genuine
        radical of `x` (`Sqrt[x]`, `Sqrt[1+x^2]`, `Cos[Sqrt[x]]`) is still
        refused and left to the algebraic routes, and the non-elementary
        decision half never runs on the admitted path. `x`-rational
        elementary-part coefficients remain a later increment;
      - `c w^(p-1) w'/Log[w]` → `c LogIntegral[w^p]` (single-li fast path) —
        subsumes `K/Log[x] → K LogIntegral[x]`, a scaled/affine argument
        (`1/Log[2x] → LogIntegral[2x]/2`) and a monomial numerator (`x/Log[x] →
        LogIntegral[x^2]`);
      - general single-log tower `gamma` over `C(x, Log[w])` (w a polynomial) →
        `v(x,Log[w]) + Sum_k d_k LogIntegral[w^k]` (Cherry 1986, MULTI-li via the
        degree-1 Sigma-decomposition over the generator `w`, matched in the tower):
        `Integrate[x^2/Log[x+1],x] = LogIntegral[(1+x)^3] - 2 LogIntegral[(1+x)^2]
        + LogIntegral[1+x]`, `Integrate[x/Log[x]^2,x] = 2 LogIntegral[x^2] -
        x^2/Log[x]`, `Integrate[x^3/Log[x^2-1],x] = LogIntegral[x^2-1]/2 +
        LogIntegral[(x^2-1)^2]/2`. Includes the **transcendental-constant
        rescaling** (a constant root `rho` of the `Log`-denominator →
        `e^rho LogIntegral[e^(-rho) w]`): `Integrate[1/(Log[x]+3),x] =
        LogIntegral[E^3 x]/E^3`, `Integrate[(Log[x]^2+3)/(Log[x]^2+3Log[x]+2),x] =
        x + 4 LogIntegral[E x]/E - 7 LogIntegral[E^2 x]/E^2`. Multi-log towers
        (reducible `w` needing a product decomposition) decline. The tower
        Sigma-decomposition is **Laurent**: a negative-degree (Laurent-polynomial)
        numerator over `Log[w]` admits negative li powers `w^k`, emitted as
        `ExpIntegralEi[k Log[w]]` (`= LogIntegral[w^k]`, the Wolfram form, which
        also avoids `LogIntegral` of an argument `<1`), and higher `Log` poles
        carry an elementary Laurent part `Sum_i b_i(x) Log[w]^i` with `b_i`
        themselves Laurent in `x`: `Integrate[1/(x^4 Log[x]),x] =
        ExpIntegralEi[-3 Log[x]]`, `Integrate[1/(x^4 Log[x]^2),x] =
        -3 ExpIntegralEi[-3 Log[x]] - 1/(x^3 Log[x])`,
        `Integrate[(x-1)^2/(x^4 Log[x]^3),x]` (three `ExpIntegralEi` terms plus a
        Laurent elementary part);
      - `c g'(x)/Log[g]` → `c LogIntegral[g]` for a **mixed exp/log kernel** `g`
        that is not a polynomial in `x` (so `Log[g]` appears as a SUM, not a
        `Log[...]` node, and the polynomial-`w` engines above cannot see it). The
        kernel is discovered from the denominator: for `f = N/D` (`Together`), the
        only elementary `g` with `f ∝ g'/Log[g]` has `Log[g] = D`, i.e. `g =
        Exp[D]`; the certificate `Together[f D / g']` must be a nonzero constant,
        and a `PowerExpand` diff-back reconfirms the branch identity. Closes
        `Integrate[(E^x (1+x))/(x + Log[x]), x] = LogIntegral[x E^x]` (here
        `Log[x E^x] = x + Log[x]`). Gated to a `D` containing a logarithm — a
        `Log`-free `D` is either `Ei`-reducible (`LogIntegral[Exp[u]] =
        ExpIntegralEi[u]`, owned by the `Ei` recognizer) or spurious, so it
        declines fast without normalizing over an `Exp[...]` generator;
      - `K Log[1 + p x]/x` → `PolyLog[2, -p x]` (dilogarithm fast path); and the
        general `R(x) Log[w]` → `Log-Log + PolyLog[2]` form (Cherry degree-2
        Sigma-decomposition, matched in the log tower with dilog arguments the
        interpolants `(x-r_i)/(r_j-r_i)` between the rational roots of the linear
        factors): `Integrate[Log[x]/(1+x),x] = Log[x] Log[1+x] + PolyLog[2,-x]`,
        `Integrate[Log[x]/(1-x),x] = PolyLog[2,1-x]`, `Log[x]/(x^2-1)`,
        `Log[2x+1]/(x+1)`. **Transcendental-constant root spacing** (spacing `!= 1`)
        now closes too — a dilog whose derivative leaves a `Log` of a positive
        constant contributes a real `Log-Log` term:
        `Integrate[Log[2+x]/x,x] = Log[2] Log[x] - PolyLog[2,-x/2]`,
        `Log[2x+3]/(x-1)`. Non-monic linear kernels (`Log[3+2x]/x`) and products
        of logs still decline;
      - fractional (Rothstein–Trager) log-part: a proper rational function of
        `theta` with squarefree denominator `prod g_i` gives `sum_i c_i Log(g_i)`,
        the constant residues `c_i` solved from `num = sum_i c_i D(g_i)(d/g_i)`
        via `SolveAlways` over `theta` and `x` (`Integrate[1/(x(1+Log[x])),x] =
        Log[1+Log[x]]`, `Integrate[E^x/(1+E^x),x] = Log[1+E^x]`); the residues are
        verified constant (an `x`-dependent `SolveAlways` pseudo-solution is
        rejected so `1/Log[c x+d]` is never mis-closed with a polynomial-coefficient
        logarithm);
      - **pure resultant Lazard–Rioboo–Trager** log-part (when the `SolveAlways`
        path above declines): an irreducible-over-Q denominator factor in `theta`
        with **algebraic (non-rational) residues** is closed by the exact resultant
        `Res_t(a - z D(d), d)` (residues) + Rioboo `LogToReal` (real `Log + ArcTan`
        form), delegated to the internal `Integrate`TranscendentalLogPart` (reuses
        the rational-LRT machinery with the monomial derivation), closing
        `Integrate[1/(x(Log[x]^2+1)),x] = ArcTan[Log[x]]`,
        `Integrate[E^x/(E^(2x)+1),x] = ArcTan[E^x]`, mixed
        `Integrate[(1+Log[x])/(x(Log[x]^2+1)),x] = ArcTan[Log[x]] + Log[1+Log[x]^2]/2`
        and higher-degree denominators; diff-back verified, so non-elementary
        siblings (`1/(Log[x]^2+1)`, `1/(x^2(Log[x]^2+1))`) decline; this same
        resultant LRT is **lifted into the tower recursion** at a logarithmic top,
        so a nested-log proper part with algebraic residues closes too —
        `Integrate[1/(x Log[x] (Log[Log[x]]^2+1)),x] = ArcTan[Log[Log[x]]]` — with
        the residues gated free of *every* lower-field variable (`{x, Log[x], …}`),
        and `1/(Log[x](Log[Log[x]]^2+1))` (x-dependent residues) declining;
    Rational-argument exponents such as `E^(1/x)` are handled (`Integrate[-E^(1/x)/
    x^2, x] = E^(1/x)`) via the `q = h/Denominator[p]` RDE ansatz.  **Multiplicatively
    commensurate merged kernels** are reduced in the tower builder — a collected
    exponential whose exponent is an integer multiple of a class primitive
    (`E^(2 E^x) = (E^(E^x))^2`) becomes a power of the primitive's tower variable
    instead of a spurious extra extension — so towers with `E^(k u)` close and the
    exp-top algebraic-residue LRT is unblocked
    (`Integrate[E^x E^(E^x)/(1+E^(2 E^x)), x] = ArcTan[E^(E^x)]`,
    `Integrate[E^x E^(2 E^x)/(1+E^(E^x)), x] = E^(E^x) - Log[1+E^(E^x)]`).  The class
    primitive is **synthesized** (`p = w_0/lcm(ratio denominators)`) rather than
    required to be a member, so *non-integer* commensurate exponents close too:
    `Integrate[1/(E^(x/2)+E^(x/3)), x]` (primitive `E^(x/6)`) and the nested
    `Integrate[D[E^(E^x/2)/(1+E^(2 E^x/3)), x], x]` (primitive `E^(E^x/6)`); only a
    class whose members are not all rational multiples of one another declines.
    The **RDE-solver degree bounds** are
    Bronstein's exact `RdeBoundDegree` (leading-degree balance,
    `deg_v(q) = deg_v(p) − deg_v(f)` where the exponential dominates), with **no
    arbitrary cap**, so an exponential-Laurent coefficient of any degree closes
    (`Integrate[(6 Log[x]^5 + 2 Log[x]^7)/x E^(Log[x]^2), x] = Log[x]^6 E^(Log[x]^2)`,
    and deg-20 …).  The **flat-tower and proper-part Hermite ansätze are likewise
    cap-free**: exact top-kernel log/exp Laurent bounds and derived inner-exp windows,
    so all top degrees close (`Integrate[Log[Log[x]]^5/(x Log[x]), x] =
    Log[Log[x]]^6/6`, `Integrate[E^x E^(6 E^x)/(1+E^(E^x)), x]`).  The
    **leading-coefficient cancellation / resonance** sub-case of `RdeBoundDegree`
    completes the SPDE degree machinery — the bound is widened monotonically to
    `max(naive, m_res)` at the Bronstein resonance integer `m_res` (detected live for the
    exponential top), correct by the same certification-and-diff-back gate.  Only
    algebraic extensions (`Sqrt`, `RootSum`) remain unimplemented, so integrands needing
    them return unevaluated.
  - `"CRCTable"` — `Integrate\`CRCTable[f, x]`.
  - `"Undefined"` — `Integrate\`Undefined[f, x]`.
  - `"Symmetry"` — origin-symmetry reduction for an interval `[-c, c]`
    (`Integrate\`Symmetry[f, {x, -c, c}]`): an odd integrand integrates to `0`,
    an even one to `2 Integrate[f, {x, 0, c}]`. The parity is proved by
    `Simplify`, and a value is claimed only when the half integral converges, so
    a divergent principal value is never reported as `0`. Under Automatic it runs
    after residue and before Newton-Leibniz.
  - `"Beta"` — Euler-Beta reduction on `[0,1]`
    (`Integrate\`Beta[f, {x, 0, 1}]`): `x^(k-1) (1-x)^(l-1) → Beta[k, l]`, with
    `Log[x]^i Log[1-x]^j` weights giving the mixed parameter derivative of
    `Beta`. Gated on `Re[k] > 0 && Re[l] > 0`.
  - `"TrigPower"` — `Sin[x]^m Cos[x]^n` over a canonical trig interval
    (`Integrate\`TrigPower[f, {x, 0, c}]`): over `[0, Pi/2]` it is
    `Beta[(m+1)/2, (n+1)/2]/2`; over `[0, Pi]`/`[0, 2Pi]` the standard parity
    multipliers apply (an odd power integrates to `0`).
  - `"NewtonLeibniz"` — the real-axis definite-integral mechanism (implicit for
    the `{x, a, b}` form); see **Definite integration** below.
  - `"LineIntegral"` — the complex contour mechanism (implicit for the
    `{x, z0, …, zn}` form); see **Complex line integration** below.
  - `"DiffUnderInt"` (alias `"DifferentiationUnderIntegral"`) — definite
    integration by differentiation under the integral sign (Feynman's trick);
    `Integrate\`DiffUnderInt[f, {x, a, b}]`. Tried last in the definite cascade
    (after residue and Newton-Leibniz). See **Differentiation under the integral
    sign** below.
  - `"SinPowerMonomial"` — `Sin[r x]^k / x^m` on `[0, Infinity)` (the ssp
    family); `Integrate\`SinPowerMonomial[f, {x, 0, Infinity}]`.
  - `"OscillatoryPower"` — Fresnel-type `Cos[b x^n]` / `Sin[b x^n]` on
    `[0, Infinity)`; `Integrate\`OscillatoryPower[f, {x, 0, Infinity}]`.
  - `"RationalLog"` — `R(x) Log[x]^n` on `[0, Infinity)` for a proper rational
    `R` with negative-real-axis poles; `Integrate\`RationalLog[f, {x, 0, Infinity}]`.
  - `"RamanujanMasterTheorem"` (alias `"Mellin"`) — half-line `∫₀^∞ x^{s-1} f(x) dx`
    by the Mellin-transform / Ramanujan Master Theorem method;
    `Integrate\`RamanujanMasterTheorem[f, {x, 0, Infinity}]`. Under Automatic it
    runs after Newton-Leibniz and before DiffUnderInt. See **Mellin / Ramanujan
    Master Theorem** below.
  The definite mechanisms name themselves only: the actual mechanism is
  chosen from the spec type, so on a definite integral any *other* method name
  is passed through to the inner indefinite integration that produces the
  antiderivative, and either definite-mechanism name on the indefinite
  `Integrate[f, x, …]` form is a no-op (stays unevaluated).
  Unknown method names emit `Integrate::method` and bubble back.
- Universal correctness predicate: `Cancel[Together[D[Integrate[f,x],x] - f]] === 0`.

**Attributes:** `Protected`.

## References

**See also:** [PolynomialQuotientRemainder](../../calculus/PolynomialQuotientRemainder/), [Apart](../../algebra/Apart/), [Log](../../elementary-functions/Log/), [ToRadicals](../../solutions-of-equations/ToRadicals/), [Root](../../solutions-of-equations/Root/), [RootSum](../../solutions-of-equations/RootSum/), [Plus](../../arithmetic/Plus/), [Sqrt](../../arithmetic/Sqrt/)

- Bronstein, "Symbolic Integration I: Transcendental Functions", 2nd ed. (Springer, 2005).
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (Kluwer, 1992), ch. 11–12.
- M. Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed. (Springer, 2005).
- M. Bronstein, "poor man's integrator" (pmint), 2004 — parallel Risch / Risch–Norman heuristic.
- B. M. Trager, "Algebraic Factoring and Rational Function Integration", ACM SYMSAC 1976.
- K. O. Geddes, S. R. Czapor, G. Labahn, *Algorithms for Computer Algebra* (Kluwer, 1992).
- K. Roach, "Symbolic Integration", 1992 (undefined-function integrands, §1.7).
- CRC Standard Mathematical Tables, 31st ed. — reference integral tables.
- Source: [`src/calculus/integrate.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/integrate.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
- Tests: [`tests/test_cherry_dilog.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_dilog.c)
- Tests: [`tests/test_cherry_ei.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_ei.c)
- Tests: [`tests/test_cherry_li.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_li.c)
- Tests: [`tests/test_cherry_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_stress.c)

## Notes & additional examples

### Notes

`Integrate[f, x]` computes the indefinite integral via a cascade: Bronstein's rational-function algorithm, then the Risch–Norman (`pmint`) heuristic, then the lazy-loaded CRC integral tables; `Method -> "<name>"` pins a single subroutine. Antiderivatives are returned without an integration constant and are not always simplified — for example `Integrate[Sin[x], x]` returns `-(1 + Cos[x])` rather than `-Cos[x]`. Definite integration is **not** supported in this build: `Integrate[x^2, {x, 0, 1}]` threads `Integrate` over the bound list and returns the garbage form `{1/3 x^3, Integrate[x^2, 0], Integrate[x^2, 1]}` instead of `1/3`. Restrict use to the indefinite, single-variable form.
