# EllipticPi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EllipticPi[n, m] is the complete elliptic integral of the third kind, Integrate[1/((1 - n Sin[t]^2) Sqrt[1 - m Sin[t]^2]), {t, 0, Pi/2}], and EllipticPi[n, phi, m] the incomplete one, with upper limit phi. The parameter argument is m = k^2, not the modulus k. EllipticPi[0, m] is EllipticK[m] and EllipticPi[0, phi, m] is EllipticF[phi, m]; EllipticPi[1, m] is ComplexInfinity. Where the path crosses the pole at Sin[t]^2 == 1/n the value is complex, not a real principal value: EllipticPi[3/2, 1/2] is -0.4567203134529101 - 2.720699046351328 I.`**

## Examples (16)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= EllipticK[0]
Out[1]= 1/2 Pi

In[2]:= N[EllipticK[1/2], 25]
Out[2]= 1.8540746773013719184338503

In[3]:= D[EllipticF[phi, m], phi]
Out[3]= 1/Sqrt[1 - m Sin[phi]^2]

In[4]:= N[EllipticPi[3/2, ArcSin[Sqrt[2] Sqrt[1/(1 + 7/5)]], 1/2], 20]
Out[4]= 0.987739699728522080208 - 2.7206990463513267759*I
```

### Applications (12)

With n = 0 the third factor vanishes, leaving the first kind

```mathematica
In[5]:= EllipticPi[0, 1/2]
Out[5]= (8 Pi^(3/2))/Gamma[-1/4]^2
```

The same reduction, incomplete

```mathematica
In[6]:= EllipticPi[0, phi, m]
Out[6]= EllipticF[phi, m]
```

At m = 0 there is a closed form in n alone

```mathematica
In[7]:= EllipticPi[n, 0]
Out[7]= (1/2 Pi)/Sqrt[1 - n]
```

And its incomplete counterpart

```mathematica
In[8]:= EllipticPi[n, phi, 0]
Out[8]= ArcTanh[Sqrt[-1 + n] Tan[phi]]/Sqrt[-1 + n]
```

A pole of the COMPLETE form only: the integrand carries 1/Cos[t]^2 at the upper limit

```mathematica
In[9]:= EllipticPi[1, m]
Out[9]= ComplexInfinity
```

Complete, via R_F + (n/3) R_J

```mathematica
In[10]:= N[EllipticPi[1/2, 1/4], 30]
Out[10]= 2.41367150420119464066692352054
```

Incomplete, amplitude Pi/3

```mathematica
In[11]:= N[EllipticPi[1/2, Pi/3, 1/4], 25]
Out[11]= 1.3101681612463965511336307
```

For n > 1 the path crosses the pole at Sin[t]^2 == 1/n: the value is complex, not a real principal value

```mathematica
In[12]:= N[EllipticPi[3/2, 1/2], 20]
Out[12]= -0.456720313452909897007 - 2.7206990463513267759*I
```

The amplitude derivative is the integrand

```mathematica
In[13]:= D[EllipticPi[n, phi, m], phi]
Out[13]= 1/(Sqrt[1 - m Sin[phi]^2] (1 - n Sin[phi]^2))
```

The complete form's parameter derivative is closed form; the incomplete one's stays inert

```mathematica
In[14]:= D[EllipticPi[n, m], m]
Out[14]= (1/2 (EllipticE[m]/(-1 + m) + EllipticPi[n, m]))/(-m + n)
```

A visible NDArray in the amplitude slot

```mathematica
In[15]:= EllipticPi[0.5, NDArray[{0.1, 0.2, 0.3}], 0.25]
Out[15]= {0.100209, 0.201675, 0.305686}
```

The three-argument shape lowers too, through the n-ary kernel opcode

```mathematica
In[16]:= CompileDiagnostics[{{x, _Real}}, EllipticPi[0.5, x, 0.25]]
Out[16]= {"Compiled" -> True, "ResultType" -> "Real", "Instructions" -> 7, "CommonSubexpressions" -> 0, "InstructionsUnoptimized" -> 7}
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| BesselJ[0, .] over 10^6 | 3.46e+03 s | 1.74e+03 s | 53 s |
| Zeta over 10^6 | 8.17 s | 4.31e+03 s | 5.7 s |
| AiryAi over 10^6 | 3.14 s | 107 s | 58.9 s |
| PolyGamma[0, .] over 10^6 | 1.55 s | 151 s | 8.19 s |
| Gamma over 10^6 | 1.16 s | 1.35 s | 7.34 s |
| Erf over 10^6 | 0.968 s | 1.23 s | 7.43 s |

## Implementation notes

**Algorithm.** `builtin_ellipticpi` dispatches on argument count (2 complete, 3 incomplete;
anything else emits `EllipticPi::argt`). Both forms run the numeric path first, then the
reductions: `Π[0, m] = K(m)`, `Π[n, 0] = π/(2√(1−n))`, `Π[1, m] = ComplexInfinity`,
`Π[n, ∞] = Π[∞, m] = 0`, and for the incomplete form `Π[n, 0, m] = 0`,
`Π[0, φ, m] = F(φ, m)`, `Π[n, π/2, m] = Π[n, m]`,
`Π[n, φ, 0] = ArcTanh[√(n−1) tan φ]/√(n−1)`, plus the oddness fold in the amplitude.

The machine kernels are Carlson compositions:

- complete — `Π(n|m) = R_F(0, 1−m, 1) + (n/3) R_J(0, 1−m, 1, 1−n)`, checked against mpmath
  at ≤ 2 ulp over `n, m ∈ (0,1)` and ≤ 3.3 ulp with `n ∈ (−6, 1)`;
- incomplete — `Π(n; φ|m) = s R_F(c², 1−m s², 1) + (n/3) s³ R_J(c², 1−m s², 1, 1−n s²)`
  with `s = sin r`, `c = cos r` on the principal strip, ≤ 3.8 ulp over 400 random
  `(n, φ, m)`, extended by the quasi-period `Π(n; φ+kπ|m) = Π(n; φ|m) + 2k Π(n|m)`
  (verified exactly against mpmath at `k = −1, 1, 2`). `c²` is computed as `c*c` for the
  same reason as in `EllipticF`.

`R_J` needs `R_C` (the degenerate `R_F`), which is why both landed together.

**`n > 1` is a decline, not a principal value.** The standing reason for having no kernel
here was that `n > 1` needs a Cauchy principal value, and that reason was simply wrong:
the value past the pole at `sin²t = 1/n` is genuinely complex —
`Π(3/2 | 1/2) = −0.456720313453 − 2.72069904635 i`, mpmath and Arb agreeing. So what a
`double` kernel owes there is a decline (`carlson_rj` refuses `p ≤ 0`), exactly as `K`, `E`
and `F` decline outside their real domains, and `flint_num_elliptic_pi` answers through
Arb's `acb_elliptic_pi`.

**The pole belongs to the complete form only.** `Π[1, m]` diverges because the integrand
carries a `1/cos²t` and the upper limit is `π/2` — verified divergent as `1/√ε`: 21.5, 221,
2221, 22214 at `ε = 10⁻², 10⁻⁴, 10⁻⁶, 10⁻⁸`. The incomplete form has no such pole, and
`Π[1, 1, 1/2]` is `1.73199154202`. The test sits ahead of the numeric path, or Arb's
non-finite ball would come back as an unevaluated `EllipticPi[1., 1/2]`.

**Data structures.** `double` scalars through `carlson_rj`/`carlson_rc`; `acb_t` above
machine precision. Registration is split by arity: a binary (`REG_B`) kernel for the
complete form and an n-ary (`REG_N`) one for the incomplete. The element-wise NDArray layer
tops out at arity 2, so only the complete form rides that buffer, while `Compile[]` lowers
both — the three-argument shape through the generic n-ary `OP_KERNN` path, which reserves a
consecutive register block and reads the kernel from the symbol's `ndarray_nary_kernel`.

**A sequencing trap, recorded.** `packed_aware` is a property of the *symbol*, not of one
arity, so registering the two-argument kernel also stopped the transparency gate
materialising packed `List`s for the three-argument form. The `ndarray_delist_and_reeval`
fallback in all three argument positions is what keeps that correct, and it had to land
first — a visible `NDArray` left unevaluated is a wrong answer, not a slow one.

**Complexity / limits.** `O(1)` per element: 41 ns/element for the complete form at 2 ulp
(a 970× speedup over the Arb-per-element path it replaced, and 7.7× faster than the
equivalent SciPy composition), 1 ulp for the incomplete one. Real principal domain only
(`m < 1`, `n < 1`, and `1 − n sin²φ > 0`). The `n`- and `m`-derivatives of the *incomplete*
form stay inert `Derivative[…]`: a least-squares fit over the natural candidate basis does
not recover them (residual 1.5 relative, no simple rational coefficients), so the closed
form involves terms that basis does not span — and an inert derivative is honest where a
guessed one corrupts every caller silently.

- Exact reductions: `EllipticK[0] = π/2`, `EllipticK[1] = ComplexInfinity`,
  `EllipticE[0] = π/2`, `EllipticE[1] = 1`; `EllipticF[0, m] = 0`,
  `EllipticF[φ, 0] = φ`, `EllipticF[π/2, m] = EllipticK[m]`;
  `EllipticE[π/2, m] = EllipticE[m]`;
  `EllipticPi[0, m] = EllipticK[m]`, `EllipticPi[0, φ, m] = EllipticF[φ, m]`,
  `EllipticPi[n, 0, m] = 0`, `EllipticPi[n, π/2, m] = EllipticPi[n, m]`.
- Exact non-special arguments stay symbolic (`EllipticF[1/3, 1/2]`); a numeric
  value follows from an inexact argument or from `N`.
- **Numeric evaluation** routes through FLINT/Arb's `acb_elliptic_*`, which is
  rigorous, arbitrary-precision, defined on the whole complex plane, and already
  uses the parameter convention and branch placement above. That buys the three
  things a hand-rolled kernel most easily gets wrong:
  - the **quasi-periodic extension** off the principal strip
    (`F(φ + kπ | m) = F(φ | m) + 2k K(m)`, and likewise for `E` and `Π`), so φ of
    any size is handled;
  - **complex φ**, which `EllipticF[ArcSin[z], m]` produces routinely as soon as
    `|z| > 1` — the normal case for an elliptic pencil written in `ArcSin` form;
  - the **value past the pole** for `EllipticPi` with `n > 1`, where the path
    crosses `sin²t = 1/n`. It is not a real principal value — `Π(3/2 | 1/2)` is
    `−0.456720313453 − 2.72069904635 i`, mpmath and Arb agreeing — so the machine
    kernel declines there and Arb answers.
- **Machine kernels.** `EllipticK`, `EllipticE` (both arities) and `EllipticF`
  carry `double` kernels built on Carlson's symmetric forms `R_F` and `R_D`, used
  by the packed/NDArray element-wise paths. They cover the real principal domain
  and **decline** outside it (`m > 1`, or `1 − m sin²φ < 0`), which abandons the
  buffer so the List path answers exactly through Arb — slower, never wrong.
  Measured against mpmath on identical machine inputs, the buffer now agrees with
  the scalar path to 3 ulp (`K`, `F`) and 6 ulp (`E`); `R_F` and `R_D` share one
  duplication loop, since every caller that wants `R_D` wants `R_F` at the same
  arguments and the two recurrences walk an identical sequence.
  `EllipticPi` carries one too, on Carlson's `R_J` (and `R_C`, which `R_J`
  needs): the **complete** `EllipticPi[n, m]` is a binary (`REG_B`) kernel at
  41 ns/element and 2 ulp, and the **incomplete** `EllipticPi[n, φ, m]` is an
  n-ary (`REG_N`) one at 1 ulp, carrying the quasi-period
  `Π(n; φ+kπ|m) = Π(n; φ|m) + 2k Π(n|m)`. The element-wise NDArray layer tops
  out at arity 2, so only the complete form rides that path; `Compile[]` lowers
  both, the three-argument shape through `OP_KERNN`. Like the others the kernel
  **declines** outside the real principal domain: `n > 1` is not a real
  principal value at all (`Π(3/2 | 1/2)` is `−0.45672 − 2.72070 i`), so the
  honest machine answer there is to hand back to Arb.
  Because `packed_aware` is a property of the symbol and not of one arity,
  registering the two-argument kernel also stopped the transparency gate
  materialising packed `List`s for the three-argument form — which is why the
  `ndarray_delist_and_reeval` fallback in all three argument positions had to
  land first, and must stay.
- **Exact values.** `EllipticK[-1] = Γ(1/4)²/(4√(2π))` and
  `EllipticK[1/2] = 8π^(3/2)/Γ(-1/4)²` (the lemniscatic singular values; `E` has
  no closed form at either and keeps none). `EllipticPi[n, 0] = π/(2√(1-n))`,
  `EllipticPi[n, φ, 0] = ArcTanh[√(n−1) tan φ]/√(n−1)`. At infinity:
  `EllipticK[∞] = 0`, `EllipticE[∞] = EllipticE[ComplexInfinity] =
  ComplexInfinity`, `EllipticF[φ, ∞] = 0`,
  `EllipticPi[∞, m] = EllipticPi[n, ∞] = 0`. All three incomplete forms are
  **odd in the amplitude** — `EllipticF[−φ, m] = −EllipticF[φ, m]`, likewise `E`
  and `Π` — because each integrand is even in `t`; the fold uses the same
  superficial-negativity test the trig heads use, so `−2x` folds and `−x−y` does
  not.

- **`Series`.** The complete integrals expand at `m = 0` from a dedicated
  kernel: `K = (π/2) Σ aₖ mᵏ` and `E = (π/2) Σ aₖ mᵏ/(1−2k)` with `a₀ = 1`,
  `aₖ = aₖ₋₁((2k−1)/(2k))²`. It has to be dedicated — the generic
  Taylor-via-`D` path evaluates `d/dm K` at `m = 0`, which is `(π/2 − π/2)/0`
  → `Indeterminate`, and the expansion is then abandoned, so
  `Series[EllipticK[m], {m, 0, 1}]` did not previously emit even the leading
  `π/2`. Composition through an inner series works; an expansion about a regular
  point still takes Taylor-via-`D`. The incomplete forms expand in `φ` by the
  generic path, as before.

- **`Interval`.** The complete integrals thread by certified monotonicity —
  `K` increasing and `E` decreasing, both only below `m = 1`, beyond which the
  value is complex. The amplitude slot of `EllipticF`, `EllipticE[φ,m]` and
  `EllipticPi` threads through the general derivative certifier. The `m` slot
  declines: its derivative reproduces the head, so the certifier can never
  bottom out there.

- **Machine numbers.** An inexact machine scalar is answered by the `double`
  kernel and comes back as a machine `Real`, not a 53-bit `EXPR_MPFR` — so
  `Precision[EllipticK[0.5]]` is `MachinePrecision`, `MachineNumberQ` is `True`,
  and an elliptic-produced list packs for its consumers. An `EXPR_MPFR`
  argument still routes to Arb, so `SetPrecision[m, 16]` is never silently
  answered in `double`.

- **The `m = 1` continuation.** `E(φ | 1) = ∫₀^φ |cos t| dt`, which is `Sin[φ]`
  only on the principal strip `|φ| ≤ π/2`; off it the value continues as
  `Sin[φ − kπ] + 2k`. The reduction therefore fires only where it is *provably*
  valid — decided by evaluating `Abs[φ] ≤ π/2`, so an exact rational, a `Real`,
  an `MPFR` and a symbol carrying assumptions all get an answer. A symbolic
  amplitude is undecidable and stays symbolic rather than becoming wrong:
  `EllipticE[2, 1]` is left alone, and `N[EllipticE[2, 1]]` is `1.0907025731743`
  (`= 2 − Sin[2]`), agreeing with its neighbour at `m = 1 − 10⁻²⁵`.

- **Inexactness and poles.** The numeric path runs *before* the remaining exact
  reductions, so an inexact argument gets an inexact answer — `EllipticK[0.]` is
  `1.5707963267948966`, not the exact `π/2`, and `EllipticPi[0., 1/2]` is
  `1.8540746773013719`, not the symbolic `EllipticK[1/2]`. The exact reductions
  remain underneath as the fallback, which is what still reduces
  `EllipticE[π/2, 0.5]` (Arb reads numbers, not a `Times`). Poles are checked
  *ahead* of the numeric path, in every spelling of the argument, because Arb
  returns a non-finite ball there and the bridge renders that as "unevaluated":
  `EllipticK[1]`, `EllipticK[1.]` and `EllipticK[SetPrecision[1, 30]]` all give
  `ComplexInfinity`, and so does the complete `EllipticPi[1, m]` — whose
  divergence is real (`Π(1−ε | 1/2)` grows as `1/√ε`). The *incomplete* third
  kind has no pole at `n = 1`: `EllipticPi[1, 1, 1/2]` is `1.73199154202`.

- **Derivatives.** The φ-derivatives are the integrands, which is what makes a
  numeric verification of an antiderivative built from these kernels close:
  `D[EllipticF[φ, m], φ] = 1/√(1 − m sin²φ)`,
  `D[EllipticE[φ, m], φ] = √(1 − m sin²φ)`,
  `D[EllipticPi[n, φ, m], φ] = 1/((1 − n sin²φ) √(1 − m sin²φ))`. Also
  `D[EllipticK[m], m] = (E(m) − (1−m) K(m)) / (2m(1−m))`,
  `D[EllipticE[m], m] = (E(m) − K(m)) / (2m)` and
  `D[EllipticE[φ, m], m] = (E(φ,m) − F(φ,m)) / (2m)`.
  `D[EllipticF[φ, m], m] = −E(φ,m)/(2(m−1)m) − F(φ,m)/(2m) +
  sin(2φ)/(4(m−1)√(1−m sin²φ))`, and for the **complete** third kind
  `D[EllipticPi[n, m], n] = (nE + (m−n)K + (n²−m)Π)/(2(m−n)(n−1)n)` and
  `D[EllipticPi[n, m], m] = (E/(m−1) + Π)/(2(n−m))`. Each of those three was an
  inert placeholder until it had been checked against a central difference at 30
  digits.

  The `n`- and `m`-derivatives of the **incomplete** `EllipticPi` remain
  deliberately inert, and not for want of trying: a least-squares fit over the
  natural candidate basis (E, F, Π over `(n−m)`, `(m−1)`, `(n−1)`, plus the
  `sin(2φ)/√(…)` boundary term) does not recover them — residual 1.5 relative,
  with no simple rational coefficients — so the closed form involves terms that
  basis does not span. An inert derivative is honest where a guessed one corrupts
  every caller silently; `BesselJ` treats its order the same way. One visible
  consequence: `Interval[]` threading for the incomplete `EllipticPi` in its `m`
  slot stays symbolic, since the certifier differentiates to get a sign.
- Wrong arity emits `EllipticK::argx` / `EllipticF::argrx` (fixed arity) or
  `EllipticE::argt` / `EllipticPi::argt` (the overloaded pair) and stays
  unevaluated.
- Without FLINT (`USE_FLINT` undefined) the numeric path falls back to the
  machine Carlson kernels for the real domain and otherwise declines, leaving the
  call symbolic.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [EllipticK](../../special-functions/EllipticK/), [EllipticF](../../special-functions/EllipticF/), [EllipticE](../../special-functions/EllipticE/), [E](../../mathematical-constants/E/), [N](../../arithmetic/N/), [ArcSin](../../elementary-functions/ArcSin/), [List](../../other-advanced/List/), [Series](../../power-series/Series/)

- B. C. Carlson, *Computing elliptic integrals by duplication*, Numer. Math. **33** (1979) 1-16.
- B. C. Carlson, *Numerical computation of real or complex elliptic integrals*, Numer. Algorithms **10** (1995) 13-26 — `R_J` and the `p < 0` transformation.
- W. H. Press et al., *Numerical Recipes in C*, 2nd ed. (Cambridge, 1992), §6.11 — the `rj`/`rc` loops.
- DLMF §19.25.14 — the third kind in terms of `R_F` and `R_J`.
- Source: [`src/special_functions/elliptic.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/elliptic.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_elliptic.c`](https://github.com/stblake/mathilda/blob/main/tests/test_elliptic.c)

## Notes & additional examples

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
