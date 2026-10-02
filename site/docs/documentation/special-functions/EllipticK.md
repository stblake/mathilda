# EllipticK

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EllipticK[m] is the complete elliptic integral of the first kind, Integrate[1/Sqrt[1 - m Sin[t]^2], {t, 0, Pi/2}]. The argument is the PARAMETER m = k^2, not the modulus k. EllipticK[0] is Pi/2 and EllipticK[1] is ComplexInfinity; exact arguments otherwise stay symbolic and inexact ones evaluate numerically at their precision.`**

## Examples (15)

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

### Applications (11)

At m = 0 the integrand is 1, so K is the quarter period Pi/2

```mathematica
In[5]:= EllipticK[0]
Out[5]= 1/2 Pi
```

The integrand loses its last denominator at t = Pi/2

```mathematica
In[6]:= EllipticK[1]
Out[6]= ComplexInfinity
```

A lemniscatic singular value: closed form in Gamma, not a decimal

```mathematica
In[7]:= EllipticK[-1]
Out[7]= (1/4 Gamma[1/4]^2)/Sqrt[2 Pi]
```

The second singular value; EllipticE has no closed form at either

```mathematica
In[8]:= EllipticK[1/2]
Out[8]= (8 Pi^(3/2))/Gamma[-1/4]^2
```

The same number to 30 digits, through Arb

```mathematica
In[9]:= N[EllipticK[1/2], 30]
Out[9]= 1.854074677301371918433850347195
```

A dedicated kernel: the generic Taylor path hits 0/0 here

```mathematica
In[10]:= Series[EllipticK[m], {m, 0, 3}]
Out[10]= 1/2 Pi + 1/8 Pi m + 9/128 Pi m^2 + 25/512 Pi m^3 + O[m]^4
```

Past the branch point at m = 1 the value is complex

```mathematica
In[11]:= N[EllipticK[9/5], 20]
Out[11]= 1.41933775128651291929 - 1.348846512193268578*I
```

Here dK/dm ~ 2^56, and all 20 digits are still right

```mathematica
In[12]:= N[EllipticK[1 - 10^-17], 20]
Out[12]= 20.9582676515692789829
```

The parameter derivative, in terms of E and K themselves

```mathematica
In[13]:= D[EllipticK[m], m]
Out[13]= (1/2 (EllipticE[m] - EllipticK[m] (1 - m)))/(m (1 - m))
```

Certified: K is increasing in m below 1

```mathematica
In[14]:= EllipticK[Interval[{1/4, 1/2}]]
Out[14]= Interval[{EllipticK[1/4], (8 Pi^(3/2))/Gamma[-1/4]^2}]
```

Real in, Real out -- so the result packs

```mathematica
In[15]:= Precision[EllipticK[0.5]]
Out[15]= MachinePrecision
```

## Implementation notes

**Algorithm.** `builtin_elliptick` tries, in this order: the pole (`m == 1` →
`ComplexInfinity`), an exact zero (`π/2`), infinity (`K(m) ~ (π − i log 16m)/(2√m)`,
so the magnitude vanishes and the answer is `0`), then the two lemniscatic singular
values `m = 1/2` and `m = −1` as closed forms in `Γ` (`ell_k_closed_half`,
`ell_k_closed_minus_one`). Only then does an inexact argument reach a number:
`elliptic_machine_k` evaluates `K(m) = R_F(0, 1−m, 1)` in `double` and the result is
returned as an `EXPR_REAL`, so `Precision[EllipticK[0.5]]` is `MachinePrecision` and the
value packs. If the kernel declines — `m ≥ 1`, where the value is a pole or genuinely
complex — `flint_num_elliptic_k` answers through Arb's `acb_elliptic_k`.

**Order is load-bearing, twice.** The pole is tested *before* the numeric path because
Arb returns a non-finite ball there and the bridge maps that to `NULL`, which would leave
`EllipticK[1.]` unevaluated instead of `ComplexInfinity`; `ell_is_one` therefore also
knows `EXPR_MPFR`, so `SetPrecision[1, 30]` and `1.` cannot disagree about being at a
pole. Conversely the numeric path runs *before* the remaining exact reductions, so an
inexact argument gets an inexact answer (`EllipticK[0.]` is `1.5707963267948966`, not the
exact `π/2`) — the reductions stay underneath as the fallback for spellings Arb cannot
read, such as an exact `Pi/2` upper limit.

**Data structures.** Plain `double` scalars: `carlson_rf` runs the duplication
`x,y,z → (x+λ)/4` with `λ = √x√y + √y√z + √z√x`, stopping when every relative deviation
from the mean is below `EC_ERRTOL_RF = 0.0025`, then adds the fifth-order tail in the
`E₂`/`E₃` invariants. The tolerance is not cosmetic: the tail is fifth order, so the error
goes as `(q/A)⁶` and the previous `0.01` bought 1e-12 where the comment claimed 1e-16.
The arbitrary-precision path is `acb_t` throughout, with the bridge's accuracy ladder
(`nb_eval1`) doubling the working precision while `acb_rel_accuracy_bits` falls short of
the digits requested.

**Complexity / limits.** Seven duplication steps at most, so the machine path is `O(1)`
with ~3 ulp error and no allocation; the `Listable`/NDArray path runs it element-wise
across cores. The machine kernel is real-domain only (`m < 1`); everything else — complex
`m`, `m > 1`, and every request above machine precision — is Arb's, which is rigorous but
single-threaded. `Series` at `m = 0` comes from a dedicated kernel rather than
Taylor-via-`D`, because `d/dm K` at `m = 0` is `(π/2 − π/2)/0`.

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

**See also:** [EllipticF](../../special-functions/EllipticF/), [EllipticE](../../special-functions/EllipticE/), [EllipticPi](../../special-functions/EllipticPi/), [E](../../mathematical-constants/E/), [N](../../arithmetic/N/), [ArcSin](../../elementary-functions/ArcSin/), [List](../../other-advanced/List/), [Series](../../power-series/Series/)

- B. C. Carlson, *Computing elliptic integrals by duplication*, Numer. Math. **33** (1979) 1-16.
- B. C. Carlson, *Numerical computation of real or complex elliptic integrals*, Numer. Algorithms **10** (1995) 13-26.
- W. H. Press et al., *Numerical Recipes in C*, 2nd ed. (Cambridge, 1992), §6.11 — the `rf`/`rd` duplication loops.
- M. Abramowitz and I. A. Stegun, *Handbook of Mathematical Functions* (Dover, 1964), ch. 17.
- DLMF §19.25 — Legendre integrals in terms of Carlson's symmetric forms.
- Source: [`src/special_functions/elliptic.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/elliptic.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_elliptic.c`](https://github.com/stblake/mathilda/blob/main/tests/test_elliptic.c)

## Notes & additional examples

### Notes

`EllipticK[m]` takes the **parameter** `m = k²`, not the modulus `k` — the Wolfram
convention, and the one place a reader of these signatures goes silently wrong.

Three argument classes are answered by three different routes. An exact `0`, `1`, `∞`,
`1/2` or `−1` reduces symbolically; any other exact argument stays symbolic, since a
decimal would be a loss of information rather than an evaluation. An inexact argument
inside the real principal domain `m < 1` goes to the `double` Carlson kernel and comes
back a machine `Real`. Everything else — `m > 1`, complex `m`, or any request above
machine precision — goes to Arb, which is rigorous and carries the branch placement.

The conditioning example is worth reading twice: `dK/dm ≈ 1/(2(1−m))`, which at
`m = 1 − 10⁻¹⁷` is about `2⁵⁶`, so an arbitrary-precision request has to be *planned* for
the loss rather than merely issued at the precision asked for. How much a function loses
is set by its conditioning, not by how its argument is spelled.
