# Symbolic summation

The [calculus tutorial](08-calculus.md) introduced `Sum` for *finite* sums —
adding `expr` as an index runs over a finite range. This tutorial takes the far
harder step of *infinite* summation: closing an entire series
`Sum[expr, {k, 1, Infinity}]` into a symbolic constant. That is one of the
oldest and deepest problems in analysis, and it is where a computer algebra
system earns its keep.

Mathilda's `Sum` is not one algorithm but a *cascade* of them, implemented
natively in C under `src/sum/` — a dispatcher (`sum.c`) that tries, in turn,
polynomial (Faulhaber), geometric, Gosper telescoping, rational–`Zeta`,
alternating–`LerchPhi`, Euler-sum, trigonometric–Fourier, and
hypergeometric summation. Each method is also exposed as a context-qualified
builtin (``Sum`Gosper``, ``Sum`Rational``, …) so you can see which one fired. This
mirrors the way `Integrate` layers its methods.

Every transcript below was produced by the actual Mathilda binary. Type the
`In[...]` lines yourself and you will get the same `Out[...]`. As always the
output form may be arranged differently from a textbook (`1/6 Pi^2` rather than
`π²/6`), but it is always mathematically correct — and, as one example will
show, occasionally *more* correct than the textbook.

## From finite to infinite

A finite sum with numeric bounds is just arithmetic:

```mathematica
In[1]:= Sum[k^2, {k, 1, 10}]
Out[1]= 385
```

Give `Sum` a *symbolic* upper bound and it returns a closed formula — the
[Faulhaber polynomials](https://en.wikipedia.org/wiki/Faulhaber%27s_formula),
which Mathilda builds from Bernoulli numbers:

```mathematica
In[1]:= Sum[k, {k, 1, n}]
Out[1]= 1/2 n (1 + n)

In[2]:= Sum[k^3, {k, 1, n}]
Out[2]= 1/4 n^2 (1 + n)^2
```

An *infinite* sum is the limit of these partial sums as `n -> Infinity`. When
that limit exists and has a name, Mathilda finds it. The rest of this tutorial
is a tour of the machinery that turns each family of series into a constant.

## Telescoping series and Gosper's algorithm

The gentlest infinite sums are *telescoping*: consecutive terms cancel, leaving
only the ends. The classic is `1/(k(k+1))`. Its partial fraction split is the
whole trick:

```mathematica
In[1]:= Apart[1/(k*(k + 1)), k]
Out[1]= 1/k - 1/(1 + k)
```

Each `-1/(k+1)` cancels the `1/k` of the next term, so the partial sum collapses
to `1 - 1/(n+1)`, and the infinite sum is `1`:

```mathematica
In[1]:= Sum[1/(k*(k + 1)), {k, 1, n}]
Out[1]= 1 - 1/(1 + n)

In[2]:= Sum[1/(k*(k + 1)), {k, 1, Infinity}]
Out[2]= 1
```

This idea goes back to [Mengoli](https://en.wikipedia.org/wiki/Pietro_Mengoli)
in the 1650s, but the general question — *given a term, is there a telescoping
antidifference at all?* — was only settled by
[**R. W. Gosper** in 1978](https://en.wikipedia.org/wiki/Gosper%27s_algorithm).
Gosper's algorithm decides, for any *hypergeometric* term `t(k)` (one whose ratio
`t(k+1)/t(k)` is a rational function of `k`), whether a hypergeometric `T(k)`
exists with `T(k+1) - T(k) = t(k)`, and constructs it when it does. Mathilda's
``Sum`Gosper`` module is a direct implementation; it is what closes the shifted
half-integer series that needs the digamma function to telescope:

```mathematica
In[1]:= Sum[1/(k*(2*k + 1)), {k, 1, Infinity}]
Out[1]= 2 - 2 Log[2]
```

A higher-order telescoping example combines partial fractions with the Basel
constant. Splitting `1/(k²(k+1)²)` produces both `1/k²` pieces (which sum to
`ζ(2)`) and `1/k` pieces (which telescope):

```mathematica
In[1]:= Apart[1/(k^2*(k + 1)^2), k]
Out[1]= 1/k^2 - 2/k + 1/(1 + k)^2 + 2/(1 + k)

In[2]:= Sum[1/(k^2*(k + 1)^2), {k, 1, Infinity}]
Out[2]= -3 + 1/3 Pi^2
```

The `-2/k` and `2/(k+1)` telescope to `2`, while the two inverse-square pieces
each contribute `ζ(2) = π²/6`; assembling them gives `π²/3 - 3`.

## Geometric and arithmetico-geometric series

A geometric series `Sum[r^k]` converges to `1/(1-r)` for `|r| < 1`; Mathilda
knows this even symbolically. Differentiating the geometric series with respect
to the ratio multiplies each term by `k`, producing *arithmetico-geometric*
sums. Mathilda handles the whole family, of which these are the simplest cases:

```mathematica
In[1]:= Sum[1/2^k, {k, 1, Infinity}]
Out[1]= 1

In[2]:= Sum[k/2^k, {k, 1, Infinity}]
Out[2]= 2

In[3]:= Sum[k^2/2^k, {k, 1, Infinity}]
Out[3]= 6
```

The pattern — `1, 2, 6` — is exactly what the derivative operator
`x d/dx` produces when applied repeatedly to `x/(1-x)` and evaluated at
`x = 1/2`. With a symbolic ratio the closed forms appear directly:

```mathematica
In[1]:= Sum[k*x^k, {k, 1, Infinity}]
Out[1]= x/(1 - x)^2
```

The exponential series is the limiting case where the "ratio" is `1/k`, giving
[Euler's number](https://en.wikipedia.org/wiki/E_(mathematical_constant)):

```mathematica
In[1]:= Sum[1/k!, {k, 0, Infinity}]
Out[1]= E
```

## The Zeta family: the Basel problem and its relatives

The most celebrated infinite sum in mathematics is the **Basel problem**, posed
by Mengoli in 1650 and left open until the 28-year-old
[Euler solved it in 1734](https://en.wikipedia.org/wiki/Basel_problem):

```mathematica
In[1]:= Sum[1/k^2, {k, 1, Infinity}]
Out[1]= 1/6 Pi^2
```

Euler went on to evaluate every even power, `Sum[1/k^(2m)]`, as a rational
multiple of `π^(2m)`. These are the even values of the
[Riemann zeta function](https://en.wikipedia.org/wiki/Riemann_zeta_function),
and Mathilda produces them the same way Euler did — from the Bernoulli numbers:

```mathematica
In[1]:= Sum[1/k^4, {k, 1, Infinity}]
Out[1]= 1/90 Pi^4

In[2]:= Sum[1/k^3, {k, 1, Infinity}]
Out[2]= Zeta[3]
```

`ζ(3)` (Apéry's constant) is left as `Zeta[3]` because — unlike the even
values — no closed form in terms of `π` is known; it was only proved
*irrational* by [Apéry in 1979](https://en.wikipedia.org/wiki/Ap%C3%A9ry%27s_theorem).

The engine that produces these is ``Sum`Rational``. Its master identity splits any
rational summand into partial fractions and evaluates each piece as a value or
derivative of the digamma/`Zeta` functions. When the denominator has
**complex-conjugate roots**, the two conjugate digamma terms recombine into a
real hyperbolic cotangent — which is how Mathilda evaluates
[Euler's cotangent series](https://en.wikipedia.org/wiki/Mittag-Leffler%27s_theorem):

```mathematica
In[1]:= Sum[1/(k^2 + 1), {k, 1, Infinity}]
Out[1]= 1/2 (-1 + Pi Coth[Pi])
```

Restricting to *odd* denominators, or attaching an alternating sign `(-1)^k`,
produces the **Dirichlet** relatives handled by ``Sum`Alternating`` (via the
Lerch transcendent and Dirichlet beta/eta functions). Two milestones of the
subject appear immediately — the sum of inverse odd squares, and
**Catalan's constant**, introduced by
[Eugène Catalan in 1865](https://en.wikipedia.org/wiki/Catalan%27s_constant)
and still not known to be irrational:

```mathematica
In[1]:= Sum[1/(2*k - 1)^2, {k, 1, Infinity}]
Out[1]= 1/8 Pi^2

In[2]:= Sum[(-1)^k/(2*k + 1)^2, {k, 0, Infinity}]
Out[2]= Catalan
```

The same alternating machinery recovers the **Leibniz–Madhava series** for
`π/4` (known to the Kerala school by 1400, rediscovered by Leibniz in 1676) and
the alternating form of `ζ(3)`, which — by separating even and odd terms —
equals `3/4` of Apéry's constant:

```mathematica
In[1]:= Sum[(-1)^(k+1)/(2*k - 1), {k, 1, Infinity}]
Out[1]= 1/4 Pi

In[2]:= Sum[(-1)^(k - 1)/k^3, {k, 1, Infinity}]
Out[2]= 3/4 Zeta[3]
```

## Logarithms, digamma, and Fourier series

Not every convergent series lands on a power of `π`. The **Mercator series**
(Nicholas Mercator, 1668) — the alternating harmonic series — is the value of
`log(1 + x)` at `x = 1`:

```mathematica
In[1]:= Sum[x^k/k, {k, 1, Infinity}]
Out[1]= -Log[1 - x]

In[2]:= Sum[(-1)^(k + 1)/k, {k, 1, Infinity}]
Out[2]= Log[2]
```

More subtle are the **Fourier series** — sums of `Sin[k]/k` and `Cos[k]/k`.
These converge only *conditionally* and depend on the specific argument. The
sawtooth-wave series `Sum[Sin[k]/k]` evaluates at `x = 1`, and the companion
cosine series gives a logarithm of a sine; Mathilda's ``Sum`Trigonometric``
module recognises both by writing the summand as the imaginary/real part of
`z^k/k` with `z = e^i` and reusing the Mercator series:

```mathematica
In[1]:= Sum[Sin[k]/k, {k, 1, Infinity}]
Out[1]= 1/2 (-1 + Pi)

In[2]:= Sum[Cos[k]/k, {k, 1, Infinity}]
Out[2]= -(Log[2] + Log[Sin[1/2]])
```

The first is the value `(π - 1)/2` of the sawtooth `(π - x)/2` at `x = 1`; the
second is `-log(2 sin(1/2))`, the real part of `-log(1 - e^i)`.

## Euler sums and multiple zeta values

Weighting each term by a **harmonic number** `H_k = 1 + 1/2 + ... + 1/k`
produces the *Euler sums*, which Euler studied in 1775 and which are now
understood as [multiple zeta values](https://en.wikipedia.org/wiki/Multiple_zeta_value).
The linear Euler sums reduce entirely to ordinary zeta values by Euler's own
reflection formula `Sum H_k/k^q = (1 + q/2) ζ(q+1) - (1/2) Σ ζ(j+1) ζ(q-j)`,
which Mathilda's ``Sum`Euler`` module implements:

```mathematica
In[1]:= Sum[HarmonicNumber[k]/k^2, {k, 1, Infinity}]
Out[1]= 2 Zeta[3]

In[2]:= Sum[HarmonicNumber[k]/k^3, {k, 1, Infinity}]
Out[2]= 1/72 Pi^4

In[3]:= Sum[HarmonicNumber[k]/k^4, {k, 1, Infinity}]
Out[3]= 3 Zeta[5] - 1/6 Pi^2 Zeta[3]
```

`In[1]` is the famous identity `Σ H_k/k² = 2ζ(3)`. `In[2]` has even weight `4`,
so it collapses all the way to a power of `π`: `π⁴/72`.

*Nonlinear* Euler sums — with `H_k²` or a generalized harmonic number
`H_k^(2) = HarmonicNumber[k, 2]` — are much harder and were a research frontier
into the 1990s. Mathilda reduces the quadratic and generalized cases through the
[Borwein–Borwein–Girgensohn](https://en.wikipedia.org/wiki/Euler_sum)
double-zeta reduction:

```mathematica
In[1]:= Sum[HarmonicNumber[k]^2/k^2, {k, 1, Infinity}]
Out[1]= 17/360 Pi^4

In[2]:= Sum[HarmonicNumber[k, 2]/k^3, {k, 1, Infinity}]
Out[2]= -9/2 Zeta[5] + 1/2 Pi^2 Zeta[3]
```

The weight-5 result in `In[2]` is `3ζ(2)ζ(3) - (9/2)ζ(5)` (recall
`ζ(2) = π²/6`, so `3ζ(2) = π²/2`). Getting a symbolic engine to *simplify* into
this form, rather than leaving an opaque multiple-zeta expression, is exactly
the kind of task these identities were built for.

## Binomial sums

Series involving the [central binomial coefficient](https://en.wikipedia.org/wiki/Central_binomial_coefficient)
`Binomial[2k, k]` are notoriously delicate — inverse central-binomial sums are
tied to the [Beta-function integral](https://en.wikipedia.org/wiki/Beta_function)
and often produce `π` divided by `√3`. Mathilda's ``Sum`Hypergeometric`` module
rewrites `Binomial[2k, k]` through `Gamma` and recognises the resulting
`2F1` hypergeometric series:

```mathematica
In[1]:= Sum[2^k/Binomial[2*k, k], {k, 1, Infinity}]
Out[1]= 1 + 1/2 Pi

In[2]:= Sum[1/(k*Binomial[2*k, k]), {k, 1, Infinity}]
Out[2]= (1/3 Pi)/Sqrt[3]
```

`In[1]` is the arcsine-family sum `1 + π/2`; `In[2]` is the elegant
`π/(3√3)`, whose closed form comes from the arcsine series evaluated at `1/2`.

## The hypergeometric frontier: machines for π

The deepest series on this tour are the ones invented specifically to *compute*
`π` — and they sit right at the edge of what closed-form summation can do. Each
is a hypergeometric series, and Mathilda expresses them exactly as such
(`HypergeometricPFQ[...]`); confirming the classical constant is then a matter of
evaluating that hypergeometric numerically to as many digits as you like.

The [**Bailey–Borwein–Plouffe** formula (1995)](https://en.wikipedia.org/wiki/Bailey%E2%80%93Borwein%E2%80%93Plouffe_formula)
is famous for letting you compute the *n*-th hexadecimal digit of `π` without
the preceding ones. Summed in full it is exactly `π`:

```mathematica
In[1]:= N[Sum[(1/16^k)*(4/(8*k + 1) - 2/(8*k + 4) - 1/(8*k + 5) - 1/(8*k + 6)), {k, 0, Infinity}], 30]
Out[1]= 3.141592653589793238462643383279

In[2]:= N[Pi, 30]
Out[2]= 3.141592653589793238462643383279
```

[**Ramanujan** (1914)](https://en.wikipedia.org/wiki/Ramanujan%E2%80%93Sato_series)
found series of astonishing efficiency. This one adds roughly eight correct
decimal places per term and converges to `9801√2/(4π)`. Mathilda sums it to a
hypergeometric value; twenty-two digits confirm the closed form:

```mathematica
In[1]:= N[Sum[(1103 + 26390*k)*(4*k)!/((k!)^4*396^(4*k)), {k, 0, Infinity}], 22]
Out[1]= 1103.0000268319745734639

In[2]:= N[9801*Sqrt[2]/(4*Pi), 22]
Out[2]= 1103.0000268319745734639
```

One series *does* close symbolically all the way to `2/π` — a
[Ramanujan-type](https://en.wikipedia.org/wiki/Ramanujan%E2%80%93Sato_series)
cubic-binomial series that Mathilda's hypergeometric machinery reduces
completely:

```mathematica
In[1]:= Sum[(4*k + 1)*Binomial[2*k, k]^3/(-64)^k, {k, 0, Infinity}]
Out[1]= 2/Pi
```

!!! note "Symbolic versus numeric"
    When a series is provably a named constant, `Sum` returns it symbolically.
    When it lands on a hypergeometric value with no simpler closed form, `Sum`
    returns that hypergeometric — and wrapping it in `N[..., d]` gives you `d`
    correct digits. This is the same division of labour you met with
    `Integrate`: an exact answer when one exists, an honest special-function
    value otherwise.

## Where to next

You have seen every major method in Mathilda's summation cascade: Faulhaber
polynomials, Gosper telescoping, the rational–`Zeta` master identity with its
hyperbolic-cotangent branch, the Dirichlet/alternating family, Fourier series,
Euler sums and multiple zeta values, binomial/hypergeometric summation, and the
`π`-machine frontier.

- The companion [infinite products tutorial](12-infinite-products.md) does for
  `Product` what this one did for `Sum` — telescoping, Weierstrass
  factorizations, Euler prime products, and the Gamma/Glaisher constants.
- When a series has no closed form at all, the
  [numerical calculus tutorial](09-numerical-calculus.md) shows how `NSum`
  computes its value to arbitrary precision.
- The [special functions tutorial](10-special-functions.md) covers `Zeta`,
  `PolyGamma`, `HarmonicNumber`, and `HypergeometricPFQ` — the constants and
  functions these sums evaluate to.

A good habit while exploring: type `?Sum` (or ``?Sum`Gosper``) at the prompt to
read the built-in help without leaving the REPL.
