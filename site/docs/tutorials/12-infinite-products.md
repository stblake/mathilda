# Infinite products

An infinite product `Product[expr, {k, 1, Infinity}]` multiplies infinitely many
factors together. Where the [summation tutorial](11-symbolic-summation.md) turned
series into constants, this one turns products into constants — and the answers
are, if anything, more surprising: rational products collapse by telescoping,
trigonometric products fall out of the Weierstrass factorization of `sin` and
`cos`, products over the *primes* reproduce the zeta function, and a handful of
exponential products distil the constants `e`, `γ` (Euler–Mascheroni), and
Glaisher's `A` out of thin air.

Mathilda evaluates products with a native C cascade under `src/product/`: a
dispatcher (`product.c`) plus one module per family — rational telescoping,
infinite-rational (Gamma), Viète cosine, Cantor double-exponential, Euler prime
products, log-sum/Glaisher, geometric, and `q`-products. Each is also a
context-qualified builtin (``Product`Viete``, ``Product`EulerPrime``, …). `Product`
is `HoldAll`.

Every transcript was produced by the actual Mathilda binary — and two of the
examples below are worked cases where the "textbook" closed form circulating
online is subtly *wrong*; Mathilda (checked with `NProduct`) gets them right.

## From finite to infinite

A finite product with a symbolic bound closes to a factorial or a telescoping
form:

```mathematica
In[1]:= Product[k, {k, 1, n}]
Out[1]= Factorial[n]

In[2]:= Product[(k + 1)/k, {k, 1, n}]
Out[2]= 1 + n
```

`In[2]` is the archetype of telescoping: `(k+1)/k` has numerator and denominator
that cancel across neighbours, leaving only `(n+1)/1`. An *infinite* product is
the limit of these partial products — and convergence is stricter than for sums:
the factors must approach `1` fast enough that `Sum[Log[factor]]` converges.

## Rational telescoping products

The workhorse family is *rational* products, where each factor is a ratio of
polynomials in `k` that factor into shifted linear pieces. The
[Wikipedia telescoping product](https://en.wikipedia.org/wiki/Telescoping_series#Products)
`1 - 1/k²` is the canonical example: writing `1 - 1/k² = ((k-1)(k+1))/k²` makes
the cancellation visible, and the partial product is `(n+1)/(2n)`, tending to
`1/2`:

```mathematica
In[1]:= Product[1 - 1/k^2, {k, 2, n}]
Out[1]= (1/2 (1 + n))/n

In[2]:= Product[1 - 1/k^2, {k, 2, Infinity}]
Out[2]= 1/2
```

The same mechanism handles the quadratic and higher-degree cousins. A quadratic
numerator `1 - 2/(k(k+1))` factors as `(k-1)(k+2)/(k(k+1))`; the difference and
sum of cubes factor the cubic case; and a symmetric shift `k(k+2)/(k+1)²`
telescopes directly:

```mathematica
In[1]:= Product[1 - 2/(k*(k + 1)), {k, 2, Infinity}]
Out[1]= 1/3

In[2]:= Product[(k*(k + 2))/(k + 1)^2, {k, 1, Infinity}]
Out[2]= 1/2
```

For the cubic product, the factorizations Mathilda relies on are visible with
`Factor`:

```mathematica
In[1]:= Factor[k^3 - 1]
Out[1]= (-1 + k) (1 + k + k^2)

In[2]:= Factor[k^3 + 1]
Out[2]= (1 + k) (1 - k + k^2)
```

The linear parts `(k-1)/(k+1)` telescope while the quadratic parts
`(k²+k+1)/(k²-k+1)` shift by one index, so the whole product collapses to `2/3`:

```mathematica
In[1]:= Product[(k^3 - 1)/(k^3 + 1), {k, 2, Infinity}]
Out[1]= 2/3
```

When the roots become *complex*, telescoping is no longer purely rational — the
answer picks up a hyperbolic function through the Gamma reflection formula
`Γ(z)Γ(1-z) = π/sin(πz)`. Mathilda's ``Product`RationalInfinite`` module writes the
product as a ratio of Gamma values at the roots and applies the reflection:

```mathematica
In[1]:= Product[(k^2 - 1)/(k^2 + 1), {k, 2, Infinity}]
Out[1]= Pi Csch[Pi]
```

That is `π/sinh(π)`: the `k²-1` part telescopes to a rational, while the
`k²+1 = (k-i)(k+i)` part contributes `Γ`-values at `±i` that reflect into
`sinh(π)`.

## Trigonometric and hyperbolic factorizations

The deepest source of infinite products is
[**Euler's factorization** (1735)](https://en.wikipedia.org/wiki/Basel_problem#Euler's_approach)
of the sine as a product over its roots — later made rigorous by the
[Weierstrass factorization theorem](https://en.wikipedia.org/wiki/Weierstrass_factorization_theorem).
Writing the two forms side by side,

> `sin(π x)/(π x) = ∏ (1 - x²/k²)` &nbsp;&nbsp; and &nbsp;&nbsp; `sinh(π x)/(π x) = ∏ (1 + x²/k²)`,

each an infinite product over `k = 1, 2, 3, ...` .
Evaluating the sine product at `x = 1/2` gives the
[**Wallis product** (John Wallis, 1656)](https://en.wikipedia.org/wiki/Wallis_product),
the first infinite product ever discovered for `π`; at `x = 1/4` it gives an
algebraic multiple of `1/π`:

```mathematica
In[1]:= Product[1 - 1/(4*k^2), {k, 1, Infinity}]
Out[1]= 2/Pi

In[2]:= Product[1 - 1/(16*k^2), {k, 1, Infinity}]
Out[2]= (2 Sqrt[2])/Pi
```

Flipping the sign selects the *hyperbolic* sine product, and restricting to odd
integers selects the hyperbolic *cosine*:

```mathematica
In[1]:= Product[1 + 1/k^2, {k, 1, Infinity}]
Out[1]= Sinh[Pi]/Pi

In[2]:= Product[1 + 1/(2*k - 1)^2, {k, 1, Infinity}]
Out[2]= Cosh[1/2 Pi]
```

A different trigonometric mechanism drives
[**Viète's formula** (1593)](https://en.wikipedia.org/wiki/Vi%C3%A8te%27s_formula) —
the oldest infinite product in all of mathematics. It comes from iterating the
half-angle identity `sin θ = 2 sin(θ/2) cos(θ/2)`, so that a product of nested
cosines telescopes against a single sine. Mathilda's ``Product`Viete`` module
recognises exactly this cosine double-angle structure:

```mathematica
In[1]:= Product[Cos[Pi/2^(k + 1)], {k, 1, Infinity}]
Out[1]= 2/Pi
```

### A frontier case: the quartic product

Push the telescoping idea to a *quartic* denominator and you reach the current
edge of Mathilda's closed-form engine. The product `(k⁴-1)/(k⁴+1)` mixes real
roots (`±1`) with the complex roots of `k⁴+1`, and Mathilda returns it
unevaluated:

```mathematica
In[1]:= Product[(k^4 - 1)/(k^4 + 1), {k, 2, Infinity}]
Out[1]= Product[(k^4 - 1)/(k^4 + 1), {k, 2, Infinity}]
```

It still has a beautiful closed form, obtained by applying the sine and sinh
factorizations simultaneously (`k⁴-1` uses both `sin` and `sinh`, `k⁴+1` uses
the fourth roots of `-1`). The result is `π sinh(π)/(cosh(√2 π) - cos(√2 π))` —
and here is a genuine caution for the reader: several online sources quote this
as `cosh(π) - cos(π)`, **dropping the `√2`**. `NProduct` settles the matter to
eighteen digits:

```mathematica
In[1]:= NProduct[(k^4 - 1)/(k^4 + 1), {k, 2, Infinity}, WorkingPrecision -> 20]
Out[1]= 0.848054049352900392127

In[2]:= N[Pi*Sinh[Pi]/(Cosh[Sqrt[2]*Pi] - Cos[Sqrt[2]*Pi]), 20]
Out[2]= 0.848054049352900392134
```

The `√2` version matches; the `cosh(π) - cos(π)` version would give `≈ 2.881`, a
different number entirely. When a symbolic engine declines to guess, `NProduct`
is how you check which "known" answer is actually correct.

## Exponential limits: e, the Euler–Mascheroni constant, and Glaisher

A third family of products has *exponential* factors and evaluates to the great
analytic constants. The cleanest is the exponential of the alternating harmonic
series (which sums to `log 2`), so the product is `e^(log(1/2)) = 1/2`:

```mathematica
In[1]:= Product[Exp[(-1)^k/k], {k, 1, Infinity}]
Out[1]= 1/2
```

The next combines a telescoping factor `(1 - 1/k)` with `e^(1/k)`. The rational
part telescopes to `1/N`, while `Σ 1/k = H_N` grows like `log N + γ`; the two
combine so that the `log N` cancels and the [**Euler–Mascheroni
constant**](https://en.wikipedia.org/wiki/Euler%27s_constant) `γ` survives:

```mathematica
In[1]:= Product[E^(1/k)*(1 - 1/k), {k, 2, Infinity}]
Out[1]= E^(-1 + EulerGamma)
```

!!! warning "Mind the sign of γ"
    The partial product is `(1/N)·e^(H_N - 1)`, and since `e^(H_N) ≈ e^γ·N`, the
    limit is `e^(γ - 1) ≈ 0.6552` — which is what Mathilda returns
    (`E^(-1 + EulerGamma)`). It is a common slip to write this as `e^(1 - γ)`;
    that would be its reciprocal, `≈ 1.526`. The telescoping fixes the sign.

The [**Somos quadratic recurrence constant**](https://en.wikipedia.org/wiki/Somos%27_quadratic_recurrence_constant)
underlies the rapidly converging `2^(k/2^k)`, whose exponents `Σ k/2^k = 2`
sum to give a clean integer:

```mathematica
In[1]:= Product[2^(k/2^k), {k, 1, Infinity}]
Out[1]= 4
```

Deeper still, `k^(1/k²)` links directly to the derivative of the zeta function.
Taking logarithms gives `Σ log(k)/k² = -ζ'(2)`, so the product is `e^(-ζ'(2))`.
Mathilda's ``Product`LogSum`` module reports the equivalent
[**Glaisher–Kinkelin**](https://en.wikipedia.org/wiki/Glaisher%E2%80%93Kinkelin_constant)
form, using the identity `ζ'(2) = (π²/6)(γ + log 2π - 12 log A)`:

```mathematica
In[1]:= Product[k^(1/k^2), {k, 1, Infinity}]
Out[1]= E^(1/6 Pi^2 (-EulerGamma + 12 Log[Glaisher] - Log[2 Pi]))
```

### A second frontier case: Stirling's product

Closely related to [**Stirling's approximation**](https://en.wikipedia.org/wiki/Stirling%27s_approximation)
of `n!` is the product `(1 + 1/k)^(k + 1/2)/e`. Like the quartic product above,
Mathilda leaves it symbolically unevaluated — but its true value is a small,
memorable rearrangement of Stirling's constant `√(2π)`:

```mathematica
In[1]:= NProduct[(1 + 1/k)^(k + 1/2)/E, {k, 1, Infinity}, WorkingPrecision -> 20]
Out[1]= 1.08443755141922754661

In[2]:= N[E/Sqrt[2*Pi], 20]
Out[2]= 1.08443755141922754661
```

The value is `e/√(2π) ≈ 1.0844`. Here too the naive reciprocal `√(2π)/e ≈ 0.922`
is the version that most often appears in problem sets — and, again, it is wrong.
A one-line Stirling estimate of the partial product `e^(-N)·∏((k+1)/k)^(k+1/2)`
pins down the correct orientation.

## Euler prime products

The most consequential infinite products in mathematics run not over the
integers but over the **primes**. Euler's 1737 discovery that

> `ζ(s) = ∏ 1/(1 - p⁻ˢ)`, the product running over all primes `p`,

is the analytic gateway to the whole theory of the distribution of primes.
Mathilda's ``Product`EulerPrime`` module recognises the pattern `1/(1 - Prime[k]^-s)`
and returns the corresponding zeta value — so the Basel constant reappears, this
time as a product over primes:

```mathematica
In[1]:= Product[1/(1 - 1/Prime[k]^2), {k, 1, Infinity}]
Out[1]= 1/6 Pi^2

In[2]:= Product[1/(1 - 1/Prime[k]^4), {k, 1, Infinity}]
Out[2]= 1/90 Pi^4

In[3]:= Product[1/(1 - 1/Prime[k]^6), {k, 1, Infinity}]
Out[3]= 1/945 Pi^6
```

Twisting the prime product by a **Dirichlet character** — here the sign
`(-1)^((p-1)/2)`, which is `+1` for primes `p ≡ 1 (mod 4)` and `-1` for
`p ≡ 3 (mod 4)` — gives the Euler product for the
[Dirichlet beta function](https://en.wikipedia.org/wiki/Dirichlet_beta_function),
the prime-side mirror of the Leibniz series for `π/4`:

```mathematica
In[1]:= Product[1/(1 - (-1)^((Prime[k] - 1)/2)/Prime[k]), {k, 2, Infinity}]
Out[1]= 1/4 Pi
```

## Geometric and Fermat-number products

A final family exploits repeated **difference of squares**. In the product
`(1 + x^(2^k))`, each factor is exactly what is needed to complete
`(1 - x)·(1 + x)(1 + x²)(1 + x⁴)... = 1 - x^(2^(N+1))`, so the whole product
telescopes to `1/(1 - x)`. Mathilda's ``Product`Cantor`` module handles these
double-exponential products; with `x = 1/3` it gives `3/2`:

```mathematica
In[1]:= Product[1 + (1/3)^(2^k), {k, 0, Infinity}]
Out[1]= 3/2
```

This is the product form of the
[binary expansion / Fermat-number identity](https://en.wikipedia.org/wiki/Fermat_number#Other_theorems_about_Fermat_numbers),
and it converges *doubly exponentially* — three factors already give five correct
digits.

## Where to next

You have now toured Mathilda's full product cascade: rational telescoping, the
Gamma-reflection route for complex roots, the sine/sinh/cosine Weierstrass
factorizations, Viète's nested cosines, the exponential products for `e`, `γ`,
and Glaisher's constant, the Euler prime products for `ζ` and Dirichlet `β`, and
the double-exponential Cantor/Fermat products — plus two honest frontier cases
where `NProduct` is the arbiter of a disputed closed form.

- The companion [symbolic summation tutorial](11-symbolic-summation.md) covers
  the parallel machinery for `Sum` — Gosper telescoping, the zeta family, Euler
  sums, and the hypergeometric `π`-machines.
- The [numerical calculus tutorial](09-numerical-calculus.md) develops `NProduct`
  (and `NSum`) in full, for products with no closed form.
- The [special functions tutorial](10-special-functions.md) covers `Gamma`,
  `Zeta`, and the constants (`EulerGamma`, `Glaisher`, `Catalan`) these products
  evaluate to.

As always, `?Product` (or ``?Product`EulerPrime``) at the prompt shows the built-in
help without leaving the REPL.
