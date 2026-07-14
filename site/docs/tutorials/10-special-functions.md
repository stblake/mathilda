# Special functions

Beyond the elementary functions (`Sin`, `Exp`, `Log`, …) lies a second tier of
**special functions** — the higher transcendental functions of analysis and
number theory. Mathilda implements a substantial library of them: the gamma and
beta functions, the digamma/polygamma family, the Riemann zeta function, the
error functions, the exponential and logarithmic integrals, the polylogarithm,
the Bernoulli and Euler numbers, and the hypergeometric family that unifies many
of the others.

Two themes run through all of them. First, on **special arguments** they reduce
to exact closed forms — `Gamma[1/2]` is `Sqrt[Pi]`, `Zeta[2]` is `π²/6` — so they
participate in symbolic algebra and calculus. Second, on **general numeric
arguments** they evaluate to machine-precision or, via `N[expr, prec]`,
arbitrary-precision (MPFR) numbers, including complex ones.

Every transcript below was produced by the actual Mathilda binary.

## The gamma function

`Gamma[z]` is Euler's gamma function `Γ(z)`, the continuous extension of the
factorial: `Γ(n) = (n-1)!`. Integer and half-integer arguments reduce exactly,
and everything else evaluates numerically:

```mathematica
In[1]:= Gamma[5]
Out[1]= 24

In[2]:= Gamma[1/2]
Out[2]= Sqrt[Pi]

In[3]:= N[Gamma[1/2]]
Out[3]= 1.77245
```

`In[1]` is `4! = 24`. `In[2]` is the classic value `Γ(1/2) = √π`; half-integer
arguments always come out as rational multiples of `Sqrt[Pi]`. `In[3]` is its
decimal value. `Gamma` also has a two-argument *incomplete* form `Gamma[a, z]`,
and the closely-related **beta function** `Beta[a, b] = Γ(a)Γ(b)/Γ(a+b)` reduces
to exact rationals when its arguments are integers:

```mathematica
In[1]:= Beta[2, 3]
Out[1]= 1/12

In[2]:= Beta[5, 3]
Out[2]= 1/105
```

### Log-gamma

`Γ` overflows quickly — `Γ(100)` is a 156-digit number. When you only need its
logarithm (for factorial ratios, asymptotics, or statistics) reach for
`LogGamma[z] = log Γ(z)`, which stays finite where `Gamma` would blow up:

```mathematica
In[1]:= N[LogGamma[100]]
Out[1]= 359.134
```

`LogGamma` is the analytic continuation of `log(Γ(z))` with a single branch cut
on the negative reals — subtly different from `Log[Gamma[z]]`, which would wrap.
Its derivative is the digamma function, which we meet next.

## Digamma and polygamma

`PolyGamma[n, z]` is the `n`-th derivative of `log Γ(z)`. The zeroth,
`PolyGamma[0, z]`, is the **digamma** function `ψ(z) = Γ′(z)/Γ(z)`; higher `n`
give the **polygamma** functions. At positive integers they reduce to exact
combinations of `EulerGamma` and powers of `Pi`:

```mathematica
In[1]:= PolyGamma[0, 1]
Out[1]= -EulerGamma

In[2]:= PolyGamma[1, 1]
Out[2]= 1/6 Pi^2

In[3]:= N[PolyGamma[0, 3]]
Out[3]= 0.922784
```

`In[1]` is `ψ(1) = -γ` (the negative of Euler's constant); `In[2]` is the
trigamma value `ψ′(1) = π²/6`. `In[3]` shows the numerical value at a general
point.

## Pochhammer symbol

`Pochhammer[a, n]` is the rising factorial `(a)ₙ = a(a+1)⋯(a+n-1) =
Γ(a+n)/Γ(a)`. For integer `n` it expands to a product of `n` linear factors —
exact for numeric `a`, polynomial for symbolic `a`:

```mathematica
In[1]:= Pochhammer[3, 4]
Out[1]= 360

In[2]:= Pochhammer[x, 3]
Out[2]= x (1 + x) (2 + x)
```

`In[1]` is `3·4·5·6 = 360`. The Pochhammer symbol is the building block of the
hypergeometric series at the end of this tutorial.

## The Riemann zeta function

`Zeta[s]` is the Riemann zeta function `ζ(s) = Σ_{k≥1} k^-s`, the central object
of analytic number theory. At even positive integers it is a rational multiple of
a power of `π`; at negative integers it is rational:

```mathematica
In[1]:= Zeta[2]
Out[1]= 1/6 Pi^2

In[2]:= Zeta[4]
Out[2]= 1/90 Pi^4

In[3]:= Zeta[-1]
Out[3]= -1/12

In[4]:= N[Zeta[3]]
Out[4]= 1.20206
```

`In[1]` is the Basel value `ζ(2) = π²/6`. `In[3]`, `ζ(-1) = -1/12`, is the value
behind the (in)famous "1 + 2 + 3 + ⋯ = -1/12" regularisation. `In[4]` is Apéry's
constant `ζ(3) ≈ 1.20206`, which has no known elementary closed form, so it stays
symbolic until you ask for `N`. The two-argument form `Zeta[s, a]` is the Hurwitz
zeta function.

The Laurent coefficients of `ζ` about `s = 1` are the **Stieltjes constants**,
available as `StieltjesGamma[n]`. The zeroth is Euler's constant:

```mathematica
In[1]:= StieltjesGamma[0]
Out[1]= EulerGamma
```

## Bernoulli and Euler numbers

`BernoulliB[n]` gives the `n`-th Bernoulli number — exact rationals that appear in
the Euler–Maclaurin formula, sums of powers, and the values of `Zeta` at even
integers. `EulerE[n]` gives the integer Euler numbers. Both also have a
polynomial form, `BernoulliB[n, x]` and `EulerE[n, x]`:

```mathematica
In[1]:= BernoulliB[10]
Out[1]= 5/66

In[2]:= BernoulliB[4, x]
Out[2]= -1/30 + x^2 - 2 x^3 + x^4

In[3]:= EulerE[4]
Out[3]= 5

In[4]:= EulerE[2, x]
Out[4]= -x + x^2
```

`In[1]` is `B₁₀ = 5/66`; the odd-index Bernoulli numbers past the first are all
zero. `In[2]` is the Bernoulli polynomial `B₄(x)`, whose coefficients are exact
rationals.

## The polylogarithm

`PolyLog[n, z]` is the polylogarithm `Liₙ(z) = Σ_{k≥1} zᵏ/kⁿ`. It generalises
the ordinary logarithm and connects back to `Zeta`:

```mathematica
In[1]:= PolyLog[1, z]
Out[1]= -Log[1 - z]

In[2]:= PolyLog[2, 1]
Out[2]= 1/6 Pi^2

In[3]:= N[PolyLog[2, 1/2]]
Out[3]= 0.582241
```

`In[1]` shows the base case `Li₁(z) = -log(1-z)`. `In[2]` is `Li₂(1) = ζ(2) =
π²/6` — at `z = 1` the polylogarithm reduces to zeta. `In[3]` is the dilogarithm
at `1/2`, a general numeric value.

## The error function family

`Erf[z] = (2/√π) ∫₀^z e^(-t²) dt` is the error function, ubiquitous in
probability and diffusion. Its relatives are the complementary error function
`Erfc[z] = 1 - Erf[z]` and the imaginary error function `Erfi[z] = -I·Erf[I z]`:

```mathematica
In[1]:= Erf[0]
Out[1]= 0

In[2]:= Erf[Infinity]
Out[2]= 1

In[3]:= N[Erf[1]]
Out[3]= 0.842701

In[4]:= N[Erfc[1]]
Out[4]= 0.157299
```

`In[1]` and `In[2]` are the exact endpoint values; note `Erf[1] + Erfc[1] = 1`,
which `In[3]` and `In[4]` confirm numerically (`0.842701 + 0.157299 = 1`). These
are entire functions — no branch cuts — and evaluate for complex arguments too.

The error function is invertible: `InverseErf[s]` solves `Erf[z] == s` for `z`
(real `s` in `[-1, 1]`), and `InverseErfc` does the same for `Erfc`:

```mathematica
In[1]:= InverseErf[0]
Out[1]= 0

In[2]:= N[InverseErf[1/2]]
Out[2]= 0.476936
```

## Exponential and logarithmic integrals

`ExpIntegralEi[z]` is the exponential integral `Ei(z)`, and `LogIntegral[z]` is
the logarithmic integral `li(z) = Ei(log z)` — the latter is the leading
approximation to the prime-counting function in the prime number theorem. Both
are non-elementary and evaluate numerically:

```mathematica
In[1]:= N[ExpIntegralEi[1]]
Out[1]= 1.89512

In[2]:= N[LogIntegral[2]]
Out[2]= 1.04516
```

## Sine and cosine integrals

`SinIntegral[z]` is the sine integral `Si(z) = Integral_0^z Sin[t]/t dt`, and
`CosIntegral[z]` is the cosine integral `Ci(z) = -Integral_z^Infinity Cos[t]/t dt`.
Both are non-elementary and evaluate numerically at machine or arbitrary precision:

```mathematica
In[1]:= N[SinIntegral[2]]
Out[1]= 1.60541

In[2]:= N[CosIntegral[2]]
Out[2]= 0.422981
```

`Si` is entire and odd, with `SinIntegral[±Infinity] = ±Pi/2`. `Ci` is different:
it has a logarithmic singularity at the origin and a branch cut along the negative
real axis, so `CosIntegral[0] = -Infinity` and a negative real argument returns a
complex value:

```mathematica
In[1]:= SinIntegral[Infinity]
Out[1]= 1/2 Pi

In[2]:= CosIntegral[0]
Out[2]= -Infinity

In[3]:= CosIntegral[-2.]
Out[3]= 0.422981 + 3.14159 I
```

The two are linked through their derivatives. Differentiating the sine integral
produces the cardinal sine `Sinc[z] = Sin[z]/z` (with the removable singularity
filled in as `Sinc[0] = 1`), while the cosine integral differentiates to `Cos[z]/z`:

```mathematica
In[1]:= D[SinIntegral[x], x]
Out[1]= Sinc[x]

In[2]:= D[CosIntegral[x], x]
Out[2]= Cos[x]/x

In[3]:= Sinc[0]
Out[3]= 1
```

## The hypergeometric family

Many of the functions above are special cases of one master function: the
generalized hypergeometric `HypergeometricPFQ[{a₁,…}, {b₁,…}, z]`, whose series
is built from Pochhammer symbols. The common low-order cases have their own
names — `Hypergeometric0F1`, `Hypergeometric1F1` (Kummer's confluent function),
and `Hypergeometric2F1` (Gauss's function). On nice parameters they collapse to
elementary functions:

```mathematica
In[1]:= Hypergeometric2F1[1, 1, 2, z]
Out[1]= -Log[1 - z]/z

In[2]:= N[Hypergeometric2F1[1, 1, 2, 1/2]]
Out[2]= 1.38629

In[3]:= N[Hypergeometric1F1[1, 2, 1]]
Out[3]= 1.71828
```

`In[1]` shows `₂F₁(1, 1; 2; z) = -log(1-z)/z` in closed form. `In[2]` evaluates it
at `z = 1/2`, giving `2 log 2 ≈ 1.38629`. `In[3]` is `₁F₁(1; 2; 1) = e - 1 ≈
1.71828`. Equivalently, `HypergeometricPFQ[{1, 1}, {2}, z]` is the same `₂F₁` as
`In[1]`.

## Where to next

This tutorial introduced Mathilda's special-function library: `Gamma`/`Beta` and
`LogGamma`, the `PolyGamma` family, `Pochhammer`, `Zeta` and `StieltjesGamma`,
`BernoulliB`/`EulerE`, `PolyLog`, the `Erf` family with its inverses, the
`ExpIntegralEi`/`LogIntegral` pair, and the unifying hypergeometric functions.

- For the full behaviour of each — exact reductions, derivatives, branch cuts,
  and numerical methods — see the
  [special functions](../documentation/special-functions/index.md) section of the
  function documentation.
- Many special values returned here (`π²/6`, `√π`, `EulerGamma`, `e - 1`) show up
  as the answers to the limits, sums, and integrals in the
  [numerical calculus tutorial](09-numerical-calculus.md).
- To revisit the guided path from the start, head back to the
  [tutorials index](index.md).

As always, `?Name` at the prompt — for example `?Zeta` — prints a function's
built-in help string without leaving the REPL.
