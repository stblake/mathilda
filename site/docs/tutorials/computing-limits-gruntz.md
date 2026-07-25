# Computing limits with Gruntz's algorithm

`Limit[expr, x -> a]` in Mathilda is not one procedure but a *cascade* of them:
rational-function shortcuts, series expansion, l'Hôpital, symmetry tricks, and —
as the last, most powerful layer — a faithful implementation of **Dominik
Gruntz's most-rapidly-varying (mrv) algorithm** (1996 ETH thesis, *On Computing
Limits in a Symbolic Manipulation System*). This tutorial is a deep dive into
that engine: what makes limits hard, the idea that makes Gruntz's method work,
and a graded battery of worked examples that push it — the same families that
form Mathilda's [`gruntz_stress_tests`](https://github.com/stblake/mathilda)
suite.

You can force the Gruntz engine directly with `Method -> "Gruntz"`, which is what
we do throughout so you always know exactly which machinery answered:

```mathematica
In[1]:= Limit[(3^x + 5^x)^(1/x), x -> Infinity, Method -> "Gruntz"]
Out[1]= 5
```

Every transcript below was produced by the actual Mathilda binary. Mathilda's
output form is sometimes arranged differently from a textbook — `E^x` instead of
`Exp[x]`, `1/E^2` instead of `E^(-2)` — but it is always mathematically correct.

---

## 1. Why limits are hard

Substituting the limit point usually gives an indeterminate form (`0/0`,
`∞ - ∞`, `1^∞`, `∞/∞`). The classical repair tools each have a failure mode that
a computer algebra system runs into constantly.

**l'Hôpital may never terminate.** Differentiating numerator and denominator can
cycle forever or make the expression *grow*. **Power series need the right
scale.** A bottom-up series expansion has to carry exact remainder terms so that
mutually cancelling pieces do not silently drop precision — the *cancellation
problem*, or "intermediate expression swell". Consider two nested exponentials
that agree to many orders:

```mathematica
In[1]:= Limit[E^x (E^(1/x - E^-x) - E^(1/x)), x -> Infinity, Method -> "Gruntz"]
Out[1]= -1
```

Naively, `E^(1/x - E^-x)` and `E^(1/x)` both tend to `1`, so the bracket is
`0`, multiplied by `E^x -> ∞`: an `∞ · 0`. The two exponentials cancel to leading
order, and the *answer* lives entirely in the term that survives the
cancellation. Gruntz's insight is to expand the **whole** function as a series in
its single most-rapidly-varying subexpression, so the surviving term is produced
directly — never reconstructed from a difference of large quantities.

---

## 2. The idea in one paragraph

Given `f(x)` as `x -> ∞`, the algorithm:

1. finds the **most rapidly varying (mrv) set** — the subexpressions that grow or
   decay fastest and share a *comparability class* (`f`, `g` are comparable when
   `ln f / ln g` tends to a nonzero finite constant);
2. picks one representative `ω` from that set with `ω -> 0⁺`, and **rewrites**
   every mrv element as a power `A · ω^c`;
3. computes the **leading term** `c₀ · ω^{e₀}` of the series of `f` in `ω`;
4. reads the limit off the leading exponent `e₀`: negative `⇒ ±∞`, positive
   `⇒ 0`, zero `⇒` recurse on the coefficient `c₀`.

The scale `ω` is discovered *on the fly* — there is never a fixed asymptotic
basis to overflow. Everything below is an application of these four steps; the
families are ordered so each one stresses a different part of the machinery.

---

## 3. Extracting the dominant base — `(Σ aᵢ^x)^{1/x}`

*(Thesis example 8.12.)* The cleanest illustration of "find the dominant scale".
For positive bases, `(Σ cᵢ aᵢ^x)^{1/x} -> max aᵢ`, because the largest base
governs the mrv class and the coefficients `cᵢ` are invisible to it.

```mathematica
In[1]:= Limit[(2^x + 3^x)^(1/x), x -> Infinity, Method -> "Gruntz"]
Out[1]= 3

In[2]:= Limit[(2^x + 3^x + 5^x + 7^x + 11^x)^(1/x), x -> Infinity, Method -> "Gruntz"]
Out[2]= 11

In[3]:= Limit[(3 2^x + 5^x)^(1/x), x -> Infinity, Method -> "Gruntz"]
Out[3]= 5
```

The engine does not care whether the dominant base is written as an integer, a
power, or a transcendental constant — only its magnitude matters. Here `e < 3`,
so `3` wins; and `e² ≈ 7.39 > 5`, so `E^2` wins:

```mathematica
In[1]:= Limit[(E^x + 3^x)^(1/x), x -> Infinity, Method -> "Gruntz"]
Out[1]= 3

In[2]:= Limit[(E^(2 x) + 5^x)^(1/x), x -> Infinity, Method -> "Gruntz"]
Out[2]= E^2

In[3]:= Limit[(Pi^x + 3^x)^(1/x), x -> Infinity, Method -> "Gruntz"]
Out[3]= Pi
```

An outer power just rides along — `(2^x + 3^x)^{p/x} -> 3^p`:

```mathematica
In[1]:= Limit[(2^x + 3^x)^(3/x), x -> Infinity, Method -> "Gruntz"]
Out[1]= 27
```

---

## 4. The cancellation problem — `E^{ax}(E^{g - c E^{-ax}} - E^{g})`

*(Thesis example 8.1.)* This is the family that defeats bottom-up series. The two
exponentials agree to leading order; the limit is exactly the coefficient of the
surviving term, `-c`, **independent** of the finite part `g` and of the decay
rate `a`:

\[
\lim_{x\to\infty} e^{ax}\!\left(e^{\,g(x) - c\,e^{-ax}} - e^{\,g(x)}\right) = -c.
\]

Sweep the coefficient `c` (including sign and fractions):

```mathematica
In[1]:= Limit[E^x (E^(1/x - 2 E^-x) - E^(1/x)), x -> Infinity, Method -> "Gruntz"]
Out[1]= -2

In[2]:= Limit[E^x (E^(1/x - E^-x/2) - E^(1/x)), x -> Infinity, Method -> "Gruntz"]
Out[2]= -1/2

In[3]:= Limit[E^x (E^(1/x + E^-x) - E^(1/x)), x -> Infinity, Method -> "Gruntz"]
Out[3]= 1
```

Now change the decay rate `a` — the balancing `E^{ax}` prefactor tracks it, and
the answer is still `-c`:

```mathematica
In[1]:= Limit[E^(2 x) (E^(1/x - 2 E^(-2 x)) - E^(1/x)), x -> Infinity, Method -> "Gruntz"]
Out[1]= -2

In[2]:= Limit[E^(3 x) (E^(1/x - 5 E^(-3 x)) - E^(1/x)), x -> Infinity, Method -> "Gruntz"]
Out[2]= -5
```

---

## 5. Nested exponential towers — `E^{E^{x + a E^{-x}}} / E^{E^x}`

*(Thesis examples 8.5–8.8.)* Climb one level. Since `E^{x + a E^{-x}} ~ E^x + a`,
the ratio tends to `E^a` — a whole real parameter recovered from deep inside a
double exponential:

```mathematica
In[1]:= Limit[E^(E^(x + 2 E^-x))/E^(E^x), x -> Infinity, Method -> "Gruntz"]
Out[1]= E^2

In[2]:= Limit[E^(E^(x - E^-x))/E^(E^x), x -> Infinity, Method -> "Gruntz"]
Out[2]= 1/E

In[3]:= Limit[E^(E^(x + E^-x/3))/E^(E^x), x -> Infinity, Method -> "Gruntz"]
Out[3]= E^(1/3)
```

If the perturbation decays *faster* than `E^{-x}` the shift vanishes and the
ratio is `1` — the engine correctly detects that the two towers share a leading
term to all visible orders:

```mathematica
In[1]:= Limit[E^(E^(x + E^(-x^2)))/E^(E^x), x -> Infinity, Method -> "Gruntz"]
Out[1]= 1
```

!!! note "An honest limit of the reach"
    Push the same balance one level deeper — `E^{E^{E^{x + E^{-x}}}}/E^{E^{E^x}}`
    — and the required cancellation exceeds what the dense-series machinery
    reaches here. Rather than risk a wrong value, Gruntz's engine returns the
    input **unevaluated**. Knowing when to abstain is part of a trustworthy
    limit engine; §11 collects these.

---

## 6. Smooth functions at a vanishing argument

*(Thesis examples 8.21/8.22.)* For a differentiable `F`,
`E^{ax}(F[b/x + c E^{-ax}] - F[b/x]) -> c · F'(0⁺)`. The mrv machinery reduces
each `F` to its series and the derivative falls out. For `Sin`, `Tan`, `Sinh`,
`Tanh`, `ArcTan` the slope at `0` is `1`, so the limit is `c`:

```mathematica
In[1]:= Limit[E^x (Sin[1/x + 2 E^-x] - Sin[1/x]), x -> Infinity, Method -> "Gruntz"]
Out[1]= 2

In[2]:= Limit[E^x (Tan[1/x + E^-x] - Tan[1/x]), x -> Infinity, Method -> "Gruntz"]
Out[2]= 1

In[3]:= Limit[E^x (ArcTan[1/x + E^-x] - ArcTan[1/x]), x -> Infinity, Method -> "Gruntz"]
Out[3]= 1
```

`Cos` is the instructive exception: its slope is `-sin(1/x) -> 0`, so the whole
thing collapses to `0`:

```mathematica
In[1]:= Limit[E^x (Cos[1/x + E^-x] - Cos[1/x]), x -> Infinity, Method -> "Gruntz"]
Out[1]= 0
```

And the construction survives being wrapped in another exponential tower (8.22):

```mathematica
In[1]:= Limit[E^(E^x) (E^Sin[1/x + 2 E^(-E^x)] - E^Sin[1/x]), x -> Infinity, Method -> "Gruntz"]
Out[1]= 2
```

---

## 7. Nested logarithms

*(Thesis examples 8.19/8.20.)* Logs live at the *bottom* of the asymptotic scale,
and ratios of them measure polynomial degree. `Log[P(x)]/Log[x] -> deg P`, and a
further log flattens the leading coefficient away:

```mathematica
In[1]:= Limit[Log[x^3 + x]/Log[x], x -> Infinity, Method -> "Gruntz"]
Out[1]= 3

In[2]:= Limit[Log[Log[x^7]]/Log[Log[x]], x -> Infinity, Method -> "Gruntz"]
Out[2]= 1

In[3]:= Limit[E^(Log[x^2 + x]/Log[x]), x -> Infinity, Method -> "Gruntz"]
Out[3]= E^2
```

The subtle members require a *leading-order cancellation of two nearly-equal
logs* — the exact difficulty of thesis 8.19:

```mathematica
In[1]:= Limit[Log[x + Log[x]] - Log[x], x -> Infinity, Method -> "Gruntz"]
Out[1]= 0

In[2]:= Limit[x (Log[x + Log[x]] - Log[x])/Log[x], x -> Infinity, Method -> "Gruntz"]
Out[2]= 1
```

---

## 8. The Hardy scale — sub-polynomial growth

*(Thesis example 8.9.)* Anything of the form `e^{o(log x)}` — a power of a
logarithm, or an exponential of a *fractional* power of a logarithm — grows more
slowly than **every** positive power of `x`. So dividing by `x^p` sends it to
`0`, and the machinery places it correctly on the scale between the constants and
the powers:

```mathematica
In[1]:= Limit[E^(Sqrt[Log[x]] Log[Log[x]]^2)/Sqrt[x], x -> Infinity, Method -> "Gruntz"]
Out[1]= 0

In[2]:= Limit[(Log[x])^100/x, x -> Infinity, Method -> "Gruntz"]
Out[2]= 0

In[3]:= Limit[E^(Log[x]^(2/3))/x, x -> Infinity, Method -> "Gruntz"]
Out[3]= 0
```

Turn any of these upside down and the same reasoning gives `Infinity`:

```mathematica
In[1]:= Limit[x/E^(Sqrt[Log[x]]), x -> Infinity, Method -> "Gruntz"]
Out[1]= Infinity
```

---

## 9. Everyday limits the same engine also nails

The mrv method is not only for exotic towers. Set on classic textbook limits it
handles them without special-casing.

**Conjugate radical differences** *(generalizing thesis example 2.5)*:

```mathematica
In[1]:= Limit[Sqrt[x^2 + 3 x] - x, x -> Infinity, Method -> "Gruntz"]
Out[1]= 3/2

In[2]:= Limit[(x^3 + x^2)^(1/3) - x, x -> Infinity, Method -> "Gruntz"]
Out[2]= 1/3

In[3]:= Limit[Sqrt[x] (Sqrt[x + 1] - Sqrt[x]), x -> Infinity, Method -> "Gruntz"]
Out[3]= 1/2
```

**Power-series limits at a finite point** *(thesis examples 2.6/2.7)* — note the
limit point is `x -> 0` and `x -> 1`:

```mathematica
In[1]:= Limit[((1 + 2 x)^3 - 1)/x, x -> 0, Method -> "Gruntz"]
Out[1]= 6

In[2]:= Limit[(Sqrt[x] - 1)/(x^(1/3) - 1), x -> 1, Method -> "Gruntz"]
Out[2]= 3/2

In[3]:= Limit[(x^4 - 1)/(x^7 - 1), x -> 1, Method -> "Gruntz"]
Out[3]= 4/7
```

**Max/Min inside a limit** *(thesis example 8.37)* resolve to the dominant
argument before the limit is taken:

```mathematica
In[1]:= Limit[x Max[1/x, 2/x, 3/x], x -> Infinity, Method -> "Gruntz"]
Out[1]= 3

In[2]:= Limit[Log[Max[E^x, E^(2 x)]]/x, x -> Infinity, Method -> "Gruntz"]
Out[2]= 2
```

---

## 10. Essential singularities — special functions

*(Thesis examples 8.23–8.34.)* Gruntz's *Phase 2* extension isolates a
semi-tractable special function `F[g]` (with `g -> ±∞`) and replaces it by its
asymptotic series, turning the essential singularity into explicit
exponentials and logs the core engine can expand.

The Gaussian tail shift `Erfc[x + c/x]/Erfc[x] -> e^{-2c}` is a razor-sharp test
of that expansion:

```mathematica
In[1]:= Limit[Erfc[x + 1/x]/Erfc[x], x -> Infinity, Method -> "Gruntz"]
Out[1]= 1/E^2

In[2]:= Limit[Erfc[x + 3/x]/Erfc[x], x -> Infinity, Method -> "Gruntz"]
Out[2]= 1/E^6

In[3]:= Limit[x Erfc[x] E^(x^2), x -> Infinity, Method -> "Gruntz"]
Out[3]= 1/Sqrt[Pi]
```

The exponential integral obeys `Ei[z] ~ e^z / z`:

```mathematica
In[1]:= Limit[x ExpIntegralEi[x] E^-x, x -> Infinity, Method -> "Gruntz"]
Out[1]= 1

In[2]:= Limit[ExpIntegralEi[2 x] E^(-2 x) x, x -> Infinity, Method -> "Gruntz"]
Out[2]= 1/2
```

The Riemann zeta tail `ζ(x) - 1 ~ 2^{-x}` — the exponential scale, not `1/x` —
gives a family of clean Dirichlet ratios `(ζ(x)-1)/(ζ(x+k)-1) -> 2^k`:

```mathematica
In[1]:= Limit[(Zeta[x] - 1) 2^x, x -> Infinity, Method -> "Gruntz"]
Out[1]= 1

In[2]:= Limit[(Zeta[x] - 1)/(Zeta[x + 3] - 1), x -> Infinity, Method -> "Gruntz"]
Out[2]= 8

In[3]:= Limit[Log[Zeta[x] - 1]/x, x -> Infinity, Method -> "Gruntz"]
Out[3]= -Log[2]
```

Digamma and log-Gamma growth laws — `ψ(x) ~ Log[x] - 1/(2x)`,
`Log[Γ(x)] ~ x Log[x]`:

```mathematica
In[1]:= Limit[x (PolyGamma[0, x] - Log[x]), x -> Infinity, Method -> "Gruntz"]
Out[1]= -1/2

In[2]:= Limit[Log[Gamma[x]]/(x Log[x]), x -> Infinity, Method -> "Gruntz"]
Out[2]= 1

In[3]:= Limit[Exp[PolyGamma[0, x]]/x, x -> Infinity, Method -> "Gruntz"]
Out[3]= 1
```

---

## 11. Knowing when to abstain

A limit engine that occasionally returns a *wrong* finite value is worse than
useless. When a case needs a cancellation or an asymptotic form beyond the
current machinery, Mathilda's Gruntz engine returns the input **unevaluated**
rather than guess. These are honest gaps, not errors:

```mathematica
In[1]:= Limit[E^(E^(E^(x + E^-x)))/E^(E^(E^x)), x -> Infinity, Method -> "Gruntz"]
Out[1]= Limit[E^(-E^E^x + E^E^(x + E^(-x))), x -> Infinity, Method -> "Gruntz"]

In[2]:= Limit[(Gamma[x + 1/Gamma[x]] - Gamma[x])/Log[x], x -> Infinity, Method -> "Gruntz"]
Out[2]= Limit[(-Gamma[x] + Gamma[x + 1/Gamma[x]])/Log[x], x -> Infinity, Method -> "Gruntz"]
```

The first is the three-level tower of §5; the second is the Gamma-*difference*
asymptotic of thesis 8.28, which needs the Stirling series expanded and
cancelled to an order the dense-polynomial series does not reach. Both are
correct *refusals*.

---

## Where this lives

- **Engine:** `src/calculus/gruntz.c` (a C port of Gruntz's mrv algorithm,
  structured after SymPy's `gruntz.py`), reached from `src/calculus/limit.c`.
- **Tests:** the graded families above are the
  [`gruntz_stress_tests`](https://github.com/stblake/mathilda) suite
  (`tests/test_gruntz_stress.c`); the exact thesis expressions are in
  `tests/test_gruntz.c`.
- **Further reading:** Dominik Gruntz, *On Computing Limits in a Symbolic
  Manipulation System*, Diss. ETH No. 11432 (1996).

The Gruntz layer also runs automatically inside plain `Limit[...]` — the
`Method -> "Gruntz"` option in this tutorial simply pins it so you can see
exactly which engine produced each result.
