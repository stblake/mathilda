# Integration by differentiation under the integral sign

`Method -> "DiffUnderInt"` evaluates a **parameter-dependent definite integral**

\[
I(p) \;=\; \int_a^b f(x, p)\, dx
\]

by the **Leibniz integral rule** — differentiating the integrand with respect to
a parameter \(p\), integrating the (simpler) result, and integrating back. It is
the technique Leibniz stated in the 1690s, that Euler wielded to evaluate hard
integrals, and that Richard Feynman made famous in *Surely You're Joking, Mr.
Feynman!* — where he calls it "a different box of tools" and puts it to work at
Los Alamos. Its informal name, **"Feynman's trick,"** dates from that book.

```mathematica
In[1]:= Integrate[(x^a - 1)/Log[x], {x, 0, 1}, Method -> "DiffUnderInt"]
Out[1]= Log[1 + a]
```

The integrand \((x^a-1)/\log x\) has no elementary antiderivative in \(x\), yet
the answer is a clean logarithm. The trick: differentiate under the integral
sign with respect to \(a\), and the integral *collapses*.

---

## 1. The mathematical idea

### The Leibniz integral rule

For an integral with a parameter \(p\) and fixed limits, differentiation passes
through the integral sign,

\[
\frac{d}{dp}\int_a^b f(x,p)\,dx \;=\; \int_a^b \frac{\partial f}{\partial p}(x,p)\,dx,
\]

provided \(f\) and \(\partial_p f\) are continuous and the differentiated
integral converges uniformly in \(p\) (a dominated-convergence condition). The
power of the identity is that \(\partial_p f\) is very often **simpler** than
\(f\) — an exponential loses its polynomial prefactor, a logarithm or arctangent
turns into a rational function — so the integral on the right is one you can do.

### The method as a first-order ODE

Write \(J(p)=\int_a^b \partial_p f\,dx\) for that simpler integral. Then \(I(p)\)
satisfies the elementary ordinary differential equation

\[
I'(p) \;=\; J(p),
\qquad\text{so}\qquad
I(p) \;=\; \int J(p)\,dp \;+\; C,
\]

and the constant \(C\) is pinned by evaluating both sides at a **base point**
\(p_0\) where \(I(p_0)\) is known exactly — usually a value of the parameter at
which the integrand vanishes identically, giving \(I(p_0)=0\).

For the opening example, with \(p=a\) and \(f=(x^a-1)/\log x\):

\[
I'(a)=\int_0^1 \frac{\partial}{\partial a}\frac{x^a-1}{\log x}\,dx
      =\int_0^1 x^a\,dx
      =\frac{1}{a+1},
\]

because \(\partial_a x^a = x^a\log x\) cancels the \(\log x\) in the denominator.
Integrating, \(I(a)=\log(a+1)+C\), and the base point \(a=0\) makes the integrand
identically zero, so \(I(0)=0\) forces \(C=0\). Hence \(\int_0^1
(x^a-1)/\log x\,dx=\log(1+a)\) — exactly `Out[1]` above. Mathilda even shows the
first step for you; the \(a\)-derivative really is that simple:

```mathematica
In[1]:= D[(x^a - 1)/Log[x], a]
Out[1]= x^a
```

### A little history

The rule is Leibniz's, and Euler used it heavily in the 18th century. Its
20th-century fame is due to Feynman, who learned it from Woods's *Advanced
Calculus* as a schoolboy and later:

> "So I got a great reputation for doing integrals, only because my box of tools
> was different from everybody else's…" — R. P. Feynman

The uniform framing "every case is the first-order ODE \(I'(p)=J(p)\)" follows
Boulnois (2023). A careful account of *when* the rule is valid — and the
conditionally-convergent pitfalls where it is not — is Conrad's expository note.

### References

1. R. P. Feynman, *Surely You're Joking, Mr. Feynman!*, W. W. Norton, 1985 —
   "A different box of tools" (the popular account).
2. F. S. Woods, *Advanced Calculus*, Ginn & Co., 1926 — the text Feynman learned
   the rule from.
3. K. Conrad, "Differentiating under the integral sign," expository notes,
   [kconrad.math.uconn.edu](https://kconrad.math.uconn.edu/blurbs/analysis/diffunderint.pdf)
   — rigorous hypotheses and worked examples (§12: conditional convergence).
4. D. Boulnois, "A unified first-order ODE approach to differentiation under the
   integral sign," arXiv:2308.09619, 2023.
5. G. B. Folland, *Real Analysis*, 2nd ed., Wiley, 1999 — the dominated-
   convergence justification (Thm 2.27).

---

## 2. How Mathilda realises it

The method runs the four steps above literally, and every step is symbolic —
there is **no numerical quadrature** anywhere in it.

1. **Pick a parameter.** The free symbols of the integrand other than \(x\) are
   the candidate parameters. Their signs are read from `Assumptions`.

2. **Differentiate and integrate.** For a candidate \(p\), form
   \(g=\operatorname{Simplify}[\partial_p f]\) and evaluate the inner definite
   integral \(J=\int_a^b g\,dx\). If it does not close to a finite closed form,
   that parameter is abandoned and the next is tried.

3. **Integrate back over the parameter** to get the antiderivative
   \(G(p)=\int J\,dp\).

4. **Fix the constant** from a base point: \(I(p)=G(p)-G(p_0)+I(p_0)\), where
   \(p_0\) is a value at which \(f\) vanishes identically (\(I(p_0)=0\)) or the
   integral is directly computable.

5. **Verify.** The result is accepted only if
   \(\operatorname{Simplify}[\,I'(p)-J\,]\equiv 0\) under the assumptions, together
   with the exact base value. The verification is a genuine proof, not a numeric
   spot-check — consistent with the project rule that `Integrate` performs no
   `NIntegrate` crosscheck.

### Why the inner integrals use dedicated families

The one delicate part is step 2: the differentiated integral \(J\) is itself a
parameter-dependent definite integral, and the general integrator is slow or
hangs on exactly the half-line exponential/trigonometric forms Feynman's trick
produces. So `DiffUnderInt` **computes those inner integrals with its own
closed-form family evaluators** and never hands them to the general engine. The
families it carries are:

| Inner family | Shape | Closed form |
|---|---|---|
| Laplace / Fourier half-line | \(\int_0^\infty x^n e^{\alpha x}\{1,\cos,\sin\}\,dx\), \(\operatorname{Re}\alpha<0\) | \(n!\,/\,(-\alpha)^{n+1}\) combinations |
| sinc / Frullani | \(\int_0^\infty (\cdots)/x\,dx\) | via \(\int_0^\infty M(s)\,ds\), \(M\) the Laplace image |
| even-rational half-line | \(\int_0^\infty P(x)/Q(x^2)\,dx\) | \(v=x^2\Rightarrow\) Beta integrals |
| general (non-even) rational half-line | \(\int_0^\infty R(s)\,ds\) | real \(\operatorname{ArcTan}\)/\(\log\) boundary values |
| Gaussian moment | \(\int_0^\infty x^n e^{-p x^2}\{1,\cos\}\,dx\) | in \(\sqrt\pi\), \(e^{-q^2/4p}\) |

The Gaussian **parameter** back-integration \(\int c\,e^{-k p^2}\,dp\) is supplied
directly as an `Erf` (the engine produces no such antiderivative). A form outside
every family is declined — the integral comes back unevaluated, fast, never a
wrong value (see [§7](#7-limitations)).

`DiffUnderInt` is tried near the end of the definite `Integrate`
[cascade](index.md), after the residue and Newton–Leibniz methods — it is the
tool for integrals that are *only* tractable through a parameter.

---

## 3. Invoking the method

Pin the method with `Method -> "DiffUnderInt"` (the long alias
`"DifferentiationUnderIntegral"` also works), or call the builtin directly as
`` Integrate`DiffUnderInt[f, {x, a, b}] ``. All three are equivalent:

```mathematica
In[1]:= Integrate[(x^a - 1)/Log[x], {x, 0, 1}, Method -> "DifferentiationUnderIntegral"]
Out[1]= Log[1 + a]

In[2]:= Integrate`DiffUnderInt[(x^a - 1)/Log[x], {x, 0, 1}]
Out[2]= Log[1 + a]
```

The integrand must contain a **free parameter** to differentiate with respect to
— an integral of numbers only is not this method's business. Supply
`Assumptions` to pin the signs of the parameters; they fix the convergence gates
and are used to clean the closed forms:

```mathematica
In[1]:= Integrate[Sin[a x]^2/x^2, {x, 0, Infinity},
          Method -> "DiffUnderInt", Assumptions -> a > 0]
Out[1]= 1/2 Pi a
```

You do not normally need to name the parameter: the method tries each free symbol
in turn and keeps the first that yields a verified closed form.

---

## 4. Worked examples

Every transcript below is reproduced from the current build. Timings (from
`// Timing`) are indicative only.

### 4.1 The logarithm trick

The archetype: a `Log[x]` in the denominator killed by \(\partial_a x^a = x^a\log
x\). With two exponents the answers just superpose (linearity of the parameter
integral):

```mathematica
In[1]:= Integrate[(x^a - 1)/Log[x], {x, 0, 1}, Method -> "DiffUnderInt"] // Timing
Out[1]= {0.005486, Log[1 + a]}

In[2]:= Integrate[(x^a - x^b)/Log[x], {x, 0, 1},
          Method -> "DiffUnderInt", Assumptions -> {a > -1, b > -1}]
Out[2]= Log[1 + a] - Log[1 + b]
```

### 4.2 Frullani integrals

A **Frullani integral** \(\int_0^\infty \frac{g(ax)-g(bx)}{x}\,dx\) is a natural
target: differentiating in \(a\) removes the \(1/x\) and leaves a standard
transform. The classic exponential and arctangent cases:

```mathematica
In[1]:= Integrate[(Exp[-a x] - Exp[-b x])/x, {x, 0, Infinity},
          Method -> "DiffUnderInt", Assumptions -> {a > 0, b > 0}]
Out[1]= -Log[a] + Log[b]

In[2]:= Integrate[(ArcTan[a x] - ArcTan[b x])/x, {x, 0, Infinity},
          Method -> "DiffUnderInt", Assumptions -> {a > 0, b > 0}]
Out[2]= 1/2 Pi (Log[a] - Log[b])
```

The first is the textbook Frullani value \(\log(b/a)\). Attaching an oscillatory
factor \(\cos(cx)\) gives a two-parameter version that still closes:

```mathematica
In[1]:= Integrate[(Exp[-a x] - Exp[-b x])/x Cos[c x], {x, 0, Infinity},
          Method -> "DiffUnderInt",
          Assumptions -> {a > 0, b > 0, Element[c, Reals]}]
Out[1]= -1/2 Log[a^2 + c^2] + 1/2 Log[b^2 + c^2]
```

i.e. \(\tfrac12\log\!\big((b^2+c^2)/(a^2+c^2)\big)\). The split-`Log` form is the
tidy real result of a Laplace-transform inner integral.

### 4.3 Sinc and Dirichlet integrals

Integrands with \(1/x^2\) or \(1/x\) reduce, after one differentiation, to a
**sinc** integral over the half-line. The two standard \(\pi a/2\) results:

```mathematica
In[1]:= Integrate[Sin[a x]^2/x^2, {x, 0, Infinity},
          Method -> "DiffUnderInt", Assumptions -> a > 0] // Timing
Out[1]= {0.063488, 1/2 Pi a}

In[2]:= Integrate[(1 - Cos[a x])/x^2, {x, 0, Infinity},
          Method -> "DiffUnderInt", Assumptions -> a > 0]
Out[2]= 1/2 Pi a
```

Adding an exponential window differentiates to \(\int_0^\infty e^{-ax}\cos(bx)\,dx\),
whose parameter integral is the arctangent:

```mathematica
In[1]:= Integrate[Exp[-a x] Sin[b x]/x, {x, 0, Infinity},
          Method -> "DiffUnderInt", Assumptions -> {a > 0, b > 0}]
Out[1]= 1/2 Pi - ArcTan[a/b]
```

which is \(\operatorname{arctan}(b/a)\) (the two differ by the identity
\(\operatorname{arctan} t+\operatorname{arctan}(1/t)=\pi/2\)).

### 4.4 Rational half-line with a logarithm

Here the parameter sits inside a logarithm; differentiating turns it into a
rational function integrated against \(1/(1+x^2)\):

```mathematica
In[1]:= Integrate[Log[1 + a^2 x^2]/(1 + x^2), {x, 0, Infinity},
          Method -> "DiffUnderInt", Assumptions -> a > 0]
Out[1]= Pi Log[1 + a]

In[2]:= Integrate[Log[1 + a^2 x^2]/(x^2 (1 + x^2)), {x, 0, Infinity},
          Method -> "DiffUnderInt", Assumptions -> a > 0] // Timing
Out[2]= {0.050109, Pi (a - Log[1 + a])}

In[3]:= Integrate[ArcTan[a x]/(x (1 + x^2)), {x, 0, Infinity},
          Method -> "DiffUnderInt", Assumptions -> a > 0]
Out[3]= 1/2 Pi Log[1 + a]
```

### 4.5 A decaying sinc — the non-even rational family

Combining an exponential decay with a \((1-\cos ax)/x^2\) window differentiates to
a **decaying** sinc \(\int_0^\infty e^{-cx}\sin(ax)/x\,dx=\operatorname{arctan}(a/c)\).
Its Laplace image \(M(s)=a/((s+c)^2+a^2)\) is *not even* in \(s\), so Mathilda
integrates it with the general non-even rational half-line evaluator, which
returns a real arctangent directly:

```mathematica
In[1]:= Integrate[Exp[-c x] (1 - Cos[a x])/x^2, {x, 0, Infinity},
          Method -> "DiffUnderInt", Assumptions -> {a > 0, c > 0}] // Timing
Out[1]= {0.653151, a ArcTan[a/c] - 1/2 c (-2 Log[c] + Log[a^2 + c^2])}
```

that is, \(a\operatorname{arctan}(a/c)-\tfrac{c}{2}\log\!\big(1+a^2/c^2\big)\).

### 4.6 A Gaussian integral

A Gaussian weight is no obstacle. Differentiating \(e^{-x^2}\sin(ax)/x\) in \(a\)
removes the \(1/x\) and leaves the **cosine moment** \(\int_0^\infty
e^{-x^2}\cos(ax)\,dx=\tfrac{\sqrt\pi}{2}e^{-a^2/4}\); integrating that back over
\(a\) produces an error function:

```mathematica
In[1]:= Integrate[Exp[-x^2] Sin[a x]/x, {x, 0, Infinity},
          Method -> "DiffUnderInt"] // Timing
Out[1]= {0.15825, 1/2 Pi Erf[1/2 a]}
```

The `Erf` comes from the Gaussian parameter back-integration \(\int
\tfrac{\sqrt\pi}{2}e^{-a^2/4}\,da\), which `DiffUnderInt` supplies in closed form
because the general integrator does not produce it.

---

## 5. Verification and the base point

Because the answer is reconstructed from a derivative, `DiffUnderInt` **proves**
each result before returning it. Two things are checked:

- **The differential equation.** \(\operatorname{Simplify}[\,I'(p)-J\,]\) must
  reduce to the literal `0` under your `Assumptions`. This is a symbolic identity
  check on elementary functions — a genuine proof, not a sampled comparison.
  (`PossibleZeroQ` is deliberately *not* used here: its shrinkage heuristic
  false-positives on decaying expressions, which could mis-accept a wrong result.)

- **The base value.** The additive constant is fixed from an *exact* \(I(p_0)\),
  most often a parameter value at which the integrand is identically zero, so
  \(I(p_0)=0\) with no computation at all.

Together these make the closed form correct by construction. The
conditionally-convergent pitfall that makes naive differentiation under the
integral sign fail (Conrad, §12) is caught automatically: a non-integrable
\(\partial_p f\) simply fails to close, and that parameter is skipped rather than
producing a spurious value.

---

## 6. The families at a glance

The inner definite integral \(J\) is matched against the closed-form families
listed in [§2](#2-how-mathilda-realises-it). Their combined reach is: half-line
integrals of `poly × exp × {1, sin, cos}` (with or without a \(1/x^k\) window),
even and non-even rational functions, and Gaussian moments — plus the
parameter-side integrals of the arctangents, logarithms, rationals, and
Gaussians those produce. That is enough to close the standard corpus of
"Feynman's trick" exercises, as [§4](#4-worked-examples) shows.

---

## 7. Limitations

- **A free parameter is required.** With no symbol to differentiate against, the
  method does not apply and returns the input unevaluated. Reach for another
  method in the [cascade](index.md) (`"NewtonLeibniz"`, `"Residue"`, `"Mellin"`).

- **Piecewise results are declined.** When the true answer is a genuine
  case-split — for example \(\int_0^\infty \frac{\sin(ax)\sin(bx)}{x^2}\,dx =
  \frac{\pi}{2}\min(a,b)\) — the differentiated inner integral is itself
  conditional, which the families do not carry. The integral is returned
  unevaluated, cleanly and fast:

  ```mathematica
  In[1]:= Integrate[Sin[a x] Sin[b x]/x^2, {x, 0, Infinity},
            Method -> "DiffUnderInt", Assumptions -> {a > 0, b > 0}]
  Out[1]= Integrate[(Sin[a x] Sin[b x])/x^2, {x, 0, Infinity},
            Assumptions -> {a > 0, b > 0}, Method -> "DiffUnderInt"]
  ```

- **Some inner integrals are out of family.** Finite-period trigonometric
  integrals over \(\{0,\pi\}\), finite-interval radicals, and half-line
  `log × rational` *base* values fall outside the current evaluators and are
  declined rather than guessed. The purely oscillatory Sine–Gaussian moment
  \(\int_0^\infty e^{-x^2}\sin(ax)\,dx\) (a Dawson/`Erfi` form) is likewise not
  carried — but, as [§4.6](#46-a-gaussian-integral) shows, the cases that
  *differentiate into* a cosine moment still close.

- **No numeric fallback.** Consistent with the rest of `Integrate`, the method
  never approximates: it either returns a proven closed form or the unevaluated
  input.

---

## See also

- [Calculus tutorial](../08-calculus.md) — `Integrate`, `D`, `Series`, `Limit`.
- [Mellin transforms](mellin.md) — the sibling method for half-line
  \(\int_0^\infty x^{s-1}f\,dx\) integrals via a transform table.
- [Residue theorem](residue.md) and [Contour / line integrals](lineintegral.md)
  — the complex-analytic methods in the cascade.
- `` ?Integrate`DiffUnderInt `` in the REPL for the built-in help string.
