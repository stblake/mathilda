# Integration by the residue theorem

`Method -> "Residue"` evaluates the classical families of **improper and periodic
real integrals** that complex analysis dispatches by *summing residues*. The idea
is Cauchy's: close the real integration path into a contour in the complex plane,
and the integral equals \(2\pi i\) times the sum of the residues the contour
encloses.

```mathematica
In[1]:= Integrate[1/(x^2 + 1), {x, -Infinity, Infinity}, Method -> "Residue"]
Out[1]= Pi
```

The integrand \(1/(x^2+1)\) has poles at \(x = \pm i\); closing into the upper
half-plane encircles the single pole \(x = i\), whose residue is \(1/(2i)\), and
\(2\pi i \cdot \tfrac{1}{2i} = \pi\).

This method is the engine behind Mathilda's `Integrate`ContourResidue` builtin.
It runs in the automatic cascade (just before the Newton–Leibniz antiderivative
path); pinning `Method -> "Residue"` forces it and skips the rest.

---

## 1. The mathematical idea

### The residue theorem

If \(f\) is meromorphic inside and on a positively-oriented simple closed contour
\(C\), with no poles on \(C\), then

\[
\oint_C f(z)\,dz \;=\; 2\pi i \sum_k \operatorname{Res}(f, z_k),
\]

the sum running over the poles \(z_k\) enclosed by \(C\). The **residue** at a
pole \(z_0\) is the coefficient of \((z - z_0)^{-1}\) in the Laurent expansion of
\(f\) about \(z_0\); Mathilda computes it symbolically with the same primitive
that backs the [`Residue`](../10-special-functions.md) builtin.

To turn a *real* integral over the line or a period into a contour integral, you
add an arc (or another edge) that closes the path, then argue that the added
piece contributes nothing — or a known amount — in the limit. Which contour you
add is what distinguishes the families below.

### Jordan's lemma and closing the contour

For a **rational** integrand that decays like \(1/z^2\) or faster, a large
semicircle of radius \(R\) in the upper half-plane contributes \(O(1/R) \to 0\),
so the line integral is \(2\pi i\) times the residues in the upper half-plane.

For a **Fourier** integrand \(R(x)\, e^{i a x}\) with \(a > 0\), *Jordan's lemma*
guarantees the upper semicircle vanishes even when \(R\) decays only like
\(1/z\), because \(|e^{i a z}| = e^{-a\,\operatorname{Im} z}\) is exponentially
small in the upper half-plane. This is what makes oscillatory integrals like
\(\int \cos(x)/(x^2+1)\,dx\) tractable.

Other integrands call for other contours: a **wedge** of angle \(2\pi/n\) for
\(1/(1+x^n)\), a **keyhole** wrapping the branch cut of \(x^p\) along
\([0,\infty)\), or a **rectangle** of height \(2\pi i\) for integrands periodic in
the imaginary direction.

### References

1. E. T. Whittaker, G. N. Watson, *A Course of Modern Analysis*, 4th ed.,
   Cambridge, 1927 — Ch. VI, the calculus of residues.
2. L. V. Ahlfors, *Complex Analysis*, 3rd ed., McGraw-Hill, 1979 — Ch. 4–5
   (residue theorem, evaluation of definite integrals).
3. G. F. Carrier, M. Krook, C. E. Pearson, *Functions of a Complex Variable:
   Theory and Technique*, McGraw-Hill, 1966 — keyhole and sector contours.
4. M. J. Ablowitz, A. S. Fokas, *Complex Variables: Introduction and
   Applications*, 2nd ed., Cambridge, 2003 — Ch. 4, Jordan's lemma.
5. NIST *Digital Library of Mathematical Functions*, §1.10(iv) (the residue
   theorem), [dlmf.nist.gov/1.10](https://dlmf.nist.gov/1.10).

---

## 2. How Mathilda realises it

The method is a set of **narrow, conjunctive recognizers** — one per standard
contour. Each recognizer (i) checks the integration bounds, (ii) reduces the
integrand to a rational function whose poles it finds with `Solve`, (iii)
classifies each pole against the family's contour (upper half-plane, unit disk,
or real axis) by the numeric sign of \(\operatorname{Im}\) or the value of
\(|\cdot|\), (iv) sums the enclosed residues, and (v) closes the sum to a scalar.

Because the residue theorem makes each family's value **correct by construction**
once its structural gates hold, there is **no numerical-quadrature crosscheck**
anywhere in the method. The only post-hoc gate is a self-consistency check:

- the closed form must reduce to a scalar (no surviving `x` or `Root`), and
- for the real-valued families, its imaginary part must vanish — a residual
  \(\operatorname{Im}\) betrays a mis-classified pole.

Any failure returns unevaluated, and the automatic cascade's Newton–Leibniz path
takes over; no existing result is weakened.

**Symbolic parameters.** When you pass `Assumptions -> ...`, the recognizers
classify parameter-dependent poles and frequencies by reading their sign at a
single generic point of the region the assumptions pin, while the residue
*arithmetic stays symbolic* — so \(\int \cos(kx)/(x^2+a^2)\,dx\) comes out as a
function of \(a\) and \(k\). A parameter the assumptions leave two-sided unbounded
is refused, since its sign is undetermined.

---

## 3. Invoking the method

Pin the method with `Method -> "Residue"`. The integral must be over the whole
line `{x, -Infinity, Infinity}`, a half-line `{x, 0, Infinity}`, or a full period
such as `{x, 0, 2 Pi}` — these are the ranges the standard contours close.

A **finite, generic** range has no closing contour, so the method declines and
returns the integral unevaluated (in the automatic cascade, an antiderivative
method would handle it instead):

```mathematica
In[1]:= Integrate[1/(x^2 + 1), {x, 0, 5}, Method -> "Residue"]
Out[1]= Integrate[1/(1 + x^2), {x, 0, 5}, Method -> "Residue"]
```

Supply `Assumptions` to pin the signs of any free parameters — see
[§5](#5-worked-examples).

---

## 4. Introductory examples

### 4.1 Rational integrands (Family A)

A rational function decaying at least as fast as \(1/z^2\) closes into the upper
half-plane; the value is \(2\pi i\) times the residues there.

```mathematica
In[1]:= Integrate[1/(x^4 + 1), {x, -Infinity, Infinity}, Method -> "Residue"]
Out[1]= Pi/Sqrt[2]

In[2]:= Integrate[1/(x^6 + 1), {x, -Infinity, Infinity}, Method -> "Residue"]
Out[2]= 2/3 Pi
```

When the integrand is **even**, the half-line integral is exactly half the
full-line value, and the method uses that shortcut:

```mathematica
In[1]:= Integrate[1/(x^2 + 1), {x, 0, Infinity}, Method -> "Residue"]
Out[1]= 1/2 Pi

In[2]:= Integrate[1/(x^4 + 1), {x, 0, Infinity}, Method -> "Residue"]
Out[2]= (1/2 Pi)/Sqrt[2]
```

### 4.2 Fourier / Jordan integrands (Family B)

For \(R(x)\) rational times \(\cos\), \(\sin\), or \(e^{i\,\cdot}\) with a positive
frequency, Jordan's lemma closes the contour upward. Mathilda evaluates the
complex integral \(J = 2\pi i \sum_{\text{UHP}} \operatorname{Res}[R(z)\,e^{i a z}]\)
and takes \(\operatorname{Re}[J]\) for a cosine kernel, \(\operatorname{Im}[J]\)
for a sine kernel:

```mathematica
In[1]:= Integrate[Cos[x]/(x^2 + 1), {x, -Infinity, Infinity}, Method -> "Residue"]
Out[1]= Pi/E

In[2]:= Integrate[x Sin[x]/(x^2 + 1), {x, -Infinity, Infinity}, Method -> "Residue"]
Out[2]= Pi/E
```

### 4.3 Trigonometric integrands over a period (Family C)

A function rational in \(\sin x\) and \(\cos x\), integrated over a full period,
maps under \(z = e^{i x}\) to a rational contour integral around the **unit
circle**; the value is \(2\pi i\) times the residues inside the disk.

```mathematica
In[1]:= Integrate[1/(2 + Cos[x]), {x, 0, 2 Pi}, Method -> "Residue"]
Out[1]= (2 Pi)/Sqrt[3]

In[2]:= Integrate[1/(5 - 4 Cos[x]), {x, 0, 2 Pi}, Method -> "Residue"]
Out[2]= 2/3 Pi
```

---

## 5. Worked examples

The set below exercises every contour the method knows — the rational and Fourier
families, a principal-value indentation, the rectangular strip, the sector wedge
(with a *symbolic* exponent), and the keyhole branch cut. Each transcript is
reproduced from the current build.

**Rational, upper half-plane.** The simplest residue integral: one pole at
\(x = i\).

```mathematica
In[1]:= Integrate[1/(x^2 + 1), {x, -Infinity, Infinity}, Method -> "Residue"]
Out[1]= Pi
```

**Rational, four poles.** The two upper-half-plane roots of \(x^4 = -1\)
contribute.

```mathematica
In[1]:= Integrate[1/(x^4 + 1), {x, -Infinity, Infinity}, Method -> "Residue"]
Out[1]= Pi/Sqrt[2]
```

**Fourier integrand with parameters.** With \(a > 0\) and \(k > 0\) the contour
closes upward; the single pole at \(x = i a\) gives the exponential decay
\(e^{-a k}\) — a function of both parameters, produced by symbolic residue
arithmetic:

```mathematica
In[1]:= Integrate[Cos[k*x]/(x^2 + a^2), {x, -Infinity, Infinity},
          Assumptions -> {a > 0, k > 0}, Method -> "Residue"]
Out[1]= (Pi E^(-a k))/a
```

**Principal value — a removable pole on the axis.** `Sin[x]/x` has a *removable*
singularity at \(x = 0\); the indented contour skirts it with a half-circle, so
the pole contributes a **half** residue \(\pi i \operatorname{Res}\). The result
is Dirichlet's integral:

```mathematica
In[1]:= Integrate[Sin[x]/x, {x, -Infinity, Infinity}, Method -> "Residue"]
Out[1]= Pi
```

**Rectangular strip — quasi-periodic integrand.** \(e^{a x}/(e^x + 1)\) is
periodic in the imaginary direction with period \(2\pi i\). A rectangle of that
height reduces (via \(w = e^x\)) to a single enclosed pole at \(w = -1\), and the
two horizontal edges differ by the constant factor \(1 - e^{2\pi i a}\). With
\(0 < a < 1\) the value is the reflection formula:

```mathematica
In[1]:= Integrate[Exp[a*x]/(Exp[x] + 1), {x, -Infinity, Infinity},
          Assumptions -> 0 < a < 1, Method -> "Residue"]
Out[1]= Pi Csc[Pi a]
```

**Sector wedge — a symbolic exponent.** \(1/(1+x^n)\) on the half-line closes into
a wedge of angle \(2\pi/n\). This is the one family that admits a **symbolic**
exponent \(n\): with the convergence assumption \(n > 1\) the value comes out in
closed form in \(n\):

```mathematica
In[1]:= Integrate[1/(1 + x^n), {x, 0, Infinity}, Assumptions -> n > 1, Method -> "Residue"]
Out[1]= (Pi Csc[Pi/n])/n
```

**Keyhole — a branch-power integrand.** \(x^{1/3}/(x^2+1)\) carries the branch cut
of \(x^{1/3}\) along \([0,\infty)\). A keyhole contour wraps that cut; the jump
across it lands the value on a cosecant (here evaluating to a surd):

```mathematica
In[1]:= Integrate[x^(1/3)/(x^2 + 1), {x, 0, Infinity}, Method -> "Residue"]
Out[1]= Pi/Sqrt[3]
```

The keyhole and sector families share one residue engine,
\(\int_0^\infty x^p R(x)\,dx\) reducing in both cases to a Mellin-type sum over
the poles of \(R\) weighted by the branch value \(x^{s-1}\); two more half-line
branch-power integrals:

```mathematica
In[1]:= Integrate[Sqrt[x]/(x^2 + 1), {x, 0, Infinity}, Method -> "Residue"]
Out[1]= Pi/Sqrt[2]

In[2]:= Integrate[1/(1 + x^3), {x, 0, Infinity}, Method -> "Residue"]
Out[2]= (2/3 Pi)/Sqrt[3]
```

---

## 6. Convergence, principal values, and divergence

- **Genuine axis poles diverge.** A pure rational integrand with a
  *non-removable* pole on the real axis has no convergent improper integral (only
  a principal value would exist, which plain `Integrate` does not compute). The
  method detects this and reports divergence rather than returning a wrong finite
  number.

- **Removable axis poles give a principal value.** A simple axis pole whose kernel
  supplies a matching zero — like `Sin[x]/x` at \(0\) — is admitted as an indented
  contour and contributes a half residue, as in [§5](#5-worked-examples).

- **The value must close to a real scalar.** For the real-valued families the
  method verifies that the summed residues produce a scalar with vanishing
  imaginary part. If a pole was mis-classified the check fails and the integral is
  left for another method — the method never returns an unverified value.

---

## 7. The contour table

| Family | Surface form | Contour | Value |
|---|---|---|---|
| A — rational | \(R(x)\) on \((-\infty, \infty)\) | upper semicircle | \(2\pi i \sum_{\text{UHP}} \operatorname{Res}\) |
| B — Fourier / Jordan | \(R(x)\,\{\cos,\sin,e^{i\cdot}\}(a x)\) on \((-\infty, \infty)\) | upper semicircle (Jordan) | \(\operatorname{Re}/\operatorname{Im}\) of \(2\pi i \sum_{\text{UHP}} \operatorname{Res}[R\,e^{i a x}]\) |
| C — trigonometric | rational in \(\sin x, \cos x\) over a period | unit circle, \(z = e^{i x}\) | \(2\pi i \sum \operatorname{Res}\) inside the disk |
| Principal value | simple removable pole on the axis (A/B) | indented axis | \(+\, \pi i \operatorname{Res}\) at the axis pole |
| Half-line (even) | even \(R(x)\) on \([0, \infty)\) | — | \(\tfrac12 \times\) the full-line value |
| Sector | \(x^m / (c + x^n)\) on \([0, \infty)\), \(n\) possibly symbolic | wedge of angle \(2\pi/n\) | \(\tfrac{\pi}{n} c^{s/n-1} \csc(\pi s/n)\), \(s = m{+}1\) |
| Keyhole / Mellin | \(x^p R(x)\) on \([0, \infty)\), \(p\) non-integer | keyhole around \([0,\infty)\) | Mellin sum over poles of \(R\), \(\times\, \tfrac{1}{1 - e^{2\pi i s}}\) |
| Rectangular strip | \(e^{c x} R(e^x)\) on \((-\infty, \infty)\) | rectangle, height \(2\pi i\) | reduced to the keyhole core by \(w = e^x\) |

All the algebraic families close their residue sums with `RootReduce`; the
Fourier family closes via the \(\operatorname{Re}\)/\(\operatorname{Im}\) of the
Jordan integral.

---

## 8. Limitations

- **Improper or periodic ranges only.** The method needs a closing contour, so a
  generic finite range is declined (see [§3](#3-invoking-the-method)).

- **Logarithmic integrands are largely out of scope.** The keyhole engine carries
  a single branch *power* \(x^p\), not the \(\log\)-weighted integrands whose
  keyhole picks up a \((\log)^2\) or a squared pole. Both of these return
  unevaluated:

  ```mathematica
  In[1]:= Integrate[Log[x]/(x^2 + 1)^2, {x, 0, Infinity}, Method -> "Residue"]
  Out[1]= Integrate[Log[x]/(1 + x^2)^2, {x, 0, Infinity}, Method -> "Residue"]

  In[2]:= Integrate[Log[1 + x^2]/(1 + x^2), {x, 0, Infinity}, Method -> "Residue"]
  Out[2]= Integrate[Log[1 + x^2]/(1 + x^2), {x, 0, Infinity}, Method -> "Residue"]
  ```

  The first pairs a logarithm with a double pole; the second has branch points at
  \(x = \pm i\) *on* the closing contour. Both are classical residue exercises, but
  outside the families this method recognises — a natural direction for a future
  extension.

- **Undetermined parameter signs are refused.** A symbolic parameter that
  `Assumptions` leaves two-sided unbounded cannot be classified against the
  contour, so the family declines rather than guess a sign.

---

## See also

- [Mellin transforms](mellin.md) — the companion half-line method; the keyhole
  family here shares its Mellin-sum core.
- [Contour / line integrals](lineintegral.md) — the path-based counterpart, which
  evaluates a contour integral segment by segment.
- [Calculus tutorial](../08-calculus.md) — `Integrate`, `D`, `Series`, `Limit`.
- `?Integrate`ContourResidue` in the REPL for the built-in help string.
