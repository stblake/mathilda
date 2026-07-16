# Integration by Mellin transforms

`Method -> "Mellin"` evaluates **half-line definite integrals**

\[
M_f(s) \;=\; \int_0^{\infty} x^{\,s-1}\, f(x)\, dx
\]

in closed form. The quantity \(M_f(s)\) is the **Mellin transform** of \(f\), and
a great many named constants and special functions are, at heart, Mellin
transforms evaluated at a particular \(s\) — the Gamma function, the Beta
function, \(\pi/2\) (Dirichlet's integral), and the whole hypergeometric family
among them.

This method is the engine behind Mathilda's
`Integrate`RamanujanMasterTheorem` builtin; `"Mellin"` and
`"RamanujanMasterTheorem"` are two names for the same method.

```mathematica
In[1]:= Integrate[x^(s - 1) Exp[-x], {x, 0, Infinity}, Method -> "Mellin"]
Out[1]= ConditionalExpression[Gamma[s], s > 0]
```

That one line *is* the definition of the Gamma function, and Mathilda returns it
together with the exact strip of convergence, `s > 0`, wrapped in a
`ConditionalExpression`.

---

## 1. The mathematical idea

### Mellin transforms

The Mellin transform \(M_f(s)=\int_0^\infty x^{s-1}f(x)\,dx\) is the
multiplicative analogue of the Laplace/Fourier transform (substitute
\(x=e^{-t}\) and it becomes a two-sided Laplace transform). It converges on a
vertical **strip** \(\alpha < \operatorname{Re} s < \beta\) in the complex
\(s\)-plane, whose edges are fixed by the growth of \(f\) at the two ends of the
ray:

- if \(f(x)\sim x^{-\alpha}\) as \(x\to 0^+\), the integral needs
  \(\operatorname{Re} s > \alpha\) for convergence at the origin;
- if \(f(x)\sim x^{-\beta}\) as \(x\to\infty\), it needs
  \(\operatorname{Re} s < \beta\) for convergence at infinity.

A closed form is only *valid inside this strip*, which is exactly why Mathilda
reports it. Outside the strip the same algebraic expression is the **analytic
continuation** of the integral, not the integral itself.

### The Ramanujan Master Theorem

For an integrand given by a power series
\(f(x)=\sum_{k\ge 0}\frac{(-1)^k}{k!}\,\varphi(k)\,x^{k}\), Ramanujan's Master
Theorem (RMT) states

\[
\int_0^\infty x^{\,s-1} f(x)\,dx \;=\; \Gamma(s)\,\varphi(-s),
\]

turning integration into **analytic interpolation of the coefficients**: read
off \(\varphi(k)\) from the series, replace \(k \mapsto -s\), and multiply by
\(\Gamma(s)\). For a generalized hypergeometric integrand
\({}_pF_q(\mathbf a;\mathbf b;-\lambda x)\) the coefficient function is
\(\varphi(k)=\prod_i (a_i)_k \big/ \prod_j (b_j)_k\), and RMT delivers the
ratio-of-Gammas transform in one step — see
[§4](#4-the-ramanujan-master-theorem-hypergeometric-integrands).

### References

1. E. C. Titchmarsh, *Introduction to the Theory of Fourier Integrals*, 2nd ed.,
   Oxford, 1948 — Mellin inversion and convergence strips (Ch. VII).
2. G. H. Hardy, *Ramanujan: Twelve Lectures on Subjects Suggested by His Life and
   Work*, Cambridge, 1940 — Ch. XI, the Master Theorem.
3. T. Amdeberhan, O. Espinosa, I. Gonzalez, M. Harrison, V. H. Moll, A. Straub,
   "Ramanujan's Master Theorem," *The Ramanujan Journal* **29** (2012) 103–120.
4. P. Flajolet, X. Gourdon, P. Dumas, "Mellin transforms and asymptotics:
   Harmonic sums," *Theoretical Computer Science* **144** (1995) 3–58.
5. F. Oberhettinger, *Tables of Mellin Transforms*, Springer, 1974.
6. NIST *Digital Library of Mathematical Functions*, §1.14 (integral transforms),
   [dlmf.nist.gov/1.14](https://dlmf.nist.gov/1.14).

---

## 2. How Mathilda realises it

The method is a **table plus a small set of operational rules**, applied
symbolically. There is no numerical quadrature anywhere in it — every answer is
correct by construction.

1. **Split into terms.** The integrand is expanded into a sum
   \(\sum_i c_i\, x^{\rho_i}\, \mathrm{kernel}_i(x)\). Each additive term is
   integrated independently; because every surviving term converges absolutely on
   its own strip, the sum of the integrals equals the integral of the sum.

2. **Peel off the power prefactor.** Each term is written as
   \(C\,x^{\rho}\,\mathrm{kernel}(x)\) with \(C\) free of \(x\). The pure power is
   absorbed into the transform variable by setting \(s = \rho + 1\): the leading
   \(x^{s-1}\) is just bookkeeping for "which point of \(M_{\mathrm{kernel}}\) do
   we want".

3. **Match a base kernel.** The single remaining factor is matched against a
   table of *proven* base Mellin transforms (see [§7](#7-the-kernel-table)). A
   linear internal scaling \(x \mapsto \lambda x\) folds into the transform as the
   \(\lambda^{-s}\) factor; a monomial substitution \(x \mapsto x^{k}\) (as in
   `Cos[Sqrt[x]]`) rescales the strip and the argument of \(\Gamma\) accordingly.

4. **Gate on the convergence strip.** Every table identity carries its strip.
   Mathilda reduces the strip against your `Assumptions` (or evaluates it
   numerically when \(s\) is a number). If the strip is proven, the bare closed
   form is returned; otherwise it is carried as a `ConditionalExpression`. This
   gate is what makes the result *unconditionally correct* rather than a formal
   manipulation.

A product of **two** transcendental kernels (say \(e^{-x}\cos x\)) is a Mellin
*convolution* — Meijer-\(G\) territory — and is deliberately out of scope; the
method returns the integral unevaluated rather than guess (see
[§8](#8-limitations)).

---

## 3. Invoking the method

Pin the method with `Method -> "Mellin"`. The integrand must be presented in the
half-line Mellin shape `x^(s-1) * kernel`, over `{x, 0, Infinity}`.

You do not have to keep `s` symbolic. If the prefactor is a concrete power, the
strip is checked numerically and a bare value comes back:

```mathematica
In[1]:= Integrate[Exp[-x], {x, 0, Infinity}, Method -> "Mellin"]
Out[1]= 1

In[2]:= Integrate[x^2 Exp[-x], {x, 0, Infinity}, Method -> "Mellin"]
Out[2]= 2
```

Here `x^2 Exp[-x]` has \(s = 3\), so the answer is \(\Gamma(3) = 2! = 2\).

When `s` is symbolic, supply the strip yourself through `Assumptions` and
Mathilda will discharge the `ConditionalExpression` for you:

```mathematica
In[1]:= Integrate[x^(s - 1) Exp[-x], {x, 0, Infinity},
          Method -> "Mellin", Assumptions -> s > 0]
Out[1]= Gamma[s]
```

Without the assumption the same call returns
`ConditionalExpression[Gamma[s], s > 0]` — the closed form annotated with the
condition under which it holds.

---

## 4. The Ramanujan Master Theorem (hypergeometric integrands)

The RMT layer handles integrands built from a generalized hypergeometric
function \({}_pF_q\) with a linear argument, provided \(p \le q + 1\). For
\({}_1F_1\) (Kummer's confluent function) the transform is a clean ratio of
Gammas:

```mathematica
In[1]:= Integrate[x^(s - 1) Hypergeometric1F1[a, c, -x], {x, 0, Infinity},
          Method -> "Mellin"]
Out[1]= ConditionalExpression[(Gamma[c] Gamma[s] Gamma[a - s])/(Gamma[a] Gamma[c - s]),
          s > 0 && s < a]
```

and for the Gauss function \({}_2F_1\) there is simply one more Gamma factor top
and bottom:

```mathematica
In[1]:= Integrate[x^(s - 1) Hypergeometric2F1[a, b, c, -x], {x, 0, Infinity},
          Method -> "Mellin"]
Out[1]= ConditionalExpression[
          (Gamma[c] Gamma[s] Gamma[a - s] Gamma[b - s])/(Gamma[a] Gamma[b] Gamma[c - s]),
          s > 0 && s < a && s < b]
```

These are the general pattern

\[
M_{{}_pF_q}(s)=\frac{\prod_j\Gamma(b_j)}{\prod_i\Gamma(a_i)}\;\Gamma(s)\;
\frac{\prod_i\Gamma(a_i-s)}{\prod_j\Gamma(b_j-s)}\;\lambda^{-s},
\]

valid on \(0 < \operatorname{Re} s < \min_i \operatorname{Re} a_i\), produced
directly by \(\varphi(-s)\).

Many special functions are hypergeometric in disguise. Mathilda rewrites a few
of them (`Erf`, `Erfc`, the lower incomplete Gamma, spherical Bessel) into
\({}_pF_q\) form first, so their Mellin transforms come out of the same machine:

```mathematica
In[1]:= Integrate[x^(s - 1) Erf[x], {x, 0, Infinity}, Method -> "Mellin"]
Out[1]= ConditionalExpression[-Gamma[1/2 (1 + s)]/(Sqrt[Pi] s),
          1/2 (1 + s) > 0 && 1/2 (1 + s) < 1/2]
```

---

## 5. Worked examples

Every transcript below is reproduced from the current build. Timings (from
`// Timing`) are indicative only.

### 5.1 The classical transforms

The Gamma and Beta functions are the two workhorse Mellin transforms. The Beta
integral is the transform of an algebraic binomial \((1+x)^{-a}\):

```mathematica
In[1]:= Integrate[x^(s - 1)/(1 + x)^a, {x, 0, Infinity}, Method -> "Mellin"]
Out[1]= ConditionalExpression[Beta[s, a - s], s > 0 && s < a]

In[2]:= Integrate[x^(s - 1)/(1 + x), {x, 0, Infinity}, Method -> "Mellin"]
Out[2]= ConditionalExpression[Beta[s, 1 - s], s > 0 && s < 1]
```

The second is the \(a=1\) special case; recall \(B(s,1-s)=\pi/\sin(\pi s)\) — the
reflection formula in integral form.

The Gaussian is handled by a monomial substitution \(x \mapsto x^2\), which halves
the argument of \(\Gamma\) and the strip:

```mathematica
In[1]:= Integrate[Exp[-x^2], {x, 0, Infinity}, Method -> "Mellin"]
Out[1]= 1/2 Sqrt[Pi]

In[2]:= Integrate[x^(s - 1) Exp[-x^2], {x, 0, Infinity}, Method -> "Mellin"]
Out[2]= ConditionalExpression[1/2 Gamma[1/2 s], s > 0]
```

### 5.2 Trigonometric kernels

The Mellin transforms of \(\cos\) and \(\sin\) are \(\Gamma(s)\cos(\pi s/2)\) and
\(\Gamma(s)\sin(\pi s/2)\) on their strips. Mathilda returns them in the
equivalent reflection form (a \(\pi\csc/\Gamma\) or \(\pi\sec/\Gamma\) ratio),
which is the numerically robust representation near the strip edges:

```mathematica
In[1]:= Integrate[x^(s - 1) Cos[x], {x, 0, Infinity}, Method -> "Mellin"]
Out[1]= ConditionalExpression[(1/2 Pi Csc[1/2 Pi s])/Gamma[1 - s], s > 0 && s < 1]
```

The Dirichlet integral \(\int_0^\infty \frac{\sin x}{x}\,dx=\frac{\pi}{2}\) falls
straight out as the \(s=1\) slice of the \(\sin x / x\) transform:

```mathematica
In[1]:= Integrate[x^(s - 1) Sin[x]/x, {x, 0, Infinity}, Method -> "Mellin"]
Out[1]= ConditionalExpression[(1/2 Pi Sec[1/2 Pi (-1 + s)])/Gamma[2 - s], s > 0 && s < 2]

In[2]:= Integrate[Sin[x]/x, {x, 0, Infinity}, Method -> "Mellin"]
Out[2]= 1/2 Pi
```

`Sin[x]/x` at \(s=1\) hits a removable \(0\cdot\infty\) coincidence in the raw
formula, which Mathilda resolves with a limit before returning \(\pi/2\).

### 5.3 The showcase set

The following integrals exercise every branch of the method — internal monomial
substitutions, the logarithm and arctangent kernels, Bessel functions, and the
hypergeometric layer. Each is returned with its exact strip.

**Exponential — the Gamma function.**

```mathematica
In[1]:= Integrate[x^(s - 1) Exp[-x], {x, 0, Infinity}, Method -> "Mellin"] // Timing
Out[1]= {0.001454, ConditionalExpression[Gamma[s], s > 0]}
```

**Algebraic binomial — the Beta function.**

```mathematica
In[1]:= Integrate[x^(s - 1)/(1 + x)^a, {x, 0, Infinity}, Method -> "Mellin"] // Timing
Out[1]= {0.001405, ConditionalExpression[Beta[s, a - s], s > 0 && s < a]}
```

**A logarithm kernel.** `Log[1 + x]/x` has \(\rho=-1\), shifting \(s\to s-1\)
inside the \(\pi\csc\) transform:

```mathematica
In[1]:= Integrate[x^(s - 1) Log[1 + x]/x, {x, 0, Infinity}, Method -> "Mellin"] // Timing
Out[1]= {0.147829, ConditionalExpression[(Pi Csc[Pi (-1 + s)])/(-1 + s), s > 0 && s < 1]}
```

**Trigonometric kernel with a square-root argument.** `Cos[Sqrt[x]]` triggers the
\(x \mapsto x^2\) monomial substitution, doubling the strip's slope
(\(2s < 1\)) and the \(\Gamma\) argument:

```mathematica
In[1]:= Integrate[x^(s - 1) Cos[Sqrt[x]], {x, 0, Infinity}, Method -> "Mellin"] // Timing
Out[1]= {0.057947, ConditionalExpression[(Pi Csc[Pi s])/Gamma[1 - 2 s], s > 0 && 2 s < 1]}
```

**Another square-root argument, with a prefactor.**

```mathematica
In[1]:= Integrate[x^(s - 1) Sin[Sqrt[x]]/Sqrt[x], {x, 0, Infinity}, Method -> "Mellin"] // Timing
Out[1]= {0.246038, ConditionalExpression[(Pi Sec[-1/2 Pi + Pi s])/Gamma[2 - 2 s], s > 0 && s < 1]}
```

**A Bessel-function kernel.** The Weber–Schafheitlin-type transform of
\(J_\nu(2\sqrt{x})\,x^{-\nu/2}\) collapses to a ratio of two Gammas:

```mathematica
In[1]:= Integrate[x^(s - 1) BesselJ[nu, 2 Sqrt[x]]/x^(nu/2), {x, 0, Infinity}, Method -> "Mellin"] // Timing
Out[1]= {0.012294, ConditionalExpression[Gamma[s]/Gamma[1 + nu - s], s > 0 && -3/2 + 2 s < nu]}
```

**An arctangent kernel.**

```mathematica
In[1]:= Integrate[x^(s - 1) ArcTan[Sqrt[x]]/Sqrt[x], {x, 0, Infinity}, Method -> "Mellin"] // Timing
Out[1]= {0.31456, ConditionalExpression[-(Pi Sec[-1/2 Pi + Pi s])/(-1 + 2 s), s > 0 && 2 s < 1]}
```

**Confluent hypergeometric — Ramanujan's Master Theorem.**

```mathematica
In[1]:= Integrate[x^(s - 1) Hypergeometric1F1[a, c, -x], {x, 0, Infinity}, Method -> "Mellin"] // Timing
Out[1]= {0.008704, ConditionalExpression[(Gamma[c] Gamma[s] Gamma[a - s])/(Gamma[a] Gamma[c - s]),
          s > 0 && s < a]}
```

**Gauss hypergeometric — Ramanujan's Master Theorem.**

```mathematica
In[1]:= Integrate[x^(s - 1) Hypergeometric2F1[a, b, c, -x], {x, 0, Infinity}, Method -> "Mellin"] // Timing
Out[1]= {0.01787, ConditionalExpression[
          (Gamma[c] Gamma[s] Gamma[a - s] Gamma[b - s])/(Gamma[a] Gamma[b] Gamma[c - s]),
          s > 0 && s < a && s < b]}
```

---

## 6. Convergence and `ConditionalExpression`

When `s` (or any other parameter) is symbolic, the closed form is only valid
inside the convergence strip, so Mathilda returns a
`ConditionalExpression[value, condition]`. The `condition` is the strip, reduced
against your `Assumptions`:

- Supply `Assumptions -> (0 < s < a)` and the whole strip is discharged — the
  bare `value` is returned.
- Supply a *partial* assumption and the residual, still-undischarged part of the
  strip is carried in the `ConditionalExpression`.
- Supply nothing and you get the full strip back, as above.

A `ConditionalExpression` collapses to its value automatically once the
condition is known to hold (for instance when the parameters later become
numbers that satisfy it). Treat the condition as an honest statement of *where*
the formula is the integral — outside the strip it is the analytic
continuation, not the integral itself.

---

## 7. The kernel table

The base transforms the method knows, each with its convergence strip
(\(s\) = Mellin variable, \(\lambda,c > 0\)):

| Kernel \(f(x)\) | Mellin transform \(M_f(s)\) | Strip |
|---|---|---|
| \(e^{c x}\) | \(\Gamma(s)\,(-c)^{-s}\) | \(0 < \operatorname{Re} s\) |
| \(e^{c x^2}\) | \(\tfrac12(-c)^{-s/2}\,\Gamma(s/2)\) | \(0 < \operatorname{Re} s\) |
| \((p+q x^{m})^{-a}\) | \(\tfrac1m\, p^{s/m-a}\, q^{-s/m}\,B(s/m,\,a-s/m)\) | \(0 < \operatorname{Re} s < m\operatorname{Re} a\) |
| \(\cos(\lambda x)\) | \(\Gamma(s)\cos(\pi s/2)\,\lambda^{-s}\) | \(0 < \operatorname{Re} s < 1\) |
| \(\sin(\lambda x)\) | \(\Gamma(s)\sin(\pi s/2)\,\lambda^{-s}\) | \(-1 < \operatorname{Re} s < 1\) |
| \(\log(1+\lambda x)\) | \(\dfrac{\pi}{s\,\sin(\pi s)}\,\lambda^{-s}\) | \(-1 < \operatorname{Re} s < 0\) |
| \(\arctan(\lambda x)\) | \(-\dfrac{\pi}{2 s\,\cos(\pi s/2)}\,\lambda^{-s}\) | \(-1 < \operatorname{Re} s < 0\) |
| \(J_\nu(\lambda x)\) | \(2^{s-1}\lambda^{-s}\,\dfrac{\Gamma((\nu+s)/2)}{\Gamma((\nu-s)/2+1)}\) | \(-\operatorname{Re}\nu < \operatorname{Re} s < 3/2\) |
| \(\operatorname{Li}_\nu(-\lambda x)\) | \(-\pi(-s)^{-\nu}\lambda^{-s}/\sin(\pi s)\) | \(-1 < \operatorname{Re} s < 0\) |
| \({}_pF_q(\mathbf a;\mathbf b;-\lambda x)\), \(p\le q{+}1\) | \(\dfrac{\prod_j\Gamma(b_j)}{\prod_i\Gamma(a_i)}\,\Gamma(s)\,\dfrac{\prod_i\Gamma(a_i-s)}{\prod_j\Gamma(b_j-s)}\,\lambda^{-s}\) | \(0<\operatorname{Re} s<\min_i\operatorname{Re} a_i\) |

Two operational rules extend the table:

- **Power prefactor** \(x^\rho\): shifts the transform argument,
  \(s \mapsto s+\rho\) (equivalently, sets \(s=\rho+1\) for a bare kernel).
- **Monomial substitution** \(x \mapsto x^{k}\): replaces \(s\) by \(s/k\) and
  rescales the strip — this is how `Cos[Sqrt[x]]`, `Sin[Sqrt[x]]/Sqrt[x]`, and
  `Exp[-x^2]` are reduced to table entries.

---

## 8. Limitations

- **One transcendental kernel per term.** A product of two transcendental
  factors is a Mellin convolution and is out of scope. Such an integral is
  returned unevaluated rather than approximated:

  ```mathematica
  In[1]:= Integrate[x^(s - 1) Exp[-x] Cos[x], {x, 0, Infinity}, Method -> "Mellin"]
  Out[1]= Integrate[Cos[x] E^(-x) x^(-1 + s), {x, 0, Infinity}, Method -> "Mellin"]
  ```

  (This particular integral does have a closed form via a contour or
  Laplace-transform argument; the Mellin *table* method simply does not carry
  convolution rules.)

- **Half-line only.** The method integrates over `{x, 0, Infinity}`. For finite
  or doubly-infinite ranges use the other methods in the
  [cascade](index.md) (`"NewtonLeibniz"`, `"Residue"`, `"LineIntegral"`, …).

- **The strip must be provable.** If neither the numeric value of `s` nor your
  `Assumptions` can settle the convergence strip, the strip is carried as a
  `ConditionalExpression`; the method never silently drops a convergence
  condition.

---

## See also

- [Calculus tutorial](../08-calculus.md) — `Integrate`, `D`, `Series`, `Limit`.
- [Special functions tutorial](../10-special-functions.md) — `Gamma`, `Beta`,
  `PolyLog`, `BesselJ`, and the hypergeometric family that appear as answers
  here.
- `?Integrate`RamanujanMasterTheorem` in the REPL for the built-in help string.
