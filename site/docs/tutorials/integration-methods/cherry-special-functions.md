# Cherry's special-function extensions

The [previous tutorial](risch-transcendental.md) ended on an uncomfortable note:
the Risch algorithm can *prove* that \(e^{x^2}\), \(e^x/x\), \(1/\log x\), and
\(\sin x/x\) have **no elementary antiderivative**. Yet every one of them is a
perfectly well-behaved integral that shows up throughout mathematics and physics.
The resolution, worked out by **G. W. Cherry** in a series of papers in the
1980s, is to *enlarge* the class of allowed answers with a fixed set of
**special functions** — the error function \(\mathrm{erf}\), the exponential
integral \(\mathrm{Ei}\), the logarithmic integral \(\mathrm{li}\), and the
dilogarithm — and to extend Risch's decision procedure to compute in that larger
class.

```mathematica
In[1]:= Integrate[Exp[x]/x, x]
Out[1]= ExpIntegralEi[x]
```

Liouville's theorem forbids an elementary answer here; Cherry's theory supplies
the *right* one. Mathilda's Cherry engines run inside the transcendental Risch
recursion, so they fire automatically from ordinary `Integrate[f, x]` — you do
not need a `Method` option. This tutorial shows what each engine can do, from a
one-line introduction to the hard cases.

---

## 1. The idea: an extended Liouville form

Cherry's generalisation of Liouville's theorem allows the antiderivative to
contain, in addition to the elementary part, a **constant-coefficient sum of one
distinguished special function** applied to field elements. For the exponential
integral, the *extended Liouville form* of an integrand \(g\,e^{f}\)
(\(f, g \in \mathbb{C}(x)\)) is

\[
\int g\,e^{f}\,dx \;=\; y\,e^{f} \;+\; \sum_i c_i\,e^{-\alpha_i}\,
\operatorname{Ei}(f + \alpha_i),
\qquad y \in \mathbb{C}(x),\; c_i, \alpha_i \text{ constant.}
\]

Differentiating and dividing out the common exponential kernel turns this into a
purely **rational** matching identity

\[
g \;=\; y' + f'\,y + \sum_i c_i \,\frac{f'}{f + \alpha_i},
\]

a linear system over the unknown constants \(c_i\) and the coefficients of \(y\),
closed by `SolveAlways`. This is Cherry's *undetermined-coefficient* solve — not
a Risch differential equation — and it is what makes the exponential-integral and
error-function cases decidable together. The logarithmic-integral and
dilogarithm cases replace this with a **\(\Sigma\)-decomposition** (Cherry 1986,
Thm 4.4) that matches *in the logarithm tower*, treating each \(\log\) kernel as
an independent variable. In every case the exact tower identity is the
certificate; the answer is additionally re-checked by a branch-safe `PowerExpand`
diff-back before it is emitted.

### References

1. G. W. Cherry, *Integration in finite terms with special functions: the error
   function*, Journal of Symbolic Computation **1** (1985), 283–302.
2. G. W. Cherry, *Integration in finite terms with special functions: the
   logarithmic integral*, SIAM Journal on Computing **15** (1986), 1–21.
3. G. W. Cherry, *An analysis of the rational exponential integral*, SIAM
   Journal on Computing **18** (1989), 893–905.
4. G. W. Cherry, *Algorithms for integrating elementary functions in terms of
   special functions*, Ph.D. thesis, University of Delaware, 1983.
5. M. Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed.,
   Springer, 2005 — the surrounding Risch framework these engines extend.
6. M. Bronstein, *Integration of elementary functions*, Journal of Symbolic
   Computation **9** (1990), 117–173.

---

## 2. The exponential integral, \(\mathrm{Ei}\)

The exponential integral \(\operatorname{Ei}(x) = \int_{-\infty}^{x}
e^{t}/t\,dt\) is the special function for rational multiples of \(e^{f}\) whose
elementary part is absent or incomplete. The introductory case is a single term:

```mathematica
In[1]:= Integrate[Exp[x]/x, x]
Out[1]= ExpIntegralEi[x]

In[2]:= Integrate[Exp[2*x]/x, x]
Out[2]= ExpIntegralEi[2 x]
```

The engine solves for the elementary part \(y\) and the \(\mathrm{Ei}\)
coefficients **simultaneously**, so it handles integrands that are part
elementary, part special:

```mathematica
In[3]:= Integrate[Exp[x]/x^2, x]
Out[3]= -E^x/x + ExpIntegralEi[x]

In[4]:= Integrate[Exp[x]/x + Exp[x], x]
Out[4]= E^x + ExpIntegralEi[x]
```

In `Out[3]` the algorithm found the rational \(y = -1/x\) *and* the residual
\(\operatorname{Ei}(x)\) from one linear system. Differentiating back,
\(\frac{d}{dx}\!\left(-e^x/x + \operatorname{Ei}(x)\right) = e^x/x^2\). ✓

---

## 3. The error function, \(\mathrm{erf}\)

The same 1989 exponential-integral engine emits the **error-function** term of
the extended Liouville form when the exponent is a perfect square. This is the
Gaussian that §3.2 of the previous tutorial proved non-elementary:

```mathematica
In[1]:= Integrate[Exp[-x^2], x]
Out[1]= 1/2 Sqrt[Pi] Erf[x]

In[2]:= Integrate[Exp[x^2], x]
Out[2]= (1/2*I) Sqrt[Pi] Erf[-I x]
```

`Out[2]` is the imaginary-error-function form
\(\tfrac{\sqrt\pi}{2}\,\mathrm{erfi}(x)\) written through `Erf`; differentiating
recovers \(e^{x^2}\) exactly. The engine that decided **False** for
`Risch`ElementaryIntegralQ[Exp[x^2], x]` is the very one that now supplies the
closed form — the decision and the construction are two faces of Cherry's theory.

---

## 4. The logarithmic integral, \(\mathrm{li}\)

The logarithmic integral \(\operatorname{li}(x) = \int_0^{x} dt/\log t\) is the
special function for the *dual* tower — integrands rational in \(x\) over a single
logarithm \(\log w\). The introductory case is again a single term:

```mathematica
In[1]:= Integrate[1/Log[x], x]
Out[1]= LogIntegral[x]
```

The prime-counting connection makes the next one memorable — the density
\(x/\log x\) of the prime number theorem integrates to \(\operatorname{li}(x^2)\),
found by Cherry's degree-1 \(\Sigma\)-decomposition over the generator \(w = x\):

```mathematica
In[2]:= Integrate[x/Log[x], x]
Out[2]= LogIntegral[x^2]
```

And the two special functions meet: \(\operatorname{li}\) and
\(\operatorname{Ei}\) are related by \(\operatorname{li}(x) =
\operatorname{Ei}(\log x)\), which the engine uses to split a mixed integrand
into an elementary piece plus an \(\mathrm{Ei}\) of a logarithm:

```mathematica
In[3]:= Integrate[1/Log[x]^2, x]
Out[3]= -x/Log[x] + ExpIntegralEi[Log[x]]
```

Differentiating back, \(\frac{d}{dx}\!\left(-x/\log x +
\operatorname{Ei}(\log x)\right) = 1/\log^2 x\). ✓

---

## 5. The sine and cosine integrals, \(\mathrm{Si}\) and \(\mathrm{Ci}\)

Because \(\sin\) and \(\cos\) are exponentials of an imaginary argument, the
oscillatory integrals \(\int \sin x/x\) and \(\int \cos x/x\) are the
exponential-integral case run over \(\mathbb{C}(i)(x)\), and come back as the
**sine and cosine integrals**:

```mathematica
In[1]:= Integrate[Sin[x]/x, x]
Out[1]= SinIntegral[x]

In[2]:= Integrate[Cos[x]/x, x]
Out[2]= CosIntegral[x]
```

The quadratic-argument oscillators land on the **Fresnel** integrals — the
optical-diffraction functions — by the same complex-exponential route:

```mathematica
In[3]:= Integrate[Sin[x^2], x]
Out[3]= Sqrt[1/2 Pi] FresnelS[Sqrt[2/Pi] x]
```

---

## 6. The dilogarithm

The hardest of Cherry's cases is the **dilogarithm**
\(\operatorname{Li}_2(z) = -\int_0^{z} \frac{\log(1-t)}{t}\,dt\), which arises
from integrands \(R(x)\log w\) that are *not* elementary. Mathilda's engine
(a degree-2 \(\Sigma\)-decomposition) searches the rational roots of the linear
factors for the dilogarithm arguments and matches in the log tower. The
introductory identities are the classical ones:

```mathematica
In[1]:= Integrate[Log[1 - x]/x, x]
Out[1]= -PolyLog[2, x]

In[2]:= Integrate[Log[1 + x]/x, x]
Out[2]= -PolyLog[2, -x]
```

(Mathilda writes the dilogarithm as `PolyLog[2, z]`.) The harder cases mix a
\(\log\cdot\log\) elementary product with the dilogarithm — the engine finds
**both** parts at once:

```mathematica
In[3]:= Integrate[Log[x]/(1 - x), x]
Out[3]= PolyLog[2, -(-1 + x)]

In[4]:= Integrate[Log[x]/(1 + x), x]
Out[4]= Log[x] Log[1 + x] + PolyLog[2, -x]
```

`Out[3]` is \(\operatorname{Li}_2(1-x)\); `Out[4]` is the celebrated
\(\log x \log(1+x) + \operatorname{Li}_2(-x)\). Differentiating `Out[4]`:
\(\frac{d}{dx}\!\left(\log x\log(1+x) + \operatorname{Li}_2(-x)\right)
= \log x/(1+x)\). ✓

---

## 7. When even special functions run out

Cherry's class is finite: \(\mathrm{erf}\), \(\mathrm{Ei}\), \(\mathrm{li}\),
\(\mathrm{Si}\), \(\mathrm{Ci}\), the Fresnel integrals, and the dilogarithm.
Integrands whose antiderivative needs something *outside* this class — a genuine
iterated exponential, say — are still returned unevaluated, honestly:

```mathematica
In[1]:= Integrate[Exp[Exp[x]], x]
Out[1]= Integrate[E^E^x, x]
```

There is no closed form for \(\int e^{e^x}\,dx\) in elementary functions *or*
Cherry's special functions, and Mathilda declines to invent one. This is the same
discipline as the Risch decision procedure: a result is emitted only behind a
proof, and silence is the correct answer when no proof exists.

!!! tip "Pairing the two tutorials"
    Use `Risch`ElementaryIntegralQ[f, x]` (previous tutorial) to *ask whether* an
    integral is elementary, and plain `Integrate[f, x]` to *get* the closed form —
    elementary when one exists, and a Cherry special function when Liouville's
    theorem rules the elementary answer out. Together they turn "can this be
    integrated?" into a question with a definite, provable answer.
