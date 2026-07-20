# The transcendental Risch algorithm

`Method -> "RischTranscendental"` is Mathilda's implementation of the **recursive
Risch algorithm** for transcendental elementary functions — the true *decision
procedure* for integration in finite terms. Given a function built from rational
operations, exponentials, and logarithms, it either returns an elementary
antiderivative or **proves that none exists**.

```mathematica
In[1]:= Integrate[(Log[x] - 1)/Log[x]^2, x, Method -> "RischTranscendental"]
Out[1]= x/Log[x]
```

That answer is not a lucky pattern match. The integrand \((\log x - 1)/\log^2 x\)
looks nothing like its antiderivative \(x/\log x\), and no table of forms would
list it. The Risch algorithm *derives* it, and — this is the key property — it is
**correct by construction**: every branch fires only behind an exact structural
certificate that already proves the closed form it emits. Unlike a heuristic
search, the result is never checked by differentiation, because the mathematics
guarantees it.

This distinguishes `"RischTranscendental"` from its cousin
[`"RischNorman"`](../08-calculus.md) (the Bronstein–Norman *parallel* Risch, a
fast heuristic ansatz): the parallel method guesses a template and solves for
coefficients; the recursive method here is the complete Risch decision
procedure, and it never falls back on guessing.

---

## 1. The mathematical idea

### Differential fields and towers

Risch's insight was to make integration an *algebraic* question. A **differential
field** is a field \(K\) with a derivation \(D\) satisfying \(D(ab) = a\,Db +
b\,Da\). The rational functions \(\mathbb{C}(x)\) with \(D = d/dx\) are the base.
Every transcendental elementary function is then reached by stacking a **tower**
of *monomial* extensions

\[
\mathbb{C}(x) = K_0 \subset K_1 \subset \cdots \subset K_n,
\qquad K_{i} = K_{i-1}(t_i),
\]

where each new generator \(t_i\) is either

- **logarithmic** over \(K_{i-1}\): \(Dt_i = Du/u\) for some \(u \in K_{i-1}\)
  (that is, \(t_i = \log u\)), or
- **exponential** over \(K_{i-1}\): \(Dt_i/t_i = Du\) for some \(u \in K_{i-1}\)
  (that is, \(t_i = e^{u}\)).

So \(\tan x\) enters as an exponential kernel (via \(e^{ix}\)), \(\log\log x\) as
a logarithm over \(\mathbb{C}(x)(\log x)\), and \(e^{x^2}\) as an exponential
whose "inner" derivative \(2x\) already lives in \(\mathbb{C}(x)\).

### Liouville's theorem

The engine that makes the problem *decidable* is a theorem of Liouville: if
\(f\) has an elementary antiderivative at all, then it has one of the very
restricted shape

\[
\int f = v_0 + \sum_{k} c_k \log v_k,
\qquad v_k \in K,\; c_k \text{ constant.}
\]

No genuinely new functions can appear — only elements of the field \(K\) you
started in, plus a sum of logarithms with *constant* coefficients. The Risch
algorithm is the constructive search for such a \(v_0\) and the \(v_k\); when the
search provably fails, the integral is **not elementary**.

### The three cases and the Risch differential equation

Working down from the top monomial \(t = t_n\), an integrand is a rational
function of \(t\) over the field below. The algorithm splits into:

- a **Hermite reduction** that removes the multiple-pole (fractional) part,
  reducing \(\deg_t\) of the denominator to a squarefree remainder;
- a **Rothstein–Trager** residue computation that produces the logarithmic part
  \(\sum c_k \log v_k\) — and, crucially, the *residue criterion* whose failure
  certifies non-elementarity;
- a **polynomial part** whose coefficients are found by recursively solving a
  **Risch differential equation** (RDE) \(Dq + f q = c\) one field-level down.

The RDE is the recursive heart. Bronstein's `Chapters 5–6` organise it into
`SPDE` (Rothstein's degree-reduction, §6.4) and a polynomial non-cancellation
solve (§6.5); Mathilda exposes each as its own builtin — see §4.

### References

1. R. H. Risch, *The problem of integration in finite terms*,
   Transactions of the AMS **139** (1969), 167–189.
2. R. H. Risch, *The solution of the problem of integration in finite terms*,
   Bulletin of the AMS **76** (1970), 605–608.
3. M. Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed.,
   Springer, 2005 — the definitive modern account; Chapters 5 (integration) and
   6 (the Risch differential equation).
4. M. Bronstein, *The transcendental Risch differential equation*,
   Journal of Symbolic Computation **9** (1990), 49–60.
5. J. Liouville, *Mémoire sur l'intégration d'une classe de fonctions
   transcendantes*, Crelle's Journal **13** (1835), 93–118 — the impossibility
   theorem behind the whole method.
6. K. O. Geddes, S. R. Czapor, G. Labahn, *Algorithms for Computer Algebra*,
   Kluwer, 1992 — Ch. 11–12, a readable algorithmic treatment.

---

## 2. Elementary answers, graded

Start the REPL with `./Mathilda` and type each `In[...]` line yourself. Every
transcript below is reproduced from the current build.

### 2.1 The exponential case

The simplest exponential tower is \(\mathbb{C}(x)(e^x)\). Integrating a
polynomial-in-\(e^x\) integrand solves a one-level RDE:

```mathematica
In[1]:= Integrate[x*Exp[x], x, Method -> "RischTranscendental"]
Out[1]= E^x (-1 + x)

In[2]:= Integrate[Exp[x]/(1 + Exp[x]), x, Method -> "RischTranscendental"]
Out[2]= Log[1 + E^x]
```

The second is a logarithmic answer \(\log(1+e^x)\) produced by the
Rothstein–Trager residue step over the exponential monomial.

### 2.2 The "exceptional" exponential — where naive methods fail

Here is the case that defeats table lookup and substitution. When the *inner*
function of an exponential is a genuine polynomial, the antiderivative — if it is
elementary — must have the form \(p(x)\,e^{g(x)}\), and Risch finds \(p\) by
solving \(Dp + g'\,p = (\text{integrand}/e^g)\):

```mathematica
In[3]:= Integrate[2*x*Exp[x^2], x, Method -> "RischTranscendental"]
Out[3]= E^x^2

In[4]:= Integrate[(2*x^2 + 1)*Exp[x^2], x, Method -> "RischTranscendental"]
Out[4]= x E^x^2
```

Notice the delicacy: \(2x\,e^{x^2}\) integrates (the RDE has solution \(p = 1\)),
and \((2x^2+1)e^{x^2}\) integrates to \(x\,e^{x^2}\) (the RDE has solution
\(p = x\)) — but \(e^{x^2}\) *alone* has **no** elementary antiderivative
(§3.2). The algorithm decides which is which by whether the RDE has a
polynomial solution.

### 2.3 The logarithmic case

The dual tower \(\mathbb{C}(x)(\log x)\) yields answers built from \(\log x\):

```mathematica
In[5]:= Integrate[1/(x*Log[x]), x, Method -> "RischTranscendental"]
Out[5]= Log[Log[x]]

In[6]:= Integrate[Log[x], x, Method -> "RischTranscendental"]
Out[6]= -x + x Log[x]
```

### 2.4 A trigonometric integrand is an exponential tower

Mathilda represents \(\tan x\) through the complex exponential \(e^{ix}\), so it
is the transcendental Risch engine — not a trig table — that returns:

```mathematica
In[7]:= Integrate[Tan[x], x, Method -> "RischTranscendental"]
Out[7]= -Log[Cos[x]]
```

---

## 3. The decision procedure: `Risch`ElementaryIntegralQ`

The same machinery answers the *yes/no* question directly.
`Risch`ElementaryIntegralQ[f, x]` returns `True` when it exhibits an elementary
closed form, and `False` when it **proves** none can exist — the constructive
form of Liouville's theorem.

### 3.1 The affirmative side

```mathematica
In[1]:= Risch`ElementaryIntegralQ[x*Exp[x], x]
Out[1]= True

In[2]:= Risch`ElementaryIntegralQ[Tan[x], x]
Out[2]= True

In[3]:= Risch`ElementaryIntegralQ[1/(x*Log[x]), x]
Out[3]= True
```

### 3.2 The negative side — famous impossibility results

This is the payoff. Each `False` below is a *theorem*, first proved (for the
Gaussian) by Liouville in the 1830s and mechanised here by the residue and
Risch-DE no-solution certificates:

```mathematica
In[4]:= Risch`ElementaryIntegralQ[Exp[x^2], x]
Out[4]= False

In[5]:= Risch`ElementaryIntegralQ[Exp[x]/x, x]
Out[5]= False

In[6]:= Risch`ElementaryIntegralQ[1/Log[x], x]
Out[6]= False

In[7]:= Risch`ElementaryIntegralQ[Sin[x]/x, x]
Out[7]= False

In[8]:= Risch`ElementaryIntegralQ[Cos[x]/x, x]
Out[8]= False
```

The Gaussian \(e^{x^2}\), the exponential integrand \(e^x/x\), the reciprocal
logarithm \(1/\log x\), and the sinc function \(\sin x / x\) are the textbook
"non-elementary" integrals. Mathilda does not merely *fail* to integrate them —
it **certifies** that no elementary antiderivative exists. (Each of these *does*
have an antiderivative in terms of a **special function** — \(\mathrm{erf}\),
\(\mathrm{Ei}\), \(\mathrm{li}\), \(\mathrm{Si}\) — which is exactly the subject
of the [next tutorial](cherry-special-functions.md).)

!!! note "Verdicts outside the tower scope"
    `ElementaryIntegralQ` answers `True`/`False` inside the transcendental-tower
    field it can model, and returns *unevaluated* when a verdict would require
    machinery outside that scope (algebraic extensions, or the dilogarithmic
    decision). An unevaluated result is an honest "I cannot decide", never a
    guess.

---

## 4. Building blocks: the Risch differential equation

The recursion bottoms out in the **Risch differential equation** \(Dy + f\,y =
g\). Mathilda exposes the solver and its two reduction stages as first-class
builtins, so you can watch the machine work.

### 4.1 `Risch`RischDE[f, g, x]`

Solves \(y' + f\,y = g\) for a rational \(y\), returning the solution or staying
unevaluated when none exists:

```mathematica
In[1]:= Risch`RischDE[1, 2*x, x]
Out[1]= -2 + 2 x

In[2]:= Risch`RischDE[0, 2*x, x]
Out[2]= x^2

In[3]:= Risch`RischDE[1/x, 1, x]
Out[3]= 1/2 x
```

Check the first: with \(y = 2x - 2\), \(y' + y = 2 + (2x - 2) = 2x\). ✓

### 4.2 `Risch`SPDE` — Rothstein's degree reduction

Rothstein's `SPDE` (Bronstein §6.4) reduces \(Dq + b\,q = c\) with a degree
bound \(n\) to a *smaller* equation, returning
\(\{b', c', m, \alpha, \beta\}\) so that any solution is
\(q = \alpha H + \beta\) with \(\deg H \le m\) and \(DH + b'H = c'\):

```mathematica
In[4]:= Risch`SPDE[1, 2*x, x^2, x, 3]
Out[4]= {2 x, x^2, 3, 1, 0}
```

### 4.3 `Risch`PolyRischDENoCancel` — the polynomial solve

The non-cancellation polynomial RDE (Bronstein §6.5) finds a polynomial \(q\)
with \(\deg q \le n\) solving \(Dq + b\,q = c\) when \(b \ne 0\):

```mathematica
In[5]:= Risch`PolyRischDENoCancel[2*x, 2*x^2 + 1, x, 2]
Out[5]= x
```

Check: with \(q = x\), \(q' + 2x\,q = 1 + 2x^2 = 2x^2 + 1\). ✓ This is exactly
the RDE that decides `Integrate[(2x^2+1) Exp[x^2], x]` in §2.2.

---

## 5. Scope and how it fits the cascade

`"RischTranscendental"` is a **transcendental** decision procedure: its field is
a tower of exponentials and logarithms over \(\mathbb{C}(x)\) (Gaussian constants
\(\mathbb{C}(i)\) are allowed). Integrands requiring **algebraic** extensions
(nested radicals, elliptic integrands) are out of scope and return unevaluated —
those are handled by other cascade stages such as
[`"GoursatAlgebraic"`](../08-calculus.md).

In the automatic dispatch, `Integrate[f, x]` reaches this engine *after* the
faster `"BronsteinRational"` and `"RischNorman"` stages, so you only pay for the
full recursion when the cheaper methods decline. Pin `Method ->
"RischTranscendental"` to force it — useful for confirming *which* engine solved
an integral, and for seeing a clean "unevaluated, because non-elementary"
outcome rather than a special-function answer from a later stage.

---

Next: **[Cherry's special-function extensions](cherry-special-functions.md)** —
what to do when Risch proves an integral is *not* elementary, but it still has a
closed form in \(\mathrm{erf}\), \(\mathrm{Ei}\), \(\mathrm{li}\), or the
dilogarithm.
