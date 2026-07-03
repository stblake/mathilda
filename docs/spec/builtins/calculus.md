# Calculus

Symbolic differentiation is implemented natively in C (`src/calculus/deriv.c`).
Earlier versions of Mathilda bootstrapped `D`, `Dt` and `Derivative`
from a rule file (`src/internal/deriv.m`); that approach was fragile
and slow -- every call walked a linear list of ~60 DownValues, re-ran
pattern matching against each, and recursed through the full rule
engine. The native implementation dispatches on the head symbol in a
single strcmp, uses an allocation-free structural walk for the
`FreeQ[f, x]` test, and builds the derivative tree directly, letting
the ordinary evaluator simplify the arithmetic afterwards. Measured
speedups range from roughly 2x on simple chain-rule calls to ~8x on
mixed-partial and higher-order derivatives; see "Performance Notes"
below.

## D

Partial derivative.
- `D[f, x]` -- derivative of `f` with respect to `x` (treating all
  other symbols as constants).
- `D[f, {x, n}]` -- the `n`-th partial derivative; `n` must be a
  non-negative integer.
- `D[f, x, y, ...]` -- mixed partial derivative, equivalent to
  `D[D[f, x], y, ...]`.
- `D[f, {x, n}, {y, m}, ...]` -- mixed multiple partial derivative.
- `D[f, {{x1, ..., xN}}]` -- gradient: the vector
  `{D[f, x1], ..., D[f, xN]}` for scalar `f`. For a list `f`, gives
  the Jacobian (one gradient row per component of `f`).
- `D[f, {{x1, ..., xN}, n}]` -- the `n`-th array derivative, equal to
  applying `D[f, {{x1, ..., xN}}]` `n` times. `n == 2` yields the
  Hessian (for scalar `f`).
- `D[f, {array1}, {array2}, ...]` -- chained array derivatives;
  semantically `First[Outer[D, {f}, array1, array2, ...]]`. Mixes
  freely with scalar specs.
- `D[f, x, NonConstants -> {y, ...}]` -- treat each listed symbol as
  an implicit function of `x`. Differentiating such a symbol produces
  the canonical unevaluated form `D[y, x, NonConstants -> {y, ...}]`
  rather than `0`, so chain-rule terms surface for implicit
  differentiation. The single-symbol shorthand
  `NonConstants -> y` canonicalises to a singleton List, and the
  option may be combined with multiple scalar/array specs (it applies
  to the same `D` call only; nested `D` expressions in the result
  preserve their own NonConstants list).

**Features**:
- `Protected`, `ReadProtected`.
- Recognises the elementary heads `Plus`, `Times`, `Power`, `Sqrt`,
  `Exp`, `Log`, `Log[b, f]`, all six trig heads and their inverses,
  all six hyperbolic heads and their inverses, and threads
  element-wise over `List`.
- For a single-argument unknown head, applies the standard chain
  rule `D[f[g], x] -> Derivative[1][f][g] * D[g, x]`.
- For multi-argument unknown heads, emits the multi-index
  Derivative form, e.g. `D[f[x, y], x] -> Derivative[1, 0][f][x, y]`.
- Recognises existing `Derivative[n1, ..., nm][f][g1, ..., gn]`
  expressions and advances the appropriate partial-index.
- Short-circuits via an allocation-free FreeQ walk so constants and
  sub-trees independent of the differentiation variable are dropped
  without further work. The short-circuit is suppressed for `Equal[...]`
  expressions and for sub-trees that mention a declared `NonConstants`
  symbol, so both implicit derivatives and relational arguments
  retain their structure.
- Distributes over `Equal`: `D[a == b, x]` becomes
  `Equal[D[a, x], D[b, x]]`. The evaluator folds `Equal[0, 0]` to
  `True`, matching Mathematica.
- Differentiates `Piecewise` clause-wise:
  `D[Piecewise[{{v1, c1}, ...}, d], x]` becomes
  `Piecewise[{{D[v1, x], c1}, ...}, D[d, x]]` — the value expressions
  (and the default) are differentiated while the conditions ride
  through unchanged. Because `Piecewise` is `HoldAll`, each value
  derivative is reduced before being placed back. This makes higher
  derivatives of `UnitStep` stable: `UnitStep'[x]` is
  `Piecewise[{{Indeterminate, x == 0}}, 0]`, and every further
  derivative reproduces the same Piecewise (with
  `D[Indeterminate, x] = Indeterminate`) instead of degrading into a
  `Derivative[1, 0][Piecewise][...]` chain-rule form.

**Examples**:
```mathematica
In[1]:= D[x^3, x]
Out[1]= 3 x^2

In[2]:= D[Sin[x^2], x]
Out[2]= 2 x Cos[x^2]

In[3]:= D[Sin[a x], {x, 3}]
Out[3]= -a^3 Cos[a x]

In[4]:= D[f[g[x]], x]
Out[4]= Derivative[1][f][g[x]] Derivative[1][g][x]

In[5]:= D[f[x, y], x]
Out[5]= Derivative[1, 0][f][x, y]

In[6]:= D[Derivative[2][f][x], x]
Out[6]= Derivative[3][f][x]

In[7]:= D[{x, x^2, Sin[x]}, x]
Out[7]= {1, 2 x, Cos[x]}

In[8]:= D[Log[b, x], x]
Out[8]= 1 / (x Log[b])

In[9]:= D[x^2 + 5 y^3, {{x, y}}]            (* gradient *)
Out[9]= {2 x, 15 y^2}

In[10]:= D[x^2 + 5 y^3, {{x, y}, 2}]        (* Hessian *)
Out[10]= {{2, 0}, {0, 30 y}}

In[11]:= D[{x^2 + y, x y}, {{x, y}}]        (* Jacobian *)
Out[11]= {{2 x, 1}, {y, x}}

In[12]:= D[x^2 + y^2 == 1, x]               (* Equal distribution *)
Out[12]= 2 x == 0

In[13]:= D[x^2 + y^2 == 1, x, NonConstants -> y]   (* implicit diff *)
Out[13]= 2 x + 2 y D[y, x, NonConstants -> {y}] == 0
```

## Dt

Total derivative.
- `Dt[f]` -- total derivative of `f`. Numeric literals and the
  distinguished constants (`Pi`, `E`, `I`, `Infinity`,
  `ComplexInfinity`, `EulerGamma`, `Catalan`, `GoldenRatio`,
  `Degree`) vanish. Unknown symbols `y` appear as `Dt[y]` factors,
  modelling implicit functional dependence.
- `Dt[f, x]`, `Dt[f, {x, n}]`, `Dt[f, x, y, ...]` -- *total* derivative
  with respect to each variable. Unlike `D`, every symbol other than the
  differentiation variable and the distinguished constants is an implicit
  function of that variable, so it contributes a `Dt[s, x]` term rather
  than vanishing. `Dt[x, x]` is `1`; a bare `Dt[s, x]` is its own normal
  form and stays unevaluated.

**Features**:
- `Protected`, `ReadProtected`.
- Shares the elementary-function derivative table with `D`; the
  only dispatch difference is the base-case handling of symbols
  (free symbols become `Dt[s, x]` factors instead of `0`).

**Examples**:
```mathematica
In[1]:= Dt[y^2 + Sin[x]]
Out[1]= 2 y Dt[y] + Cos[x] Dt[x]

In[2]:= Dt[Pi + 3 + x y]
Out[2]= x Dt[y] + y Dt[x]

In[3]:= Dt[x^n, x]
Out[3]= x^(-1 + n) (n + x Log[x] Dt[n, x])

In[4]:= Dt[a x + b, x]
Out[4]= a + x Dt[a, x] + Dt[b, x]

In[5]:= Dt[x^2, {x, 2}]
Out[5]= 2
```

## Derivative

Higher-order derivative operator; represents the symbolic object
`Derivative[n1, ..., nk][f]` (the mixed `(n1, ..., nk)`-th partial
of a `k`-argument function `f`).

**Features**:
- `Protected`, `ReadProtected`.
- Acts primarily as a tag carried through the differentiation
  pipeline: `D` and `Dt` produce `Derivative[...]` heads for
  unknown functions and advance their indices when differentiating
  an expression that already contains one.
- **Operator-form reduction**: `Derivative[n1, ..., nm][f]` is reduced
  at evaluation time when `f` has DownValues that match. The evaluator
  synthesises `Function[{t1, ..., tm}, f[t1, ..., tm]]` (with the
  DownValue rule expanded into the body) and differentiates that pure
  function via the existing pure-function pipeline. This makes
  `f'[x]` (i.e. `Derivative[1][f][x]`) compute the derivative of a
  user-defined `f`. When `f` has no matching DownValue the form is
  left unevaluated, matching Mathematica.

**Examples**:
```mathematica
In[1]:= Derivative[2][f][x]
Out[1]= Derivative[2][f][x]

In[2]:= D[%, x]
Out[2]= Derivative[3][f][x]

In[3]:= D[f[x, y], y]
Out[3]= Derivative[0, 1][f][x, y]

In[4]:= f[x_] := x^5 + 6 x^3
In[5]:= f'[x]
Out[5]= 18 x^2 + 5 x^4

In[6]:= f'[5]
Out[6]= 3575

In[7]:= g[x_, y_] := x^2 y^3
In[8]:= Derivative[1, 1][g][a, b]
Out[8]= 6 a b^2
```

## Integrate (rational-function integration, Phase 1-8d)

`Integrate[f, x]` is the public entry point for the rational-function
integrator implemented in `src/calculus/integrate.c` (System dispatcher) and
`src/calculus/intrat.c` (algorithm package).  Phase 1 of the
`IntegrateRational.m` port (see `plans/INTEGRATE_PLAN.md`) closes the
following classes of integrand:

- **Polynomials in `x`** — term-by-term integration via
  `Integrate`IntegratePolynomial`: `a x^n -> a x^(n+1)/(n+1)` for
  `n != -1`, `a/x -> a Log[x]`.
- **Improper rational functions** — split into polynomial part +
  proper rational via `PolynomialQuotientRemainder`.
- **Repeated roots** — Mack's linear Hermite reduction
  (`Integrate`HermiteReduce`) extracts the rational part `g` such
  that `f = D[g, x] + h` with `Denominator[h]` squarefree
  (Bronstein, *Symbolic Integration I*, p. 44).
- **Derivative-recognition fast path** — when the residual `h` has
  the form `c * D'/D^k` with `c` free of `x` and `k >= 1`, emits the
  closed form `c Log[D]` (k=1) or `-c/((k-1) D^(k-1))` (k>=2).
- **Per-summand Apart loop** — Phase 5 splits the squarefree-
  denominator residual via `Apart` and dispatches each summand
  independently through the derivative-recognition / LRT /
  LogToReal stack, after `ExtractConstants` factors out scalar
  prefactors.  The integral closes only when every piece closes.
- **Lazard-Rioboo-Trager log part with bounded-Solve closure** —
  Phase 2 / 4 run `Integrate`IntRationalLogPart` on the squarefree-
  denominator residual; every resulting `(Q_i, S_i)` pair is then
  dispatched through `Integrate`LogToReal`.  Closure scope:
  - Linear factor `t - α`: contributes `α Log[S_i(α, x)]`.
  - Quadratic factor with positive discriminant: two real roots,
    two `Log` terms.
  - Quadratic factor with negative discriminant (the ArcTan family):
    complex conjugate pair `u ± I v`, contributes
    `u Log[A^2 + B^2] + v LogToAtan[A, B, x]` (Bronstein, *Symbolic
    Integration I*, p. 63 — Rioboo's complex-to-real conversion).
  - Higher-degree irreducible-over-Q factors fall through to
    `Integrate`NaiveLogPart` (Phase 8b/c) — never to a
    speculative-extension retry.  Algebraic-extension closure is
    deferred to `solve.c` / `ToRadicals`.
- **Phase 8b — NaiveLogPart RootSum fallback** —
  `Integrate`NaiveLogPart[f, x]` returns the held-symbolic
  `RootSum[Function[t, d(t)], Function[t, a(t) Log(x - t) / d'(t)]]`
  representation of the log part, mirroring the Lagrange / Bronstein
  closed form `int a/d dx = sum_α a(α) Log(x - α) / d'(α)`.  This is
  universal — every proper rational integrand admits this form, with
  derivatives flowing through the body via the `D[RootSum, x]` rule
  in `src/calculus/deriv.c`.  See `src/root.c` for the held `Root` and
  `RootSum` constructs (HoldAll + Protected).  Direct port of
  `IntegrateRational.m:1116-1124`.
- **Phase 8c — NaiveLogPart wired as universal LogToReal fallback** —
  the closure preference becomes
    1. `IntRationalLogPart -> LogToReal` (real elementary form),
    2. linear-Q closer (`Log + ArcTanh` combination),
    3. `NaiveLogPart` (held `RootSum` form),
  with step 3 universal so `try_lrt_close` never returns NULL for a
  well-formed proper rational input.
- **Phase 8d-bonus — radical-form root expansion** —
  the tail of `NaiveLogPart` runs `expand_simple_rootsum`: when
  `d(t)` is one of the polynomial shapes whose roots have an explicit
  radical-formula closed form (linear / quadratic / biquadratic), the
  `RootSum` is expanded in place to a `Plus` of `body` evaluated at
  each root via the quadratic formula and `Sqrt`.  Higher-degree
  factors stay in held `RootSum` form until `solve.c` is implemented.
- **Stability fixes surfaced by the corpus runs** — the two corpora
  exercised numeric and recursive paths the unit tests had never
  reached, and ran into two distinct child-process crash classes:
  - **Plus/Times BigInt-Rational** (`src/plus.c`, `src/times.c`):
    `add_numbers` and `multiply_numbers` had fast paths only for
    `Rational[Integer, Integer]` and returned `NULL` on
    `Rational[BigInt, ...]` operands (which arose from the
    intermediate resultant computations).  The callers in
    `builtin_plus` / `builtin_times` then dereferenced the NULL
    via `is_overflow()`.  Both helpers now fall through to a
    generic GMP rational add/multiply that recognises any
    combination of `Integer / BigInt / Rational[Integer-or-BigInt,
    Integer-or-BigInt]`, and the callers defensively re-stash the
    operand on a NULL return.  Locked in by
    `tests/test_bigint.c::test_plus_rational_with_bigint_parts`.
  - **`is_zero_poly` recursion bound** (`src/poly/poly.c`):
    `is_zero_poly`'s deep path recurses on each coefficient
    returned by `CoefficientList(expanded, vars[0])`, which
    normally strips one variable per descent.  When `vars[0]` is
    an algebraic constant like `Sqrt[5]` and the polynomial mixes
    several radicals (`Sqrt[5]`, `Sqrt[21]`, `Sqrt[105]`, …),
    `CoefficientList` does not actually strip the variable and
    the recursion sees the same expression again, overflowing the
    C stack.  Threading a `depth` argument and bailing out
    conservatively (returning `false`) at depth 32 keeps the
    function correct on every genuine polynomial while preventing
    SIGSEGV on opaque algebraic-coefficient inputs.

Combined effect on the `IntegrateRationalTests.m` corpus: the
remaining child crashes from earlier runs are eliminated, while
`diff_nonzero` stays at the documented baseline of 1.
- **Phase 6 LogToArcTanh post-processing** — pairs of
  `c Log[A] + c Log[B]` collapse to `c Log[A B]`; sign-paired
  `c Log[A] - c Log[B]` go to `c Log[A/B]` or
  `2 c ArcTanh[(A + B)/(B - A)]` when the simplified argument is
  rational in x.  Beautifies the output without affecting
  differentiation.
- **Phase 7 ArcTan/ArcTanh sign normalisation and option parsing** —
  the final pass strips a leading minus sign from `ArcTan` and
  `ArcTanh` arguments (`ArcTan[-arg] -> -ArcTan[arg]`).  The package
  entry `Integrate`BronsteinRational[f, x, opts...]` accepts trailing
  `Rule` options:
  - `"PFD" -> True | False` (default True) — toggles the per-summand
    Apart loop;
  - `"LogToArcTan" -> True | False` (default True) — toggles the
    Phase 6 post-processing;
  - `"Radicals" -> True | False` (default False) — reserved for
    future `ToRadicals` integration;
  - `Extension -> alpha` — reserved for future use; currently
    advisory.
  Options are stripped before dispatch — Phase 7 keeps them advisory
  while the algorithmic path uses the Mathematica defaults.
- **Post-method normalisation of radical antiderivatives
  (`intsimp_finalize`)** — closed antiderivatives carrying a
  `Power[f(x), p/q]` of the integration variable pass through one final
  cleanup at the `builtin_integrate` chokepoint
  (`src/calculus/intsimp.c`); non-radical and partially-closed results
  are left exactly as the method produced them. For radical results:
  x-free reciprocal/product powers over a positive base are flattened
  (`(1/a)^(2/3) -> a^(-2/3)`, branch-safe because only x-free,
  positive-base subexpressions are `PowerExpand`-ed); the algebraic part
  is recombined with its `Numerator` / `Denominator` `Expand`-ed so
  integer-power compounds collapse in x
  (`6 a^4 - 12 a^3 (a+b x) + 6 a^2 (a+b x)^2 -> 6 a^2 b^2 x^2`, gated to
  fire only when it strictly shrinks); and inverse-trig arguments are
  distributed and sign-normalised. This is why
  `Integrate[1/(x^3 (a+b x)^(1/3)), x]` returns the compact
  `(4 (a+b x)^(5/3) - 7 a (a+b x)^(2/3))/(6 a^2 x^2) + ...` form.

Inputs whose denominator is real-irreducible-over-Q with degree > 4
or non-biquadratic degree 4 close in held `RootSum` form rather than
unevaluated — see Phase 8b / 8c above.  Phase 8a-bonus emits
`Integrate::inexact` and bubbles back unevaluated when the integrand
contains a real-valued constant that cannot be coerced to an exact
rational by `Rationalize` (e.g. `N[Pi]`); `4.5` and `3.125`-style
exact-rational floats are silently coerced and processed.

**Inexact integrand handling**:
```mathematica
In[?]:= Integrate[1/(4.5 x^4 - 3.125), x]
        (* Rationalize coerces 4.5 -> 9/2, 3.125 -> 25/8, then
           closes via the standard pipeline. *)

In[?]:= Integrate[1/(4.5 x^4 - N[Pi]), x]
        Integrate::inexact: Integrand contains inexact numbers...
Out[?]= Integrate[1/(-3.14159 + 4.5 x^4), x]  (* unevaluated *)
```

**Comprehensive corpus runner** (Phase 8d/8e):
`tests/test_intrat_corpus.c` loads any `IntegrateRational` test
corpus file at runtime via Mathilda's own `Get[]`, runs every
`{integrand, var, ...}` entry through the rational integrator using
fork-per-case isolation, and classifies the result via differential
check.  Two CTest entries are wired:
- `intrat_corpus_tests` — RUBI corpus `IntegrateRationalTests.m`
  (113 reference cases with optimal antiderivatives).
- `intrat_inline_corpus_tests` — `IntegrateRationalInlineCases.m`
  generated from the inline `(*IntegrateRational[...]*)` test cells
  in `IntegrateRational.m` via `tools/extract_inline_cases.py`
  (96 cases).

Per-case strict 10 s wall-clock cap (child SIGALRMs at 10 s; the
parent SIGKILLs at the same deadline).  The corpus runner is the
public progress dashboard for the rational integrator: each phase
that lands new closure machinery moves cases out of TIMEOUT /
RootSum into `diff_zero` (closed in real elementary form, with
`D[result, x] - integrand` reducing to 0).  The
`CORPUS_DIFF_NONZERO_BASELINE` constant in the test source is the
high-water mark of known-broken cases — improvements drive it
monotonically down.

**Features**:
- `Protected`, `Listable`.
- Eleven-stage dispatch cascade (`DerivativeDivides`, `LinearRadicals`,
  `QuadraticRadicals` and `LinearRatioRadicals` added 2026-06-06; `Weierstrass`
  added 2026-06-09; `ChebychevAlgebraic` and `GoursatAlgebraic` added 2026-06-29):
  `Integrate[f, x]` (Method -> Automatic, default) tries each subroutine in
  order and returns the first non-`NULL` result:
  1. `Integrate\`Undefined[f, x]` — when `f` contains an undefined-function
     derivative (e.g. `f'[x]`); see below.
  2. `Integrate\`BronsteinRational[f, x]` — when `PolynomialQ[f, x] ||
     rationalQ[f, x]`.
  3. `Integrate\`LinearRadicals[f, x]` — rational functions of `x` and radicals
     `(a x + b)^(m/n)` of one shared linear argument; rationalised by
     `u = (a x + b)^(1/n)`.
  4. `Integrate\`QuadraticRadicals[f, x]` — rational functions of `x` and square
     roots `(a x^2 + b x + c)^(m/2)` of one shared quadratic argument;
     rationalised by a single real-valued Euler substitution.
  5. `Integrate\`LinearRatioRadicals[f, x]` — rational functions of `x` and
     radicals `((a x + b)/(c x + d))^(m/n)` of one shared linear-fractional
     argument; rationalised by `u = ((a x + b)/(c x + d))^(1/n)`.
  6. `Integrate\`ChebychevAlgebraic[f, x]` — Chebychev binomial differentials
     `x^p (a x^r + b)^q` (`p, q, r` rational, `a, b` free of `x`).  Elementary
     iff one of `q`, `(p+1)/r`, `q+(p+1)/r` is an integer (Chebychev's theorem),
     with substitutions `x = u^N` (Type I), `u^s = a x^r + b` (Type II), or
     `u = x^r` then `t^s = (a u + b)/u` (Type III) that rationalise `f`.
     Recognition is a single structural scan, so it runs ahead of
     `DerivativeDivides`'s Eliminate/Solve search.  Non-elementary binomials
     return `NULL` (the cascade falls through to later methods).
  7. `Integrate\`GoursatAlgebraic[f, x]` — pseudo-elliptic integrands
     `F(x) R(x)^q` (`F` rational, `R` a polynomial) with `q` any rational of
     reduced denominator `2`, `3`, or `4` by Goursat's algorithm and its
     cube-/fourth-root generalisations (Blake 2026).  The exponent is split
     `R^q = R^k R^(-p)` with radical order `p in {1/2, 1/3, 2/3, 1/4, 3/4}` and
     the integer `R^k` absorbed into `F`, so positive-power radicals such as
     `(1-x^3)^(1/3)/x` are handled, not only radicals already in a denominator.
     A Mobius automorphism cycling the roots of `R` splits the integrand into
     eigencomponents that descend to genus-0 curves when the elementarity
     criterion holds (`p=1/2`: Klein four-group `V4`, trivial projection
     vanishes; `p=1/3,2/3`: order-3 cycle; `p=1/4,3/4`: order-4 cycle on
     harmonic roots).  For `p=1/2` with `R` a cubic carrying the `t^3-1` higher
     symmetry, when `V4` declines a Section-4 (Goursat 1887) period-3 reduction
     is tried: an order-3 Mobius `S` fixes one ramification point and cycles the
     other three, and the integral is elementary when `F` is a non-trivial
     period-3 character `F(S) = Exp[2 Pi I/3] F` (so `(x-1)/((x+2) Sqrt[x^3-1])`
     integrates).  For `p=1/3` with `R` a cubic, when the order-3 eigendescent is
     obstructed (`H1 != 0`) but `F` has a pole at a non-branch point, a
     constructive third-kind logarithmic-part reduction is tried instead: the
     antiderivative is
     `C Sum_{j=0..2} Exp[4 Pi I j/3] Log[R^(1/3) - Exp[2 Pi I j/3] (lead R)^(1/3) x]`,
     so the parametric family
     `(2-(k+1)x)/((1-(k+1)x) (x(1-x)(1-k x))^(1/3))` integrates over `Q(k)`.
     The rational reductions are integrated recursively and
     back-substituted.  Obstructed (genuinely elliptic) integrands,
     non-harmonic quartics, and the cross-character A4 cases (Section 5, e.g.
     `t/((t^3+8) Sqrt[t^3-1])`) return `NULL`.  A differentiate-back guard rejects the
     rare cases where the eigenspace zero-test misfires on deeply nested radical
     roots; the whole attempt runs under a CPU-time budget so a cyclotomic-root
     `R` with an unlucky cofactor (where algebraic-number `Together`/`Cancel`
     blows up) declines rather than hanging the cascade.  Uses
     `Solve[..., Cubics -> True, Quartics -> True]` (the Ferrari quartic solver,
     added 2026-06-29).  The Boolean global `Integrate\`GoursatDebug` (default
     `False`; added 2026-06-30) traces the descent to stderr when set `True`:
     whether the integrand matches the `F(x) R(x)^(-p)` form (and the recognised
     `F`, `R`, `p`), which involution / eigenspace criterion is tested and
     whether it holds (`V4` trivial projection, the order-3/order-4
     ω-eigencomponents, the period-3 trivial projection per fixed point), and the
     differentiate-back verdict.  The flag is latched once at the outermost call,
     so recursive genus-0 reductions share it and indent by depth.
     Before returning (square-root case), the antiderivative is normalised so its
     one radical is the single generator `Sqrt[R]`: any term left carrying `Sqrt[R]`
     *split* across a numerator/denominator as a product of factor roots
     (`Sqrt[x] Sqrt[(1-x)(1-k^2 x)] ...`) has its `x`-dependent half-power factors
     recombined into one radicand and reduced over `R` (to a rational when it is a
     perfect square, else `rational*Sqrt[R]`), removing the spurious branch point and
     keeping a downstream `Simplify` to a single radical generator.  Constant radicals
     are left intact, and the rewrite is kept only if the differentiate-back guard
     still passes.
     A graded battery of worked examples (every exponent `p`, both numerator and
     denominator radicals, and every involution equation, with the negative
     controls that decline) is collected in
     [`GOURSAT_EXERCISES.md`](../../../GOURSAT_EXERCISES.md) and mirrored as the
     `test_graded` ladder in `tests/test_integrate_goursat.c`.
  8. `Integrate\`Weierstrass[f, x]` — rational functions of the trig kernels
     `Sin/Cos/Tan/Cot/Sec/Csc[x]` (or hyperbolic `Sinh/Cosh/.../Csch[x]`) with a
     kernel in a denominator; continuous `Tan[x/2]` / `Tanh[x/2]` substitution
     (Jeffrey & Rich 1994).  Runs ahead of `DerivativeDivides`: it is
     domain-specific, deterministic, correct by construction, and yields a real,
     continuous antiderivative rather than a complex-logarithm form.
  9. `Integrate\`DerivativeDivides[f, x]` — substitution `u(x)`; in the
     cascade the quiet, branch-correct **direct quotient** strategy only.
  10. `Integrate\`RischNorman[f, x]` — Bronstein pmint, all integrands.
  11. `Integrate\`CRCTable[f, x]` — CRC integral table lookup (lazy-loaded
     from `src/internal/CRCMathTablesIntegrals.m` on first call).
  If every stage gives up the call bubbles back unevaluated.
- `Method -> "<name>"` option (3rd argument) bypasses the cascade and
  dispatches strictly to a single subroutine, with no fallback:
  - `"Automatic"` — default cascade above.
  - `"BronsteinRational"` — `Integrate\`BronsteinRational[f, x]`.
  - `"DerivativeDivides"` — `Integrate\`DerivativeDivides[f, x]` (direct **and**
    the more thorough Eliminate/Solve branch search).  The list form
    `Method -> {"DerivativeDivides", "Substitution" -> u}` **pins** the kernel
    `u(x)`: instead of collecting and trialing every `x`-dependent subexpression,
    only that one substitution is attempted (still both strategies, still strict
    — no fallback to other kernels if `u` does not close the integral). E.g.
    `Integrate[Sqrt[x]/(1 + Sqrt[x]), x, Method -> {"DerivativeDivides", "Substitution" -> Sqrt[x]}]`.
  - `"LinearRadicals"` — `Integrate\`LinearRadicals[f, x]`.
  - `"QuadraticRadicals"` — `Integrate\`QuadraticRadicals[f, x]`.
  - `"LinearRatioRadicals"` — `Integrate\`LinearRatioRadicals[f, x]`.
  - `"ChebychevAlgebraic"` — `Integrate\`ChebychevAlgebraic[f, x]` (Chebychev
    binomial differential `x^p (a x^r + b)^q`).
  - `"GoursatAlgebraic"` — `Integrate\`GoursatAlgebraic[f, x]` (pseudo-elliptic
    `F/R^p`, `p` in `{1/2, 1/3, 2/3, 1/4, 3/4}`, via Mobius eigendescent).
  - `"Weierstrass"` — `Integrate\`Weierstrass[f, x]` (no denominator gate: applies
    to any rational function of the trig/hyperbolic kernels of `x`, including
    polynomial trig).
  - `"RischNorman"` — `Integrate\`RischNorman[f, x]`.
  - `"CRCTable"` — `Integrate\`CRCTable[f, x]`.
  - `"Undefined"` — `Integrate\`Undefined[f, x]`.
  Unknown method names emit `Integrate::method` and bubble back.
- Universal correctness predicate: `Cancel[Together[D[Integrate[f,x],x] - f]] === 0`.

**Examples**:
```mathematica
In[1]:= Integrate[3 + 5 x + 2 x^2, x]
Out[1]= 3 x + 5/2 x^2 + 2/3 x^3

In[2]:= Integrate[2 x/(x^2 + 1), x]
Out[2]= Log[1 + x^2]

In[3]:= Integrate[1/(x - a)^2, x]
Out[3]= -1/(-a + x)

In[4]:= Integrate[(2x+3)/(x^2+3x+5)^2, x]
Out[4]= -1/(5 + 3 x + x^2)

In[5]:= Integrate[1/((x-1)(x-2)(x-3)), x]          (* Phase 2 LRT closes this *)
Out[5]= 1/2 Log[-3 + x] - Log[-2 + x] + 1/2 Log[-1 + x]

In[6]:= Integrate[1/(x^2 + 1), x]                  (* Phase 4 LogToReal *)
Out[6]= ArcTan[x]

In[7]:= Integrate[1/(x^4 + x^2 + 1), x]            (* two quadratic factors *)
Out[7]= 1/6 Sqrt[3] ArcTan[(-1 + 2 x)/Sqrt[3]] +
        1/6 Sqrt[3] ArcTan[(1 + 2 x)/Sqrt[3]] +
        1/4 Log[1 + x + x^2] - 1/4 Log[1 - x + x^2]

In[8]:= Integrate[Sin[x], x, Method -> "RischNorman"]  (* strict, no fallback *)
Out[8]= -Cos[x]

In[9]:= Integrate[x^3, x, Method -> "BronsteinRational"]
Out[9]= 1/4 x^4
```

The `Integrate`` package also exposes the lower-level helpers
`Integrate`HermiteReduce`, `Integrate`IntegratePolynomial`,
`Integrate`BronsteinRational` (the explicit form),
`Integrate`IntRationalLogPart` (Phase 2's LRT computation),
`Integrate`RischNorman` (Bronstein pmint), `Integrate`LinearRadicals`
(linear-radical substitution), `Integrate`QuadraticRadicals`
(quadratic-radical Euler substitution), `Integrate`LinearRatioRadicals`
(linear-fractional / Möbius radical substitution), `Integrate`Weierstrass`
(continuous `Tan[x/2]` / `Tanh[x/2]` substitution), `Integrate`CRCTable`
(table lookup), `Integrate`Undefined` (unknown-function integration,
Roach §1.7), and the unit-test helpers `Integrate`Helpers`Content`,
`...`Primitive`, `...`Monic`, `...`LeadingCoefficient`,
`...`SquareFree`, `...`ExtractConstants`, `...`ApartList`.  All are
`Protected`; the BronsteinRational helpers additionally have
`ReadProtected`.

### Integrate`CRCTable

`Integrate`CRCTable[f, x]` looks `f` up in the CRC Standard Mathematical
Tables (31st ed., 600+ formulas).  The rules live in
`src/internal/CRCMathTablesIntegrals.m`, internally on the head
`IntegrateTable`; `Integrate`CRCTable` is a thin wrapper around it.
The .m file is `Get`-loaded on the first invocation of the CRCTable
stage rather than at startup, so sessions that never call `Integrate`
pay nothing for the table.

Every recursive rule in the table carries an `IntegerQ[index] && index
> base` (or `< base`) guard sufficient to guarantee termination via
first-principles analysis of the reduction direction.  Without these
guards rules such as Formula 49 (the `1/(x^2 - c^2)^n` reduction)
would diverge on negative or non-integer `n`.  As defence-in-depth,
the C dispatcher caps CRC-rule recursion depth at 256 levels and
emits `Integrate`CRCTable::depth` rather than locking up on any rule
that escapes the audit.

The table currently fires on only a small subset of inputs because
Mathilda's pattern matcher does not yet fully support `/;`-guarded
multi-argument rules; this is a separate issue tracked under the
matcher work.

### Integrate`Undefined

`Integrate`Undefined[f, x]` integrates expressions that are rational in
**undefined functions** `u[x]` and their derivatives, following Kelly
Roach, "Indefinite and Definite Integration" (1992), §1.7 ("Undefined
Functions").  Each undefined function value `u[g]` and its derivative
tower `u'[g], u''[g], …` is treated as a differential-field generator;
the integrand is reduced by recognising integration-by-parts /
total-derivative structure in the top generator.  A single inner call to
the rational integrator over a substituted generator symbol subsumes
Roach's polynomial, fraction, and log parts (so `1/u -> Log[u]`,
`1/u^2 -> -1/u`, and `a/(a^2+u'^2) -> ArcTan` all fall out).  Composite
arguments are handled via the chain rule (`g' = D[g, x]`).  A logarithm of
an unknown-function expression, `Log[eta]`, is itself recognised as a
transcendental generator and reduced by parts
(`Integrate[C L + D, x] = G L - Integrate[G L' - D, x]` with
`G = Integrate[C, x]`, `L' = eta'/eta`), with self-referential resolution
for perfect powers (`Integrate[Log[f] f'/f, x] = Log[f]^2/2`).  The stage
is gated to only run when `f` contains an undefined-function derivative or
such a logarithm; genuinely non-elementary integrands (e.g. `f'[x] g'[x]`,
`f'[x]^2`) are left unevaluated, with a cycle guard preventing the by-parts
recursion from looping.  `Protected`, `ReadProtected`.

Known limitations: transcendental generators other than `Log` (e.g.
`ArcTan[eta]`, `Exp[eta]` with `eta` containing an unknown function) are
not yet recognised; nested unknown arguments `f[g[x]]` (with `g` itself
undefined) are deferred.

```mathematica
In[1]:= Integrate[x f'[x] + f[x], x]
Out[1]= x f[x]

In[2]:= Integrate[f'[x] g[x] + f[x] g'[x], x]
Out[2]= f[x] g[x]

In[3]:= Integrate[f'[x]/f[x], x]
Out[3]= Log[f[x]]

In[4]:= Integrate[(f'[x] g'[x] - f[x] g''[x])/(f[x]^2 + g'[x]^2), x]
Out[4]= -ArcTan[Derivative[1][g][x]/f[x]]

In[5]:= Integrate[2 x f'[x^2], x]            (* composite argument *)
Out[5]= f[x^2]

In[6]:= Integrate[(f[x] - x f[x] + f[x] Log[x f[x]] + x f'[x])/f[x], x]
Out[6]= -1/2 x^2 + x Log[x f[x]]            (* Log[eta] generator *)

In[7]:= Integrate[Log[f[x]] f'[x]/f[x], x]
Out[7]= 1/2 Log[f[x]]^2                     (* self-referential by-parts *)
```

### Integrate`DerivativeDivides

`Integrate`DerivativeDivides[f, x]` integrates **by substitution** — the
classical "derivative-divides" method (Moses' SIN, Maxima's `diffdiv`).  It
recognises an integrand of the shape `f(x) = c · h(u(x)) · u'(x)`, reduces the
problem to `Integrate[h[u], u]`, and back-substitutes `u -> u(x)`.
Implemented in `src/calculus/integrate_derivdivides.c`.

Candidate kernels `u(x)` are every distinct subexpression of `f` that depends
on `x`, except `x` and `f` themselves (the C analogue of
`Cases[Union[Level[f, {0, Infinity}]], e_ /; !FreeQ[e, x]]`).  Two
complementary strategies are tried per kernel, each followed by a
**verification gate** `PossibleZeroQ[D[result, x] - f] === True` (numeric
sampling; as of 2026-06-09 the former rigorous `Simplify` confirm was dropped —
it cost ~1.1 s of `trigrat` normalization on the radical-trig cases for no gain
over the sampler):

1. **Direct quotient** — `q = Cancel[Together[f / D[u(x), x]]]`, then
   `q /. u(x) -> u`; accepted when free of `x`.  Cheap, emits no diagnostics,
   handles transcendental compositions (`x Exp[x^2]`), and **selects the
   correct radical branch inherently** (no squaring).  Tried first.

2. **Eliminate / Solve** — builds the differential relation
   `Eliminate[{Dt[y] == f Dt[x], u == u(x), Dt[u == u(x)]}, {x, Dt[x]}]`,
   solves for `Dt[y]`, and reduces each branch with `Factor //@ `,
   `PowerExpand` and `Cancel[. / Dt[u]]`.  It closes integrands the direct
   strategy cannot — including **radical substitutions** like
   `u = Sqrt[Tan[x]]` (reducing `Integrate[Sqrt[Tan[x]], x]` to the rational
   integral `2 u^2/(1 + u^4)`), and cases where `Cancel` canonicalises
   `1/Cos[x]` to `Sec[x]`, breaking the syntactic substitution.  Because the
   algebraisation can square radicals and invert functions, `Solve` returns
   several branches; the verification gate is what **selects the branch that
   differentiates back to `f`**.  This strategy runs in the Automatic cascade as
   well as under the explicit method (its Eliminate `::ifun` / `::alg`
   diagnostics are muted while the integrator drives it).  Because it is
   heavyweight (~0.1–1 s per kernel), as of 2026-06-09 it runs **only on the
   outermost integrand** — reduced sub-integrals are finished by the direct
   strategy and the rest of the cascade.

The reduced integral re-enters the full `Integrate`, so substitutions compose.
Three guards keep the recursion finite and cheap: an **integrand memo** that
short-circuits any integrand (canonicalised by renaming the integration variable
to a fixed sentinel) already attempted in the current top-level descent — this
breaks circular substitution chains and collapses overlapping subproblems that
would otherwise fan out exponentially (e.g. `Integrate[x Sin[x^2], x]`); the
**outermost-only** restriction on the Eliminate/Solve search above; and a hard
depth backstop (8) with per-call fresh substitution symbols.  Strict: returns
unevaluated when no substitution closes the integral.  `Protected`,
`ReadProtected`.

Known limitations: kernels must appear **literally** in `f` (so `Tan[x]`,
which Mathilda keeps atomic rather than `Sin[x]/Cos[x]`, exposes no `Cos[x]`
kernel — such integrands are handled by RischNorman instead); the reduced
integral must itself close under the other methods.

```mathematica
In[1]:= Integrate[Sin[x] Sqrt[1 - Cos[x]], x]      (* direct, correct branch *)
Out[1]= 2/3 (1 - Cos[x])^(3/2)

In[2]:= Integrate[1/(x Log[x]), x]
Out[2]= Log[Log[x]]

In[3]:= Integrate[Sin[x]/Cos[x]^2, x, Method -> "DerivativeDivides"]
Out[3]= Sec[x]                                      (* via Eliminate/Solve *)

In[4]:= Integrate`DerivativeDivides[2 x Exp[x^2], x]
Out[4]= E^x^2

In[5]:= Integrate[Sqrt[Tan[x]], x]                  (* u = Sqrt[Tan[x]] *)
Out[5]= ArcTan[-1 + Sqrt[2] Sqrt[Tan[x]]]/Sqrt[2] + ArcTan[1 + Sqrt[2] Sqrt[Tan[x]]]/Sqrt[2]
          - Log[1 + Tan[x] + Sqrt[2] Sqrt[Tan[x]]]/(2 Sqrt[2])
          + Log[1 + Tan[x] - Sqrt[2] Sqrt[Tan[x]]]/(2 Sqrt[2])

In[6]:= Integrate[Sqrt[Cot[x]], x]                  (* u = Sqrt[Cot[x]] *)
Out[6]= -ArcTan[-1 + Sqrt[2] Sqrt[Cot[x]]]/Sqrt[2] - ArcTan[1 + Sqrt[2] Sqrt[Cot[x]]]/Sqrt[2]
          - Log[1 + Cot[x] - Sqrt[2] Sqrt[Cot[x]]]/(2 Sqrt[2])
          + Log[1 + Cot[x] + Sqrt[2] Sqrt[Cot[x]]]/(2 Sqrt[2])
```

### Integrate`LinearRadicals

`Integrate`LinearRadicals[f, x]` integrates a **rational function of `x` and
radicals `(a x + b)^(m/n)` that share one linear argument** by the classical
rationalising substitution `u = (a x + b)^(1/n)`, where `n = LCM` of the radical
denominators.  Implemented in `src/calculus/integrate_linrad.c`.

It first scans `f` for every `Power[base, p/q]` with `q > 1`: each `base` must be
a degree-1 polynomial in `x` and all must be the **same** linear form `a x + b`
(distinct bases — e.g. `Sqrt[x] + Sqrt[x+1]` — are rejected, as are radicals of a
non-linear argument such as `Sqrt[x^2+1]`).  With `x = (u^n - b)/a` and
`dx = (n/a) u^(n-1) du` the integrand becomes

```
Integrate[f, x] = (n/a) Integrate[ R[(u^n - b)/a, u^M1, u^M2, ...] u^(n-1), u ],
                  Mk = mk n / nk,
```

a **rational function of `u`** that `Integrate`BronsteinRational` closes exactly.
The radical substitution reuses `poly_subst_radical_to_gen` (shared with the
algebraic-factoring path); after the substituted integrand is confirmed rational
in `{x, u}`, the reduced integral re-enters the full `Integrate` and the result
is back-substituted `u -> (a x + b)^(1/n)`.  Unlike the heuristic methods, the
result is **not** put through a `Simplify[D[result, x] - f] === 0` gate: the
substitution is an exact bijection that introduces no branch issues, so the
antiderivative is correct by construction once the rational sub-integral closes.
(Skipping the gate also avoids a prohibitively expensive — and on integrands with
symbolic parameters, non-terminating — `PossibleZeroQ`/`Simplify` over the
symbolic-radical residue.)  A depth guard (8) and per-call fresh substitution
symbols keep the recursion finite and collision-free.  Strict: returns
unevaluated when `f` is not of this form or the reduced integral does not close.
`Protected`, `ReadProtected`.

```mathematica
In[1]:= Integrate[1/Sqrt[x + 1], x]
Out[1]= 2 Sqrt[1 + x]

In[2]:= Integrate[Sqrt[x]/(1 + Sqrt[x]), x]
Out[2]= -2 Sqrt[x] + x + 2 Log[1 + Sqrt[x]]

In[3]:= Integrate[1/(1 + x^(1/3)), x, Method -> "LinearRadicals"]
Out[3]= -3 x^(1/3) + 3/2 x^(2/3) + 3 Log[1 + x^(1/3)]

In[4]:= Integrate[Sqrt[2 x + 3]/x, x, Method -> "LinearRadicals"]
Out[4]= 2 Sqrt[3 + 2 x] - 2 Sqrt[3] ArcTanh[Sqrt[3 + 2 x]/Sqrt[3]]
```

### Integrate`QuadraticRadicals

`Integrate`QuadraticRadicals[f, x]` integrates a **rational function of `x` and
square roots `(a x^2 + b x + c)^(m/2)` that share one quadratic argument** by a
classical **Euler substitution**.  Implemented in
`src/calculus/integrate_quadrad.c`.

It scans `f` for every `Power[base, p/q]` whose `x`-dependent `base` is a
degree-2 polynomial: each must be a **square root** (`q == 2`) and all must be
the **same** radicand `rad = a x^2 + b x + c` (a cube root such as
`(x^2+1)^(1/3)`, a radical of a cubic such as `Sqrt[x^3+1]`, or distinct
radicands are rejected; a purely linear radical belongs to
`Integrate`LinearRadicals`).  Substituting the radicals to a fresh symbol `y`
(`rad^(p/2) -> y^p`, via `poly_subst_radical_to_gen`) must leave a rational
function of `{x, y}`.

To keep the antiderivative **real-valued**, exactly **one** substitution is
chosen by the sign of the leading coefficient `a` (and, when `a < 0`, of the
discriminant `b^2 - 4 a c`) — the routine does *not* try all three Euler forms:

| Condition | Euler substitution | Image `y = Sqrt[rad]` |
|-----------|--------------------|-----------------------|
| `a > 0` | first | `Sqrt[a] x + u` |
| `a < 0` and `b^2 - 4 a c > 0` | third | `(x - alpha) u`, with `alpha` a real root |
| `a` symbolic | first (best-effort `a > 0` branch) | `Sqrt[a] x + u` |

For numeric `a < 0` a real radicand requires real roots, so the third
substitution subsumes the second (the `c > 0` form).  Each branch yields
`x = X(u)`, `dx = X'(u) du` and the radical image above; the rationalised
integrand `f dx /. {y -> ..., x -> X}` re-enters the full `Integrate` and the
result is back-substituted `u -> U(x)`.  As with `Integrate`LinearRadicals`, the
Euler substitution is an exact bijection on the relevant domain, so the result
is **not** put through a `Simplify[D[result, x] - f] === 0` gate — it is correct
by construction once the rational sub-integral closes.  A depth guard (8) and
fresh per-call substitution symbols keep the recursion finite.  Strict: returns
unevaluated when `f` is not of this form or the reduced integral does not close.
`Protected`, `ReadProtected`.

```mathematica
In[1]:= Integrate[1/Sqrt[x^2 + 1], x]
Out[1]= -Log[-x + Sqrt[1 + x^2]]

In[2]:= Integrate[1/Sqrt[1 - x^2], x]
Out[2]= -2 ArcTan[Sqrt[1 - x^2]/(-1 + x)]

In[3]:= Integrate[1/Sqrt[a x^2 + 1], x]               (* symbolic leading coeff *)
Out[3]= -Log[-Sqrt[a] x + Sqrt[1 + a x^2]]/Sqrt[a]

In[4]:= Integrate[1/Sqrt[x^2 - 1], x, Method -> "QuadraticRadicals"]
Out[4]= -Log[-x + Sqrt[-1 + x^2]]
```

### Integrate`LinearRatioRadicals

`Integrate`LinearRatioRadicals[f, x]` integrates a **rational function of `x`
and radicals `((a x + b)/(c x + d))^(m/n)` that share one linear-fractional
(Möbius) argument** by a rationalising substitution.  Implemented in
`src/calculus/integrate_linratiorad.c`.

It scans `f` for every `Power[base, p/q]` (`|q| > 1`) whose `base` depends on
`x`; all such bases must be **structurally identical** and `n = LCM` of the
radical denominators (distinct bases are rejected).  Unlike
`Integrate`LinearRadicals` / `Integrate`QuadraticRadicals`, the base is **not**
required to be a polynomial — it is the ratio
`Times[a x + b, Power[c x + d, -1]]`, which is exactly what partitions this
method from those two (their scans demand a polynomial base, so a ratio radical
falls through to here).  The shared base is canonicalised with
`Cancel[Together[.]]` and the Möbius coefficients read off its numerator
(`a, b`) and denominator (`c, d`); the denominator must be genuinely linear
(degree 1 — a constant denominator is `Integrate`LinearRadicals`' job) and the
map non-degenerate (`a d - b c != 0`).

With `m = n` the substitution `u = ((a x + b)/(c x + d))^(1/m)` inverts the
Möbius map,

```
x  = (d u^m - b)/(a - c u^m),
dx = m (a d - b c) u^(m-1)/(a - c u^m)^2 du,
```

and the radicals rewrite to `u` (`r^(p/q) -> u^(p m/q)`, via
`poly_subst_radical_to_gen`).  The result must be rational in `{x, u}`; the
rationalised integrand re-enters the full `Integrate` and the antiderivative is
back-substituted `u -> ((a x + b)/(c x + d))^(1/m)`.  As with the other radical
stages, the Möbius substitution is an exact bijection on the relevant domain, so
the result is **not** put through a `Simplify[D[result, x] - f] === 0` gate — it
is correct by construction once the rational sub-integral closes.  A depth guard
(8) and fresh per-call substitution symbols keep the recursion finite.  Strict:
returns unevaluated when `f` is not of this form or the reduced integral does
not close.  `Protected`, `ReadProtected`.

```mathematica
In[1]:= Integrate[1/Sqrt[(x + 1)/(x - 1)], x]
Out[1]= 2 Sqrt[(1 + x)/(-1 + x)]/(-1 + (1 + x)/(-1 + x)) -
        2 ArcTanh[Sqrt[(1 + x)/(-1 + x)]]

In[2]:= Integrate[1/Sqrt[(2 x + 1)/(x + 3)], x, Method -> "LinearRatioRadicals"]
Out[2]= -5 Sqrt[(1 + 2 x)/(3 + x)]/(-4 + 2 (1 + 2 x)/(3 + x)) +
        5/2 ArcTanh[Sqrt[(1 + 2 x)/(3 + x)]/Sqrt[2]]/Sqrt[2]
```

### Integrate`Weierstrass

`Integrate`Weierstrass[f, x]` integrates a **rational function of the
trigonometric kernels** `Sin/Cos/Tan/Cot/Sec/Csc[x]` — or the **hyperbolic
kernels** `Sinh/Cosh/Tanh/Coth/Sech/Csch[x]` — by the continuous Weierstrass
substitution of Jeffrey & Rich (*The Evaluation of Trigonometric Integrals
Avoiding Spurious Discontinuities*, ACM TOMS 20(1), 1994).  Added 2026-06-09.

Algorithm: substitute `u = Tan[x/2]` (`u = Tanh[x/2]` for hyperbolic), turning
`f` into a rational function of `u`; integrate that (recursing through
`Integrate`BronsteinRational`); back-substitute; and — for the trigonometric
case — add the secular correction `K Floor[(x - b)/p]` (`b = Pi`, `p = 2 Pi`)
that removes the spurious jump discontinuities the classical substitution
introduces at the poles of `Tan[x/2]` (odd multiples of `Pi`).  The jump `K` is
the difference of the one-sided limits of the `u`-antiderivative at `±Infinity`;
if that limit diverges (a *genuine* singularity of the integrand) no correction
is applied.  A `TrigExpand` pre-pass reduces multiple/sum-angle arguments
(`Cos[2 x]`, `Cosh[x] Cosh[2 x]`, ...) to kernels of the bare variable.

**Hyperbolic case needs no `Floor` correction**: `Tanh[x/2]` is a smooth,
strictly monotone bijection `R -> (-1, 1)` with no poles, so the substitution
introduces no spurious discontinuity — the back-substituted antiderivative is
already continuous (genuine singularities such as `Cosh[x] = 2` are real poles
and are correctly left in place).

The substitution is an exact identity, the rational sub-integral is closed by a
verified integrator, and `Floor' = 0` almost everywhere, so the result is
**correct by construction** — no differentiate-back gate is applied (and the
`Floor` term would defeat symbolic `D` anyway).  In the Automatic cascade only
genuine rational integrands (a kernel in a denominator) are intercepted, so
polynomial trig such as `Integrate[Sin[x], x]` keeps its cleaner table form; the
explicit `Method -> "Weierstrass"` has no such gate.  Strict: returns unevaluated
when `f` is not a rational function of the trig/hyperbolic kernels of `x` (e.g.
`x` outside a kernel, a kernel of a nonlinear argument, mixed trig + hyperbolic,
or a radical of a kernel).  `Protected`, `ReadProtected`.

```mathematica
In[1]:= Integrate[3/(5 - 4 Cos[x]), x]               (* paper eq. (10) *)
Out[1]= 2 Pi Floor[(-Pi + x)/(2 Pi)] + 2 ArcTan[3 Tan[x/2]]

In[2]:= Integrate[1/(2 + Cos[x]), x]                 (* continuous form *)
Out[2]= 2 Pi Floor[(-Pi + x)/(2 Pi)]/Sqrt[3] + 2 ArcTan[Tan[x/2]/Sqrt[3]]/Sqrt[3]

In[3]:= Integrate[1/(2 + Cosh[x]), x]                (* hyperbolic: no Floor *)
Out[3]= 2 ArcTanh[Tanh[x/2]/Sqrt[3]]/Sqrt[3]
```

### InterpolatingFunction integrands

`Integrate[InterpolatingFunction[...][x], x]` returns the antiderivative as a
fresh applied `InterpolatingFunction[...][x]`, mirroring how
`D[InterpolatingFunction[...][x], x]` differentiates such objects.
Differentiation only bumps the object's derivative-order annotation; the
per-window evaluation kernels evaluate derivatives of order `>= 0` only and
cannot produce an antiderivative, so integration instead builds a genuinely new
interpolant (`src/calculus/integrate_interp.c`):

1. Read the grid x-coordinates from the object's stored table.
2. Sample the original function values `y_i = ifun[x_i]`.
3. Accumulate the antiderivative node values `F_0 = 0`,
   `F_i = F_{i-1} + Integrate_{x_{i-1}}^{x_i} ifun` by 5-point Gauss-Legendre
   quadrature (exact through degree 9 — i.e. the default/Spline/Hermite
   piecewise-cubic interpolants; very high explicit `InterpolationOrder` panels
   incur a small quadrature error).
4. Build a Hermite `InterpolatingFunction` through `{{x_i}, F_i, y_i}` — because
   the antiderivative's exact derivative is the original function, supplying
   `F'(x_i) = y_i` makes `D[Integrate[ifun[x], x], x]` round-trip to `ifun[x]`.

Only the 1-D, direct case (the applied argument is the integration variable
itself) is reduced; `Integrate[ifun[g[x]], x]` is not generally expressible as
an `InterpolatingFunction` and is left to the cascade above. The construction
also handles the derivative-annotated objects produced by `D` (integrating the
sampled derivative recovers the lower-order antiderivative). Computations use
machine doubles.

## Integrate`IntRationalLogPart

The Lazard-Rioboo-Trager logarithmic-part computation
(Bronstein, *Symbolic Integration I*, p. 51).

- `Integrate`IntRationalLogPart[A/D, x, t]` returns a list of
  `{Q_i(t), S_i(t, x)}` pairs encoding the integral as
  `Σ_i RootSum[Q_i(t)==0, t Log[S_i(t, x)]]`.
- `Integrate`IntRationalLogPart[A/D, x, t, RootSum -> True]` returns
  the symbolic `Sum[RootSum[Function[t, Q_i], Function[t, t Log[S_i]]]]`
  form directly.

The denominator `D` is assumed squarefree; the typical caller is
`Integrate`BronsteinRational` after Hermite reduction.

```mathematica
In[1]:= Integrate`IntRationalLogPart[1/(x^4+1), x, t]
Out[1]= {{1 + 256 t^4, 4 t + x}}

In[2]:= Integrate`IntRationalLogPart[1/((x-1)(x-2)), x, t]
Out[2]= {{1 - t^2, -3/2 - 1/2 t + x}}

In[3]:= Integrate`IntRationalLogPart[1/(x^2+1), x, t, RootSum -> True]
Out[3]= RootSum[Function[t, 1 + 4 t^2], Function[t, t Log[2 t + x]]]
```

## PolynomialQuotientRemainder

Single-pass companion to `PolynomialQuotient` / `PolynomialRemainder`.

- `PolynomialQuotientRemainder[p, q, x]` returns `{Quotient,
  Remainder}` such that `p == Quotient*q + Remainder` and
  `deg(Remainder, x) < deg(q, x)`.
- Accepts an optional `Extension -> alpha` rule (default `None`) for
  division over `Q(alpha)[x]` rather than the rational coefficient
  field.

```mathematica
In[1]:= PolynomialQuotientRemainder[x^3 + x + 1, x^2 + 1, x]
Out[1]= {x, 1}

In[2]:= PolynomialQuotientRemainder[x^2 - 2, x - Sqrt[2], x, Extension -> Sqrt[2]]
Out[2]= {Sqrt[2] + x, 0}
```

## SubresultantPolynomialRemainders

Polynomial-remainder chain in `K(coeffs)[x]`, used by the
Lazard-Rioboo-Trager log-part computation in the
`Integrate`` package (Phase 2 of `plans/INTEGRATE_PLAN.md`).

- `SubresultantPolynomialRemainders[a, b, x]` gives the chain
  `{a, b, R_2, R_3, ...}` obtained by iterating pseudo-remainder
  until a constant or zero remainder is reached.

The chain is correct modulo content scaling, which downstream
consumers strip via the `primitive[]` operation; this is the
property the LRT algorithm depends on (it consumes the degree of
each chain element and the primitive part of each, both of which are
content-invariant).

```mathematica
In[1]:= SubresultantPolynomialRemainders[x^4 + 1, 2 x^3, x]
Out[1]= {1 + x^4, 2 x^3, 2}
```


## Sum

Definite, indefinite and symbolic summation.  Implemented natively in C
under `src/sum/`: a dispatcher (`sum.c`) plus one file per algorithm,
each algorithm also exposed as a context-qualified builtin
(`Sum`Polynomial`, `Sum`Geometric`, `Sum`Gosper`) — mirroring how
`Integrate` exposes `Integrate`BronsteinRational`.  `Sum` is `HoldAll`.

### Forms

- `Sum[f, {i, imax}]` — sum of `f` for `i` from `1` to `imax`.
- `Sum[f, {i, imin, imax}]` — `i` from `imin` to `imax`.
- `Sum[f, {i, imin, imax, di}]` — step `di`.
- `Sum[f, {i, {i1, i2, ...}}]` — successive values from an explicit list.
- `Sum[f, {i, ...}, {j, ...}, ...]` — multiple (nested) sums; the
  outermost iterator is given first and inner iterators may appear in
  the bounds of outer ones.
- `Sum[f, i]` — the indefinite sum (antidifference): the `F` with
  `DifferenceDelta[F, i] == f`.

For a finite, unit-step, integer range with a summable body, `Sum` first
tries the closed-form method cascade below — its cost is independent of
the span width, so `Sum[i^2, {i, 1, 100000}]` is as cheap as the
symbolic-bound form rather than 100000 term evaluations.  When no closed
form applies (a non-summable body, e.g. `Sum[Prime[i], {i, 1, 6}]`), or
the step is not `1`, or the body iterates an explicit list, the sum falls
back to direct term-by-term expansion (empty ranges fold to `0`).
Symbolic bounds and the indefinite form go straight to the cascade; if
none applies the `Sum[...]` is returned unevaluated.  The cascade has no
step-aware closed form, so a symbolic (or otherwise non-expandable) range
with a non-unit step — e.g. `Sum[i, {i, 1, n, 2}]` — is returned
unevaluated rather than reduced to the (wrong) unit-step result.
`Method -> "Polynomial" | "Geometric" | "Gosper" | "Rational" | "Euler" |
"Alternating"` forces a single algorithm (strict, no fallback), and now also
takes effect on finite unit-step integer ranges.

```mathematica
In[1]:= Sum[i^2, {i, 1, 100}]
Out[1]= 338350
In[2]:= Sum[i^2, {i, 1, n}]
Out[2]= 1/6 n (1 + n) (1 + 2 n)
In[3]:= Sum[f[i, j], {i, 1, 3}, {j, 1, i}]
Out[3]= f[1, 1] + f[2, 1] + f[2, 2] + f[3, 1] + f[3, 2] + f[3, 3]
```

### Sum`Polynomial

Closed-form summation of polynomials via Newton forward differences in
the falling-factorial basis (no Bernoulli numbers).  `Sum`Polynomial[f, i]`
gives the antidifference; `Sum`Polynomial[f, i, imin, imax]` the definite
sum `F(imax+1) - F(imin)`.

```mathematica
In[1]:= Sum[i^3, {i, 1, n}]
Out[1]= 1/4 n^2 (1 + n)^2
In[2]:= Sum[i^2, i]
Out[2]= 1/6 i (-1 + i) (-1 + 2 i)
```

### Sum`Geometric

Summation of `p(i) r^i` with `p` polynomial in `i` and `r` free of `i`
(the ratio of the geometric factors, combined across several such
factors).  A factor is recognised as geometric by the ratio test —
`g(i+1)/g(i)` free of `i` — so every surface form of the same
exponential is caught alike (`r^i`, `r^(-i)`, `(r^i)^(-1)`, `2^(-k)`
held as `(2^k)^(-1)`, `(1/2)^k`, `r^(a i + b)`, …).  The antidifference
has the form `q(i) r^i`; `q` solves `r q(i+1) - q(i) = p(i)` by
undetermined coefficients.

For a finite range the definite sum is `F(imax+1) - F(imin)`.  For an
**infinite** upper limit the sum converges iff `|r| < 1` (a polynomial
times a decaying exponential), in which case the boundary term
`q(n+1) r^(n+1) -> 0` and the sum collapses to `-F(imin)` — an exact
rational.  An undecidable ratio (symbolic `r`) or `|r| >= 1` (including
the divergent `r = -1` edge) is left unevaluated rather than assigned a
false analytic continuation.

```mathematica
In[1]:= Sum[a^i, i]
Out[1]= a^i/(-1 + a)
In[2]:= Sum[q1^i q2^i, i]
Out[2]= (q1 q2)^i/(-1 + q1 q2)
In[3]:= Sum[k^2/2^k, {k, 0, Infinity}]
Out[3]= 6
In[4]:= Sum[k^3/2^k, {k, 0, Infinity}]
Out[4]= 26
```

### Sum`Gosper

Gosper's algorithm for indefinite summation of a hypergeometric term
`t(i)` (one whose ratio `t(i+1)/t(i)` is rational): term-ratio test →
Gosper–Petkovšek normal form (dispersion + gcd peeling) → degree-bounded
key equation `a(i) x(i+1) - b(i-1) x(i) = c(i)` → antidifference
`F = (b(i-1)/c(i)) x(i) t(i)`.  Returns unevaluated when `t` is not a
hypergeometric term or is not Gosper-summable.

```mathematica
In[1]:= Sum[k k!, k]
Out[1]= k!
In[2]:= Sum[k k!, {k, 1, n}]
Out[2]= -1 + (1 + n)!
In[3]:= Sum[1/(i (i + 1)), {i, 1, n}]
Out[3]= 1 - 1/(1 + n)
```

### Sum`Hypergeometric

Infinite sums `Sum[t(k), {k, imin, Infinity}]` (concrete integer `imin`)
whose summand is a hypergeometric term — `t(k+1)/t(k)` rational in `k`. Any
`Binomial` in the summand is first expanded to factorials (whose ratios the
simplifier reduces), so binomial terms such as `2^k/Binomial[2k, k]` are
recognised. The reindexed term ratio is factored into monic linear factors over
the rationals;
the numerator roots give the upper parameters, the denominator roots minus the
canonical `(m+1)` factorial factor give the lower parameters, and the
leading-coefficient ratio is the argument `z`. The result
`t(imin) HypergeometricPFQ[{a}, {b}, z]` is reduced to a closed form by
[`HypergeometricPFQ`](special-functions.md) when possible. Summed only in the
convergent regime (`p<=q` for all `z`; `p==q+1` only for numeric `|z|<1`);
divergent series are left unevaluated rather than returned as the analytic
continuation.

```mathematica
In[1]:= Sum[z^k/k!, {k, 0, Infinity}]
Out[1]= E^z
In[2]:= Sum[x^k, {k, 0, Infinity}]
Out[2]= 1/(1 - x)
In[3]:= Sum[z^k/(2 k)!, {k, 0, Infinity}]
Out[3]= Cosh[2 Sqrt[1/4 z]]
In[4]:= Sum[2^k/Binomial[2 k, k], {k, 1, Infinity}]
Out[4]= 1 + Pi/2
In[5]:= Sum[1/Binomial[2 k, k], {k, 0, Infinity}]
Out[5]= 4/3 + 2 Pi/(9 Sqrt[3])
```

### Sum`Rational

Infinite sums `Sum[p(i)/q(i), {i, imin, Infinity}]` (concrete integer `imin`)
of a rational function of the index, in the convergent regime
`deg q >= deg p + 2`. Decompose into linear partial fractions over the
denominator's splitting field and sum each pole term `c/(i-rho)^k` with the
master identity (`n = i - rho`):

- `k >= 2`: `Zeta[k, imin - rho]` (a Hurwitz Zeta; reduces `1/i^s` to `Zeta[s]`,
  even `s` to powers of `Pi`).
- `k == 1`: `-c (PolyGamma[0, imin - rho] + EulerGamma)` (individually divergent,
  but the residues sum to zero so the combination converges).

Rational poles are found by `Apart[f, i]` over the rationals. An irreducible
quadratic factor with **complex-conjugate** roots (`disc < 0`) is handled in
real partial-fraction form: its constant-numerator (symmetric) part collapses to
`Coth`, while its linear-numerator (antisymmetric) part is the conjugate digamma
sum. A factor with **real radical** roots (`disc >= 0`) is split over the
auto-detected field extension (`Apart[f, i, Extension -> {I, Sqrt[...], ...}]`),
giving `PolyGamma` at radical/complex arguments. Symbolic-coefficient
denominators and finite/indefinite rational sums are out of scope (left
unevaluated). Dispatched in the `Automatic` cascade ahead of
`Sum`Hypergeometric`; `Method -> "Rational"` forces it.

```mathematica
In[1]:= Sum[1/i^2, {i, 1, Infinity}]
Out[1]= 1/6 Pi^2
In[2]:= Sum[1/(i^2 (i^2 + 1)), {i, 1, Infinity}]
Out[2]= 1/6 (3 + Pi^2 - 3 Pi Coth[Pi])
In[3]:= Sum[1/(i (i^2 + 1)), {i, 1, Infinity}]
Out[3]= 1/2 (2 EulerGamma + PolyGamma[0, 1 - I] + PolyGamma[0, 1 + I])
```

### Sum`Alternating

Infinite **alternating rational sums** `Sum[sigma (-1)^k R(k), {k, imin,
Infinity}]` (concrete integer `imin`, `R` rational). The sign factor is
recognised in any form `(-1)^(a k + b)` (`a` odd), with `sigma = (-1)^b`. `R` is
partial-fractioned over the rationals; each linear pole term `c0/(b1 k + b0)^m`
contributes `c0 b1^{-m} (-1)^imin LerchPhi[-1, m, imin - rho]` (`rho = -b0/b1`),
since `Sum_{k>=imin} (-1)^k/(k+a)^m = (-1)^imin LerchPhi[-1, m, imin+a]`. The
Lerch transcendent `Phi(-1, m, .)` reduces to elementary constants — `Log[2]`
and `Pi` at `m = 1` (digamma / Gauss digamma), and `Catalan` / rational·`Pi^m`
at half-integer arguments (Dirichlet beta). Convergence requires `deg q >= deg p
+ 1`. Irreducible-quadratic (complex) poles, radical poles, and divergent inputs
are out of scope and left unevaluated. Dispatched in the `Automatic` cascade
ahead of `Sum`Rational` (so the `(-1)^k` is peeled first); `Method ->
"Alternating"` forces it.

```mathematica
In[1]:= Sum[(-1)^(k + 1)/k, {k, 1, Infinity}]
Out[1]= Log[2]
In[2]:= Sum[(-1)^k/(2 k + 1), {k, 1, Infinity}]
Out[2]= 1/4 (-4 + Pi)
In[3]:= Sum[(-1)^k/(2 k + 1)^2, {k, 0, Infinity}]
Out[3]= Catalan
```

### Sum`Euler

Infinite **linear Euler sums** `Sum[HarmonicNumber[k, p]/k^q, {k, 1, Infinity}]`
(and the order-1 form `HarmonicNumber[k]`), times an optional constant, reduced
to Riemann zeta values. Writing `sigma_h(p, q) = Sum_{k>=1} H_k^{(p)}/k^q`
(converges for `q >= 2`):

- **order 1** (`p = 1`), Euler's formula for any `q >= 2`:
  `Sum H_k/k^q = (1 + q/2) Zeta[q+1] - (1/2) Sum_{j=1}^{q-2} Zeta[j+1] Zeta[q-j]`.
- **diagonal** (`p = q`), from the reflection `sigma_h(p,q) + sigma_h(q,p) =
  Zeta[p] Zeta[q] + Zeta[p+q]`: `Sum H_k^{(p)}/k^p = (Zeta[p]^2 + Zeta[2p])/2`.
- **non-diagonal, odd weight** (`p != q`, `p, q >= 2`, `p + q` odd), via the
  double-zeta reduction `sigma_h(p,q) = Z(q,p) + Zeta[p+q]` with the
  Borwein–Borwein–Girgensohn closed form for `Z(s,t)` (odd `s`), and the
  reflection for even outer `q`. E.g. `Sum H_k^{(2)}/k^3 = 3 Zeta[2] Zeta[3] -
  9/2 Zeta[5]`. Even weight has no zeta reduction and stays unevaluated.
- **quadratic** `Sum H_k^2/k^q` for `q = 2..5` (weight `<= 7`), from the rigorous
  per-weight reduction: `q=2 -> 17/4 Zeta[4]`, `q=3 -> 7/2 Zeta[5] -
  Zeta[2] Zeta[3]`, etc. At weight 8 (`q >= 6`) an irreducible multiple zeta
  value appears and the sum stays unevaluated.

Even zeta values collapse to powers of `Pi` automatically. Higher nonlinear
sums, divergent `q < 1`, and finite/indefinite forms are out of scope and left
unevaluated (never a fabricated value). Dispatched in the `Automatic` cascade
after `Sum`Rational` and before `Sum`Hypergeometric`; `Method -> "Euler"` forces
it.

```mathematica
In[1]:= Sum[HarmonicNumber[k]/k^2, {k, 1, Infinity}]
Out[1]= 2 Zeta[3]
In[2]:= Sum[HarmonicNumber[k, 2]/k^3, {k, 1, Infinity}]
Out[2]= 3 Zeta[2] Zeta[3] - 9/2 Zeta[5]
In[3]:= Sum[HarmonicNumber[k]^2/k^2, {k, 1, Infinity}]
Out[3]= 17 Pi^4/360
```

### Sum`Trigonometric

Infinite **Fourier-type sums** `Sum[T(k)/k^s, {k, imin, Infinity}]` where `T(k)`
is a trigonometric polynomial in `k` (products/powers of `Sin`/`Cos` of
arguments linear in `k`) and `s` is a positive integer. `T` is linearised with
`TrigReduce` into single-angle terms `c_j {Sin|Cos}[a_j k + phi_j]`, each mapped
via `Sum_{k>=1} Sin[a k]/k^s = Im PolyLog[s, e^{i a}]` and `Sum Cos[a k]/k^s =
Re PolyLog[s, e^{i a}]`; constant terms give `c Zeta[s]` (for `s >= 2`). At
`s = 1` with no phase and numeric `0 < a < 2 Pi` the result collapses to
elementary form `Sum Sin[a k]/k = (Pi - a)/2`, `Sum Cos[a k]/k = -Log[2
Sin[a/2]]`. For `s >= 2` (or a phase / symbolic / out-of-range coefficient) the
`Im`/`Re` of a `PolyLog` is returned, exactly as Wolfram leaves it. A divergent
DC term (`s = 1`) leaves the sum unevaluated. Dispatched after `Sum`Alternating`;
`Method -> "Trigonometric"` forces it.

```mathematica
In[1]:= Sum[Sin[k]/k, {k, 1, Infinity}]
Out[1]= 1/2 (-1 + Pi)
In[2]:= Sum[Cos[k]/k, {k, 1, Infinity}]
Out[2]= -(Log[2] + Log[Sin[1/2]])
In[3]:= Sum[Sin[k]/k^2, {k, 1, Infinity}]
Out[3]= Im[PolyLog[2, E^I]]
```

### Sum`LogZeta

Log-weighted zeta series `Sum[c Log[k]/k^s, {k, 1, Infinity}] = -c Zeta'[s]`.
Since there is no symbolic `Zeta'`, an elementary closed form is emitted only for
`s == 2`, via the Glaisher bridge
`-Zeta'[2] = (Pi^2/6)(12 Log[Glaisher] - EulerGamma - Log[2 Pi])`; other `s`
stays unevaluated.  Unblocks `Product[k^(1/k^2)]` through `Product`LogSum`.

### Sum`LogRational

Convergent sums of a rational function plus a log of a rational function,
`Sum[P(k) + Log R(k)]`, where the pieces may individually diverge but combine
(`Sum[1/k + Log[(k-1)/k], {k, 2, Infinity}] = EulerGamma - 1`).  Uses matched
digamma / `LogGamma` asymptotics: with `P` decomposed by `Apart` and `R`'s roots
over `Q`, the value is `sum_{m>=2} c Zeta[m, imin-rho] - sum_{m==1} c PolyGamma[0,
imin-rho] - sum_num e LogGamma[imin-a] + sum_den e LogGamma[imin-b]`, returned
under the log-divergence-cancellation conditions (else unevaluated).

## DifferenceDelta

`DifferenceDelta[f, i]` gives the forward difference
`(f /. i -> i+1) - f`, the discrete analogue of `D` and the left inverse
of indefinite `Sum`.

```mathematica
In[1]:= DifferenceDelta[i^2, i]
Out[1]= 1 + 2 i
In[2]:= DifferenceDelta[Sum[k k!, k], k]
Out[2]= -k! + (1 + k)!
```


## Product

Definite, indefinite, symbolic and convergent-infinite products — the
multiplicative analogue of `Sum`.  Implemented natively in C under
`src/product/`: a dispatcher (`product.c`) plus one file per algorithm,
each exposed as a context-qualified builtin (`Product`Telescoping`,
`Product`Rational`, `Product`Geometric`, `Product`QProduct`,
`Product`Special`, `Product`Infinite`).  `Product` is `HoldAll`.

### Forms

- `Product[f, {i, imax}]` — product of `f` for `i` from `1` to `imax`.
- `Product[f, {i, imin, imax}]` — `i` from `imin` to `imax`.
- `Product[f, {i, imin, imax, di}]` — step `di`.
- `Product[f, {i, {i1, i2, ...}}]` — successive values from a list.
- `Product[f, {i, ...}, {j, ...}, ...]` — multiple (nested) products.
- `Product[f, i]` — the indefinite product (anti-quotient).

A finite, unit-step, integer range with a closed-form body first tries
the method cascade below (cost independent of span width); otherwise the
product expands term-by-term.  **Empty or reversed ranges fold to `1`**
(the multiplicative identity), not `0`.  Symbolic bounds, the indefinite
form, and `imax == Infinity` go to the cascade; if none applies the
`Product[...]` is returned unevaluated.  The `Automatic` cascade runs
Telescoping → Rational → Geometric → QProduct → Special → Infinite, so
the cleanest closed form wins; `Method -> "Telescoping" | "Rational" |
"Geometric" | "QProduct"` forces a single algorithm.  Other options:
`VerifyConvergence` (default `True`; a divergent infinite product prints
`Product::div` and stays unevaluated).  `N[Product[...]]` routes to
`NProduct`.

```mathematica
In[1]:= Product[k, {k, 1, n}]
Out[1]= Factorial[n]
In[2]:= Product[k + a, {k, 1, n}]
Out[2]= Pochhammer[1 + a, n]
In[3]:= Product[2^k, {k, 1, n}]
Out[3]= 2^(1/2 n (1 + n))
In[4]:= Product[1 + 1/k^2, {k, 1, Infinity}]
Out[4]= Sinh[Pi]/Pi
In[5]:= Product[1 - a q^k, {k, 0, n - 1}]
Out[5]= QPochhammer[a, q, n]
In[6]:= Product[k^k, {k, 1, n}]
Out[6]= Hyperfactorial[n]
```

### Product`Telescoping

Rational products whose anti-quotient is itself rational (Gamma-free),
via integer-spaced root chains: `Product[1 - 1/k^2, {k, 2, n}]`
→ `(1 + n)/(2 n)`.  Falls through when the anti-quotient is not rational.

### Product`Rational

The workhorse: any rational `f` whose numerator and denominator factor
into linear factors over `Q`, output in `Pochhammer` / `Factorial`
(`Product[k, {k, 1, n}]` → `n!`).  Falls through on irreducible
quadratic-or-higher factors.

### Product`Geometric

Factors `base^(p(i))` with `base` free of `i`, routed through `Sum`:
`Product[base^p(i)] = base^Sum[p(i), {i, imin, imax}]`, multiplying any
rational cofactor via `Product`Rational`.  The `base^...` factor is never
passed to `Together`/`Factor` (which loop on a symbolic exponent).  An
**infinite** bound is allowed for a pure `base^(summable exponent)` product
(no rational cofactor), the exponent closed by the shipped `Sum`:
`Product[2^(k/2^k), {k, 1, Infinity}]` → `4` (since `Sum[k/2^k] = 2`).

### Product`QProduct

Factors linear in `q^i` (`q` free of `i`) in terms of `QPochhammer`:
`Product[1 - a q^i, {i, imin, imax}] = QPochhammer[a q^imin, q, imax-imin+1]`.

### Product`Special

Named special-function products: `Product[i^i, {i, 1, n}]`
→ `Hyperfactorial[n]`, `Product[Gamma[i], {i, 1, n-1}]` → `BarnesG[n]`.

### Product`Infinite

`imax == Infinity`: a rational convergence gate (equal degree, leading
and next-to-leading coefficients) separates convergent from divergent
products; convergent ones are evaluated as the limit of the finite closed
form, and the Weierstrass family `Product[1 + c/k^2, {k, 1, Infinity}]`
→ `Sinh[Pi Sqrt[c]]/(Pi Sqrt[c])` is recognised directly.  A divergent
product prints `Product::div` and stays unevaluated.

### Product`RationalInfinite

Convergent infinite rational products with **complex-conjugate roots**, via the
Gamma canonical form.  Numerator and denominator are factored over `Q` with
irreducible quadratics resolved by the quadratic formula; under the convergence
conditions (leading constant `1`, balanced degrees, balanced root sums) the
product equals `Product_j Gamma(a - beta_j)^{n_j} / Product_i Gamma(a - alpha_i)^{m_i}`.
The resulting Gamma product is reduced in C (integer shift, cancellation of
equal arguments, conjugate pairs `Gamma[1+ib]Gamma[1-ib] = Pi b/Sinh[Pi b]` and
`Gamma[1/2+ib]Gamma[1/2-ib] = Pi/Cosh[Pi b]`, real reflection) and the candidate
is numerically checked against the raw Gamma product.  Fires only when a genuine
complex root is present (all-real cases stay with `Product`Infinite`); an
irreducible cubic-or-higher factor bails.

```mathematica
In[1]:= Product[(k^2 - 1)/(k^2 + 1), {k, 2, Infinity}]
Out[1]= Pi Csch[Pi]
In[2]:= Product[(k^3 - 1)/(k^3 + 1), {k, 2, Infinity}]
Out[2]= 2/3
```

### Product`Cantor

Double-exponential (Cantor) telescoping `Product[1 + x^(2^i), {i, imin, Infinity}]`
`= 1/(1 - x^(2^imin))` for `|x| < 1` (the exponent must double each step).
`Product[1 + (1/3)^(2^k), {k, 0, Infinity}]` → `3/2`.

### Product`Viete

Viète-type cosine products `Product[Cos[a(i)], {i, imin, Infinity}]` whose angle
halves each step (`a(i+1) == a(i)/2`), giving `Sin[2 a(imin)]/(2 a(imin))`.
`Product[Cos[Pi/2^(k+1)], {k, 1, Infinity}]` → `2/Pi`;
`Product[Cos[x/2^k], {k, 1, Infinity}]` → `Sin[x]/x`.

### Product`EulerPrime

Euler products over the primes: bodies depending on the index only through
`Prime[i]`.  `Product[1/(1 - Prime[i]^-s)]` → `Zeta[s]` (even `s` closes:
`Pi^2/6`, `Pi^4/90`, …), and the `chi_4` product
`Product[1/(1 - (-1)^((Prime[i]-1)/2)/Prime[i]), {i, 2, Infinity}]` → `Pi/4`
(Dirichlet `L(1, chi_4)`).  A divergent (`s <= 1`) product stays unevaluated.

### Product`LogSum

The general Exp/log-sum bridge `Product[f] = Exp[Sum[Log[f]]]` (runs last).
Engages on symbolic-power bodies; `PowerExpand[Log[f]]` is summed by the shipped
`Sum` and re-exponentiated when it closes:
`Product[Exp[(-1)^k/k], {k, 1, Infinity}]` → `1/2`,
`Product[k^(1/k^2), {k, 1, Infinity}]` → `Exp[-Zeta'[2]]` (Glaisher form, via
`Sum`LogZeta`), `Product[E^(1/k)(1 - 1/k), {k, 2, Infinity}]` → `E^(EulerGamma-1)`
(via `Sum`LogRational`).


## FindRoot

Iterative numerical root finder.  Implemented natively in C in
`src/findroot.c`.  Has `HoldAll, Protected` attributes and uses
`Block`-style local binding of the search variables so the user's
global symbol table is not perturbed during iteration.

### Forms

- `FindRoot[f, {x, x0}]` -- Newton from a single start.
- `FindRoot[lhs == rhs, {x, x0}]` -- equation form (normalised to
  `lhs - rhs`).
- `FindRoot[f, {x, x0, x1}]` -- secant from two starts (no derivative
  needed).
- `FindRoot[f, {x, xstart, xmin, xmax}]` -- Brent's method on the
  bracket `[xmin, xmax]`.
- `FindRoot[{f1, ..., fn}, {{x, x0}, {y, y0}, ...}]` -- multivariate
  Newton with a symbolically computed Jacobian, solved per step via
  in-place Gaussian elimination with partial pivoting.

Each `f_i` may be either an expression (treated as `f_i == 0`) or an
explicit equation `lhs_i == rhs_i`.

### Output

`{ var -> value }` for the scalar case, or
`{ var1 -> v1, ..., varN -> vN }` for a system -- the same `Rule`-list
form that `Solve` returns.

### Method dispatch

| Spec form                                  | Default method |
|--------------------------------------------|----------------|
| `{x, x0, xmin, xmax}` (4-element bracket)  | Brent          |
| `{x, x0, x1}` (two start points)           | Secant         |
| `{x, x0}` (single start), scalar           | Newton         |
| `{{x,x0},{y,y0},...}` (system)             | Newton         |

The `Method` option overrides this when explicit (`"Newton"`,
`"Secant"`, or `"Brent"`).  `Method -> "Brent"` accepts either a
4-element bracket spec or a 2-start spec used as `[a, b]`.

When the symbolic derivative `D[f, x]` cannot be reduced to a
numericalizable form — `D` returns an inert `Derivative[..]` head, as for
`Zeta'`, `HurwitzZeta`'s order derivative, `PolyGamma`, etc. — the Newton
iteration falls back to a **central finite-difference** derivative
(`(f(x+h) − f(x−h))/(2h)`) instead of failing with `nlnum`. This covers all
four scalar Newton kernels (real / complex × machine / MPFR), so e.g.
`FindRoot[Zeta[s] == 1.05, {s, 5}]` → `{s -> 4.61297}` and
`FindRoot[2^(-s)(Zeta[s] - Zeta[s, 3/2]) == 0, {s, 20 I}]` →
`{s -> 0.948962 + 20.3778 I}`.

### Complex roots

`FindRoot` automatically engages a complex Newton inner loop when the
starting value evaluates to a `Complex` with non-zero imaginary part.
At machine precision this uses `double complex` arithmetic; at
arbitrary precision the iteration runs on an `(mpfr_t re, mpfr_t im)`
pair and the codebase's existing MPFR-complex transcendental paths
(in `numeric.c`, `trig.c`, `logexp.c`) supply `Sin`, `Exp`, `Zeta`,
etc.

### Options

| Option              | Default        | Effect |
|---------------------|----------------|--------|
| `Method`            | `Automatic`    | `"Newton"`, `"Secant"`, `"Brent"`, or `Automatic`. |
| `WorkingPrecision`  | `MachinePrecision` | `MachinePrecision`, or a digit count (>= ~16 routes through MPFR). |
| `MaxIterations`     | `100`          | Iteration limit. |
| `AccuracyGoal`      | `Automatic`    | Digit count `n` ⇒ stop when `|f| < 10^{-n}`. `Infinity` disables this criterion. `Automatic` resolves to `WorkingPrecision/2`. |
| `PrecisionGoal`     | `Automatic`    | Digit count `n` ⇒ stop when `|step| < |x| * 10^{-n}`.  Same defaults / specials as `AccuracyGoal`. |
| `DampingFactor`     | `1`            | Multiplier on the Newton step.  Useful for repeated roots (set to the multiplicity for quadratic convergence). |
| `Jacobian`          | `Automatic`    | A symbolic expression (scalar form) or `n*n` nested list (system form) that overrides the call to `D[]`. |
| `StepMonitor`       | `None`         | A held expression evaluated after each step. |
| `EvaluationMonitor` | `None`         | A held expression evaluated each time `f` (or any of its derivatives) is evaluated. |

### Diagnostics (stderr)

| Tag         | Triggered when |
|-------------|----------------|
| `FindRoot::argt`    | Wrong arg count. |
| `FindRoot::ivar`    | Malformed variable spec. |
| `FindRoot::nlnum`   | `f` or its derivative did not evaluate to a number. |
| `FindRoot::cvmit`   | `MaxIterations` exhausted (the last iterate is still returned). |
| `FindRoot::noconv`  | Divergence (non-finite step). |
| `FindRoot::brnoth`  | Brent given non-bracketing endpoints. |
| `FindRoot::badopt`  | Unknown option name or invalid value. |
| `FindRoot::badmeth` | Unknown `Method` value. |
| `FindRoot::vecvar`  | Vector-valued variable spec (deferred to a follow-up). |
| `FindRoot::dsing`   | Derivative vanished and `|f|` was not already below tolerance. |

### Examples

```mathematica
In[1]:= FindRoot[Sin[x] + Exp[x], {x, 0}]
Out[1]= {x -> -0.588533}

In[2]:= FindRoot[Cos[x] == x, {x, 0}]
Out[2]= {x -> 0.739085}

In[3]:= FindRoot[{y == Exp[x], x + y == 2}, {{x, 1}, {y, 1}}]
Out[3]= {x -> 0.442854, y -> 1.55715}

In[4]:= FindRoot[Sin[x], {x, 3}, WorkingPrecision -> 50]
Out[4]= {x -> 3.1415926535897932384626433832795028841971693993751}

In[5]:= FindRoot[(Cos[z + I] - 2) (z + 2), {z, 1.0 + 0.1 I}]
Out[5]= {z -> -3.1302*10^-27 + 0.316958 I}

In[6]:= FindRoot[Cos[x] - x, {x, 0, 1}, Method -> "Brent"]
Out[6]= {x -> 0.739085}

In[7]:= FindRoot[x^2 - 2, {x, 1.0, 2.0}, Method -> "Secant"]
Out[7]= {x -> 1.41421}

In[8]:= FindRoot[(x - 1)^3, {x, 0.5}, DampingFactor -> 3]
Out[8]= {x -> 1.0}
```

## FindMinimum / FindMaximum

Iterative local optimisation.  Implemented natively in C in
`src/findmin.c`.  Both have `HoldAll, Protected` attributes and use
`Block`-style local binding of the search variables.  `FindMaximum[f, ...]`
is a thin wrapper around `FindMinimum[-f, ...]` that negates the
objective value before returning.

### Forms

- `FindMinimum[f, {x, x0}]` -- 1D from a single start (default Brent).
- `FindMinimum[f, {x, x0, x1}]` -- bracket Brent on `[x0, x1]`.
- `FindMinimum[f, {x, xstart, xmin, xmax}]` -- bracket Brent on `[xmin, xmax]`.
- `FindMinimum[f, {{x, x0}, {y, y0}, ...}]` -- n-D from a user start.
- `FindMinimum[f, {x, y, ...}]` -- n-D auto-start at 1 for each variable
  (matches Mathematica; avoids the common saddle-at-origin trap for
  oscillatory objectives like `Sin[x] Sin[2 y]` whose gradient vanishes
  at the origin).
- `FindMinimum[{f, cons}, vars]` -- constrained minimisation.
- `FindMaximum[...]` -- same forms; maximises `f` (equivalent to negating
  the objective and the f-value of the result).

### Output

`{f_min, {var1 -> v1, ..., varN -> vN}}` -- a 2-element list whose first
element is the function value and whose second is the rule list for the
optimising variable assignments.

### Method dispatch

| Spec                  | Default method |
|-----------------------|----------------|
| n = 1                 | Brent          |
| n >= 2                | QuasiNewton (BFGS) |
| `{x, x0, x1}` (1D)    | Brent (bracket) |

Methods overridable via `Method -> "Brent" | "Newton" | "QuasiNewton"
| "ConjugateGradient"`.  Brent is golden-section search with parabolic
interpolation (derivative-free), QuasiNewton is BFGS with Armijo
backtracking line search, ConjugateGradient is Polak-Ribière+ with
restart, Newton uses the symbolic Hessian (via repeated `D[]`) with
modified-Cholesky safeguarding.  Gradients default to a symbolic
gradient (`D[f, x_i]` per variable) with a central-difference fallback;
override via `Gradient -> {dfdx1, dfdx2, ...}`.

### Constraints

Inside the `{f, cons}` form, `cons` is a boolean tree of comparisons:

- `Less / LessEqual / Greater / GreaterEqual` between a bare iteration
  variable and a numeric constant become **box constraints** (enforced
  by per-step projection).
- Other inequalities (`g(x) <= 0`) and equalities (`h(x) == 0`) feed a
  **quadratic-penalty** wrapper around the inner solver.  The outer μ
  schedule starts at 1 and multiplies by 10 each round until feasible
  (max 9 rounds, μ up to 10^8).  The inner BFGS/CG/Newton iterations
  drive the *augmented* objective `f + μ·Σ_k max(0, g_k)^2 + μ·Σ_j h_j^2`
  using a matching *augmented* gradient `∇f + 2μ·Σ_k (active) g_k ∇g_k +
  2μ·Σ_j h_j ∇h_j` — the gradient of each constraint expression is
  computed symbolically at setup and falls back to central differences
  per-constraint when symbolic differentiation fails.
- `Or[...]`, `Element[...]`, `x ∈ Integers` and the rest of the
  Mathematica constraint surface are not yet implemented -- they emit
  `FindMinimum::nimpl`.

The penalty wrapper converges from feasible *or* infeasible starting
points on smooth nonlinear constraints (linear/quadratic inequalities,
linear/quadratic equalities, intersections thereof).  At very high μ
the inner solver may exit early on line-search exhaustion; the outer
loop's feasibility check is authoritative, and only emits a diagnostic
(`FindMinimum::infeas`) when the final iterate genuinely fails the
constraint tolerance (1e-12).

### Options

| Option              | Default        | Effect |
|---------------------|----------------|--------|
| `Method`            | `Automatic`    | `"Brent"`, `"QuasiNewton"`, `"ConjugateGradient"`, `"Newton"`, or `Automatic`. |
| `WorkingPrecision`  | `MachinePrecision` | `MachinePrecision`, or a digit count (>= ~16 routes through MPFR).  Lifts the 1D `Brent` and n-D `QuasiNewton` iterations into MPFR at the requested precision so the result `{f_min, {x -> ...}}` carries MPFR leaves with that many digits.  Explicit `Method -> "Newton"` or `"ConjugateGradient"` at MPFR currently falls back to `QuasiNewton` with a `FindMinimum::nimpl` diagnostic; general (non-box) constraints at MPFR fall back to machine precision similarly. |
| `MaxIterations`     | `500`          | Iteration limit on the inner loop. |
| `AccuracyGoal`      | `Automatic`    | Digit count `n` ⇒ stop when `\|grad\| < 10^{-n}`. `Infinity` disables. `Automatic` resolves to `WorkingPrecision/2`. |
| `PrecisionGoal`     | `Automatic`    | Digit count `n` ⇒ stop when `\|step\| < \|x\| * 10^{-n}`. |
| `Gradient`          | `Automatic`    | Explicit `{dfdx1, ..., dfdxN}` overrides the symbolic gradient. |
| `StepMonitor`       | `None`         | A held expression evaluated after each step. |
| `EvaluationMonitor` | `None`         | A held expression evaluated each time `f` (or any partial) is evaluated. |

### Diagnostics (stderr)

| Tag                  | Triggered when |
|----------------------|----------------|
| `FindMinimum::argt`    | Wrong arg count. |
| `FindMinimum::ivar`    | Malformed variable spec. |
| `FindMinimum::vecvar`  | Vector-valued variable (deferred). |
| `FindMinimum::badmeth` | Unknown `Method` value. |
| `FindMinimum::badopt`  | Unknown option name or invalid value. |
| `FindMinimum::nimpl`   | Method, constraint shape, or domain restriction not yet supported. |
| `FindMinimum::nlnum`   | `f`, gradient, or constraint did not evaluate to a number. |
| `FindMinimum::cvmit`   | Inner-loop `MaxIterations` exhausted. |
| `FindMinimum::lstol`   | Line search could not find a sufficient decrease. |
| `FindMinimum::dsing`   | Hessian non-positive-definite (Newton). |
| `FindMinimum::infeas`  | Penalty outer loop could not satisfy all constraints. |

### Examples

```mathematica
In[1]:= FindMinimum[(x - 3)^2 + 1, {x, 0}]
Out[1]= {1.0, {x -> 3.0}}

In[2]:= FindMinimum[x Cos[x], {x, 2}]
Out[2]= {-3.28837, {x -> 3.42562}}

In[3]:= FindMinimum[x Cos[x], {x, 7, 1, 15}]
Out[3]= {-9.47729, {x -> 9.52933}}

In[4]:= FindMinimum[Sin[x] Sin[2 y], {{x, 2}, {y, 2}}]
Out[4]= {-1.0, {x -> 1.5708, y -> 2.35619}}

In[5]:= FindMinimum[{x Cos[x], 1 <= x && x <= 15}, {x, 7}]
Out[5]= {-9.47729, {x -> 9.52933}}

(* Chained `lo <= x <= hi` parses as an Inequality[...] node; FindMinimum
   walks its (value, op, value) triples and registers each as a box bound. *)
In[5b]:= FindMaximum[{x Cos[x], 1 <= x <= 15}, {x, 7}]
Out[5b]= {6.36096, {x -> 6.4373}}

In[6]:= FindMinimum[(1-x)^2 + 100 (y-x^2)^2, {{x, 0}, {y, 0}}]
Out[6]= {0.0, {x -> 1.0, y -> 1.0}}

In[7]:= FindMaximum[Cos[x], {x, 0}]
Out[7]= {1.0, {x -> 0.0}}

In[8]:= FindMinimum[(x - 3)^2, {x, 0}, Method -> "ConjugateGradient"]
Out[8]= {0.0, {x -> 3.0}}

(* Arbitrary precision via WorkingPrecision: the 1D Brent and n-D BFGS
   iterations both run in MPFR at the requested precision and the
   returned `{f_min, {x -> ...}}` carries MPFR leaves with that many
   digits. *)
In[9]:= FindMinimum[(x - Pi)^2, {x, 0}, WorkingPrecision -> 50]
Out[9]= {0.0, {x -> 3.1415926535897932384626433832795028841971693993751}}

In[10]:= FindMinimum[x Cos[x], {x, 2}, WorkingPrecision -> 80]
Out[10]= {-3.2883713955908964865125964571235283975158511553846230554230811211040811736596049,
         {x -> 3.42561845948172814647771386218545617769640923939753965919739613085112431446169}}
```
