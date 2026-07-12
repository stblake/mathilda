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

## Limit

`Limit[f, x -> a]` finds the limiting value of `f` as `x` approaches `a`.
`Limit[f, {x1 -> a1, ...}]` is the iterated (rightmost-first) limit and
`Limit[f, {x1, ...} -> {a1, ...}]` the joint multivariate limit.

**Features**:
- `HoldAll`, `Protected`, `ReadProtected`.
- Options: `Direction -> Automatic`, `Assumptions -> Automatic`,
  `Method -> Automatic`.
- **Direction** selects the approach: `Reals`/`"TwoSided"` (default),
  `"FromAbove"` (or `-1`), `"FromBelow"` (or `+1`), or `Complexes`.
- **Method** restricts the internal strategy cascade to a single named
  group. `Automatic` runs the full cascade in order; a named method runs
  only that group and leaves `Limit` unevaluated if it does not apply
  (an unrecognised name emits `Limit::method` and is likewise left
  unevaluated). The restriction applies only to the outermost call —
  recursive sub-limits always use the full cascade. Methods:
  - `"Substitution"` — continuity / direct substitution, `Abs` kink
    resolution, atom-substitution and one-sided probes.
  - `"RationalFunction"` — leading-degree comparison for `P(x)/Q(x)`.
  - `"Series"` — Taylor / Laurent / Puiseux leading-term expansion.
  - `"LHospital"` — L'Hospital's rule with growth guardrails.
  - `"Asymptotic"` — dominant-term / `Log` / exponential reductions at
    infinity, including `f^g` via `Exp[g Log f]`, and the compose-at-infinity
    rule: for `f[g(x)]` whose inner argument diverges to `±Infinity`, apply the
    builtin's own value at Infinity (`Erf[Infinity] = 1`, `Tanh[Infinity] = 1`,
    `ArcTan[Infinity] = Pi/2`, `Gamma[Infinity] = Infinity`, …). Functions that
    do not self-evaluate there (oscillatory `Sin`, `Cos`) fall through and yield
    `Indeterminate`.
  - `"Bounded"` — squeeze envelope and bounded-oscillation `Interval`. Also
    covers a bounded base raised to a divergent positive power (`exp -> +Infinity`):
    with `B >= |base|` the pointwise magnitude bound, the limit is `0` when
    `Limit[B]` lies in `[0, 1)` (e.g. `(Sin[1/x]/2)^(1/x^2) -> 0` and the
    shrinking-bound `(x Sin[1/x]/2)^(1/x^2) -> 0` at `x -> 0`), and `Infinity`
    when the base is bounded below by a constant `> 1` and positive (e.g.
    `(2 + Sin[1/x]/2)^(1/x^2) -> Infinity`).
- May return a finite value, `Infinity`, `-Infinity`, `ComplexInfinity`,
  `Indeterminate`, an `Interval[{lo, hi}]`, or the original expression
  unevaluated when the limit cannot be determined.

**Examples**:
```mathematica
In[1]:= Limit[Sin[x]/x, x -> 0]
Out[1]= 1

In[2]:= Limit[(x^2 - 1)/(x - 1), x -> 1]
Out[2]= 2

In[3]:= Limit[(1 + 1/x)^x, x -> Infinity]
Out[3]= E

In[4]:= Limit[1/x, x -> 0]
Out[4]= ComplexInfinity

In[5]:= Limit[1/x, x -> 0, Direction -> "FromAbove"]
Out[5]= Infinity

In[6]:= Limit[x^2 + y^2, {x, y} -> {1, 2}]
Out[6]= 5

In[7]:= Limit[Sin[x]/x, x -> 0, Method -> "Series"]
Out[7]= 1

In[8]:= Limit[(2 x^2 + 1)/(x^2 + x), x -> Infinity, Method -> "RationalFunction"]
Out[8]= 2

In[9]:= Limit[Sin[x]/x, x -> 0, Method -> "RationalFunction"]
Out[9]= Limit[Sin[x]/x, x -> 0, Method -> "RationalFunction"]
```

## Residue

`Residue[f, {z, z0}]` gives the residue of `f` at the isolated singularity
`z = z0` — the coefficient of `(z - z0)^-1` in the Laurent expansion of `f`
(`src/calculus/residue.c`). Attribute: `Protected`. The numerical counterpart is
[`NResidue`](../builtins/numerical-calculus.md); the symbolic engine here is
exact but needs `f` to admit a Laurent series at `z0`.

**Method.** The residue is read straight out of `Series[f, {z, z0, 0}]`: an
order-0 expansion spans exponent `-1` for most integrands, and its `-1`
coefficient is the residue. When the pole's leading behaviour is carried by an
unknown function (e.g. `f[z]/z^5`), the series engine truncates relative to that
function's depth, so `Residue` raises the expansion order until the `-1`
coefficient is among the explicit terms (bounded, to stay safe at essential
singularities). A fractional-power (Puiseux) expansion, `den > 1`, signals a
**branch point**, where the residue is undefined — the call is left unevaluated.

For a **rational** integrand the expansion is taken about `z0` with an *expanded*
denominator (`z -> z0 + w`, then `Expand`): this collapses the radical arithmetic
in the denominator's constant term (`Sqrt[3]^2 -> 3`, …) so a pole whose location
is a **sum of radicals** — e.g. `z0 = -2 + Sqrt[3]`, a root of `1 + 4 z + z^2` —
is detected (a naïve `Series` would evaluate `Denominator(z0)` to a non-simplified
nonzero form and wrongly report residue `0`). Transcendental / special-function
integrands keep the direct expansion, which uses the engine's own knowledge of
the Laurent series at `z0` (e.g. `Zeta` at `1`).

```
In[1]:= Residue[1/z, {z, 0}]
Out[1]= 1

In[2]:= Residue[1/z^2, {z, 0}]
Out[2]= 0

In[3]:= Residue[1/Sin[z]^5, {z, 0}]
Out[3]= 3/8

In[4]:= Residue[(z + 1)/(z^2 (z - 2)), {z, 0}]      (* order-2 pole *)
Out[4]= -3/4

In[5]:= Residue[1/(z^2 + 1), {z, I}]                (* complex pole *)
Out[5]= -I/2

In[6]:= Residue[x^3/(x^4 - 2), {x, 2^(1/4)}]        (* algebraic pole *)
Out[6]= 1/4

In[7]:= Residue[f[z]/z^5, {z, 0}]                   (* unknown numerator *)
Out[7]= 1/24 Derivative[4][f][0]

In[8]:= Residue[Zeta[z]/(z - 1)^10, {z, 1}]
Out[8]= -StieltjesGamma[9]/362880

In[9]:= Residue[1/Sqrt[z], {z, 0}]                  (* branch point *)
Out[9]= Residue[1/Sqrt[z], {z, 0}]
```

Consistency with Cauchy's theorem: for a contour enclosing only the pole,
`NIntegrate[f, {z, ...loop...}]/(2 Pi I)` agrees with `Residue`; e.g.
`Residue[1/Sin[z]^7, {z, 0}]` is `5/16`, matching `NResidue[1/Sin[z]^7, {z, 0}]`
≈ `0.3125`.

## Integrate (rational-function integration, Phase 1-8d)

`Integrate[f, x]` is the public entry point for the rational-function
integrator implemented in `src/calculus/integrate.c` (System dispatcher) and
`src/calculus/intrat.c` (algorithm package).  Phase 1 of the
`IntegrateRational.m` port (see `plans/INTEGRATE_PLAN.md`) closes the
following classes of integrand:

A **list integrand** threads element-wise: `Integrate[{f1, ..., fn}, spec...]`
returns `{Integrate[f1, spec...], ..., Integrate[fn, spec...]}` for both the
indefinite `x` and the definite `{x, a, b}` / contour spec forms.  (`Integrate`
is deliberately not `Listable`, which would wrongly also thread over the range
spec; the integrand-only threading is handled explicitly.)

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
  10. `Integrate\`RischNorman[f, x]` — Bronstein pmint (parallel Risch), all
     integrands.
  11. `Integrate\`RischTranscendental[f, x]` — the **recursive** transcendental
     Risch algorithm; runs after RischNorman and only adds
     closed forms the earlier stages missed.  Correct by construction (no
     differentiation check).  Handles logarithmic polynomials and the
     special-function cases below (Erf, ExpIntegralEi, LogIntegral, PolyLog).
     Its Risch differential equation is solved by Bronstein's rational one-step
     (SPDE) reduction (polynomial-gcd time, no undetermined-coefficient
     blow-up), closing high-degree `R(x) e^x` forms and — via the
     exponential special-denominator Laurent ansatz and an exact
     tower-variable verification — nested exp/log towers such as
     `∫((e^x−x²+2x)/(x²(x+e^x)²)) e^((x²−1)/x + 1/(x+e^x)) dx = e^(−x+1/(e^x+x)+(x²−1)/x)`
     and `∫(1/(x log(1+e^x)) − …) dx = log(x)/log(1+e^x)`.  Arithmetic warnings
     from transient internal singular expressions are muted (as in Mathematica).
  12. `Integrate\`CRCTable[f, x]` — CRC integral table lookup (lazy-loaded
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
  - `"RischNorman"` — `Integrate\`RischNorman[f, x]` (parallel Risch / pmint).
  - `"RischTranscendental"` — `Integrate\`RischTranscendental[f, x]`, the recursive
    transcendental Risch algorithm (`src/calculus/integrate_risch_transcendental.c`).
    A decision procedure over a differential transcendental tower, distinct
    from the parallel-Risch heuristic `"RischNorman"`.  Every case is correct
    by construction — it fires only behind an exact structural certificate, so
    the result is not checked by differentiation.  Cases:
      - rational: delegated to `Integrate\`BronsteinRational`;
      - logarithmic polynomial `P(x, Log[u])`: the recursive primitive-
        polynomial coefficient matching, with a limited-integration oracle that
        folds a would-be new logarithm back into the tower (e.g.
        `Integrate[Log[2 x + 3], x]`, `Integrate[Log[x]/x, x]`);
      - exponential (Laurent) polynomial `sum_i p_i(x) E^(i u)`, `u` polynomial
        in `x`, `i` positive or negative: the powers of `E^u` decouple and each
        `i != 0` term solves the Risch differential equation
        `q_i' + i u' q_i = p_i` by a polynomial ansatz
        (`Integrate[x E^x, x] = (x-1) E^x`, `Integrate[x E^(x^2), x] = E^(x^2)/2`,
        `Integrate[(E^x+E^(-x))/2, x] = Sinh[x]`);
      - Hermite reduction for a repeated pole of `theta = Log[u]` or
        `theta = E^u` (the latter when `D` is coprime to `theta`):
        `Q = H(theta)/Hden(theta) + sum_j c_j Log(g_j)` with
        `Hden = gcd(D, dD/dtheta)`, solved by `SolveAlways` over `theta` and `x`
        (`Integrate[1/(x (1+Log[x])^2), x] = -1/(1+Log[x])`,
        `Integrate[E^x/(1+E^x)^2, x] = -1/(1+E^x)`);
      - a coupled hyperexponential case (a unified ansatz
        `Q = sum_i w_i(x) E^(i u) + H(E^u)/Hden(E^u) + sum_j c_j Log(g_j)` solved
        by `SolveAlways` over `theta` and `x`) that closes mixed
        polynomial-plus-log exponentials such as
        `Integrate[1/(1 + E^x), x] = x - Log[1 + E^x]`, and — with the Hermite
        term `H/Hden` fused in (the `theta`-coprime denominator split into its
        repeated part `Hden = gcd(Dtil, dDtil/dtheta)` and squarefree radical) —
        the repeated / `theta = 0` exponential poles
        `Integrate[1/(1 + E^x)^2, x] = x + 1/(1 + E^x) - Log[1 + E^x]`,
        `Integrate[1/(E^x (1 + E^x)^2), x]`, `Integrate[1/(1 + E^x)^3, x]`;
      - a multi-kernel **sum-of-exponentials** case: an integrand that
        exponentializes to a sum `sum_k p_k(x) E^(W_k)` of NON-commensurate
        exponentials `E^(W_k)` (e.g. the `(1 ± I) x` pair from `E^x Sin[x]`)
        decouples — each term integrates by its own Risch DE
        `q_k' + W_k' q_k = p_k` — closing `Integrate[E^x Sin[x], x]`,
        `Integrate[x E^x Sin[x], x]`, `Integrate[E^(2x) Cos[3x], x]`, ... (via
        the complex exponentials the answer is left in an I-laden `Cosh`/`Sinh`
        form, a `Simplify` opportunity; the diff-back is exactly `0`);
      - nested **logarithmic** and **exponential tower** cases: a rational
        function of a chain of nested logarithms (`Log[x]`, `Log[Log[x]]`, ...) or
        nested exponentials (`E^x`, `E^(E^x)`, ...) is integrated over the tower
        derivation `D = d/dx + sum_i Dt_i d/dt_i` by one unified `SolveAlways`
        ansatz over all tower variables, closing
        `Integrate[1/(x Log[x] Log[Log[x]]), x] = Log[Log[Log[x]]]`,
        `Integrate[Log[Log[x]]/(x Log[x]), x] = Log[Log[x]]^2/2`,
        `Integrate[E^x E^(E^x), x] = E^(E^x)`, ...  The tower cases are
        bounded-ansatz searches and so are diff-back verified (a non-elementary
        integrand declines rather than returning a spurious form); a whole-tower
        rationality gate and a single-kernel nesting gate keep the other cases
        from ever certifying a wrong nested answer;
      - the **genuine one-extension-at-a-time recursion** (Bronstein ch. 5),
        which the flat tower ansatz above cannot express: it builds
        the ordered differential tower (structure-theorem triangularity check) and
        peels one kernel at a time, integrating the polynomial/Laurent part in the
        top kernel *coefficient by coefficient* — each coefficient integral is an
        integration in the LOWER field (the recursion), bottoming out in `C(x)`.
        This closes **mixed exp/log towers** and **rational lower-field
        coefficients** the single-kind flat cases decline, e.g.
        `Integrate[E^x/x + E^x Log[x], x] = E^x Log[x]`,
        `Integrate[Log[1+E^x] + x E^x/(1+E^x), x] = x Log[1+E^x]`,
        `Integrate[1/(x^2 Log[x]) - Log[Log[x]]/x^2, x] = Log[Log[x]]/x`; diff-back
        verified, so non-elementary mixed integrands (`E^x Log[x]`) decline.  A
        **proper rational part** at a logarithmic top level (a genuine `t_n`-pole)
        is integrated by tower **Hermite reduction + Rothstein–Trager**
        (`H(t)/Hden + sum_j c_j Log(g_j)`, `Hden = gcd(den, d den/dt_n)`, constant
        residues), closing `Integrate[1/(x Log[x] (1+Log[Log[x]])^2), x] =
        -1/(1+Log[Log[x]])` and `Integrate[1/(x(1+Log[x])) + E^x, x] = E^x +
        Log[1+Log[x]]` (a non-constant residue declines).  For an **exponential**
        top level the Laurent and log parts couple, so a proper `t_n`-pole is closed
        by a unified coupled-hyperexponential ansatz over the tower derivation, e.g.
        `Integrate[(2 Log[x]/x) E^(Log[x]^2)/(1+E^(Log[x]^2)), x] =
        Log[1+E^(Log[x]^2)]`.  The exponential Laurent step solves a **field Risch
        differential equation** `q_i' + i w_n' q_i = p_i` in the lower field for each
        power; a `q_i` that is rational there — with an arbitrary denominator, by the
        RDE denominator theorem `q_i = h/Denominator[p_i]` with `h` a bounded
        polynomial — is found, closing
        `Integrate[(2/x - 1/(x Log[x]^2)) E^(Log[x]^2), x] = E^(Log[x]^2)/Log[x]`
        (`q = 1/Log[x]`) and
        `Integrate[(2 Log[x]/(x(1+Log[x])) - 1/(x(1+Log[x])^2)) E^(Log[x]^2), x] =
        E^(Log[x]^2)/(1+Log[x])` (`q = 1/(1+Log[x])`, non-monomial denominator).
        A structural pre-pass re-splits an **evaluator-merged exponential monomial**
        `E^(a+b) -> E^a E^b` (the evaluator folds `E^x E^(E^x)` into `E^(x+E^x)`,
        whose exponent carries the foreign kernel `E^x` and breaks tower
        independence), restoring the independent basis `{E^x, E^(E^x)}` and closing
        `Integrate[E^x E^(E^x)/(1+E^(E^x)), x] = Log[1+E^(E^x)]` (the non-elementary
        `E^(E^x)/(1+E^(E^x))` still declines);
      - a trig/hyperbolic front-end (`TrigToExp` -> exponential machinery ->
        `ExpToTrig`) that closes `Sin`, `Cos`, `Sinh`, `Cosh`, `Sin[x]^2`,
        `Sin[x] Cos[x]`, `Tan`, `Tanh`, ...; through the complex substitution
        `Tan`/`Tanh` come out in a correct but I-laden form (e.g.
        `I x - Log[1 + E^(2 I x)] = -Log[Cos[x]]`) that no current simplifier
        reduces to real closed form (a `Simplify` improvement opportunity);
      - `K E^(a x^2 + b x + c)` (`a != 0`) → `Erf`/`Erfi`;
      - `(M E^(a x + b))/(c x + d)` → `ExpIntegralEi` — the `E^v` kernel is
        extracted directly, so a negative leading coefficient (`E^(-x)/x →
        ExpIntegralEi[-x]`) and a nonzero exponent constant close uniformly;
      - `c w^(p-1) w'/Log[w]` → `c LogIntegral[w^p]` — subsumes `K/Log[x] → K
        LogIntegral[x]` and adds a scaled/affine argument (`1/Log[2x] →
        LogIntegral[2x]/2`) and a monomial numerator (`x/Log[x] →
        LogIntegral[x^2]`);
      - `K Log[1 + p x]/x` → `PolyLog[2, -p x]` (dilogarithm);
      - fractional (Rothstein–Trager) log-part: a proper rational function of
        `theta` with squarefree denominator `prod g_i` gives `sum_i c_i Log(g_i)`,
        the constant residues `c_i` solved from `num = sum_i c_i D(g_i)(d/g_i)`
        via `SolveAlways` over `theta` and `x` (`Integrate[1/(x(1+Log[x])),x] =
        Log[1+Log[x]]`, `Integrate[E^x/(1+E^x),x] = Log[1+E^x]`); the residues are
        verified constant (an `x`-dependent `SolveAlways` pseudo-solution is
        rejected so `1/Log[c x+d]` is never mis-closed with a polynomial-coefficient
        logarithm);
      - **pure resultant Lazard–Rioboo–Trager** log-part (when the `SolveAlways`
        path above declines): an irreducible-over-Q denominator factor in `theta`
        with **algebraic (non-rational) residues** is closed by the exact resultant
        `Res_t(a - z D(d), d)` (residues) + Rioboo `LogToReal` (real `Log + ArcTan`
        form), delegated to the internal `Integrate`TranscendentalLogPart` (reuses
        the rational-LRT machinery with the monomial derivation), closing
        `Integrate[1/(x(Log[x]^2+1)),x] = ArcTan[Log[x]]`,
        `Integrate[E^x/(E^(2x)+1),x] = ArcTan[E^x]`, mixed
        `Integrate[(1+Log[x])/(x(Log[x]^2+1)),x] = ArcTan[Log[x]] + Log[1+Log[x]^2]/2`
        and higher-degree denominators; diff-back verified, so non-elementary
        siblings (`1/(Log[x]^2+1)`, `1/(x^2(Log[x]^2+1))`) decline; this same
        resultant LRT is **lifted into the tower recursion** at a logarithmic top,
        so a nested-log proper part with algebraic residues closes too —
        `Integrate[1/(x Log[x] (Log[Log[x]]^2+1)),x] = ArcTan[Log[Log[x]]]` — with
        the residues gated free of *every* lower-field variable (`{x, Log[x], …}`),
        and `1/(Log[x](Log[Log[x]]^2+1))` (x-dependent residues) declining;
    Rational-argument exponents such as `E^(1/x)` are handled (`Integrate[-E^(1/x)/
    x^2, x] = E^(1/x)`) via the `q = h/Denominator[p]` RDE ansatz.  **Multiplicatively
    commensurate merged kernels** are reduced in the tower builder — a collected
    exponential whose exponent is an integer multiple of a class primitive
    (`E^(2 E^x) = (E^(E^x))^2`) becomes a power of the primitive's tower variable
    instead of a spurious extra extension — so towers with `E^(k u)` close and the
    exp-top algebraic-residue LRT is unblocked
    (`Integrate[E^x E^(E^x)/(1+E^(2 E^x)), x] = ArcTan[E^(E^x)]`,
    `Integrate[E^x E^(2 E^x)/(1+E^(E^x)), x] = E^(E^x) - Log[1+E^(E^x)]`).  The class
    primitive is **synthesized** (`p = w_0/lcm(ratio denominators)`) rather than
    required to be a member, so *non-integer* commensurate exponents close too:
    `Integrate[1/(E^(x/2)+E^(x/3)), x]` (primitive `E^(x/6)`) and the nested
    `Integrate[D[E^(E^x/2)/(1+E^(2 E^x/3)), x], x]` (primitive `E^(E^x/6)`); only a
    class whose members are not all rational multiples of one another declines.
    The **RDE-solver degree bounds** are
    Bronstein's exact `RdeBoundDegree` (leading-degree balance,
    `deg_v(q) = deg_v(p) − deg_v(f)` where the exponential dominates), with **no
    arbitrary cap**, so an exponential-Laurent coefficient of any degree closes
    (`Integrate[(6 Log[x]^5 + 2 Log[x]^7)/x E^(Log[x]^2), x] = Log[x]^6 E^(Log[x]^2)`,
    and deg-20 …).  The **flat-tower and proper-part Hermite ansätze are likewise
    cap-free**: exact top-kernel log/exp Laurent bounds and derived inner-exp windows,
    so all top degrees close (`Integrate[Log[Log[x]]^5/(x Log[x]), x] =
    Log[Log[x]]^6/6`, `Integrate[E^x E^(6 E^x)/(1+E^(E^x)), x]`).  The
    **leading-coefficient cancellation / resonance** sub-case of `RdeBoundDegree`
    completes the SPDE degree machinery — the bound is widened monotonically to
    `max(naive, m_res)` at the Bronstein resonance integer `m_res` (detected live for the
    exponential top), correct by the same certification-and-diff-back gate.  Only
    algebraic extensions (`Sqrt`, `RootSum`) remain unimplemented, so integrands needing
    them return unevaluated.
  - `"CRCTable"` — `Integrate\`CRCTable[f, x]`.
  - `"Undefined"` — `Integrate\`Undefined[f, x]`.
  - `"Symmetry"` — origin-symmetry reduction for an interval `[-c, c]`
    (`Integrate\`Symmetry[f, {x, -c, c}]`): an odd integrand integrates to `0`,
    an even one to `2 Integrate[f, {x, 0, c}]`. The parity is proved by
    `Simplify`, and a value is claimed only when the half integral converges, so
    a divergent principal value is never reported as `0`. Under Automatic it runs
    after residue and before Newton-Leibniz.
  - `"Beta"` — Euler-Beta reduction on `[0,1]`
    (`Integrate\`Beta[f, {x, 0, 1}]`): `x^(k-1) (1-x)^(l-1) → Beta[k, l]`, with
    `Log[x]^i Log[1-x]^j` weights giving the mixed parameter derivative of
    `Beta`. Gated on `Re[k] > 0 && Re[l] > 0`.
  - `"TrigPower"` — `Sin[x]^m Cos[x]^n` over a canonical trig interval
    (`Integrate\`TrigPower[f, {x, 0, c}]`): over `[0, Pi/2]` it is
    `Beta[(m+1)/2, (n+1)/2]/2`; over `[0, Pi]`/`[0, 2Pi]` the standard parity
    multipliers apply (an odd power integrates to `0`).
  - `"NewtonLeibniz"` — the real-axis definite-integral mechanism (implicit for
    the `{x, a, b}` form); see **Definite integration** below.
  - `"LineIntegral"` — the complex contour mechanism (implicit for the
    `{x, z0, …, zn}` form); see **Complex line integration** below.
  - `"DiffUnderInt"` (alias `"DifferentiationUnderIntegral"`) — definite
    integration by differentiation under the integral sign (Feynman's trick);
    `Integrate\`DiffUnderInt[f, {x, a, b}]`. Tried last in the definite cascade
    (after residue and Newton-Leibniz). See **Differentiation under the integral
    sign** below.
  - `"SinPowerMonomial"` — `Sin[r x]^k / x^m` on `[0, Infinity)` (the ssp
    family); `Integrate\`SinPowerMonomial[f, {x, 0, Infinity}]`.
  - `"OscillatoryPower"` — Fresnel-type `Cos[b x^n]` / `Sin[b x^n]` on
    `[0, Infinity)`; `Integrate\`OscillatoryPower[f, {x, 0, Infinity}]`.
  - `"RationalLog"` — `R(x) Log[x]^n` on `[0, Infinity)` for a proper rational
    `R` with negative-real-axis poles; `Integrate\`RationalLog[f, {x, 0, Infinity}]`.
  - `"RamanujanMasterTheorem"` (alias `"Mellin"`) — half-line `∫₀^∞ x^{s-1} f(x) dx`
    by the Mellin-transform / Ramanujan Master Theorem method;
    `Integrate\`RamanujanMasterTheorem[f, {x, 0, Infinity}]`. Under Automatic it
    runs after Newton-Leibniz and before DiffUnderInt. See **Mellin / Ramanujan
    Master Theorem** below.
  The definite mechanisms name themselves only: the actual mechanism is
  chosen from the spec type, so on a definite integral any *other* method name
  is passed through to the inner indefinite integration that produces the
  antiderivative, and either definite-mechanism name on the indefinite
  `Integrate[f, x, …]` form is a no-op (stays unevaluated).
  Unknown method names emit `Integrate::method` and bubble back.
- Universal correctness predicate: `Cancel[Together[D[Integrate[f,x],x] - f]] === 0`.

#### Differentiation under the integral sign (`Integrate\`DiffUnderInt`)

For a parameter-dependent definite integral `I(p) = Integrate[f(x,p), {x,a,b}]`,
this method (Leibniz rule / "Feynman's trick") differentiates the integrand with
respect to a free parameter `p`, evaluates the resulting simpler definite
integral `J(p) = Integrate[D[f,p], {x,a,b}]`, integrates `J(p)` back over the
parameter, and fixes the constant of integration from an **exact** base value
`I(p0)` (a `p` where `f` vanishes identically, or reduces to a directly-
integrable form). Every case is the first-order ODE `I'(p) = J(p)`
(Boulnois 2023). Verification is symbolic and correct-by-construction
(`Simplify[D[I,p] - J] === 0` plus the exact base) — there is **no** numeric
crosscheck. Assumptions (`a > 0`, …) are honoured and used to clean the closed
forms.

Because the general integrator is slow/hangs on the parameter-dependent inner
integrals Feynman's trick produces, `DiffUnderInt` evaluates the standard
families itself with closed-form formulas: the **Laplace/Fourier half-line**
`∫₀^∞ xⁿ e^{-p x}{1,cos,sin} dx`, the **sinc/Frullani** `∫₀^∞ …/x dx`, the
**even-rational half-line** `∫₀^∞ P(x)/Q(x²) dx`, the **general (non-even)
rational half-line** `∫₀^∞ R(s) ds` (real `ArcTan`/`Log` boundary values — this
is what closes a *decaying* sinc such as `∫₀^∞ e^{-p x} Sin[q x]/x dx = ArcTan[q/p]`
whose Laplace image is non-even), and the **Gaussian moment** family
`∫₀^∞ xⁿ e^{-p x²}{1,cos} dx` in `Sqrt[Pi]`/`e^{-q²/4p}`. The Gaussian
parameter back-integration `∫ c e^{-k p²} dp` is supplied directly as an `Erf`
(the engine does not produce it). Forms outside these families (finite-period
trig, piecewise/`Min`-`Max` results, the Sin-Gaussian Dawson/Erfi moment) are
declined — the integral is returned unevaluated, fast, never a wrong value. The
two rational half-line families gate on the inner integrand being a **rational
function of `x`**: a differentiated exp-geometric/Mellin form (still carrying
`e^x` or `x^{s-1}`) is not, and feeding it to `Apart[·, x]` otherwise drives a
non-terminating rewrite — so such forms are declined up front and left for the
Ramanujan/Mellin method.

Worked examples that close:
`Integrate[(x^a-1)/Log[x], {x,0,1}]` → `Log[1+a]`;
`Integrate[Exp[-a x] Sin[b x]/x, {x,0,Infinity}, Assumptions->a>0]` → `ArcTan[b/a]`;
`Integrate[Sin[a x]^2/x^2, {x,0,Infinity}, Assumptions->a>0]` → `π a/2`;
`Integrate[Log[1+a^2 x^2]/(1+x^2), {x,0,Infinity}, Assumptions->a>0]` → `π Log[1+a]`;
`Integrate[Exp[-c x](1-Cos[a x])/x^2, {x,0,Infinity}, Assumptions->{a>0,c>0}]` → `a ArcTan[a/c] − (c/2) Log[1+a²/c²]`;
`Integrate[Exp[-x^2] Sin[a x]/x, {x,0,Infinity}]` → `(π/2) Erf[a/2]`.

#### Mellin / Ramanujan Master Theorem (`Integrate\`RamanujanMasterTheorem`)

The series/transform-based mechanism for half-line integrals
`∫₀^∞ x^{s-1} f(x) dx` of a *transcendental* `f` (the class residue and FTC do
not close). By Ramanujan's Master Theorem, if `f(x) = Σ (-1)^k φ(k) x^k / k!`
then the Mellin transform is `Γ(s) φ(-s)` on the fundamental strip. The
integrand is Expanded and, term by term, decomposed into `C · x^ρ · f(λx)`; the
kernel `f` is matched against a table of proven base Mellin transforms and the
power prefactor sets `s = ρ + 1`:

| kernel `f` | `∫₀^∞ x^{s-1} f dx` | strip |
|------------|--------------------|-------|
| `Exp[c x]`, `Re c<0` | `Γ(s) (-c)^{-s}` | `0<Re s` |
| `Exp[c x^2]`, `Re c<0` | `½ (-c)^{-s/2} Γ(s/2)` | `0<Re s` |
| `(p + q x^m)^{-a}`, `p,q>0` | `(1/m) p^{s/m-a} q^{-s/m} B(s/m, a-s/m)` | `0<Re s<m Re a` |
| `Cos[λ x]`, `λ>0` | `π / (2 Sin(πs/2) Γ(1-s)) λ^{-s}` | `0<Re s<1` |
| `Sin[λ x]`, `λ>0` | `π / (2 Cos(πs/2) Γ(1-s)) λ^{-s}` | `-1<Re s<1` |
| `Log[1 + λ x]`, `λ>0` | `π / (s Sin(π s)) λ^{-s}` | `-1<Re s<0` |
| `ArcTan[λ x]`, `λ>0` | `-π / (2 s Cos(πs/2)) λ^{-s}` | `-1<Re s<0` |
| `BesselJ[ν, λ x]`, `λ>0` | `2^{s-1} λ^{-s} Γ((ν+s)/2)/Γ((ν-s)/2+1)` | `-Re ν<Re s<3/2` |
| `pFq[{a}; {b}; -λ x]`, `λ>0` | `(∏Γ(b_j)/∏Γ(a_i)) Γ(s) (∏Γ(a_i-s)/∏Γ(b_j-s)) λ^{-s}` | `0<Re s<min Re a_i` |
| `PolyLog[ν, -λ x]`, `λ>0` | `π (-s)^{-ν} λ^{-s} / Sin(π s)` | `-1<Re s<0` |
| `1/(e^{c x}+γ)`, `c>0`, `-1≤γ≤1` | `c^{-s} Γ(s) (-1/γ) PolyLog(s, -γ)` | `0<Re s` (`1<Re s` if `γ=-1`) |

The last row is the **exponential-geometric** kernel of the statistical-mechanics
integrals: expanding `1/(e^{cx}+γ) = (-1/γ) Σ_{j≥1} (-γ)^j e^{-jcx}` and
integrating term by term lands on `PolyLog`. Its two headline specialisations are
`γ=-1` **Bose–Einstein** `∫₀^∞ x^{s-1}/(e^{cx}-1) = c^{-s} Γ(s) ζ(s)` (Planck /
Debye; the denominator zero at `x=0` tightens the strip to `Re s>1`) and `γ=+1`
**Fermi–Dirac** `∫₀^∞ x^{s-1}/(e^{cx}+1) = c^{-s} Γ(s) η(s)` (emitted as
`-Γ(s) PolyLog(s,-1)`, which stays finite at `s=1` where `(1-2^{1-s})ζ(s)` would
be `0·∞`). A **symbolic fugacity** is admitted too — the general Bose integral
`∫₀^∞ x^{s-1}/(z^{-1} e^x - 1) dx = Γ(s) PolyLog(s, z)` closes for a symbolic `z`
whenever the `Assumptions` confine `γ' = -z` to `(-1, 1]`. The built-in
assumption engine only discharges syntactic matches (it proves neither `1/z>0`
nor `-1≤-z≤1` from `0<z<1`), so the `-1<γ'≤1` gate is decided by a small **sound
interval-bound prover** over the parameter box read off the `Assumptions`:
interval arithmetic yields a guaranteed enclosure, so the gate never accepts an
inadmissible fugacity (an unbounded or out-of-range `z` simply declines).

Four operational layers extend the table:

- **Monomial substitution** `g(x^k)` (`k≠1`) via `y = x^k`:
  `∫ x^{s-1} g(x^k) = (1/k) ∫ y^{s/k-1} g(y)`, so `Sin[√x]`, `BesselJ[ν,2√x]`,
  `ArcTan[√x]`, `Cos[x²]` reduce to the linear table at `s/k`.
- **Hypergeometric reduction** (applied before Expand, so a cancellation kernel
  is never split): `Erf[u] → u·₁F₁`, `Γ[a]-Γ[a,x] → x^a/a·₁F₁` (lower incomplete
  gamma), and the product `BesselJ[ν,·]² → ₁F₂` (a Mellin convolution closed via
  the `J²` identity rather than a Barnes integral).
- **Parametric differentiation** for `Log[1+λx]^n (1+λx)^{-w₀}`:
  `M = (-1)^n ∂ⁿ_w[λ^{-s} B(s, w-s)]|_{w=w₀}`, strip `-n<Re s<w₀`.
- **`Log[x]^k` weight** (a bare `Log[x]`, distinct from the `Log[1+λx]` kernel):
  since `∂_s x^{s-1} = x^{s-1} Log x`, a `Log[x]^k` factor is the `k`-th
  `s`-derivative of the base transform `M_R(s)`. The open strip carries unchanged
  (`Log^k` is dominated by `x^{±ε}`), so e.g. `∫₀^∞ Log[x]/(1+x²) dx = 0` and
  `∫₀^∞ x Log[x] e^{-x} dx = Γ'(2) = 1-γ`.
- The `pFq` transform is the master kernel — `1F1`, `2F1`, `3F2`, … close
  uniformly (`Hypergeometric1F1`/`2F1` are stored as `HypergeometricPFQ`).

A **Frullani pre-pass** (run on the whole integrand before Expand, since each
half is individually divergent) recognises `(f(a x)-f(b x))/x` and returns
`(f(0⁺)-f(∞)) Log(b/a)` (`a,b>0`): the scale ratio is read structurally and the
pairing verified by the exact identity `(t₁ /. x→ρx)+t₂ = 0`; the boundary values
are the finite limits of `f`. So `∫₀^∞ (e^{-2x}-e^{-5x})/x dx = Log(5/2)` and
`∫₀^∞ (ArcTan(5x)-ArcTan(2x))/x dx = (π/2) Log(5/2)`.

Each application is **gated on its convergence strip** — checked by `Simplify`
(numerically for a numeric `s`, or against the supplied `Assumptions` for a
symbolic `s`), so every result is correct by construction. **When the
assumptions do not prove the strip, the value is returned as a
`ConditionalExpression[value, strip]`** (matching Wolfram) — it collapses to the
bare value once the strip is proved and to `Undefined` if it is refuted. A
provably-violated strip declines. Verification is symbolic; there is **no**
numeric crosscheck (the trig/PolyLog transforms use reflection-formula forms
regular at `s=0`, so e.g. `∫₀^∞ Sin[x]/x dx = π/2` falls out with no limit). A
sum is integrated term by term (each term must converge on its own). Out of
scope — products of three or more transcendental kernels, finite intervals, and
two-sided reductions — return unevaluated, never a wrong value.

Worked examples that close:
`Integrate[Exp[-x^2], {x,0,Infinity}]` → `√π/2`;
`Integrate[x^(s-1) Exp[-x], {x,0,Infinity}]` → `ConditionalExpression[Γ[s], s>0]`;
`Integrate[x^(s-1) BesselJ[ν,2√x]/x^(ν/2), {x,0,Infinity}]` → `Γ[s]/Γ[1+ν-s]` (Ramanujan's canonical example);
`Integrate[x^(s-1) (Γ[a]-Γ[a,x])/x^a, {x,0,Infinity}]` → `Γ[s]/(a-s)`;
`Integrate[x^(s-1) Hypergeometric2F1[a,b,c,-x], {x,0,Infinity}]` → `Γ[c]Γ[s]Γ[a-s]Γ[b-s]/(Γ[a]Γ[b]Γ[c-s])`;
`Integrate[Sin[x]/x, {x,0,Infinity}]` → `π/2`;
`Integrate[BesselJ[0, x], {x,0,Infinity}]` → `1`;
`Integrate[x^3/(Exp[x]-1), {x,0,Infinity}]` → `π⁴/15` (Debye);
`Integrate[x^3/(Exp[x]+1), {x,0,Infinity}]` → `7π⁴/120` (Fermi–Dirac);
`Integrate[x^(s-1)/(Exp[x]-1), {x,0,Infinity}, Assumptions→s>1]` → `Γ[s] ζ[s]`;
`Integrate[(Exp[-2x]-Exp[-5x])/x, {x,0,Infinity}]` → `Log[5/2]` (Frullani);
`Integrate[Log[x]/(1+x^2), {x,0,Infinity}]` → `0`.

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

### Definite integration (Newton-Leibniz)

`Integrate[f, {x, xmin, xmax}]` gives the definite integral by the
**fundamental theorem of calculus** (`src/calculus/integrate_newton_leibniz.c`,
exposed explicitly as `Integrate\`NewtonLeibniz[f, {x, xmin, xmax}]` and as
`Method -> "NewtonLeibniz"`).  The mechanism:

1. Antidifferentiate `f` with the ordinary indefinite cascade to get `F`.  If
   no closed antiderivative exists, the definite integral is left unevaluated
   — never assigned a wrong value.
2. Locate the real poles of the **integrand** `f` strictly inside
   `(xmin, xmax)` via `Integrate\`SingularPoints` (roots of
   `Denominator[Together[f]]`).  A continuous `f` has an antiderivative that is
   pole-free where `f` is, so only `f`'s own poles bound the segments.
3. Split `[xmin, xmax]` at those poles and form the telescoping sum
   `Σ (F(p_{i+1}⁻) − F(p_i⁺))`, evaluating each boundary through the `Limit`
   engine: a plain limit at infinite endpoints, a one-sided limit
   (`Direction -> "FromBelow"/"FromAbove"`) at interior poles, and direct
   substitution at ordinary finite endpoints.  Improper integrals thereby
   acquire their correct finite value; a genuinely divergent integral emits
   `Integrate::idiv` and is left unevaluated (matching Mathematica).

`Integrate[f, {x, a, b}, {y, c, d}, …]` gives the iterated multiple integral,
reduced innermost-first (the last spec is the inner integral), so an inner
bound may depend on an outer variable.

`Integrate\`SingularPoints[expr, {x, a, b}]` returns the sorted list of the
interior real poles used for the split — exposed for inspection and reuse.

**Cauchy principal value.** `Integrate[f, {x, a, b}, PrincipalValue -> True]`
computes the Cauchy principal value across an interior pole. It is defined only
for **odd-order** poles (the integrand changes sign across the pole, so the
symmetric one-sided divergences cancel); the value is then `Re[F(b) − F(a)]`
with the real-branch antiderivative — the branch-cut crossing at each pole
contributes a purely imaginary `I Pi` (residue) that `Re` removes. An
**even-order** interior pole has no principal value and emits `Integrate::idiv`.
With no interior pole the option is a no-op. Examples: `∫₋₁¹ dx/x = 0`,
`∫₀³ dx/(x−2) = −Log 2`, `∫₀² x/(x²−1) dx = ½ Log 3`.

**Continuous branch-form antiderivatives.** Many continuous periodic integrands
(e.g. `1/(2 + Cos[x])`) antidifferentiate to a Weierstrass branch form that is
*already continuous* — an `ArcTan[Tan[x/2] …]` term whose jump is exactly
cancelled by an accompanying `Floor` step — so `F(b) − F(a)` by direct
substitution is the correct value with no interior split.  To stay correct when
an antiderivative could instead carry a genuine *uncorrected* branch jump, a
result built from a step head (`Floor`, `Sign`, …) or an inverse-trig node over
a pole-bearing trig head is accepted only when a numeric `NIntegrate`
cross-check agrees; otherwise the integral is left unevaluated rather than
risking a wrong value.

**Examples**:
```mathematica
In[1]:= Integrate[x^2, {x, 0, 1}]
Out[1]= 1/3

In[2]:= Integrate[1/(1 + x^2), {x, 0, Infinity}]
Out[2]= 1/2 Pi

In[3]:= Integrate[1/Sqrt[x], {x, 0, 1}]            (* improper, convergent *)
Out[3]= 2

In[4]:= Integrate[1/(2 + Cos[x]), {x, 0, 2 Pi}]    (* continuous branch form *)
Out[4]= (2 Pi)/Sqrt[3]

In[5]:= Integrate[1/x, {x, -1, 1}]                 (* divergent *)
Integrate::idiv: Integral of 1/x does not converge on {-1, 1}.
Out[5]= Integrate[1/x, {x, -1, 1}]

In[6]:= Integrate[x y, {x, 0, 1}, {y, 0, 1}]       (* iterated *)
Out[6]= 1/4

In[7]:= Integrate`SingularPoints[1/((x - 1)(x - 2)), {x, 0, 3}]
Out[7]= {1, 2}
```

### Complex line / contour integration

When a `{x, a, b}` spec has a **non-real endpoint**, or when a spec lists more
than two points `{x, z0, z1, …, zn}`, `Integrate` evaluates the **contour
integral** of `f` along the straight segments `z0 → z1 → … → zn` in the complex
plane (`src/calculus/integrate_line.c`, exposed explicitly as
`Integrate\`LineIntegral[f, {x, z0, …, zn}]`).  Real-endpoint two-point specs
still go through the real-axis Newton-Leibniz path above.

Each segment `a → b` is parametrised by a **real** parameter,
`γ(t) = a + t (b − a)`, `t ∈ [0, 1]`, which reduces the complex problem to the
real machinery already in place:

1. Antidifferentiate `f` in `x` to get `F` (bail if unknown).
2. **On-path singularities** become real roots `t* ∈ (0, 1)` of
   `Denominator[Together[f(γ(t))]]` — a singularity strictly on the contour
   makes the integral divergent (`Integrate::idiv`, left unevaluated).
   `Integrate\`PathSingularPoints[f, {x, z0, …, zn}]` returns those points.
3. The segment value is the continuous change of `F` along the segment, with
   endpoint values taken as **real one-sided limits in `t`** when substitution
   is singular (so a complex-ray approach is a real one-sided limit the `Limit`
   engine can take).  For rational integrands whose antiderivative is a sum of
   logarithms / inverse-tangents of **affine** arguments, the branch-correct
   value is recovered by combining each into a single principal `Log` of a ratio
   `Log[(u(b))/(u(a))]` — exact because a straight segment subtends an angle
   `< π` at any point off it, so closed-contour residues come out exactly.
4. Every segment value is numerically cross-checked against a complex quadrature
   of `f(γ(t)) γ'(t)`; an uncorrectable branch crossing leaves the integral
   unevaluated rather than returning a wrong branch.

**Examples**:
```mathematica
In[1]:= Integrate[1/x, {x, 1 - I, 2 + 3 I}]
Out[1]= Log[2 + 3 I] - Log[1 - I]

In[2]:= Integrate[z^2, {z, 0, 1 + I}]
Out[2]= -2/3 + 2/3 I

In[3]:= Integrate[1/Sqrt[z], {z, 0, 1 + I}]        (* branch-point endpoint *)
Out[3]= 2 Sqrt[1 + I]

In[4]:= Integrate`PathSingularPoints[1/z, {z, -1 - I, 1 + I}]
Out[4]= {0}

In[5]:= Integrate[1/z, {z, -1 - I, 1 + I}]         (* pole on the path *)
Integrate::idiv: Integral of 1/z does not converge on the contour {-1 - I, 1 + I}.
Out[5]= Integrate[1/z, {z, -1 - I, 1 + I}]

In[6]:= Chop[N[Integrate[1/z, {z, 1, I, -1, -I, 1}]]]   (* CCW loop about 0 *)
Out[6]= 0.0 + 6.28319 I                                 (* = 2 Pi I *)
```

### Contour / residue-theorem definite integration

For the classical families of **improper** and **periodic** real integrals that
complex analysis dispatches by summing residues over the poles enclosed by a
standard contour, `Integrate[f, {x, a, b}]` runs a residue-theorem method
**before** Newton-Leibniz (also reachable as `Method -> "Residue"` or
`Integrate`ContourResidue[f, {x, a, b}]`).  Each answer is
**correct-by-construction**: once a family's structural gates hold, the residue
theorem gives the exact value, so there is *no* numeric quadrature crosscheck.
The only post-hoc gate is a self-consistency check that the residue sum closed to
a scalar (no surviving `x`/`Root`) and — for the real-valued families — that its
imaginary part vanishes (a residual `Im` would betray a mis-classified pole); a
failure returns unevaluated and the Newton-Leibniz path takes over.  Four
recognizers:

- **Rational on `(-∞, ∞)`** — `f = P/Q` with `deg Q ≥ deg P + 2` and **no real
  pole**: value `= 2 π i · Σ Res` over the poles in the upper half-plane.
- **Fourier / Jordan on `(-∞, ∞)`** — `f = R(x) · K` with
  `K ∈ {Cos[a x], Sin[a x], Exp[I a x]}`, `R` rational with `deg`-drop `≥ 1`,
  `a` a nonzero real: with `J = 2 π i · Σ_UHP Res[R Exp[I a x]]` (for `a > 0`),
  `∫ R Cos = Re[J]` and `∫ R Sin = Im[J]`. The bare complex-exponential kernel is
  recognised in both spellings — `Exp[I a x]` and the evaluator-normalised
  `Power[E, I a x]` — and returns `J` directly (`∫ R Exp[I a x] = J`); for
  negative `a` the lower half-plane is closed. Real-exponent `Exp[a x]` is *not*
  a Fourier kernel (it routes to the rectangular family). `Cos`/`Sin` are closed
  via `Re`/`Im` of the single decaying exponential, not split into two.
- **Rational-in-`{Sin, Cos}` over a full period** `(0, 2π)` or `(-π, π)` — via
  `z = Exp[I x]` on the unit circle: value `= 2 π i · Σ Res` over the poles
  inside the unit disk.
- **Removable axis singularity (Fourier)** — a **simple** real-axis pole of `R`
  at which the kernel vanishes (so `f = R·K` is analytic there, e.g. `Sin[x]/x`
  at `0`) contributes a **half residue** `π i · Res`, the indented-contour value,
  giving `∫ Sin[x]/x = π`.  A *genuine* axis pole (kernel nonzero there, e.g.
  `Cos[x]/x`) makes the ordinary integral diverge and returns unevaluated —
  plain `Integrate` does not compute a principal value.

A one-line **half-line** add-on covers even integrands:
`∫₀^∞ f = ½ ∫₋∞^∞ f`.  Three further recognizers handle branch-cut and
symbolic-exponent contours:

- **Keyhole / Mellin on `(0, ∞)`** — a branch power times a rational function,
  `f = x^p R(x)` with `p` non-integer (so `s = p + 1 ∉ ℤ`): value
  `= −π · Σ_k Res[ z^{s-1} R(z), z_k ] · e^{−iπs} / sin(π s)` over the poles
  `z_k` of `R`, with `z^{s-1}` on the branch `arg ∈ (0, 2π)`.  The residue is of
  the **full** integrand `z^{s-1} R(z)` — for a simple pole this equals
  `z_k^{s-1} Res(R, z_k)`, but for an order-≥2 pole it does not (a pure double
  pole has `Res(R)=0`), so poles of any order are handled by shifting `w=z−z_k`
  and expanding the analytic factor `(1+w/z_k)^{s-1}`.  The `e^{−iπs}/sin(π s)`
  prefactor is the exact reduction of the keyhole jump `1/(1 − e^{2πi s})`,
  landing a numeric `s` on an algebraic multiple of `π` (e.g.
  `∫₀^∞ x^{1/3}/(x²+1) = π/√3`, `∫₀^∞ √x/(1+x)² = π/2`).  Requires
  `0 < Re(s) < deg Q − deg P`.
- **Sector on `(0, ∞)`** — `f = x^m/(c + x^n)` with the exponent `n` possibly a
  **symbolic parameter**: the wedge of angle `2π/n` gives
  `(π/n) c^{s/n − 1} csc(π s/n)`, `s = m + 1`.  This is the one family admitting a
  symbolic `n` (the keyhole cannot enumerate `n` poles), powering
  `Integrate[1/(1 + x^n), {x, 0, ∞}, Assumptions -> n > 1] = (π/n) csc(π/n)`.
- **Rectangular / quasi-periodic on `(-∞, ∞)`** — `f = Exp[c x] R(Exp[x])`
  (period `2πi`): reduced to the keyhole core by `w = Exp[x]`
  (`∫_{-∞}^∞ f dx = ∫₀^∞ f(Log w)/w dw = ∫₀^∞ w^{c-1} R(w) dw`), so e.g.
  `Integrate[Exp[a x]/(Exp[x]+1), {x, -∞, ∞}, Assumptions -> 0 < a < 1] = π csc(π a)`.

**Assumptions and symbolic parameters.**  An `Integrate[f, {x, a, b},
Assumptions -> …]` option lets the residue families evaluate integrals whose
parameters are symbolic (`a > 0`, `0 < a < 1`, `n > 1`, …).  The recognizers
classify a parameter-dependent pole or kernel frequency by reading its sign at a
single generic point of the region the assumptions pin (a sign-consistent
instantiation), while the residue arithmetic stays fully symbolic — so the
closed form is still **correct by construction, with no numeric crosscheck**.
Convergence/applicability gates (`n > m + 1`, `0 < s < deg`-drop, `c > 0`) are
verified against the assumption-**guaranteed** interval bounds, not the sample
point, so an under-constrained problem (e.g. only `n > 0` for the sector family)
is refused rather than guessed; a parameter the assumptions leave two-sided
unbounded is likewise refused.  Radical pole locations are `PowerExpand`-cleaned
under all-positive parameters (`Sqrt[-4 a²] → 2 I a`) so a rational answer closes
to `π/a` rather than a `Sqrt[-4 a²]` surface, and real-parameter conjugation
uses `I → −I` (the symbolic `Conjugate` would not reduce).

```
In[0a]:= Integrate[Cos[k x]/(x^2 + a^2), {x, -Infinity, Infinity},
           Assumptions -> {a > 0, k > 0}]
Out[0a]= (Pi E^(-a k))/a

In[0b]:= Integrate[x^(1/3)/(x^2 + 1), {x, 0, Infinity}]     (* keyhole/Mellin *)
Out[0b]= Pi/Sqrt[3]

In[0c]:= Integrate[1/(1 + x^n), {x, 0, Infinity}, Assumptions -> n > 1]  (* sector *)
Out[0c]= (Pi Csc[Pi/n])/n

In[0d]:= Integrate[Exp[a x]/(Exp[x] + 1), {x, -Infinity, Infinity},
           Assumptions -> 0 < a < 1]                         (* rectangular *)
Out[0d]= Pi Csc[Pi a]

In[1]:= Integrate[1/(1 + x^4), {x, -Infinity, Infinity}]
Out[1]= Pi/Sqrt[2]

In[2]:= Integrate[1/(1 + x^2)^2, {x, -Infinity, Infinity}]   (* order-2 pole *)
Out[2]= 1/2 Pi

In[3]:= Integrate[Cos[x]/(1 + x^2), {x, -Infinity, Infinity}]
Out[3]= Pi/E

In[4]:= Integrate[1/(2 + Cos[x]), {x, 0, 2 Pi}]
Out[4]= (2 Pi)/Sqrt[3]

In[5]:= Integrate[Sin[x]/x, {x, -Infinity, Infinity}]        (* principal value *)
Out[5]= Pi

In[6]:= Integrate[1/(1 + x^4), {x, 0, Infinity}]             (* even half-line *)
Out[6]= (1/2 Pi)/Sqrt[2]
```

The residue sums are closed with `RootReduce` (algebraic families) or the
`Conjugate` identities for `Re[J]`/`Im[J]` (the Fourier family); a value that
still contains an unreduced `Root`, or a surviving imaginary part, is treated as
a mis-fire and returned unevaluated.  Negative controls such as
`Integrate[1/(1 + x^3), {x, -Infinity, Infinity}]` and
`Integrate[1/(x^2 - 1), {x, -Infinity, Infinity}, Method -> "Residue"]`
(genuine real-axis poles),
`Integrate[Cos[x]/x, {x, -Infinity, Infinity}, Method -> "Residue"]`
(genuine axis pole, kernel nonzero),
`Integrate[1/Sqrt[1 + x^4], {x, -Infinity, Infinity}]` (branch point, not
rational), and `Integrate[1/(2 + Cos[x]), {x, 0, Pi}]` (not a full period) all
stay unevaluated.  The keyhole/Mellin, sector and rectangular families
described above extend the reach to branch-cut and symbolic-exponent contours;
a log-keyhole (`∫₀^∞ Log[x] R(x)`) with symbolic on-circle poles remains out of
scope (it needs assumption-aware `Arg`/`Log` branch reasoning that Mathilda does
not yet have).

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

The table covers the inverse-trig and inverse-hyperbolic families
(Formulas 427–464 and their `427h`–`464h` analogs) among many others.
Rules whose denominator carries a *squared* coefficient
(`Sqrt[1 ∓ a^2 x^2]`, `(1 ± a^2 x^2)`) bind that coefficient linearly as
`a_` via `c_ + a_. x_^2` and recover the linear coefficient as `Sqrt[±a]`
on the RHS (with a `Condition` linking it to the numerator coefficient) —
the matcher does not invert a square, so `a_^2` in a pattern will not
bind. Some more elaborate multi-argument `/;`-guarded rules still do not
fire; this is a separate issue tracked under the matcher work.

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

## Integrate`TranscendentalLogPart

The transcendental (recursive-Risch) Lazard-Rioboo-Trager log part —
the monomial-extension generalisation of `Integrate`IntRationalLogPart`,
consumed by `Integrate`RischTranscendental` for frac-case integrands with
algebraic residues.

- `Integrate`TranscendentalLogPart[a, d, tau, z, Dd, g]` integrates the
  proper, squarefree-denominator fraction `a(tau)/d(tau)` of a
  transcendental monomial `tau` (= `Log[u]` or `E^u`), where `Dd` is the
  monomial derivation `D(d)`, `z` is the residue variable, and `g` is the
  residue-constant gate.  The residues are the roots of
  `Res_tau(a - z Dd, d)`; each squarefree factor with its log argument is
  converted to real `Log + ArcTan` form by Rioboo's `LogToReal`.  The gate
  `g` is either a **symbol** (single-kernel extension over `Q(x)`, residues
  required free of `x`) or a **List of symbols** (a tower proper part,
  residues required free of `x` and every lower tower variable
  `t_0, …, t_{L-1}` — i.e. constants of the whole tower derivation).
  Returns the antiderivative as a function of `tau` (the caller substitutes
  `tau -> kernel`), or unevaluated when the integral is not elementary in
  this form (a residue depending on a gate variable) / a factor exceeds
  `LogToReal`'s bounded scope.

```mathematica
(* 1/(x (Log[x]^2+1)): tau = Log[x], D(d) = (1+t)^2 *)
In[1]:= Integrate`TranscendentalLogPart[1, x + x t^2, t, z, 1 + 2 t + t^2, x]
Out[1]= ArcTan[t]
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
