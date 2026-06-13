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
- `Dt[f, x]`, `Dt[f, {x, n}]`, `Dt[f, x, y, ...]` -- delegate to
  `D[...]` for partial-derivative behaviour.

**Features**:
- `Protected`, `ReadProtected`.
- Shares the elementary-function derivative table with `D`; the
  only dispatch difference is the base-case handling of symbols.

**Examples**:
```mathematica
In[1]:= Dt[y^2 + Sin[x]]
Out[1]= 2 y Dt[y] + Cos[x] Dt[x]

In[2]:= Dt[Pi + 3 + x y]
Out[2]= x Dt[y] + y Dt[x]

In[3]:= Dt[y^2, x]
Out[3]= 0

In[4]:= Dt[x^2, {x, 2}]
Out[4]= 2
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
- Nine-stage dispatch cascade (`DerivativeDivides`, `LinearRadicals`,
  `QuadraticRadicals` and `LinearRatioRadicals` added 2026-06-06; `Weierstrass`
  added 2026-06-09):
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
  6. `Integrate\`Weierstrass[f, x]` — rational functions of the trig kernels
     `Sin/Cos/Tan/Cot/Sec/Csc[x]` (or hyperbolic `Sinh/Cosh/.../Csch[x]`) with a
     kernel in a denominator; continuous `Tan[x/2]` / `Tanh[x/2]` substitution
     (Jeffrey & Rich 1994).  Runs ahead of `DerivativeDivides`: it is
     domain-specific, deterministic, correct by construction, and yields a real,
     continuous antiderivative rather than a complex-logarithm form.
  7. `Integrate\`DerivativeDivides[f, x]` — substitution `u(x)`; in the
     cascade the quiet, branch-correct **direct quotient** strategy only.
  8. `Integrate\`RischNorman[f, x]` — Bronstein pmint, all integrands.
  9. `Integrate\`CRCTable[f, x]` — CRC integral table lookup (lazy-loaded
     from `src/internal/CRCMathTablesIntegrals.m` on first call).
  If every stage gives up the call bubbles back unevaluated.
- `Method -> "<name>"` option (3rd argument) bypasses the cascade and
  dispatches strictly to a single subroutine, with no fallback:
  - `"Automatic"` — default cascade above.
  - `"BronsteinRational"` — `Integrate\`BronsteinRational[f, x]`.
  - `"DerivativeDivides"` — `Integrate\`DerivativeDivides[f, x]` (direct **and**
    the more thorough Eliminate/Solve branch search).
  - `"LinearRadicals"` — `Integrate\`LinearRadicals[f, x]`.
  - `"QuadraticRadicals"` — `Integrate\`QuadraticRadicals[f, x]`.
  - `"LinearRatioRadicals"` — `Integrate\`LinearRatioRadicals[f, x]`.
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
`Method -> "Polynomial" | "Geometric" | "Gosper"` forces a single
algorithm (strict, no fallback), and now also takes effect on finite
unit-step integer ranges.

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
(the base of the exponential factors, combined across several `r^i`
factors).  The antidifference has the form `q(i) r^i`; `q` solves
`r q(i+1) - q(i) = p(i)` by undetermined coefficients.

```mathematica
In[1]:= Sum[a^i, i]
Out[1]= a^i/(-1 + a)
In[2]:= Sum[q1^i q2^i, i]
Out[2]= (q1 q2)^i/(-1 + q1 q2)
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
whose summand is a hypergeometric term — `t(k+1)/t(k)` rational in `k`. The
reindexed term ratio is factored into monic linear factors over the rationals;
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
```

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


## NResidue

Numerical residue by contour integration.  `NResidue[expr, {z, z0}]`
estimates the residue of `expr` at `z = z0` -- the coefficient of
`(z - z0)^-1` in the Laurent expansion -- by integrating around a small
circle in the complex plane.  Implemented natively in C in
`src/numerical_calculus/`: the reusable periodic-trapezoidal contour core
lives in `quadrature.{c,h}`, the builtin in `nresidue.{c,h}`.  Attribute:
`Protected`.

Unlike the symbolic `Residue` (which needs a power series at `z0`),
`NResidue` works for **essential singularities** such as `Exp[1/x]` and
`Sin[1/x]`.  It cannot tell a tiny spurious residual from a true zero --
`Chop` the result when needed -- and returns an incorrect value if the
contour encloses another singularity or crosses a branch cut.

### Method

The residue equals the Cauchy integral over the circle of radius `r`:

```
Res(f, z0) = (1/2 pi i) oint f dz = (r/N) sum_{k=0}^{N-1} f(z0 + r e^{i th_k}) e^{i th_k},  th_k = 2 pi k / N.
```

The integrand is 2*pi*-periodic and analytic in *theta*, so the periodic
trapezoidal rule converges geometrically (Trefethen & Weideman, *SIAM
Review* 2014).  The engine doubles `N` (reusing samples) until
`|S_{2N} - S_N|` meets the precision goal, applies Aitken/Shanks
extrapolation to the doubling sequence, and runs at machine precision or,
when `WorkingPrecision` requests it, in MPFR complex arithmetic.

### Forms

- `NResidue[expr, {z, z0}]` -- residue of `expr` near `z = z0`.
- `NResidue[{e1, e2, ...}, {z, z0}]` -- threads element-wise over the
  first argument (manual threading; `NResidue` is deliberately **not**
  `Listable`, so the `{z, z0}` spec is never split).

### Options

| Option | Default | Meaning |
|--------|---------|---------|
| `Radius`           | `1/100`         | Contour radius, or `Automatic` for an adaptive (Fornberg/Bornemann-style) search that favours fast-converging radii. |
| `WorkingPrecision` | `MachinePrecision` | Machine doubles, or MPFR at the requested decimal precision. |
| `PrecisionGoal`    | `WorkingPrecision - 2` | Target accuracy (decimal digits). |
| `MaxRecursion`     | `10`            | Maximum number of `N`-doublings. |
| `Method`           | `Trapezoidal`   | Only the trapezoidal rule is implemented. |

### Beyond Mathematica's NResidue

- `Radius -> Automatic` removes the manual radius tuning the reference
  implementation requires (its own `Exp[1/x]` example fails until you
  guess `Radius -> 1`).
- A relative-jump / decay-rate diagnostic flags a branch cut that
  **crosses** the contour (`NResidue::bcut`) instead of silently returning
  garbage.  A cut lying entirely *inside* the disk (e.g.
  `Sqrt[x-1] Sqrt[x+1]` on `Radius -> 2`, where the integrand is analytic
  on the circle) is undetectable on the contour and returns the same value
  as Mathematica.
- Aitken/Shanks extrapolation plus a reported error estimate.

### Diagnostics (stderr)

| Tag                | Triggered when |
|--------------------|----------------|
| `NResidue::ivar`   | Second argument is not a `{z, z0}` list, or the variable is not a symbol. |
| `NResidue::nnum`   | `z0` is not numeric, or `expr` did not evaluate to a number on the contour. |
| `NResidue::ncvi`   | Did not converge to the precision goal within `MaxRecursion` (the best estimate is still returned). |
| `NResidue::bcut`   | The integrand appears non-analytic on the contour (a branch-cut crossing); the result is unreliable. |
| `NResidue::badopt` / `NResidue::badmeth` | Invalid option value / unsupported `Method`. |

### Examples

```mathematica
In[1]:= NResidue[1/x, {x, 0}]
Out[1]= 1.

In[2]:= NResidue[Sin[1/(10 x)], {x, 0}] // Chop
Out[2]= 0.1

In[3]:= NResidue[1/(1.7 - 2.7 z + z^2), {z, 1.}] // Chop
Out[3]= -1.42857

In[4]:= NResidue[Exp[1/x], {x, 0}, Radius -> 1] // Chop
Out[4]= 1.

In[5]:= NResidue[{Exp[1/x], Sin[1/x], Cos[1/x]}, {x, 0}, Radius -> 1] // Chop
Out[5]= {1., 1., 0}

In[6]:= NResidue[1/x + 1/(x + 0.005), {x, 0}, Radius -> 0.001] // Chop
Out[6]= 1.

In[7]:= NResidue[Exp[1/x], {x, 0}, Radius -> Automatic] // Chop
Out[7]= 1.

In[8]:= 10! NResidue[Zeta[x]/x^11, {x, 0}, Radius -> 1/2, WorkingPrecision -> 30]
Out[8]= -3.6287999994567658842202915*10^6   (= Derivative[10][Zeta][0])
```


## ND

Numerical derivative.  `ND[expr, x, x0]` gives a numerical approximation to
the derivative of `expr` with respect to `x` at `x = x0`; `ND[expr, {x, n},
x0]` gives the `n`-th derivative.  Implemented natively in C in
`src/numerical_calculus/nderiv.{c,h}`.  Attribute: `Protected`.  Like
`NResidue`, `ND` cannot tell a tiny spurious residual from a true zero --
`Chop` when needed.

### Methods

**`Method -> EulerSum`** (default).  Richardson (Romberg/Neville)
extrapolation of the `n`-th **forward** finite difference taken along the
complex direction `Scale = s`:

```
D(h) = (1/(s h)^n) sum_{k=0}^n (-1)^{n-k} C(n,k) f(x0 + k s h),   h_i = 2^-i,
T(i,0) = D(h_i),   T(i,j) = T(i,j-1) + (T(i,j-1) - T(i-1,j-1)) / (2^j - 1),
result = T(Terms-1, Terms-1).
```

The forward (one-sided) stencil along `s` is what produces directional and
one-sided derivatives -- the left/right derivatives of `Abs`, and complex
directions such as `Scale -> 1 + I`.  Because a forward difference has an
error expansion in *all* powers of `h`, the Richardson denominator is
`2^j - 1` (not the `4^j - 1` of a central stencil).  EulerSum samples only
along `s`, so it works for **non-analytic** `expr` (e.g. `Re[Cos[I y]]`).  It
requires an integer order `n >= 0` and, for high-order derivatives, fights
subtractive cancellation via higher `WorkingPrecision` and more `Terms`.

**`Method -> NIntegrate`**.  Cauchy's integral formula, evaluated by reusing
`NResidue`:

```
f^(n)(x0) = n! Res_{z=x0} f(z)/(z - x0)^(n+1)
          = Gamma(n+1) * NResidue[expr/(x-x0)^(n+1), {x, x0}, Radius -> Scale].
```

`Scale` is the contour radius (default `1`).  `Gamma(n+1)` in place of `n!`
lets the order be **fractional or complex** (e.g. `ND[x, {x, -1/2}, 1]` =
`4/(3 Sqrt[Pi])`).  This method requires `expr` to be **analytic** near `x0`
and silently returns the derivative of the analytic continuation otherwise
(so `Re[Cos[I y]]` gives a wrong answer -- use EulerSum there).

### Forms and options

- `ND[expr, x, x0]` -- first derivative at `x0` (equivalent to `{x, 1}`).
- `ND[expr, {x, n}, x0]` -- `n`-th derivative.
- `ND[{e1, e2, ...}, x, x0]` -- threads element-wise over the first argument
  (manual threading; `ND` is deliberately **not** `Listable`, which would
  split the `{x, n}` spec).

| Option | Default | Meaning |
|--------|---------|---------|
| `Method` | `EulerSum` | `EulerSum` or `NIntegrate`. |
| `Scale` | `1` | EulerSum: step size / complex direction.  NIntegrate: contour radius. |
| `Terms` | `7` | EulerSum extrapolation depth. |
| `WorkingPrecision` | `MachinePrecision` | machine `double` or, for a digit count, MPFR. |
| `PrecisionGoal` | `Automatic` | target accuracy (NIntegrate). |
| `MaxRecursion` | `10` | max contour refinements (NIntegrate). |

| Message | When |
|---------|------|
| `ND::ivar` | The second argument is not `x` or `{x, n}` with `x` a symbol. |
| `ND::nnum` | `x0`/`Scale` is not numeric, or `expr` did not evaluate to a number at a sample point. |
| `ND::ord`  | `Method -> EulerSum` was given a non-integer order (use `NIntegrate`). |
| `ND::badscl` | `Scale` did not evaluate to a nonzero number. |
| `ND::badopt` / `ND::badmeth` | Invalid option value / unsupported `Method`. |

### Examples

```mathematica
In[1]:= ND[Exp[x], x, 1]
Out[1]= 2.71828

In[2]:= ND[Cos[x]^3, {x, 2}, 0]
Out[2]= -3.

In[3]:= ND[Sin[x], x, Pi I]
Out[3]= 11.592 + 1.32527*10^-10 I

In[4]:= ND[{Exp[x], Sin[x]}, x, 1]
Out[4]= {2.71828, 0.540302}

In[5]:= ND[Re[Cos[I y]], y, 1]          (* non-analytic: use EulerSum *)
Out[5]= 1.1752

In[6]:= ND[Abs[x], {x, 1}, 0, Scale -> 1 + I]
Out[6]= 0.707107 - 0.707107 I

In[7]:= ND[Sin[100 x], x, 0, Scale -> 1/100]
Out[7]= 100.

In[8]:= ND[Exp[x^2], {x, 4}, 0, Method -> NIntegrate]
Out[8]= 12. - 3.3723*10^-15 I

In[9]:= ND[x, {x, -1/2}, 1, Method -> NIntegrate]   (* = 4/(3 Sqrt[Pi]) *)
Out[9]= 0.752253

In[10]:= ND[Sin[x^2], {x, 3}, 1, Terms -> 20, WorkingPrecision -> 40]
Out[10]= -14.420070264639875819037588981065446865125
```


## NSeries

Numerical Taylor/Laurent series.  `NSeries[f, {x, x0, n}]` gives a numerical
approximation to the series expansion of `f` about `x = x0`, including the
terms `(x - x0)^-n` through `(x - x0)^n`, as a `SeriesData` object.
Implemented natively in C in `src/numerical_calculus/nseries.{c,h}`.
Attribute: `Protected`.

Unlike the symbolic `Series` (which needs a power series at `x0`), `NSeries`
needs only to **sample `f` numerically**, so it works for functions whose
coefficients have no closed form and for **Laurent expansions about essential
singularities** such as `Sin[x + 1/x]`.  It cannot tell a tiny spurious
residual from a true zero -- `Chop` the result when needed -- and returns an
incorrect value if the disk centred at `x0` contains a branch cut of `f`.
For a Laurent result, the `SeriesData` neglects higher-order poles.

### Method

`f` is sampled at `N` equispaced points on a circle of radius `r` centred at
`x0`, and a discrete Fourier transform of the samples recovers the Laurent
coefficients via Cauchy's integral formula (Lyness & Sande, 1971;
Bornemann, *FoCM* 2011):

```
z_j = x0 + r e^{2 pi i j / N},  j = 0 .. N-1,
c_k = (1/N) sum_j f(z_j) e^{-2 pi i j k / N},     (DFT of the samples)
a_e = c_{e mod N} * r^{-e}        for e = -n .. n  (coeff of (x-x0)^e).
```

The upper-half DFT bins (`k = N - m`) supply the **negative**-power
coefficients, so one transform yields both the principal part and the
analytic part; this is exact when `f` is analytic on an annulus containing
the circle.  The sample count is a power of two with an oversampling margin,
`N = 2^{ceil(log2 n) + 2}`, which pushes the leading aliased term below the
round-off floor.  A direct `O(N^2)` DFT is used (no FFT dependency): `N` is
small and each sample requires a full symbolic evaluation of `f`, which
dominates the runtime.  The same path serves machine (`double _Complex`) and
arbitrary-precision (MPFR) computations.

### Result

`SeriesData[x, x0, {a_-n, ..., a_n}, -n, n+1, 1]` -- the coefficient list runs
from exponent `-n` upward, with an `O[(x-x0)^{n+1}]` term.  Coefficients are
real (`Real`/MPFR) when their imaginary part is exactly zero, else
`Complex[re, im]`.

### Options

| Option | Default | Meaning |
|--------|---------|---------|
| `Radius`           | `1`             | Radius of the sampled circle; picks the annulus within which a Laurent series converges. |
| `WorkingPrecision` | `MachinePrecision` | Machine doubles, or MPFR at the requested decimal precision (also shrinks spurious imaginary residuals). |

### Diagnostics (stderr)

| Tag              | Triggered when |
|------------------|----------------|
| `NSeries::ivar`  | Second argument is not a `{x, x0, n}` list, the variable is not a symbol, or `n` is not a non-negative integer. |
| `NSeries::nnum`  | `x0` is not numeric, or `f` did not evaluate to a number on the contour. |
| `NSeries::badopt`| Invalid option value or unrecognised option. |

### Examples

```mathematica
In[1]:= NSeries[Exp[x], {x, 0, 5}] // Chop
Out[1]= 1. + x + 0.5 x^2 + 0.166667 x^3 + 0.0416667 x^4 + 0.00833333 x^5 + O[x]^6

In[2]:= NSeries[Exp[x], {x, I, 5}] // Chop
Out[2]= (0.540302 + 0.841471 I) + (0.540302 + 0.841471 I) (x - I) + ... + O[x - I]^6

In[3]:= NSeries[Sin[x + 1/x], {x, 0, 10}] // Chop
Out[3]= 2.49234*10^-6/x^9 - 0.000174944/x^7 + ... + 0.576725/x + 0.576725 x - ... + O[x]^11

In[4]:= NSeries[1/((1 + x) (3 + x)), {x, 0, 10}, Radius -> 5] // Chop
Out[4]= 9841./x^10 - 3280./x^9 + 1093./x^8 - 364./x^7 + 121./x^6 - 40./x^5 + 13./x^4 - 4./x^3 + 1/x^2 + O[x]^11

In[5]:= NSeries[Exp[x], {x, 0, 5}, WorkingPrecision -> 30] // Chop
Out[5]= 1. + x + 0.5 x^2 + 0.16666666666666666... x^3 + ... + O[x]^6
```

### Notes

- Unlike Mathematica's `NSeries` (whose documented count is
  `2^{ceil(log2 n) + 1}`), Mathilda oversamples by `+2` for a wider
  anti-aliasing margin.
- `Series[Sin[x + 1/x], {x, 0, 10}]` returns unevaluated (no power series at
  the essential singularity); `NSeries` recovers the Laurent expansion.


## NLimit

Numerical limit.  `NLimit[expr, z -> z0]` numerically finds the limiting value
of `expr` as `z` approaches `z0`.  Implemented natively in C in
`src/numerical_calculus/nlimit.{c,h}`.  Attribute: `Protected`.  Like the other
numerical-calculus builtins, `NLimit` cannot tell a tiny spurious residual from
a true zero -- `Chop` the result when needed.

`NLimit` constructs a geometric sequence of sample points approaching `z0` and
recovers the limit by sequence acceleration:

- **Finite `z0`:** samples `z_k = z0 - d * Scale * 2^-k`, where `d` is the
  (unit) `Direction` vector.  The points lie on the `-d` side of `z0`, so one
  moves *along* `d` to reach it.
- **Infinite `z0`** (`Infinity`, `-Infinity`, `I Infinity`,
  `DirectedInfinity[d]`, `ComplexInfinity`): samples march outward on the
  point's ray from the origin, `z_k = u * Scale * 2^k`.

### Methods

- **`EulerSum`** (default) -- Richardson / Romberg extrapolation of the sample
  sequence, using the all-powers denominator `2^j - 1` (the same convention as
  `ND`'s `EulerSum`).  Best for smooth power-series approaches; depth is set by
  `Terms`.
- **`SequenceLimit`** -- Wynn's epsilon algorithm (iterated Shanks transform).
  Exact in one step for a geometric / exponential tail; the number of
  iterations is `WynnDegree` (which needs at least `2(WynnDegree + 1)` terms).
  The estimate is read from the `ε_{2·WynnDegree}` column at the entry that best
  agrees with its neighbour (avoiding the roundoff-amplified bottom corner).

### Options

| Option | Default | Meaning |
|--------|---------|---------|
| `Method` | `EulerSum` | `EulerSum` or `SequenceLimit`. |
| `WorkingPrecision` | `MachinePrecision` | `MachinePrecision`, or digits → MPFR. |
| `Direction` | `Automatic` (≡ `-1`) | complex approach vector for finite `z0`. |
| `Scale` | `1` | initial step (finite) / distance from origin (infinite). |
| `Terms` | `7` | number of sample points / extrapolation depth. |
| `WynnDegree` | `1` | `SequenceLimit` iterations. |

`Direction -> Automatic` (`-1`) approaches a finite point from larger values;
`Direction -> 1` from smaller; complex rays such as `-Exp[225 Degree I]` select
an arbitrary direction in the complex plane (essential for path-dependent
limits and branch cuts).

### Robustness

The last two extrapolates are compared against the *sample scale* (the largest
`|S_k|`).  A divergent or non-settling sequence -- e.g. a power-law approach to
infinity like `NLimit[1/x, x -> 0]` -- yields `NLimit::noise` and the form is
returned unevaluated.  A bounded but only roughly-resolved limit is returned
rather than refused, matching Mathematica.

### Beyond / unlike Mathematica's NLimit

- `EulerSum` is implemented as Richardson/Romberg extrapolation (consistent
  with `ND`), not Euler series summation; results agree on smooth cases.
- Some sequences that make Mathematica's `EulerSum` fail (e.g.
  `NLimit[Tanh[x], x -> Infinity]` at the default `Terms`) succeed here because
  Richardson does not divide by a vanishing Euler ratio.

### Messages

| Message | When |
|---------|------|
| `NLimit::noise` | Cannot recognise a limiting value (divergent / noisy sequence). |
| `NLimit::notnum` | `expr` not numerical at a sample point, or the point/Direction/Scale is not numerical. |
| `NLimit::ndterm` | Not enough `Terms` for the chosen `Method`. |
| `NLimit::badopt` / `NLimit::badmeth` | Invalid option value / unsupported `Method`. |

### Examples

```mathematica
In[1]:= NLimit[Sin[x]/x, x -> 0]
Out[1]= 1.

In[2]:= NLimit[(1 + 1/n)^n, n -> Infinity]
Out[2]= 2.71828

In[3]:= NLimit[(1 + I/x)^x, x -> Infinity]
Out[3]= 0.540302 + 0.841471 I

In[4]:= NLimit[Tanh[Pi x]/(1 + x^2), x -> I] // Chop
Out[4]= -1.5708 I

In[5]:= NLimit[(10^x - 1)/x, x -> 0, Terms -> 10, Method -> SequenceLimit]
Out[5]= 2.30259                                   (= Log[10])

In[6]:= NLimit[z + Conjugate[z]/z, z -> 0, Direction -> -I] // Chop
Out[6]= -1.

In[7]:= NLimit[Tan[z], z -> Infinity I, Method -> SequenceLimit] // Chop
Out[7]= 1. I

In[8]:= NLimit[(2^x - 1)/x, x -> 0, WorkingPrecision -> 30, Terms -> 14]
Out[8]= 0.693147180559945309417232121458   (= Log[2])

In[9]:= NLimit[1/x, x -> 0]
        NLimit::noise: Cannot recognize a limiting value. ...
Out[9]= NLimit[1/x, x -> 0]
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
