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
- Three-stage dispatch cascade (2026-05): `Integrate[f, x]` (Method ->
  Automatic, default) tries each subroutine in order and returns the
  first non-`NULL` result:
  1. `Integrate\`BronsteinRational[f, x]` — when `PolynomialQ[f, x] ||
     rationalQ[f, x]`.
  2. `Integrate\`RischNorman[f, x]` — Bronstein pmint, all integrands.
  3. `Integrate\`CRCTable[f, x]` — CRC integral table lookup (lazy-loaded
     from `src/internal/CRCMathTablesIntegrals.m` on first call).
  If every stage gives up the call bubbles back unevaluated.
- `Method -> "<name>"` option (3rd argument) bypasses the cascade and
  dispatches strictly to a single subroutine, with no fallback:
  - `"Automatic"` — default cascade above.
  - `"BronsteinRational"` — `Integrate\`BronsteinRational[f, x]`.
  - `"RischNorman"` — `Integrate\`RischNorman[f, x]`.
  - `"CRCTable"` — `Integrate\`CRCTable[f, x]`.
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
`Integrate`RischNorman` (Bronstein pmint), `Integrate`CRCTable`
(table lookup), and the unit-test helpers `Integrate`Helpers`Content`,
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
- `FindMinimum[f, {x, y, ...}]` -- n-D auto-start at zero.
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
