# Computing Limits in Mathilda

## Introduction

`Limit` is one of the few builtins in Mathilda where correctness is not
local: you cannot decide the answer by looking at one node of the
expression tree. `Sin[x]/x` at `x -> 0` is indeterminate by
substitution, finite in value, and undetectable without some theory of
growth rates. A computer algebra system that ships `Limit` ships an
opinion about how all of calculus fits together.

This post is a tour of that opinion in Mathilda. We will:

1. set the scene with a couple of easy cases,
2. walk through the dispatcher layer by layer, using actual entries
   from `tests/test_limit.c` as worked examples,
3. compare the resulting architecture with Maxima's `limit.lisp`.

The implementation lives in `src/limit.c`
(~1640 lines) and the test battery in `tests/test_limit.c`
(16 test groups, ~70 assertions).

## Two easy examples

```text
In[1]:= Limit[x^2, x -> 2]
Out[1]= 4

In[2]:= Limit[Sin[x]/x, x -> 0]
Out[2]= 1
```

The first is continuous substitution: plug `x = 2` into `x^2`. The
second is everything `Limit` is allowed to do that substitution is not.

Mathilda treats the two cases with different machinery. The first falls
out of a fast path that checks whether the expression is free of the
variable at the target point. The second triggers a Taylor series
expansion, reads off the leading term, and hands it back. Every other
limit in the test suite is one of those two with an additional
preprocessing pass in front, or a different structural recognizer.

## The dispatcher

`compute_limit` in `src/limit.c` is a cascade of specialists.
Each is a function of signature

```c
static Expr* layer_xxx(Expr* f, LimitCtx* ctx);
```

that either returns a freshly allocated result (done, stop cascading)
or `NULL` (pass to the next layer). The order encodes the priority we
want at each shape: cheap and exact first, speculative and expensive
last.

```
  rewrite_reciprocal_trig   (pre-pass: Csc/Sec/Cot -> Sin/Cos)
      |
      v
  layer1_fast_paths           (free-of-x, variable-itself,
                               continuous substitution)
  layer_arctan_infinity       (ArcTan/ArcCot of divergent arg)
  layer3_rational             (polynomial / polynomial)
  layer_log_of_finite         (Log[g] with finite inner limit)
  layer_log_merge             (Log[..] + polynomial at infinity)
  layer_plus_termwise         (Plus at infinity with finite terms)
  layer5_log_reduction        (f^g indeterminate forms)
  layer_bounded_envelope      (squeeze for Sin / Cos / ... decay)
  layer2_series               (Series-based, the workhorse)
  layer5_lhospital            (L'Hospital with growth guard)
  layer6_bounded              (Interval[{-1,1}] for Sin[1/x])
```

The rest of this post walks the layers in the same order, with the
test cases each one is responsible for.

---

### Reciprocal-trig rewrite (pre-pass)

```c
/* Normalize reciprocal trig up-front. */
Expr* rewritten = rewrite_reciprocal_trig(f_in);
Expr* f = simp(rewritten);
```

Mathilda's evaluator aggressively folds `0 * ComplexInfinity -> 0`. That
is the correct algebraic identity when the `ComplexInfinity` is really
a pole of fixed order, but it is a trap for symbolic limit computation:
`Csc[x]` at `x = 0` evaluates to `ComplexInfinity`, so `x Csc[x]` would
fold to `0` by substitution before a limit routine has a chance to
think.

We avoid the trap by rewriting `Csc`, `Sec`, `Cot`, `Csch`, `Sech`,
`Coth` into their `Sin` / `Cos` / `Sinh` / `Cosh` forms **before** the
dispatcher runs. Eighteen entries in `test_reciprocal_trig` exercise
this path, including:

```text
Limit[x Csc[x], x -> 0]                  = 1
Limit[Cot[4 x] Tan[3 x], x -> 0]         = 3/4
Limit[Csc[Pi x] Log[x], x -> 1]          = -1/Pi
Limit[(x - ArcSin[x]) Csc[x]^3, x -> 0]  = -1/6
Limit[Sec[2 x] (1 - Tan[x]), x -> Pi/4]  = 1
```

After the rewrite, `x Csc[x]` becomes `x/Sin[x]`, which the Series
layer dispatches to `1` in one step.

### Layer 1 -- Fast paths and continuous substitution

Three cheap checks, in order:

1. **Free of `x`.** If the expression does not contain the variable, it
   is its own limit. `Limit[7, x -> 0]` is `7`; `Limit[a, x -> 5]` is
   `a`.
2. **The variable itself.** `Limit[x, x -> 5]` is `5`.
3. **Continuous substitution.** `Together`-normalize `f` into `N/D`,
   evaluate `D(point)`; if non-zero and non-divergent, plug `point`
   into both `N` and `D` and divide.

```c
Expr* tog = simp(mk_fn1("Together", expr_copy(f)));
Expr* num = simp(mk_fn1("Numerator",   expr_copy(tog)));
Expr* den = simp(mk_fn1("Denominator", expr_copy(tog)));
Expr* den_at = subst_eval(den, ctx->x, ctx->point);
bool den_bad = is_lit_zero(den_at) || is_divergent(den_at);
```

The `is_lit_zero` check is enough because Mathilda's `Power` evaluator
folds `Power[0, b] -> 0` for any positive `b` and
`Power[0, b] -> ComplexInfinity` for negative `b`. `Sqrt[x-1]/x` at
`x = 1` substitutes to `Sqrt[0]/1`, which collapses to `0/1 = 0`
automatically. Earlier versions of the limit layer carried a
`reduces_to_zero` predicate to handle `Sqrt[0]`; moving the fold into
`builtin_power` removed the need.

Worked examples (from `test_fast_paths`, `test_power_zero_folding`):

```text
Limit[7, x -> 0]                        = 7          (free-of-x)
Limit[x, x -> 5]                        = 5          (variable itself)
Limit[x^2, x -> 2]                      = 4          (substitution)
Limit[a x^2 + 4 x + 7, x -> 1]          = 11 + a     (substitution)
Limit[Sqrt[x - 1]/x, x -> 1]            = 0          (Power[0, 1/2] fold)
```

Multivariate limits `Limit[f, {x, y} -> {a, b}]` run through a
specialised variant of the continuous-substitution check; iterated
limits `Limit[f, {x -> a, y -> b}]` recurse right-to-left.

### Layer -- ArcTan / ArcCot at infinity

Series cannot see through `ArcTan[h(x)]` when `h(x) -> oo`. We
intercept that shape early: compute `Limit[h(x)]`, and if it diverges
to `+/- oo`, return `+/- Pi/2` (for `ArcTan`) or `0 / Pi` (for
`ArcCot`).

```text
Limit[ArcTan[x^2 - x^4], x -> Infinity] = -Pi/2
```

Internally: `Limit[x^2 - x^4] = -Infinity` via the rational-function
layer, so the outer `ArcTan` limit is `-Pi/2`.

### Layer 3 -- Rational functions

For `P(x) / Q(x)` with both sides polynomial in `x`, we extract degrees
and leading coefficients via `CoefficientList` and decide in one shot:

- At `+/- Infinity`: compare degrees, use the leading-coefficient
  ratio when they are equal, and pick a sign of infinity from the
  parity of `deg(P) - deg(Q)` when unequal.
- At a finite point: substitute, handle `0/0` by deferring to a later
  layer, handle `nonzero/0` as a directional infinity.

```text
Limit[x/(x + 1), x -> Infinity]         = 1
Limit[3 x^2/(x^2 - 2), x -> -Infinity]  = 3
Limit[x/(x - 2 a), x -> a]              = -1
```

### Layer -- Log of a finite expression

If `f = Log[g]` and the inner limit of `g` is a finite non-zero `c`,
return `Log[c]`. One-liner, but load-bearing: it turns

```text
Log[(2 + E^x) Exp[-x]]  ->  Log[1 + 2 Exp[-x]]  ->  Log[1]  =  0
```

at `x -> oo`, which is the finishing move for the `layer_log_merge`
rewrite just below.

### Layer -- Log merge at infinity

Individually divergent summands can cancel into a finite `Log`. The
canonical shape is

```text
Limit[-x + Log[2 + E^x], x -> Infinity]
```

where `-x -> -oo` and `Log[2 + E^x] -> +oo`. Term-by-term the answer
is `-oo + oo`, which no finite-term dispatcher will resolve. The trick
is to fuse the sum under a single logarithm:

```
 Sum(Log[g_i]) + Sum(h_j)  ->  Log[ Prod(g_i) * Prod(Exp[h_j]) ]
```

The rewrite uses `Exp[h_j]` as a fictitious kernel -- we are just
moving everything inside `Log`, counting on `Expand` to simplify the
argument afterwards. In the example:

```
Log[(2 + E^x) * Exp[-x]]
  = Log[2 Exp[-x] + 1]
  -> Log[1] = 0                (as x -> Infinity)
```

Termination is guarded by a structural-change check
(`expr_eq(new_f, f)` bails out); recursion is bounded by `ctx->depth`.

### Layer -- Plus termwise at infinity

If every summand of a top-level `Plus` has a finite limit individually,
add them up. This is what the previous two layers feed into when they
need a clean "finish". It is also what rescues harmless shapes that
only look divergent because of a dominant term hiding a cancellation
that `Together` did not catch.

### Layer 5.3 -- Log-reduction for exponential indeterminate forms

`f^g` with `f -> 1` and `g -> oo` is the classical `1^oo`
indeterminate form. We use the identity

```
 f^g = Exp[g Log[f]]
```

and compute the inner limit. The layer fires whenever any `Power`
factor at the top level of `f` has the limit variable in its exponent.

```c
static Expr* log_contribution(Expr* factor) {
    if (has_head(factor, "Power") && factor->data.function.arg_count == 2) {
        /* Log[a^b] = b Log[a]. */
        return mk_times(expr_copy(factor->data.function.args[1]),
                        mk_fn1("Log", expr_copy(factor->data.function.args[0])));
    }
    return mk_fn1("Log", expr_copy(factor));
}
```

A guard keeps the layer from firing on shapes like
`(Cos[x]-1)/(Exp[x^2]-1)` where `Exp` is buried inside a `Plus` -- those
are determinate `0/0` forms best left for the Series layer.

Worked examples (from `test_log_reduction` and `test_exp_indeterminate`):

```text
Limit[(1 + A x)^(1/x), x -> 0]                 = E^A
Limit[(1 + a/x)^(b x), x -> Infinity]          = E^(a b)
Limit[Cos[x]^Cot[x], x -> 0]                   = 1
Limit[(1 + Sin[a x])^Cot[b x], x -> 0]         = E^(a/b)
Limit[x^Tan[(Pi x)/2], x -> 1]                 = E^(-2/Pi)
Limit[((-3+2x)/(5+2x))^(1+2x), x -> Infinity]  = 1/E^8
Limit[E^(-x) x^20, x -> Infinity]              = 0
```

The last one is worth unpacking. After log-reduction we are asking for
`Limit[-x + 20 Log[x], x -> Infinity]`, which the `layer_log_merge`
pass rewrites as `Log[Exp[-x] * x^20]`. `Exp[-x] * x^20` on its own has
limit `0` (this happens through the bounded-envelope and series
routines), so `Log[...] = -Infinity`, and `Exp[-Infinity] = 0`.

### Layer -- Bounded envelope (squeeze theorem)

`Sin[t^2]` has no Taylor expansion at infinity, so if the series layer
were to run first it would flail. We intercept shapes where a bounded
oscillatory head (currently `Sin`, `Cos`, and relatives) is multiplied
by a factor that decays, producing `|f|` times a bound on the
oscillator's magnitude. If the bound's limit is `0`, the whole
expression limit is `0`.

```text
Limit[Sin[t^2]/t^2, t -> Infinity]          = 0
Limit[(1 - Cos[x])/x, x -> Infinity]        = 0
Limit[x^2 Sin[1/x], x -> 0]                 = 0
Limit[(x Sin[x])/(5 + x^2), x -> Infinity]  = 0
```

The layer is currently gated to `+Infinity` because the
magnitude-upper-bound utility returns `x` (not `Abs[x]`) for the bare
variable, which is only an upper bound when `x` is positive.
Negative-infinity cases would need a preparatory `x -> -y` substitution.

### Layer 2 -- Series-based evaluation

This is the workhorse. We call `Series[f, {x, x0, k}]` with increasing
order `k` (starting at 4, doubling up to 64) and read off the leading
term from the resulting `SeriesData`. The three outcomes:

| leading exponent | limit |
|-----------------:|:------|
| `> 0`            | `0`   |
| `= 0`            | leading coefficient |
| `< 0`            | directional infinity (sign from the leading coefficient and the direction) |

Worked example: `Limit[Sin[x]/x, x -> 0]`.

```
Series[Sin[x]/x, {x, 0, 4}]
  = 1 - x^2/6 + x^4/120 + O[x]^5
```

The leading exponent is `0`; the leading coefficient is `1`; done.

Another: `Limit[(Tan[x] - x)/x^3, x -> 0]`.

```
Series[(Tan[x] - x)/x^3, {x, 0, 4}]
  = 1/3 + (2 x^2)/15 + O[x]^3
```

The leading exponent is `0`; the leading coefficient is `1/3`.

At `-Infinity` we first substitute `x -> -y` so Series can use its
`+Infinity` kernel. Tests covered by this layer:

```text
Limit[Sin[x]/x, x -> 0]                  = 1
Limit[(Cos[x]-1)/(Exp[x^2]-1), x -> 0]   = -1/2
Limit[(Tan[x]-x)/x^3, x -> 0]            = 1/3
Limit[Sin[2 x]/Sin[x], x -> Pi]          = -2
Limit[(x^a - 1)/a, a -> 0]               = Log[x]
Limit[Sqrt[x^2 + a x] - Sqrt[x^2 + b x], x -> Infinity] = (a - b)/2
```

The last, `Sqrt[x^2 + a x] - Sqrt[x^2 + b x]`, is an "RP-form"
(radical-polynomial) family where Series at infinity expands in
`1/x`:

```
Series[Sqrt[x^2 + a x] - Sqrt[x^2 + b x], {x, Infinity, 2}]
  = (a - b)/2 + (b^2 - a^2)/(8 x) + O[1/x]^2
```

Leading term `(a - b)/2`, exponent `0`: that is the answer.

### Layer 5.1 -- L'Hospital with a growth guard

If Series cannot see through the expression (typical symptoms: leading
coefficients all literally zero up to a large order, or Series
returning itself), we try L'Hospital. We insist on a `N/D` split
(non-trivial denominator), and we insist on the numerator and
denominator being simultaneously zero (or simultaneously divergent) at
the point; otherwise the rule is not applicable.

A growth guard aborts if the leaf count of `N' / D'` grows in three
consecutive iterations without the indeterminacy resolving. This is
the escape hatch for pathological inputs -- L'Hospital on the wrong
shape can generate expressions that blow up in size faster than they
simplify.

### Layer 6 -- Bounded-oscillation interval returns

`Limit[Sin[1/x], x -> 0]` has no value in the ordinary sense; we
return `Interval[{-1, 1}]`. The layer fires for `Sin` or `Cos` of an
expression that diverges at the target point.

```text
Limit[Sin[1/x], x -> 0]       = Interval[{-1, 1}]
Limit[Sin[x], x -> Infinity]  = Interval[{-1, 1}]
```

---

## Comparison with Maxima

Maxima's limit implementation lives in
[`src/limit.lisp`](https://github.com/calyau/maxima/blob/master/src/limit.lisp)
(a ~4300-line Common Lisp file) with an auxiliary `src/gruntz.lisp`
for the Gruntz algorithm on real-valued limits at infinity. The shapes
of the two systems are worth putting side by side.

### Architectural similarities

- **Dispatcher on expression shape.** Maxima's `limit` entry point
  dispatches to per-head handlers (`simplim%sin`, `simplim%log`,
  `simplim%atan`, ...) via the `SIMPLIM%` property on symbols; Mathilda
  cascades through `layer_*` functions. Both amount to a dictionary
  keyed by the head of the expression.
- **Series/Taylor as the primary tool for indeterminate forms.**
  Maxima has a dedicated `taylim` / `tlimit` mode that replaces the
  default limit algorithm with "take a Taylor series and read off the
  leading term." Mathilda's `layer2_series` is exactly this, hooked in
  as the main engine rather than an alternative mode.
- **L'Hospital as a fallback with a growth guard.** Maxima uses
  `lhospitallim` (default 4) as the maximum number of L'Hospital
  applications; Mathilda uses an analogous iteration cap plus a leaf-
  count growth detector.
- **Special-case logarithmic reduction for `f^g`.** Both rewrite `f^g`
  as `Exp[g Log[f]]` when the form is indeterminate. Maxima calls this
  in `exp-lim`; Mathilda calls it `layer5_log_reduction`.
- **Bounded-oscillation / squeeze handling for `Sin`, `Cos`.** Maxima
  returns `ind` (for "indeterminate but bounded") when it can only say
  the limit is in an interval; Mathilda returns `Interval[{-1, 1}]`.
- **Finite-point 0/0 and infinity/infinity both use the same
  pipeline.** Neither system has a separate algorithm for limits at
  `Infinity`; we substitute the point into the series kernel.

### Architectural differences

- **No Gruntz algorithm.** Maxima's `gruntz.lisp` implements the
  Gruntz "MrvLimit" algorithm for limits of real-valued exp-log
  functions at infinity. It handles arbitrarily-nested iterated
  exponentials like `Exp[Exp[x] - Exp[x + 1/x]]` that our Series-plus-
  log-merge combination cannot reach. This is the biggest capability
  gap.
- **Sign/assumption machinery.** Maxima uses its global `asksign` and
  assumption database to resolve sign ambiguities during limit
  computation ("Is `a` positive, negative, or zero?"). Mathilda has no
  assumption database; our sign decisions are all literal
  (`literal_sign` checks structural sign only), and we conservatively
  return `DirectedInfinity[c]` or leave the expression unevaluated
  when the sign of a leading coefficient is symbolic.
- **Direction convention.** Maxima's third argument is a plus/minus
  symbol; Mathilda accepts both the Mathematica strings
  (`"FromAbove"`, `"FromBelow"`) and Mathematica integer direction
  codes (`-1` meaning "approaching from above", `+1` from below).
  Both are tested.
- **Log-merge preprocessing.** Maxima does not have a direct analogue
  of Mathilda's `layer_log_merge`; it would typically handle
  `-x + Log[2 + E^x]` through the Gruntz layer's exp-log rewriting.
  Our merge pass is a simpler, targeted rewrite for the finite
  combinations of `Plus`, `Log`, and polynomial summands.
- **Radical simplification.** Maxima has `radcan` as a user-callable
  simplifier; Mathilda runs a narrow version of radical fusion inside
  `builtin_times` itself, so that `Sqrt[6]/Sqrt[2]` becomes `Sqrt[3]`
  system-wide -- not only when it appears as `Limit` output.
- **Scope.** Maxima's `limit.lisp` is a mature file with dozens of
  special-function handlers (`%bessel_j`, `%gamma`, `%zeta`, etc.).
  Mathilda currently implements only the elementary-function
  subsection.

### Capability snapshot

| Shape                              | Mathilda        | Maxima |
|:-----------------------------------|:---------------|:-------|
| `P(x)/Q(x)` at `+/- Infinity`      | yes (Layer 3)  | yes    |
| Taylor-leading-term at a finite point | yes (Layer 2) | yes    |
| `1^Infinity`, `0^0`, `Infinity^0`  | yes (Layer 5.3)| yes    |
| L'Hospital fallback                | yes (Layer 5.1)| yes    |
| `Sin[1/x]` -> `Interval`           | yes (Layer 6)  | yes (`ind`) |
| Squeeze for `Sin[h(x)]/p(x)`       | yes            | partial|
| Gruntz nested exp-log at infinity  | no             | yes (`gruntz.lisp`) |
| Assumption-driven sign resolution  | no             | yes    |
| Special functions (Gamma, Bessel)  | no             | yes    |

## Recent additions: nine work packages

After the initial dispatcher stabilised, a single REPL report turned
up roughly thirty more failing shapes. Rather than patch each in
isolation, we grouped them into nine thematic work packages. Each
landed as an independent module with its own `test_wp*` in
`tests/test_limit.c`, so a future regression points at exactly one
work area.

Below is a tour of each package with the motivating example it was
written to fix.

### WP-8 -- Numeric-point fast path

A plain numeric limit point where the expression is analytic should
resolve by substitution, regardless of how large the expression is.
The pre-batch dispatcher would sometimes ride through Series or
L'Hospital on such inputs and spin indefinitely on heavy-but-analytic
shapes.

`layer1_fast_paths` now calls `try_numeric_point_substitution`, which
plugs the numeric point into the `Together`-normalised form after
two gates: the denominator must not vanish at the point, and no
`Power[base, x-dependent exponent]` with a divergent exponent may
remain (otherwise we'd short-cut an `1^oo` / `0^0` / `oo^0` seed).

```text
Limit[Log[1 - (Log[Exp[z]/z - 1] + Log[z])/z]/z, z -> 100]
  -> Log[1 + (-Log[100] - Log[-1 + E^100/100])/100]/100
```

Previously a hang; now a one-shot substitution. The guards keep
`Limit[Sin[x]/x, x -> 0]` (denominator vanishes) and
`Limit[(1 + A x)^(1/x), x -> 0]` (divergent exponent) on the Series
and log-reduction paths where they belong.

### WP-2 -- b^x kernel: series inversion with symbolic-Log coefficients

`Limit[(25^x - 5^x)/(4^x - 2^x), x -> 0]` used to hang inside
`so_inv`. Each iteration of the Laurent-series inverse multiplied
growing polynomials in `Log[a]`, `Log[b]`, `Log[2]`, `Log[5]`, and
Mathilda's evaluator does not canonicalise those into a unique form,
so the iteration count is unbounded.

The fix is a data-dependent iteration cap in `so_inv`, sized by the
total leaf count of the input coefficients:

```c
if (total_leaves <= 4)        N_cap = 64;
else if (total_leaves <= 20)  N_cap = 16;
else if (total_leaves <= 60)  N_cap = 6;
else                           N_cap = 4;
```

Trivial numeric inputs retain the full order; heavy symbolic inputs
still produce a valid leading-term Laurent -- which is all `Limit`
needs. A single-`simp`-per-iteration rewrite in the loop cuts the
per-step cost by roughly 3x on top.

```text
Limit[(25^x - 5^x)/(4^x - 2^x), x -> 0] = (Log[25] - Log[5])/(Log[4] - Log[2])
Limit[(3^x - 2^x)/(5^x - 2^x), x -> 0]  = (Log[3] - Log[2])/(Log[5] - Log[2])
Limit[(a^x - 1)/x, x -> 0]              = Log[a]
```

### WP-4 -- Pole sign-disagreement classifier

Mathematica's `Direction -> Reals` asks for the real-line answer at
a pole. If the two one-sided limits disagree in sign (odd-order real
pole, `Sqrt`-type branch singularity, `Tan` at `Pi/2`) the real-line
limit does not exist, so the answer must be `Indeterminate`.
`Direction -> Complexes` asks for the radial / all-direction answer:
any isolated pole maps to `ComplexInfinity` regardless of parity.

A new `LIMIT_DIR_REALS` tag (distinct from the implicit two-sided
default) threads this through the layer cascade:

```text
Limit[Tan[x], x -> Pi/2]                         = ComplexInfinity
Limit[Tan[x], x -> Pi/2, Direction -> Reals]     = Indeterminate
Limit[Tan[x], x -> Pi/2, Direction -> Complexes] = ComplexInfinity
Limit[(x - 3)/(Sqrt[x] - 3), x -> 9, Direction -> Reals] = Indeterminate
Limit[Sin[Tan[x]], x -> Pi/2, Direction -> Reals]        = Indeterminate
Limit[1/(x-2)^3, x -> 2, Direction -> Reals]             = Indeterminate
```

The implicit two-sided default keeps the old `ComplexInfinity`
fall-back for unqualified rational-function limits: this matches
our internal tests and the older Mathematica convention.

### WP-1 -- Multivariate path-dependence analyzer

Joint limits `{x1, ..., xn} -> {a1, ..., an}` now run a polar (2D)
or spherical (3D) substitution, then a direction-sampling
cross-check along the axes and main diagonals. Three outcomes:

| polar r-limit | behavior |
|:--------------|:---------|
| free of angles | return it |
| depends on angles | sample paths; agree -> common value, disagree -> `Indeterminate` |
| diverges | sample paths to classify |

The origin fast path also scans for any reciprocal subterm whose base
vanishes at the joint point, which catches `0/0` buried inside
`ArcTan`, `Sin`, etc., that Mathilda's scalar evaluator would
otherwise fold to `0`.

```text
Limit[Tan[x y]/(x y), {x, y} -> {0, 0}]               = 1
Limit[y/(x + y), {x, y} -> {0, 0}]                    = Indeterminate
Limit[(x^3 + y^3)/(x^2 + y^2), {x, y} -> {0, 0}]      = 0
Limit[ArcTan[y^2/(x^2 + x^3)], {x, y} -> {0, 0}]      = Indeterminate
Limit[(x z)/(x^2 + y^2 + z^2), {x, y, z} -> {0, 0, 0}] = Indeterminate
Limit[ArcTan[y/x], {x, y} -> {Infinity, Infinity}]    = Indeterminate
Limit[x y/(x^2 + y^2 + 1), {x, y} -> {0, 0}]          = 0
```

### WP-3 -- Series of `x^a` at a nonzero finite point

Two coupled bugs conspired to make `Series[x^a, {x, 1, 2}]` collapse
to `0`:

1. `try_factor_power_prefactor` pulled `x^a` out of any top-level
   `Times` regardless of `x0`. At `x0 = 1` the remaining factors
   expanded in isolation and evaluated to `0`.
2. `try_apart_preprocess` ran partial fractions on expressions with
   `Power[base, non-rational]`, and our `Apart` collapses such inputs
   to `0`.

The fix gates prefactor extraction to `x0 = 0` or `Infinity` only,
and teaches `try_apart_preprocess` to refuse non-rational powers.
Together these correctly expand `x^a` at `x0 = 1` via the
`(1 + u)^alpha` binomial kernel:

```text
Normal[Series[x^a, {x, 1, 2}]]
  = 1 + a (x - 1) + (1/2) a (a - 1) (x - 1)^2

Limit[3 (x^a - a x + a - 1)/(x - 1)^2, x -> 1] = (3/2) a (a - 1)
```

### WP-6 -- Sinh / Cosh / Tanh exponentialisation at infinity

A pre-Series tree walk (`rewrite_hyperbolic_to_exp`) replaces

```
Sinh[z] -> (E^z - E^-z)/2
Cosh[z] -> (E^z + E^-z)/2
Tanh[z] -> (E^z - E^-z)/(E^z + E^-z)
```

whenever the limit point is `+/- Infinity`, and hands the result to
`Expand`. After expansion the term-wise `Plus` layer folds
cancelling `Exp[k x]` summands that the Sinh/Cosh series kernels
themselves couldn't see:

```text
Limit[(1 + Sinh[x])/Exp[x], x -> Infinity]  = 1/2
Limit[Cosh[x]/Exp[x], x -> Infinity]        = 1/2
Limit[Cosh[x]/Exp[-x], x -> -Infinity]      = 1/2
```

The rewrite is gated to infinite limit points only: at a finite
point the Taylor expansion of `Sinh`/`Cosh` is already fine.

### WP-5 -- Sign-aware envelope / dominant term at infinity

`layer_plus_termwise` gained a `growth_exponent_upper` helper: a
polynomial-degree upper bound that

- treats `Sin`, `Cos`, `Tanh`, `ArcTan`, `ArcCot` as `0` (bounded),
- multiplies across `Times`,
- takes max over `Plus`,
- multiplies by `n` in `Power[base, n]` for nonnegative integer `n`,
- treats `Power[base, negative_int]` with a growing base as `0`
  (reciprocal of something diverging is bounded).

When one `Plus` summand has strictly larger growth than every other
and that summand's limit is `+/- Infinity`, it wins:

```text
Limit[x^2 + x Sin[x^2], x -> Infinity]  = Infinity     (x^2 absorbs x*Sin)
Limit[x + Sin[x], x -> Infinity]        = Infinity     (x absorbs Sin)
Limit[x - Cos[x], x -> Infinity]        = Infinity
Limit[x + Sin[x], x -> -Infinity]       = -Infinity
Limit[Sin[x^2] + Cos[x], x -> Infinity] = (unevaluated)  (no dominant term)
```

### WP-9 -- Complex-direction branch cuts

A new `LIMIT_DIR_IMAGINARY` tag handles numeric-imaginary
directions (`I`, `k I`, `Complex[0, k > 0]`). The analytic layers
compute the principal-branch result as usual; a post-pass in
`builtin_limit` then

- conjugates the imaginary part via `ReplaceAll[I -> -I]` when the
  direction was `I` (we're landing below the cut), and
- returns `Indeterminate` for `Direction -> Complexes` when the
  principal result carries an imaginary part (the point sits on a
  branch cut and the radial approach is two-valued).

Any signed-infinity pole under `Complexes` collapses to
`ComplexInfinity`. As a bonus, `Limit[{f1, f2, ...}, ...]` now
threads over a top-level `List` in its first argument:

```text
Limit[Sqrt[z], z -> -1, Direction -> Reals]     = I
Limit[Sqrt[z], z -> -1, Direction -> I]         = -I
Limit[Sqrt[z], z -> -1, Direction -> Complexes] = Indeterminate
Limit[Log[z], z -> -1, Direction -> Reals]      = I Pi
Limit[Log[z], z -> -1, Direction -> I]          = -I Pi
Limit[Log[z], z -> -1, Direction -> Complexes]  = Indeterminate
Limit[{Sqrt[z], Log[z]}, z -> -1, Direction -> I] = {-I, -I Pi}
Limit[1/(x - 2), x -> 2, Direction -> I]        = ComplexInfinity
```

### WP-7 -- Gruntz-lite for `Log[sum]`

A partial, targeted Gruntz-style rewrite: at `+Infinity`, if one
summand of `Log[dom + rest_1 + rest_2 + ...]` strictly dominates
the others (by the WP-5 `growth_exponent_upper` ranking), rewrite
as

```
Log[dom + rest...] -> Log[dom] + Log[1 + rest/dom]
```

and recurse. This unblocks the iterated-log family where the
factored-out `Log[dom]` cancels against an outer term:

```text
Limit[Log[x + Log[x]] - Log[x], x -> Infinity]     = 0
Limit[Log[x^2 + x] - 2 Log[x], x -> Infinity]      = 0
```

This is **not** full Gruntz. Full Hardy-field comparative
asymptotics (stacked `Exp[Exp[-x/(1 + Exp[-x])]]`, multi-level
log-exp dominance, `Sin[x] + Log[x-a]/Log[E^x-E^a]` as `x -> a`,
Gruntz's canonical `Exp[Log[Log[x + Exp[Log[x] Log[Log[x]]]]]/...]
-> E` nested-log case) remain future work. What lands here is the
machinery that the full Gruntz algorithm needs -- the
`growth_exponent_upper` ranking, the `Log[sum]` rewrite, the
hyperbolic exponentialisation -- so the next step is plugging in a
proper MRV (most-rapidly-varying) set construction rather than a
whole new layer.

---

### Capability snapshot, post-batch

| Shape                                       | Status |
|:--------------------------------------------|:-------|
| Heavy analytic expression at numeric point  | WP-8   |
| `(b1^x - b2^x)/(b3^x - b4^x)` at 0          | WP-2   |
| `Direction -> Reals` / `Complexes`          | WP-4   |
| Polar/spherical multivariate at origin      | WP-1   |
| `0/0` buried inside `ArcTan`, etc.          | WP-1   |
| `Series[x^a, {x, nonzero, ...}]`            | WP-3   |
| `Sinh`/`Cosh`/`Tanh` at infinity            | WP-6   |
| Dominant polynomial absorbs bounded oscillator | WP-5 |
| `Sqrt`/`Log` at negative real, `Direction -> I` | WP-9 |
| `Limit[{f1, f2, ...}, ...]` list threading  | WP-9   |
| Simple `Log[dom + rest]` Gruntz rewrite     | WP-7   |
| `Abs[g(x)]/h(x)` at a kink of `g`           | Abs-rewrite |
| Full Gruntz MRV (stacked exp, iterated log) | future work |

Each row maps to its own `test_wp*` in `tests/test_limit.c`. The
test battery now stands at 24 groups and ~130 assertions, all
passing under Valgrind with zero leaks.

---

### Layer -- `Abs[g(x)]` direction-aware rewrite

`Abs` is continuous everywhere on the real line but
non-differentiable at the kink. Both Series and L'Hospital fail
silently on shapes like `Abs[x]/x` near `x = 0`: Series cannot
expand `Abs[x]` (the leading coefficient is `Derivative[1][Abs][0]`,
itself an unevaluated symbolic value) and L'Hospital differentiates
`Abs[x]` into the same `Derivative[1][Abs][x]`, then accepts
`Derivative[1][Abs][0] / 1` as a "clean" answer because the
classifier only checks for zero, infinity, and `Indeterminate` --
nothing tells it that the numerator is a non-numeric symbolic head.

The fix is a dedicated layer that runs *before* Series and
L'Hospital. It does three things, each cheap enough to short-circuit
when it cannot make progress:

1.  **Continuity fast path.** If `f` is exactly `Abs[g(x)]` and
    `Limit[g, x -> a]` resolves, return `Abs` of that value
    (with `Abs[Infinity] = Abs[-Infinity] = Abs[ComplexInfinity] =
    Infinity` folded). This handles `Abs[1/x]` near 0, `Abs[Tan[x]]`
    near `Pi/2`, and `Abs[Log[x]]` near `0+` -- shapes where the
    inner limit diverges and the structural rewrite path would
    produce an intermediate `-g` that downstream layers can't
    classify (`Limit[-Log[x], x -> 0+]` is itself a separate gap).

2.  **One-sided structural rewrite.** For `Direction -> "FromAbove"`
    or `"FromBelow"`, walk `f` and replace each x-dependent
    `Abs[g]` subterm by `g` or `-g` according to the sign of `g`
    approaching the limit point in that direction:

    -   If `g(a)` evaluates to a non-zero literal (or to `±Infinity`),
        the sign is immediate.
    -   If `g(a) = 0`, the leading-order derivative test gives the
        sign: if `g^(k)(a)` is the first non-zero derivative,
        `g(a + t) = g^(k)(a)/k! · t^k + O(t^(k+1))`, so

            FromAbove  (t > 0):  sign(g) = sign(g^(k)(a))
            FromBelow  (t < 0):  sign(g) = sign(g^(k)(a)) · (-1)^k

    -   If `g(a)` is itself a `ComplexInfinity` / `DirectedInfinity`
        / unsigned divergent, recurse with `compute_limit(g, ...)`
        on the same direction and use the sign of that result.

    The rewritten expression is then handed back to `compute_limit`,
    which drops through the rest of the pipeline without any `Abs`
    in the way.

3.  **Two-sided splitting.** When `Direction` is `TwoSided` (or
    `Reals`), the sign of `g` may differ on either side of the
    kink, so the structural rewrite cannot pick a single
    replacement. The layer probes both one-sided limits via
    `compute_limit` and either returns the common value (when they
    agree) or `Indeterminate` (when they disagree). Two-sided
    splitting is only triggered when `f` actually contains an
    x-dependent `Abs` -- the existing
    `layer_onesided_disagree` is unchanged for the rest of the
    surface.

Concrete examples now resolving correctly:

```
Limit[Abs[x]/x, x -> 0, Direction -> -1]              ==> 1
Limit[Abs[x]/x, x -> 0, Direction -> 1]               ==> -1
Limit[Abs[x]/x, x -> 0]                               ==> Indeterminate
Limit[Abs[x^2 - 1]/(x - 1), x -> 1, Direction -> -1]  ==> 2
Limit[Abs[Sin[x] - x]/x^3, x -> 0, Direction -> -1]   ==> 1/6
Limit[Abs[Cos[x] - 1]/x^2, x -> 0]                    ==> 1/2
Limit[Abs[E^(2x) - 1]/x, x -> 0, Direction -> -1]     ==> 2
Limit[Abs[1/x], x -> 0]                               ==> Infinity
Limit[Abs[Log[x]], x -> 0, Direction -> -1]           ==> Infinity
Limit[Abs[Tan[x]], x -> Pi/2]                         ==> Infinity
Limit[Abs[x]/x, x -> Infinity]                        ==> 1
Limit[Abs[x]/x, x -> -Infinity]                       ==> -1
```

Coverage lives in `test_abs_directional` in
`tests/test_limit.c`.

## Conclusion

The layered structure in `src/limit.c` is not a natural decomposition
of a single theory; it is a decomposition of the *failure modes* of
the two tools we actually have -- substitution and series expansion.
Every layer in between exists because something in the test suite
showed us a shape where substitution collapses wrong and series
cannot see through the kernel. The pre-pass for reciprocal trig, the
log-merge at infinity, the bounded envelope, the radical canonicalizer,
the nine work packages of batches 4 through 12 -- each one is a patch
over an observed failure, not a general theory.

The near-term path forward is still the Gruntz algorithm for
exp-log limits at infinity. WP-5, WP-6, and WP-7 laid the three
pieces of machinery it needs -- `growth_exponent_upper` as a cheap
MRV proxy, hyperbolic exponentialisation, and a `Log[sum]` rewrite
-- so the remaining step is a proper most-rapidly-varying set
construction rather than a fresh architectural layer. After that,
an assumption database would unlock the symbolic-leading-coefficient
cases that currently bail out with `DirectedInfinity[c]`, and the
sign-from-outer-context iterated ArcTan limits that WP-1 cannot yet
resolve.

The tests pass, the memory is clean under Valgrind, and the ~130
cases in `tests/test_limit.c` are a regression line for every change
to the file.
