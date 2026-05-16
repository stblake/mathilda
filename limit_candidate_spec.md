# Limits in Mathilda

## 0. Scope and Goals

Mathilda needs a `Limit` built-in that matches Mathematica's user-facing
semantics while borrowing the proven algorithmic pipeline from MACSYMA's
`DELIMITER` and Maxima's modern `limit.lisp`. The implementation must be
a native C module consistent with Mathilda's expression model, memory
ownership rules, and evaluation pipeline (see `SPEC.md`).

### 0.1 User-facing Interface

The `Limit` built-in must accept three calling forms:

```mathematica
Limit[f, x -> x*]                               (* single variable *)
Limit[f, {x1 -> x1*, ..., xn -> xn*}]           (* iterated / nested *)
Limit[f, {x1, ..., xn} -> {x1*, ..., xn*}]      (* multivariate *)
```

Plus an optional `Direction` option:

| Setting | Meaning |
| :--- | :--- |
| `Reals` / `"TwoSided"` | From both real directions (default on the real line) |
| `"FromAbove"` / `-1` | Approach from larger values ($x \to x^{*+}$) |
| `"FromBelow"` / `+1` | Approach from smaller values ($x \to x^{*-}$) |
| `Complexes` | From all complex directions |
| `Exp[I θ]` | Along the ray in direction θ |
| `{dir1, ..., dirn}` | Per-variable directions (multivariate form) |

**Note on sign convention.** Mathematica's sign convention is the
*opposite* of MACSYMA's. In Mathematica, `Direction -> -1` means
*FromAbove* (approach from $+$, i.e. $0^+$), and `Direction -> +1`
means *FromBelow* (approach from $-$, i.e. $0^-$). The dispatcher
must translate `Direction` into the internal `ZERO_PLUS`/`ZERO_MINUS`
tags **after** this sign flip. This is a common source of bugs when
porting from Maxima-style code.

### 0.2 Return Values

`Limit` must be able to return any of the following:

| Return | Condition |
| :--- | :--- |
| A finite value or $\pm\infty$ | Limit exists |
| `Indeterminate` | Provably no limit (e.g. `Limit[Sin[1/x], x->0]` when asked for a single value; `Limit[Log[x]/x, x->0]`) |
| `Interval[{a, b}]` | Bounded but no single limit; returned only when bounds can be established (smallest bounds not guaranteed) |
| The original `Limit[...]` expression | No limit could be determined — remain unevaluated |

---

## 1. Core Data Structures and Global State

Before evaluating limits, the CAS must represent the specialized
mathematical concepts of infinity, boundedness, and directionality.

* **Directional Limit Nodes.** Extend the internal AST vocabulary to
  handle directional limits inherently. These are *internal* tokens used
  during limit computation, not necessarily user-visible symbols:
    * `ZERO_PLUS` ($0^+$) and `ZERO_MINUS` ($0^-$).
    * `INF` (positive real infinity) and `MINF` (negative real infinity).
    * `INFINITY` (complex/directionless infinity).
    * `UND` (undefined due to discontinuity — maps to user-facing
      `Indeterminate`).
    * `IND` (indefinite but bounded, e.g. $\lim_{x\to\infty} \sin x$ —
      maps to user-facing `Interval[...]` when bounds are known, else to
      `Indeterminate` or unevaluated).
* **The Complex Flag.** Maintain a context-passed state boolean `CPLX`.
  If the limit point is `INFINITY` or the expression involves `%I`
  ($\sqrt{-1}$), set `CPLX = TRUE`. When true, one-sided approaches are
  invalid and `INF`/`MINF` fold into `INFINITY`.
* **Limit Context Structure.** A `LimitCtx` struct passed through the
  pipeline carrying: the variable, the target point, the direction ray
  (an `Expr*` of the form `Exp[I θ]` or one of the directional tokens),
  the `CPLX` flag, and a recursion-depth counter for stack safety.

### 1.1 Existing Infrastructure to Leverage

The following Mathilda facilities already exist and should be used by
the limit evaluator rather than reimplemented:

* `D` / `Derivative` (native C, per recent commit `8081f55`) — L'Hospital.
* `Series` (per commits `f416915`, `5c50026`) — asymptotic expansion.
* `Simplify`, `Together`, `Cancel`, `Apart` — pre-normalization.
* `Expand`, `Collect` — polynomial manipulation.
* `PolynomialQ`, `Numerator`, `Denominator`, `Coefficient` — rational
  function dispatch.

---

## 2. The Evaluation Pipeline

The main entry point is a C function mirroring the signature

```c
Expr* builtin_limit(Expr* res);
```

which conforms to Mathilda's built-in contract (ownership of `res`, `NULL`
to leave unevaluated, else a freshly allocated result). Internally it
routes through the six layers below. Each layer either returns a result
(short-circuiting the rest) or passes control down.

### Layer 0 — Interface Normalization

Responsible for translating the three Mathematica calling forms and the
`Direction` option into a canonical internal representation before any
mathematics is attempted.

1. **Parse the second argument:**
    * `x -> a` → single-variable limit.
    * `{x1 -> a1, ..., xn -> an}` → iterated limit: fold from
      *innermost* (last in the list) outward. Note that iterated limit
      order is semantically significant; `Limit[x/y, {x->0, y->0}]` ≠
      `Limit[x/y, {y->0, x->0}]`.
    * `{x1, ..., xn} -> {a1, ..., an}` → multivariate limit (see
      Layer 5a).
2. **Parse the `Direction` option.** Normalize to a direction ray
   expressed as an `Expr*`:
    * `Reals`, `"TwoSided"`, or omitted on the real line → two-sided.
    * `"FromAbove"` or `-1` → `ZERO_PLUS` internal tag.
    * `"FromBelow"` or `+1` → `ZERO_MINUS` internal tag.
    * `Complexes` → set `CPLX = TRUE`.
    * `Exp[I θ]` → complex ray; the limit reduces to a single-variable
      limit after the substitution $x = x^* + t\,e^{i\theta}$, $t \to 0^+$.
      This is the general framework that subsumes all directional cases.
    * `{dir1, ..., dirn}` → zip with variable list for multivariate.
3. **Produce a `LimitCtx`.** Pass it, along with the (possibly
   substituted) expression, to Layer 1.

### Layer 1 — Fast Paths

Cheap structural checks that resolve the majority of trivial cases.

* If the expression is independent of the limit variable, return it
  unchanged.
* If the expression *is* the variable, return the limit point.
* **Continuous substitution:** if the expression is a sum/product/power
  of known-continuous functions and the substituted value is finite and
  well-defined (no division by zero, no log of zero, etc.), return the
  substituted value. This handles e.g. `Limit[Sin[x], x->Pi/2] = 1`
  and `Limit[a x^2 + 4x + 7, x->1] = 11 + a` in one step.

### Layer 2 — Series-Based Evaluation *(promoted to first-class)*

Most of the classical DELIMITER test cases — `(Cos[x]-1)/(Exp[x^2]-1)`,
`(1+A x)^(1/x)`, `(Tan[x]-x)/x^3`, etc. — are solved most robustly by
Taylor / Laurent / Puiseux series. Mathilda already has `Series`, so we
elevate series expansion above L'Hospital in the dispatch order. This
avoids the L'Hospital stack-blow guardrails entirely for most cases.

1. If the limit point is `INF`/`MINF`, substitute $x = 1/t$ and reduce
   to $t \to 0^+$ or $t \to 0^-$.
2. Call `Series[expr, {x, point, k}]` starting with `k = 4` and
   doubling up to a configurable cap (e.g. `k = 32`) until a nonzero
   leading term is found or series cannot be computed.
3. Read off the limit from the leading term:
    * Leading exponent $> 0$ → limit is $0$.
    * Leading exponent $= 0$ → limit is the leading coefficient.
    * Leading exponent $< 0$ → limit is $\pm\infty$ (sign from
      coefficient; use direction if odd exponent).
4. If `Series` fails or returns a non-convergent form, fall through.

### Layer 3 — Rational Functions

If the (pre-normalized via `Together`) expression matches $P(x)/Q(x)$
with both polynomial in $x$:

* **Finite limit point $L$:** if $Q(L) \neq 0$, return $P(L)/Q(L)$. If
  $P(L) = Q(L) = 0$, factor out $(x - L)$ by polynomial division and
  recurse. If only $Q(L) = 0$, the limit is $\pm\infty$ or `INFINITY`
  depending on the sign of $P(L)$ and the multiplicity / direction.
* **Infinite limit point:** extract leading coefficients and degrees.
  Return $0$, the ratio of leading coefficients, or $\pm\infty$ based on
  strict degree comparison.

### Layer 4 — RP-Forms (Radicals of Polynomials)

RP-forms are expressions composed of polynomials and rational powers of
polynomials.

* If the limit point is `MINF`, substitute $x = -y$ and evaluate at
  `INF`.
* **Asymptotic dominance:** extract the highest exponent of $x$ in
  numerator and denominator via an `expo()` helper.
* Replace polynomials under radicals with their leading terms. If the
  dominant exponents don't match, expand radicals in a Laurent series
  around infinity (reusing Layer 2) to produce a rational function and
  recurse into Layer 3.

### Layer 5 — General Heuristics

For expressions that fail the previous layers, extract $N(x)$ and $D(x)$
via `Numerator`/`Denominator` after `Together`.

1. **L'Hospital's Rule.** Apply to indeterminate forms $0/0$ or
   $\infty/\infty$ using the existing `D` built-in.
    * *Critical guardrail.* Uncontrolled recursion will blow the C
      stack. Count the distinct nonrational components in the AST via a
      `LeafCount`-style helper. If this complexity metric *strictly
      grows* for three consecutive L'Hospital derivations, abort and
      fall through. Also abort if the `LimitCtx` recursion depth exceeds
      a configurable cap.
2. **Order of Infinity (STRENGTH routine).** Used when limits evaluate
   to $\infty - \infty$ or when L'Hospital stalls.
    * Assign hierarchical strengths: $\text{EXP} \gg \text{POLY} \gg
      \text{LOG} \gg \text{CONST}$.
    * Traverse `Times` products and take the maximum strength. If the
      numerator's strength strictly dominates, the limit is
      $\pm\infty$; if the denominator's, the limit is $0$.
3. **Logarithmic Reduction.** Transform exponential indeterminate forms
   ($1^\infty$, $0^0$, $\infty^0$) via the identity
   $f(x)^{g(x)} = e^{g(x)\,\log f(x)}$ and recurse.

### Layer 5a — Multivariate Limits

Multivariate limits are fundamentally harder than single-variable limits
because path-independence must be verified, not merely assumed. The
initial implementation should cover a modest but useful subset and
cleanly return unevaluated for the rest:

1. **Jointly continuous fast path.** If the substitution of all
   variables simultaneously gives a finite, well-defined value and the
   expression is a sum/product/composition of jointly continuous
   functions on a neighborhood of the target, return that value.
   Handles e.g. `Limit[Tan[x y]/(x y), {x,y}->{0,0}] = 1` after one
   `Series` in the product $u = xy$.
2. **Polar substitution (2D only).** For targets of the form
   `{x, y} -> {0, 0}`, try $x = r\cos\theta$, $y = r\sin\theta$ and
   take the univariate limit $r \to 0^+$. If the result is independent
   of $\theta$, that is the limit; if it depends on $\theta$, the
   multivariate limit does not exist → return `Indeterminate`.
3. **Path-dependence check (heuristic).** Try a small set of linear
   paths $y = m x$ for several $m$; if the univariate limits disagree,
   return `Indeterminate`.
4. **Fallback.** Return the expression unevaluated. Document clearly in
   the docstring that the multivariate form is a best-effort
   implementation.

### Layer 6 — Bound Analysis for `Interval` Returns

If none of the layers above produces a value, attempt to establish
finite bounds so that an `Interval[{lo, hi}]` can be returned:

1. Decompose the expression into bounded and unbounded factors.
2. If every unbounded factor has a finite limit and each bounded factor
   has known range bounds (e.g. `Sin`, `Cos` ∈ `[-1, 1]`; `ArcTan` ∈
   `(-Pi/2, Pi/2)`), propagate the bounds to the whole expression.
3. Return `Interval[{lo, hi}]` if bounds are finite; otherwise
   `Indeterminate` if the expression is known-bounded but bounds cannot
   be computed; otherwise leave unevaluated.

This layer justifies e.g. `Limit[Sin[1/x], x->0] = Interval[{-1, 1}]`.

---

## 3. C Implementation Notes

* **Module location:** `src/limit.c` / `src/limit.h`, initialized from
  `core.c` via a new `limit_init()` call, following the standard
  Mathilda module pattern.
* **Attributes:** `Limit` should be `ATTR_PROTECTED`. It should
  *probably not* be `ATTR_HOLDALL` — the first argument (the function
  expression) should be evaluated normally — but the rule `x -> a` must
  not be prematurely evaluated if `x` has an OwnValue in scope. Consider
  `ATTR_HOLDREST` on the second argument, with explicit evaluation
  inside the built-in after substitution handling.
* **Memory ownership:** the built-in takes ownership of `res`. Every
  intermediate `Series` / `D` / `Together` call returns a fresh `Expr*`
  that must be freed. Use the standard NULL-out-before-free pattern
  (SPEC §4.2) for stealing arguments.
* **Recursion safety:** the `LimitCtx` must carry a depth counter;
  abort into "unevaluated" rather than segfaulting on pathological
  input.
* **Direction ray representation:** carry the direction as an `Expr*`
  inside `LimitCtx` (either a canonical tag symbol or an `Exp[I θ]`
  expression). This lets Layer 0's `Exp[I θ]` substitution path share
  code with the classical directional cases.

---

## 4. Test Suite

The test suite must cover all four return-value categories, all three
calling forms, all `Direction` settings, and symbolic limit points.
Add the tests to `tests/test_limit.c` using the existing
`assert_eval_eq` helper.

### 4.1 Definite Numeric Limits (DELIMITER Reference)

| Input | Expected | Notes |
| :--- | :--- | :--- |
| `Limit[x^Log[1/x], x->Infinity]` | `0` | Logarithmic reduction |
| `Limit[(Cos[x]-1)/(Exp[x^2]-1), x->0]` | `-1/2` | Series |
| `Limit[(1+A x)^(1/x), x->0]` | `E^A` | Log reduction + series |
| `Limit[x^2 Sqrt[4 x^4 + 5] - 2 x^4, x->Infinity]` | `5/4` | RP-form |
| `Limit[1/x - 1/Sin[x], x->0, Direction->-1]` | `0` | FromAbove |
| `Limit[E^(x Sqrt[x^2+1]) - E^(x^2), x->Infinity]` | `Infinity` | Strength |
| `Limit[E^x Log[x] / (E^(x^3+1) Sqrt[Log[x^4+x+1]]), x->Infinity]` | `0` | Strength |
| `Limit[1/(x^3 - 6 x^2 + 11 x - 6), x->2, Direction->+1]` | `Infinity` | FromBelow; note Mathematica sign |
| `Limit[(x Sqrt[x+5] + 1)/(Sqrt[4 x^3 + 1] + x), x->Infinity]` | `1/2` | RP-form |
| `Limit[Tan[x]/Log[Cos[x]], x->Pi/2, Direction->+1]` | `-Infinity` | FromBelow |
| `Limit[D[(z - I Pi/2) z (z - 2 Pi I)/(Sinh[z] - I), z], z -> I Pi/2]` | `-2 I Pi` | Complex |

### 4.2 Mathematica Interface Examples

| Input | Expected |
| :--- | :--- |
| `Limit[Sin[x]/x, x->0]` | `1` |
| `Limit[(1 + Sinh[x])/Exp[x], x->Infinity]` | `1/2` |
| `Limit[x/(x - 2 a), x->a]` | `-1` |
| `Limit[3 x^2/(x^2 - 2), x->-Infinity]` | `3` |
| `Limit[Sqrt[x - 1]/x, x->1]` | `0` |
| `Limit[Sin[2 x]/Sin[x], x->Pi]` | `-2` |
| `Limit[x Log[x^2], x->0]` | `0` |
| `Limit[(x^a - 1)/a, a->0]` | `Log[x]` |
| `Limit[Sqrt[x^2 + a x] - Sqrt[x^2 + b x], x->Infinity]` | `(a - b)/2` |
| `Limit[ArcTan[x^3 - x], x->-Infinity]` | `-Pi/2` |
| `Limit[(1 + a/x)^(b x), x->Infinity]` | `E^(a b)` |

### 4.3 Directional Limits (Sign Convention Tests)

| Input | Expected |
| :--- | :--- |
| `Limit[1/x, x->0, Direction->"FromBelow"]` | `-Infinity` |
| `Limit[1/x, x->0, Direction->1]` | `-Infinity` |
| `Limit[1/x, x->0, Direction->"FromAbove"]` | `Infinity` |
| `Limit[1/x, x->0, Direction->-1]` | `Infinity` |
| `Limit[Tan[x], x->Pi/2, Direction->"FromBelow"]` | `Infinity` |
| `Limit[Tan[x], x->Pi/2, Direction->"FromAbove"]` | `-Infinity` |

### 4.4 `Indeterminate` Returns

| Input | Expected |
| :--- | :--- |
| `Limit[Sin[1/x], x->0]` | `Interval[{-1, 1}]` (Layer 6) |
| `Limit[Log[x]/x, x->0]` | `Indeterminate` |
| `Limit[Sin[Tan[x]], x->Pi/2, Direction->Reals]` | `Indeterminate` |
| `Limit[Sin[x], x->Infinity]` | `Indeterminate` |

### 4.5 `Interval` Returns (Bounded, No Limit)

| Input | Expected |
| :--- | :--- |
| `Limit[Sin[1/x], x->0]` | `Interval[{-1, 1}]` |

### 4.6 Unevaluated Returns

| Input | Expected |
| :--- | :--- |
| `Limit[x f[x], x->0]` | `Limit[x f[x], x->0]` (unknown `f`) |

### 4.7 Multivariate Examples

| Input | Expected |
| :--- | :--- |
| `Limit[Tan[x y]/(x y), {x, y} -> {0, 0}]` | `1` |
| `Limit[x y/(x^2 + y^2), {x, y} -> {0, 0}]` | `Indeterminate` (path-dependent) |

### 4.8 Iterated (Nested) Examples

| Input | Expected |
| :--- | :--- |
| `Limit[x/(x + y), {x -> 0, y -> 0}]` | `0` (inner `y->0` first gives `1`; wait — see below) |

*Test-writer note:* iterated limits apply innermost (rightmost) first.
Cross-check each expected value against Mathematica before committing.

---

## 5. Incremental Delivery Plan

A suggested ordering so that a useful `Limit` ships early and grows:

1. **Milestone 1 — Layers 0, 1, 3.** Interface normalization, fast
   paths, and rational functions. Covers all polynomial/rational
   examples.
2. **Milestone 2 — Layer 2.** Series-based evaluation. Unlocks most
   transcendental cases cheaply.
3. **Milestone 3 — Layer 5.1 + 5.3.** L'Hospital and logarithmic
   reduction, with the complexity-growth guardrail.
4. **Milestone 4 — Layers 4, 5.2.** RP-forms and STRENGTH. Handles the
   hardest asymptotic cases.
5. **Milestone 5 — Layer 6.** Bound analysis / `Interval` returns.
6. **Milestone 6 — Layer 5a.** Multivariate and iterated forms.

Each milestone ships with passing tests from the corresponding
sections of §4 and updates to `Mathilda_spec.md`.
