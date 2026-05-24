# PossibleZeroQ — Hybrid Symbolic/Numeric Zero‑Recognition Plan

Author: design proposal, 2026-05-24
Status: planning — not yet implemented

---

## 1. Goal

Implement the Mathematica built-in `PossibleZeroQ[expr]`:

> *Gives `True` if basic symbolic and numerical methods suggest that `expr` has
> value zero, and `False` otherwise. The general problem is undecidable
> (Richardson 1968); `PossibleZeroQ` is a quick but not always accurate test.*

Attributes: `Listable`, `Protected`.

Required semantics from the reference examples:

| Input | Result | Reason |
|-------|--------|--------|
| `E^(I Pi/4) - (-1)^(1/4)` | `True` | Closed numeric expression that simplifies to 0 |
| `(x+1)(x-1) - x^2 + 1` | `True` | Polynomial identity over `Z[x]` |
| `(E+Pi)^2 - E^2 - Pi^2 - 2 E Pi` | `True` | Numeric, exact cancellation |
| `E^Pi - Pi^E` | `False` | Numeric, ~0.6814 |
| `2^(2 I) - 2^(-2 I) - 2 I Sin[Log[4]]` | `True` | Closed numeric complex identity |
| `1/x + 1/y - (x+y)/(x y)` | `True` | Rational identity after `Together`/`Cancel` |
| `Sqrt[x^2] - x` | `False` | Not identically zero (false for x<0) |

Out of scope (Phase 2 / future): expressions involving non-implemented heads
(e.g. `Erf`, hypergeometrics). The reference example using `Erf` triggers
Mathematica's `ztest1::` warning even there — see §10.

---

## 2. Theoretical background and literature

Zero-recognition is the central oracle assumed by symbolic integration
(Risch 1969, Bronstein 2005), simplification, polynomial-GCD and
differential-equation solving. Without it, every algorithm above is forced
into heuristics.

### 2.1 Decidability landscape

- **Richardson (1968)**, *Some unsolvable problems involving elementary
  functions of a real variable*, JSL 33(4). Proves identity-testing is
  undecidable for the class of expressions built from `Q`, `pi`, `ln 2`, `x`,
  `+`, `−`, `×`, `sin`, `|·|` and function composition. So `PossibleZeroQ`
  is necessarily heuristic.
- **Caviness (1970)** strengthens Richardson by removing `|·|` if `e^x` is
  added.
- **Wang (1974)**, *The undecidability of the existence of zeros of real
  elementary functions*, JACM 21(4) — analogous result for zero-existence
  over an interval.

These results say: *no algorithm* can decide identity to zero on the full
elementary class. Hence the Mathematica name choice: `PossibleZeroQ`.

### 2.2 Decidable sub-cases we should fully solve

- **Polynomial identity testing over Q (or any computable field)** is
  decidable in time linear in the dense representation after `Expand`. Reduces
  to "all coefficients of the canonical form are 0". Already implemented in
  Mathilda as `is_zero_poly` at `src/poly/poly.c:1066`.
- **Rational-function identity testing over Q(x_1,…,x_n)**: decidable via
  `Together`/`Cancel` → reduce to polynomial case. Already wired in
  `intrat_canonic` at `src/calculus/intrat.c:106` (`Cancel[Together[e]]`).
- **Algebraic-number identity testing**: decidable; reduces to working in a
  primitive element representation. Mathilda has partial support
  (`src/poly/algfac.c`); leveraging it is a Phase 2 item.

### 2.3 Heuristic numeric methods we will adopt

- **Schwartz (1980)**, *Fast probabilistic algorithms for verification of
  polynomial identities*, JACM 27(4) and **Zippel (1979)**, *Probabilistic
  algorithms for sparse polynomials*, EUROSAM. The Schwartz–Zippel lemma:
  *a non-zero polynomial of total degree d over a field F evaluates to 0 at
  a uniformly random point in S^n ⊆ F^n with probability ≤ d/|S|.* So
  evaluating at k random points in a set of size ≥ Kd gives a one-sided
  Monte-Carlo identity test with error ≤ (d/Kd)^k = K^{-k}. For typical
  degrees ≤ 100 and |S| ≥ 10⁹, four samples give false-positive ≤ 10⁻³⁰.
- **DeMillo & Lipton (1978)**, *A probabilistic remark on algebraic program
  testing*, IPL 7(4). Earlier statement of the same lemma.
- **Stoutemyer (1989)**, *Crimes and misdemeanors in the computer algebra
  trade*, Notices AMS — discusses how transcendental zero-recognition is
  routinely heuristic in CASes, and warns about catastrophic cancellation
  when using floating point.
- **Bronstein (2005)**, *Symbolic Integration I — Transcendental
  Functions*, 2nd ed., Springer, §3.1. Describes the "structure theorem"
  approach for elementary towers (Risch, Trager): test transcendence of
  monomials, then test 0-equivalence in `K(t_1,…,t_n)`. The fall-back to
  numeric evaluation at multiple random points is explicitly endorsed.
- **Bailey, Borwein & Plouffe (1997)** and **PSLQ (Ferguson, Bailey 1992)** —
  integer-relation algorithms that detect "is this real number actually a
  rational combination of known constants?" Used internally by Mathematica's
  `PossibleZeroQ` for hard cases. *Out of scope for v1 but mentioned for
  Phase 3.*
- **Aslaksen (1996)**, *Can your computer do complex analysis?*, Math.
  Intelligencer 18(3) — survey of CAS failures around `Sqrt`, `Log` branch
  cuts. Justifies our `Sqrt[x^2] − x ≠ 0` answer.
- **Hubert (1997)** and **Aubry, Lazard & Moreno-Maza (1999)** — equality
  testing in differential and difference algebra via regular chains. Out of
  scope; mentioned only to keep `PossibleZeroQ` extensible.

### 2.4 Practical strategy this plan adopts

The proven, widely deployed recipe (Maple `testzero`, Mathematica
`PossibleZeroQ`, SymPy `is_zero`, FriCAS `zero?`):

1. Cheap structural checks (literal `0`, `Complex[0,0]`, …) → constant time.
2. Rational normalization (`Together` ∘ `Cancel` ∘ `Expand`) →
   `is_zero_poly`. Decides every rational identity over `Q`.
3. **Numeric path** when expression has no free symbols: evaluate at machine
   precision, and if the result is "near zero" bump to high-precision MPFR
   and confirm. Cancellation tolerance is set proportional to the *operand
   magnitudes*, not absolute, to handle expressions like
   `(E+Pi)^2 - E^2 - Pi^2 - 2 E Pi`.
4. **Schwartz–Zippel** path when free symbols are present and the expression
   is not a rational function: substitute random rational values, recurse
   into the numeric path, repeat with independent samples.

Probabilistic guarantee (one-sided): false-positives ≤ 10⁻²⁰ under standard
sampling parameters. False-negatives (returning `False` for something
actually zero) only happen when the symbolic engine fails *and* the numeric
sample falls on a "branch-cut surprise" — handled by §6.5.

---

## 3. Architecture

A new module `src/zero_test.c` / `src/zero_test.h`. Single public entry:

```c
/* zero_test.h */
#ifndef MATHILDA_ZERO_TEST_H
#define MATHILDA_ZERO_TEST_H

#include "expr.h"

/* Public C API — usable from other builtins (Equal, Simplify, integrators). */
typedef enum {
    ZERO_TEST_FALSE = 0,   /* proved (or strongly believed) non-zero          */
    ZERO_TEST_TRUE  = 1,   /* proved zero (rational) or strongly believed zero */
    ZERO_TEST_UNKNOWN = 2  /* could not decide (heuristic exhausted)          */
} ZeroTestResult;

ZeroTestResult zero_test_decide(const Expr* e);

/* The Mathilda builtin wrapper — registered via symtab_add_builtin. */
Expr* builtin_possible_zero_q(Expr* res);

void zero_test_init(void);
#endif
```

`PossibleZeroQ[expr]` collapses `UNKNOWN` to `True` (matching Mathematica's
`PossibleZeroQ::ztest1` behaviour — when unsure, *assume zero* and emit a
message). This is the documented behaviour of the reference example with
`Erf[...]`. We emit the message via the standard Mathilda message channel
(see §9). The richer three-valued API is exported for internal callers
(notably `Equal`, `Simplify`, integration) that need to distinguish.

### 3.1 Pipeline overview

```
                 ┌──────────────────────────────────────┐
PossibleZeroQ ─▶ │ Stage 0: trivially zero / non-zero?  │ ─▶ True / False
                 └──────────────────────────────────────┘
                                 │ unknown
                                 ▼
                 ┌──────────────────────────────────────┐
                 │ Stage 1: rational normalization      │ ─▶ True if is_zero_poly
                 │   Together → Cancel → Expand         │
                 │   → is_zero_poly                     │
                 └──────────────────────────────────────┘
                                 │ unknown
                                 ▼
                 ┌──────────────────────────────────────┐
        no free  │ Stage 2: numeric closed-form         │ ─▶ True / False
        symbols  │   numericalize @ 53 bits             │
       ────────▶ │   if |z| < ε · scale, retry @ 200    │
                 │   bits, then 500 bits; final verdict │
                 └──────────────────────────────────────┘
                                 │
                  has free symbols
                                 ▼
                 ┌──────────────────────────────────────┐
                 │ Stage 3: Schwartz–Zippel             │ ─▶ True / UNKNOWN→True
                 │   k random rational substitutions    │     False on any miss
                 │   each ⇒ recurse into Stage 2        │
                 └──────────────────────────────────────┘
```

### 3.2 Why this order

Stages are ordered by **average cost**, *not* by power. Stage 1 is O(n log n)
in expression size and *decides* polynomial identities exactly — never
descend to floating point if symbolic cancellation suffices. Stage 2 is
~µs at machine precision; only the rare "near-zero" case pays the MPFR
re-evaluation cost. Stage 3 is invoked only when free symbols make
Stages 0/1 unable to decide.

---

## 4. Stage details

### 4.1 Stage 0 — structural shortcuts

Decide in constant or shallow-linear time. No allocation.

| Condition | Verdict |
|-----------|---------|
| `EXPR_INTEGER` value 0 | True |
| `EXPR_REAL` value `0.0` or `-0.0` | True |
| `EXPR_BIGINT` with `mpz_sgn(z) == 0` | True |
| `EXPR_MPFR` with `mpfr_zero_p(x)` | True |
| `Complex[0, 0]` (after evaluation) | True |
| `Times[…, 0, …]` head | True (defensive — should already be folded) |
| `EXPR_STRING`, `Hold[…]`, `Symbol` that is not numeric constant | False |
| Non-numeric symbol (e.g. bare `x`) with no OwnValue | False |
| Numeric-constant symbol (`Pi`, `E`, `EulerGamma`, …) | False |
| `EXPR_INTEGER`, `EXPR_REAL`, `EXPR_BIGINT`, `EXPR_MPFR` non-zero | False |

Cost: O(1) plus one symbol-table lookup for symbol cases. No allocations.

### 4.2 Stage 1 — rational normalization

```c
Expr* canon = call_internal(Cancel ∘ Together)(e);
if (is_zero_poly(canon)) { expr_free(canon); return TRUE; }
/* Try also pure Expand for cases without denominators. */
Expr* expanded = expr_expand(e);                     /* src/expand.c:172 */
bool z = is_zero_poly(expanded);
expr_free(expanded);
if (z) { expr_free(canon); return TRUE; }
```

`is_zero_poly` (`src/poly/poly.c:1066`) already does the right thing for
multivariate polynomials over `Z`: it expands, collects variables, extracts
coefficients via `internal_coefficientlist`, and recursively checks every
coefficient. Depth bound 32 prevents pathological recursion; if exceeded,
the verdict is "not proved zero" (we treat that as continue, not False).

**`Together`/`Cancel` already exist** in `src/rat.c:1438,1450` with memoization
via `FactorMemo`, so this is a near-free reuse.

**Caveat — `Sqrt[x^2]`**. `Expand` will not turn `Sqrt[x^2] - x` into `0`,
and that is correct: the identity is false on the principal branch when
`x < 0`. So Stage 1 returns "not zero" here, which is what we want.

Allocations are released with `expr_free` before each early return — the
NULL-out-before-free pattern is not needed (we never reuse fragments).

### 4.3 Stage 2 — numeric, adaptive precision

Used when **`expr_has_free_symbols(e) == false`**. Free symbol = `EXPR_SYMBOL`
with no `OwnValue` *and* not a registered numeric constant
(`Pi`, `E`, `EulerGamma`, `Catalan`, `GoldenRatio`, `Degree`,
`I`, `Infinity`, `True`, `False`). Implemented via `collect_symbols_in`
(`src/sort.c:208`) plus a small allow-list check.

Algorithm:

```c
static const int PRECISION_LADDER[] = { 53, 200, 500, 1000 };   /* bits */
static const double TOLERANCE_FACTOR = 8.0;                     /* slack */

Expr* magnitude_estimate(const Expr* e);  /* sum of |leaf| over leaves */

for (int i = 0; i < 4; ++i) {
    NumericSpec spec = numeric_spec_from_bits(PRECISION_LADDER[i]);
    Expr* z = numericalize(e, spec);
    if (!z) return UNKNOWN;                /* numericalize failed */
    /* Normalize Complex[re, im] → real magnitude. */
    double mag = numeric_abs_double(z);    /* may use mpfr_get_d */
    double scale = numeric_magnitude_scale(e, spec);  /* see below */
    expr_free(z);

    double tol = scale * TOLERANCE_FACTOR * pow(2.0, -PRECISION_LADDER[i] + 6);
    if (mag > tol) return FALSE;           /* definitively non-zero */
    if (mag == 0.0)            return TRUE;
    /* mag is small but nonzero — possibly catastrophic cancellation; retry */
}
return TRUE;          /* survived all 4 ladder rungs ⇒ assume zero */
```

**Scale**: `numeric_magnitude_scale(e, spec)` computes `Σ |numericalize(leaf_i)|`
over numeric leaves of `e` (replacing free constants with their numeric
value, walking only the numeric skeleton). This is the same quantity used
by IEEE-style cancellation analysis — if your operands have magnitude
`~10`, a residual of `10·2⁻⁴⁷` is plausible round-off, but `10·2⁻³⁰` is
not. Without this scale factor the absolute test `|z| < 10⁻¹⁰` is wrong
both ways: it would falsely declare `1e-15` zero, and falsely declare
`1e-100 - 1e-100` non-zero after MPFR refinement.

**Why a ladder?** Catastrophic cancellation can hide an exact zero behind
spurious rounding noise at 53 bits. Bumping to 200 bits and re-evaluating
distinguishes two regimes:
- *Actually zero.* Residual shrinks geometrically with precision (it is
  bounded by `2^{-p} · scale`). Each rung divides `mag` by ~`2^{147}`,
  so a true zero certainly survives every rung.
- *Cancellation that looks like zero at low precision.* The residual
  stabilises near its true non-zero value once precision exceeds the
  cancellation depth. A `False` verdict at rung 2 or 3 is reliable.

Citing **Higham (2002)**, *Accuracy and Stability of Numerical Algorithms*,
2nd ed., SIAM, §1.7 — the standard treatment of how doubling working
precision reveals true vs. apparent zeros.

**Complex contagion.** `numericalize_function` (`src/numeric.c:479`)
already returns `Complex[re, im]` for complex-valued numeric expressions.
`numeric_abs_double` must accept that shape: `|re + i·im| = hypot(re, im)`,
using `mpfr_hypot` in MPFR mode.

**`I` and `(-1)^(1/4)`.** Both go through `Power`/`Complex` normalization
in the evaluator before reaching us; our reference test
`E^(I Pi/4) - (-1)^(1/4)` will evaluate both halves to
`Complex[Sqrt[2]/2, Sqrt[2]/2]`, subtract, and *Stage 1* already returns
True via `Together`/`Cancel` on rational combinations of `Sqrt[2]`.
Stage 2 is the safety net.

### 4.4 Stage 3 — Schwartz–Zippel for symbolic non-rational expressions

Reached when free symbols are present and Stages 0/1 did not decide.
Examples: `Sqrt[Sin[x]^2 + Cos[x]^2] - 1`, `Log[E^x] - x` (latter is also
false in general on the complex plane — see §6.5).

```c
#define ZT_NUM_SAMPLES        4         /* independent trials                */
#define ZT_SAMPLE_RANGE       (1 << 20) /* uniform on [-2^20, 2^20] ∩ Q     */
#define ZT_SAMPLE_DENOM_MAX   (1 << 16) /* random denominator bound          */

bool any_failure = false;
for (int trial = 0; trial < ZT_NUM_SAMPLES; ++trial) {
    Expr* sub = substitute_symbols_with_random_rationals(e);
    ZeroTestResult r = numeric_decide(sub);  /* = Stage 2 */
    expr_free(sub);
    if (r == ZERO_TEST_FALSE)   { any_failure = true; break; }
    if (r == ZERO_TEST_UNKNOWN) { any_failure = true; break; }
    /* r == TRUE → continue, need k consecutive hits */
}
return any_failure ? FALSE : TRUE;
```

By Schwartz–Zippel, false-positive probability after `k=4` samples on a
set of size `2·2²⁰·2¹⁶ ≈ 2³⁷` is bounded by `(d/|S|)^k`. For `d ≤ 1000`
(very generous bound on combined polynomial degree), that's `2^{-108}`.
For arbitrary elementary expressions the lemma does not strictly apply,
but in practice the same exponential decay holds — see Bronstein 2005
§3.1 for the standard argument.

Sampling implementation reuses `random.c`:
- `RandomInteger[{-2^20, 2^20}]` for numerator (via `internal_randominteger`)
- `RandomInteger[{1, 2^16}]` for denominator
- Combine via `expr_new_function("Rational", …)` and re-evaluate

Branch-cut traps: pure imaginary sample points are *also drawn* with
probability 1/2, so we sample uniformly from `Q[i]` rather than `Q`. This
catches `Sqrt[x^2] - x = 0` (false at `x = -2`) the same way it catches
identities valid only on `R`.

### 4.5 Trip-wire / fall-through

If *all* stages decline:

- Stage 1 didn't reduce to a known polynomial form.
- Stage 2/3 returned `UNKNOWN` because `numericalize` produced an
  unevaluated symbol (e.g. `Erf[...]`).

We emit message `PossibleZeroQ::ztest1` ("Unable to decide whether
numeric quantity X is equal to zero. Assuming it is.") and return `True`.
This matches the Mathematica behaviour quoted in the prompt.

---

## 5. File changes

| File | Change |
|------|--------|
| `src/zero_test.h` | **new** — public header (API in §3) |
| `src/zero_test.c` | **new** — pipeline implementation |
| `src/sym_names.h` | add `extern const char* SYM_PossibleZeroQ;` |
| `src/sym_names.c` | intern the name; add entry to bootstrap list |
| `src/core.c` | call `zero_test_init()` from `core_init()` |
| `src/info.c` | add `symtab_set_docstring("PossibleZeroQ", "…")` |
| `src/attr.c` | (only if attribute set is module-managed) set Listable+Protected |
| `tests/test_zero_test.c` | **new** — unit tests (§7) |
| `tests/CMakeLists.txt` | register `zero_test_tests` executable |
| `docs/spec/builtins/expression-information.md` | add `PossibleZeroQ` entry |
| `docs/spec/changelog/2026-05.md` | summary under "Added" |
| `Mathilda_spec.md` | no change (overview file, no new top-level section) |

`zero_test_init()` body:

```c
void zero_test_init(void) {
    symtab_add_builtin("PossibleZeroQ", builtin_possible_zero_q);
    symtab_get_def("PossibleZeroQ")->attributes
        |= ATTR_LISTABLE | ATTR_PROTECTED;
}
```

Docstring (info.c) — terse per project rule, no examples:

> "PossibleZeroQ[expr] gives True if symbolic and numerical methods suggest
> that expr has value zero, and False otherwise. The general problem of
> deciding whether an expression is zero is undecidable; PossibleZeroQ is
> a quick but not always accurate test."

---

## 6. Correctness & robustness considerations

### 6.1 Memory management

All builtin allocations follow the §4 SPEC.md ownership contract:

- `builtin_possible_zero_q(res)` consumes `res` indirectly via the evaluator
  and returns either `True`, `False`, or `NULL`. Result expressions are
  freshly allocated.
- Internal helpers consume their input only when documented. `zero_test_decide`
  takes `const Expr*` — it never frees its input, only intermediates.
- `numericalize`, `expr_expand`, `internal_together`, `internal_cancel` are
  documented as returning fresh ownership; we pair every call with
  `expr_free` on the path out.
- Tests run under `valgrind --leak-check=full --error-exitcode=1` in CI.

### 6.2 Listable threading

`ATTR_LISTABLE` is handled by `apply_listable` in `src/eval.c:300` *before*
our builtin is invoked. We see only scalar calls, e.g.
`PossibleZeroQ[{0, x-x, Pi-Pi}]` arrives as three independent
`PossibleZeroQ[…]` calls. No threading code in the builtin itself.

### 6.3 Holding semantics

`PossibleZeroQ` does **not** have Hold attributes — the evaluator
fully evaluates the argument before we see it. So when the user writes
`PossibleZeroQ[1/x + 1/y - (x+y)/(x*y)]`, the evaluator first attempts to
simplify the addition; if the result still contains the algebraic form, we
get whatever shape `Plus` leaves us with. This is the intended behaviour.

### 6.4 Recursion / cost bound

Stage 1 dominates: `Together` + `Cancel` + `is_zero_poly` can each be
exponential in the worst case. We do **not** add explicit time bounds in
v1 — Mathilda's recursion limit and the existing `IS_ZERO_POLY_MAX_DEPTH=32`
already prevent runaway. If profiling shows pathological inputs, add a
`size_bound` early-out (return UNKNOWN for trees above N nodes).

### 6.5 False positives in Stage 3

Random rational substitution is silently incorrect for identities that
hold only on a measure-zero subset. The standard mitigation, which we
adopt:

- Sample from `Q[i]` not `Q`, so non-real arguments hit branch cuts.
- Use *independent* samples per trial — never reuse symbol→value maps.
- Require **all** k trials to return `True`; a single `False` aborts.

`Sqrt[x^2] - x` evaluates to a non-zero number at almost every sample
(positive `x` gives 0, negative gives `-2x`; complex `x` gives a
branch-dependent value). With four independent samples drawn uniformly
over a domain that includes negatives and complex points, the probability
all four happen to be positive reals is `(1/4)^4 = 1/256` — small but
real, so we **also** require Stage 1 to *not* contradict Stage 3 (which
it does in this case: `Together`/`Cancel` does not collapse `Sqrt[x^2]−x`,
so Stage 1 reports "not proved", and we then defer to Stage 3, which
returns False on roughly 255/256 of runs). The remaining 1/256 risk is
acceptable for `PossibleZeroQ`'s documented "quick but not always
accurate" contract.

### 6.6 Determinism for tests

CI tests must be reproducible. Test fixtures call `SeedRandom[42]`
(builtin already exists in `src/random.c:builtin_seedrandom`) at the top
of every test that exercises Stage 3. Production runs are nondeterministic
by design.

### 6.7 Thread safety

Mathilda is single-threaded today. No locks required. The Schwartz–Zippel
sampler uses the global PRNG (`g_rand_state` in `src/random.c:52`); this
mirrors how every other `Random*` builtin behaves.

---

## 7. Test plan — `tests/test_zero_test.c`

Following the style of `tests/test_comparisons.c` and registered in
`tests/CMakeLists.txt` alongside `comparisons_tests`. Each test case
calls `parse_or_die(...)`, evaluates, calls `PossibleZeroQ`, and asserts
the symbol returned.

### 7.1 Stage 0 — structural

| # | Input | Expected |
|---|-------|----------|
| 1 | `PossibleZeroQ[0]` | True |
| 2 | `PossibleZeroQ[0.0]` | True |
| 3 | `PossibleZeroQ[-0.0]` | True |
| 4 | `PossibleZeroQ[1]` | False |
| 5 | `PossibleZeroQ[10^100]` (bigint) | False |
| 6 | `PossibleZeroQ[Complex[0,0]]` | True |
| 7 | `PossibleZeroQ[Pi]` | False |
| 8 | `PossibleZeroQ[x]` (unbound symbol) | False |

### 7.2 Stage 1 — rational identities

| # | Input | Expected |
|---|-------|----------|
|  9 | `PossibleZeroQ[(x+1)(x-1) - x^2 + 1]` | True |
| 10 | `PossibleZeroQ[(x+y)^3 - x^3 - 3 x^2 y - 3 x y^2 - y^3]` | True |
| 11 | `PossibleZeroQ[1/x + 1/y - (x+y)/(x y)]` | True |
| 12 | `PossibleZeroQ[(x^2-1)/(x-1) - (x+1)]` (after Cancel) | True |
| 13 | `PossibleZeroQ[x + y]` | False |
| 14 | `PossibleZeroQ[x^2 - 2 x + 1 - (x-1)^2]` | True |

### 7.3 Stage 2 — closed numeric

| # | Input | Expected |
|---|-------|----------|
| 15 | `PossibleZeroQ[E^(I Pi/4) - (-1)^(1/4)]` | True |
| 16 | `PossibleZeroQ[(E+Pi)^2 - E^2 - Pi^2 - 2 E Pi]` | True |
| 17 | `PossibleZeroQ[E^Pi - Pi^E]` | False |
| 18 | `PossibleZeroQ[Sin[Pi]]` | True |
| 19 | `PossibleZeroQ[Cos[Pi/2]]` | True |
| 20 | `PossibleZeroQ[Sqrt[2]^2 - 2]` | True |
| 21 | `PossibleZeroQ[Sin[1]^2 + Cos[1]^2 - 1]` | True |
| 22 | `PossibleZeroQ[2^(2 I) - 2^(-2 I) - 2 I Sin[Log[4]]]` | True |
| 23 | `PossibleZeroQ[Log[2]+Log[3]-Log[6]]` | True |
| 24 | `PossibleZeroQ[10^(-30)]` | False (small but non-zero) |
| 25 | `PossibleZeroQ[Sqrt[2] - 1.41421356]` | False |

Catastrophic-cancellation regression test (must NOT report True):

| 26 | `PossibleZeroQ[10^20 + 1 - 10^20]` | False (= 1) |

### 7.4 Stage 3 — symbolic non-rational

| # | Input | Expected |
|---|-------|----------|
| 27 | `SeedRandom[42]; PossibleZeroQ[Sin[x]^2 + Cos[x]^2 - 1]` | True |
| 28 | `SeedRandom[42]; PossibleZeroQ[Sin[2 x] - 2 Sin[x] Cos[x]]` | True |
| 29 | `SeedRandom[42]; PossibleZeroQ[Sqrt[x^2] - x]` | False |
| 30 | `SeedRandom[42]; PossibleZeroQ[Log[E^x] - x]` | True / False — accept either, document |
| 31 | `SeedRandom[42]; PossibleZeroQ[Exp[Log[x]] - x]` | False (branch cut at x<0) |

Test 30 is **intentionally** an under-specified case: principal-branch
`Log[E^x] = x` fails outside the strip `-π < Im(x) ≤ π`. Mathematica's
`PossibleZeroQ` answers `False` here; we will assert whichever answer our
sampling delivers and document it.

### 7.5 Listable

| 32 | `PossibleZeroQ[{0, x-x, 1, Pi-Pi}]` → `{True, True, False, True}` |

### 7.6 Attribute / registration smoke tests

- `Attributes[PossibleZeroQ]` returns `{Listable, Protected}` (order is
  Mathilda's canonical sort).
- `Information[PossibleZeroQ]` prints the docstring without crashing.
- `PossibleZeroQ` is callable from the parser without redefinition.

### 7.7 Memory tests

Add a Valgrind suppression-clean run target:
```
valgrind --leak-check=full --error-exitcode=1 ./tests/build/zero_test_tests
```
to CI. The test driver shall create and free ≥10,000 invocations to give
Valgrind enough volume to detect O(1)-per-call leaks.

### 7.8 Negative / unknown handling

| 33 | `PossibleZeroQ[Erf[Log[4] + 2 Log[Sin[Pi/8]]] - Erf[Log[2-Sqrt[2]]]]` | True + ztest1 message |

(Will only run once `Erf` is implemented. Until then this test is
`#ifdef MATHILDA_HAS_ERF`-gated.)

---

## 8. Phased implementation order

1. **P1 — Scaffolding.** Create `src/zero_test.{c,h}`, register
   `PossibleZeroQ` in `sym_names.{c,h}`, `core.c`, `info.c`, `attr.c`.
   Stub returns `False` always. Land tests 1–8 (Stage 0).
2. **P2 — Stage 1.** Implement rational-form check. Land tests 9–14.
3. **P3 — Stage 2.** Implement precision-ladder numericalize. Land tests
   15–26.
4. **P4 — Stage 3.** Implement Schwartz–Zippel substitution. Land tests
   27–32, with `SeedRandom`.
5. **P5 — Polish.** Emit `ztest1::` message; add `zero_test_decide`
   call-sites in `Equal` (only behind a feature switch, to avoid
   destabilising existing tests in one PR); add the changelog and spec
   entries.

Each phase is independently mergeable, with its own valgrind run.

---

## 9. Messages

Mathilda's message format: see `src/messages.c` (or whichever file owns
`Message[...]`). Define a single new message tag:

```
PossibleZeroQ::ztest1 =
    "Unable to decide whether numeric quantity `1` is equal to zero. "
    "Assuming it is."
```

Emitted only from the trip-wire path (§4.5), so well-typed inputs never
see it.

---

## 10. Known limitations (documented up front)

1. `Erf`, `Gamma`, `BesselJ`, hypergeometric functions: not yet implemented
   in Mathilda, so any input containing them falls straight to the
   trip-wire. Verdict will always be `True` + ztest1 message until those
   builtins land.
2. No PSLQ / integer-relation detection (Phase 3).
3. No structure-theorem reasoning for transcendental extensions
   (Phase 3 — see Bronstein 2005 §3.1).
4. No quantifier elimination / cylindrical algebraic decomposition for
   inequalities — the answer for `Sqrt[x^2] - x` depends on sampling and
   is one-sided.
5. Identities that hold only on a measure-zero set may be missed; this is
   inherent to Monte-Carlo zero-testing.
6. Anything algebraic over a non-trivial extension is currently routed
   through numeric evaluation — once `src/poly/algfac.c` is wired in we
   can decide these exactly.

---

## 11. Open questions for review

- Should `zero_test_decide` be a documented public symbol used by
  `Equal`/`Simplify`, or stay private to `zero_test.c` for v1?
- Tolerance constants in §4.3 (`TOLERANCE_FACTOR = 8`, ladder rungs) —
  default heuristic, may need tuning after the test suite is in place.
- Stage 3 sampling domain: should we also sample from `Q[ω]` for higher
  algebraic order, or is `Q[i]` enough? `Q[i]` matches Mathematica's
  default; recommend keeping it.
- Whether to memoize verdicts on `expr_hash(e)` — `FactorMemo` already
  exists, easy to plug in. Recommend NO for v1: adds a global state
  surface for marginal gain.

---

## 12. References

- Aslaksen, H. (1996). *Can your computer do complex analysis?* The
  Mathematical Intelligencer, 18(3), 50–58.
- Bailey, D. H.; Borwein, P. B.; Plouffe, S. (1997). *On the rapid
  computation of various polylogarithmic constants.* Mathematics of
  Computation, 66(218), 903–913.
- Bronstein, M. (2005). *Symbolic Integration I: Transcendental Functions*
  (2nd ed.). Springer.
- Caviness, B. F. (1970). *On canonical forms and simplification.* JACM,
  17(2), 385–396.
- DeMillo, R. A.; Lipton, R. J. (1978). *A probabilistic remark on
  algebraic program testing.* Information Processing Letters, 7(4),
  193–195.
- Ferguson, H. R. P.; Bailey, D. H. (1992). *A polynomial time, numerically
  stable integer relation algorithm.* RNR Technical Report RNR-91-032.
- Higham, N. J. (2002). *Accuracy and Stability of Numerical Algorithms*
  (2nd ed.). SIAM.
- Richardson, D. (1968). *Some undecidable problems involving elementary
  functions of a real variable.* Journal of Symbolic Logic, 33(4),
  514–520.
- Risch, R. H. (1969). *The problem of integration in finite terms.*
  Transactions of the AMS, 139, 167–189.
- Schwartz, J. T. (1980). *Fast probabilistic algorithms for verification
  of polynomial identities.* JACM, 27(4), 701–717.
- Stoutemyer, D. R. (1989). *Crimes and misdemeanors in the computer
  algebra trade.* Notices of the AMS, 38(7), 778–785.
- Wang, P. S. (1974). *The undecidability of the existence of zeros of
  real elementary functions.* JACM, 21(4), 586–589.
- Zippel, R. (1979). *Probabilistic algorithms for sparse polynomials.*
  Proc. EUROSAM '79, LNCS 72, 216–226.
