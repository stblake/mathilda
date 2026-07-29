# Implementation Plan: `NProduct` (numerical product)

> Self-contained build sheet. All research, design decisions, exact edit sites
> (file + line anchors), full source, and the test suite are captured here so the
> work can be executed in one pass without re-deriving anything.

## 1. Goal & scope

Add a builtin `NProduct[f, {i, imin, imax (, di)}, opts]` to Mathilda's
`src/numerical_calculus/` subsystem — the multiplicative analogue of the
existing `NSum`. It must:

- give numerical approximations to finite **and** infinite products, with an
  optional step `di`;
- support multidimensional products `NProduct[f, {i,…}, {j,…}, …]`, including an
  inner bound that depends on an outer index;
- support complex, alternating, and oscillatory products;
- honour the options `Method`, `WorkingPrecision`, `AccuracyGoal`,
  `PrecisionGoal`, `VerifyConvergence`, `NProductFactors`,
  `NProductExtraFactors`, `WynnDegree`, `EvaluationMonitor`;
- carry attributes `HoldAll, Protected` (the index is `Block`-localised);
- evaluate efficiently at machine precision and arbitrary (MPFR) precision;
- introduce **zero** memory leaks (valgrind clean).

## 2. Key design decision — reduce to `NSum`, don't reimplement

The authoritative source (Keiper, *The N functions of Mathematica*, 1992, §4.2,
p. G.15) states verbatim:

> "For `NProduct[]` the method `Integrate` means that the product is evaluated
> as `Exp[NSum[Log[...], ...]]` where the method for the `NSum` is `Integrate`."

`mpmath.nprod` does the same (`exp(nsum(lambda n: log(f(n)), interval))` at +10
guard bits). Every serious implementation reduces a product to **the sum of the
logarithms, then exponentiates**:

```
∏_{i=imin}^{imax} f(i)  =  exp( Σ_{i=imin}^{imax} log f(i) )
```

Mathilda's `NSum` (`src/numerical_calculus/nsum.c`) is already a complete, tested
engine implementing Euler–Maclaurin (`Integrate`), Wynn-ε / SequenceLimit,
Cohen–Villegas–Zagier (alternating), automatic method selection, divergence
detection, MPFR working precision, multidimensional recursion, and tail-difference
handling of large finite ranges. Reusing it is the **simplest, most faithful,
lowest-risk** design (CLAUDE.md: "Simplicity First"; "C for performance, rules
for mathematics").

### Why the reduction is provably correct (incl. the complex branch-cut worry)

For a *finite* range, with principal `Log`, the identity is **exact** — no
branch-cut error:

```
exp( Σ Log f(i) ) = ∏ exp(Log f(i)) = ∏ f(i)
```

because `exp(Log z) = z` for any branch, and the principal logs' imaginary parts
sum to a multiple of 2π that `exp` wraps away. (This is also why Mathematica/
Mathilda legitimately show a harmless `+0. I` residue on products of real
negative factors — we match that.)

For an *infinite* product, convergence **requires** `f(i) → 1`, hence
`Log f(i) → 0` on the principal branch (continuous near 1). So the tail `NSum`
extrapolates / Euler–Maclaurins is smooth, and all derivatives of `Log f` used
by Euler–Maclaurin are `f'/f` — branch-independent. **No separate
phase-unwrapping code is needed.**

### Bonus: faithful failure modes & convergence semantics

- `VerifyConvergence` (Keiper: "factors are tested for convergence to one") is
  inherited for free — `NSum` tests `Log f(i) → 0`, i.e. `f(i) → 1`.
- `∏(1 + 1/k)` diverges but Mathematica does **not** detect it (ratio test on
  `Log(1+1/k) ≈ 1/k` → 1). Our reduction reproduces this exactly.
- `∏ 2^i`-style divergence → inner `NSum` returns `ComplexInfinity`; we
  special-case it to return `ComplexInfinity` (rather than `Exp[ComplexInfinity]`).

## 3. Algorithm (after option mapping)

```
builtin_nprod(res):
  if not EXPR_FUNCTION or argc < 2: return NULL
  peel trailing NProduct-option Rules            (np_is_option_arg)
  pos_end < 2  -> return NULL                     (need body + >=1 spec)
  body  = args[0]
  for i in 1..pos_end-1: if !np_is_spec(args[i]): return NULL
  parse NProduct options -> NpOpts                (borrowed rhs ptrs; §4)

  // Multidimensional: outer body is an inner NProduct over remaining specs.
  if nspecs >= 2:
      eff_body = NProduct[ body, args[2..pos_end-1], <ORIGINAL option Rules> ]
  else:
      eff_body = copy(body)

  // --- The reduction ---
  logbody = Log[ eff_body ]
  nsum    = NSum[ logbody, copy(spec0), <MAPPED NSum option Rules incl. guard prec> ]
  s       = eval_and_free(nsum)

  if s is ComplexInfinity symbol:             free s; return ComplexInfinity
  if s is not number-like and not an infinity: free s; return NULL   // stay unevaluated
  result = eval_and_free( Exp[s] )            // s ownership moves into Exp node
  if opts.mpfr and result is number-like:     result = eval_and_free( N[result, wdigits] )
  return result
```

Finite-vs-infinite, small-vs-large, method selection, complex/alternating, and
large-finite tail-difference all happen **inside** the delegated `NSum` call.

### Guard precision

`exp` amplifies absolute error in the exponent into relative error in the
product (mpmath adds 10 bits). So, only on the **MPFR path** (user
`WorkingPrecision = p` digits with `p > machine`): run inner `NSum` at
`p + NP_GUARD_DIGITS` (=10) digits, then round the final `Exp[...]` back to `p`
digits with `N[result, p]`. **Machine path**: no guard available; accept machine
(matches Mathematica's machine-precision `NProduct`). Goals (`AccuracyGoal`,
`PrecisionGoal`) pass through unchanged — they are targets, not the guard.

## 4. Options & their mapping to `NSum`

`NProduct` parses its **own** option names, then emits `NSum`-compatible option
`Rule[...]` nodes onto the inner `NSum[]`:

| `NProduct` option       | Default             | Mapped `NSum` option |
|-------------------------|---------------------|----------------------|
| `NProductFactors`       | 15 (NSum default)   | `NSumTerms`          |
| `NProductExtraFactors`  | auto                | `NSumExtraTerms`     |
| `Method` (`"EulerMaclaurin"`, `"WynnEpsilon"`, `Integrate`, `SequenceLimit`, `Automatic`) | `Automatic` | `Method` (NSum already maps these strings/symbols, see `nsum.c:329-356`) |
| `WynnDegree`            | 1                   | `WynnDegree`         |
| `VerifyConvergence`     | `True`              | `VerifyConvergence`  |
| `WorkingPrecision`      | `MachinePrecision`  | `WorkingPrecision` (+ guard on MPFR path) |
| `AccuracyGoal`          | `Infinity`          | `AccuracyGoal`       |
| `PrecisionGoal`         | `Automatic`         | `PrecisionGoal`      |
| `EvaluationMonitor`     | `None`              | dropped (NSum ignores anyway) |

Unknown option ⇒ the Rule is not recognised ⇒ treated as a positional arg ⇒
fails the spec check ⇒ `return NULL` (stay unevaluated), matching `NSum`.

## 5. Exact edit sites (verified against current tree)

**Create 3 files:** `src/numerical_calculus/nprod.h`, `src/numerical_calculus/nprod.c`,
`tests/test_nprod.c` (full source in §7–§9).

**Edit `src/core.c`** — after the `nsum_init();` block at lines **553–554**, insert:
```c
    void nprod_init(void);
    nprod_init();
```

**Edit `src/sym_names.h`** — after line **270** (`extern const char* SYM_NSumExtraTerms;`):
```c
extern const char* SYM_NProduct;
extern const char* SYM_NProductFactors;
extern const char* SYM_NProductExtraFactors;
```

**Edit `src/sym_names.c`** — two sites:
- after line **255** (`const char* SYM_NSumExtraTerms = NULL;`):
```c
const char* SYM_NProduct = NULL;
const char* SYM_NProductFactors = NULL;
const char* SYM_NProductExtraFactors = NULL;
```
- in `sym_names_init`, after line **660** (`SYM_NSumExtraTerms = intern_symbol("NSumExtraTerms");`):
```c
    SYM_NProduct                   = intern_symbol("NProduct");
    SYM_NProductFactors            = intern_symbol("NProductFactors");
    SYM_NProductExtraFactors       = intern_symbol("NProductExtraFactors");
```

**Edit `src/info.c`** — after the `NSum` docstring block (ends line **271**), insert the
`symtab_set_docstring("NProduct", …)` from §10 (terse — no examples, per repo rule).

**Edit `tests/CMakeLists.txt`** — two sites:
- in `COMMON_SRC`, after line **284** (`../src/numerical_calculus/nsum.c`):
```cmake
    ../src/numerical_calculus/nprod.c
```
- after the `nsum_tests` block (lines **364–365**):
```cmake
add_executable(nprod_tests test_nprod.c ${COMMON_SRC})
target_include_directories(nprod_tests PRIVATE ../src ../src/numerical_calculus)
```

**Edit `docs/spec/builtins/calculus.md`** — add a `## NProduct` section beside `## NSum`
(content in §11).

**Edit `docs/spec/changelog/2026-06-08.md`** — add the feature note from §11
(Monday-of-ISO-week file for 2026-06-13 is `2026-06-08.md`).

**Build caveats (from project memory):**
- A new `src/*.c` MUST be in `COMMON_SRC` or *every* `*_tests` binary fails to
  link `_nprod_init`.
- The top-level `makefile` auto-discovers `src/**.c`, so `nprod.c` compiles into
  `./Mathilda` automatically.
- Do **not** `expr_free(res)` in the builtin — the evaluator owns `res`.

## 6. Confirmed API facts (so the source below is correct)

- Interned symbols all exist: `SYM_Log, SYM_Exp, SYM_ComplexInfinity,
  SYM_DirectedInfinity, SYM_Indeterminate, SYM_Infinity, SYM_Rule,
  SYM_RuleDelayed, SYM_List, SYM_Method, SYM_WorkingPrecision, SYM_MachinePrecision,
  SYM_AccuracyGoal, SYM_PrecisionGoal, SYM_VerifyConvergence, SYM_WynnDegree,
  SYM_EvaluationMonitor, SYM_Compiled, SYM_NSumTerms, SYM_NSumExtraTerms,
  SYM_Complex, SYM_Rational, SYM_True, SYM_False`. New: `SYM_NProduct,
  SYM_NProductFactors, SYM_NProductExtraFactors`.
- `expr_new_function(head, args, n)` **copies the args buffer** but **adopts**
  the head and each child `Expr*` (so a stack array is fine; the children are
  owned by the new node). Confirmed by `nsum.c:1135-1141`.
- `eval_and_free(e)` (in `eval.h:71`) evaluates `e`, freeing `e`, returns the
  result. `evaluate(e)` does not free `e`.
- A registered builtin is dispatched by head name, so `expr_new_symbol("NSum")`
  etc. resolve correctly (NSum already does this for its multidim recursion).
- `NSum`/`Exp`/`N`/`Log` builtins all handle `EXPR_MPFR` and complex.
- `NUMERIC_MACHINE_PRECISION_DIGITS` (≈15.95) is the machine/MPFR cutoff
  (`numeric.h:72`).

## 7. `src/numerical_calculus/nprod.h`

```c
/*
 * nprod.h — NProduct[f, {i, imin, imax (, di)}, opts] and multidimensional
 * products.
 *
 * Numerical approximation of a (finite or infinite) product.  Mirrors
 * Mathematica's NProduct: HoldAll, Block-localises the product index, and — per
 * Keiper 1992 — is evaluated as Exp[NSum[Log[f], ...]], reusing the full NSum
 * engine (Euler-Maclaurin / Wynn epsilon / Cohen-Villegas-Zagier, MPFR,
 * divergence detection, multidimensional recursion).
 *
 * Options: Method, WorkingPrecision, NProductFactors, NProductExtraFactors,
 * WynnDegree, VerifyConvergence, AccuracyGoal, PrecisionGoal (EvaluationMonitor
 * accepted and ignored).  Attributes: HoldAll, Protected.
 */
#ifndef MATHILDA_NPROD_H
#define MATHILDA_NPROD_H

#include "expr.h"

Expr* builtin_nprod(Expr* res);
void  nprod_init(void);

#endif /* MATHILDA_NPROD_H */
```

## 8. `src/numerical_calculus/nprod.c`

```c
/*
 * nprod.c — NProduct[f, {i, imin, imax (, di)}, opts]   (see nprod.h)
 *
 * Strategy
 * --------
 * Per Keiper 1992 ("The N functions of Mathematica", G.15), a numerical product
 * is evaluated as the exponential of a numerical sum of logarithms:
 *
 *     Prod_{i=imin}^{imax} f(i)  =  Exp[ NSum[ Log[f(i)], {i, imin, imax} ] ].
 *
 * This is exact for finite ranges (Exp inverts Log on any branch and the
 * principal-log phases that wind by multiples of 2*pi are unwrapped by Exp) and
 * correct for convergent infinite products (factors -> 1 => Log f -> 0 is smooth
 * on the principal branch, and Euler-Maclaurin uses only branch-independent
 * derivatives f'/f).  We therefore delegate every hard part — method selection,
 * Euler-Maclaurin, Wynn epsilon, Cohen-Villegas-Zagier, MPFR working precision,
 * large-finite tail differences, and divergence detection — to the existing,
 * tested NSum engine, and only:
 *   - parse NProduct's own option names and map them onto NSum's;
 *   - add guard digits on the arbitrary-precision path (Exp amplifies the
 *     absolute error of the exponent into relative error of the product);
 *   - handle multidimensional products by recursion (inner NProduct as body);
 *   - special-case divergence (NSum -> ComplexInfinity).
 *
 * Memory: receives `res` owned by the evaluator; returns a fresh Expr* on
 * success or NULL (unevaluated).  Never frees `res`.
 */

#include "nprod.h"

#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "attr.h"
#include "eval.h"
#include "numeric.h"
#include "sym_names.h"
#include "symtab.h"

/* Extra digits carried on the MPFR path before the final Exp/round. */
#define NP_GUARD_DIGITS 10.0

/* ------------------------------------------------------------------ *
 *  Diagnostics                                                        *
 * ------------------------------------------------------------------ */

static void np_warn(const char* tag, const char* fmt, ...) {
    va_list ap;
    fprintf(stderr, "NProduct::%s: ", tag);
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
}

/* ------------------------------------------------------------------ *
 *  Small helpers                                                      *
 * ------------------------------------------------------------------ */

/* Iterator spec: List with 2..4 args whose first is a symbol (mirrors NSum). */
static bool np_is_spec(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL
        || e->data.function.head->data.symbol != SYM_List) return false;
    size_t n = e->data.function.arg_count;
    return n >= 2 && n <= 4 && e->data.function.args[0]->type == EXPR_SYMBOL;
}

static bool np_is_known_option(const char* s) {
    return s == SYM_Method || s == SYM_WorkingPrecision
        || s == SYM_NProductFactors || s == SYM_NProductExtraFactors
        || s == SYM_WynnDegree || s == SYM_VerifyConvergence
        || s == SYM_AccuracyGoal || s == SYM_PrecisionGoal
        || s == SYM_Compiled || s == SYM_EvaluationMonitor;
}

static bool np_is_option_arg(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* lhs = e->data.function.args[0];
    return lhs->type == EXPR_SYMBOL && np_is_known_option(lhs->data.symbol);
}

/* WorkingPrecision values are small integer/real literals in practice. */
static bool np_to_double(Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real;            return true; }
    return false;
}

static bool np_is_numberlike(Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_REAL:
        case EXPR_BIGINT:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        case EXPR_FUNCTION: {
            Expr* h = e->data.function.head;
            return h->type == EXPR_SYMBOL
                && (h->data.symbol == SYM_Complex || h->data.symbol == SYM_Rational);
        }
        default:
            return false;
    }
}

static bool np_is_infinity(Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL)
        return e->data.symbol == SYM_ComplexInfinity
            || e->data.symbol == SYM_Indeterminate
            || e->data.symbol == SYM_Infinity;
    if (e->type == EXPR_FUNCTION) {
        Expr* h = e->data.function.head;
        return h->type == EXPR_SYMBOL && h->data.symbol == SYM_DirectedInfinity;
    }
    return false;
}

/* Build Rule[Symbol(name), rhs] adopting rhs. */
static Expr* np_rule(const char* name, Expr* rhs) {
    Expr* a[2];
    a[0] = expr_new_symbol(name);
    a[1] = rhs;
    return expr_new_function(expr_new_symbol("Rule"), a, 2);
}

/* Wrap head[child] adopting child. */
static Expr* np_apply1(const char* head, Expr* child) {
    Expr* a[1];
    a[0] = child;
    return expr_new_function(expr_new_symbol(head), a, 1);
}

/* ------------------------------------------------------------------ *
 *  Options                                                            *
 * ------------------------------------------------------------------ */

typedef struct {
    Expr* factors;       /* NProductFactors rhs (borrowed) or NULL      */
    Expr* extra_factors; /* NProductExtraFactors rhs or NULL            */
    Expr* method;        /* Method rhs or NULL                          */
    Expr* wynn;          /* WynnDegree rhs or NULL                      */
    Expr* verify;        /* VerifyConvergence rhs or NULL               */
    Expr* wprec;         /* WorkingPrecision rhs or NULL                */
    Expr* accgoal;       /* AccuracyGoal rhs or NULL                    */
    Expr* precgoal;      /* PrecisionGoal rhs or NULL                   */
    bool   mpfr;         /* WorkingPrecision selects arbitrary precision */
    double wdigits;      /* requested working-precision digits (mpfr)   */
} NpOpts;

static bool np_apply_option(Expr* rule, NpOpts* o) {
    Expr* lhs = rule->data.function.args[0];
    Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol;

    if (name == SYM_NProductFactors)      { o->factors = rhs;       return true; }
    if (name == SYM_NProductExtraFactors) { o->extra_factors = rhs; return true; }
    if (name == SYM_Method)               { o->method = rhs;        return true; }
    if (name == SYM_WynnDegree)           { o->wynn = rhs;          return true; }
    if (name == SYM_VerifyConvergence)    { o->verify = rhs;        return true; }
    if (name == SYM_AccuracyGoal)         { o->accgoal = rhs;       return true; }
    if (name == SYM_PrecisionGoal)        { o->precgoal = rhs;      return true; }
    if (name == SYM_WorkingPrecision) {
        o->wprec = rhs;
        o->mpfr = false; o->wdigits = 0.0;
        if (!(rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_MachinePrecision)) {
            double d;
            if (np_to_double(rhs, &d) && d > NUMERIC_MACHINE_PRECISION_DIGITS) {
                o->mpfr = true; o->wdigits = d;
            }
        }
        return true;
    }
    /* Compiled / EvaluationMonitor: accepted, ignored. */
    if (name == SYM_Compiled || name == SYM_EvaluationMonitor) return true;
    return false;
}

/* ------------------------------------------------------------------ *
 *  Entry point                                                        *
 * ------------------------------------------------------------------ */

Expr* builtin_nprod(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;
    Expr** args = res->data.function.args;

    /* Peel trailing options. */
    size_t pos_end = argc;
    while (pos_end > 0 && np_is_option_arg(args[pos_end - 1])) pos_end--;
    if (pos_end < 2) return NULL;                       /* need f + >=1 spec */

    Expr* body = args[0];
    for (size_t i = 1; i < pos_end; i++)
        if (!np_is_spec(args[i])) return NULL;

    NpOpts o; memset(&o, 0, sizeof o);
    for (size_t i = pos_end; i < argc; i++)
        if (!np_apply_option(args[i], &o)) return NULL;

    size_t nspecs = pos_end - 1;
    Expr* spec0 = args[1];

    /* Effective body of the outer product: for >1 spec it is an inner NProduct
     * over the remaining specs (HoldAll + index localisation lets a dependent
     * inner bound such as {k,1,n} see the outer index). The inner NProduct
     * receives the ORIGINAL options (it maps them itself). */
    Expr* eff_body;
    if (nspecs >= 2) {
        size_t n = 1 /*body*/ + (nspecs - 1) /*specs[2..]*/ + (argc - pos_end) /*opts*/;
        Expr** v = malloc(sizeof(Expr*) * (n ? n : 1));
        if (!v) return NULL;
        size_t w = 0;
        v[w++] = expr_copy(body);
        for (size_t i = 2; i < pos_end; i++) v[w++] = expr_copy(args[i]);
        for (size_t i = pos_end; i < argc; i++) v[w++] = expr_copy(args[i]);
        eff_body = expr_new_function(expr_new_symbol("NProduct"), v, w);
        free(v);
    } else {
        eff_body = expr_copy(body);
    }

    /* Build the inner NSum[Log[eff_body], spec0, <mapped options>]. */
    Expr* logbody = np_apply1("Log", eff_body);

    Expr** nv = malloc(sizeof(Expr*) * 12);
    if (!nv) { expr_free(logbody); return NULL; }
    size_t nw = 0;
    nv[nw++] = logbody;                 /* adopts eff_body                  */
    nv[nw++] = expr_copy(spec0);
    if (o.factors)       nv[nw++] = np_rule("NSumTerms",      expr_copy(o.factors));
    if (o.extra_factors) nv[nw++] = np_rule("NSumExtraTerms", expr_copy(o.extra_factors));
    if (o.method)        nv[nw++] = np_rule("Method",         expr_copy(o.method));
    if (o.wynn)          nv[nw++] = np_rule("WynnDegree",     expr_copy(o.wynn));
    if (o.verify)        nv[nw++] = np_rule("VerifyConvergence", expr_copy(o.verify));
    if (o.accgoal)       nv[nw++] = np_rule("AccuracyGoal",   expr_copy(o.accgoal));
    if (o.precgoal)      nv[nw++] = np_rule("PrecisionGoal",  expr_copy(o.precgoal));
    if (o.mpfr) {
        long inner = (long)(o.wdigits + NP_GUARD_DIGITS + 0.5);
        nv[nw++] = np_rule("WorkingPrecision", expr_new_integer(inner));
    } else if (o.wprec) {
        nv[nw++] = np_rule("WorkingPrecision", expr_copy(o.wprec));
    }

    Expr* nsum = expr_new_function(expr_new_symbol("NSum"), nv, nw);
    free(nv);

    Expr* s = eval_and_free(nsum);      /* the numeric Log-sum, or special  */

    /* Divergent product: NSum reports the divergent log-sum as ComplexInfinity. */
    if (s && s->type == EXPR_SYMBOL && s->data.symbol == SYM_ComplexInfinity) {
        expr_free(s);
        return expr_new_symbol("ComplexInfinity");
    }
    /* Could not reduce to a number (symbolic body / bound): stay unevaluated. */
    if (!s || (!np_is_numberlike(s) && !np_is_infinity(s))) {
        expr_free(s);
        return NULL;
    }

    Expr* result = eval_and_free(np_apply1("Exp", s));   /* s -> Exp node    */

    /* Round back to the requested working precision on the MPFR path. */
    if (o.mpfr && np_is_numberlike(result)) {
        Expr* a[2];
        a[0] = result;
        a[1] = expr_new_integer((long)(o.wdigits + 0.5));
        result = eval_and_free(expr_new_function(expr_new_symbol("N"), a, 2));
    }

    (void)np_warn;   /* reserved for future divergence/convergence messages  */
    return result;
}

/* ------------------------------------------------------------------ *
 *  Registration                                                       *
 * ------------------------------------------------------------------ */

void nprod_init(void) {
    symtab_add_builtin("NProduct", builtin_nprod);
    /* HoldAll: the factor and the iterator specs must not be pre-evaluated; the
     * index is Block-localised inside the delegated NSum. Not Listable. */
    symtab_get_def("NProduct")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
}
```

> **Field-name check before compiling:** confirm `Expr.data.integer`,
> `Expr.data.real`, `EXPR_BIGINT`, and (under `USE_MPFR`) `EXPR_MPFR` names
> against `expr.h` — `nsum.c:75-122` uses the same accessors, so copy from there
> if any differ. `np_to_double` deliberately only reads INTEGER/REAL (a
> WorkingPrecision literal), so it does not depend on the bigint/mpfr field
> names.

## 9. `tests/test_nprod.c`

Mirror `test_nsum.c` exactly (same `eval_str`, `close_to`, `ASSERT_CLOSE`,
`main` scaffolding — compares **inside** the language, never parses printed
reals). Cases:

```c
/* 1. Infinite, Euler-Maclaurin (default). Keiper: ~ Sinh[Pi]/Pi. */
ASSERT_CLOSE("NProduct[1+1/i^2,{i,1,Infinity}]", "Sinh[Pi]/Pi", 1e-7);
ASSERT_CLOSE("NProduct[1+1/n^2,{n,1,Infinity}]", "Sinh[Pi]/Pi", 1e-7);

/* 2. Finite (positive factors -> real, no spurious imaginary). */
ASSERT_CLOSE("NProduct[1+(-1)^i/i^2,{i,100,10000}]", "1.00005", 1e-4);
ASSERT_CLOSE("NProduct[i,{i,1,5}]", "120", 1e-9);          /* exact 5! */
ASSERT_CLOSE("NProduct[1+1/i,{i,1,9}]", "10", 1e-9);       /* telescopes to 10 */

/* 3. Step di, and the equivalent reindexing. */
ASSERT_CLOSE("NProduct[1+1/2^i,{i,0,Infinity,2}]", "2.71182", 1e-4);
ASSERT_CLOSE("NProduct[1+1/2^(2j),{j,0,Infinity}]", "2.71182", 1e-4);

/* 4. Explicit Wynn epsilon method. */
ASSERT_CLOSE("NProduct[1+1/(n+1)^2,{n,0,Infinity},Method->\"WynnEpsilon\"]",
             "Sinh[Pi]/Pi", 1e-2);

/* 5. Multidimensional, independent and dependent inner bounds. */
ASSERT_CLOSE("NProduct[1+(-1)^n(2/n)^k/k^2,{n,2,Infinity},{k,1,Infinity}]",
             "0.812331", 1e-4);
ASSERT_CLOSE("NProduct[1+(-1)^n(2/n)^k/k^2,{n,2,Infinity},{k,1,n}]",
             "0.564097", 1e-4);

/* 6. Complex infinite product. Keiper: 1.43706 + 1.07945 I. */
ASSERT_CLOSE("NProduct[1+E^(I n 2/3)/n^2,{n,1,Infinity}]",
             "1.43706 + 1.07945 I", 1e-4);

/* 7. Goals tighten the error vs the closed form. */
ASSERT_CLOSE("NProduct[1+1/n^2,{n,1,Infinity},AccuracyGoal->8,PrecisionGoal->8]",
             "Sinh[Pi]/Pi", 1e-9);

/* 8. NProductFactors improves the peaked case. */
ASSERT_CLOSE("NProduct[1+1/(1+(k-20)^2),{k,0,Infinity},NProductFactors->30]",
             "12.9013", 1e-3);

/* 9. Sin product identity. */
ASSERT_CLOSE("Block[{x=1/2}, x NProduct[1-x^2/(Pi^2 k^2),{k,1,Infinity}]]",
             "Sin[0.5]", 1e-5);

/* 10. Arbitrary precision (per project rule: numeric builtins need an N[...,p] case). */
ASSERT_CLOSE("NProduct[E^(-1/(2n))(1+1/(2n)),{n,1,Infinity},WorkingPrecision->25]",
             "0.8455012816335291822442405", 1e-12);

/* 11. Edge / unevaluated / attributes. */
//  - symbolic upper bound stays unevaluated (strstr "NProduct[")
//  - 1-arg NProduct[f] stays unevaluated
//  - non-list spec NProduct[f,i] stays unevaluated
//  - unknown option NProduct[1+1/i^2,{i,1,Infinity},Bogus->1] stays unevaluated
//  - MemberQ[Attributes[NProduct], Protected] / HoldAll  == True
//  - memory loop (parse/evaluate/free over the inputs above)
```

> Verify each expected constant numerically when first run (e.g. evaluate
> `Sinh[Pi]/Pi`, the multidim targets, and the WP-25 value in the REPL) and
> adjust tolerances if the engine's accuracy differs; treat any case needing a
> looser tolerance than its NSum analogue as a signal to investigate, not to
> paper over.

## 10. Docstring (`src/info.c`, terse — no examples)

```c
    symtab_set_docstring("NProduct",
        "NProduct[f, {i, imin, imax}]\n"
        "\tgives a numerical approximation to the product of f for i from imin "
        "to imax.\n\n"
        "NProduct[f, {i, imin, imax, di}] uses step di. imax may be Infinity. "
        "NProduct[f, {i, ...}, {j, ...}, ...] evaluates a multidimensional "
        "product (an inner bound may depend on an outer index). The index is "
        "localised (HoldAll). Evaluated as Exp[NSum[Log[f], ...]], so the NSum "
        "engine (Euler-Maclaurin for monotone factors, Wynn's epsilon "
        "otherwise) and its convergence test carry over. Machine or arbitrary "
        "precision via WorkingPrecision.\n\n"
        "Options: Method (Automatic | EulerMaclaurin | WynnEpsilon), "
        "WorkingPrecision (default MachinePrecision), NProductFactors (leading "
        "factors taken explicitly, default 15), NProductExtraFactors, "
        "WynnDegree, VerifyConvergence (default True; a divergent product gives "
        "ComplexInfinity), AccuracyGoal, PrecisionGoal.");
```

## 11. Docs & changelog

**`docs/spec/builtins/calculus.md`** — add beside `## NSum`:

```markdown
## NProduct

Numerical multiplication.  `NProduct[f, {i, imin, imax}]` gives a numerical
approximation to the product of `f` for `i` from `imin` to `imax` (which may be
`Infinity`); `NProduct[f, {i, imin, imax, di}]` uses step `di`, and
`NProduct[f, {i,…}, {j,…}, …]` is multidimensional (an inner bound may depend on
an outer index).  Implemented in `src/numerical_calculus/nprod.{c,h}`.
Attributes: `HoldAll, Protected`.

Per Keiper (1992) the product is evaluated as `Exp[NSum[Log[f], …]]`, reusing
the full NSum engine: Euler-Maclaurin (`Method -> "EulerMaclaurin"`, default for
monotone factors), Wynn epsilon (`Method -> "WynnEpsilon"`), automatic method
selection, MPFR working precision, large-finite tail differences, and the
convergence test (factors are checked for `-> 1`; a divergent product such as
`∏(1+2^i)` returns `ComplexInfinity`).  Options: `Method`, `WorkingPrecision`,
`NProductFactors` (leading factors taken explicitly, default 15),
`NProductExtraFactors`, `WynnDegree`, `VerifyConvergence`, `AccuracyGoal`,
`PrecisionGoal`.  On the arbitrary-precision path NProduct carries guard digits
because `Exp` turns the exponent's absolute error into the product's relative
error.

Like Mathematica, NProduct can miss the divergence of slowly diverging products
(e.g. `∏(1+1/k)`, whose log-sum is the harmonic series) and may leave a harmless
`+0. I` residue on products of real negative factors.
```

**`docs/spec/changelog/2026-06-08.md`** — append:

```markdown
## Feature: `NProduct` — numerical product (2026-06-13)

New builtin `NProduct[f, {i, imin, imax (, di)}]` giving a numerical
approximation to a finite or infinite product, plus multidimensional
`NProduct[f, {i,…}, {j,…}, …]` (an inner bound may depend on an outer index).
Sixth member of the **numerical calculus** group; implemented in
`src/numerical_calculus/nprod.{c,h}`.  Following Keiper (1992) it is evaluated as
`Exp[NSum[Log[f], …]]`, delegating method selection, Euler-Maclaurin / Wynn
epsilon / Cohen-Villegas-Zagier, MPFR working precision, large-finite tail
differences, and divergence detection to the existing NSum engine, with guard
digits added on the arbitrary-precision path.  Attributes `HoldAll, Protected`
(the index is `Block`-localised inside NSum); options `Method`,
`WorkingPrecision`, `NProductFactors`, `NProductExtraFactors`, `WynnDegree`,
`VerifyConvergence`, `AccuracyGoal`, `PrecisionGoal` (`EvaluationMonitor`
accepted and ignored).  Symbols `NProduct`, `NProductFactors`,
`NProductExtraFactors` interned in `sym_names.{c,h}`; registered from `core_init`
via `nprod_init`; docstring in `info.c`; documented in
`docs/spec/builtins/calculus.md`; unit tests in `tests/test_nprod.c`.
```

## 12. Build & verify (run **foreground only**, one process at a time)

> The previous session leaked `./Mathilda` processes and left 23h+ `until`-loop
> poller shells that exhausted application memory. Do **not** background builds
> or REPL sessions; do **not** write until-loop pollers. Build in the
> foreground, let it finish, then run the one test binary.

```bash
make -j$(sysctl -n hw.ncpu)                  # builds ./Mathilda (auto-discovers nprod.c)
cd tests && mkdir -p build && cd build
cmake .. >/dev/null && make nprod_tests nsum_tests -j$(sysctl -n hw.ncpu)
./nprod_tests                                # scoped: only the new + NSum regression
./nsum_tests
```

Memory hygiene (short, foreground; never background a valgrind on the REPL):
```bash
printf 'NProduct[1+1/i^2,{i,1,Infinity}]\nNProduct[1+E^(I n 2/3)/n^2,{n,1,Infinity}]\nQuit[]\n' \
  | valgrind --leak-check=full ./Mathilda 2>/tmp/np.val
# diff against the known macOS baseline noise (~12.8 kB / 400 blocks from dyld/
# libobjc/Accelerate init); grep the stacks for Mathilda src to find real leaks.
```

## 13. Risks & mitigations

- **`exp` error amplification** → guard digits on the MPFR path (§3), `N`-round back.
- **Divergent products** (`∏(1+1/k)`) — undetectable by the ratio test; we match
  Mathematica (no false guarantee). Documented as a known limitation.
- **Zero / negative finite factors** — the reduction stays correct (`+0. I`
  residue only, matching Mathematica). An exact-zero factor makes the product 0
  via `Log[0]` → `-∞`/`ComplexInfinity` handling.
- **Leaks** — every constructed `Expr*` (eff_body, logbody, NSum wrapper, Exp
  wrapper, `N` wrapper, the malloc'd arg arrays) is paired with
  `eval_and_free`/`expr_free`/`free`; follow the `builtin_nsum` ownership
  pattern. Do **not** `expr_free(res)`.

## 14. Deferred (out of scope)

- Native direct-multiply fast path for small finite products (avoids the log/exp
  round-trip and the `+0. I` residue). Not needed for correctness.
- `EvaluationMonitor` actually firing per factor evaluation.
- A dedicated convergence/divergence **message** (the `np_warn` hook is in place
  but unused).
```
