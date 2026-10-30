/* integrate.c
 *
 * `Integrate[f, x]` System` dispatcher.  Cascades through three
 * stages and supports an explicit `Method -> "..."` option:
 *
 *   1. Integrate`BronsteinRational      — polynomial / rational integrands
 *   2. Integrate`RischTranscendental    — recursive transcendental Risch
 *   3. Integrate`CRCTable               — CRC integral table (lazy-loaded)
 *   4. Integrate`ParallelMixedTower     — parallel Risch-Norman over a mixed
 *                                         radical tower (last resort; .m)
 *
 * Method values: "Automatic" (default, full cascade), "BronsteinRational",
 * "RischTranscendental", "CRCTable", "ParallelMixedTower" (strict passthrough,
 * no fallback).  RischNorman / RischNormanBlake were removed in v0.163 --
 * ParallelMixedTower subsumes both.
 *
 * The CRC table is large and most sessions never need it, so its
 * .m file is Get-loaded on first invocation of try_crctable() rather
 * than at startup.  See LAZY_LOAD_CRC below.
 */

/* clock_gettime / CLOCK_MONOTONIC (used by ip_mono_seconds, the cascade
 * profiler's clock) are POSIX, not C99.  glibc hides them under -std=c99 unless
 * a feature-test macro is defined BEFORE the first include; placed below the
 * first #include it has no effect at all (see SPEC.md §10). */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "integrate.h"
#include "mtime.h"   /* mth_mono_seconds */
#include "core.h"    /* tc_register_unwind_counter */
#include "integrate_interp.h"
#include "integrate_unknown.h"
#include "integrate_derivdivides.h"
#include "integrate_linrad.h"
#include "integrate_quadrad.h"
#include "integrate_linratiorad.h"
#include "integrate_chebychev.h"
#include "integrate_goursat.h"
#include "integrate_fresnel.h"
#include "integrate_gammapower.h"
#include "integrate_logbyparts.h"
#include "integrate_jeffrey.h"
#include "integrate_linarg.h"
#include "integrate_newton_leibniz.h"
#include "integrate_symmetry.h"
#include "integrate_piecewise.h"
#include "integrate_beta.h"
#include "integrate_line.h"
#include "integrate_residue.h"
#include "integrate_dirac.h"
#include "integrate_diffunderint.h"
#include "integrate_ramanujan.h"
#include "integrate_intrep.h"
#include "integrate_series.h"
#include "integrate_king.h"
#include "integrate_feynparam.h"
#include "intrat.h"
#include "integrate_risch_transcendental.h"
#include "risch_canonical.h"
#include "risch_structure.h"
#include "risch_hypertangent.h"
#include "risch_coupled.h"
#include "risch_hermite.h"
#include "intsimp.h"
#include "common.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"
#include "attr.h"
#include "internal.h"
#include "loadmodule.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "message.h"   /* mth_msg_suppress_push/pop: quiet ParallelMixedTower's internal probes */
#include "print.h"     /* expr_to_string: naming the integrand in the Integrate::nonelem diagnostic */
#include "series.h"
#include "functionexpand.h" /* fe_reduce_definite_result: MeijerG / PolyLog closed forms */

#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>    /* clock_gettime / clock: the cascade profiler's clock */

/* True iff `e` is the symbol `True`.  The PolynomialQ / rationalQ
 * predicates we call below return either True or False. */
static bool is_true_symbol(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_True;
}

/* True iff `f` contains `x` as a subexpression (structural FreeQ negation). */
static bool depends_on_var(const Expr* f, const Expr* x) {
    if (expr_eq((Expr*)f, (Expr*)x)) return true;
    if (f->type == EXPR_FUNCTION) {
        if (depends_on_var(f->data.function.head, x)) return true;
        for (size_t i = 0; i < f->data.function.arg_count; i++)
            if (depends_on_var(f->data.function.args[i], x)) return true;
    }
    return false;
}

/* True iff `f` contains a Power[b, e] whose base involves `x` and whose
 * exponent is not a concrete number -- a symbolic-exponent power such as x^k,
 * x^(k-1) or (1-x)^l.  Such an integrand is provably not a rational function of
 * x, yet feeding it to the Together / PolynomialGCD rationality probe (or the
 * derivative-divides quotient normalisation) drives a pseudo-remainder-sequence
 * blow-up over the symbolic powers as independent generators.  Callers use it
 * to decline those Together-backed paths structurally, up front. */
static bool has_symbolic_power_in(const Expr* f, const Expr* x) {
    if (!f || f->type != EXPR_FUNCTION) return false;
    if (head_is((Expr*)f, SYM_Power) && f->data.function.arg_count == 2) {
        const Expr* b = f->data.function.args[0];
        const Expr* e = f->data.function.args[1];
        if (depends_on_var(b, x) && !expr_is_numeric_like(e)) return true;
    }
    if (has_symbolic_power_in(f->data.function.head, x)) return true;
    for (size_t i = 0; i < f->data.function.arg_count; i++)
        if (has_symbolic_power_in(f->data.function.args[i], x)) return true;
    return false;
}

/* Test whether `f` is a polynomial in `x`.  Calls the existing
 * PolynomialQ builtin so we get the same definition the rest of
 * Mathilda uses. */
static bool is_polynomial_in(Expr* f, Expr* x) {
    Expr* args[2] = { expr_copy(f), expr_copy(x) };
    Expr* call = internal_polynomialq(args, 2);
    Expr* val  = evaluate(call);
    expr_free(call);
    bool ok = is_true_symbol(val);
    expr_free(val);
    return ok;
}

/* Test whether `f` is a rational function in `x`: Together[f] must
 * have a non-trivial denominator that is itself polynomial in `x`,
 * with a polynomial numerator. */
static bool is_rational_in(Expr* f, Expr* x) {
    Expr* together = internal_together((Expr*[]){expr_copy(f)}, 1);
    Expr* combined = evaluate(together);
    expr_free(together);
    if (!combined) return false;

    Expr* num = internal_numerator((Expr*[]){expr_copy(combined)}, 1);
    Expr* den = internal_denominator((Expr*[]){expr_copy(combined)}, 1);
    expr_free(combined);
    Expr* num_v = evaluate(num);
    Expr* den_v = evaluate(den);
    expr_free(num); expr_free(den);

    bool ok = false;
    if (num_v && den_v) {
        bool num_is_poly = is_polynomial_in(num_v, x);
        bool den_is_poly = is_polynomial_in(den_v, x);
        bool den_is_unit = (den_v->type == EXPR_INTEGER && den_v->data.integer == 1);
        ok = num_is_poly && den_is_poly && !den_is_unit;
    }
    if (num_v) expr_free(num_v);
    if (den_v) expr_free(den_v);
    return ok;
}

/* True iff `result` is the unevaluated call `head_name[...]`. */
static bool result_is_unresolved(const Expr* result, const char* head_name) {
    return head_is(result, intern_symbol(head_name));
}

/* True iff any subexpression of `e` is an unevaluated call `head_name[...]`.
 * Needed for the CRC pipeline: the inner table rules can fire partially
 * (e.g. extracting a `-1` via `IntegrateTable[a_ f_, x_] := a IntegrateTable[f, x]`)
 * and leave a stray `IntegrateTable[...]` nested inside a Times / Plus.
 * The top-level head_is check would miss those, so we walk the tree. */
static bool result_contains_head(const Expr* e, const char* head_name) {
    if (!e) return false;
    if (head_is(e, intern_symbol(head_name))) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (result_contains_head(e->data.function.head, head_name)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (result_contains_head(e->data.function.args[i], head_name)) return true;
    }
    return false;
}

/* Helper: build and evaluate `head_name[f, x]`, freeing the call
 * expression.  Returns whatever evaluate() produces. */
static Expr* call_stage(const char* head_name, Expr* f, Expr* x) {
    Expr* call = expr_new_function(
        expr_new_symbol(head_name),
        (Expr*[]){ expr_copy(f), expr_copy(x) }, 2);
    Expr* result = evaluate(call);
    expr_free(call);
    return result;
}

static void   integrate_arm_budget(double seconds);

/* NOTE -- a stage-boundary `TimeConstrained[Integrate`DerivativeDivides[f,x], t]`
 * wrapper was built here and REMOVED. It did make TimeConstraint a hard bound
 * (3 s measured 3.00 s, 1 s measured 1.00 s), but it changed the stage's RESULT:
 * `Sqrt[Sin[x]]/(Sin[x]^2 + 1)`, which the recursion gate in
 * integrate_derivdivides.c closes in 0.066 s, became a 3 s DECLINE under the
 * wrapper. Something about reaching the stage through its registered head inside
 * TimeConstrained -- the fail-memo keying, or the argument-evaluation order --
 * loses the close, so the wrapper cost an answer to buy a bound. Not reinstated
 * without understanding that; the recursion gate removes the need in this case
 * (20 s abort -> 0.066 s close) and TimeConstraint remains the between-stages
 * backstop described with g_intg_deadline below. */

/* Stage 0: undefined-function integrator (Roach 1992, §1.7).  Handles
 * integrands rational in unknown functions u[x] and their derivatives,
 * e.g. Integrate[x f'[x] + f[x], x] -> x f[x].  Cheaply gated: returns
 * NULL immediately unless the integrand contains an undefined-function
 * derivative, so ordinary integrands skip it. */
static Expr* try_undefined(Expr* f, Expr* x) {
    return integrate_unknown_try(f, x);
}

/* Stage 1: rational-function integrator.  Returns the antiderivative on
 * success, NULL when the integrand is non-rational or pmint gives up. */
static Expr* try_rational(Expr* f, Expr* x) {
    /* A symbolic-exponent power of x (x^k, (1-x)^l, ...) is never rational in
     * x; decline structurally so the Together-based rationality probe below is
     * never handed an integrand that sends PolynomialGCD into a blow-up. */
    if (has_symbolic_power_in(f, x)) return NULL;
    if (!is_polynomial_in(f, x) && !is_rational_in(f, x)) return NULL;
    Expr* result = call_stage("Integrate`BronsteinRational", f, x);
    if (!result) return NULL;
    if (result_is_unresolved(result, "Integrate`BronsteinRational")) {
        expr_free(result);
        return NULL;
    }
    return result;
}

/* True for a pseudo-elliptic F/R^p integrand (Sqrt/odd-root of a degree>=3
 * polynomial) -- Goursat's domain.  Defined with the linearity helpers below;
 * forward-declared here for the derivative-divides gate. */
static bool has_pseudoelliptic_radical(const Expr* e, const Expr* x);

/* Stage 1b: derivative-divides substitution.  Recognises integrands of the
 * shape c h(u(x)) u'(x) and reduces to Integrate[h[u], u].  Runs the fast,
 * branch-correct direct-quotient strategy first and then the more thorough
 * Eliminate/Solve branch-search (which closes radical substitutions such as
 * u = Sqrt[Tan[x]]); the latter's Eliminate diagnostics are muted while the
 * integrator drives it (see integrate_derivdivides.c).
 *
 * A pseudo-elliptic F/R^p integrand (R a degree>=3 polynomial) is Goursat's
 * domain, not derivative-divides': the Eliminate/Solve branch-search grinds
 * for seconds on such a radicand (e.g. (t-1)/((t+2) Sqrt[t^3-1])) and closes
 * nothing.  For those, run the cheap direct-quotient strategy only
 * (integrate_derivdivides_try) -- it still catches the genuine folds like
 * t^2/(t^3-1)^(1/3) via Pass 1 -- and leave the intact integrand for the
 * Goursat stage.  Non-pseudo-elliptic integrands keep the full search. */
static Expr* try_derivdivides(Expr* f, Expr* x) {
    if (has_pseudoelliptic_radical(f, x))
        return integrate_derivdivides_try(f, x);
    return integrate_derivdivides_full(f, x);
}

/* Stage 1c: linear-radical substitution.  Recognises a rational function of x
 * and radicals (a x + b)^(m/n) of one shared linear argument, rationalises via
 * u = (a x + b)^(1/n), integrates the rational result, and verifies. */
static Expr* try_linrad(Expr* f, Expr* x) {
    return integrate_linrad_try(f, x);
}

/* Stage 1d: quadratic-radical (Euler) substitution.  Recognises a rational
 * function of x and square roots (a x^2 + b x + c)^(m/2) of one shared quadratic
 * argument, applies a single real-valued Euler substitution, integrates the
 * rational result, and back-substitutes. */
static Expr* try_quadrad(Expr* f, Expr* x) {
    return integrate_quadrad_try(f, x);
}

/* Stage 1e: linear-ratio-radical (Möbius) substitution.  Recognises a rational
 * function of x and radicals ((a x + b)/(c x + d))^(m/n) of one shared
 * linear-fractional argument, rationalises via u = ((a x + b)/(c x + d))^(1/n),
 * integrates the rational result, and back-substitutes. */
static Expr* try_linratiorad(Expr* f, Expr* x) {
    return integrate_linratiorad_try(f, x);
}

/* Stage 1e2: Chebychev binomial differential.  Recognises an integrand
 * x^p (a x^r + b)^q with p, q, r rational and a, b free of x, and -- when one
 * of q, (p+1)/r, q+(p+1)/r is an integer (Chebychev's theorem) -- applies the
 * matching rationalising substitution (x = u^N, u^s = a x^r + b, or
 * u = x^r then t^s = (a u + b)/u), integrates the rational result, and
 * back-substitutes.  Recognition is a single structural scan, so it is cheap
 * enough to run ahead of the Eliminate/Solve search in derivative-divides.
 * Non-elementary binomials return NULL, so the cascade falls through to the
 * later methods (which may one day carry special-function representations). */
static Expr* try_chebychev(Expr* f, Expr* x) {
    return integrate_chebychev_try(f, x);
}

/* Stage 1f: continuous Weierstrass substitution (Jeffrey & Rich 1994).
 * Recognises a rational function of the trig kernels Sin/Cos/Tan/Cot/Sec/Csc[x]
 * (or the hyperbolic kernels) carrying a kernel in a denominator, substitutes
 * u = Tan[x/2] (Tanh[x/2] for hyperbolic), integrates the resulting rational
 * function of u, back-substitutes, and -- for the trig case -- adds the
 * K Floor[(x - b)/p] secular term that removes the spurious discontinuities at
 * the poles of Tan[x/2].  Runs ahead of Risch-Norman so genuine rational-trig
 * integrands get the real, continuous antiderivative rather than pmint's
 * complex-logarithm form. */
/* c Sec[x]^n / c Csc[x]^n, n >= 3 an integer, c free of x: the reduction formula
 *     Integrate[Sec^n] = Sec^(n-2) Tan/(n-1) + (n-2)/(n-1) Integrate[Sec^(n-2)]
 *     Integrate[Csc^n] = -Csc^(n-2) Cot/(n-1) + (n-2)/(n-1) Integrate[Csc^(n-2)]
 * recursing down to Integrate[Sec] / Integrate[Sec^2] (and Csc), which close.
 * These were the one trig-power family the automatic cascade could not do: a POSITIVE
 * power of a reciprocal kernel carries no explicit negative exponent, so the
 * Weierstrass gate (kernel in a denominator) passed it on and DerivativeDivides spun
 * past any budget -- Integrate[Sec[x]^4, x] hung, and with it the Bernoulli
 * y' == (Tan[x] + Sec[x] y^3) y of §2.2.48-4702.  Every Tan/Cot power and mixed
 * Sec^m Tan^k product already closed instantly and is not matched here. */
static Expr* try_reciprocal_trig_power(Expr* f, Expr* x) {
    if (x->type != EXPR_SYMBOL) return NULL;
    Expr* coeff = expr_new_integer(1);
    const Expr* kp = NULL;                       /* the single Power[Sec|Csc[x], n] */
    size_t nf = head_is(f, SYM_Times) ? f->data.function.arg_count : 1;
    for (size_t i = 0; i < nf; i++) {
        Expr* fi = head_is(f, SYM_Times) ? f->data.function.args[i] : f;
        if (!depends_on_var(fi, x)) {
            coeff = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ coeff, expr_copy(fi) }, 2));
        } else if (!kp) kp = fi;
        else { expr_free(coeff); return NULL; }
    }
    if (!kp || !head_is((Expr*)kp, SYM_Power) || kp->data.function.arg_count != 2) {
        expr_free(coeff); return NULL;
    }
    const Expr* k = kp->data.function.args[0];
    const Expr* ne = kp->data.function.args[1];
    if (ne->type != EXPR_INTEGER || ne->data.integer < 3 || k->type != EXPR_FUNCTION
        || k->data.function.arg_count != 1 || !expr_eq(k->data.function.args[0], x)
        || k->data.function.head->type != EXPR_SYMBOL) { expr_free(coeff); return NULL; }
    const char* kh = k->data.function.head->data.symbol.name;
    if (kh != SYM_Sec && kh != SYM_Csc) { expr_free(coeff); return NULL; }
    int64_t n = ne->data.integer;
    bool sec = (kh == SYM_Sec);
    /* Sec^(n-2) Tan / (n-1)   or   -Csc^(n-2) Cot / (n-1) */
    Expr* kpow = expr_new_function(expr_new_symbol(SYM_Power),
                     (Expr*[]){ expr_copy((Expr*)k), expr_new_integer(n - 2) }, 2);
    Expr* partner = expr_new_function(expr_new_symbol(sec ? SYM_Tan : SYM_Cot),
                        (Expr*[]){ expr_copy(x) }, 1);
    Expr* lead = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
                     expr_new_function(expr_new_symbol(SYM_Rational),
                         (Expr*[]){ expr_new_integer(sec ? 1 : -1), expr_new_integer(n - 1) }, 2),
                     kpow, partner }, 3);
    /* (n-2)/(n-1) Integrate[K^(n-2), x] -- must close (it bottoms out at n = 1, 2) */
    Expr* sub = eval_and_free(expr_new_function(expr_new_symbol(SYM_Integrate), (Expr*[]){
                    expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ expr_copy((Expr*)k), expr_new_integer(n - 2) }, 2),
                    expr_copy(x) }, 2));
    if (!sub || result_contains_head(sub, SYM_Integrate)) {
        if (sub) expr_free(sub);
        expr_free(lead); expr_free(coeff); return NULL;
    }
    Expr* rest = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
                     expr_new_function(expr_new_symbol(SYM_Rational),
                         (Expr*[]){ expr_new_integer(n - 2), expr_new_integer(n - 1) }, 2),
                     sub }, 2);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
               coeff, expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){ lead, rest }, 2) }, 2));
}

static Expr* try_weierstrass(Expr* f, Expr* x) {
    Expr* r = try_reciprocal_trig_power(f, x);
    if (r) return r;
    return integrate_jeffrey_try(f, x);
}

/* Stage 1f-: linear-argument substitution, run just before Weierstrass.  A
 * rational trig/hyperbolic integrand whose kernels all share one non-trivial
 * linear argument a*x+b is reduced to the bare-argument integral via u = a*x+b,
 * which the recursive cascade (Weierstrass itself, cleanly) then closes:
 * Integrate[Sec[3x]^2, x] = Tan[3x]/3 rather than the degree-12 Tan[x/2] mess the
 * always-x/2 Weierstrass substitution would produce.  Declines on everything
 * else, so bare-argument and polynomial-trig integrands are untouched. */
static Expr* try_linarg(Expr* f, Expr* x) {
    return integrate_linarg_try(f, x);
}

/* Does the integrand contain a circular-trig kernel anywhere? */
static bool integrand_has_trig(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        const char* n = h->data.symbol.name;
        if (n == SYM_Sin || n == SYM_Cos || n == SYM_Tan ||
            n == SYM_Cot || n == SYM_Sec || n == SYM_Csc) return true;
    }
    if (integrand_has_trig(h)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (integrand_has_trig(e->data.function.args[i])) return true;
    return false;
}

/* Every circular-trig kernel in `e` has an argument that depends on the
 * integration variable `x` (i.e. `e` is a trig expression in x, not one carrying
 * a foreign trig kernel such as Sec[y] treated as a parameter).  This keeps the
 * Expand[TrigReduce[·]] linearization off multi-variable trig blobs, whose
 * Expand blows up (the Lie first-integral integrands). */
static bool trig_all_in_x(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return true;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && e->data.function.arg_count == 1) {
        const char* n = h->data.symbol.name;
        if ((n == SYM_Sin || n == SYM_Cos || n == SYM_Tan ||
             n == SYM_Cot || n == SYM_Sec || n == SYM_Csc) &&
            !depends_on_var(e->data.function.args[0], x))
            return false;
    }
    if (!trig_all_in_x(h, x)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!trig_all_in_x(e->data.function.args[i], x)) return false;
    return true;
}

/* Node count, capped: returns >= limit as soon as the budget is exceeded so a
 * large integrand is rejected cheaply. */
static long node_count_capped(const Expr* e, long limit) {
    if (!e || limit <= 0) return 0;
    long c = 1;
    if (e->type == EXPR_FUNCTION) {
        c += node_count_capped(e->data.function.head, limit - c);
        for (size_t i = 0; i < e->data.function.arg_count && c < limit; i++)
            c += node_count_capped(e->data.function.args[i], limit - c);
    }
    return c;
}

/* Build head[arg] and evaluate it (borrowing arg; the caller keeps ownership). */
static Expr* apply1_eval(const char* head, const Expr* arg) {
    Expr* call = expr_new_function(expr_new_symbol(head),
                                   (Expr*[]){ expr_copy((Expr*)arg) }, 1);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* Stage 1f0: linearity.  Split a Plus integrand and integrate term-by-term,
 * committing ONLY if every term closes to an elementary antiderivative (no
 * residual Integrate[...]); otherwise decline so the whole-integrand cascade
 * (Weierstrass, Risch, ...) still runs and can capture a sum whose pieces are
 * individually non-elementary.  Placed after the rational integrator (so
 * rational sums keep their combined partial-fraction form) and before
 * Weierstrass, which otherwise routes a trig polynomial such as
 * Sin[2 x] + Sin[3 x] through the tan-half-angle substitution and produces a
 * divergent rational form — or loops (the y' + Tan[x] y == q integrating-factor
 * hang, and the constant-coefficient forcing integrals).
 *
 * A trig integrand that is NOT already a Plus (a product/power such as
 * Cos[x]^2 or Cos[x](3 Sin[x] + K)) is first linearised by Expand[TrigReduce[·]]
 * into a sum of Sin[k x]/Cos[k x] monomials, which the same term-by-term split
 * then integrates cleanly.  Recursion handles nested terms; TrigReduce strictly
 * lowers trig degree, so the descent terminates. */
/* True when `f` is a product carrying a sum as a direct factor — the shape
 * `c (g + h)` that linearity can distribute into `c g + c h`. */
static bool times_has_plus_factor(const Expr* f) {
    if (!head_is((Expr*)f, SYM_Times)) return false;
    for (size_t i = 0; i < f->data.function.arg_count; i++)
        if (head_is(f->data.function.args[i], SYM_Plus)) return true;
    return false;
}

/* Exponent[poly, x] as a long, or -1 when it is not a non-negative machine
 * integer (poly not polynomial in x, or Exponent left unevaluated).  Borrows
 * both arguments; the caller keeps ownership. */
static long poly_degree_in(const Expr* poly, const Expr* x) {
    Expr* call = expr_new_function(expr_new_symbol("Exponent"),
                                   (Expr*[]){ expr_copy((Expr*)poly),
                                              expr_copy((Expr*)x) }, 2);
    Expr* r = evaluate(call);
    expr_free(call);
    long d = -1;
    if (r && r->type == EXPR_INTEGER && r->data.integer >= 0)
        d = (long)r->data.integer;
    if (r) expr_free(r);
    return d;
}

/* True if `e` contains a radical Power[base, r], r a non-integer rational
 * (head Rational), whose base is a polynomial of degree >= 3 in x -- the
 * pseudo-elliptic F(x)/R(x)^p shape (Goursat's domain, p in
 * {1/2,1/3,2/3,1/4,3/4}).  Distributing such an integrand over a Plus shreds
 * the shared radical into pieces that are individually non-elementary
 * (elliptic), each of which then grinds through the whole-integrand cascade
 * -- ParallelMixedTower especially -- before linearity's all-terms-close
 * guard rejects the split.  Radicals of a degree<=2 radicand (Sqrt[quadratic]
 * -> ArcSin/ArcTanh/Log) yield elementary pieces and split cleanly, so they
 * are deliberately NOT flagged. */
static bool has_pseudoelliptic_radical(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_is((Expr*)e, SYM_Power) && e->data.function.arg_count == 2
        && head_is(e->data.function.args[1], SYM_Rational)
        && poly_degree_in(e->data.function.args[0], x) >= 3)
        return true;
    if (has_pseudoelliptic_radical(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_pseudoelliptic_radical(e->data.function.args[i], x)) return true;
    return false;
}

/* A circular-trig kernel sits in a DENOMINATOR of `e`: a Sec/Csc/Cot/Tan kernel
 * (each an implicit reciprocal of Cos/Sin), or a Sin/Cos kernel raised to a
 * negative numeric power.  That is exactly the rational-trig domain the
 * Weierstrass (Jeffrey-Rich) stage closes directly and fast, and TrigReduce
 * CANNOT clear a trig denominator -- so the Expand[TrigReduce] linearity branch
 * below must not pre-expand such an integrand: its term-by-term recursion
 * re-enters the same rational-trig form and spins (the cold Integrate[Sec[x]
 * Tan[x]^2, x] / `y' + Tan[x] y == q` 12 s hang, closed in 0.03 s by Weierstrass
 * once linearity declines).  Only ever consulted when trig_all_in_x already
 * holds, so every kernel here is a kernel in x. */
static bool trig_in_denominator(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && e->data.function.arg_count == 1) {
        const char* n = h->data.symbol.name;
        if (n == SYM_Sec || n == SYM_Csc || n == SYM_Cot || n == SYM_Tan) return true;
    }
    if (h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Power
        && e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* ex   = e->data.function.args[1];
        if (base->type == EXPR_FUNCTION && base->data.function.head->type == EXPR_SYMBOL) {
            const char* bn = base->data.function.head->data.symbol.name;
            if ((bn == SYM_Sin || bn == SYM_Cos)
                && ((ex->type == EXPR_INTEGER && ex->data.integer < 0)
                    || (ex->type == EXPR_REAL && ex->data.real < 0.0))) return true;
        }
    }
    if (trig_in_denominator(h)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (trig_in_denominator(e->data.function.args[i])) return true;
    return false;
}

/* `e` is a RATIONAL FUNCTION of circular-trig kernels of the integration variable
 * and numeric constants -- built only from numbers, the constants Pi/E/EulerGamma/…,
 * the trig kernels Sin/Cos/Tan/Cot/Sec/Csc[·], and the operators +, *, and integer
 * powers.  This is exactly the Jeffrey–Rich (Weierstrass) domain.  A trig kernel is
 * a leaf here: its argument's x-dependence is the caller's `trig_all_in_x`
 * business, so we do NOT recurse into it (that keeps Sin[2x], Cos[x/2] rational).
 * The point is to EXCLUDE the integrands that merely CONTAIN a trig kernel inside a
 * larger structure the Weierstrass stage cannot take but Linearity can -- an Exp /
 * Log / bare-parameter / bare-x body such as an integrating-factor expression
 * `Y^2 E^(Log[1 + 2 Y Tan[x]] + …)` (§2.2.35-3461): its `Power[Y, 2]` base is a
 * bare non-constant symbol, so this returns false and the integrand stays with
 * Linearity, which closes it (Weierstrass does not).  A bare `x` likewise fails
 * (e.g. `x Sin[x]` is not a rational-trig integrand). */
static bool is_rational_trig_in_x(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT) return true;
    if (e->type == EXPR_SYMBOL) {
        const char* n = e->data.symbol.name;         /* numeric constants only */
        return n == SYM_Pi || n == SYM_E || n == SYM_EulerGamma
            || n == SYM_GoldenRatio || n == SYM_Catalan || n == SYM_Degree || n == SYM_I;
    }
    if (e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL) return false;
    const char* n = h->data.symbol.name;
    if (e->data.function.arg_count == 1
        && (n == SYM_Sin || n == SYM_Cos || n == SYM_Tan ||
            n == SYM_Cot || n == SYM_Sec || n == SYM_Csc)) return true;   /* trig leaf */
    if (n == SYM_Rational) return true;                                   /* numeric coeff */
    if (n == SYM_Plus || n == SYM_Times) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!is_rational_trig_in_x(e->data.function.args[i])) return false;
        return true;
    }
    if (n == SYM_Power && e->data.function.arg_count == 2) {
        const Expr* b  = e->data.function.args[0];
        const Expr* ex = e->data.function.args[1];
        return ex->type == EXPR_INTEGER && is_rational_trig_in_x(b);      /* integer power of a rational-trig */
    }
    return false;
}

static Expr* try_linearity(Expr* f, Expr* x) {
    Expr* target = NULL;                 /* a Plus we own (or borrow from f) */
    Expr* owned  = NULL;                 /* non-NULL when we built `target` */
    if (head_is(f, SYM_Plus) && f->data.function.arg_count >= 2) {
        target = f;
    } else if (integrand_has_trig(f) && trig_all_in_x(f, x)
               && !(is_rational_trig_in_x(f) && trig_in_denominator(f))
               && node_count_capped(f, 256) < 256) {
        Expr* tr = apply1_eval("TrigReduce", f);
        Expr* g  = tr ? apply1_eval("Expand", tr) : NULL;
        if (tr) expr_free(tr);
        if (g && head_is(g, SYM_Plus) && g->data.function.arg_count >= 2) {
            target = g; owned = g;
        } else { if (g) expr_free(g); return NULL; }
    } else if (times_has_plus_factor(f) && !has_pseudoelliptic_radical(f, x)
               && node_count_capped(f, 256) < 256) {
        /* Distribute a product over a sum factor: Integrate is linear, so
         * c (g + h) -> c g + c h, integrated term-by-term.  This is what combines
         * an exponential product such as E^{-9x}(a E^{2x} - b E^{2x}) into the
         * sum a E^{-7x} - b E^{-7x}; the whole-integrand substitution path
         * (DerivativeDivides) otherwise routes the un-distributed product through
         * a fractional-power base substitution (u = E^{2x}), which is slow and,
         * for funcapp/symbolic coefficients, returns a branch-WRONG antiderivative
         * (a spurious (-1)^(1/d) factor).  try_rational has already consumed any
         * polynomial/rational Times-Plus, so this fires only on transcendental /
         * radical products.  The commit-only-if-every-term-closes guard below
         * leaves a product whose distributed pieces are individually
         * non-elementary on the whole-integrand cascade, so no sum that is
         * elementary only as a whole is lost.  The has_pseudoelliptic_radical
         * gate above additionally skips a product sharing a Sqrt/odd-root of a
         * degree>=3 polynomial (Goursat's pseudo-elliptic F/R^p domain):
         * distributing there manufactures individually non-elementary elliptic
         * pieces that each grind through PMT before the guard rejects the
         * split, so the integrand is left intact for the Goursat stage. */
        Expr* g = apply1_eval("Expand", f);
        if (g && head_is(g, SYM_Plus) && g->data.function.arg_count >= 2) {
            target = g; owned = g;
        } else { if (g) expr_free(g); return NULL; }
    } else {
        return NULL;
    }

    size_t n = target->data.function.arg_count;
    Expr** parts = malloc(n * sizeof(Expr*));
    size_t got = 0;
    bool ok = true;
    for (size_t i = 0; i < n; i++) {
        Expr* gi = call_stage("Integrate", target->data.function.args[i], x);
        if (!gi || result_contains_head(gi, "Integrate")) {
            if (gi) expr_free(gi);
            ok = false;
            break;
        }
        parts[got++] = gi;
    }
    if (!ok) {
        for (size_t i = 0; i < got; i++) expr_free(parts[i]);
        free(parts);
        if (owned) expr_free(owned);
        return NULL;
    }
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), parts, n);
    free(parts);
    if (owned) expr_free(owned);
    Expr* r = evaluate(sum);   /* evaluate BORROWS sum */
    expr_free(sum);
    return r;
}

/* Stage 1e3: Goursat's pseudo-elliptic algorithm and its cube-/fourth-root
 * generalisations.  Recognises F(x)/R(x)^p (p in {1/2,1/3,2/3,1/4,3/4}); when
 * a Mobius automorphism of the roots of R makes the integrand pseudo-elliptic,
 * descends to genus-0 curves and integrates in closed form.  Deterministic and
 * correct by construction, so it runs ahead of the Eliminate/Solve search and
 * Risch-Norman.  Non-elementary (obstructed) integrands return NULL. */
static Expr* try_goursat(Expr* f, Expr* x) {
    return integrate_goursat_try(f, x);
}

/* K Sin[a x^2+b x+c] / K Cos[...] -> FresnelS/FresnelC (complete the square). */
static Expr* try_fresnel(Expr* f, Expr* x) {
    return integrate_fresnel_try(f, x);
}

/* K x^p E^(a x^m), p SYMBOLIC -> the incomplete Gamma[(p+1)/m, -a x^m].  Gated to
 * a symbolic exponent so a numeric one keeps its elementary / Erf answer from the
 * later stages; see integrate_gammapower.c for why it is worth a stage of its own
 * (the general stages SEARCH for 12.9 s before declining on this shape). */
static Expr* try_gammapower(Expr* f, Expr* x) {
    return integrate_gammapower_try(f, x);
}

/* c Log[g] K, K carrying a trig/hyperbolic kernel of x -> one integration by
 * parts (Log[g] INT K - INT K'-log-weight), closing the Log*trig family to a
 * Log*trig + SinIntegral/CosIntegral form instead of the ~3.8-12 s
 * ParallelMixedSpecial search; see integrate_logbyparts.c.  Diff-back verified,
 * so it declines (never emits a wrong form) on a mis-match. */
static Expr* try_logbyparts(Expr* f, Expr* x) {
    return integrate_logbyparts_try(f, x);
}

/* Stage: polylogarithm ladders (src/internal/PolyLogLadder.m).
 *
 * Two exact constructions, loaded lazily on first use:
 *   - R(x) Log[a x + b]^k (rational R, poles anywhere -- algebraic or complex
 *     roots of any multiplicity -- and any log power, one linear log kernel per
 *     term): partial fractions in the kernel variable W = a x + b closed by the
 *     polylogarithm ladder INT Log[W]^k/(W - c) = Log[W]^k Log[1 - W/c] +
 *     k Sum (-1)^t (k-1)!/(k-1-t)! Log[W]^(k-1-t) PolyLog[t+2, W/c];
 *   - a polynomial in x times a rational function of trig/hyperbolic functions
 *     of linear arguments: TrigToExp in one exponential kernel T = E^(base x),
 *     partial fractions in T, and the exponential ladder for x^n/(T - rho)
 *     (rho-derivatives of it for repeated roots).
 * Both used to be missing: Integrate[Log[x]/(1 + x^2), x],
 * Integrate[Log[x]^2/(1 + x), x] and Integrate[x^2 Csc[x], x] all came back
 * unevaluated (Mathar corpus, the largest failure class).  The package accepts
 * a result only after a differentiation check; anything outside the class
 * returns $Failed and the cascade continues. */
static unsigned pll_load_gen = 0;   /* symtab_generation() of the attempt */
static bool pll_load_succeeded = false;

static Expr* try_polylogladder(Expr* f, Expr* x) {
    if (pll_load_gen != symtab_generation()) {
        pll_load_gen = symtab_generation();
        mth_msg_suppress_push();
        pll_load_succeeded = mathilda_load_module("PolyLogLadder.m") != 0;
        mth_msg_suppress_pop();
        if (!pll_load_succeeded)
            mth_message("Integrate`PolyLogLadder", "nofile",
                "cannot locate src/internal/PolyLogLadder.m on disk.");
    }
    if (!pll_load_succeeded) return NULL;
    mth_msg_suppress_push();
    Expr* result = call_stage("PolyLogLadder`Integrate", f, x);
    mth_msg_suppress_pop();
    if (!result) return NULL;
    /* Any atom is a decline: the package answers $Failed, and a symbol read in
     * while loading a module can carry a context prefix (`...`$Failed), so the
     * name is not compared at all -- an antiderivative is never a bare symbol
     * here (the stage only accepts integrands that depend on x). */
    if (result_is_unresolved(result, "PolyLogLadder`Integrate") ||
        result->type != EXPR_FUNCTION ||
        result_contains_head(result, "Integrate") ||
        result_contains_head(result, "Return") || result_contains_head(result, "Throw")) {
        expr_free(result);
        return NULL;
    }
    return result;
}

/* Stage: holonomic ansatz (src/internal/Holonomic.m).
 *
 * x^alpha L(x) times monomials in Bessel J/Y/I/K, Airy and Sin/Cos kernels of
 * linear argument.  Each kernel family spans a 2-dimensional space closed
 * under d/dx with Laurent-polynomial coefficients, so a product of kernels
 * lives in a finite-dimensional differential module; the package seeks an
 * antiderivative in that module with Laurent-polynomial coefficients by
 * solving one linear system (Integrate[BesselJ[1, x]^4/x, x],
 * Integrate[BesselJ[nu, x]^2/x^4, x], Integrate[x^(1 + mu + nu) BesselJ[mu, x]
 * BesselY[nu, x], x] used to come back unevaluated).  Exact by construction
 * and re-checked symbolically; anything else is $Failed and the cascade
 * continues.  Gated on a Bessel or Airy head, so it costs nothing elsewhere. */
static unsigned hol_load_gen = 0;   /* symtab_generation() of the attempt */
static bool hol_load_succeeded = false;

static bool hol_has_kernel(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        const char* n = h->data.symbol.name;
        if (strcmp(n, "BesselJ") == 0 || strcmp(n, "BesselY") == 0 ||
            strcmp(n, "BesselI") == 0 || strcmp(n, "BesselK") == 0 ||
            strcmp(n, "AiryAi") == 0 || strcmp(n, "AiryBi") == 0 ||
            strcmp(n, "AiryAiPrime") == 0 || strcmp(n, "AiryBiPrime") == 0)
            return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (hol_has_kernel(e->data.function.args[i])) return true;
    return false;
}

static Expr* try_holonomic(Expr* f, Expr* x) {
    if (!hol_has_kernel(f)) return NULL;
    if (hol_load_gen != symtab_generation()) {
        hol_load_gen = symtab_generation();
        mth_msg_suppress_push();
        hol_load_succeeded = mathilda_load_module("Holonomic.m") != 0;
        mth_msg_suppress_pop();
        if (!hol_load_succeeded)
            mth_message("Integrate`Holonomic", "nofile",
                "cannot locate src/internal/Holonomic.m on disk.");
    }
    if (!hol_load_succeeded) return NULL;
    mth_msg_suppress_push();
    Expr* result = call_stage("Holonomic`Integrate", f, x);
    mth_msg_suppress_pop();
    if (!result) return NULL;
    /* any atom is a decline (the package's $Failed), as for PolyLogLadder */
    if (result_is_unresolved(result, "Holonomic`Integrate") ||
        result->type != EXPR_FUNCTION ||
        result_contains_head(result, "Integrate") ||
        result_contains_head(result, "Return") || result_contains_head(result, "Throw")) {
        expr_free(result);
        return NULL;
    }
    return result;
}

/* Stage: Bessel hypergeometric closed forms (src/internal/BesselHyper.m).
 *
 * c x^m Z1(b x^p) [Z2(b x^p)], Z in {J, Y, I, K}: the product's power series
 * (the Cauchy product of two Bessel series sums by Chu-Vandermonde to a 2F3)
 * integrated term by term into a HypergeometricPFQ, or -- for integer and
 * symbolic orders of Y / K, where the series route divides by Sin[nu Pi] --
 * its Mellin-Barnes MeijerG integrated under the contour.  Integer order
 * differences reduce to diagonal products by the recurrences and the
 * Wronskian.  Runs after Holonomic so an antiderivative that closes in Bessel
 * functions keeps that form; this covers the ones that do not
 * (Integrate[BesselJ[1/3, x]^2, x], Integrate[BesselK[0, x], x], ...).
 * Exact by construction; anything outside the class is $Failed. */
static unsigned bh_load_gen = 0;   /* symtab_generation() of the attempt */
static bool bh_load_succeeded = false;

static Expr* try_besselhyper(Expr* f, Expr* x) {
    if (!hol_has_kernel(f)) return NULL;
    if (bh_load_gen != symtab_generation()) {
        bh_load_gen = symtab_generation();
        mth_msg_suppress_push();
        bh_load_succeeded = mathilda_load_module("BesselHyper.m") != 0;
        mth_msg_suppress_pop();
        if (!bh_load_succeeded)
            mth_message("Integrate`BesselHypergeometric", "nofile",
                "cannot locate src/internal/BesselHyper.m on disk.");
    }
    if (!bh_load_succeeded) return NULL;
    mth_msg_suppress_push();
    Expr* result = call_stage("BesselHyper`Integrate", f, x);
    mth_msg_suppress_pop();
    if (!result) return NULL;
    if (result_is_unresolved(result, "BesselHyper`Integrate") ||
        result->type != EXPR_FUNCTION ||
        result_contains_head(result, "Integrate") ||
        result_contains_head(result, "Return") || result_contains_head(result, "Throw")) {
        expr_free(result);
        return NULL;
    }
    return result;
}

/* Stage 2: recursive transcendental Risch integrator.
 * A decision procedure over a differential transcendental tower, with
 * rational / logarithmic / exponential / special-function cases, each
 * correct by construction (it fires only behind an exact certificate, so
 * it emits no wrong closed forms and needs no differentiation check).  The
 * former parallel-Risch (pmint) heuristics Integrate`RischNorman and
 * Integrate`RischNormanBlake were removed in v0.163: ParallelMixedTower (a
 * later stage) subsumes both -- the transcendental case and the simple radical. */
static Expr* try_rischtranscendental(Expr* f, Expr* x) {
    Expr* result = call_stage("Integrate`RischTranscendental", f, x);
    if (!result) return NULL;
    if (result_is_unresolved(result, "Integrate`RischTranscendental")) {
        expr_free(result);
        return NULL;
    }
    return result;
}

/* Stage 3: CRC integral table.  Loaded from disk on first invocation
 * (lazy) so unaffected sessions pay nothing for it. */
#define MAX_CRC_DEPTH 256
static int crc_depth = 0;
static unsigned crc_load_gen = 0;   /* symtab_generation() of the attempt */
static bool crc_load_succeeded  = false;

/* Load the CRC integral table lazily. Path resolution (mathilda_resolve_
 * internal) is CWD-independent and matches what LoadModule uses, so this works
 * from a relocated/installed binary as well as from the repo root or tests/. */
static void crc_lazy_load(void) {
    if (crc_load_gen == symtab_generation()) return;
    crc_load_gen = symtab_generation();

    char path[2048];
    if (!mathilda_resolve_internal("CRCMathTablesIntegrals.m",
                                   path, sizeof(path))) {
        mth_message("Integrate`CRCTable", "nofile",
            "cannot locate src/internal/CRCMathTablesIntegrals.m on disk.");
        return;
    }

    int opened = 0;
    Expr* res = mathilda_run_file(path, &opened);
    bool failed = res && res->type == EXPR_SYMBOL
                      && strcmp(res->data.symbol.name, "$Failed") == 0;
    if (res) expr_free(res);
    if (opened && !failed) crc_load_succeeded = true;
}

static Expr* try_crctable(Expr* f, Expr* x) {
    crc_lazy_load();
    if (!crc_load_succeeded) return NULL;

    if (crc_depth >= MAX_CRC_DEPTH) {
        mth_message("Integrate`CRCTable", "depth",
            "rule recursion exceeded %d levels; "
            "the table has a divergent rule for this integrand.",
            MAX_CRC_DEPTH);
        return NULL;
    }
    crc_depth++;
    Expr* result = call_stage("Integrate`CRCTable", f, x);
    crc_depth--;
    if (!result) return NULL;

    /* "Failed" sentinels: either the public head is unresolved (no
     * rule matched the CRCTable wrapper) or the internal IntegrateTable
     * head leaked through (table lookup found no matching formula).
     * The IntegrateTable scan is deep — a partial-match rule like
     * `IntegrateTable[a_ f_, x_] /; FreeQ[a, x] := a IntegrateTable[f, x]`
     * can leave `IntegrateTable[...]` nested under Times/Plus when the
     * inner call has no matching formula; that must also count as
     * unresolved so Integrate falls back to its own unevaluated form. */
    if (result_is_unresolved(result, "Integrate`CRCTable") ||
        result_contains_head(result, "IntegrateTable")) {
        expr_free(result);
        return NULL;
    }
    return result;
}

/* Stage: ParallelMixedTower.  A parallel (Risch-Norman) integrator over a
 * simple radical in a mixed transcendental tower (S. Blake, Parallel
 * Integration over Simple Radical Extensions II), implemented as a Wolfram-
 * language package in src/internal/mixed/ParallelMixed.m and loaded lazily on
 * first use exactly like the CRC table above -- sessions that never reach this
 * cascade stage pay nothing for it. */
static unsigned pmt_load_gen = 0;   /* symtab_generation() of the attempt */
static bool pmt_load_succeeded  = false;
/* Async-jump defer (crash-safety inside TimeConstrained) + re-entry guard; see
 * core.c.  tc_async_region_active is non-zero while a heavy method is already
 * running, so the package's internal Integrate[] on a sub-problem declines here
 * rather than re-entering the method and looping. */
void tc_async_defer_push(void);
void tc_async_defer_pop(void);
int  tc_async_region_active(void);

static void pmt_lazy_load(void) {
    if (pmt_load_gen == symtab_generation()) return;
    pmt_load_gen = symtab_generation();
    /* mathilda_load_module rather than a bare mathilda_run_file: it carries a
     * load-once memo keyed on the relpath, which ParallelMixedSpecial.m also
     * loads Part II through.  Without the shared memo a session that reaches
     * both methods parsed and evaluated this 217 KB file TWICE (~1 s each),
     * re-running EndPackage[] and the logrewrite load -- harmless semantically
     * (the definitions are just re-Set) but a visible cold-start cost and a
     * source of duplicated Print output. */
    if (mathilda_load_module("mixed/ParallelMixed.m")) {
        pmt_load_succeeded = true;
        return;
    }
    mth_message("Integrate`ParallelMixedTower", "nofile",
        "cannot locate src/internal/mixed/ParallelMixed.m on disk.");
}

/* True iff `r` is a {"not elementary", ...} list -- the ParallelMixedTower
 * package's RIGOROUS non-elementarity certificate (paper Algorithm 4, the three
 * guarded NotElementary exits: a non-constant residue exactly decided in the
 * residue field, a proved non-torsion residue divisor, or a verified
 * residue-free holomorphic remainder insolvable under the exact proved bounds).
 * A {"failed", ...} list, by contrast, is an inconclusive give-up that proves
 * nothing -- the classical parallel method's failure -- and must NOT warn.  The
 * package draws that boundary internally; here we only read the leading tag. */
static bool pmt_is_nonelementary_certificate(const Expr* r) {
    if (!r || r->type != EXPR_FUNCTION) return false;
    Expr* h = r->data.function.head;
    if (!h || h->type != EXPR_SYMBOL || h->data.symbol.name != SYM_List) return false;
    if (r->data.function.arg_count < 1) return false;
    const Expr* a0 = r->data.function.args[0];
    return a0 && a0->type == EXPR_STRING
        && strcmp(a0->data.string, "not elementary") == 0;
}

/* The Integrate`ParallelMixedTower method symbol.  A thin C builtin so the
 * qualified-symbol surface Integrate`ParallelMixedTower[f, x] is callable from a
 * cold session without paying the ~1s package load until it is actually reached:
 * the builtin loads the .m worker lazily and delegates to it, returning the
 * worker's value unchanged -- the antiderivative, $Failed, or a List
 * certificate/failure.
 *
 * The delegation runs inside an async-jump-deferred region
 * (tc_async_defer_push/pop): the package is malloc-heavy, and under an outer
 * TimeConstrained the SIGPROF handler could otherwise siglongjmp out of a
 * regular malloc (an Expr node, a FLINT structure) and crash.  Deferred, the
 * deadline is still enforced at safe points (tc_check_deadline every rewrite
 * step, tc_alloc_safepoint after each GMP alloc), so a runaway is still aborted
 * -- just cleanly.  tc_async_region_active doubles as the single re-entry guard
 * through which both the direct call and the cascade stage pass.  Both the defer
 * count and the message-suppression depth are save/restored by tc_run_guarded,
 * so a timeout unwinding through here leaves neither stuck.
 *
 * Time budgets.  The package carries three wall-clock backstops
 * ($ParallelMixedTimeBudget 45 s for the solve, $ParallelMixedVerifyBudget 15 s
 * for the numeric verification, $RewriteBudget 30 s for the real-form rewrite of
 * the logarithms).  They exist for the Automatic cascade, which runs this stage
 * on integrands nobody chose for it and must get out of the way.  An explicit
 * request -- Method -> "ParallelMixedTower" or the qualified symbol -- is a
 * deliberate request for the full decision, so all three are Blocked to
 * Infinity: the method runs until it answers, declines on its own terms, or the
 * caller's TimeConstrained fires.  (An exhausted verify budget ACCEPTED the
 * antiderivative unverified as the principal branch; unbounded, every explicit
 * answer is verified to completion.)  The cascade arms pmt_in_cascade, consumed
 * here exactly like pms_in_cascade (see its comment), so an unwind out of a
 * cascade call cannot hand the backstops to the next explicit call or vice versa. */
static int pmt_in_cascade = 0;

static Expr* builtin_integrate_pmt(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    if (tc_async_region_active()) return NULL;   /* no self re-entry / nested heavy method */
    pmt_lazy_load();
    if (!pmt_load_succeeded) return NULL;

    int cascade = pmt_in_cascade;
    pmt_in_cascade = 0;            /* consume; see pms_in_cascade */

    tc_async_defer_push();
    /* Suppress the package's own diagnostics (BuildTower::function on an
     * unsupported generator, etc.): as a method probe it must not spew to the
     * user.  The fired-counter that the package's internal Check reads still
     * advances while suppressed, so its failure detection is unaffected. */
    mth_msg_suppress_push();
    bool overflow_before = eval_overflow_fired();
    Expr* result;
    if (cascade) {
        result = call_stage("ParallelMixed`ParallelIntegrateMixed",
                            res->data.function.args[0], res->data.function.args[1]);
    } else {
        /* Block[{$ParallelMixedTimeBudget = Infinity, $ParallelMixedVerifyBudget
         * = Infinity, $RewriteBudget = Infinity}, worker[f, x]] -- Block so the
         * backstops are restored on every exit path, including a TimeConstrained
         * unwind; fully qualified, since a misspelt context would bind a
         * different symbol and silently leave the backstop in force. */
        static const char* const budgets[] = {
            "ParallelMixed`$ParallelMixedTimeBudget",
            "ParallelMixed`$ParallelMixedVerifyBudget",
            "ParallelMixed`$RewriteBudget" };
        Expr* binds[3];
        for (int i = 0; i < 3; i++)
            binds[i] = expr_new_function(expr_new_symbol(SYM_Set),
                (Expr*[]){ expr_new_symbol(budgets[i]), expr_new_symbol(SYM_Infinity) }, 2);
        Expr* call = expr_new_function(
            expr_new_symbol("ParallelMixed`ParallelIntegrateMixed"),
            (Expr*[]){ expr_copy(res->data.function.args[0]),
                       expr_copy(res->data.function.args[1]) }, 2);
        Expr* locals = expr_new_function(expr_new_symbol(SYM_List), binds, 3);
        Expr* blk = expr_new_function(expr_new_symbol(SYM_Block),
                                      (Expr*[]){ locals, call }, 2);
        result = evaluate(blk);
        expr_free(blk);
    }
    mth_msg_suppress_pop();
    tc_async_defer_pop();
    /* A $RecursionLimit overflow inside the package does not abort: it leaves
     * every enclosing loop in a one-rewrite degraded mode, in which the
     * package's own decline branches can be skipped.  Nothing produced in that
     * state is evidence of anything -- neither an antiderivative nor a
     * certificate -- so it becomes a failure tuple (a List: the cascade and
     * the Method stage treat it as a decline). */
    if (!overflow_before && eval_overflow_fired()) {
        if (result) expr_free(result);
        result = expr_new_function(expr_new_symbol(SYM_List),
            (Expr*[]){ expr_new_string("failed"),
                       expr_new_string("internal error: $RecursionLimit exceeded") }, 2);
    }
    /* A {"not elementary", ...} certificate is the package's PROOF that no
     * elementary antiderivative exists -- report it Mathematica-style, exactly
     * as RischTranscendental does for its own field decision.  Emitted after the
     * message-mute region (the diagnostic is about the user's integrand, not a
     * package probe) and after the async-defer pop.  The result itself is
     * returned unchanged: the qualified-symbol surface still exposes the raw
     * certificate, and the cascade/Method stage (try_parallelmixedtower) still
     * discards the List to leave Integrate[...] unevaluated. */
    if (pmt_is_nonelementary_certificate(result))
        integrate_announce_nonelementary(res->data.function.args[0],
                                         res->data.function.args[1]);
    return result;   /* expression, $Failed, or a List certificate/failure */
}

/* True if e is (structurally) a rational function -- built only from numbers,
 * symbols, Plus, Times, Rational/Complex, and Power with an integer exponent.  A
 * non-integer power (a radical) or any other head (Log, Exp, Sin, ...) makes it
 * non-rational.  ParallelMixedTower's value is the radical/transcendental-tower
 * case; a pure rational function is BronsteinRational's job (an earlier stage),
 * so the cascade must NOT spend this stage's search on one -- notably a
 * PARAMETRIC rational that BronsteinRational declines, which would otherwise
 * burn the full iPIM attempt before failing anyway. */
static bool pmt_is_rational_structure(const Expr* e) {
    if (!e) return true;
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_BIGINT: case EXPR_REAL: case EXPR_SYMBOL:
            return true;
        case EXPR_FUNCTION: {
            if (e->data.function.head->type != EXPR_SYMBOL) return false;
            const char* h = e->data.function.head->data.symbol.name;
            size_t n = e->data.function.arg_count;
            if (h == SYM_Plus || h == SYM_Times ||
                h == SYM_Rational || h == SYM_Complex) {
                for (size_t i = 0; i < n; i++)
                    if (!pmt_is_rational_structure(e->data.function.args[i]))
                        return false;
                return true;
            }
            if (h == SYM_Power && n == 2) {
                if (e->data.function.args[1]->type != EXPR_INTEGER) return false;
                return pmt_is_rational_structure(e->data.function.args[0]);
            }
            return false;   /* Log, Exp, Sin, Sqrt(=Power^rational), ... */
        }
        default: return false;   /* String, NDArray, ... */
    }
}

/* True iff `e` contains the symbol `x` anywhere (structural). */
static bool pmt_depends(const Expr* e, const Expr* x) {
    if (!e) return false;
    if (expr_eq((Expr*)e, (Expr*)x)) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (pmt_depends(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (pmt_depends(e->data.function.args[i], x)) return true;
    return false;
}

/* ParallelMixedTower decides ELEMENTARY integrands: a tower of exp/log over the
 * rational functions, with trig/hyperbolic (and their inverses) and radicals
 * rewritten into it.  A kernel outside that class applied to an x-dependent
 * argument -- LegendreP[n, m, Cos[x]], an unknown g[Cos[x]] -- is not modelled:
 * the package used to treat it as a CONSTANT after substituting a tower variable
 * for its argument, returning a confident wrong antiderivative with the tower
 * variable leaked into it (Mathar M13537: (2 g[r25])/(n + n Tan[...]^2) for
 * Integrate[g[Cos[t]] Sin[n (z - t)], t]).  Such an integrand is outside the
 * stage's class, so it declines. */
static bool pmt_has_foreign_kernel(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (!pmt_depends(e, x)) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return true;
    static const char* const ELEM[] = {
        "Plus", "Times", "Power", "Sqrt", "Exp", "Log", "Rational", "Complex",
        "Sin", "Cos", "Tan", "Cot", "Sec", "Csc",
        "Sinh", "Cosh", "Tanh", "Coth", "Sech", "Csch",
        "ArcSin", "ArcCos", "ArcTan", "ArcCot", "ArcSec", "ArcCsc",
        "ArcSinh", "ArcCosh", "ArcTanh", "ArcCoth", "ArcSech", "ArcCsch", NULL };
    const char* h = e->data.function.head->data.symbol.name;
    bool ok = false;
    for (int i = 0; ELEM[i]; i++) if (strcmp(h, ELEM[i]) == 0) { ok = true; break; }
    if (!ok) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (pmt_has_foreign_kernel(e->data.function.args[i], x)) return true;
    return false;
}

/* True iff `r` mentions a user-level symbol that occurs in neither `f` nor `x`
 * -- a variable leaked out of an internal substitution.  System symbols (Pi, E,
 * I, builtin heads) are exempt, and so are the bound variables of a Function
 * (RootSum[Function[z, ...], ...]); a Module local (Temporary) is NOT exempt
 * anywhere else, since an attribute alone once let the unassigned
 * ParallelMixed`Private`r$NNN through as the "antiderivative". */
#define PMT_MAX_BOUND 32
static bool pmt_leaks_symbol_b(const Expr* r, const Expr* f, const Expr* x,
                               const Expr** bound, int nb) {
    if (!r) return false;
    if (r->type == EXPR_SYMBOL) {
        for (int i = 0; i < nb; i++)
            if (expr_eq((Expr*)r, (Expr*)bound[i])) return false;
        if (pmt_depends(f, r) || expr_eq((Expr*)r, (Expr*)x)) return false;
        SymbolDef* d = symtab_lookup(r->data.symbol.name);
        if (d && (d->module_temp || (d->attributes & ATTR_TEMPORARY))) return true;
        if (d && (d->is_system || d->builtin_func || d->attributes)) return false;
        return true;
    }
    if (r->type != EXPR_FUNCTION) return false;
    const Expr* h = r->data.function.head;
    if (h->type == EXPR_SYMBOL && strcmp(h->data.symbol.name, "Function") == 0 &&
        r->data.function.arg_count >= 2) {
        /* Function[v, body] / Function[{v1, ...}, body]: v is bound in body */
        const Expr* v = r->data.function.args[0];
        const Expr* nbound[PMT_MAX_BOUND];
        int n = nb < PMT_MAX_BOUND ? nb : PMT_MAX_BOUND;
        for (int i = 0; i < n; i++) nbound[i] = bound[i];
        if (v->type == EXPR_SYMBOL && n < PMT_MAX_BOUND) nbound[n++] = v;
        else if (v->type == EXPR_FUNCTION)
            for (size_t i = 0; i < v->data.function.arg_count && n < PMT_MAX_BOUND; i++)
                if (v->data.function.args[i]->type == EXPR_SYMBOL)
                    nbound[n++] = v->data.function.args[i];
        for (size_t i = 1; i < r->data.function.arg_count; i++)
            if (pmt_leaks_symbol_b(r->data.function.args[i], f, x, nbound, n)) return true;
        return false;
    }
    if (pmt_leaks_symbol_b(h, f, x, bound, nb)) return true;
    for (size_t i = 0; i < r->data.function.arg_count; i++)
        if (pmt_leaks_symbol_b(r->data.function.args[i], f, x, bound, nb)) return true;
    return false;
}

static bool pmt_leaks_symbol(const Expr* r, const Expr* f, const Expr* x) {
    return pmt_leaks_symbol_b(r, f, x, NULL, 0);
}

/* Automatic-cascade / Method -> "ParallelMixedTower" stage: evaluate the method
 * symbol and decline (continue the cascade) unless it returned an
 * antiderivative.  ParallelIntegrateMixed signals "no result" with $Failed (a
 * tower-build failure) or a List (a certificate {"not elementary", ...} or a
 * failure tuple {"failed", ...}); an antiderivative is never a bare List.  An
 * unresolved head means the .m worker never loaded. */
static Expr* try_parallelmixedtower(Expr* f, Expr* x, bool cascade) {
    /* Skip pure rational functions -- BronsteinRational owns them (see above). */
    if (pmt_is_rational_structure(f)) return NULL;
    if (pmt_has_foreign_kernel(f, x)) return NULL;   /* outside the elementary class */
    pmt_in_cascade = cascade;   /* consumed by builtin_integrate_pmt: backstop budgets */
    Expr* result = call_stage("Integrate`ParallelMixedTower", f, x);
    pmt_in_cascade = 0;         /* the builtin may have declined before consuming it */
    if (!result) return NULL;

    bool decline =
        result_is_unresolved(result, "Integrate`ParallelMixedTower") ||
        (result->type == EXPR_SYMBOL &&
         strcmp(result->data.symbol.name, "$Failed") == 0) ||
        (result->type == EXPR_FUNCTION &&
         result->data.function.head->type == EXPR_SYMBOL &&
         strcmp(result->data.function.head->data.symbol.name, "List") == 0) ||
        pmt_leaks_symbol(result, f, x) ||
        /* the worker's own call left unevaluated (its package definitions are
         * missing) is not an antiderivative */
        result_contains_head(result, "ParallelMixed`ParallelIntegrateMixed");
    if (decline) {
        expr_free(result);
        return NULL;
    }
    return result;
}

/* ------------------------------------------------------------------------- *
 * Methods: Trager and TragerKauers.
 * ------------------------------------------------------------------------- *
 * Trager's algorithm for algebraic functions over a simple radical extension
 * (Trager 1984; Integrate[f, x, Method -> "Trager"]) and the same pipeline with
 * Kauers' Groebner-basis heuristic for the logarithmic part (Kauers 2008;
 * Method -> "TragerKauers").  Both are written in Mathilda's own language in
 * src/internal/trager/ and loaded on first use.  The worker
 * Trager`TragerIntegrate[f, x, engine] returns the antiderivative (verified
 * exactly before it is returned), {"not elementary", reason} -- a proof -- or
 * {"failed", reason}.  Neither method is part of the Automatic cascade. */
static unsigned trager_load_gen = 0;   /* symtab_generation() of the attempt */
static bool trager_load_succeeded = false;

static void trager_lazy_load(void) {
    if (trager_load_gen == symtab_generation()) return;
    trager_load_gen = symtab_generation();
    if (mathilda_load_module("trager/Trager.m")) {
        trager_load_succeeded = true;
        return;
    }
    mth_message("Integrate`Trager", "nofile",
        "cannot locate src/internal/trager/Trager.m on disk.");
}

/* Shared body of the two method symbols: `engine` is "Trager" or "Kauers".
 * Same async-defer / message-suppression discipline as the ParallelMixedTower
 * symbol (see builtin_integrate_pmt); a {"not elementary", ...} certificate is
 * announced as Integrate::nonelem after both regions are popped. */
static Expr* integrate_trager_common(Expr* res, const char* engine) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    if (tc_async_region_active()) return NULL;   /* no nested heavy method */
    trager_lazy_load();
    if (!trager_load_succeeded) return NULL;

    Expr* args[3] = { expr_copy(res->data.function.args[0]),
                      expr_copy(res->data.function.args[1]),
                      expr_new_string(engine) };
    Expr* call = expr_new_function(expr_new_symbol("Trager`TragerIntegrate"), args, 3);
    tc_async_defer_push();
    mth_msg_suppress_push();
    Expr* result = evaluate(call);
    mth_msg_suppress_pop();
    tc_async_defer_pop();
    expr_free(call);
    if (pmt_is_nonelementary_certificate(result))
        integrate_announce_nonelementary(res->data.function.args[0],
                                         res->data.function.args[1]);
    return result;   /* antiderivative, or a {"not elementary"|"failed", reason} List */
}

static Expr* builtin_integrate_trager(Expr* res) {
    return integrate_trager_common(res, "Trager");
}

static Expr* builtin_integrate_trager_kauers(Expr* res) {
    return integrate_trager_common(res, "Kauers");
}

/* Method -> "Trager" / "TragerKauers" stage: decline (leave Integrate
 * unevaluated) unless the method symbol returned an antiderivative. */
static Expr* try_trager(Expr* f, Expr* x, const char* head) {
    Expr* result = call_stage(head, f, x);
    if (!result) return NULL;
    bool decline =
        result_is_unresolved(result, head) ||
        (result->type == EXPR_FUNCTION &&
         result->data.function.head->type == EXPR_SYMBOL &&
         result->data.function.head->data.symbol.name == SYM_List) ||
        result_contains_head(result, "Trager`TragerIntegrate");
    if (decline) {
        expr_free(result);
        return NULL;
    }
    return result;
}

/* ------------------------------------------------------------------------- *
 * Stage: ParallelMixedSpecial.
 * ------------------------------------------------------------------------- *
 * The SPECIAL-FUNCTION stage of the same parallel (Risch-Norman) method: the
 * ansatz of ParallelMixedTower extended with special-function kernels, so the
 * answer may contain ExpIntegralEi, LogIntegral, Erf / Erfc / Erfi, the
 * incomplete Gamma[s, z], SinIntegral / CosIntegral, and the elliptic
 * EllipticF / EllipticE / EllipticPi of a pencil.  Implemented as a Wolfram-
 * language module in src/internal/ParallelMixedSpecial.m, which is loaded INTO
 * ParallelMixed`Private` (the logrewrite.m pattern) and so reaches Part II's
 * internals by their short names -- see the file's own header.
 *
 * Contract of the worker.  IntegrateSurfaceSpecial[f, x, opts] answers with the
 * PAIR {answer, verified} on success, or a status List: {"failed", reason, ...}
 * (inconclusive), {"not elementary", ...} (a proof), {"not in class", ...} (the
 * integrand is outside the class the stage decides).  This builtin unwraps a
 * VERIFIED pair to the bare antiderivative and otherwise hands back exactly what
 * the package said -- an unverified pair included.  Nothing is invented here, so
 * the qualified-symbol surface still shows the package's own words, and
 * `Head[r] === List` remains the caller's test for "no answer". */
static unsigned pms_load_gen = 0;   /* symtab_generation() of the attempt */
static bool pms_load_succeeded = false;

static void pms_lazy_load(void) {
    if (pms_load_gen == symtab_generation()) return;
    pms_load_gen = symtab_generation();
    /* mathilda_load_module, not mathilda_run_file: the module has its own
     * load-once memo keyed on the relpath, and the special module loads Part II
     * through the same mechanism.  Going through the memo is what stops the
     * 217 KB Part II file being parsed and evaluated TWICE (~1 s each) when both
     * methods are reached in one session -- see the note on pmt_lazy_load. */
    /* Muted across the load so a build without the file reports ONE message --
     * the method-level one below, which names the method -- instead of that plus
     * the loader's own LoadModule::nofile. Load-time chatter from the package is
     * a probe's noise here, the same judgement builtin_integrate_pmt makes for
     * the per-call diagnostics; a parser error still prints, being exempt from
     * the funnel. */
    mth_msg_suppress_push();
    int ok = mathilda_load_module("ParallelMixedSpecial.m");
    mth_msg_suppress_pop();
    if (ok) {
        pms_load_succeeded = true;
        return;
    }
    mth_message("Integrate`ParallelMixedSpecial", "nofile",
        "cannot locate src/internal/ParallelMixedSpecial.m on disk.");
}

/* {answer, verified} with verified === True -> the answer; else NULL. */
static Expr* pms_verified_answer(Expr* r) {
    if (!r || r->type != EXPR_FUNCTION) return NULL;
    Expr* h = r->data.function.head;
    if (!h || h->type != EXPR_SYMBOL || h->data.symbol.name != SYM_List) return NULL;
    if (r->data.function.arg_count != 2) return NULL;
    const Expr* v = r->data.function.args[1];
    if (!v || v->type != EXPR_SYMBOL || v->data.symbol.name != SYM_True) return NULL;
    Expr* ans = r->data.function.args[0];
    r->data.function.args[0] = NULL;          /* move out; caller frees the wrapper */
    return ans;
}

/* Wall-clock budget, in seconds, that the AUTOMATIC CASCADE allows the special
 * stage -- scoping down the package's own $SpecialTimeBudget (45 s, the parity
 * value with Part II) for the duration of one cascade call.
 *
 * Why the cascade needs its own, much smaller, number.  The stage runs LAST, on
 * every integrand no earlier stage closed, and Algorithm S6's kernel search
 * escalates through six {split, retry} configurations: on an integrand with two
 * transcendental generators and no kernel that can close it, the ansatz grows
 * until the linear solve dominates.  Measured: Integrate[Sin[x^2 + Log[x]] Cos[x],
 * {x, 0, 1}] -- a plain user-level call, and one of the cases in
 * tests/test_integrate_newton_leibniz.c -- ran past 120 s before this gate, and
 * x^3 Sin[x] Log[x]^2 costs 12.1 s.  An explicit Method -> "ParallelMixedSpecial"
 * (or the qualified symbol) is a deliberate request for that search and keeps the
 * full 45 s; the cascade is not, and must get out of the way.
 *
 * Choosing the value (measured at v0.240, $SpecialTimeBudget raised to 120 so
 * the budget is not what is being timed -- the table is in
 * MATHILDA_PARALLEL_MIXED_SPECIAL_PLAN.md).  Every integrand this stage CLOSES
 * closes fast: 0.019 s (x E^x^2) to 0.187 s (Sin[x]/x), with the elliptic
 * pencils at 0.07-0.08 s.  The slowest work that is not a close is the
 * Log[x] Sin[x] / Cos[x] Log[x] pair at 3.19 / 3.34 s -- a special answer found
 * but withheld by strict mode for want of a certificate.  So 10 s leaves every
 * close a ~50x margin AND keeps that pair entirely inside the budget, which is
 * what stops the certificate work of the plan (the route from 247 to 305) from
 * needing this number changed to take effect.  It remains 4.5x below the
 * package default, which is the point.
 *
 * The stress corpus is unaffected either way: its harness drives
 * IntegrateSurfaceSpecial directly and pins both budgets to the 120 s cap. */
#define PMS_CASCADE_BUDGET_SECONDS 10

/* Armed by try_parallelmixedspecial when it is running the stage as part of the
 * Automatic cascade rather than on an explicit Method -> / qualified-symbol
 * request.  A plain int, not a counter: the stage cannot re-enter itself
 * (tc_async_region_active bars that), so there is no nesting to track.
 *
 * CONSUMED, not merely read, by the builtin below: an outer TimeConstrained can
 * siglongjmp out of the middle of the stage, which would otherwise leave the
 * flag armed and silently give the cascade's tight budget to the next EXPLICIT
 * Method -> call in the session.  Clearing it at the point of use makes that
 * unwind harmless without a save/restore dance. */
static int pms_in_cascade = 0;

static Expr* builtin_integrate_pms(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    if (tc_async_region_active()) return NULL;   /* no self re-entry / nested heavy method */
    pms_lazy_load();
    if (!pms_load_succeeded) return NULL;

    int cascade = pms_in_cascade;
    pms_in_cascade = 0;            /* consume; see the flag's own comment */

    tc_async_defer_push();
    mth_msg_suppress_push();
    Expr* raw;
    if (cascade) {
        /* Block[{$SpecialTimeBudget = n}, worker[f, x]].
         * Block, not Set: dynamic scope restores the package's own default on
         * every exit path, including the TimeConstrained unwind (v0.238 made
         * Block restore OwnValues and DownValues across that unwind), so an
         * explicit Method -> call later in the same session still gets 45 s.
         *
         * Fully qualified: the budget is exported from the ParallelMixed`
         * public section (mixed/ParallelMixed.m).  Before v0.523 the resolver
         * left every `$name` bare (issue #104), so this Blocked the bare name;
         * a wrong spelling here binds a DIFFERENT symbol and silently leaves
         * the 45 s default in force. */
        Expr* call = expr_new_function(
            expr_new_symbol("ParallelMixed`Private`IntegrateSurfaceSpecial"),
            (Expr*[]){ expr_copy(res->data.function.args[0]),
                       expr_copy(res->data.function.args[1]) }, 2);
        Expr* bind = expr_new_function(
            expr_new_symbol(SYM_Set),
            (Expr*[]){ expr_new_symbol("ParallelMixed`$SpecialTimeBudget"),
                       expr_new_integer(PMS_CASCADE_BUDGET_SECONDS) }, 2);
        Expr* locals = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ bind }, 1);
        Expr* blk = expr_new_function(expr_new_symbol(SYM_Block),
                                      (Expr*[]){ locals, call }, 2);
        raw = evaluate(blk);
        expr_free(blk);
    } else {
        raw = call_stage("ParallelMixed`Private`IntegrateSurfaceSpecial",
                         res->data.function.args[0], res->data.function.args[1]);
    }
    mth_msg_suppress_pop();
    tc_async_defer_pop();
    if (!raw) return NULL;

    Expr* ans = pms_verified_answer(raw);
    if (ans) { expr_free(raw); return ans; }

    /* Not a verified pair.  A {"not elementary", ...} certificate is a PROOF, so
     * it warns exactly as ParallelMixedTower's does; {"failed", ...} and
     * {"not in class", ...} prove nothing and stay silent. */
    if (pmt_is_nonelementary_certificate(raw))
        integrate_announce_nonelementary(res->data.function.args[0],
                                         res->data.function.args[1]);
    return raw;
}

/* Automatic-cascade / Method -> "ParallelMixedSpecial" stage.  Declines unless
 * the method returned an antiderivative, i.e. unless the builtin above unwrapped
 * a verified pair -- so a List of any shape, $Failed, and an unresolved head all
 * continue the cascade.
 *
 * Deliberately NOT gated on has_pseudoelliptic_radical, the gate that keeps
 * ParallelMixedTower off a genus>0 F/R^p curve.  That gate is right for the
 * elementary-only stage, which grinds on such a curve and closes nothing -- and
 * exactly wrong here, because an elliptic pencil IS this stage's subject:
 * 1/Sqrt[x^3 - x] is pseudo-elliptic by that test and closes as
 * -Sqrt[2] EllipticF[ArcSin[Sqrt[2]/Sqrt[1 + x]], 1/2].  The rational-structure
 * skip is kept: a pure rational function is BronsteinRational's. */
static Expr* try_parallelmixedspecial(Expr* f, Expr* x, bool cascade) {
    if (pmt_is_rational_structure(f)) return NULL;
    pms_in_cascade = cascade ? 1 : 0;   /* budget scope; see PMS_CASCADE_BUDGET_SECONDS */
    Expr* result = call_stage("Integrate`ParallelMixedSpecial", f, x);
    pms_in_cascade = 0;                 /* belt-and-braces: the builtin consumed it */
    if (!result) return NULL;

    /* An antiderivative never mentions a symbol absent from the integrand: the
     * pms_leaks check rejects a partially-evaluated package body that escaped
     * with Module locals (fp$NNN, T$NNN) -- Mathar M8041's
     * {ParallelIntegrateSpecial[fp$58989, ...], AllTrue[None, False &]} reached
     * the user that way -- and $Aborted is a timeout, not an answer. */
    bool decline =
        result_is_unresolved(result, "Integrate`ParallelMixedSpecial") ||
        (result->type == EXPR_SYMBOL &&
         (strcmp(result->data.symbol.name, "$Failed") == 0 ||
          strcmp(result->data.symbol.name, "$Aborted") == 0)) ||
        (result->type == EXPR_FUNCTION &&
         result->data.function.head->type == EXPR_SYMBOL &&
         result->data.function.head->data.symbol.name == SYM_List) ||
        pmt_leaks_symbol(result, f, x) ||
        /* the worker's own call left unevaluated (its package definitions are
         * missing) is not an antiderivative */
        result_contains_head(result, "ParallelMixed`Private`IntegrateSurfaceSpecial");
    if (decline) {
        expr_free(result);
        return NULL;
    }
    return result;
}

/* Method-option parsing.  Mirrors the canonical SYM_Method / SYM_Rule
 * idiom (see src/list.c:1480-1491 and src/facint.c:1276-1310). */
typedef enum {
    METHOD_AUTOMATIC = 0,
    METHOD_RATIONAL,
    METHOD_DERIVATIVE_DIVIDES,
    METHOD_LINEAR_RADICALS,
    METHOD_QUADRATIC_RADICALS,
    METHOD_LINEAR_RATIO_RADICALS,
    METHOD_CHEBYCHEV,
    METHOD_GOURSAT,
    METHOD_WEIERSTRASS,
    METHOD_RISCH_TRANSCENDENTAL,
    METHOD_CRCTABLE,
    METHOD_POLYLOG_LADDER,
    METHOD_HOLONOMIC,
    METHOD_BESSEL_HYPERGEOMETRIC,
    METHOD_PARALLEL_MIXED_TOWER,
    METHOD_PARALLEL_MIXED_SPECIAL,
    METHOD_TRAGER,           /* Trager's algorithm, simple radicals (src/internal/trager) */
    METHOD_TRAGER_KAUERS,    /* the same with Kauers' log-part heuristic */
    METHOD_UNDEFINED,
    METHOD_NEWTON_LEIBNIZ,   /* definite-only: selects the real-axis FTC mechanism */
    METHOD_LINE_INTEGRAL,    /* definite-only: selects the complex contour mechanism */
    METHOD_RESIDUE,          /* definite-only: selects the residue-theorem mechanism */
    METHOD_DIFF_UNDER_INT,   /* definite-only: differentiation under the integral sign */
    METHOD_RAMANUJAN,        /* definite-only: Mellin / Ramanujan Master Theorem (half-line) */
    METHOD_SYMMETRY,         /* definite-only: origin-symmetry (odd -> 0, even -> 2 half) */
    METHOD_BETA,             /* definite-only: Euler-Beta on [0,1] */
    METHOD_TRIG_POWER,       /* definite-only: Sin^m Cos^n over [0,Pi/2]/[0,Pi]/[0,2Pi] */
    METHOD_SINPOW_MONO,      /* definite-only: Sin[r x]^k / x^m on [0,Inf) */
    METHOD_OSC_POWER,        /* definite-only: Cos/Sin[b x^n] on [0,Inf) */
    METHOD_RATIONAL_LOG,     /* definite-only: R(x) Log[x]^n on [0,Inf) */
    METHOD_INTEGRAL_REP,     /* definite-only: special-function integral reps (half-line) */
    METHOD_MEIJER_CONV,      /* definite-only: Meijer G convolution of two kernels (half-line) */
    METHOD_PIECEWISE,        /* definite-only: Abs/Sign/Heaviside/UnitStep/Min/Max/Piecewise */
    METHOD_SERIES,           /* definite-only: series term-by-term on [0,1] (R(x) Log[x]^k) */
    METHOD_KING,             /* definite-only: general King reflection x -> a+b-x */
    METHOD_FEYNMAN_PARAM,    /* definite-only: parameter/Frullani-log (Sum c x^e)/Log[x] on [0,1] */
    METHOD_INVALID
} IntegrateMethod;

/* Map a method-name string to its enum, or METHOD_INVALID if unrecognised. */
static IntegrateMethod method_from_string(const char* s) {
    if (strcmp(s, "Automatic")   == 0) return METHOD_AUTOMATIC;
    if (strcmp(s, "BronsteinRational") == 0) return METHOD_RATIONAL;
    if (strcmp(s, "DerivativeDivides") == 0) return METHOD_DERIVATIVE_DIVIDES;
    if (strcmp(s, "LinearRadicals") == 0) return METHOD_LINEAR_RADICALS;
    if (strcmp(s, "QuadraticRadicals") == 0) return METHOD_QUADRATIC_RADICALS;
    if (strcmp(s, "LinearRatioRadicals") == 0) return METHOD_LINEAR_RATIO_RADICALS;
    if (strcmp(s, "ChebychevAlgebraic") == 0) return METHOD_CHEBYCHEV;
    if (strcmp(s, "GoursatAlgebraic") == 0) return METHOD_GOURSAT;
    if (strcmp(s, "Weierstrass") == 0) return METHOD_WEIERSTRASS;
    if (strcmp(s, "RischTranscendental") == 0) return METHOD_RISCH_TRANSCENDENTAL;
    if (strcmp(s, "CRCTable")    == 0) return METHOD_CRCTABLE;
    if (strcmp(s, "PolyLogLadder") == 0) return METHOD_POLYLOG_LADDER;
    if (strcmp(s, "Holonomic") == 0) return METHOD_HOLONOMIC;
    if (strcmp(s, "BesselHypergeometric") == 0) return METHOD_BESSEL_HYPERGEOMETRIC;
    if (strcmp(s, "ParallelMixedTower") == 0) return METHOD_PARALLEL_MIXED_TOWER;
    if (strcmp(s, "ParallelMixedSpecial") == 0) return METHOD_PARALLEL_MIXED_SPECIAL;
    if (strcmp(s, "Trager") == 0) return METHOD_TRAGER;
    if (strcmp(s, "TragerKauers") == 0) return METHOD_TRAGER_KAUERS;
    if (strcmp(s, "Undefined")   == 0) return METHOD_UNDEFINED;
    if (strcmp(s, "NewtonLeibniz") == 0) return METHOD_NEWTON_LEIBNIZ;
    if (strcmp(s, "LineIntegral") == 0) return METHOD_LINE_INTEGRAL;
    if (strcmp(s, "Residue") == 0 || strcmp(s, "ContourResidue") == 0)
        return METHOD_RESIDUE;
    if (strcmp(s, "DiffUnderInt") == 0 ||
        strcmp(s, "DifferentiationUnderIntegral") == 0)
        return METHOD_DIFF_UNDER_INT;
    if (strcmp(s, "Symmetry") == 0) return METHOD_SYMMETRY;
    if (strcmp(s, "Beta") == 0) return METHOD_BETA;
    if (strcmp(s, "TrigPower") == 0) return METHOD_TRIG_POWER;
    if (strcmp(s, "SinPowerMonomial") == 0) return METHOD_SINPOW_MONO;
    if (strcmp(s, "OscillatoryPower") == 0) return METHOD_OSC_POWER;
    if (strcmp(s, "RationalLog") == 0) return METHOD_RATIONAL_LOG;
    if (strcmp(s, "RamanujanMasterTheorem") == 0 || strcmp(s, "Mellin") == 0)
        return METHOD_RAMANUJAN;
    if (strcmp(s, "IntegralRepresentation") == 0) return METHOD_INTEGRAL_REP;
    if (strcmp(s, "MeijerConvolution") == 0) return METHOD_MEIJER_CONV;
    if (strcmp(s, "Piecewise") == 0) return METHOD_PIECEWISE;
    if (strcmp(s, "SeriesExpansion") == 0 || strcmp(s, "SeriesTermByTerm") == 0)
        return METHOD_SERIES;
    if (strcmp(s, "King") == 0 || strcmp(s, "Reflection") == 0) return METHOD_KING;
    if (strcmp(s, "FeynmanParameter") == 0) return METHOD_FEYNMAN_PARAM;
    return METHOD_INVALID;
}

/* Parse the `Method -> ...` option.  The RHS may be:
 *   - the symbol `Automatic` or a bare method-name string;
 *   - a list `{"<method>", subopt -> val, ...}` carrying method sub-options
 *     (mirrors the FactorInteger list-form in src/facint.c:1173).  The only
 *     recognised sub-option is `"Substitution" -> u`, valid only for
 *     "DerivativeDivides": on success `*out_sub` is set to an owned copy of `u`
 *     (the caller must free it).
 * `*out_sub` is always initialised (to NULL) before any other work. */
static IntegrateMethod parse_method_option(Expr* opt, Expr** out_sub) {
    *out_sub = NULL;
    if (opt->type != EXPR_FUNCTION) return METHOD_INVALID;
    if (opt->data.function.head->type != EXPR_SYMBOL) return METHOD_INVALID;
    const char* hd = opt->data.function.head->data.symbol.name;
    if ((hd != SYM_Rule && hd != SYM_RuleDelayed) ||
        opt->data.function.arg_count != 2) return METHOD_INVALID;
    Expr* lhs = opt->data.function.args[0];
    Expr* rhs = opt->data.function.args[1];
    if (lhs->type != EXPR_SYMBOL || lhs->data.symbol.name != SYM_Method)
        return METHOD_INVALID;

    /* Accept either a string ("Automatic") or the symbol Automatic. */
    if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic)
        return METHOD_AUTOMATIC;

    /* List form: {"<method>", subopt -> val, ...}. */
    if (rhs->type == EXPR_FUNCTION &&
        rhs->data.function.head->type == EXPR_SYMBOL &&
        rhs->data.function.head->data.symbol.name == SYM_List &&
        rhs->data.function.arg_count >= 1) {
        Expr* mname = rhs->data.function.args[0];
        if (mname->type != EXPR_STRING) return METHOD_INVALID;
        IntegrateMethod m = method_from_string(mname->data.string);
        if (m == METHOD_INVALID) return METHOD_INVALID;
        for (size_t i = 1; i < rhs->data.function.arg_count; i++) {
            Expr* so = rhs->data.function.args[i];
            if (so->type != EXPR_FUNCTION ||
                so->data.function.head->type != EXPR_SYMBOL) goto bad_subopt;
            const char* sh = so->data.function.head->data.symbol.name;
            if ((sh != SYM_Rule && sh != SYM_RuleDelayed) ||
                so->data.function.arg_count != 2) goto bad_subopt;
            Expr* skey = so->data.function.args[0];
            Expr* sval = so->data.function.args[1];
            if (skey->type == EXPR_STRING &&
                strcmp(skey->data.string, "Substitution") == 0 &&
                m == METHOD_DERIVATIVE_DIVIDES) {
                if (*out_sub) expr_free(*out_sub);
                *out_sub = expr_copy(sval);
                continue;
            }
            goto bad_subopt;
        }
        return m;
    bad_subopt:
        if (*out_sub) { expr_free(*out_sub); *out_sub = NULL; }
        return METHOD_INVALID;
    }

    if (rhs->type != EXPR_STRING) return METHOD_INVALID;
    return method_from_string(rhs->data.string);
}

/* Parse a definite integral's trailing `Method -> "..."` option.  `*name` is set
 * to a borrowed method-name string to pass through to the inner indefinite
 * Integrate (NULL for Automatic and for the definite mechanisms NewtonLeibniz /
 * LineIntegral / Residue, which select the mechanism itself).  `*mech` is set to
 * the chosen definite mechanism so the dispatcher knows whether to engage the
 * residue method and/or fall through to Newton-Leibniz.  Returns false on a
 * malformed / unrecognised option. */
static bool definite_parse_method(Expr* opt, const char** name,
                                  IntegrateMethod* mech) {
    *name = NULL;
    *mech = METHOD_AUTOMATIC;
    if (opt->type != EXPR_FUNCTION || opt->data.function.arg_count != 2)
        return false;
    if (opt->data.function.head->type != EXPR_SYMBOL) return false;
    const char* hd = opt->data.function.head->data.symbol.name;
    if (hd != SYM_Rule && hd != SYM_RuleDelayed) return false;
    Expr* lhs = opt->data.function.args[0];
    Expr* rhs = opt->data.function.args[1];
    if (lhs->type != EXPR_SYMBOL || lhs->data.symbol.name != SYM_Method) return false;
    if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic) return true;
    if (rhs->type != EXPR_STRING) return false;
    IntegrateMethod m = method_from_string(rhs->data.string);
    if (m == METHOD_INVALID) return false;
    *mech = m;
    /* NewtonLeibniz / LineIntegral / Residue / DiffUnderInt name the definite
     * mechanism itself; the actual mechanism is chosen from the spec type, so
     * they pass NULL through to the inner indefinite Integrate. */
    if (m == METHOD_AUTOMATIC || m == METHOD_NEWTON_LEIBNIZ ||
        m == METHOD_LINE_INTEGRAL || m == METHOD_RESIDUE ||
        m == METHOD_DIFF_UNDER_INT || m == METHOD_RAMANUJAN ||
        m == METHOD_SYMMETRY || m == METHOD_BETA ||
        m == METHOD_TRIG_POWER || m == METHOD_SINPOW_MONO ||
        m == METHOD_OSC_POWER || m == METHOD_RATIONAL_LOG ||
        m == METHOD_INTEGRAL_REP || m == METHOD_MEIJER_CONV ||
        m == METHOD_PIECEWISE) return true;
    *name = rhs->data.string;   /* borrowed: valid while `res` is alive */
    return true;
}

/* True iff `opt` is a Rule/RuleDelayed whose LHS is the symbol `sym`. */

/* The default TimeConstraint for the Automatic cascade, in seconds.
 *
 * Chosen by measurement, not by taste. Profiling the cascade from the inside
 * (MATHILDA_INTEGRATE_PROFILE=1) over the integrands it fails to close shows
 * `try_derivdivides` is very nearly the whole cost -- 25.05 s of a 25.06 s
 * cascade on `Sqrt[Sin[x]]/(Sin[x]^2 + 1)`, against 0.003 s for the ten stages
 * ahead of it. Timing the stage on the integrands it DOES close bounds the price
 * of cutting it off: median 0.40 s, p90 1.72 s, and only one close in 45 runs
 * longer than 3 s -- and that one (an ArcCos/Log mixed radical) is closed by
 * ParallelMixedTower three stages later in 0.21 s. So 3 s costs no answers and
 * caps the cases that currently do not terminate. Override per call with
 * `TimeConstraint -> t`, or `Infinity` for the old unbounded behaviour. */
#define INTEGRATE_DEFAULT_TIME_CONSTRAINT 3.0

/* An option whose LHS is the named (non-interned) symbol, e.g. "TimeConstraint".
 * Mirrors is_rule_with_lhs in src/simp/simp_builtins.c, which is how the same
 * option is matched on Simplify / FullSimplify / Refine. */
static bool integrate_option_is(const Expr* opt, const char* name) {
    if (!opt || opt->type != EXPR_FUNCTION || opt->data.function.arg_count != 2)
        return false;
    const Expr* h = opt->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return false;
    if (h->data.symbol.name != SYM_Rule && h->data.symbol.name != SYM_RuleDelayed)
        return false;
    const Expr* lhs = opt->data.function.args[0];
    return lhs->type == EXPR_SYMBOL && strcmp(lhs->data.symbol.name, name) == 0;
}

/* TimeConstraint value -> seconds. Infinity (symbol or DirectedInfinity) and any
 * non-positive or unreadable value mean unbounded, matching
 * simp_parse_time_budget. A {t, ...} list takes its first element. */
static double integrate_parse_time_constraint(const Expr* e) {
    if (!e) return HUGE_VAL;
    if (e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_Infinity)
        return HUGE_VAL;
    if (e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == SYM_DirectedInfinity) return HUGE_VAL;
        if (h == SYM_List && e->data.function.arg_count >= 1)
            return integrate_parse_time_constraint(e->data.function.args[0]);
    }
    double v;
    if (common_machine_real_value((Expr*)e, &v) && v > 0.0) return v;
    return HUGE_VAL;
}

static bool option_lhs_is(Expr* opt, const char* sym) {
    if (opt->type != EXPR_FUNCTION || opt->data.function.arg_count != 2) return false;
    if (opt->data.function.head->type != EXPR_SYMBOL) return false;
    const char* hd = opt->data.function.head->data.symbol.name;
    if (hd != SYM_Rule && hd != SYM_RuleDelayed) return false;
    Expr* lhs = opt->data.function.args[0];
    return lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == sym;
}

/* True if e carries an arbitrary/undefined function -- a head symbol with no
 * builtin and no DownValues (e.g. the forcing f in a Green's-function
 * convolution Integrate[K(t,s) f(s), {s,0,t}]), or an inert Derivative.  The
 * improper/parametric definite methods (residue, Ramanujan, differentiation
 * under the integral, ...) cannot close such an integrand and spend seconds
 * churning before declining, so the definite driver skips them for it. */
static bool def_has_undefined_function(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    /* An x-free subtree is a constant w.r.t. the integration: every definite
     * method folds it into an opaque coefficient, so an undefined head buried in
     * it (a parameter List like {a,b}, an opaque constant g[a]) cannot obstruct
     * the method or trigger the churn this guard exists to avoid.  Only an
     * x-DEPENDENT arbitrary function (the forcing f[x] in a Green's-function
     * convolution) must still be flagged -- so test depends-on-x first. */
    if (!depends_on_var(e, x)) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        if (h->data.symbol.name == SYM_Derivative) return true;
        SymbolDef* d = symtab_lookup(h->data.symbol.name);
        if (!d || (!d->builtin_func && !d->down_values)) return true;
    } else if (def_has_undefined_function(h, x)) {
        return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (def_has_undefined_function(e->data.function.args[i], x)) return true;
    return false;
}

/* Definite / iterated integration Integrate[f, {x,a,b}, {y,c,d}, ..., opts].
 * Reduces innermost-first (the last spec is the inner integral) so an inner
 * bound may depend on an outer variable.  Returns NULL (unevaluated) if any
 * stage cannot be evaluated. */
/* Definite stage: Meijer G convolution (src/internal/MeijerConvolution.m),
 * loaded lazily.  Half-line (0, Infinity) only; any atom result (the package's
 * $Failed) or a surviving Integrate is a decline. */
static unsigned mcv_load_gen = 0;   /* symtab_generation() of the attempt */
static bool mcv_load_succeeded = false;

static Expr* try_meijer_convolution(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions) {
    if (!(a->type == EXPR_INTEGER && a->data.integer == 0)) return NULL;
    if (!(b->type == EXPR_SYMBOL && strcmp(b->data.symbol.name, "Infinity") == 0)) return NULL;
    if (mcv_load_gen != symtab_generation()) {
        mcv_load_gen = symtab_generation();
        mth_msg_suppress_push();
        mcv_load_succeeded = mathilda_load_module("MeijerConvolution.m") != 0;
        mth_msg_suppress_pop();
        if (!mcv_load_succeeded)
            mth_message("Integrate`MeijerConvolution", "nofile",
                "cannot locate src/internal/MeijerConvolution.m on disk.");
    }
    if (!mcv_load_succeeded) return NULL;
    Expr* call = expr_new_function(expr_new_symbol("MeijerConvolution`Run"),
        (Expr*[]){ expr_copy(f), expr_copy(x),
                   assumptions ? expr_copy(assumptions) : expr_new_symbol("True") }, 3);
    mth_msg_suppress_push();
    Expr* r = evaluate(call);
    mth_msg_suppress_pop();
    expr_free(call);
    if (!r) return NULL;
    /* An atomic answer is a value (0, Pi, Catalan, ...) unless it is the
     * stage's $Failed: the MeijerG it builds can evaluate all the way to a
     * number (the J1(x) Y1(3x) integral is exactly 0). */
    if ((r->type == EXPR_SYMBOL && (r->data.symbol.name == SYM_DollarFailed ||
                                    r->data.symbol.name == SYM_Indeterminate ||
                                    r->data.symbol.name == SYM_ComplexInfinity)) ||
        result_contains_head(r, "Integrate") ||
        result_contains_head(r, "MeijerConvolution`Run") ||
        result_contains_head(r, "Throw")) {
        expr_free(r);
        return NULL;
    }
    return r;
}

static Expr* integrate_definite(Expr* res) {
    size_t argc = res->data.function.arg_count;
    Expr* f = res->data.function.args[0];

    /* Leading run of range / contour specs, then at most one trailing Method
     * option.  A contour spec `{x, z0, ..., zn}` (polyline, or any non-real
     * endpoint) is a superset of the real 3-element range spec. */
    size_t nspecs = 0;
    while (1 + nspecs < argc &&
           integrate_line_is_contour_spec(res->data.function.args[1 + nspecs]))
        nspecs++;
    if (nspecs == 0) return NULL;

    const char* method = NULL;
    IntegrateMethod mech = METHOD_AUTOMATIC;
    Expr* assumptions = NULL;   /* borrowed from res: the Assumptions option value */
    bool principal_value = false;
    size_t tail = 1 + nspecs;
    /* Trailing options, in any order: `Method -> ...` and/or `Assumptions -> ...`.
     * Any other trailing argument (unknown option, stray expression) is rejected
     * with the Integrate::method diagnostic, leaving the call unevaluated. */
    for (size_t t = tail; t < argc; t++) {
        Expr* opt = res->data.function.args[t];
        if (option_lhs_is(opt, SYM_Assumptions)) {
            assumptions = opt->data.function.args[1];   /* borrowed */
            continue;
        }
        /* TimeConstraint -> t: same option as on the indefinite form.  The
         * definite path integrates indefinitely first, so the budget reaches the
         * kernel search through that call. */
        if (integrate_option_is(opt, "TimeConstraint")) {
            integrate_arm_budget(
                integrate_parse_time_constraint(opt->data.function.args[1]));
            continue;
        }
        /* PrincipalValue -> True / False (not an interned system symbol; match
         * by name).  Selects the Cauchy principal value for interior poles. */
        if (opt->type == EXPR_FUNCTION && opt->data.function.arg_count == 2 &&
            opt->data.function.head->type == EXPR_SYMBOL &&
            (opt->data.function.head->data.symbol.name == SYM_Rule ||
             opt->data.function.head->data.symbol.name == SYM_RuleDelayed) &&
            opt->data.function.args[0]->type == EXPR_SYMBOL &&
            strcmp(opt->data.function.args[0]->data.symbol.name, "PrincipalValue") == 0) {
            Expr* rhs = opt->data.function.args[1];
            principal_value = (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_True);
            continue;
        }
        if (option_lhs_is(opt, SYM_Method)) {
            if (definite_parse_method(opt, &method, &mech)) continue;
        }
        static uint64_t last_warned_hash = 0;
        uint64_t h = expr_hash(res);
        if (h != last_warned_hash) {
            mth_message("Integrate", "method", "Method option value is not a recognised integration method name.");
            last_warned_hash = h;
        }
        return NULL;
    }

    /* Fold from the innermost (last) spec outward. */
    Expr* cur = expr_copy(f);
    for (size_t k = nspecs; k >= 1; k--) {
        Expr* spec = res->data.function.args[1 + (k - 1)];
        Expr* r;
        if (integrate_line_spec_is_complex(spec)) {
            /* Complex line / contour integral (non-real endpoint or polyline). */
            r = integrate_line_from_spec(cur, spec, method);
        } else {
            /* Real-axis definite integral.  Under Automatic (or an explicit
             * Method -> "Residue") try the residue-theorem method first: for
             * improper / periodic forms it yields cleaner, NIntegrate-verified
             * closed forms.  On NULL it falls through to Newton-Leibniz, except
             * when the user pinned "Residue" (strict: no fallback). */
            Expr* x = spec->data.function.args[0];
            Expr* a = spec->data.function.args[1];
            Expr* b = spec->data.function.args[2];
            r = NULL;
            bool diverges = false;
            /* Exact identities, ahead of every method: an empty interval
             * integrates to 0, and a DiracDelta[linear] integrand sifts to a
             * point value (impulse / Green's-function forcing).  DiracDelta is
             * otherwise inert, so these are the sole source of a value for such
             * inputs and never override an existing result. */
            if (expr_eq(a, b)) r = expr_new_integer(0);
            if (!r) r = integrate_dirac_try(cur, x, a, b);
            /* Whole-line Log[Abs[affine]]/(quadratic) -- taken BEFORE piecewise,
             * which would otherwise split at the Abs breakpoint into two
             * half-infinite pieces that grind in the parametric cascade. */
            if (!r && mech == METHOD_AUTOMATIC)
                r = integrate_logabs_quad_try(cur, x, a, b, assumptions);
            /* Piecewise integrands (Abs/Sign/HeavisideTheta/UnitStep/Min/Max/
             * Piecewise in x): split the interval at the interior breakpoints,
             * reduce each sub-interval to its single active branch, integrate,
             * and sum.  Runs early -- ahead of residue/FTC, none of which model
             * a branch switch -- and declines cheaply (no target head) for every
             * other integrand. */
            if (!r && (mech == METHOD_AUTOMATIC || mech == METHOD_PIECEWISE))
                r = integrate_piecewise_try(cur, x, a, b, assumptions);
            /* An integrand carrying an arbitrary/undefined function cannot be
             * closed by the improper/parametric methods; under Automatic they
             * are skipped for it (they churn for seconds before declining),
             * leaving the integral in its correct unevaluated form.  Pinned
             * methods are honoured regardless. */
            bool has_undef = def_has_undefined_function(cur, x);
            if (!r && !has_undef && (mech == METHOD_AUTOMATIC || mech == METHOD_RESIDUE))
                r = integrate_residue_try(cur, x, a, b, assumptions, &diverges,
                                          principal_value);
            /* The residue method conclusively found a pole on the integration
             * contour: the integral does not converge.  Emit Integrate::idiv
             * and stop -- do NOT fall through to Newton-Leibniz (which would
             * re-derive the same divergence and emit a duplicate warning). */
            if (!r && diverges) {
                integrate_emit_idiv(cur, a, b);
                expr_free(cur);
                return NULL;
            }
            /* Origin-symmetry reduction, after residue (which keeps its clean
             * closed forms for the families it owns) but before the general FTC:
             * an odd integrand over a symmetric interval is 0 with no
             * antiderivative, and an even one halves the interval -- often
             * unlocking the half-line Ramanujan/Mellin path for a non-rational
             * even integrand that residue does not own.  It only claims a value
             * once the half converges, so a divergent principal value is never
             * reported as 0, and NULL always falls through unchanged. */
            if (!r && !has_undef && (mech == METHOD_AUTOMATIC || mech == METHOD_SYMMETRY))
                r = integrate_symmetry_try(cur, x, a, b, assumptions, principal_value);
            /* Deterministic, correct-by-construction transform reductions run
             * BEFORE the Newton-Leibniz FTC search (confidence-first ordering,
             * the rule the Weierstrass and Symmetry placements already follow):
             * each is gated by a cheap structural shape test and closes exactly
             * the half-line / canonical-interval family it owns, so a
             * symbolic-exponent integrand (x^(s-1) f, Sin^a Cos^b, ...) yields
             * its closed form at once -- instead of first letting FTC grind a
             * doomed elementary-antiderivative search (Risch / ParallelMixed*)
             * for tens of seconds and emit a spurious Integrate::nonelem.  Each
             * declines (NULL, cheaply: the interval/shape check is first) on a
             * shape it does not own, so FTC below still owns the integer-power
             * and genuinely elementary cases. */
            /* Euler-Beta: x^(k-1)(1-x)^(l-1) [Log weights] on [0,1]; then the
             * Euler->Beta*2F1 generalization with one extra linear factor
             * (alpha+beta x)^e. */
            /* Half-line x^(c-1) Log[x]^k / (1 +/- x^d): d^k/ds^k of the base Csc/Cot
             * Mellin at s = c.  Tightly gated to [0,Inf) and that exact shape.  Runs
             * after residue (its convergence prover needs the assumption context
             * residue sets up) but before the [0,1] Beta/series block. */
            if (!r && !has_undef && mech == METHOD_AUTOMATIC)
                r = integrate_mellin_logpow_try(cur, x, a, b, principal_value, assumptions);
            /* Log-weighted rational x^(c-1) Log[x]^n/(1-x^m) on [0,1] -- the l=0
             * limit the Beta stage cannot represent; runs first so it also
             * avoids the spurious Beta[c,0]=ComplexInfinity Infinity::indet. */
            if (!r && !has_undef && (mech == METHOD_AUTOMATIC || mech == METHOD_BETA))
                r = integrate_ratlog01_try(cur, x, a, b, assumptions);
            /* Series term-by-term on [0,1]: R(x) Log[x]^k over a root-of-unity
             * denominator 1 +/- x^d.  The alternating (1 + x^d) sibling of
             * ratlog01 above (which owns the 1 - x^m HurwitzZeta case): expand
             * the geometric series, integrate term-by-term, resum with Sum[]
             * (-> Catalan / DirichletBeta / Zeta).  Runs in the pre-FTC [0,1]
             * block -- these integrands have non-elementary (PolyLog)
             * antiderivatives, so letting Newton-Leibniz grind them would only
             * emit a spurious Integrate::nonelem; the cyclotomic probe declines
             * non-matches instantly and the value is NIntegrate-cross-checked. */
            if (!r && !has_undef && (mech == METHOD_AUTOMATIC || mech == METHOD_SERIES))
                r = integrate_series_try(cur, x, a, b, assumptions);
            /* Parameter / Frullani-logarithm family on [0,1]:
             * (Sum c_i x^(e_i))/Log[x] -> Sum c_i Log[e_i+1] (Sum c_i = 0),
             * i.e. differentiation under the integral sign w.r.t. the exponent
             * with the inner x-power integral taken directly as 1/(e+1).  Cheap
             * structural gate (needs a 1/Log[x] weight); closes e.g.
             * (x^a-1)/Log[x] -> Log[a+1] that the FTC/Beta cascade leaves
             * unevaluated.  Correct-by-construction + NIntegrate-cross-checked. */
            if (!r && !has_undef && (mech == METHOD_AUTOMATIC || mech == METHOD_FEYNMAN_PARAM))
                r = integrate_feynparam_try(cur, x, a, b, assumptions);
            if (!r && !has_undef && (mech == METHOD_AUTOMATIC || mech == METHOD_BETA))
                r = integrate_beta_try(cur, x, a, b, assumptions);
            if (!r && !has_undef && (mech == METHOD_AUTOMATIC || mech == METHOD_BETA))
                r = integrate_euler_2f1_try(cur, x, a, b, assumptions);
            /* Half-Beta: x^(c-1)(1-x^2)^e Log[x]^i on [0,1] (substitution t=x^2),
             * the quadratic-radical shape the plain Beta stage cannot match. */
            if (!r && !has_undef && (mech == METHOD_AUTOMATIC || mech == METHOD_BETA))
                r = integrate_halfbeta_try(cur, x, a, b, assumptions);
            /* Sin^m Cos^n over a canonical trig interval ([0,Pi/2], [0,Pi], [0,2Pi]). */
            if (!r && !has_undef && (mech == METHOD_AUTOMATIC || mech == METHOD_TRIG_POWER))
                r = integrate_trigpower_try(cur, x, a, b, assumptions);
            /* Log[Sin[x]] / Log[Cos[x]] over a canonical [0, Pi/2] or [0, Pi]
             * interval (the classical log-sine values -(Pi/2)Log2 / -Pi Log2). */
            if (!r && !has_undef && mech == METHOD_AUTOMATIC)
                r = integrate_logtrig_try(cur, x, a, b, assumptions);
            /* Poisson log integral Integrate[Log[A + B Cos[x]], {x,0,Pi}]. */
            if (!r && !has_undef && mech == METHOD_AUTOMATIC)
                r = integrate_poissonlog_try(cur, x, a, b, assumptions);
            /* King's-property reflection: Integrate[x g(x), {x,0,c}] with g
             * invariant under x -> c-x (e.g. x Sin[x]/(1+Cos[x]^2) on [0,Pi]). */
            if (!r && !has_undef && mech == METHOD_AUTOMATIC)
                r = integrate_kingreflect_try(cur, x, a, b, assumptions);
            /* General King reflection x -> a+b-x: when f(x)+f(a+b-x) collapses
             * to a constant c (e.g. 1/(1+Tan[x]^Sqrt[2]) on [0,Pi/2], c=1) the
             * value is c (b-a)/2.  Runs before FTC so these no-elementary-
             * antiderivative integrands are closed without a spurious
             * Integrate::nonelem; constancy is numerically screened and the
             * value NIntegrate-cross-checked, so a mis-detection declines. */
            if (!r && !has_undef && (mech == METHOD_AUTOMATIC || mech == METHOD_KING))
                r = integrate_king_try(cur, x, a, b, assumptions);
            /* Raabe's integral: Integrate[Log[Gamma[x]], {x, t, t+1}]. */
            if (!r && !has_undef && mech == METHOD_AUTOMATIC)
                r = integrate_raabe_try(cur, x, a, b, assumptions);
            /* Mellin / Ramanujan Master Theorem: half-line ∫₀^∞ x^(s-1) f(x) dx
             * of a transcendental f (Gaussian moments, Gamma/Bessel/trig
             * transforms) that residue does not own.  Under Automatic it runs
             * before FTC and DiffUnderInt; the pinned mechanism has no fallback. */
            if (!r && !has_undef && (mech == METHOD_AUTOMATIC || mech == METHOD_RAMANUJAN ||
                       mech == METHOD_OSC_POWER))
                r = integrate_ramanujan_try(cur, x, a, b, assumptions);
            /* Principal value of x^(c-1)/(1-x^d) (on-axis pole at x=1; the d-th
             * roots of unity defeat the residue root-finder on a symbolic d). */
            if (!r && !has_undef && principal_value && mech == METHOD_AUTOMATIC)
                r = integrate_mellinpv_oneminus_try(cur, x, a, b, principal_value, assumptions);
            /* Integral-representation recognizer: a half-line integrand that IS a
             * classical integral representation of a special function
             * (Laplace-Bessel, Bessel-K, Airy) — no elementary antiderivative, so
             * the Mellin/residue/FTC families do not close it, but the closed form
             * is immediate once the shape matches and its convergence gate is
             * proved.  After Ramanujan (the general Mellin gets first crack), before
             * FTC/DiffUnderInt. */
            if (!r && !has_undef && (mech == METHOD_AUTOMATIC || mech == METHOD_INTEGRAL_REP))
                r = integrate_intrep_try(cur, x, a, b, assumptions);
            /* Meijer G convolution: C x^alpha K1(x) K2(x) on (0, Infinity) for
             * Bessel J/Y/K, exponential, Gaussian, trig and (a^2+x^2)^-r kernels
             * at arbitrary scales (src/internal/MeijerConvolution.m).  After
             * the table stages (they own their closed forms), before FTC. */
            if (!r && !has_undef && (mech == METHOD_AUTOMATIC || mech == METHOD_MEIJER_CONV))
                r = try_meijer_convolution(cur, x, a, b, assumptions);
            /* Newton-Leibniz (FTC) unless the user pinned Residue, the
             * parameter-differentiation mechanism, the Ramanujan/Mellin
             * mechanism, the symmetry mechanism, or a Beta mechanism.  Runs
             * after the deterministic reductions above so it only sees the
             * integrand classes none of them owns. */
            if (!r && mech != METHOD_RESIDUE && mech != METHOD_DIFF_UNDER_INT &&
                mech != METHOD_RAMANUJAN && mech != METHOD_SYMMETRY &&
                mech != METHOD_BETA && mech != METHOD_TRIG_POWER &&
                mech != METHOD_SINPOW_MONO && mech != METHOD_OSC_POWER &&
                mech != METHOD_RATIONAL_LOG && mech != METHOD_INTEGRAL_REP &&
                mech != METHOD_MEIJER_CONV &&
                mech != METHOD_PIECEWISE && mech != METHOD_SERIES &&
                mech != METHOD_KING && mech != METHOD_FEYNMAN_PARAM)
                r = integrate_newton_leibniz_try_pv(cur, x, a, b, method,
                                                    principal_value, assumptions);
            /* Sin[r x]^k / x^m half-line (ssp) and R(x) Log[x]^n (log*rat):
             * under Automatic these are pre-passes inside Ramanujan; the pinned
             * methods route directly here. */
            if (!r && mech == METHOD_SINPOW_MONO)
                r = integrate_sinpowmono_try(cur, x, a, b, assumptions);
            if (!r && mech == METHOD_RATIONAL_LOG)
                r = integrate_ratlogpow_try(cur, x, a, b, assumptions);
            /* Differentiation under the integral sign: last resort under
             * Automatic (parameter-dependent improper/periodic integrals that
             * residue and FTC cannot close), or the pinned mechanism.  Skipped
             * under PrincipalValue: if the PV-capable stages (residue, N-L)
             * already declined, this ordinary-integral method would compute a
             * value ACROSS the interior pole (a divergent / complex artifact),
             * so a PV request that they cannot close is left unevaluated. */
            if (!r && !has_undef && !principal_value &&
                (mech == METHOD_AUTOMATIC || mech == METHOD_DIFF_UNDER_INT))
                r = integrate_diffunderint_try(cur, x, a, b, assumptions);
            /* A real-axis definite integral must evaluate to a value free of
             * its integration variable: a result still containing x is a leaked
             * antiderivative (a stage returned F(x) instead of F(b)-F(a)), which
             * is always wrong -- discard it and leave the integral unevaluated
             * rather than return an expression in the dummy variable. */
            if (r && depends_on_var(r, x)) { expr_free(r); r = NULL; }
        }
        expr_free(cur);
        if (!r) return NULL;
        cur = r;
    }
    /* Special-function closed forms.  A method stage can close the integral in
     * a representation the user should not see: MeijerConvolution answers in
     * MeijerG, Newton-Leibniz leaves PolyLog/Log values at roots of unity
     * (E^(I Pi/4), I, ...) whose conjugate pairs cancel to Catalan / Zeta.
     * fe_reduce_definite_result rewrites those through FunctionExpand's rules
     * and keeps a rewrite only when it is a genuine closed form. */
    if (cur) cur = fe_reduce_definite_result(cur, assumptions);
    /* Assumption-aware result cleanup.  When the call carried an explicit
     * Assumptions option, Refine the folded result so a definite integral
     * returns the reduced closed form a method stage left equal-but-unsimplified
     * (Sqrt[(a^2+b^2)/a^2] -> Sqrt[a^2+b^2]/a, Sqrt[p^2] -> p, Abs/Sign under a
     * proved sign, ...).  Gated on `assumptions`, so no assumption-free call
     * changes form.  On any failure the original result is kept. */
    if (cur && assumptions) {
        Expr* call = expr_new_function(expr_new_symbol("Refine"),
            (Expr*[]){ expr_copy(cur), expr_copy(assumptions) }, 2);
        Expr* refined = evaluate(call);
        expr_free(call);
        /* Keep the refinement unless it INTRODUCED a complex value -- Refine can
         * strip a Re[] that was protecting realness (Re[CosIntegral[a] -
         * CosIntegral[-a]] = 0 -> the complex CosIntegral[a] - CosIntegral[-a]),
         * which would turn a correct real principal value into a wrong one. */
        if (refined) {
            if (nl_numeric_is_complex(refined) && !nl_numeric_is_complex(cur))
                expr_free(refined);                 /* reject: refinement went complex */
            else { expr_free(cur); cur = refined; }
        }
    }
    return cur;
}

/* Nesting depth of the method cascade (see integrate.h). */
int g_integrate_depth = 0;

/* ------------------------------------------------------------------------
 * TimeConstraint: the cascade's wall-clock budget.
 *
 * Owned HERE, not by any one stage, and checked between stages at EVERY depth.
 * That placement is forced. Bounding a single stage's own loops cannot bound it,
 * because the expensive work is a recursive cascade DESCENT the stage initiates:
 * try_direct_kernel -> Integrate[qu, u] re-enters the switch, and the stages that
 * then run (try_linearity above all -- 43.5% of the measured decline cost, and
 * itself a per-term recursion) know nothing about the caller's deadline.
 * Measured: with the budget inside derivdivides only, `TimeConstraint -> 3` let a
 * case run 30 s because one nested descent ran 27 s uninterrupted, while `-> 1`
 * stopped at 1.9 s purely because the deadline had already passed before that
 * descent began. Checked between stages at every depth, every unit of work is
 * bounded by construction.
 *
 * Armed once per top-level evaluation (eval_toplevel_id()), not per
 * builtin_integrate entry: the evaluator's fixed point re-enters several times
 * for one user-level Integrate, and re-arming would multiply the budget by the
 * number of passes.
 * ---------------------------------------------------------------------- */
static double   g_intg_deadline   = HUGE_VAL;
static uint64_t g_intg_budget_ep  = 0;
static bool     g_intg_budget_set = false;

static double intg_mono_seconds(void);     /* defined with the profiler below */

static void integrate_arm_budget(double seconds) {
    uint64_t ep = eval_toplevel_id();
    if (g_intg_budget_set && ep == g_intg_budget_ep) return;
    g_intg_budget_ep = ep;
    g_intg_budget_set = true;
    g_intg_deadline = (seconds > 0.0 && seconds < HUGE_VAL)
                          ? intg_mono_seconds() + seconds : HUGE_VAL;
}

/* Kept for the in-stage bounds to consult via integrate.h; the cascade itself
 * deliberately does NOT short-circuit on it -- see the note on IP_STAGE. */
MATHILDA_MAYBE_UNUSED static bool integrate_budget_spent(void) {
    return g_intg_deadline < HUGE_VAL && intg_mono_seconds() >= g_intg_deadline;
}


/* ------------------------------------------------------------------------
 * MATHILDA_INTEGRATE_PROFILE=1 -- per-stage cost of the Automatic cascade.
 *
 * WHAT THIS EXISTS TO ANSWER. A cascade stage is charged for the inputs it
 * DECLINES, not the ones it closes: a stage that closes is the fast path by
 * construction, while every integrand that reaches it and cannot be closed pays
 * the stage's full search. So "is the cascade ordered wrong?" is a question
 * about decline cost by position, and that cannot be answered from outside --
 * the pinned `Integrate`<Method>` surface measures something DIFFERENT, because
 * some cascade call sites are deliberately cheaper than the explicit method:
 * `try_derivdivides` routes a PSEUDO-ELLIPTIC integrand to the direct-quotient
 * path only, and `try_parallelmixedtower` is skipped for that shape entirely. So
 * a pinned-surface measurement is representative for most integrands and wrong
 * for those, and only an inside measurement can tell you which you have. (Here
 * the two agree: on Sqrt[Sin[x]]/(Sin[x]^2 + 1), which is not pseudo-elliptic
 * and so takes the full path, the pinned surface reports 15.0 s against its 15 s
 * cap and the cascade 25.05 s against a 25 s cap -- the same stage, saturating
 * whatever bound it is given.)
 *
 * Costs nothing when off: one cached int test per stage, no clock call. Same
 * family of knob as MATHILDA_PACK_DIAG / MATHILDA_NO_PACK.
 * ---------------------------------------------------------------------- */
typedef enum {
    IP_UNDEFINED = 0, IP_STEP, IP_RATIONAL, IP_LINRAD, IP_QUADRAD, IP_LINRATIORAD,
    IP_CHEBYCHEV, IP_FRESNEL, IP_GAMMAPOWER, IP_LINEARITY, IP_LINARG,
    IP_WEIERSTRASS, IP_DERIVDIVIDES, IP_RISCHTRANS, IP_CRCTABLE, IP_POLYLOG,
    IP_HOLONOMIC, IP_BESSELHYPER, IP_LOGBYPARTS, IP_PMT, IP_GOURSAT, IP_PMS,
    IP_NSTAGES
} IntegrateProfStage;

static const char* const ip_stage_name[IP_NSTAGES] = {
    "01 Undefined", "01b UnitStep", "02 BronsteinRational", "03 LinearRadicals",
    "04 QuadraticRadicals", "05 LinearRatioRadicals", "06 ChebychevAlgebraic",
    "07 Fresnel", "08 GammaPower", "09 Linearity", "10 LinearArgSubst",
    "11 Weierstrass", "12 DerivativeDivides", "13 RischTranscendental",
    "14 CRCTable", "14b PolyLogLadder", "14c Holonomic", "14d BesselHypergeometric", "15 LogByParts",
    "16 ParallelMixedTower", "17 GoursatAlgebraic", "18 ParallelMixedSpecial"
};
static long   g_ip_entries = 0;   /* cascade passes over one integrand */
static double g_ip_secs[IP_NSTAGES];
static long   g_ip_calls[IP_NSTAGES];
static long   g_ip_closes[IP_NSTAGES];
static long   g_ip_aborted[IP_NSTAGES];
static int    g_ip_on = -1;              /* -1 = env not yet read */

/* Monotonic seconds. CLOCK_MONOTONIC is immune to NTP steps; falls back to the
 * CPU clock when the syscall is unavailable. Mirrors simp_mono_seconds() in
 * src/simp/simp_util.c and dt_wall_seconds() in datetime.c -- kept local rather
 * than reaching into a simp-internal header. */
static double ip_mono_seconds(void);
static double intg_mono_seconds(void) { return ip_mono_seconds(); }

static double ip_mono_seconds(void) {
    return mth_mono_seconds();   /* portable monotonic clock; see src/mtime.c */
}

static void ip_flush_inflight(void);

static void ip_dump(void) {
    ip_flush_inflight();
    double tot = 0.0;
    for (int i = 0; i < IP_NSTAGES; i++) tot += g_ip_secs[i];
    fprintf(stderr, "\n=== Integrate cascade profile (outermost frame only) ===\n");
    fprintf(stderr, "  %-24s %7s %7s %7s %10s %9s %7s\n",
            "stage", "calls", "closes", "abort", "total(s)", "mean(s)", "share");
    for (int i = 0; i < IP_NSTAGES; i++) {
        if (!g_ip_calls[i]) continue;
        fprintf(stderr, "  %-24s %7ld %7ld %7ld %10.3f %9.4f %6.1f%%\n",
                ip_stage_name[i], g_ip_calls[i], g_ip_closes[i], g_ip_aborted[i],
                g_ip_secs[i],
                g_ip_secs[i] / (double)g_ip_calls[i],
                tot > 0.0 ? 100.0 * g_ip_secs[i] / tot : 0.0);
    }
    fprintf(stderr, "  %-24s %7ld %7s %7s %10.3f\n", "TOTAL", g_ip_entries, "", "", tot);
    fprintf(stderr, "  (TOTAL calls = cascade passes over the SAME integrand:"
                    " the evaluator fixed point re-enters)\n");
    fprintf(stderr, "========================================================\n");
}

static bool ip_enabled(void) {
    if (g_ip_on < 0) {
        const char* s = getenv("MATHILDA_INTEGRATE_PROFILE");
        g_ip_on = (s && *s && strcmp(s, "0") != 0) ? 1 : 0;
        if (g_ip_on) atexit(ip_dump);
    }
    return g_ip_on == 1;
}

/* Accounted only at the OUTERMOST cascade frame.  A stage that recurses into
 * Integrate would otherwise be charged its callees' time as well, and the
 * per-stage totals would sum past the wall clock -- which is exactly what the
 * first cut of this profiler did (36.7 s of stage time under a 25 s bound).
 * Note the cascade can still run several times for ONE user-level Integrate:
 * the evaluator's fixed point re-evaluates the expression, so g_ip_entries
 * counts those passes and the per-pass mean is total/entries, not total/calls. */
/* The stage currently in flight at the outermost frame, and when it started.
 * REQUIRED, not bookkeeping: a stage killed by the caller's TimeConstrained is
 * unwound by siglongjmp, so its ip_account() never runs and the most expensive
 * stage on every aborting integrand -- the only ones whose cost we care about --
 * would be invisible. ip_flush_inflight() charges it instead, called when the
 * next outermost pass begins and again at exit. */
static int    g_ip_cur    = -1;
static double g_ip_cur_t0 = 0.0;

static void ip_flush_inflight(void) {
    if (g_ip_cur < 0) return;
    g_ip_secs[g_ip_cur] += ip_mono_seconds() - g_ip_cur_t0;
    g_ip_calls[g_ip_cur]++;              /* an abandoned stage closed nothing */
    g_ip_aborted[g_ip_cur]++;
    g_ip_cur = -1;
}

static void ip_account(int stage, double t0, bool closed) {
    if (g_integrate_depth != 1) return;
    g_ip_secs[stage] += ip_mono_seconds() - t0;
    g_ip_calls[stage]++;
    if (closed) g_ip_closes[stage]++;
}

/* NOTE -- this macro briefly also skipped a stage once the TimeConstraint was
 * spent, i.e. it TERMINATED the cascade. Measured over 150 integrands that cost
 * EIGHT closes (125 -> 117) and is the wrong semantics: the stages after the
 * expensive one are cheap (PMT 0.27 s, PMS 0.34 s mean) and are often the ones
 * that close the integrand, so abandoning them to save time on a stage that
 * already failed trades an answer for nothing. A budget belongs AROUND one
 * stage's own search, after which the cascade continues -- which needs the hard
 * bound discussed in call_stage_bounded's removal note.
 *
 * Time one cascade stage. Only the OUTERMOST frame is accounted: a stage that
 * recurses into Integrate would otherwise be charged its callees' time too and
 * the shares would not sum. */
#define IP_STAGE(stage, callexpr)                                            \
    do {                                                                     \
        if (!result) {                                                       \
            if (ip_enabled()) {                                              \
                double _ip_t0 = ip_mono_seconds();                           \
                if (g_integrate_depth == 1) {                                \
                    g_ip_cur = (stage); g_ip_cur_t0 = _ip_t0;                \
                }                                                            \
                result = (callexpr);                                         \
                if (g_integrate_depth == 1) g_ip_cur = -1;                   \
                ip_account((stage), _ip_t0, result != NULL);                 \
            } else {                                                         \
                result = (callexpr);                                         \
            }                                                                \
        }                                                                    \
    } while (0)

/* Speculative-integration suppression counter (see integrate.h). */
int g_integrate_quiet = 0;

/* Heavy-tail-stage suppression counter (see integrate.h).  Raised by LogByParts
 * around its recursive sub-integrals so a non-closing reduction declines promptly
 * instead of paying ParallelMixedTower/Special twice. */
int g_integrate_no_special = 0;

/* De-dup flag for the user-facing Integrate::nonelem diagnostic within ONE
 * top-level Integrate cascade.  Two cascade stages can each PROVE the same
 * integrand non-elementary -- RischTranscendental (line ~1220) and, right after
 * it, ParallelMixedTower (line ~1226) -- so without this guard
 * Integrate[Exp[x^2], x] would print the identical warning twice.  Reset to 0
 * when a fresh top-level cascade begins (g_integrate_depth crosses 0 -> 1);
 * set to 1 by the first stage that speaks.  The direct qualified-symbol surface
 * runs at depth 0 with only one method active and does not consult it. */
static int g_integrate_nonelem_announced = 0;

/* Deferred payload for that diagnostic inside a cascade.  A stage that PROVES
 * non-elementarity is not the last word: a LATER stage may still produce a
 * non-elementary antiderivative in closed form, and since v0.241 it routinely
 * does -- Integrate[Log[x] Sin[x], x] is proved non-elementary by
 * RischTranscendental and then ANSWERED CosIntegral[x] - Cos[x] Log[x] by
 * ParallelMixedSpecial, three stages later.  Warning and then answering is not
 * incoherent (the answer is indeed not elementary) but it is noise no other CAS
 * emits, so the message is held until the cascade is done and printed only if
 * nothing answered.  Held as strings, not Expr*, so no ownership crosses the
 * cascade.  intg_nonelem_pending_clear() drops a held message without printing;
 * it also runs on cascade ENTRY, so a message stranded by a TimeConstrained
 * siglongjmp cannot surface against an unrelated later integrand. */
static char* g_integrate_nonelem_pending_f = NULL;
static char* g_integrate_nonelem_pending_x = NULL;

static void intg_nonelem_pending_clear(void) {
    free(g_integrate_nonelem_pending_f); g_integrate_nonelem_pending_f = NULL;
    free(g_integrate_nonelem_pending_x); g_integrate_nonelem_pending_x = NULL;
}

/* Nesting depth of the DEFINITE driver (integrate_definite_scoped).  A definite
 * integral first tries an indefinite antiderivative (Newton-Leibniz), and that
 * probe runs a complete top-level indefinite cascade (depth 0 -> 1 -> 0).  If
 * the probe proves the integrand non-elementary, the cascade's exit would flush
 * Integrate::nonelem -- even though a later DEFINITE stage (differentiation
 * under the integral sign, Mellin, contour residues, ...) may still answer, as
 * Integrate[Log[1 + a^2 x^2]/(1 + x^2), {x, 0, Infinity}, Assumptions -> a > 0]
 * = Pi Log[1 + a] does.  That spurious message is not cosmetic: it marks the
 * evaluation as failed, so Check[] returns its failure branch for a call that
 * succeeded.  While this is > 0 the indefinite exit therefore HANDS its held
 * message up to the definite-level payload below instead of printing it; the
 * outermost definite frame prints it only if the whole definite call declined. */
static int   g_integrate_definite_depth = 0;
static char* g_integrate_definite_pending_f = NULL;
static char* g_integrate_definite_pending_x = NULL;

static void intg_definite_pending_clear(void) {
    free(g_integrate_definite_pending_f); g_integrate_definite_pending_f = NULL;
    free(g_integrate_definite_pending_x); g_integrate_definite_pending_x = NULL;
}

/* Adopt (fs, xs) as the definite-level held message, keeping the FIRST one a
 * definite call produced (the same once-per-call contract as the cascade's
 * de-dup flag).  Takes ownership of both strings. */
static void intg_definite_pending_adopt(char* fs, char* xs) {
    if (g_integrate_definite_pending_f || !fs || !xs) { free(fs); free(xs); return; }
    g_integrate_definite_pending_f = fs;
    g_integrate_definite_pending_x = xs;
}

static void intg_nonelem_pending_flush(void) {
    if (g_integrate_definite_depth > 0) {
        /* Inside a definite integral: defer to its outermost frame. */
        intg_definite_pending_adopt(g_integrate_nonelem_pending_f,
                                    g_integrate_nonelem_pending_x);
        g_integrate_nonelem_pending_f = NULL;
        g_integrate_nonelem_pending_x = NULL;
        return;
    }
    if (g_integrate_nonelem_pending_f && g_integrate_nonelem_pending_x)
        mth_message("Integrate", "nonelem",
                    "The integrand %s has no antiderivative elementary in %s.",
                    g_integrate_nonelem_pending_f, g_integrate_nonelem_pending_x);
    intg_nonelem_pending_clear();
}

static Expr* integrate_definite(Expr* res);

/* integrate_definite, bracketed by the definite-level nonelem deferral (see
 * g_integrate_definite_depth).  Only the outermost definite frame decides: a
 * nested definite sub-integral's decline is not the user's answer. */
static Expr* integrate_definite_scoped(Expr* res) {
    if (g_integrate_definite_depth == 0)
        intg_definite_pending_clear();   /* drop anything a longjmp stranded */
    g_integrate_definite_depth++;
    Expr* r = integrate_definite(res);
    g_integrate_definite_depth--;
    if (g_integrate_definite_depth == 0) {
        char* fs = g_integrate_definite_pending_f;
        char* xs = g_integrate_definite_pending_x;
        if (!r && fs && xs)
            mth_message("Integrate", "nonelem",
                        "The integrand %s has no antiderivative elementary in %s.",
                        fs, xs);
        intg_definite_pending_clear();
    }
    return r;
}

/* Emit the user-facing "no elementary antiderivative" diagnostic for the
 * ORIGINAL integrand, shared by every method that can PROVE non-elementarity
 * (currently RischTranscendental and ParallelMixedTower).  The caller is
 * responsible for having established the proof (a RischTranscendental field
 * decision, or a {"not elementary", ...} certificate from the ParallelMixedTower
 * package); this routine only handles the gating and the message so the two
 * methods stay byte-identical and cannot double-print.
 *
 * Gating mirrors the historical RischTranscendental site: silent for internal
 * recursion (depth > 1, where naming a gensym variable would be meaningless),
 * silent under speculative integration (g_integrate_quiet, e.g. the Kovacic
 * solver / DSolve), and at most once per top-level cascade (depth == 1).  Inside
 * a definite integral the decision is deferred to its outermost frame (see
 * g_integrate_definite_depth).  Routed through the mth_message funnel, so Quiet[]
 * silences it and Check[] catches it. */
void integrate_announce_nonelementary(Expr* f, Expr* x) {
    if (g_integrate_depth > 1) return;    /* internal recursion: stay silent   */
    if (g_integrate_quiet != 0) return;   /* speculative / DSolve suppression   */
    char* fs = expr_to_string(f);
    char* xs = expr_to_string(x);
    if (g_integrate_depth == 1 && fs && xs) {  /* cascade: at most one stage speaks */
        if (g_integrate_nonelem_announced) { free(fs); free(xs); return; }
        g_integrate_nonelem_announced = 1;
        /* HOLD it: a later stage may still answer.  Flushed at the cascade exit
         * only when nothing did (see the comment on the pending payload). */
        intg_nonelem_pending_clear();
        g_integrate_nonelem_pending_f = fs;   /* ownership moves to the payload */
        g_integrate_nonelem_pending_x = xs;
        return;
    }
    /* depth 0 -- the direct qualified-symbol surface, one method, no cascade to
     * change its mind: say it now -- unless a definite integral is probing it,
     * in which case its outermost frame decides (see g_integrate_definite_depth). */
    if (g_integrate_definite_depth > 0) { intg_definite_pending_adopt(fs, xs); return; }
    mth_message("Integrate", "nonelem", "The integrand %s has no antiderivative elementary in %s.", fs ? fs : "?", xs ? xs : "?");
    free(fs); free(xs);
}

/* Fail-memo: the indefinite integrands the method cascade has FAILED on
 * (returned NULL for) so far in the CURRENT top-level evaluation.
 *
 * The evaluator's fixed-point loop re-invokes builtin_integrate a second time
 * on an integrand whose only change since the first attempt was Orderless
 * canonicalisation of its factors -- e.g. Integrate[Log[x] Log[Log[x]]/x^2, x]
 * enters the full cascade twice with a byte-identical integrand. That second
 * pass re-runs the entire (multi-second, Eliminate/Solve-driven) search and
 * re-emits every diagnostic, for a result the first pass already proved is
 * NULL. Since the cascade is deterministic, a re-entry with the same
 * (integrand, variable, method) inside the SAME top-level evaluation must yield
 * the same NULL, so we short-circuit it.
 *
 * The table is a small fixed array keyed on the top-level evaluation id
 * (eval_toplevel_id): when a new user command starts, `intg_fail_epoch` no
 * longer matches and the table is emptied -- so it self-invalidates every
 * command (where a changed definition could legitimately make the same integral
 * solvable) with no explicit reset, and needs no dynamic allocation. Several
 * distinct failing integrals in one command (a List of Integrate calls, which
 * itself re-evaluates) each get their own slot, so none re-emits when the outer
 * structure takes a second pass. Beyond the capacity the memo simply stops
 * helping (extra integrals re-run) -- a graceful degradation, never wrong. */
#define INTG_FAIL_SLOTS 32
static uint64_t intg_fail_epoch = 0;   /* eval_toplevel_id the table belongs to */
static int      intg_fail_count = 0;   /* live entries in intg_fail_tab          */
static struct { uint64_t hf, hx; int method; } intg_fail_tab[INTG_FAIL_SLOTS];

/* Drop the table if we have crossed into a new top-level evaluation. */
static void intg_fail_sync_epoch(uint64_t tid) {
    if (tid != intg_fail_epoch) { intg_fail_epoch = tid; intg_fail_count = 0; }
}
static bool intg_fail_seen(uint64_t hf, uint64_t hx, int method) {
    for (int i = 0; i < intg_fail_count; i++)
        if (intg_fail_tab[i].hf == hf && intg_fail_tab[i].hx == hx
            && intg_fail_tab[i].method == method)
            return true;
    return false;
}
static void intg_fail_record(uint64_t hf, uint64_t hx, int method) {
    if (intg_fail_seen(hf, hx, method)) return;      /* already noted */
    if (intg_fail_count >= INTG_FAIL_SLOTS) return;  /* full: degrade gracefully */
    intg_fail_tab[intg_fail_count].hf = hf;
    intg_fail_tab[intg_fail_count].hx = hx;
    intg_fail_tab[intg_fail_count].method = method;
    intg_fail_count++;
}

Expr* builtin_integrate(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    Expr* f = res->data.function.args[0];

    /* Thread over a list integrand: Integrate[{f1, ..., fn}, spec...] gives
     * {Integrate[f1, spec...], ..., Integrate[fn, spec...]}.  Integrate is
     * deliberately NOT Listable (that would also thread over a `{x, a, b}`
     * range spec), so the integrand-only threading is handled explicitly here.
     * We build the List of sub-Integrate calls and let the evaluator recurse
     * into each element (which stays symbolic if a piece cannot be done). */
    if (f->type == EXPR_FUNCTION &&
        f->data.function.head->type == EXPR_SYMBOL &&
        f->data.function.head->data.symbol.name == SYM_List) {
        size_t n = f->data.function.arg_count;
        Expr** elems = malloc(n * sizeof(Expr*));
        for (size_t i = 0; i < n; i++) {
            Expr** sub = malloc(argc * sizeof(Expr*));
            sub[0] = expr_copy(f->data.function.args[i]);
            for (size_t j = 1; j < argc; j++)
                sub[j] = expr_copy(res->data.function.args[j]);
            elems[i] = expr_new_function(expr_new_symbol(SYM_Integrate), sub, argc);
            free(sub);
        }
        Expr* out = expr_new_function(expr_new_symbol(SYM_List), elems, n);
        free(elems);
        return out;
    }

    /* Definite form: the second argument is a `{x, a, b}` range spec or a
     * `{x, z0, ..., zn}` complex line/contour spec.  Handles the iterated
     * multi-spec form too. */
    if (integrate_line_is_contour_spec(res->data.function.args[1]))
        return integrate_definite_scoped(res);

    /* Trailing options are parsed as a loop below, so more than one is allowed
     * (`Integrate[f, x, Method -> "...", TimeConstraint -> 3]`).  The loop
     * rejects anything that is neither a recognised option nor a method name. */

    Expr* x = res->data.function.args[1];

    /* The integration variable must be a (single) symbol so we can
     * compute derivatives, partial fractions, ... in it. */
    if (x->type != EXPR_SYMBOL) return NULL;

    /* Applied InterpolatingFunction objects integrate to a fresh antiderivative
     * InterpolatingFunction, mirroring how D differentiates them.  This is
     * handled before the inexact-rationalisation scan below: the object embeds
     * Real sample data that must NOT be force-rationalised. */
    {
        Expr* interp = integrate_interp(f, x);
        if (interp) return interp;
    }

    /* SeriesData integrates term-by-term to a new SeriesData. Handle this
     * before the inexact-rationalisation scan and method cascade, none of
     * which understand the SeriesData head. */
    if (f->type == EXPR_FUNCTION &&
        f->data.function.head->type == EXPR_SYMBOL &&
        f->data.function.head->data.symbol.name == SYM_SeriesData &&
        f->data.function.arg_count == 6) {
        Expr* r = series_integrate(f, x);
        if (r) return r;   /* residue/unsupported -> fall through, stays unevaluated */
    }

    /* Parse the optional Method -> "..." option.  `method_sub` captures the
     * user-pinned substitution from the
     * `{"DerivativeDivides", "Substitution" -> u}` list form (NULL otherwise);
     * it is freed below alongside `coerced`. */
    IntegrateMethod method = METHOD_AUTOMATIC;
    Expr* method_sub = NULL;
    /* An `Assumptions -> ...` option is accepted on the indefinite form for
     * surface compatibility (it constrains parameter domains); the indefinite
     * cascade does not yet consume it, so it is simply skipped here rather than
     * mis-parsed as a Method value. */
    /* TimeConstraint -> t bounds the kernel search of the one stage that is, by
     * measurement, nearly the whole cost of a cascade that fails to close (see
     * integrate_derivdivides.h).  Default 3 s; Infinity restores the unbounded
     * search.  A trailing option may appear in any order and any combination,
     * which is why this is a loop rather than the single-option test it replaced. */
    double time_constraint = INTEGRATE_DEFAULT_TIME_CONSTRAINT;
    for (size_t t = 2; t < argc; t++) {
        Expr* opt = res->data.function.args[t];
        if (option_lhs_is(opt, SYM_Assumptions)) continue;   /* accepted, not yet consumed */
        if (integrate_option_is(opt, "TimeConstraint")) {
            time_constraint = integrate_parse_time_constraint(opt->data.function.args[1]);
            continue;
        }
        method = parse_method_option(opt, &method_sub);
        if (method == METHOD_INVALID) {
            static uint64_t last_warned_hash = 0;
            uint64_t h = expr_hash(res);
            if (h != last_warned_hash) {
                mth_message("Integrate", "method", "Method option value is not one of \"Automatic\", \"BronsteinRational\", \"DerivativeDivides\", \"LinearRadicals\", \"QuadraticRadicals\", \"LinearRatioRadicals\", \"ChebychevAlgebraic\", \"GoursatAlgebraic\", \"Weierstrass\", \"RischTranscendental\", \"CRCTable\", \"ParallelMixedTower\", \"ParallelMixedSpecial\", \"Trager\", \"TragerKauers\", \"NewtonLeibniz\", \"LineIntegral\".");
                last_warned_hash = h;
            }
            return NULL;
        }
    }

    /* Fail-memo short-circuit: if the deterministic cascade already returned
     * NULL for this exact (integrand, variable, method) earlier in THIS
     * top-level evaluation, it will do so again -- skip the redundant second
     * search and its duplicate diagnostics.  Keyed on the incoming integrand
     * `f`: both fixed-point passes see the same Orderless-canonicalised arg, so
     * their hashes match even though the eval clock churned between them.  See
     * the intg_fail_* declaration above. */
    uint64_t intg_hf = expr_hash(f);
    uint64_t intg_hx = expr_hash(x);
    intg_fail_sync_epoch(eval_toplevel_id());
    if (intg_fail_seen(intg_hf, intg_hx, (int)method)) {
        if (method_sub) expr_free(method_sub);
        return NULL;
    }

    /* Inexact integrands: route through the shared preprocessor in
     * common.c -- if `f` contains any inexact leaf, force-rationalise
     * it (with bit-exact ½-ulp fallback so transcendental floats like
     * N[Pi] still become rationals), integrate exactly, then
     * numericalise the result so the user observes inexact-in /
     * inexact-out semantics.  Same contract as Solve.
     *
     * The scan also captures the *minimum* precision (in bits) across
     * every inexact leaf, which is then used both as the rationalisation
     * tolerance and as the precision of the final numericalised result. */
    CommonInexactInfo inexact = common_scan_inexact(f);
    Expr* coerced = NULL;
    Expr* effective_f = f;
    /* For an inexact integrand the transcendental-Risch a^u debasing must stay OFF:
     * the float base is rationalised to a large exact fraction below, and debasing
     * it would fabricate an exact Log-based closed form that numericalises to a
     * degraded exponent (2.71828^x -> 2.71828^(0.999999 x)).  Inexact-in stays
     * inexact-out / unevaluated instead.  Restored before returning. */
    bool debase_saved = true;
    if (inexact.has_inexact) {
        debase_saved = rt_transcendental_set_debase(false);
        coerced = common_rationalize_input(f, inexact.min_bits);
        if (!coerced) { rt_transcendental_set_debase(debase_saved); return NULL; }
        effective_f = coerced;
    }

    /* Enter the cascade one level deeper.  A sub-integral spawned by an internal
     * substitution (DerivativeDivides, ...) re-enters here at depth >= 2, which
     * is how the RischTranscendental stage tells the user's integrand from an
     * internal recursion variable when deciding whether to name it in a
     * nonelem diagnostic. */
    /* Arm the TimeConstraint for this descent, at the OUTERMOST frame only: the
     * deadline is absolute, so every nested re-entry from a stage's own recursion
     * inherits the same cutoff -- which is the point, since that fan-out is what
     * is being bounded.  An explicit `Method -> ...` is a deliberate request for
     * the full search and is left unbounded unless the caller passed
     * TimeConstraint themselves. */
    if (g_integrate_depth == 0 && method == METHOD_AUTOMATIC)
        integrate_arm_budget(time_constraint);

    g_integrate_depth++;
    /* Fresh top-level cascade: clear the Integrate::nonelem de-dup flag so the
     * first stage that proves non-elementarity (RischTranscendental or, after
     * it, ParallelMixedTower) speaks exactly once.  A sub-integral re-enters at
     * depth >= 2 and must NOT clear it. */
    if (g_integrate_depth == 1) {
        g_integrate_nonelem_announced = 0;
        intg_nonelem_pending_clear();   /* drop anything a longjmp stranded */
    }

    Expr* result = NULL;
    switch (method) {
        case METHOD_AUTOMATIC:
            /* Depth 1 only: a NESTED pass (a stage recursing into Integrate --
             * the cascade switch runs at every depth) would otherwise flush the
             * outermost stage still in flight and charge it only the few ms
             * before the recursion began, losing the rest of its cost. */
            if (ip_enabled() && g_integrate_depth == 1) {
                ip_flush_inflight(); g_ip_entries++;
            }
            IP_STAGE(IP_UNDEFINED, try_undefined(effective_f, x));
            /* A step factor with a linear argument, g(x) UnitStep[x - a]: the
             * continuous antiderivative (G(x) - G(a)) UnitStep[x - a].  Cheap
             * structural gate (a top-level step factor), and first because no
             * later stage integrates a step at all -- they would only grind. */
            IP_STAGE(IP_STEP,         integrate_step_indefinite_try(effective_f, x));
            IP_STAGE(IP_RATIONAL,     try_rational(effective_f, x));
            IP_STAGE(IP_LINRAD,       try_linrad(effective_f, x));
            IP_STAGE(IP_QUADRAD,      try_quadrad(effective_f, x));
            IP_STAGE(IP_LINRATIORAD,  try_linratiorad(effective_f, x));
            /* Chebychev binomial differentials: a fast, deterministic
             * rationalising substitution that closes (correct by construction),
             * so it runs ahead of the Eliminate/Solve search and Risch-Norman. */
            IP_STAGE(IP_CHEBYCHEV,    try_chebychev(effective_f, x));
            /* Fresnel: K Sin/Cos of a quadratic -> FresnelS/FresnelC by completing
             * the square (the trig sibling of the Gaussian -> Erf recognizer),
             * deterministic and diff-back verified. */
            IP_STAGE(IP_FRESNEL,      try_fresnel(effective_f, x));
            /* Symbolic power times an exponential -> incomplete Gamma.  Before the
             * general stages because they SEARCH on this shape (12.9 s to decline on
             * x^n E^(-x)); gated to a symbolic exponent, so nothing with a numeric
             * one is diverted from its elementary / Erf answer below. */
            IP_STAGE(IP_GAMMAPOWER,   try_gammapower(effective_f, x));
            /* Weierstrass before derivative-divides: it is a domain-specific,
             * deterministic algorithm for rational trig/hyperbolic integrands
             * that is guaranteed to close (and verified by construction), so it
             * runs ahead of the more expensive Eliminate/Solve substitution
             * search and ahead of Risch-Norman's complex-logarithm forms. */
            /* Linearity before Weierstrass: term-by-term over a Plus so trig
             * polynomials (Sin[2x]+Sin[3x], constant-coeff forcing sums) get the
             * clean additive antiderivative instead of a divergent/looping
             * tan-half-angle form.  Declines unless every term closes. */
            IP_STAGE(IP_LINEARITY,    try_linearity(effective_f, x));
            /* Linear-argument substitution before Weierstrass: pull a shared
             * non-trivial linear argument a*x+b out of a rational trig integrand
             * (u = a*x+b) so the bare-argument form closes cleanly, instead of
             * Weierstrass exploding Sec[3x]^2 into a degree-12 Tan[x/2] rational.
             * Declines unless every x sits in a kernel of one such argument with
             * a kernel in a denominator, so other integrands are untouched. */
            IP_STAGE(IP_LINARG,       try_linarg(effective_f, x));
            IP_STAGE(IP_WEIERSTRASS,  try_weierstrass(effective_f, x));
            IP_STAGE(IP_DERIVDIVIDES, try_derivdivides(effective_f, x));
            /* Recursive transcendental Risch: correct by construction, adding
             * closed forms the earlier stages missed (logarithmic polynomials,
             * Gaussians -> Erf, exp/x -> ExpIntegralEi, 1/Log -> LogIntegral,
             * Log[1+ax]/x -> PolyLog).  The former RischNorman / RischNormanBlake
             * pmint heuristics that ran here were removed in v0.163 -- the
             * ParallelMixedTower stage below subsumes both. */
            IP_STAGE(IP_RISCHTRANS,   try_rischtranscendental(effective_f, x));
            IP_STAGE(IP_CRCTABLE,     try_crctable(effective_f, x));
            /* Log[g] times a trig/hyperbolic kernel -> one integration by parts.
             * After the cheap elementary stages (the `if (!result)` short-circuit
             * means it sees only integrands those declined) and before the
             * grinding tower / special-function tail, so the Log*trig family --
             * which otherwise costs ParallelMixedSpecial 3.8-12 s -- closes to a
             * Log*trig + SinIntegral/CosIntegral form.  Diff-back verified;
             * declines (cascade continues) on any mis-match. */
            /* Polylogarithm ladders: R(x) Log[a x + b]^k and polynomial x trig
             * rationals -> PolyLog closed forms.  Before LogByParts and the
             * grinding tower tail: cheap structural gate, exact construction,
             * diff-back accepted. */
            IP_STAGE(IP_POLYLOG,      try_polylogladder(effective_f, x));
            /* Bessel / Airy kernel products: one linear system in the kernels'
             * differential module, before the tower tail grinds on them. */
            IP_STAGE(IP_HOLONOMIC,    try_holonomic(effective_f, x));
            /* ...and the hypergeometric / MeijerG closed form when no
             * Bessel-product antiderivative exists. */
            IP_STAGE(IP_BESSELHYPER,  try_besselhyper(effective_f, x));
            IP_STAGE(IP_LOGBYPARTS,   try_logbyparts(effective_f, x));
            /* The parallel Risch-Norman integrator over a mixed tower.  After
             * the CRC table so a tabled integral is answered cleanly and
             * cheaply rather than paying this stage's package load and search
             * on every case the table already covers.  Skipped for a pseudo-
             * elliptic F/R^p integrand (Goursat's exclusive domain): PMT's
             * tower search over that genus>0 algebraic curve only grinds --
             * tens of seconds before declining (e.g. 45s on
             * (t^4+2t^3-4)/(t^2 Sqrt[(t^2-1)(t^2-4)])) -- and closes nothing,
             * so the integrand goes straight to the Goursat stage that does.
             * The explicit Method -> "ParallelMixedTower" is unaffected.
             * Also skipped when g_integrate_no_special is raised (a speculative
             * sub-integral of LogByParts): a reduction that does not close at a
             * cheap stage must decline promptly, not grind the tower. */
            if (!g_integrate_no_special && !has_pseudoelliptic_radical(effective_f, x))
                IP_STAGE(IP_PMT,      try_parallelmixedtower(effective_f, x, true));
            /* Goursat pseudo-elliptic (and cube-/fourth-root) reductions to
             * genus-0 curves.  Runs last, after ParallelMixedTower: it is a
             * specialist for pseudo-elliptic F/R^p (p in {1/2,1/3,2/3,1/4,3/4})
             * that the general methods above do not close, so every earlier
             * stage gets first crack at any integrand it also handles.  (The
             * grind-prone search stages -- linearity's product-split,
             * derivative-divides' Eliminate/Solve, and PMT's tower search --
             * are gated off above for exactly this F/R^p shape, on which they
             * only spin before declining, so such an integrand reaches this
             * stage promptly.)  Still correct-by-construction (an internal
             * differentiate-back guard), so a decline here leaves the integral
             * unevaluated rather than wrong.  Skipped under g_integrate_no_special
             * (a LogByParts speculative sub-integral). */
            if (!g_integrate_no_special)
                IP_STAGE(IP_GOURSAT,  try_goursat(effective_f, x));
            /* ParallelMixedSpecial: the special-function stage, and the LAST
             * stage of the cascade.
             *
             * Last for two reasons. It is the only stage that may answer with a
             * non-elementary function (ExpIntegralEi, Erf, the elliptic family,
             * ...), and an ELEMENTARY antiderivative is always the better answer
             * -- so every stage that could produce one, Goursat's pseudo-elliptic
             * reductions included, gets first crack. And it is a general search,
             * so it costs least where it is reached only by integrands nothing
             * cheaper could close.
             *
             * Note it is NOT gated on has_pseudoelliptic_radical, unlike
             * ParallelMixedTower just above: that gate keeps the elementary-only
             * stage off a genus>0 F/R^p curve, where it grinds and closes
             * nothing, and an elliptic pencil is precisely what THIS stage
             * closes. try_parallelmixedspecial keeps the rational-structure skip.
             *
             * Complete answers only: the cascade calls the surface that returns
             * an antiderivative or declines, never the partial mode, so plain
             * Integrate[f, x] can no more return a half-solved
             * `answer + Inactive[Integrate][remainder, x]` than it could before.
             * Skipped under g_integrate_no_special (a LogByParts speculative
             * sub-integral): the whole point of the flag is to not pay this
             * search twice. */
            if (!g_integrate_no_special)
                IP_STAGE(IP_PMS,      try_parallelmixedspecial(effective_f, x, true));
            break;
        case METHOD_RATIONAL:
            result = try_rational(effective_f, x);
            break;
        case METHOD_DERIVATIVE_DIVIDES:
            /* Explicit method runs the thorough Eliminate/Solve search in
             * addition to the direct quotient (the Automatic cascade uses the
             * direct strategy only).  A user-pinned `"Substitution" -> u`
             * restricts the search to that single kernel. */
            result = method_sub
                   ? integrate_derivdivides_with_sub(effective_f, x, method_sub)
                   : integrate_derivdivides_full(effective_f, x);
            break;
        case METHOD_LINEAR_RADICALS:
            result = try_linrad(effective_f, x);
            break;
        case METHOD_QUADRATIC_RADICALS:
            result = try_quadrad(effective_f, x);
            break;
        case METHOD_LINEAR_RATIO_RADICALS:
            result = try_linratiorad(effective_f, x);
            break;
        case METHOD_CHEBYCHEV:
            result = try_chebychev(effective_f, x);
            break;
        case METHOD_GOURSAT:
            result = try_goursat(effective_f, x);
            break;
        case METHOD_WEIERSTRASS:
            result = integrate_jeffrey_full(effective_f, x);
            break;
        case METHOD_RISCH_TRANSCENDENTAL:
            result = try_rischtranscendental(effective_f, x);
            break;
        case METHOD_CRCTABLE:
            result = try_crctable(effective_f, x);
            break;
        case METHOD_POLYLOG_LADDER:
            result = try_polylogladder(effective_f, x);
            break;
        case METHOD_HOLONOMIC:
            result = try_holonomic(effective_f, x);
            break;
        case METHOD_BESSEL_HYPERGEOMETRIC:
            result = try_besselhyper(effective_f, x);
            break;
        case METHOD_PARALLEL_MIXED_TOWER:
            /* Explicit request: no internal time budgets (builtin_integrate_pmt). */
            result = try_parallelmixedtower(effective_f, x, false);
            break;
        case METHOD_PARALLEL_MIXED_SPECIAL:
            /* Explicit request: the package's full $SpecialTimeBudget, not the
             * cascade's tight one (see PMS_CASCADE_BUDGET_SECONDS). */
            result = try_parallelmixedspecial(effective_f, x, false);
            break;
        case METHOD_TRAGER:
            result = try_trager(effective_f, x, "Integrate`Trager");
            break;
        case METHOD_TRAGER_KAUERS:
            result = try_trager(effective_f, x, "Integrate`TragerKauers");
            break;
        case METHOD_UNDEFINED:
            result = try_undefined(effective_f, x);
            break;
        case METHOD_NEWTON_LEIBNIZ:
        case METHOD_LINE_INTEGRAL:
        case METHOD_RESIDUE:
        case METHOD_DIFF_UNDER_INT:
        case METHOD_RAMANUJAN:
        case METHOD_SYMMETRY:
        case METHOD_BETA:
        case METHOD_TRIG_POWER:
        case METHOD_SINPOW_MONO:
        case METHOD_OSC_POWER:
        case METHOD_RATIONAL_LOG:
        case METHOD_INTEGRAL_REP:
        case METHOD_MEIJER_CONV:
        case METHOD_PIECEWISE:
        case METHOD_SERIES:
        case METHOD_KING:
        case METHOD_FEYNMAN_PARAM:
            /* Definite-only mechanisms.  Meaningless on the indefinite form
             * Integrate[f, x, Method -> "NewtonLeibniz" / "LineIntegral" /
             * "Residue" / "DiffUnderInt" / "RamanujanMasterTheorem" /
             * "IntegralRepresentation"]; leave unevaluated. */
            break;
        case METHOD_INVALID:
            break;  /* unreachable: handled above */
    }

    /* Single post-method normalisation chokepoint: flatten x-free
     * reciprocal/product powers across every method's output, and — for
     * radical-bearing antiderivatives — recombine the algebraic part and
     * tidy inverse-trig arguments.  Runs on the symbolic result before the
     * inexact path numericalises it.  See intsimp_finalize. */
    if (result) result = intsimp_finalize(result, x);   /* takes ownership */

    if (coerced) expr_free(coerced);
    if (method_sub) expr_free(method_sub);

    if (inexact.has_inexact && result) {
        Expr* numeric = common_numericalize_result(result, inexact.min_bits);
        expr_free(result);
        result = numeric;
    }

    rt_transcendental_set_debase(debase_saved);   /* no-op unless the inexact path disabled it */

    g_integrate_depth--;

    /* The cascade is over: a stage that proved non-elementarity gets to speak
     * only if no LATER stage produced a closed form.  ParallelMixedSpecial now
     * answers a good many integrands RischTranscendental had already proved
     * non-elementary, and warning about an answer we are about to print is
     * noise. */
    if (g_integrate_depth == 0) {
        if (result) intg_nonelem_pending_clear(); else intg_nonelem_pending_flush();
    }

    /* Record a failure so the fixed-point loop's redundant re-entry on the
     * identical integrand (same top-level evaluation) short-circuits above.
     * Only the outermost user call arms the memo: an inner sub-integral's
     * failure is expected and must not suppress an unrelated outer retry. */
    if (!result && g_integrate_depth == 0)
        intg_fail_record(intg_hf, intg_hx, (int)method);

    return result;
}

void integrate_init(void) {
    symtab_add_builtin("Integrate", builtin_integrate);
    /* A TimeConstrained timeout that siglongjmps out of a definite integral
     * must not leave the definite depth raised, or every later Integrate::nonelem
     * would be deferred forever. */
    tc_register_unwind_counter(&g_integrate_definite_depth);

    /* The ParallelMixedTower method symbol (lazy-loads src/internal/mixed/
     * ParallelMixed.m on first use).  Registered here so the qualified-symbol
     * surface Integrate`ParallelMixedTower[f, x] works from a cold session. */
    symtab_add_builtin("Integrate`ParallelMixedTower", builtin_integrate_pmt);
    symtab_get_def("Integrate`ParallelMixedTower")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Integrate`ParallelMixedTower",
        "Integrate`ParallelMixedTower[f, x] integrates f with respect to x by the "
        "parallel (Risch-Norman) method over a simple radical in a mixed "
        "transcendental tower (S. Blake, Parallel Integration over Simple Radical "
        "Extensions II). It returns the antiderivative, or $Failed / a certificate "
        "list that Integrate's cascade treats as a decline. When the method PROVES "
        "the integrand has no elementary antiderivative it issues Integrate::nonelem, "
        "as RischTranscendental does. Called directly (or via Method -> "
        "\"ParallelMixedTower\") it has no internal time budget; the Automatic "
        "cascade bounds it by $ParallelMixedTimeBudget.");
    /* The Trager / TragerKauers method symbols (lazy-load
     * src/internal/trager/Trager.m). */
    symtab_add_builtin("Integrate`Trager", builtin_integrate_trager);
    symtab_get_def("Integrate`Trager")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Integrate`Trager",
        "Integrate`Trager[f, x] integrates f, a rational function of x and powers of "
        "a single radical r(x)^(1/n) with rational coefficients, by Trager's algorithm "
        "for simple radical extensions (Trager 1984): Hermite reduction on an integral "
        "basis, Rothstein-Trager residues, residue divisors as ideals, principal "
        "generators by normalisation at infinity, and torsion orders by reduction "
        "modulo two primes. It returns the antiderivative (verified exactly), "
        "{\"not elementary\", reason} -- a proof, also announced as Integrate::nonelem "
        "-- or {\"failed\", reason}.");
    symtab_add_builtin("Integrate`TragerKauers", builtin_integrate_trager_kauers);
    symtab_get_def("Integrate`TragerKauers")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Integrate`TragerKauers",
        "Integrate`TragerKauers[f, x] is Integrate`Trager with the logarithmic part "
        "found by Kauers' Groebner-basis heuristic (Kauers 2008): for each irreducible "
        "factor of the residue polynomial, the least power of the residue ideal that is "
        "principal modulo the curve gives the logand. Faster on many inputs but "
        "incomplete: it returns {\"failed\", reason} where the heuristic does not "
        "close. Non-elementarity certificates come only from the shared Hermite and "
        "residue stages.");
    /* The ParallelMixedSpecial method symbol (lazy-loads
     * src/internal/ParallelMixedSpecial.m, which loads Part II itself). */
    symtab_add_builtin("Integrate`ParallelMixedSpecial", builtin_integrate_pms);
    symtab_get_def("Integrate`ParallelMixedSpecial")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Integrate`ParallelMixedSpecial",
        "Integrate`ParallelMixedSpecial[f, x] integrates f with respect to x by the "
        "SPECIAL-FUNCTION stage of the parallel (Risch-Norman) method over a mixed "
        "tower (S. Blake): the ansatz of Integrate`ParallelMixedTower extended with "
        "special-function kernels, so the answer may contain ExpIntegralEi, "
        "LogIntegral, Erf, Erfc, Erfi, the incomplete Gamma[s, z], SinIntegral, "
        "CosIntegral, or the elliptic EllipticF, EllipticE and EllipticPi of a "
        "pencil. It returns the antiderivative, or a certificate / failure list "
        "that Integrate's cascade treats as a decline. When the method PROVES the "
        "integrand has no elementary antiderivative it issues Integrate::nonelem.");

    /* NOT Listable: Listable would thread over a definite-integral range
     * spec `{x, a, b}` element-wise (producing garbage like
     * {Integrate[f,x], Integrate[f,a], Integrate[f,b]}).  The `{x,a,b}` form
     * is recognised explicitly at the top of builtin_integrate and routed to
     * the definite (Newton-Leibniz) path. */
    symtab_get_def("Integrate")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Integrate",
        "Integrate[f, x] gives the indefinite integral of f with respect to x.\n"
        "Integrate[f, {x, xmin, xmax}] gives the definite integral by the\n"
        "fundamental theorem of calculus (Method -> \"NewtonLeibniz\").\n"
        "Integrate[f, {x, xmin, xmax}, {y, ymin, ymax}, ...] gives the iterated\n"
        "multiple integral (innermost/last spec integrated first; inner bounds\n"
        "may depend on outer variables).  See also Integrate`SingularPoints.\n"
        "Integrate[f, x, Method -> \"<name>\"] dispatches directly to a single\n"
        "subroutine, bypassing the default cascade.  Accepted method names:\n"
        "  \"Automatic\"          — the full dispatch cascade (default)\n"
        "  \"BronsteinRational\"  — Integrate`BronsteinRational (polynomial / rational)\n"
        "  \"DerivativeDivides\"  — Integrate`DerivativeDivides (substitution u(x); direct + Eliminate/Solve)\n"
        "  \"LinearRadicals\"     — Integrate`LinearRadicals (rationalise radicals of a x + b)\n"
        "  \"QuadraticRadicals\"  — Integrate`QuadraticRadicals (Euler substitution for Sqrt[a x^2 + b x + c])\n"
        "  \"LinearRatioRadicals\" — Integrate`LinearRatioRadicals (rationalise radicals of (a x + b)/(c x + d))\n"
        "  \"ChebychevAlgebraic\" — Integrate`ChebychevAlgebraic (binomial x^p (a x^r + b)^q via Chebychev's theorem)\n"
        "  \"GoursatAlgebraic\"   — Integrate`GoursatAlgebraic (pseudo-elliptic F/R^p, p in {1/2,1/3,2/3,1/4,3/4}, via Mobius eigendescent)\n"
        "  \"Weierstrass\"        — Integrate`Weierstrass (continuous tan(x/2) / tanh(x/2) substitution)\n"
        "  \"RischTranscendental\"       — Integrate`RischTranscendental (recursive transcendental Risch; correct by construction)\n"
        "  \"CRCTable\"           — Integrate`CRCTable (lazy-loaded CRC integral table)\n"
        "  \"ParallelMixedTower\" — Integrate`ParallelMixedTower (parallel Risch-Norman over a simple radical in a mixed transcendental tower; Blake II)\n"
        "  \"ParallelMixedSpecial\" — Integrate`ParallelMixedSpecial (the same over a mixed tower with special-function kernels: Ei, li, erf, incomplete Gamma, Si, Ci, elliptic F/E/Pi)\n"
        "  \"Trager\"             — Integrate`Trager (Trager 1984: complete algorithm for one radical r(x)^(1/n) over Q;\n"
        "                          decides elementarity, torsion by reduction mod p; not in the Automatic cascade)\n"
        "  \"TragerKauers\"       — Integrate`TragerKauers (the same with Kauers' 2008 Groebner heuristic for the log part)\n"
        "  \"Undefined\"          — Integrate`Undefined (unknown functions u[x], u'[x]; Roach §1.7)\n"
        "  \"NewtonLeibniz\"       — real definite integrals via F(b)-F(a) (implicit for the {x,a,b} form)\n"
        "  \"LineIntegral\"        — complex contour integrals (implicit for the {x,z0,...,zn} form)\n"
        "  \"Residue\"             — improper/periodic real definite integrals by the residue theorem\n"
        "                          (rational/Fourier on (-Inf,Inf), rational-in-Sin/Cos over a period,\n"
        "                          principal values, even half-lines); tried before NewtonLeibniz under Automatic\n"
        "  \"DiffUnderInt\"         — parameter-dependent definite integrals by differentiation under the\n"
        "  (\"DifferentiationUnderIntegral\") integral sign (Feynman's trick): Integrate`DiffUnderInt;\n"
        "                          Laplace/Fourier, sinc, and even-rational half-line families;\n"
        "                          tried after Residue and NewtonLeibniz in the definite cascade\n"
        "  \"RamanujanMasterTheorem\" — half-line Int_0^Inf x^(s-1) f(x) dx by the Mellin transform /\n"
        "  (\"Mellin\")              Ramanujan Master Theorem: Integrate`RamanujanMasterTheorem;\n"
        "                          exp/Gaussian/algebraic/Cos/Sin/ArcTan/Log/BesselJ/pFq/PolyLog\n"
        "                          kernels (monomial x^k substitution; Erf, incomplete Gamma, BesselJ^2\n"
        "                          reduced to pFq); also the exp-geometric kernel 1/(E^(c x)+g)\n"
        "                          (Bose-Einstein / Fermi-Dirac -> Gamma*PolyLog), a Frullani pre-pass\n"
        "                          (f(a x)-f(b x))/x -> (f(0)-f(Inf)) Log[b/a], and a Log[x]^k weight;\n"
        "                          strip-gated, yielding a ConditionalExpression when\n"
        "                          Assumptions do not prove convergence; after NewtonLeibniz under Automatic\n"
        "Method -> {\"DerivativeDivides\", \"Substitution\" -> u} pins the kernel u(x),\n"
        "trialing only that substitution.\n"
        "Named methods are strict: failure returns unevaluated, with no fallback.\n"
        "Integrate[f, x, TimeConstraint -> t] gives the cascade a wall-clock budget\n"
        "of t seconds (default 3; Infinity for no budget).  It is checked BETWEEN\n"
        "stages, at every recursion depth, and inside the substitution stage's\n"
        "kernel and branch loops and its Eliminate/Solve calls -- so it stops the\n"
        "cascade progressing once spent, but it is not yet a hard bound: a single\n"
        "long-running evaluation inside one stage can overrun it.  It matters only\n"
        "for integrands Integrate cannot close; measured from inside the cascade\n"
        "(MATHILDA_INTEGRATE_PROFILE=1), two stages account for 90.8% of that cost\n"
        "-- the linearity split (43.5%) and the substitution search (47.3%) -- while\n"
        "the eight cheap exact stages ahead of them are together under 0.1%.  An\n"
        "explicit Method -> name is a deliberate request for the full search and\n"
        "carries no budget unless you pass TimeConstraint yourself.\n"
        "The CRCTable rules are loaded from disk on first use only.\n"
        "An applied 1-D InterpolatingFunction integrates to its antiderivative\n"
        "InterpolatingFunction (mirroring D).");

    /* Initialise the Integrate` package: HermiteReduce, IntegratePolynomial,
     * helpers, and the explicit `Integrate`BronsteinRational` entry. */
    intrat_init();

    /* Cherry 1986 Thm 4.4 Sigma-decomposition: Integrate`SigmaDecomposition
     * (also the engine behind the logarithmic-integral non-existence decision). */
    void cherry_sigma_init(void);
    cherry_sigma_init();

    /* Direct debuggable surfaces for every Cherry engine: Integrate`Cherry`Ei /
     * ExpMultiterm / Li / Dilog / DilogExp. */
    void cherry_builtins_init(void);
    cherry_builtins_init();

    /* Undefined-function integrator (Roach §1.7): Integrate`Undefined. */
    integrate_unknown_init();

    /* Derivative-divides substitution: Integrate`DerivativeDivides. */
    integrate_derivdivides_init();

    /* Linear-radical rationalising substitution: Integrate`LinearRadicals. */
    integrate_linrad_init();

    /* Quadratic-radical (Euler) substitution: Integrate`QuadraticRadicals. */
    integrate_quadrad_init();

    /* Linear-ratio-radical (Möbius) substitution: Integrate`LinearRatioRadicals. */
    integrate_linratiorad_init();

    /* Chebychev binomial differential: Integrate`ChebychevAlgebraic. */
    integrate_chebychev_init();

    /* Goursat pseudo-elliptic / cube-/fourth-root: Integrate`GoursatAlgebraic. */
    integrate_goursat_init();

    /* Continuous Weierstrass substitution (Jeffrey & Rich 1994):
     * Integrate`Weierstrass. */
    integrate_jeffrey_init();

    /* Definite integration by the fundamental theorem of calculus:
     * Integrate`NewtonLeibniz and the pole detector Integrate`SingularPoints. */
    integrate_newton_leibniz_init();

    /* Origin-symmetry reduction: Integrate`Symmetry (odd -> 0, even -> 2 half).
     * Tried first in the definite dispatch. */
    integrate_symmetry_init();

    /* Piecewise integrands: Integrate`Piecewise (Abs/Sign/HeavisideTheta/
     * UnitStep/Min/Max/Piecewise); split at interior breakpoints. */
    integrate_piecewise_init();

    /* Euler-Beta reductions: Integrate`Beta ([0,1]) and Integrate`TrigPower
     * (Sin^m Cos^n over a canonical trig interval). */
    integrate_beta_init();

    /* Complex line / contour integration: Integrate`LineIntegral and the
     * on-path singularity detector Integrate`PathSingularPoints. */
    integrate_line_init();

    /* Definite integration by the residue theorem: Integrate`ContourResidue.
     * Engaged before Newton-Leibniz for a single real spec under Automatic. */
    integrate_residue_init();

    /* Definite integration by differentiation under the integral sign
     * (Feynman's trick): Integrate`DiffUnderInt.  Last resort in the definite
     * cascade, after residue and Newton-Leibniz. */
    integrate_diffunderint_init();
    integrate_ramanujan_init();
    integrate_intrep_init();
    integrate_series_init();
    integrate_king_init();
    integrate_feynparam_init();

    /* Recursive transcendental Risch integrator:
     * Integrate`RischTranscendental.  Correct by construction. */
    integrate_risch_transcendental_init();

    /* Bronstein differential-field foundation (Symbolic Integration I, Ch. 3):
     * splitting factorization / canonical representation primitives, exposed as
     * the Risch` internal builtins and reused by the recursive integrator. */
    integrate_risch_canonical_init();

    /* Bronstein structure theorems (Symbolic Integration I, §9.3): the Q-linear
     * decision of whether a proposed new exp/log monomial is algebraically
     * independent of the tower or reducible. */
    integrate_risch_structure_init();

    /* Bronstein hypertangent case (Symbolic Integration I, §5.10): direct real
     * integration over a tangent monomial (polynomial part). */
    integrate_risch_hypertangent_init();

    /* Bronstein coupled differential systems (Symbolic Integration I, Ch. 8) +
     * the hypertangent reduced case (§5.10): CoupledDESystem via the reduction
     * to a Risch DE over C(i)(x), and pole-peeling over a tangent monomial. */
    integrate_risch_coupled_init();

    /* Bronstein Hermite reduction (Symbolic Integration I, §5.3): the literal
     * quadratic HermiteReduce over the tower derivation, f = D[g] + h + r. */
    integrate_risch_hermite_init();

    /* Pre-register Integrate`CRCTable so ?Integrate`CRCTable shows a
     * docstring even before the lazy load fires.  The actual rule
     * definition is added when CRCMathTablesIntegrals.m is Get-loaded
     * on first call — that DownValue install must NOT be blocked by
     * Protected, so we deliberately leave attributes empty here.  The
     * lazy load itself promotes Integrate`CRCTable to Protected via
     * SetAttributes inside CRCMathTablesIntegrals.m. */
    symtab_get_def("Integrate`CRCTable");
    symtab_set_docstring("Integrate`CRCTable",
        "Integrate`CRCTable[f, x] looks up f in the CRC Standard Mathematical\n"
        "Tables (31st ed.) integral table. Returns the antiderivative on a hit,\n"
        "or the unevaluated form if no rule applies. The table file\n"
        "(src/internal/CRCMathTablesIntegrals.m) is loaded from disk on first\n"
        "use only.");
}
