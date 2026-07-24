/* intrat.c — rational-function integration package (`Integrate``).
 *
 * Phase 1 of plans/INTEGRATE_PLAN.md.  Implements the scaffolding,
 * polynomial integration, Mack's linear Hermite reduction, the
 * derivative-recognition fast path (`c*D'/D^k`), and the four
 * one-line helpers (Content / Primitive / Monic / LeadingCoefficient)
 * used by the later phases.  Everything is exposed under the
 * `Integrate`` context so each stage is REPL-testable in isolation.
 *
 * Algorithm references:
 *   - Bronstein, *Symbolic Integration I*, Springer 2nd ed., 2004.
 *   - Trager, *Algebraic Factoring and Rational Function
 *     Integration*, ACM SYMSAC 1976.
 * The Mathematica baseline being ported is `IntegrateRational.m` in
 * the repo root — line numbers in the comments below refer to it.
 *
 * Memory contract: every public BuiltinFunc follows the Mathilda
 * convention — caller (evaluator) owns `res`; the function returns a
 * freshly-allocated Expr* on success or NULL on failure.
 */

#include "intrat.h"
#include "intrat_internal.h"
#include "intsimp.h"
#include "common.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "internal.h"
#include "poly.h"
#include "expand.h"
#include "sym_names.h"
#include "print.h"
#include "options.h"
#include "rationalize.h"
#include "root.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

/* ------------------------------------------------------------------ */
/* Trace scaffolding (matches Mathematica `debugprint`).               */
/*                                                                    */
/* When `Integrate``$Verbose` evaluates to True the trace helper       */
/* prints an `IN ...` line on entry / `OUT ...` line on exit of every  */
/* public helper.  Cheap when off (single OwnValue lookup).            */
/* ------------------------------------------------------------------ */

static bool intrat_trace_enabled(void) {
    SymbolDef* def = symtab_lookup("Integrate`$Verbose");
    if (!def || !def->own_values) return false;
    /* The OwnValue is `Integrate`$Verbose -> True/False`.  Walk it. */
    Rule* r = def->own_values;
    while (r) {
        if (r->replacement && r->replacement->type == EXPR_SYMBOL
            && r->replacement->data.symbol.name == SYM_True) return true;
        r = r->next;
    }
    return false;
}

static void intrat_trace(const char* fn, const char* dir, Expr* e) {
    if (!intrat_trace_enabled()) return;
    char* s = expr_to_string(e);
    fprintf(stderr, "[Integrate`] %s %s %s\n", fn, dir, s ? s : "(null)");
    if (s) free(s);
}

/* ------------------------------------------------------------------ */
/* Small wrappers around Mathilda internal_* / evaluate idioms.        */
/*                                                                    */
/* All of these hand the result back as a fresh Expr* (the caller     */
/* owns it).  None of them free their input arguments — every helper  */
/* makes its own copies.                                              */
/* ------------------------------------------------------------------ */

static Expr* call_eval(const char* head, Expr** args, size_t argc) {
    Expr* call = expr_new_function(expr_new_symbol(head), args, argc);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* Expand then free the input.  expr_expand borrows its argument and
 * returns a fresh tree, so any caller that owns a freshly-built temp
 * (e.g. internal_times(...) on the spot) must free it after the
 * expand or the temp leaks.  This wraps the common idiom. */
static inline Expr* expand_and_free(Expr* e) {
    if (!e) return NULL;
    Expr* r = expr_expand(e);
    expr_free(e);
    return r;
}

/* canonic[expr] = Cancel[Together[expr]].  Phase 1 keeps this
 * stripped down — the full Mathematica `canonic` adds a RootReduce /
 * ToRadicals tail that we'll bring online in later phases.
 *
 * The full Mathematica canonic also switches to
 * `Extension -> Automatic` only when the input contains a radical
 * (intsimp_has_radical); the Mathilda port stays plain-rational for
 * Phase 1 since our Cancel/Together don't yet honour Extension. */
Expr* intrat_canonic(Expr* e) {
    Expr* tg = internal_together((Expr*[]){expr_copy(e)}, 1);
    Expr* result = internal_cancel((Expr*[]){tg}, 1);
    return result;
}

static Expr* intrat_d(Expr* e, Expr* x) {
    Expr* args[2] = {expr_copy(e), expr_copy(x)};
    return call_eval("D", args, 2);
}

static Expr* intrat_numerator(Expr* e) {
    return internal_numerator((Expr*[]){expr_copy(e)}, 1);
}

Expr* intrat_denominator(Expr* e) {
    return internal_denominator((Expr*[]){expr_copy(e)}, 1);
}

static Expr* intrat_polyq(Expr* a, Expr* b, Expr* x) {
    return internal_polynomialquotient(
        (Expr*[]){expr_copy(a), expr_copy(b), expr_copy(x)}, 3);
}

/* PolynomialRemainder[a, b, var] over the named variable.  Used by
 * the Phase 2 LRT pipeline to reduce mod Q_i(t). */
static Expr* intrat_polyrem_t(Expr* a, Expr* b, Expr* var) {
    return internal_polynomialremainder(
        (Expr*[]){expr_copy(a), expr_copy(b), expr_copy(var)}, 3);
}

/* PolynomialQuotientRemainder, returning the {q, r} components separately. */
static bool intrat_polyqr(Expr* a, Expr* b, Expr* x, Expr** out_q, Expr** out_r) {
    Expr* args[3] = {expr_copy(a), expr_copy(b), expr_copy(x)};
    Expr* call = expr_new_function(
        expr_new_symbol(SYM_PolynomialQuotientRemainder), args, 3);
    Expr* result = evaluate(call);
    expr_free(call);
    if (!result || result->type != EXPR_FUNCTION
        || result->data.function.head->type != EXPR_SYMBOL
        || result->data.function.head->data.symbol.name != SYM_List
        || result->data.function.arg_count != 2) {
        if (result) expr_free(result);
        *out_q = NULL; *out_r = NULL;
        return false;
    }
    *out_q = expr_copy(result->data.function.args[0]);
    *out_r = expr_copy(result->data.function.args[1]);
    expr_free(result);
    return true;
}

/* exquo[a, b, x] = quotient asserting zero remainder.  Bronstein's
 * "exact quotient" — IntegrateRational.m:1238-1248.  On a non-zero
 * remainder we emit a soft warning (matching the reference) and
 * return the quotient anyway; downstream consumers will detect the
 * symbolic remainder during their own canonicalisation. */
static Expr* intrat_exquo(Expr* a, Expr* b, Expr* x) {
    Expr *q = NULL, *r = NULL;
    if (!intrat_polyqr(a, b, x, &q, &r)) return NULL;
    if (r && !is_zero_poly(r)) {
        fprintf(stderr,
            "Integrate`exquo::remainder failed to produce a zero remainder.\n");
    }
    if (r) expr_free(r);
    return q;
}

static bool intrat_polyq_test(Expr* expr, Expr* var) {
    Expr* args[2] = {expr_copy(expr), expr_copy(var)};
    Expr* call = internal_polynomialq(args, 2);
    Expr* val  = evaluate(call);
    expr_free(call);
    bool ok = (val && val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_True);
    if (val) expr_free(val);
    return ok;
}

bool intrat_freeq_test(Expr* expr, Expr* var) {
    Expr* args[2] = {expr_copy(expr), expr_copy(var)};
    Expr* call = internal_freeq(args, 2);
    Expr* val  = evaluate(call);
    expr_free(call);
    bool ok = (val && val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_True);
    if (val) expr_free(val);
    return ok;
}

/* Monic GCD in K(coeffs)[x].  Mathematica's HermiteReduce baseline
 * uses First[PolynomialExtendedGCD[a, b, x]] which returns the GCD
 * normalised so leading_coefficient(g, x) == 1.  Using the content-
 * normalised PolynomialGCD instead breaks the algorithm: when D has
 * non-trivial integer content (e.g. d = 4 x^2 (5-4x)^2) the dbar
 * term silently scales by content(d), the dstar / dbarhat / dtil
 * triplet downstream gets out of sync, and the f = D[g, x] + h
 * identity fails by a constant factor.
 *
 * Implementation: compute PolynomialGCD, then divide by the leading
 * coefficient in x.  PolynomialExtendedGCD itself is correct here but
 * orders of magnitude slower on hard inputs because it also tracks
 * Bezout coefficients we never use.
 *
 * Result is expanded — every downstream consumer feeds the gcd into
 * get_degree_poly, which sees through `Power[var, n]` only, not
 * `Power[base, n]` for multi-term `base`. */
static Expr* intrat_polygcd_monic(Expr* a, Expr* b, Expr* x) {
    Expr* g = internal_polynomialgcd(
        (Expr*[]){expr_copy(a), expr_copy(b)}, 2);
    Expr* g_e = expr_expand(g);
    expr_free(g);
    if (!g_e || is_zero_poly(g_e)) return g_e;
    int deg = get_degree_poly(g_e, x);
    if (deg <= 0) {
        /* Degree-0 gcd: normalise to 1 to keep the algorithm scale-
         * invariant (and so PolynomialQuotient[d, gcd] returns d/lc
         * faithfully). */
        expr_free(g_e);
        return expr_new_integer(1);
    }
    Expr* lc = get_coeff(g_e, x, deg);
    if (!lc) return g_e;
    /* If lc is already 1, nothing to do. */
    if (lc->type == EXPR_INTEGER && lc->data.integer == 1) {
        expr_free(lc);
        return g_e;
    }
    Expr* div = internal_divide((Expr*[]){g_e, lc}, 2);
    Expr* expanded = expr_expand(div);
    expr_free(div);
    return expanded;
}

/* lc[p, x] = leading coefficient of p in x.  When deg=0 the polynomial
 * is the constant itself. */
static Expr* intrat_lc(Expr* p, Expr* x) {
    int deg = get_degree_poly(p, x);
    if (deg <= 0) return expr_copy(p);
    return get_coeff(p, x, deg);
}

/* content[p, x] = PolynomialGCD @@ CoefficientList[p, x]. */
static Expr* intrat_content(Expr* p, Expr* x) {
    Expr* clist = internal_coefficientlist(
        (Expr*[]){expr_copy(p), expr_copy(x)}, 2);
    Expr* cl = evaluate(clist);
    expr_free(clist);
    if (!cl || cl->type != EXPR_FUNCTION
        || cl->data.function.head->type != EXPR_SYMBOL
        || cl->data.function.head->data.symbol.name != SYM_List) {
        if (cl) expr_free(cl);
        return expr_new_integer(1);
    }
    size_t n = cl->data.function.arg_count;
    if (n == 0) { expr_free(cl); return expr_new_integer(0); }
    /* Fold via PolynomialGCD.  The single-variable form
     * PolynomialGCD[a, b, c, ...] handles the multi-arg form. */
    Expr** args = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) args[i] = expr_copy(cl->data.function.args[i]);
    expr_free(cl);
    Expr* g = internal_polynomialgcd(args, n);
    /* internal_polynomialgcd (via internal_call_impl/expr_new_function)
     * memcpy's the slot pointers into its own backing store; the caller-
     * owned `args` buffer is no longer needed and must be released. */
    free(args);
    /* When coefficients carry parametric fractions (e.g. 1/b^3 factors
     * after Apart on a parametric integrand), our PolynomialGCD bails
     * and returns the call held.  Treat that as content == 1 — still
     * a correct (non-strict) primitive, and crucially keeps held GCD
     * calls from polluting downstream S polynomials. */
    if (head_is(g, SYM_PolynomialGCD)) {
        expr_free(g);
        g = expr_new_integer(1);
    }
    return g;
}

/* primitive[p, x] = p / content[p, x], expanded.  When p is free of x
 * the content equals p, so the result is 1; we shortcut to `p` to
 * match Mathematica's primitive[p, x] /; FreeQ[p, x] := p case. */
static Expr* intrat_primitive(Expr* p, Expr* x) {
    if (intrat_freeq_test(p, x)) return expr_copy(p);
    Expr* c = intrat_content(p, x);
    if (!c || (c->type == EXPR_INTEGER && c->data.integer == 1)) {
        if (c) expr_free(c);
        return expr_copy(p);
    }
    Expr* div = internal_divide(
        (Expr*[]){expr_copy(p), c}, 2);
    Expr* expanded = expr_expand(div);
    expr_free(div);
    return expanded;
}

/* monic[p, x] = p / lc[p, x], expanded.  When p is free of x we
 * follow Mathematica's convention and return 1. */
static Expr* intrat_monic_internal(Expr* p, Expr* x) {
    if (intrat_freeq_test(p, x)) return expr_new_integer(1);
    Expr* lc = intrat_lc(p, x);
    if (!lc) return expr_copy(p);
    if (lc->type == EXPR_INTEGER && lc->data.integer == 1) {
        expr_free(lc);
        return expr_expand(p);
    }
    Expr* div = internal_divide((Expr*[]){expr_copy(p), lc}, 2);
    Expr* result = expr_expand(div);
    expr_free(div);
    /* Cancel any rational coefficients that appear after the divide. */
    Expr* cancelled = internal_cancel((Expr*[]){result}, 1);
    return cancelled;
}

/* ------------------------------------------------------------------ */
/* ExtendedEuclidean (Diophantine).                                    */
/*                                                                    */
/* Given a, b, c with c in the ideal (a, b), return r, s' such that    */
/*    r*a + s'*b = c   and   deg(r) < deg(b).                          */
/*                                                                    */
/* Direct port of IntegrateRational.m:1203-1210.                       */
/* ------------------------------------------------------------------ */

static bool intrat_extended_euclidean(Expr* a, Expr* b, Expr* c, Expr* x,
                                      Expr** out_r, Expr** out_sp) {
    *out_r = *out_sp = NULL;

    Expr* eg_args[3] = { expr_copy(a), expr_copy(b), expr_copy(x) };
    Expr* eg = internal_polynomialextendedgcd(eg_args, 3);
    Expr* eg_eval = evaluate(eg);
    expr_free(eg);
    if (!eg_eval || eg_eval->type != EXPR_FUNCTION
        || eg_eval->data.function.arg_count != 2) {
        if (eg_eval) expr_free(eg_eval);
        return false;
    }
    Expr* g = eg_eval->data.function.args[0];
    Expr* coefs = eg_eval->data.function.args[1];
    if (coefs->type != EXPR_FUNCTION
        || coefs->data.function.head->type != EXPR_SYMBOL
        || coefs->data.function.head->data.symbol.name != SYM_List
        || coefs->data.function.arg_count != 2) {
        expr_free(eg_eval);
        return false;
    }

    Expr* g_copy = expr_copy(g);
    Expr* t_coef = expr_copy(coefs->data.function.args[0]);
    Expr* s_coef = expr_copy(coefs->data.function.args[1]);
    expr_free(eg_eval);

    /* q = c / g via PolynomialQuotient (g divides c). */
    Expr* q = intrat_polyq(c, g_copy, x);
    expr_free(g_copy);
    if (!q) {
        expr_free(t_coef); expr_free(s_coef);
        return false;
    }

    /* t = q*t_coef, s = q*s_coef, both expanded. */
    Expr* qt_raw = internal_times((Expr*[]){expr_copy(q), t_coef}, 2);
    Expr* t_new  = expr_expand(qt_raw); expr_free(qt_raw);
    Expr* qs_raw = internal_times((Expr*[]){expr_copy(q), s_coef}, 2);
    Expr* s_new  = expr_expand(qs_raw); expr_free(qs_raw);
    expr_free(q);

    /* {q2, r} = PolynomialQuotientRemainder[t_new, b, x] */
    Expr *q2 = NULL, *r = NULL;
    if (!intrat_polyqr(t_new, b, x, &q2, &r)) {
        expr_free(t_new); expr_free(s_new);
        return false;
    }
    expr_free(t_new);

    /* sp = Together[s_new + q2*a]. */
    Expr* q2a = internal_times((Expr*[]){q2, expr_copy(a)}, 2);
    Expr* sum = internal_plus((Expr*[]){s_new, q2a}, 2);
    Expr* sp = internal_together((Expr*[]){sum}, 1);

    *out_r = r;
    *out_sp = sp;
    return true;
}

/* AXIOM-style primitivePart(p, q) generalised for the LRT pipeline.
 * Mirrors SUBRESP.primitivePart in intrf.spad.pamphlet:48-52:
 *   r := extEuc(lc(p, var_p), q, 1, var_q).coef1     // cofactor mod q
 *   if r carries unsimplified radicals, Simplify it (otherwise the
 *     downstream PolynomialRemainder step on r*p over an algebraic
 *     extension fails to terminate -- the radical-coefficient
 *     coefficients don't reduce in the underlying algebraic field
 *     and the pseudo-division loop diverges)
 *   rem := PolynomialRemainder(r * p, q, var_q)
 *   return PrimitivePart(rem, var_p)
 *
 * Returns owned Expr* on success, NULL on extEuc / division failure.
 * Borrows p, q, var_p, var_q (caller retains ownership of all four).
 *
 * This consolidates two parallel patches in intrat_int_rational_log_part
 * (the i = deg_d case and the general i case).  Applying the radical-
 * Simplify guard uniformly is defensive -- costs nothing on rational-
 * coefficient inputs because intsimp_has_radical short-circuits. */
static Expr* intrat_primitive_part_mod(Expr* p, Expr* q,
                                       Expr* var_p, Expr* var_q) {
    Expr* lc_p = intrat_lc(p, var_p);
    Expr* one = expr_new_integer(1);
    Expr *r = NULL, *_unused = NULL;
    bool ok = intrat_extended_euclidean(lc_p, q, one, var_q, &r, &_unused);
    expr_free(lc_p);
    expr_free(one);
    if (!ok) return NULL;
    if (_unused) expr_free(_unused);

    r = intsimp_simplify_if_radical(r);

    Expr* rp_raw = internal_times(
        (Expr*[]){r, expr_copy(p)}, 2);
    Expr* rp = expr_expand(rp_raw); expr_free(rp_raw);
    Expr* rem = intrat_polyrem_t(rp, q, var_q);
    expr_free(rp);
    if (!rem) return NULL;
    /* S is the argument of a Log, so it is only defined up to a constant
     * multiple.  We take the integer-content-primitive part (rather than the
     * var_p-monic form): clearing rational denominators in the K(var_q)
     * coefficient field keeps the downstream real-log expansion (LogToReal /
     * LogToAtan over the residue-root algebraic field) on integer-coefficient
     * polynomials, which is markedly faster and avoids the fractional-coefficient
     * blow-up that stalls high-degree cases such as 1/(x^5 + 1). */
    Expr* out = intrat_primitive(rem, var_p);
    expr_free(rem);
    return out;
}

/* ------------------------------------------------------------------ */
/* HermiteReduce (Mack's linear version).                              */
/*                                                                    */
/* Direct port of IntegrateRational.m:1303-1323.  Returns {g, h} with  */
/* f == D[g, x] + h and Denominator[h] squarefree.                     */
/* ------------------------------------------------------------------ */

static Expr* intrat_hermite_reduce(Expr* f, Expr* x) {
    intrat_trace("HermiteReduce", "IN", f);

    /* Numerator / Denominator may return unexpanded results — every
     * downstream `get_degree_poly` call needs the polynomial in
     * expanded form to see the leading term, so we expand once
     * upfront. */
    Expr* a_raw = intrat_numerator(f);
    Expr* a = expr_expand(a_raw); expr_free(a_raw);
    Expr* d_raw = intrat_denominator(f);
    Expr* d = expr_expand(d_raw); expr_free(d_raw);
    Expr* g = expr_new_integer(0);

    /* dbar = gcd(d, D[d, x]) — monic gcd via PolynomialExtendedGCD. */
    Expr* dprime = intrat_d(d, x);
    Expr* dprime_e = expr_expand(dprime); expr_free(dprime);
    Expr* dbar = intrat_polygcd_monic(d, dprime_e, x);
    expr_free(dprime_e);
    if (!dbar) { expr_free(a); expr_free(d); expr_free(g); return NULL; }
    /* dstar = d / dbar */
    Expr* dstar = intrat_exquo(d, dbar, x);
    expr_free(d);

    while (true) {
        int deg_dbar = get_degree_poly(dbar, x);
        if (deg_dbar <= 0) break;

        /* dbartwo = monic-gcd(dbar, D[dbar, x]) */
        Expr* dbar_prime = intrat_d(dbar, x);
        Expr* dbar_prime_e = expr_expand(dbar_prime); expr_free(dbar_prime);
        Expr* dbartwo = intrat_polygcd_monic(dbar, dbar_prime_e, x);
        expr_free(dbar_prime_e);
        if (!dbartwo) {
            expr_free(a); expr_free(dbar); expr_free(dstar); expr_free(g);
            return NULL;
        }

        /* dbarhat = dbar / dbartwo */
        Expr* dbarhat = intrat_exquo(dbar, dbartwo, x);

        /* dtil = -exquo(dstar * D[dbar, x], dbar, x) */
        Expr* dbar_prime2 = intrat_d(dbar, x);
        Expr* dt_raw = internal_times((Expr*[]){expr_copy(dstar), dbar_prime2}, 2);
        Expr* numerator_dtil = expr_expand(dt_raw); expr_free(dt_raw);
        Expr* exq = intrat_exquo(numerator_dtil, dbar, x);
        expr_free(numerator_dtil);
        Expr* dtil_raw = internal_times((Expr*[]){expr_new_integer(-1), exq}, 2);
        Expr* dtil = expr_expand(dtil_raw); expr_free(dtil_raw);

        /* {b_coef, c_coef} = ExtendedEuclidean(dtil, dbarhat, a, x) */
        Expr *b_coef = NULL, *c_coef = NULL;
        if (!intrat_extended_euclidean(dtil, dbarhat, a, x, &b_coef, &c_coef)) {
            expr_free(dtil); expr_free(dbarhat); expr_free(dbartwo);
            expr_free(a); expr_free(dbar); expr_free(dstar); expr_free(g);
            return NULL;
        }

        /* a := c - exquo(D[b, x] * dstar, dbarhat, x) */
        Expr* b_prime = intrat_d(b_coef, x);
        Expr* bp_raw = internal_times((Expr*[]){b_prime, expr_copy(dstar)}, 2);
        Expr* bp_dstar = expr_expand(bp_raw); expr_free(bp_raw);
        Expr* corr = intrat_exquo(bp_dstar, dbarhat, x);
        expr_free(bp_dstar);
        Expr* neg_corr = internal_times(
            (Expr*[]){expr_new_integer(-1), corr}, 2);
        Expr* new_a = internal_plus((Expr*[]){c_coef, neg_corr}, 2);
        new_a = eval_and_free(new_a);
        expr_free(a);
        a = new_a;

        /* g := g + canonic(b/dbar) */
        Expr* bd = internal_divide(
            (Expr*[]){b_coef, expr_copy(dbar)}, 2);
        Expr* bd_can = intrat_canonic(bd);
        expr_free(bd);
        Expr* new_g = internal_plus((Expr*[]){g, bd_can}, 2);
        new_g = eval_and_free(new_g);
        g = new_g;

        /* dbar := dbartwo */
        expr_free(dbar); expr_free(dbarhat); expr_free(dtil);
        dbar = dbartwo;
    }
    expr_free(dbar);

    /* h = canonic(a/dstar) */
    Expr* h_raw = internal_divide(
        (Expr*[]){a, dstar}, 2);
    Expr* h = intrat_canonic(h_raw);
    expr_free(h_raw);

    Expr* result = expr_new_function(expr_new_symbol(SYM_List),
        (Expr*[]){g, h}, 2);
    intrat_trace("HermiteReduce", "OUT", result);
    return result;
}

/* ------------------------------------------------------------------ */
/* IntegratePolynomial — IntegrateRational.m:1153-1162.                */
/* ------------------------------------------------------------------ */

static Expr* intrat_integrate_polynomial(Expr* poly, Expr* x) {
    intrat_trace("IntegratePolynomial", "IN", poly);

    /* CoefficientList[Expand[poly], x] gives {a_0, a_1, ..., a_n}. */
    Expr* expanded = expr_expand(poly);
    Expr* cl_call = internal_coefficientlist(
        (Expr*[]){expanded, expr_copy(x)}, 2);
    Expr* cl = evaluate(cl_call);
    expr_free(cl_call);
    if (!cl || cl->type != EXPR_FUNCTION
        || cl->data.function.head->type != EXPR_SYMBOL
        || cl->data.function.head->data.symbol.name != SYM_List) {
        if (cl) expr_free(cl);
        /* Fall back to the trivial constant case: a polynomial that is
         * a pure constant has CoefficientList -> {const}. */
        return NULL;
    }
    size_t n = cl->data.function.arg_count;
    if (n == 0) { expr_free(cl); return expr_new_integer(0); }

    /* Build Σ a_k x^(k+1) / (k+1). */
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t k = 0; k < n; k++) {
        Expr* a = cl->data.function.args[k];
        Expr* x_pow = internal_power(
            (Expr*[]){expr_copy(x), expr_new_integer((int64_t)(k + 1))}, 2);
        Expr* div = internal_divide(
            (Expr*[]){x_pow, expr_new_integer((int64_t)(k + 1))}, 2);
        terms[k] = internal_times(
            (Expr*[]){expr_copy(a), div}, 2);
    }
    expr_free(cl);
    Expr* sum = internal_plus(terms, n);
    free(terms);
    Expr* sum_eval = eval_and_free(sum);

    intrat_trace("IntegratePolynomial", "OUT", sum_eval);
    return sum_eval;
}

/* ------------------------------------------------------------------ */
/* Pre-Hermite derivative-recognition fast path.                        */
/*                                                                    */
/* If f = num/den with FactorSquareFree[den] = c0 * pol^k for a single */
/* squarefree pol, and num == c * D[pol, x] for some constant c, then  */
/*   k == 1: ∫ f dx = (c/c0) Log[pol]                                  */
/*   k >= 2: ∫ f dx = -(c/c0) / ((k-1) * pol^(k-1))                    */
/*                                                                    */
/* Returns NULL if the structure does not match — the caller falls     */
/* back to the standard Hermite + LRT pipeline.                        */
/* ------------------------------------------------------------------ */

/* Walk a FactorSquareFree[Times[..., pol^k, ...]] result and decide
 * whether it has the shape `c0 * pol^k` for a single non-constant
 * squarefree polynomial `pol` in `x`.  Returns true on match and
 * fills `*c0_out`, `*pol_out`, `*k_out`. */
static bool extract_factored_denominator(Expr* fs, Expr* x,
                                         Expr** c0_out, Expr** pol_out,
                                         int* k_out) {
    Expr* fs_expanded = fs;
    /* Times-headed: fs = c0 * factor1 * factor2 * ... */
    if (fs_expanded->type == EXPR_FUNCTION
        && fs_expanded->data.function.head->type == EXPR_SYMBOL
        && fs_expanded->data.function.head->data.symbol.name == SYM_Times) {
        Expr* c0 = expr_new_integer(1);
        Expr* pol = NULL;
        int k = 0;
        for (size_t i = 0; i < fs_expanded->data.function.arg_count; i++) {
            Expr* arg = fs_expanded->data.function.args[i];
            if (intrat_freeq_test(arg, x)) {
                Expr* prod = internal_times(
                    (Expr*[]){c0, expr_copy(arg)}, 2);
                c0 = eval_and_free(prod);
                continue;
            }
            /* Non-constant factor: must be unique and Power[pol, k] form. */
            if (pol) { expr_free(c0); return false; }
            if (arg->type == EXPR_FUNCTION
                && arg->data.function.head->type == EXPR_SYMBOL
                && arg->data.function.head->data.symbol.name == SYM_Power
                && arg->data.function.arg_count == 2
                && arg->data.function.args[1]->type == EXPR_INTEGER
                && arg->data.function.args[1]->data.integer >= 1) {
                pol = expr_copy(arg->data.function.args[0]);
                k = (int)arg->data.function.args[1]->data.integer;
            } else {
                pol = expr_copy(arg);
                k = 1;
            }
        }
        if (!pol) { expr_free(c0); return false; }
        *c0_out = c0; *pol_out = pol; *k_out = k;
        return true;
    }
    /* Single Power[pol, k]. */
    if (fs_expanded->type == EXPR_FUNCTION
        && fs_expanded->data.function.head->type == EXPR_SYMBOL
        && fs_expanded->data.function.head->data.symbol.name == SYM_Power
        && fs_expanded->data.function.arg_count == 2
        && fs_expanded->data.function.args[1]->type == EXPR_INTEGER
        && fs_expanded->data.function.args[1]->data.integer >= 1
        && !intrat_freeq_test(fs_expanded->data.function.args[0], x)) {
        *c0_out = expr_new_integer(1);
        *pol_out = expr_copy(fs_expanded->data.function.args[0]);
        *k_out = (int)fs_expanded->data.function.args[1]->data.integer;
        return true;
    }
    /* Bare polynomial. */
    if (!intrat_freeq_test(fs_expanded, x)) {
        *c0_out = expr_new_integer(1);
        *pol_out = expr_copy(fs_expanded);
        *k_out = 1;
        return true;
    }
    return false;
}

/* Returns the closed-form integral if the derivative-recognition
 * pattern matches; NULL otherwise. */
static Expr* intrat_derivative_recognition(Expr* num, Expr* den, Expr* x) {
    /* Factor the denominator into (constant) * (single squarefree pol)^k. */
    Expr* fs = internal_factorsquarefree(
        (Expr*[]){expr_copy(den)}, 1);
    Expr* fs_eval = evaluate(fs);
    expr_free(fs);
    if (!fs_eval) return NULL;

    Expr *c0 = NULL, *pol = NULL;
    int k = 0;
    bool ok = extract_factored_denominator(fs_eval, x, &c0, &pol, &k);
    expr_free(fs_eval);
    if (!ok) return NULL;

    /* Compute c = num / (c0 * D[pol, x]).  If FreeQ[c, x] we have the
     * derivative-recognition form. */
    Expr* dpol = intrat_d(pol, x);
    Expr* c0_dpol_raw = internal_times((Expr*[]){expr_copy(c0), dpol}, 2);
    Expr* c0_dpol = expr_expand(c0_dpol_raw); expr_free(c0_dpol_raw);
    Expr* ratio = internal_divide(
        (Expr*[]){expr_copy(num), c0_dpol}, 2);
    Expr* c_const = intrat_canonic(ratio);
    expr_free(ratio);

    if (!intrat_freeq_test(c_const, x)) {
        expr_free(c0); expr_free(pol); expr_free(c_const);
        return NULL;
    }

    Expr* result;
    if (k == 1) {
        /* ∫ c * D[pol]/pol dx = c * Log[pol] */
        Expr* log_pol = expr_new_function(expr_new_symbol(SYM_Log),
            (Expr*[]){pol}, 1);
        result = internal_times(
            (Expr*[]){c_const, log_pol}, 2);
    } else {
        /* ∫ c * D[pol]/pol^k dx = -c / ((k-1) * pol^(k-1)) */
        Expr* pol_pow = internal_power(
            (Expr*[]){pol, expr_new_integer((int64_t)(k - 1))}, 2);
        Expr* denom = internal_times(
            (Expr*[]){expr_new_integer((int64_t)(k - 1)), pol_pow}, 2);
        Expr* div = internal_divide(
            (Expr*[]){c_const, denom}, 2);
        result = internal_times(
            (Expr*[]){expr_new_integer(-1), div}, 2);
    }
    expr_free(c0);
    return eval_and_free(result);
}

/* ====================================================================
 * Phase 2 — Lazard-Rioboo-Trager log part
 * ====================================================================
 *
 * Implements `IntRationalLogPart[A/D, x, t]` (Bronstein, Symbolic
 * Integration I, p. 51) for square-free denominators D in K[x],
 * along with the supporting helpers SquareFree, ExtractConstants,
 * ApartList, and a structural pre-detector for the cyclotomic /
 * n-th-root denominator family that appears repeatedly in Bronstein
 * and the BronsteinRational test suite.
 *
 * Output form (RootSum -> False, the default for Phase 2's structural
 * tests):
 *   list of {Q_i(t), S_i(t, x)} pairs encoding
 *      ∫ A/D dx == Σ_i RootSum[Q_i(t)==0, t Log[S_i(t, x)]].
 *
 * Phase 4's LogToReal will consume the same (Q_i, S_i) pairs and
 * convert each one into real Log + ArcTan form when Q_i factors
 * tractably.  Until then, the Phase 2 builtin returns the symbolic
 * RootSum sum (RootSum -> True) or the raw pair list (RootSum ->
 * False).
 */

/* ----- SquareFree: {poly, multiplicity} list ----- */

/* Bucket FactorSquareFree[p]'s output into a list of (poly,
 * multiplicity) pairs, indexed densely by multiplicity from 1..max.
 * If the input has no factor of multiplicity i, the i-th entry is
 * {1, i}.  Mirrors the BronsteinRational baseline's `SquareFree[]`
 * helper at IntegrateRational.m:1474-1487.
 *
 * Implementation note: Mathilda's FactorSquareFree returns the
 * factorisation as `Times[c, Power[p1, k1], Power[p2, k2], ...]`
 * (or a single Power[]/atomic factor).  We walk that structure,
 * grouping factors by multiplicity. */
static Expr* intrat_squarefree_list(Expr* p) {
    /* Run FactorSquareFree first. */
    Expr* fs = internal_factorsquarefree((Expr*[]){expr_copy(p)}, 1);
    Expr* fs_eval = evaluate(fs);
    expr_free(fs);
    if (!fs_eval) return NULL;

    /* Walk the factor structure, accumulating (factor, mult) pairs. */
    typedef struct { Expr* poly; int mult; } Pair;
    Pair* pairs = (Pair*)malloc(sizeof(Pair) * 8);
    size_t pcap = 8, pcount = 0;

#define PUSH_PAIR(POLY, MULT) do { \
        if (pcount >= pcap) { pcap *= 2; \
            pairs = (Pair*)realloc(pairs, sizeof(Pair) * pcap); } \
        pairs[pcount].poly = (POLY); \
        pairs[pcount].mult = (MULT); \
        pcount++; \
    } while (0)

    if (fs_eval->type == EXPR_FUNCTION
        && fs_eval->data.function.head->type == EXPR_SYMBOL
        && fs_eval->data.function.head->data.symbol.name == SYM_Times) {
        for (size_t i = 0; i < fs_eval->data.function.arg_count; i++) {
            Expr* arg = fs_eval->data.function.args[i];
            if (arg->type == EXPR_FUNCTION
                && arg->data.function.head->type == EXPR_SYMBOL
                && arg->data.function.head->data.symbol.name == SYM_Power
                && arg->data.function.arg_count == 2
                && arg->data.function.args[1]->type == EXPR_INTEGER) {
                int k = (int)arg->data.function.args[1]->data.integer;
                if (k >= 1) {
                    PUSH_PAIR(expr_copy(arg->data.function.args[0]), k);
                    continue;
                }
            }
            PUSH_PAIR(expr_copy(arg), 1);
        }
    } else if (fs_eval->type == EXPR_FUNCTION
        && fs_eval->data.function.head->type == EXPR_SYMBOL
        && fs_eval->data.function.head->data.symbol.name == SYM_Power
        && fs_eval->data.function.arg_count == 2
        && fs_eval->data.function.args[1]->type == EXPR_INTEGER) {
        int k = (int)fs_eval->data.function.args[1]->data.integer;
        if (k >= 1) PUSH_PAIR(expr_copy(fs_eval->data.function.args[0]), k);
        else PUSH_PAIR(expr_copy(fs_eval), 1);
    } else {
        /* Single bare factor — multiplicity 1. */
        PUSH_PAIR(expr_copy(fs_eval), 1);
    }
    expr_free(fs_eval);

#undef PUSH_PAIR

    /* Determine max multiplicity. */
    int max_mult = 0;
    for (size_t i = 0; i < pcount; i++) {
        if (pairs[i].mult > max_mult) max_mult = pairs[i].mult;
    }
    if (max_mult == 0) {
        /* Constant input — return {{c, 1}}. */
        Expr* pair = expr_new_function(expr_new_symbol(SYM_List),
            (Expr*[]){expr_copy(p), expr_new_integer(1)}, 2);
        Expr* list = expr_new_function(expr_new_symbol(SYM_List),
            (Expr*[]){pair}, 1);
        for (size_t i = 0; i < pcount; i++) expr_free(pairs[i].poly);
        free(pairs);
        return list;
    }

    /* Build output list of length max_mult.  Index i (1-based) has
     * the product of all factors of multiplicity i, or 1 if none. */
    Expr** out_args = (Expr**)malloc(sizeof(Expr*) * max_mult);
    for (int m = 1; m <= max_mult; m++) {
        /* Collect all pairs with this multiplicity. */
        Expr** factors = NULL;
        size_t fcount = 0;
        for (size_t i = 0; i < pcount; i++) {
            if (pairs[i].mult == m) {
                factors = (Expr**)realloc(factors, sizeof(Expr*) * (fcount + 1));
                factors[fcount++] = expr_copy(pairs[i].poly);
            }
        }
        Expr* poly_for_m;
        if (fcount == 0) {
            poly_for_m = expr_new_integer(1);
        } else if (fcount == 1) {
            poly_for_m = expr_expand(factors[0]);
            expr_free(factors[0]);
        } else {
            Expr* prod = internal_times(factors, fcount);
            poly_for_m = expr_expand(prod);
            expr_free(prod);
        }
        if (factors) free(factors);
        out_args[m - 1] = expr_new_function(expr_new_symbol(SYM_List),
            (Expr*[]){poly_for_m, expr_new_integer(m)}, 2);
    }
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), out_args, max_mult);
    free(out_args);

    for (size_t i = 0; i < pcount; i++) expr_free(pairs[i].poly);
    free(pairs);
    return list;
}

/* Extract the i-th element of a Mathematica-style List expression.
 * Returns a fresh copy.  i is 1-based to match Mathematica.  Returns
 * NULL if `lst` is not a List or i is out of range. */
static Expr* list_get(Expr* lst, size_t i) {
    if (!lst || lst->type != EXPR_FUNCTION
        || lst->data.function.head->type != EXPR_SYMBOL
        || lst->data.function.head->data.symbol.name != SYM_List
        || i == 0 || i > lst->data.function.arg_count) return NULL;
    return expr_copy(lst->data.function.args[i - 1]);
}

/* ----- ExtractConstants: pull constant prefactors out of f ----- */

/* Split p into (numeric prefactor, x-bearing residue).  Walk a Times
 * head; every factor that is FreeQ[_, x] absorbs into the prefactor,
 * everything else stays in the residue.  When the input is a single
 * factor we either return (1, p) or (p, 1).  Mirrors the part of
 * IntegrateRational.m:1013-1023 that wraps FactorSquareFreeList +
 * Power@@first, but we work directly on Times args to avoid an extra
 * factorisation pass. */
static void intrat_split_constant(Expr* p, Expr* x,
                                  Expr** out_const, Expr** out_residue) {
    if (intrat_freeq_test(p, x)) {
        *out_const = expr_copy(p);
        *out_residue = expr_new_integer(1);
        return;
    }
    if (p->type != EXPR_FUNCTION
        || p->data.function.head->type != EXPR_SYMBOL
        || p->data.function.head->data.symbol.name != SYM_Times) {
        *out_const = expr_new_integer(1);
        *out_residue = expr_copy(p);
        return;
    }
    /* Times-headed: split args. */
    size_t n = p->data.function.arg_count;
    Expr** consts = (Expr**)malloc(sizeof(Expr*) * n);
    Expr** resid  = (Expr**)malloc(sizeof(Expr*) * n);
    size_t cn = 0, rn = 0;
    for (size_t i = 0; i < n; i++) {
        if (intrat_freeq_test(p->data.function.args[i], x)) {
            consts[cn++] = expr_copy(p->data.function.args[i]);
        } else {
            resid[rn++] = expr_copy(p->data.function.args[i]);
        }
    }
    *out_const  = (cn == 0) ? expr_new_integer(1)
                : (cn == 1) ? consts[0]
                : eval_and_free(internal_times(consts, cn));
    *out_residue = (rn == 0) ? expr_new_integer(1)
                : (rn == 1) ? resid[0]
                : eval_and_free(internal_times(resid, rn));
    free(consts); free(resid);
}

/* extractConstants[f, x] = {numconst/denconst, num_residue/den_residue}.
 * Returns a List[const, simplified] expression. */
static Expr* intrat_extract_constants(Expr* f, Expr* x) {
    Expr* num = intrat_numerator(f);
    Expr* den = intrat_denominator(f);
    Expr* num_e = eval_and_free(num);
    Expr* den_e = eval_and_free(den);
    /* Apply FactorSquareFree once to surface the constant prefactor. */
    Expr* num_fs = internal_factorsquarefree((Expr*[]){num_e}, 1);
    num_fs = eval_and_free(num_fs);
    Expr* den_fs = internal_factorsquarefree((Expr*[]){den_e}, 1);
    den_fs = eval_and_free(den_fs);

    Expr *nc = NULL, *nr = NULL, *dc = NULL, *dr = NULL;
    intrat_split_constant(num_fs, x, &nc, &nr);
    intrat_split_constant(den_fs, x, &dc, &dr);
    expr_free(num_fs); expr_free(den_fs);

    Expr* const_part = internal_divide(
        (Expr*[]){nc, dc}, 2);
    const_part = eval_and_free(const_part);
    Expr* residue = internal_divide(
        (Expr*[]){nr, dr}, 2);
    residue = eval_and_free(residue);

    return expr_new_function(expr_new_symbol(SYM_List),
        (Expr*[]){const_part, residue}, 2);
}

/* ----- ApartList: partial-fraction expansion as a list ----- */

/* apartList[f, x] = list of partial-fraction terms.  Phase 2 keeps
 * the lightweight implementation: call Apart[f, x], and split the
 * resulting Plus into a list of summands.  When Apart returns a
 * single non-Plus term we wrap it in a singleton list. */
static Expr* intrat_apart_list(Expr* f, Expr* x, const Expr* alpha) {
    /* Build Apart[f, x, Extension -> alpha?]. */
    size_t argc = (alpha ? 3 : 2);
    Expr** args = (Expr**)malloc(sizeof(Expr*) * argc);
    args[0] = expr_copy(f);
    args[1] = expr_copy(x);
    if (alpha) {
        args[2] = expr_new_function(expr_new_symbol(SYM_Rule),
            (Expr*[]){expr_new_symbol(SYM_Extension), expr_copy((Expr*)alpha)}, 2);
    }
    Expr* call = expr_new_function(expr_new_symbol(SYM_Apart), args, argc);
    Expr* ap = evaluate(call);
    expr_free(call);
    free(args);

    if (!ap) return NULL;

    if (ap->type == EXPR_FUNCTION
        && ap->data.function.head->type == EXPR_SYMBOL
        && ap->data.function.head->data.symbol.name == SYM_Plus) {
        size_t n = ap->data.function.arg_count;
        Expr** terms = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) terms[i] = expr_copy(ap->data.function.args[i]);
        expr_free(ap);
        /* expr_new_function memcpy's the args buffer into its own
         * backing store, so the caller's `terms` buffer must be freed. */
        Expr* list = expr_new_function(expr_new_symbol(SYM_List), terms, n);
        free(terms);
        return list;
    }
    return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ap}, 1);
}

/* ----- Cyclotomic / nth-root pre-detection ----- */

/* nthrootQ[p, x]: True iff p has the form `a x^n + b` with a, b free
 * of x.  Caught by Mathematica's MatchQ pattern at
 * IntegrateRational.m:1093.  The integer n is returned via *n_out.
 * Reserved for Phase 4's LogToReal closer (n-th-root-of-unity case);
 * unused in Phase 2 but kept here so the helper does not need to be
 * resurrected when Phase 4 lands. */
MATHILDA_MAYBE_UNUSED
static bool intrat_nthroot_q(Expr* p, Expr* x, int* n_out, Expr** a_out, Expr** b_out) {
    /* Run Collect[p, x] to gather the polynomial in canonical form. */
    Expr* coll = expr_new_function(expr_new_symbol(SYM_Collect),
        (Expr*[]){expr_copy(p), expr_copy(x)}, 2);
    Expr* c = evaluate(coll);
    expr_free(coll);
    if (!c) return false;

    /* Walk a Plus-headed result. */
    Expr* a = NULL;
    Expr* b = NULL;
    int n = 0;

    Expr** terms;
    size_t nterms;
    if (c->type == EXPR_FUNCTION
        && c->data.function.head->type == EXPR_SYMBOL
        && c->data.function.head->data.symbol.name == SYM_Plus) {
        terms = c->data.function.args;
        nterms = c->data.function.arg_count;
    } else {
        terms = &c;
        nterms = 1;
    }

    for (size_t i = 0; i < nterms; i++) {
        Expr* term = terms[i];
        if (intrat_freeq_test(term, x)) {
            if (b) { expr_free(b); expr_free(c); if (a) expr_free(a); return false; }
            b = expr_copy(term);
            continue;
        }
        /* term must be c * x^n with c free of x (or x^n alone, c=1). */
        Expr* coeff = NULL;
        int deg = 0;
        if (term->type == EXPR_FUNCTION
            && term->data.function.head->type == EXPR_SYMBOL
            && term->data.function.head->data.symbol.name == SYM_Power
            && term->data.function.arg_count == 2
            && expr_eq(term->data.function.args[0], x)
            && term->data.function.args[1]->type == EXPR_INTEGER) {
            deg = (int)term->data.function.args[1]->data.integer;
            coeff = expr_new_integer(1);
        } else if (expr_eq(term, x)) {
            deg = 1; coeff = expr_new_integer(1);
        } else if (term->type == EXPR_FUNCTION
            && term->data.function.head->type == EXPR_SYMBOL
            && term->data.function.head->data.symbol.name == SYM_Times) {
            /* coeff * x^n form. */
            Expr** factors = term->data.function.args;
            size_t nf = term->data.function.arg_count;
            Expr** ccoll = (Expr**)malloc(sizeof(Expr*) * nf);
            size_t ccn = 0;
            for (size_t k = 0; k < nf; k++) {
                Expr* f = factors[k];
                if (intrat_freeq_test(f, x)) {
                    ccoll[ccn++] = expr_copy(f);
                } else if (deg == 0 && f->type == EXPR_FUNCTION
                    && f->data.function.head->type == EXPR_SYMBOL
                    && f->data.function.head->data.symbol.name == SYM_Power
                    && f->data.function.arg_count == 2
                    && expr_eq(f->data.function.args[0], x)
                    && f->data.function.args[1]->type == EXPR_INTEGER) {
                    deg = (int)f->data.function.args[1]->data.integer;
                } else if (deg == 0 && expr_eq(f, x)) {
                    deg = 1;
                } else {
                    deg = -1; break;
                }
            }
            if (deg <= 0) {
                for (size_t kk = 0; kk < ccn; kk++) expr_free(ccoll[kk]);
                free(ccoll);
                if (a) expr_free(a);
                if (b) expr_free(b);
                expr_free(c);
                return false;
            }
            coeff = (ccn == 0) ? expr_new_integer(1)
                  : (ccn == 1) ? ccoll[0]
                  : eval_and_free(internal_times(ccoll, ccn));
            free(ccoll);
        } else {
            if (a) expr_free(a);
            if (b) expr_free(b);
            expr_free(c);
            return false;
        }
        if (a) { expr_free(a); expr_free(b ? b : NULL); expr_free(c);
                 expr_free(coeff); return false; }
        a = coeff;
        n = deg;
    }

    expr_free(c);
    if (!a || n <= 1) {
        if (a) expr_free(a);
        if (b) expr_free(b);
        return false;
    }
    *n_out = n;
    *a_out = a;
    *b_out = b ? b : expr_new_integer(0);
    return true;
}

/* ----- IntRationalLogPart: the LRT log-part computation ----- */

/* Find the unique element of the prs chain whose degree in x equals
 * `target_deg`.  Returns a fresh copy or NULL if none exists. */
static Expr* find_prs_at_degree(Expr* prs, Expr* x, int target_deg) {
    if (prs->type != EXPR_FUNCTION
        || prs->data.function.head->type != EXPR_SYMBOL
        || prs->data.function.head->data.symbol.name != SYM_List) return NULL;
    for (size_t i = 0; i < prs->data.function.arg_count; i++) {
        Expr* el = prs->data.function.args[i];
        Expr* el_e = expr_expand(el);
        int d = get_degree_poly(el_e, x);
        if (d == target_deg) {
            return el_e;
        }
        expr_free(el_e);
    }
    return NULL;
}

/* ------------------------------------------------------------------ */
/* Simple-root expansion (Phase 8d-bonus).                             */
/* ------------------------------------------------------------------ */

/* Tiny tree builders used only by the radical-form root extraction.
 * Each takes ownership of every Expr* arg and returns a fresh tree
 * the caller owns.  All paths route through eval_and_free at the
 * end so the outer Plus/Times can fold trivial identities. */
static Expr* nlp_neg(Expr* a) {
    return internal_times((Expr*[]){ expr_new_integer(-1), a }, 2);
}
static Expr* nlp_sqrt(Expr* a) {
    /* Power[a, 1/2] — Mathilda's Power evaluator turns this into the
     * canonical Sqrt[..] surface form during printing. */
    return expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ a,
                   expr_new_function(expr_new_symbol(SYM_Rational),
                       (Expr*[]){ expr_new_integer(1),
                                  expr_new_integer(2) }, 2) },
        2);
}
/* Forward declarations for helpers defined later in this file but
 * referenced by the palindromic-quartic expander immediately below. */
static Expr* subst_t(Expr* e, Expr* t, Expr* sub_expr);
static void split_re_im(Expr* p, Expr** re_out, Expr** im_out);

static Expr* nlp_sub_body(Expr* body, Expr* bvar, Expr* value) {
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){ expr_copy(bvar), expr_copy(value) }, 2);
    Expr* sub = internal_replace_all(
        (Expr*[]){ expr_copy(body), rule }, 2);
    return eval_and_free(sub);
}

/* True iff every variable in `poly` is structurally equal to `bvar`.
 * Phase 8d-bonus only fires on non-parametric polynomials: when poly
 * has parameters (e.g. b + a t^4), the radical-formula substitution
 * produces forms whose D[..., x] - integrand reduces to 0 only under
 * a stronger algebraic-extension simplifier than Mathilda currently
 * has.  Bailing out for parametric cases preserves the held RootSum
 * form, which the user can verify formally and which avoids
 * apparent-regression noise on the corpus. */
static bool poly_only_uses(Expr* poly, Expr* bvar) {
    Expr* vars_call = expr_new_function(expr_new_symbol(SYM_Variables),
        (Expr*[]){ expr_copy(poly) }, 1);
    Expr* vars = evaluate(vars_call);
    expr_free(vars_call);
    if (!vars || vars->type != EXPR_FUNCTION
        || vars->data.function.head->type != EXPR_SYMBOL
        || vars->data.function.head->data.symbol.name != SYM_List) {
        if (vars) expr_free(vars);
        return false;
    }
    bool ok = true;
    for (size_t i = 0; i < vars->data.function.arg_count; i++) {
        if (!expr_eq(vars->data.function.args[i], bvar)) { ok = false; break; }
    }
    expr_free(vars);
    return ok;
}

/* Build the real form of
 *   RootSum[ Function[t, d_t],
 *            Function[t, a_t · Log[x − t] / dd_t] ]
 * for a palindromic quartic d_t = c0 t^4 + c1 t^3 + c2 t^2 + c1 t + c0
 * (c0 == c4, c1 == c3).
 *
 * Strategy: factor d_t over R via the u = t + 1/t substitution
 *   c0 u^2 + c1 u + (c2 − 2 c0) = 0  ⇒  u = u_±,
 * giving two complex-conjugate root pairs of d_t:
 *   α_k = u_k/2 + i v_k,   v_k = Sqrt[(4 − u_k^2)/4].
 * Each conjugate pair contributes
 *   2 Re[ (a(α_k)/d'(α_k)) · Log[x − α_k] ]
 *     = C_re · Log[(x − u_k/2)^2 + v_k^2]
 *     + 2 C_im · ArcTan[v_k/(x − u_k/2)]
 * where C = a(α_k)/d'(α_k) = C_re + i C_im.
 *
 * Avoids LogToReal/LogToAtan entirely (those bog down on the nested-
 * radical S polynomials from the scaled-palindromic LRT path).
 * Returns NULL when d_t isn't palindromic, when u_± aren't real, or
 * when any sign-gate is undecidable. */
static Expr* expand_palindromic_quartic_real(
    Expr* a_t, Expr* d_t, Expr* dd_t, Expr* t, Expr* x)
{
    /* Detect palindromic quartic in t. */
    int deg = get_degree_poly(d_t, t);
    if (deg != 4) return NULL;
    Expr* c0 = get_coeff(d_t, t, 0);
    Expr* c1 = get_coeff(d_t, t, 1);
    Expr* c2 = get_coeff(d_t, t, 2);
    Expr* c3 = get_coeff(d_t, t, 3);
    Expr* c4 = get_coeff(d_t, t, 4);

    Expr* c0mc4 = internal_subtract(
        (Expr*[]){expr_copy(c0), expr_copy(c4)}, 2);
    Expr* c1mc3 = internal_subtract(
        (Expr*[]){expr_copy(c1), expr_copy(c3)}, 2);
    bool palindromic = intsimp_zero_q(c0mc4) && intsimp_zero_q(c1mc3);
    expr_free(c0mc4); expr_free(c1mc3);

    bool c0_nonzero = !(c0->type == EXPR_INTEGER && c0->data.integer == 0);
    if (!palindromic || !c0_nonzero) {
        expr_free(c0); expr_free(c1); expr_free(c2);
        expr_free(c3); expr_free(c4);
        return NULL;
    }

    /* Outer disc Du = c1^2 − 4 c0 (c2 − 2 c0).  Soundness gate:
     * Du >= 0 so u_± are real. */
    Expr* c1sq = expand_and_free(internal_power(
        (Expr*[]){expr_copy(c1), expr_new_integer(2)}, 2));
    Expr* two_c0 = expand_and_free(internal_times(
        (Expr*[]){expr_new_integer(2), expr_copy(c0)}, 2));
    Expr* c2m2c0 = eval_and_free(internal_subtract(
        (Expr*[]){expr_copy(c2), two_c0}, 2));
    Expr* four_term = expand_and_free(internal_times(
        (Expr*[]){expr_new_integer(4),
                  expr_copy(c0),
                  c2m2c0}, 3));
    Expr* neg_four_term = expand_and_free(internal_times(
        (Expr*[]){expr_new_integer(-1), four_term}, 2));
    Expr* Du = eval_and_free(internal_plus(
        (Expr*[]){c1sq, neg_four_term}, 2));
    Expr* Du_exp = expr_expand(Du); expr_free(Du);
    int Du_sign = intsimp_sign_pos_assumption(Du_exp);
    if (Du_sign == 0) Du_sign = intsimp_numeric_sign(Du_exp);
    if (Du_sign < 0) {
        expr_free(Du_exp);
        expr_free(c0); expr_free(c1); expr_free(c2);
        expr_free(c3); expr_free(c4);
        return NULL;
    }
    Expr* sqrt_Du = intsimp_pos_sqrt(Du_exp);
    expr_free(Du_exp);

    Expr* neg_c1 = expand_and_free(internal_times(
        (Expr*[]){expr_new_integer(-1), expr_copy(c1)}, 2));
    Expr* two_c0_b = expand_and_free(internal_times(
        (Expr*[]){expr_new_integer(2), expr_copy(c0)}, 2));

    Expr* u_num_p = eval_and_free(internal_plus(
        (Expr*[]){expr_copy(neg_c1), expr_copy(sqrt_Du)}, 2));
    Expr* u_num_m = eval_and_free(internal_subtract(
        (Expr*[]){neg_c1, sqrt_Du}, 2));
    Expr* u_plus = eval_and_free(internal_divide(
        (Expr*[]){u_num_p, expr_copy(two_c0_b)}, 2));
    Expr* u_minus = eval_and_free(internal_divide(
        (Expr*[]){u_num_m, two_c0_b}, 2));

    expr_free(c0); expr_free(c1); expr_free(c2);
    expr_free(c3); expr_free(c4);

    /* For each u, compute the conjugate-pair real contribution. */
    Expr* us[2] = {u_plus, u_minus};
    Expr* contributions[2] = {NULL, NULL};
    bool ok = true;

    for (int k = 0; k < 2 && ok; k++) {
        Expr* u = us[k];

        /* v_sq = (4 - u^2)/4.  Must be > 0 for v_im real. */
        Expr* u_sq = expand_and_free(internal_power(
            (Expr*[]){expr_copy(u), expr_new_integer(2)}, 2));
        Expr* neg_u_sq = expand_and_free(internal_times(
            (Expr*[]){expr_new_integer(-1), u_sq}, 2));
        Expr* four_m_usq = eval_and_free(internal_plus(
            (Expr*[]){expr_new_integer(4), neg_u_sq}, 2));
        Expr* v_sq = eval_and_free(internal_divide(
            (Expr*[]){four_m_usq, expr_new_integer(4)}, 2));
        Expr* v_sq_exp = expr_expand(v_sq); expr_free(v_sq);

        int v_sq_sign = intsimp_sign_pos_assumption(v_sq_exp);
        if (v_sq_sign == 0) v_sq_sign = intsimp_numeric_sign(v_sq_exp);
        if (v_sq_sign <= 0) {
            expr_free(v_sq_exp);
            ok = false; break;
        }

        Expr* v_im = intsimp_pos_sqrt(v_sq_exp);

        Expr* u_half = eval_and_free(internal_divide(
            (Expr*[]){expr_copy(u), expr_new_integer(2)}, 2));

        /* α = u_half + I · v_im. */
        Expr* I_unit = expr_new_function(expr_new_symbol(SYM_Complex),
            (Expr*[]){expr_new_integer(0), expr_new_integer(1)}, 2);
        Expr* I_v = eval_and_free(internal_times(
            (Expr*[]){I_unit, expr_copy(v_im)}, 2));
        Expr* alpha = eval_and_free(internal_plus(
            (Expr*[]){expr_copy(u_half), I_v}, 2));

        /* a(α) and d'(α) — split each into real + imag. */
        Expr* a_at_alpha = subst_t(a_t, t, alpha);
        Expr* dp_at_alpha = subst_t(dd_t, t, alpha);
        expr_free(alpha);

        Expr *Aa = NULL, *Ba = NULL;
        split_re_im(a_at_alpha, &Aa, &Ba);
        expr_free(a_at_alpha);
        Expr *Ad = NULL, *Bd = NULL;
        split_re_im(dp_at_alpha, &Ad, &Bd);
        expr_free(dp_at_alpha);

        /* c = (Aa + i Ba) / (Ad + i Bd)
         *   = (Aa Ad + Ba Bd + i(Ba Ad − Aa Bd)) / (Ad^2 + Bd^2). */
        Expr* AaAd = expand_and_free(internal_times(
            (Expr*[]){expr_copy(Aa), expr_copy(Ad)}, 2));
        Expr* BaBd = expand_and_free(internal_times(
            (Expr*[]){expr_copy(Ba), expr_copy(Bd)}, 2));
        Expr* re_num = eval_and_free(internal_plus(
            (Expr*[]){AaAd, BaBd}, 2));
        Expr* BaAd = expand_and_free(internal_times(
            (Expr*[]){expr_copy(Ba), expr_copy(Ad)}, 2));
        Expr* AaBd = expand_and_free(internal_times(
            (Expr*[]){expr_copy(Aa), expr_copy(Bd)}, 2));
        Expr* im_num = eval_and_free(internal_subtract(
            (Expr*[]){BaAd, AaBd}, 2));
        Expr* Ad_sq = expand_and_free(internal_power(
            (Expr*[]){Ad, expr_new_integer(2)}, 2));
        Expr* Bd_sq = expand_and_free(internal_power(
            (Expr*[]){Bd, expr_new_integer(2)}, 2));
        Expr* denom = eval_and_free(internal_plus(
            (Expr*[]){Ad_sq, Bd_sq}, 2));
        Expr* C_re = eval_and_free(internal_divide(
            (Expr*[]){re_num, expr_copy(denom)}, 2));
        Expr* C_im = eval_and_free(internal_divide(
            (Expr*[]){im_num, denom}, 2));
        expr_free(Aa); expr_free(Ba);

        /* L = (x − u_half)^2 + v_sq.
         * T = v_im / (x − u_half). */
        Expr* neg_uh = expand_and_free(internal_times(
            (Expr*[]){expr_new_integer(-1), expr_copy(u_half)}, 2));
        Expr* x_m_uh = eval_and_free(internal_plus(
            (Expr*[]){expr_copy(x), neg_uh}, 2));
        Expr* x_m_uh_sq = expand_and_free(internal_power(
            (Expr*[]){expr_copy(x_m_uh), expr_new_integer(2)}, 2));
        Expr* L = eval_and_free(internal_plus(
            (Expr*[]){x_m_uh_sq, v_sq_exp}, 2));
        Expr* T = eval_and_free(internal_divide(
            (Expr*[]){expr_copy(v_im), expr_copy(x_m_uh)}, 2));
        expr_free(x_m_uh); expr_free(v_im); expr_free(u_half);

        Expr* log_L = expr_new_function(expr_new_symbol(SYM_Log),
            (Expr*[]){L}, 1);
        Expr* atan_T = expr_new_function(expr_new_symbol(SYM_ArcTan),
            (Expr*[]){T}, 1);

        Expr* term_log = eval_and_free(internal_times(
            (Expr*[]){expr_copy(C_re), log_L}, 2));
        Expr* term_atan = eval_and_free(internal_times(
            (Expr*[]){expr_new_integer(2),
                      expr_copy(C_im), atan_T}, 3));
        expr_free(C_re); expr_free(C_im);

        contributions[k] = eval_and_free(internal_plus(
            (Expr*[]){term_log, term_atan}, 2));
    }

    expr_free(u_plus); expr_free(u_minus);

    if (!ok) {
        if (contributions[0]) expr_free(contributions[0]);
        if (contributions[1]) expr_free(contributions[1]);
        return NULL;
    }

    Expr* total = eval_and_free(internal_plus(contributions, 2));
    return total;
}

/* When `poly` (a univariate polynomial in `bvar`) is one of the
 * structural shapes whose roots have a closed-form radical formula,
 * expand `RootSum[Function[bvar, poly], Function[bvar, body]]` into
 * an explicit Plus[body[α₁], body[α₂], ...] and return it.  Returns
 * NULL when the polynomial does not match any of the supported
 * shapes — currently:
 *
 *   degree 1                                         (1 root)
 *   degree 2                                         (2 roots)
 *   degree 4 with only constant / x^2 / x^4 terms    (4 roots)
 *
 * Gated on poly being free of parameters other than bvar — see
 * poly_only_uses above for the rationale.
 *
 * Once Solve / ToRadicals lands, this routine is superseded by a
 * call to Solve[poly == 0, bvar, Reals] followed by the same body
 * substitution loop. */
static Expr* expand_simple_rootsum(Expr* poly, Expr* bvar, Expr* body) {
    if (!poly_only_uses(poly, bvar)) return NULL;

    Expr* cl_call = expr_new_function(expr_new_symbol(SYM_CoefficientList),
        (Expr*[]){ expr_copy(poly), expr_copy(bvar) }, 2);
    Expr* cl = evaluate(cl_call);
    expr_free(cl_call);
    if (!cl || cl->type != EXPR_FUNCTION
        || cl->data.function.head->type != EXPR_SYMBOL
        || cl->data.function.head->data.symbol.name != SYM_List) {
        if (cl) expr_free(cl);
        return NULL;
    }

    size_t n = cl->data.function.arg_count;
    Expr* result = NULL;

    if (n == 2) {
        /* c0 + c1 t = 0  =>  t = -c0 / c1. */
        Expr* c0 = expr_copy(cl->data.function.args[0]);
        Expr* c1 = expr_copy(cl->data.function.args[1]);
        Expr* root_raw = internal_divide(
            (Expr*[]){ nlp_neg(c0), c1 }, 2);
        Expr* root = eval_and_free(root_raw);
        result = nlp_sub_body(body, bvar, root);
        expr_free(root);
    } else if (n == 3) {
        /* c0 + c1 t + c2 t^2 = 0
         * t = (-c1 ± Sqrt[c1² - 4 c0 c2]) / (2 c2). */
        Expr* c0 = expr_copy(cl->data.function.args[0]);
        Expr* c1 = expr_copy(cl->data.function.args[1]);
        Expr* c2 = expr_copy(cl->data.function.args[2]);

        Expr* disc_raw = internal_plus((Expr*[]){
            internal_times((Expr*[]){ expr_copy(c1), expr_copy(c1) }, 2),
            nlp_neg(internal_times((Expr*[]){
                expr_new_integer(4),
                expr_copy(c0), expr_copy(c2) }, 3))
        }, 2);
        Expr* disc = eval_and_free(disc_raw);
        Expr* sqd  = eval_and_free(nlp_sqrt(disc));

        Expr* two_c2 = eval_and_free(internal_times(
            (Expr*[]){ expr_new_integer(2), expr_copy(c2) }, 2));

        Expr* num_p = eval_and_free(internal_plus((Expr*[]){
            nlp_neg(expr_copy(c1)), expr_copy(sqd) }, 2));
        Expr* num_m = eval_and_free(internal_plus((Expr*[]){
            nlp_neg(expr_copy(c1)), nlp_neg(expr_copy(sqd)) }, 2));
        expr_free(sqd); expr_free(c0); expr_free(c1); expr_free(c2);

        Expr* root_p = eval_and_free(internal_divide(
            (Expr*[]){ num_p, expr_copy(two_c2) }, 2));
        Expr* root_m = eval_and_free(internal_divide(
            (Expr*[]){ num_m, two_c2 }, 2));

        Expr* t1 = nlp_sub_body(body, bvar, root_p);
        Expr* t2 = nlp_sub_body(body, bvar, root_m);
        expr_free(root_p); expr_free(root_m);

        result = eval_and_free(internal_plus((Expr*[]){ t1, t2 }, 2));
    } else if (n == 5) {
        /* Biquadratic: only c0, c2, c4 are non-zero.
         *   c0 + c2 t² + c4 t⁴ = 0
         *   u = (-c2 ± Sqrt[c2² - 4 c0 c4]) / (2 c4),  t = ±Sqrt[u]. */
        Expr* c1 = cl->data.function.args[1];
        Expr* c3 = cl->data.function.args[3];
        bool zero_c1 = (c1->type == EXPR_INTEGER && c1->data.integer == 0);
        bool zero_c3 = (c3->type == EXPR_INTEGER && c3->data.integer == 0);
        if (zero_c1 && zero_c3) {
            Expr* c0 = expr_copy(cl->data.function.args[0]);
            Expr* c2 = expr_copy(cl->data.function.args[2]);
            Expr* c4 = expr_copy(cl->data.function.args[4]);

            Expr* disc_raw = internal_plus((Expr*[]){
                internal_times((Expr*[]){ expr_copy(c2), expr_copy(c2) }, 2),
                nlp_neg(internal_times((Expr*[]){
                    expr_new_integer(4),
                    expr_copy(c0), expr_copy(c4) }, 3))
            }, 2);
            Expr* disc = eval_and_free(disc_raw);
            Expr* sqd  = eval_and_free(nlp_sqrt(disc));

            Expr* two_c4 = eval_and_free(internal_times(
                (Expr*[]){ expr_new_integer(2), expr_copy(c4) }, 2));

            Expr* up_num = eval_and_free(internal_plus((Expr*[]){
                nlp_neg(expr_copy(c2)), expr_copy(sqd) }, 2));
            Expr* um_num = eval_and_free(internal_plus((Expr*[]){
                nlp_neg(expr_copy(c2)), nlp_neg(expr_copy(sqd)) }, 2));
            expr_free(sqd); expr_free(c0); expr_free(c2); expr_free(c4);

            Expr* u_p = eval_and_free(internal_divide(
                (Expr*[]){ up_num, expr_copy(two_c4) }, 2));
            Expr* u_m = eval_and_free(internal_divide(
                (Expr*[]){ um_num, two_c4 }, 2));

            Expr* sp = eval_and_free(nlp_sqrt(u_p));
            Expr* sm = eval_and_free(nlp_sqrt(u_m));

            Expr* t1 = nlp_sub_body(body, bvar, sp);
            Expr* t2 = nlp_sub_body(body, bvar, eval_and_free(nlp_neg(expr_copy(sp))));
            Expr* t3 = nlp_sub_body(body, bvar, sm);
            Expr* t4 = nlp_sub_body(body, bvar, eval_and_free(nlp_neg(expr_copy(sm))));
            expr_free(sp); expr_free(sm);

            result = eval_and_free(internal_plus(
                (Expr*[]){ t1, t2, t3, t4 }, 4));
        }
    }

    expr_free(cl);
    return result;
}

/* ------------------------------------------------------------------ */
/* NaiveLogPart — RootSum fallback for the log part.                   */
/* ------------------------------------------------------------------ */

/* Direct port of IntegrateRational.m:1116-1124.
 *
 *   NaiveLogPart[f, x] = Σ_{α : d(α) = 0}  a(α) Log(x - α) / d'(α)
 *
 * for f = a/d.  Returned in held-symbolic RootSum form
 *
 *   RootSum[ Function[t, d(t)],
 *            Function[t, a(t) Log(x - t) / d'(t)] ]
 *
 * with `t` a fresh bound variable distinct from the outer integration
 * variable x.  The result is the universal fallback when LogToReal
 * cannot close the log part to a real elementary expression — every
 * proper rational integrand admits this representation, and the
 * D[RootSum] rule wired in src/deriv.c differentiates it back to f.
 *
 * This routine deliberately does NOT call ToRadicals on the
 * biquadratic case (cf. IntegrateRational.m:1121).  Algebraic-extension
 * radical closure is the responsibility of Solve / ToRadicals which we
 * have not yet implemented; until then RootSum is the honest answer. */
static Expr* intrat_naive_log_part(Expr* f, Expr* x) {
    intrat_trace("NaiveLogPart", "IN", f);

    /* a/d split.  We canonicalise f via Together so a and d are
     * polynomial in x even when f arrived as a sum of partial-fraction
     * pieces or with rational coefficients. */
    Expr* tg = internal_together((Expr*[]){expr_copy(f)}, 1);
    Expr* combined = evaluate(tg);
    expr_free(tg);
    if (!combined) return NULL;

    Expr* a = eval_and_free(intrat_numerator(combined));
    Expr* d = eval_and_free(intrat_denominator(combined));
    expr_free(combined);

    /* d'(x), expanded so polynomial divisions downstream see it in
     * normal form. */
    Expr* dprime = intrat_d(d, x);
    Expr* dd = expr_expand(dprime); expr_free(dprime);

    /* Fresh bound variable.  The Private subcontext keeps it well clear
     * of any user symbol the caller might have set globally — using a
     * bare `t` would be unsafe because ReplaceAll evaluates the rule's
     * RHS, so a user-bound t = 5 would corrupt the substitution. */
    Expr* t = expr_new_symbol("Integrate`Private`t");

    /* Substitute t for x in a, d, dd via ReplaceAll. */
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){ expr_copy(x), expr_copy(t) }, 2);

    Expr* a_t = eval_and_free(internal_replace_all(
        (Expr*[]){ expr_copy(a), expr_copy(rule) }, 2));
    Expr* d_t = eval_and_free(internal_replace_all(
        (Expr*[]){ expr_copy(d), expr_copy(rule) }, 2));
    Expr* dd_t = eval_and_free(internal_replace_all(
        (Expr*[]){ expr_copy(dd), expr_copy(rule) }, 2));
    expr_free(rule); expr_free(a); expr_free(d); expr_free(dd);

    /* Phase D2 — palindromic-quartic real-form expansion.  Builds the
     * conjugate-pair contributions directly, sidestepping the LogToReal
     * / LogToAtan chain that overwhelms on nested-radical S
     * polynomials.  Try this before constructing the held RootSum body
     * so we never produce a Function-headed RootSum for inputs whose
     * denominator is a palindromic quartic. */
    {
        Expr* palindromic_real = expand_palindromic_quartic_real(
            a_t, d_t, dd_t, t, x);
        if (palindromic_real) {
            expr_free(t); expr_free(a_t); expr_free(d_t); expr_free(dd_t);
            intrat_trace("NaiveLogPart", "OUT (palindromic-real)",
                         palindromic_real);
            return palindromic_real;
        }
    }

    /* Build the body: a_t * Log[x - t] / dd_t. */
    Expr* neg_t = internal_times((Expr*[]){
        expr_new_integer(-1),
        expr_copy(t)
    }, 2);
    Expr* x_minus_t = internal_plus(
        (Expr*[]){ expr_copy(x), neg_t }, 2);
    Expr* x_minus_t_e = eval_and_free(x_minus_t);

    Expr* log_term = expr_new_function(expr_new_symbol(SYM_Log),
        (Expr*[]){ x_minus_t_e }, 1);

    Expr* numer = internal_times(
        (Expr*[]){ a_t, log_term }, 2);
    Expr* numer_e = eval_and_free(numer);

    Expr* body = internal_divide(
        (Expr*[]){ numer_e, dd_t }, 2);
    Expr* body_e = eval_and_free(body);

    /* Phase 8d-bonus finishing pass: if the polynomial d(t) has a
     * closed-form radical-formula root set (linear, quadratic, or
     * biquadratic), expand the RootSum to an explicit Plus over the
     * 1/2/4 roots and evaluate the body at each.  Otherwise fall
     * back to the held RootSum form.  Once Solve / ToRadicals lands
     * the simple-shape detector is replaced by Solve. */
    Expr* expanded = expand_simple_rootsum(d_t, t, body_e);
    if (expanded) {
        expr_free(t); expr_free(d_t); expr_free(body_e);
        intrat_trace("NaiveLogPart", "OUT (expanded)", expanded);
        return expanded;
    }

    /* Wrap in RootSum[Function[t, d_t], Function[t, body_e]]. */
    Expr* result = root_make_rootsum(t, d_t, body_e);
    intrat_trace("NaiveLogPart", "OUT", result);
    return result;
}

/* Divide p by the polynomial GCD of its NONZERO coefficients viewed as a
 * polynomial in cvar (the content of p in cvar).  internal PolynomialGCD
 * bails to a held call when handed a zero argument (PolynomialGCD[x^2, 0,
 * 4 x^2] returns 1, not x^2), so zero coefficients are filtered first.
 * Returns Cancel[p / content] when the content is a non-trivial polynomial,
 * else a copy of p.  Borrows p; caller owns the result.  Used by the
 * transcendental log-part path to strip the x-only content that the
 * denominator-clearing scaling introduces into the Rothstein-Trager
 * resultant, exposing the (constant-residue) polynomial in the residue
 * variable. */
static Expr* intrat_strip_content_var(Expr* p, Expr* cvar) {
    Expr* cl = eval_and_free(internal_coefficientlist(
        (Expr*[]){expr_copy(p), expr_copy(cvar)}, 2));
    if (!cl || cl->type != EXPR_FUNCTION
        || cl->data.function.head->type != EXPR_SYMBOL
        || cl->data.function.head->data.symbol.name != SYM_List
        || cl->data.function.arg_count == 0) {
        if (cl) expr_free(cl);
        return expr_copy(p);
    }
    size_t n = cl->data.function.arg_count;
    Expr** nz = (Expr**)malloc(sizeof(Expr*) * n);
    size_t nnz = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* c = cl->data.function.args[i];
        if (!is_zero_poly(c)) nz[nnz++] = expr_copy(c);
    }
    expr_free(cl);
    if (nnz == 0) { free(nz); return expr_copy(p); }
    Expr* content = (nnz == 1) ? nz[0] : internal_polynomialgcd(nz, nnz);
    free(nz);
    if (head_is(content, SYM_PolynomialGCD)
        || (content->type == EXPR_INTEGER
            && (content->data.integer == 1 || content->data.integer == -1))) {
        expr_free(content);
        return expr_copy(p);
    }
    Expr* div = internal_divide((Expr*[]){expr_copy(p), content}, 2);
    return eval_and_free(internal_cancel((Expr*[]){div}, 1));
}

/* True iff polynomial p is free of every variable named in `xgate` (a single
 * symbol, or a List of symbols) — i.e. the roots of p are constants of the
 * derivation.  Used by the residue criterion's κ_D split. */
static bool intrat_freeq_all(Expr* p, Expr* xgate) {
    if (xgate->type == EXPR_FUNCTION
        && xgate->data.function.head->type == EXPR_SYMBOL
        && xgate->data.function.head->data.symbol.name == SYM_List) {
        for (size_t gi = 0; gi < xgate->data.function.arg_count; gi++)
            if (!intrat_freeq_test(p, xgate->data.function.args[gi]))
                return false;
        return true;
    }
    return intrat_freeq_test(p, xgate);
}

/* Recursively walk a Factor[] result, routing each irreducible factor (base of
 * a Power, or a bare polynomial in `t`) into the constant-residue product *Qs
 * (free of the gate vars) or the non-constant product *Qn.  Numeric / t-free
 * factors are dropped.  *Qs, *Qn are running products (owned), each starting at
 * integer 1; internal_times adopts them, so they are consumed and replaced. */
static void intrat_route_factor(Expr* f, Expr* t, Expr* xgate,
                                Expr** Qs, Expr** Qn) {
    if (f->type == EXPR_FUNCTION
        && f->data.function.head->type == EXPR_SYMBOL) {
        const char* h = f->data.function.head->data.symbol.name;
        if (h == SYM_Times) {
            for (size_t k = 0; k < f->data.function.arg_count; k++)
                intrat_route_factor(f->data.function.args[k], t, xgate, Qs, Qn);
            return;
        }
        if (h == SYM_Power) {   /* squarefree bucket => exponent 1: route base */
            intrat_route_factor(f->data.function.args[0], t, xgate, Qs, Qn);
            return;
        }
    }
    if (get_degree_poly(f, t) <= 0) return;   /* numeric / t-free: drop */
    Expr** dst = intrat_freeq_all(f, xgate) ? Qs : Qn;
    *dst = expand_and_free(internal_times((Expr*[]){*dst, expr_copy(f)}, 2));
}

/* Split a squarefree bucket polynomial Qi (in `t`) into its constant-residue
 * sub-product (returned) and, multiplied into *Qn_acc, its non-constant-residue
 * sub-product.  Returns the constant product (owned; integer 1 if none). */
static Expr* intrat_bucket_const_part(Expr* Qi, Expr* t, Expr* xgate,
                                      Expr** Qn_acc) {
    Expr* factored = eval_and_free(internal_factor((Expr*[]){expr_copy(Qi)}, 1));
    Expr* Qs = expr_new_integer(1);
    if (factored) { intrat_route_factor(factored, t, xgate, &Qs, Qn_acc); expr_free(factored); }
    return Qs;
}

/* IntRationalLogPart core. Returns either the {Q, S} pair list
 * (root_sum=false) or the symbolic RootSum sum (root_sum=true).
 * On algorithmic failure returns NULL.
 *
 * `a`, `d` are the (owned) numerator / denominator in the main variable `x`
 * with squarefree `d`; `t` is the residue variable.  Two optional hooks
 * generalise the pure-rational computation to a transcendental monomial
 * extension (the recursive-Risch Lazard-Rioboo-Trager log part):
 *   - `dprime_override != NULL` injects the monomial derivation `D(d)` in
 *     place of the ordinary `d/dx` (here `x` is the monomial variable);
 *   - `xgate != NULL` requires the residues (roots of the resultant) to be
 *     constants of the derivation (free of `xgate`).  Bronstein's residue
 *     criterion (Thm 5.6.1) splits the resultant r = r_s · r_n by whether each
 *     factor's residues are constant: the constant part r_s yields elementary
 *     logs, the non-constant part r_n does not.  When `remainder_out != NULL`
 *     and the split is genuinely mixed, the elementary r_s logs are returned AND
 *     *remainder_out is set to the (owned) r_n rational-function remainder to be
 *     reported unintegrated (partial log part).  When `remainder_out == NULL`
 *     (caller refuses partials) a mixed resultant declines (NULL), as does the
 *     all-non-constant case; the all-constant case is unchanged.
 * With both NULL this is the original pure-rational log part.
 *
 * Direct port of IntegrateRational.m:751-801. */
static Expr* intrat_log_part_core(Expr* a, Expr* d, Expr* x, Expr* t,
                                  bool root_sum, Expr* dprime_override,
                                  Expr* xgate, Expr** remainder_out, bool decide) {
    if (remainder_out) *remainder_out = NULL;
    /* prs = SubresultantPolynomialRemainders[d, a - t*D(d), x] */
    Expr* dprime_e;
    if (dprime_override) {
        dprime_e = expr_expand(dprime_override);
    } else {
        Expr* dprime = intrat_d(d, x);
        dprime_e = expr_expand(dprime); expr_free(dprime);
    }
    Expr* tdprime = internal_times(
        (Expr*[]){expr_copy(t), dprime_e}, 2);
    Expr* atdp_raw = internal_subtract(
        (Expr*[]){expr_copy(a), tdprime}, 2);
    Expr* atdp;
    if (dprime_override) {
        /* the monomial derivation D(d) can carry x-denominators (e.g. the
         * u'/u of a logarithmic kernel); clear them so the resultant / PRS
         * operate on genuine polynomials.  The x-only denominator scales the
         * resultant by a unit that the x-content strip below removes. */
        Expr* tg = eval_and_free(internal_together((Expr*[]){atdp_raw}, 1));
        Expr* nn = intrat_numerator(tg);
        atdp = expr_expand(nn); expr_free(nn); expr_free(tg);
    } else {
        atdp = expr_expand(atdp_raw); expr_free(atdp_raw);
    }

    Expr* prs_args[3] = { expr_copy(d), expr_copy(atdp), expr_copy(x) };
    Expr* prs_call = expr_new_function(
        expr_new_symbol(SYM_SubresultantPolynomialRemainders), prs_args, 3);
    Expr* prs = evaluate(prs_call);
    expr_free(prs_call);
    if (!prs) {
        expr_free(a); expr_free(d); expr_free(atdp); return NULL;
    }

    /* resultant = primitive[Resultant[d, a - t*D[d,x], x], t]
     *
     * NOTE: this Resultant[] call is NOT redundant with the PRS chain
     * computed above, despite both being subresultant-PRS-flavoured
     * computations.  Mathilda's SubresultantPolynomialRemainders is a
     * pseudo-remainder chain (poly.c:1637-1642), not the Lazard-scaled
     * subresultant chain.  The chain's primitive-part-in-t and
     * degree-in-x are content-invariant and so are correct substrates
     * for the S_i lookups below (i >= 1); but the chain's degree-0-in-x
     * entry is NOT in general the resultant, even after primitive in
     * t.  Concrete witness: x/(1 + x^8) has resultant (1 + 4096 t^4)^2
     * but PRS-bottom primitive in t is t (1 + 4096 t^4), structurally
     * different. */
    /* The residue map is a/D(d) = P(x)/W(x) where atdp = P - t*W is linear in t.
     * Snapshot it before the Resultant call adopts `atdp` — the κ_D split's
     * remainder reconstruction resultant needs it (partial log part). */
    Expr* atdp_keep = xgate ? expr_copy(atdp) : NULL;
    Expr* res_args[3] = { expr_copy(d), atdp, expr_copy(x) };
    Expr* res_call = expr_new_function(
        expr_new_symbol(SYM_Resultant), res_args, 3);
    Expr* resultant_raw = evaluate(res_call);
    expr_free(res_call);
    if (!resultant_raw) {
        expr_free(a); expr_free(d); expr_free(prs);
        if (atdp_keep) expr_free(atdp_keep);
        return NULL;
    }
    Expr* resultant = intrat_primitive(resultant_raw, t);
    expr_free(resultant_raw);

    /* κ_D split state (partial log part).  In `partial_mode` the resultant has a
     * genuine mix of constant- and non-constant-residue factors: emit only the
     * constant-residue (r_s) logs and report the r_n part via *remainder_out. */
    bool partial_mode = false;
    Expr* remainder_local = NULL;

    if (xgate) {
        /* Strip the lower-field content the denominator-clearing scaling leaves,
         * then classify the residues (roots of the resultant) by whether they
         * are constants of the derivation (free of every gate variable).  A
         * single symbol gates a single-kernel extension; a `List` gates a tower
         * proper part.  Bronstein's residue criterion (Thm 5.6.1) keeps the
         * elementary constant-residue logs even when non-constant residues are
         * present; the fully-constant case is unchanged. */
        Expr* stripped = intrat_strip_content_var(resultant, t);
        expr_free(resultant);
        resultant = stripped;
        bool free_all = intrat_freeq_all(resultant, xgate);
        if (!free_all && decide) {
            /* Decision mode: ANY non-constant residue proves this simple part has no
             * elementary integral (Bronstein residue criterion, Thm 5.6.1(ii)).  We
             * only need the VERDICT, not the elementary constant-residue logs or a
             * reconstructed remainder — return an inert marker the recursive-Risch
             * caller (rt_field_lrt_logpart) decodes into a NONELEMENTARY verdict. */
            expr_free(resultant); expr_free(a); expr_free(d); expr_free(prs);
            if (atdp_keep) expr_free(atdp_keep);
            return expr_new_symbol("Integrate`$NonConstantResidue");
        }
        if (!free_all) {
            /* Mixed or fully non-constant.  A caller that cannot represent a
             * partial result (remainder_out == NULL — e.g. the exponential
             * coupled path) declines exactly as before.  Otherwise defer to the
             * per-factor partition after the squarefree split below. */
            if (!remainder_out) {
                expr_free(resultant); expr_free(a); expr_free(d); expr_free(prs);
                if (atdp_keep) expr_free(atdp_keep);
                return NULL;
            }
            partial_mode = true;
        }
    }

    /* Q = SquareFree[resultant]: list of {poly_i, mult_i}. */
    Expr* Q = intrat_squarefree_list(resultant);
    expr_free(resultant);
    if (!Q) {
        expr_free(a); expr_free(d); expr_free(prs);
        if (atdp_keep) expr_free(atdp_keep);
        return NULL;
    }

    int deg_d = get_degree_poly(d, x);
    size_t nQ = Q->data.function.arg_count;

    /* Partial κ_D partition: per bucket, the constant-residue sub-product
     * Qconst[i] (used in the output instead of the full bucket) and the global
     * non-constant product Q_n_poly (whose residues drive the remainder).  Only
     * built in partial_mode; NULL entries mean "no constant residues in bucket". */
    Expr** Qconst = NULL;
    if (partial_mode) {
        Qconst = (Expr**)malloc(sizeof(Expr*) * nQ);
        Expr* Q_n_poly = expr_new_integer(1);
        bool any_const = false;
        for (size_t i_idx = 0; i_idx < nQ; i_idx++) {
            Expr* Qi = list_get(Q->data.function.args[i_idx], 1);
            if (!Qi || get_degree_poly(Qi, t) <= 0) { if (Qi) expr_free(Qi); Qconst[i_idx] = NULL; continue; }
            Qconst[i_idx] = intrat_bucket_const_part(Qi, t, xgate, &Q_n_poly);
            expr_free(Qi);
            if (Qconst[i_idx] && get_degree_poly(Qconst[i_idx], t) > 0) any_const = true;
        }
        if (!any_const) {
            /* No elementary logs to emit — decline (also prevents an infinite
             * Integrate[remainder] regress).  Decision mode never reaches here: it
             * short-circuits at the non-constant-residue detection above. */
            for (size_t i_idx = 0; i_idx < nQ; i_idx++) if (Qconst[i_idx]) expr_free(Qconst[i_idx]);
            free(Qconst); expr_free(Q_n_poly); expr_free(Q);
            expr_free(a); expr_free(d); expr_free(prs);
            if (atdp_keep) expr_free(atdp_keep);
            return NULL;
        }
        /* Remainder = the d_n-part of a/d, where d_n = gcd(d, Res_t(Q_n, atdp))
         * collects the poles whose residue is a root of Q_n (non-constant).
         * atdp = P - t*W is linear in t, so Res_t(Q_n, atdp) = W^deg * Q_n(P/W)
         * (a poly in x); its gcd with d picks exactly those poles. */
        Expr* resn_call = expr_new_function(expr_new_symbol(SYM_Resultant),
            (Expr*[]){ Q_n_poly, atdp_keep, expr_copy(t) }, 3);   /* adopts both */
        atdp_keep = NULL;
        Expr* Res_n = eval_and_free(resn_call);
        Expr* d_n = Res_n ? intrat_polygcd_monic(d, Res_n, x) : NULL;
        if (Res_n) expr_free(Res_n);
        Expr* d_s = d_n ? intrat_exquo(d, d_n, x) : NULL;
        Expr* rr = NULL; Expr* ss = NULL;
        bool ok = d_n && d_s
            && intrat_extended_euclidean(d_n, d_s, a, x, &rr, &ss);  /* rr d_n + ss d_s = a */
        if (rr) expr_free(rr);
        if (ok) {
            /* a/(d_n d_s) = rr/d_s + ss/d_n; the non-constant (d_n) part is ss/d_n. */
            remainder_local = eval_and_free(internal_together((Expr*[]){
                internal_divide((Expr*[]){ ss, expr_copy(d_n) }, 2) }, 1));
        } else if (ss) { expr_free(ss); }
        if (d_n) expr_free(d_n);
        if (d_s) expr_free(d_s);
        if (!remainder_local) {
            /* Could not reconstruct the remainder — decline rather than emit a
             * partial result whose residual identity we cannot certify. */
            for (size_t i_idx = 0; i_idx < nQ; i_idx++) if (Qconst[i_idx]) expr_free(Qconst[i_idx]);
            free(Qconst); expr_free(Q);
            expr_free(a); expr_free(d); expr_free(prs);
            return NULL;
        }
    }
    if (atdp_keep) { expr_free(atdp_keep); atdp_keep = NULL; }

    /* Compute S[i] for each i.  Where Q[[i,1]] has non-trivial t-
     * degree, we derive a polynomial in (t, x). */
    Expr** S = (Expr**)malloc(sizeof(Expr*) * nQ);
    for (size_t i = 0; i < nQ; i++) S[i] = NULL;

    for (size_t i_idx = 0; i_idx < nQ; i_idx++) {
        int i = (int)(i_idx + 1);
        Expr* Qi_pair = Q->data.function.args[i_idx];
        Expr* Qi = list_get(Qi_pair, 1);
        if (!Qi) continue;

        if (get_degree_poly(Qi, t) <= 0) {
            expr_free(Qi);
            continue;
        }

        if (i == deg_d) {
            /* Special case (Geddes-Czapor-Labahn correction
             * IntegrateRational.m:765-768): use d itself with lc[d, x]
             * as the leading-coefficient witness in the Diophantine. */
            S[i_idx] = intrat_primitive_part_mod(d, Qi, x, t);
        } else {
            Expr* s_raw = find_prs_at_degree(prs, x, i);
            if (!s_raw) { expr_free(Qi); continue; }
            Expr* s = intrat_primitive(s_raw, x);
            expr_free(s_raw);

            /* A = SquareFree[primitive[lc[s, x], t]]. */
            Expr* lc_s = intrat_lc(s, x);
            Expr* lc_s_prim = intrat_primitive(lc_s, t);
            expr_free(lc_s);
            Expr* A = intrat_squarefree_list(lc_s_prim);
            expr_free(lc_s_prim);
            if (!A) { expr_free(Qi); expr_free(s); continue; }

            /* Divide common factors between A_j and Qi out of s. */
            for (size_t j = 0; j < A->data.function.arg_count; j++) {
                Expr* A_j_pair = A->data.function.args[j];
                Expr* A_jpoly = list_get(A_j_pair, 1);
                if (!A_jpoly) continue;
                Expr* g_args[3] = {A_jpoly, expr_copy(Qi), expr_copy(t)};
                /* Use PolynomialExtendedGCD-derived monic gcd to
                 * avoid the content-scaling issue we hit in Phase 1. */
                Expr* g = intrat_polygcd_monic(g_args[0], g_args[1], g_args[2]);
                expr_free(A_jpoly); expr_free(g_args[1]); expr_free(g_args[2]);
                if (!g || (g->type == EXPR_INTEGER && g->data.integer == 1)) {
                    if (g) expr_free(g);
                    continue;
                }
                /* s := exquo[s, g^j, t] (j is 1-based index = j+1). */
                int jpow = (int)(j + 1);
                Expr* gj = expand_and_free(internal_power(
                    (Expr*[]){g, expr_new_integer((int64_t)jpow)}, 2));
                Expr* sn = intrat_exquo(s, gj, t);
                expr_free(gj);
                if (sn) { expr_free(s); s = sn; }
            }
            expr_free(A);

            /* s := PrimitivePart_mod(s, Qi) -- see helper definition. */
            S[i_idx] = intrat_primitive_part_mod(s, Qi, x, t);
            expr_free(s);
        }
        expr_free(Qi);
    }
    expr_free(prs);
    expr_free(a); expr_free(d);

    /* Build output. */
    if (root_sum) {
        /* Sum of RootSum heads. */
        Expr** terms = (Expr**)malloc(sizeof(Expr*) * nQ);
        size_t nterms = 0;
        for (size_t i_idx = 0; i_idx < nQ; i_idx++) {
            if (!S[i_idx]) continue;
            /* In partial_mode restrict the RootSum to the constant-residue
             * sub-product of the bucket (S[i] mod Qi agrees on its sub-roots). */
            Expr* Qi;
            if (partial_mode) {
                if (!Qconst[i_idx] || get_degree_poly(Qconst[i_idx], t) <= 0) continue;
                Qi = expr_copy(Qconst[i_idx]);
            } else {
                Qi = list_get(Q->data.function.args[i_idx], 1);
                if (!Qi) continue;
                if (get_degree_poly(Qi, t) <= 0) { expr_free(Qi); continue; }
            }

            /* RootSum[Function[t, Qi], Function[t, t Log[S[i]]]].  Canonicalise
             * the held Function bodies up front (Function is HoldAll, so the
             * evaluator will not normalise them afterwards): sort the resultant
             * factor Qi and fold the t*Log[S] product so the RootSum prints in
             * normal form. */
            Qi = eval_and_free(Qi);
            Expr* func1 = expr_new_function(expr_new_symbol(SYM_Function),
                (Expr*[]){expr_copy(t), Qi}, 2);
            Expr* logS = expr_new_function(expr_new_symbol(SYM_Log),
                (Expr*[]){expr_copy(S[i_idx])}, 1);
            Expr* tlog = eval_and_free(internal_times(
                (Expr*[]){expr_copy(t), logS}, 2));
            Expr* func2 = expr_new_function(expr_new_symbol(SYM_Function),
                (Expr*[]){expr_copy(t), tlog}, 2);
            terms[nterms++] = expr_new_function(expr_new_symbol(SYM_RootSum),
                (Expr*[]){func1, func2}, 2);
        }
        for (size_t i_idx = 0; i_idx < nQ; i_idx++) if (S[i_idx]) expr_free(S[i_idx]);
        free(S);
        expr_free(Q);
        if (partial_mode) {
            for (size_t i_idx = 0; i_idx < nQ; i_idx++) if (Qconst[i_idx]) expr_free(Qconst[i_idx]);
            free(Qconst);
            *remainder_out = remainder_local;
        }
        Expr* result;
        if (nterms == 0) { free(terms); result = expr_new_integer(0); }
        else if (nterms == 1) { result = terms[0]; free(terms); }
        else { result = internal_plus(terms, nterms); free(terms); result = eval_and_free(result); }
        intrat_trace("IntRationalLogPart", "OUT", result);
        return result;
    }

    /* RootSum -> False: list of {Q_i, S_i} pairs. */
    Expr** pairs = (Expr**)malloc(sizeof(Expr*) * nQ);
    size_t npairs = 0;
    for (size_t i_idx = 0; i_idx < nQ; i_idx++) {
        if (!S[i_idx]) continue;
        Expr* Qi;
        if (partial_mode) {
            if (!Qconst[i_idx] || get_degree_poly(Qconst[i_idx], t) <= 0) continue;
            Qi = expr_copy(Qconst[i_idx]);
        } else {
            Qi = list_get(Q->data.function.args[i_idx], 1);
            if (!Qi) continue;
            if (get_degree_poly(Qi, t) <= 0) { expr_free(Qi); continue; }
        }
        pairs[npairs++] = expr_new_function(expr_new_symbol(SYM_List),
            (Expr*[]){Qi, expr_copy(S[i_idx])}, 2);
    }
    for (size_t i_idx = 0; i_idx < nQ; i_idx++) if (S[i_idx]) expr_free(S[i_idx]);
    free(S);
    expr_free(Q);
    if (partial_mode) {
        for (size_t i_idx = 0; i_idx < nQ; i_idx++) if (Qconst[i_idx]) expr_free(Qconst[i_idx]);
        free(Qconst);
        *remainder_out = remainder_local;
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), pairs, npairs);
    free(pairs);
    intrat_trace("IntRationalLogPart", "OUT", result);
    return result;
}

/* Pure-rational log part: Together-canonicalise f, split into a/d, and run
 * the core with the ordinary derivation (no transcendental hooks). */
static Expr* intrat_int_rational_log_part(Expr* f, Expr* x, Expr* t,
                                          bool root_sum) {
    intrat_trace("IntRationalLogPart", "IN", f);
    /* Together-canonicalise f so the numerator and denominator are clean
     * polynomials with no internal fractions in parametric coefficients. */
    Expr* fc = eval_and_free(internal_together((Expr*[]){expr_copy(f)}, 1));
    Expr* a = eval_and_free(intrat_numerator(fc));
    Expr* d_raw = intrat_denominator(fc);
    Expr* d = expr_expand(d_raw); expr_free(d_raw);
    expr_free(fc);
    return intrat_log_part_core(a, d, x, t, root_sum, NULL, NULL, NULL, false);
}

/* Defined below (Phase 4 consumer); forward-declared for the transcendental
 * log part which hands its {Q_i, S_i} pairs straight to Rioboo's LogToReal. */
static Expr* intrat_log_to_real_pairs(Expr* pair_list, Expr* x, Expr* t);

/* Transcendental (recursive-Risch) Lazard-Rioboo-Trager log part.
 *
 * Integrates the PROPER, SQUAREFREE-denominator fraction  a(tau)/d(tau)  of a
 * single transcendental monomial extension tau (= Log[u] or E^u over Q(x)),
 * whose log part carries ALGEBRAIC (non-rational) residues — the class the
 * single-constant-per-factor SolveAlways ansatz cannot express.  `Dd` is the
 * monomial derivation D(d) = d/dx(d) + (D tau) d/dtau(d), computed by the
 * caller.  The residues are the roots of the Rothstein-Trager resultant
 * Res_tau(a - z Dd, d) (z = the residue variable); each squarefree factor
 * Q_i(z) with its log argument S_i(z, tau) is converted to real Log + ArcTan
 * form by Rioboo's LogToReal.  Returns the antiderivative as a function of
 * `tau` (the caller substitutes tau -> kernel), or NULL when the integral is
 * not elementary in this form / a factor exceeds LogToReal's bounded scope.
 * `a`, `d`, `Dd`, `tau`, `z`, `xreal` are all borrowed.
 *
 * When `remainder_out != NULL`, the residue criterion may return a PARTIAL log
 * part (Bronstein Thm 5.6.1): the elementary constant-residue logs, with the
 * non-constant-residue rational remainder written to *remainder_out (owned, a
 * function of `tau`) for the caller to report unintegrated.  *remainder_out is
 * NULL when the whole log part is elementary. */
static Expr* intrat_transcendental_log_part(Expr* a, Expr* d, Expr* tau,
                                            Expr* z, Expr* Dd, Expr* xreal,
                                            Expr** remainder_out, bool decide) {
    if (remainder_out) *remainder_out = NULL;
    Expr* rem = NULL;
    Expr* pairs = intrat_log_part_core(expr_copy(a), expr_copy(d), tau, z,
                                       /*root_sum=*/false, Dd, xreal, &rem, decide);
    if (!pairs) { if (rem) expr_free(rem); return NULL; }
    /* Decision mode short-circuit: a non-constant residue yields the inert marker
     * Integrate`$NonConstantResidue (no pair list) — propagate it up verbatim for
     * the recursive-Risch caller to decode into a NONELEMENTARY verdict. */
    if (pairs->type == EXPR_SYMBOL
        && strcmp(pairs->data.symbol.name, "Integrate`$NonConstantResidue") == 0) {
        if (rem) expr_free(rem);
        return pairs;
    }
    Expr* real = intrat_log_to_real_pairs(pairs, tau, z);
    expr_free(pairs);
    if (!real) { if (rem) expr_free(rem); return NULL; }  /* LogToReal declined */
    if (remainder_out) *remainder_out = rem;
    else if (rem) expr_free(rem);   /* caller refuses partial: shouldn't happen */
    return real;
}

/* ====================================================================
 * Phase 3 — LogToAtan (Rioboo's recursive conversion).
 * ====================================================================
 *
 * Given A, B in K[x] with B != 0, return a sum of ArcTan of polynomials
 * in K[x] such that
 *    D[f, x] = D[I Log[(A + I B)/(A - I B)], x].
 *
 * Direct port of IntegrateRational.m:1529-1561.  Used by Phase 4's
 * LogToReal: each pair (A, B) of complex Re/Im polynomial parts feeds
 * LogToAtan to produce the corresponding real ArcTan summand.
 */

static Expr* intrat_log_to_atan(Expr* a, Expr* b, Expr* x) {
    intrat_trace("LogToAtan", "IN", a);

    /* If b is non-zero and free of x, return 0 (when a free of x) or
     * 2 ArcTan[a/b] otherwise. */
    bool b_free = intrat_freeq_test(b, x);
    if (!is_zero_poly(b) && b_free) {
        if (intrat_freeq_test(a, x)) return expr_new_integer(0);
        Expr* ratio = internal_divide((Expr*[]){expr_copy(a), expr_copy(b)}, 2);
        Expr* canon = intrat_canonic(ratio); expr_free(ratio);
        Expr* arctan = expr_new_function(expr_new_symbol(SYM_ArcTan),
            (Expr*[]){canon}, 1);
        Expr* res = internal_times((Expr*[]){expr_new_integer(2), arctan}, 2);
        return eval_and_free(res);
    }

    /* A := Collect[a, x] (we use intrat_canonic for the closest
     * Mathilda equivalent — it folds rational coefficients to canonical
     * form without RootReduce).  The .m baseline applies simproot to
     * each coefficient, but Phase 3 deliberately stays minimal so the
     * tests can still match against the un-simplified output. */
    Expr* A = expr_expand(a);
    Expr* B = expr_expand(b);

    /* If A is divisible by B in K[x], LogToAtan reduces to a single
     * 2 ArcTan[A/B] term (the recursion would otherwise terminate in
     * one step). */
    if (!is_zero_poly(B)) {
        Expr* rem = intrat_polyrem_t(A, B, x);
        if (rem && is_zero_poly(rem)) {
            expr_free(rem);
            Expr* q = intrat_polyq(A, B, x);
            if (intrat_freeq_test(q, x)) {
                expr_free(q); expr_free(A); expr_free(B);
                return expr_new_integer(0);
            }
            Expr* arctan = expr_new_function(expr_new_symbol(SYM_ArcTan),
                (Expr*[]){q}, 1);
            Expr* res = internal_times((Expr*[]){expr_new_integer(2), arctan}, 2);
            expr_free(A); expr_free(B);
            return eval_and_free(res);
        }
        if (rem) expr_free(rem);
    }

    /* If deg(A) < deg(B), swap (recurse with -B, A). */
    int dA = get_degree_poly(A, x);
    int dB = get_degree_poly(B, x);
    if (dA < dB) {
        Expr* negB = expand_and_free(internal_times(
            (Expr*[]){expr_new_integer(-1), expr_copy(B)}, 2));
        Expr* recurse = intrat_log_to_atan(negB, A, x);
        expr_free(negB); expr_free(A); expr_free(B);
        return recurse;
    }

    /* {g, {d, c}} = PolynomialExtendedGCD[B, -A, x]   so B*d - A*c = g. */
    Expr* negA = expand_and_free(internal_times(
        (Expr*[]){expr_new_integer(-1), expr_copy(A)}, 2));
    Expr* eg_args[3] = { expr_copy(B), negA, expr_copy(x) };
    Expr* eg = internal_polynomialextendedgcd(eg_args, 3);
    Expr* eg_eval = evaluate(eg);
    expr_free(eg);
    if (!eg_eval || eg_eval->type != EXPR_FUNCTION
        || eg_eval->data.function.arg_count != 2) {
        if (eg_eval) expr_free(eg_eval);
        expr_free(A); expr_free(B);
        return NULL;
    }
    Expr* g_poly = expr_copy(eg_eval->data.function.args[0]);
    Expr* coefs = eg_eval->data.function.args[1];
    if (coefs->type != EXPR_FUNCTION || coefs->data.function.arg_count != 2) {
        expr_free(eg_eval); expr_free(g_poly);
        expr_free(A); expr_free(B); return NULL;
    }
    Expr* d_coef = expr_copy(coefs->data.function.args[0]);
    Expr* c_coef = expr_copy(coefs->data.function.args[1]);
    expr_free(eg_eval);

    /* v = collectnum[(A*d + B*c)/g, x] // Together. */
    Expr* Ad_raw = internal_times(
        (Expr*[]){expr_copy(A), expr_copy(d_coef)}, 2);
    Expr* Ad = expr_expand(Ad_raw); expr_free(Ad_raw);
    Expr* Bc_raw = internal_times(
        (Expr*[]){expr_copy(B), expr_copy(c_coef)}, 2);
    Expr* Bc = expr_expand(Bc_raw); expr_free(Bc_raw);
    Expr* sum_raw = internal_plus(
        (Expr*[]){Ad, Bc}, 2);
    Expr* div = internal_divide(
        (Expr*[]){sum_raw, g_poly}, 2);
    Expr* v = intrat_canonic(div); expr_free(div);

    expr_free(A); expr_free(B);

    /* If v is free of x, return LogToAtan[d, c, x] (no ArcTan
     * contribution at this level). */
    if (intrat_freeq_test(v, x)) {
        expr_free(v);
        Expr* recurse = intrat_log_to_atan(d_coef, c_coef, x);
        expr_free(d_coef); expr_free(c_coef);
        return recurse;
    }

    /* Else: 2 ArcTan[v] + LogToAtan[d, c, x]. */
    Expr* arctan = expr_new_function(expr_new_symbol(SYM_ArcTan),
        (Expr*[]){v}, 1);
    Expr* head_term = internal_times(
        (Expr*[]){expr_new_integer(2), arctan}, 2);
    Expr* recurse = intrat_log_to_atan(d_coef, c_coef, x);
    expr_free(d_coef); expr_free(c_coef);
    if (!recurse) recurse = expr_new_integer(0);
    Expr* sum = internal_plus(
        (Expr*[]){head_term, recurse}, 2);
    Expr* result = eval_and_free(sum);
    intrat_trace("LogToAtan", "OUT", result);
    return result;
}

/* ----- Shared helpers used by both the Phase 2 linear-Q closer and
 *       the Phase 4 LogToReal dispatcher. ----- */

/* Extract the (t - α) form of a linear polynomial in `t` and return
 * the rational root α.  Returns NULL when the input is not linear in
 * t.  Fresh Expr*. */
static Expr* extract_linear_root(Expr* poly, Expr* t) {
    int deg = get_degree_poly(poly, t);
    if (deg != 1) return NULL;
    Expr* a_coef = get_coeff(poly, t, 1);
    Expr* b_coef = get_coeff(poly, t, 0);
    Expr* neg_b = internal_times(
        (Expr*[]){expr_new_integer(-1), b_coef}, 2);
    Expr* alpha = internal_divide((Expr*[]){neg_b, a_coef}, 2);
    return eval_and_free(alpha);
}

/* Produce a single Log term: α * Log[S(α, x)]. */
static Expr* build_log_term(Expr* alpha, Expr* S, Expr* t) {
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){expr_copy(t), expr_copy(alpha)}, 2);
    Expr* substituted = expr_new_function(expr_new_symbol(SYM_ReplaceAll),
        (Expr*[]){expr_copy(S), rule}, 2);
    Expr* Sat = evaluate(substituted);
    expr_free(substituted);
    Expr* logS = expr_new_function(expr_new_symbol(SYM_Log),
        (Expr*[]){Sat}, 1);
    Expr* term = internal_times(
        (Expr*[]){expr_copy(alpha), logS}, 2);
    return eval_and_free(term);
}

/* ====================================================================
 * Phase 4 — LogToReal (Rioboo's complex-to-real conversion).
 * ====================================================================
 *
 * Given r in K[t] and s in K[t,x], produce a real function f satisfying
 *   D[f, x] == D[RootSum[r(α)==0, α Log[s(α, x)]], x].
 *
 * Phase 4 implements the bounded-Solve subset called out in the plan:
 * each factor of r over Q is dispatched as
 *   linear (t - α): contribute α Log[s(α, x)];
 *   quadratic a t^2 + b t + c with disc >= 0:
 *       two real roots (-b ± sqrt(disc))/(2a) treated as two linear cases;
 *   quadratic a t^2 + b t + c with disc < 0:
 *       complex conjugate pair u ± I v with u = -b/(2a),
 *       v = sqrt(-disc)/(2a); contribute
 *           u Log[A^2 + B^2] + v LogToAtan[A, B, x]
 *       where A + I B = s(u + I v, x) split into real / imaginary parts.
 *
 * Higher-degree factors (or symbolic discriminants we can't decide)
 * fall through to the symbolic RootSum form so the result still
 * differentiates correctly via the residue formula.
 */

/* Walk an expanded `term` and split it into (real_part, imag_part).
 * Recognises Complex[a, b] either as a bare factor or as a factor
 * inside Times.  Anything else is treated as fully real. */
static void split_term_re_im(Expr* term, Expr** re_out, Expr** im_out) {
    if (term->type == EXPR_FUNCTION
        && term->data.function.head->type == EXPR_SYMBOL
        && term->data.function.head->data.symbol.name == SYM_Complex
        && term->data.function.arg_count == 2) {
        *re_out = expr_copy(term->data.function.args[0]);
        *im_out = expr_copy(term->data.function.args[1]);
        return;
    }
    if (term->type == EXPR_FUNCTION
        && term->data.function.head->type == EXPR_SYMBOL
        && term->data.function.head->data.symbol.name == SYM_Times) {
        /* Look for a Complex factor among the Times args. */
        Expr* complex_factor = NULL;
        size_t complex_idx = 0;
        size_t n = term->data.function.arg_count;
        for (size_t i = 0; i < n; i++) {
            Expr* arg = term->data.function.args[i];
            if (arg->type == EXPR_FUNCTION
                && arg->data.function.head->type == EXPR_SYMBOL
                && arg->data.function.head->data.symbol.name == SYM_Complex) {
                complex_factor = arg; complex_idx = i; break;
            }
        }
        if (!complex_factor) {
            *re_out = expr_copy(term);
            *im_out = expr_new_integer(0);
            return;
        }
        /* Build "rest" by dropping the complex factor. */
        size_t k = (n == 0) ? 0 : n - 1;
        Expr* rest;
        if (k == 0) rest = expr_new_integer(1);
        else if (k == 1) {
            rest = expr_copy(term->data.function.args[complex_idx == 0 ? 1 : 0]);
        } else {
            Expr** rest_args = (Expr**)malloc(sizeof(Expr*) * k);
            size_t kk = 0;
            for (size_t i = 0; i < n; i++) {
                if (i == complex_idx) continue;
                rest_args[kk++] = expr_copy(term->data.function.args[i]);
            }
            rest = eval_and_free(internal_times(rest_args, kk));
            free(rest_args);
        }
        Expr* a_part = expr_copy(complex_factor->data.function.args[0]);
        Expr* b_part = expr_copy(complex_factor->data.function.args[1]);
        *re_out = eval_and_free(internal_times(
            (Expr*[]){a_part, expr_copy(rest)}, 2));
        *im_out = eval_and_free(internal_times(
            (Expr*[]){b_part, rest}, 2));
        return;
    }
    /* No Complex anywhere — purely real. */
    *re_out = expr_copy(term);
    *im_out = expr_new_integer(0);
}

/* Expand p, then walk its Plus head, splitting each summand. */
static void split_re_im(Expr* p, Expr** re_out, Expr** im_out) {
    Expr* expanded = expr_expand(p);
    if (!(expanded->type == EXPR_FUNCTION
        && expanded->data.function.head->type == EXPR_SYMBOL
        && expanded->data.function.head->data.symbol.name == SYM_Plus)) {
        split_term_re_im(expanded, re_out, im_out);
        expr_free(expanded);
        return;
    }
    size_t n = expanded->data.function.arg_count;
    Expr** re_terms = (Expr**)malloc(sizeof(Expr*) * n);
    Expr** im_terms = (Expr**)malloc(sizeof(Expr*) * n);
    size_t nre = 0, nim = 0;
    for (size_t i = 0; i < n; i++) {
        Expr *re_t = NULL, *im_t = NULL;
        split_term_re_im(expanded->data.function.args[i], &re_t, &im_t);
        if (re_t && !is_zero_poly(re_t)) re_terms[nre++] = re_t;
        else if (re_t) expr_free(re_t);
        if (im_t && !is_zero_poly(im_t)) im_terms[nim++] = im_t;
        else if (im_t) expr_free(im_t);
    }
    *re_out = (nre == 0) ? expr_new_integer(0)
            : (nre == 1) ? re_terms[0]
            : eval_and_free(internal_plus(re_terms, nre));
    *im_out = (nim == 0) ? expr_new_integer(0)
            : (nim == 1) ? im_terms[0]
            : eval_and_free(internal_plus(im_terms, nim));
    free(re_terms); free(im_terms);
    expr_free(expanded);
}

/* Substitute t -> sub_expr inside e (full ReplaceAll). */
static Expr* subst_t(Expr* e, Expr* t, Expr* sub_expr) {
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){expr_copy(t), expr_copy(sub_expr)}, 2);
    Expr* call = expr_new_function(expr_new_symbol(SYM_ReplaceAll),
        (Expr*[]){expr_copy(e), rule}, 2);
    Expr* result = evaluate(call);
    expr_free(call);
    return result;
}

/* Sign tests, canonic-zero, and the positive-symbol Sqrt simplifier
 * used by logtoreal_dispatch live in intsimp.c — see intsimp.h. */

/* Try to dispatch a single quadratic factor `Q(t) = a t^2 + b t + c`
 * to its real-form contribution.  Returns NULL on failure (caller
 * falls back to symbolic RootSum). */
static Expr* logtoreal_quadratic(Expr* a, Expr* b, Expr* c,
                                 Expr* s, Expr* x, Expr* t) {
    /* discriminant = b^2 - 4 a c. */
    Expr* b2 = expand_and_free(internal_power((Expr*[]){expr_copy(b), expr_new_integer(2)}, 2));
    Expr* fourac = expand_and_free(
        internal_times((Expr*[]){expr_new_integer(4), expr_copy(a), expr_copy(c)}, 3));
    Expr* neg_fourac = internal_times(
        (Expr*[]){expr_new_integer(-1), fourac}, 2);
    Expr* disc = eval_and_free(internal_plus((Expr*[]){b2, neg_fourac}, 2));

    /* Decide the sign of disc.  Rational cases give a definitive
     * sign; symbolic cases are decided under the positive-symbol
     * assumption (matching IntegrateRational.m's `Refine[#>0]&` filter
     * that gates the parametric path).  When the sign cannot be
     * determined (parametric quadratic with mixed-sign discriminant
     * such as 1/(a x^2 + b x + c)), default to the negative-disc /
     * ArcTan branch with a held Sqrt[-disc].  That formula is the
     * principal antiderivative for both signs of the discriminant
     * (ArcTan continues analytically through imaginary argument),
     * matching Mathematica's parametric-quadratic output. */
    int disc_sign = intsimp_sign_pos_assumption(disc);

    Expr* two_a = eval_and_free(internal_times(
        (Expr*[]){expr_new_integer(2), expr_copy(a)}, 2));

    if (disc_sign > 0) {
        /* Two real roots.  α± = (-b ± Sqrt[disc])/(2 a). */
        Expr* sqrt_d = expr_new_function(expr_new_symbol(SYM_Sqrt),
            (Expr*[]){expr_copy(disc)}, 1);
        sqrt_d = eval_and_free(sqrt_d);
        Expr* neg_b = internal_times(
            (Expr*[]){expr_new_integer(-1), expr_copy(b)}, 2);
        neg_b = eval_and_free(neg_b);
        Expr* num_plus = eval_and_free(internal_plus(
            (Expr*[]){expr_copy(neg_b), expr_copy(sqrt_d)}, 2));
        Expr* num_minus = eval_and_free(internal_subtract(
            (Expr*[]){neg_b, sqrt_d}, 2));
        Expr* alpha_plus = eval_and_free(internal_divide(
            (Expr*[]){num_plus, expr_copy(two_a)}, 2));
        Expr* alpha_minus = eval_and_free(internal_divide(
            (Expr*[]){num_minus, two_a}, 2));

        /* Each contributes α * Log[s(α, x)]. */
        Expr* s_at_plus  = subst_t(s, t, alpha_plus);
        Expr* s_at_minus = subst_t(s, t, alpha_minus);
        Expr* log_plus = expr_new_function(expr_new_symbol(SYM_Log),
            (Expr*[]){s_at_plus}, 1);
        Expr* log_minus = expr_new_function(expr_new_symbol(SYM_Log),
            (Expr*[]){s_at_minus}, 1);
        Expr* term_plus  = internal_times((Expr*[]){alpha_plus, log_plus}, 2);
        Expr* term_minus = internal_times((Expr*[]){alpha_minus, log_minus}, 2);
        expr_free(disc);
        Expr* sum = internal_plus(
            (Expr*[]){term_plus, term_minus}, 2);
        return eval_and_free(sum);
    }

    /* Complex conjugate pair: u = -b/(2a), v = Sqrt[-disc]/(2a). */
    Expr* neg_b2 = eval_and_free(internal_times(
        (Expr*[]){expr_new_integer(-1), expr_copy(b)}, 2));
    Expr* u_root = eval_and_free(internal_divide(
        (Expr*[]){neg_b2, expr_copy(two_a)}, 2));
    /* Distribute -1 across `disc` before passing to the Sqrt
     * simplifier.  Without an Expand, Mathilda keeps `-1 * disc` as
     * Times[-1, Plus[...]]; intsimp_pos_sqrt then walks the Times
     * factors and pulls Sqrt[-1] = I out of the `-1`, ending in an
     * imaginary v_root that collapses logtoreal_quadratic to 0. */
    Expr* neg_disc_raw = expand_and_free(internal_times(
        (Expr*[]){expr_new_integer(-1), disc}, 2));
    /* Try the positive-symbol Sqrt simplifier first (e.g.
     * Sqrt[4 a^2 b^2] -> 2 a b); fall back to the held Sqrt[..] when
     * it can't extract a clean radical form. */
    Expr* sqrt_neg_disc = intsimp_pos_sqrt(neg_disc_raw);
    expr_free(neg_disc_raw);
    Expr* v_root = eval_and_free(internal_divide(
        (Expr*[]){sqrt_neg_disc, two_a}, 2));

    /* Substitute t -> u_root + I*v_root in s and split into A + I B. */
    Expr* iv = internal_times((Expr*[]){
        expr_new_function(expr_new_symbol(SYM_Complex),
            (Expr*[]){expr_new_integer(0), expr_new_integer(1)}, 2),
        expr_copy(v_root)
    }, 2);
    iv = eval_and_free(iv);
    Expr* substituted_t = eval_and_free(internal_plus(
        (Expr*[]){expr_copy(u_root), iv}, 2));
    Expr* s_complex = subst_t(s, t, substituted_t);
    expr_free(substituted_t);

    Expr *A = NULL, *B = NULL;
    split_re_im(s_complex, &A, &B);
    expr_free(s_complex);

    /* Term 1: u_root * Log[A^2 + B^2]. */
    Expr* A2 = expand_and_free(internal_power(
        (Expr*[]){expr_copy(A), expr_new_integer(2)}, 2));
    Expr* B2 = expand_and_free(internal_power(
        (Expr*[]){expr_copy(B), expr_new_integer(2)}, 2));
    Expr* mod2 = eval_and_free(internal_plus((Expr*[]){A2, B2}, 2));
    Expr* mod2_can = intrat_canonic(mod2); expr_free(mod2);
    Expr* logmod = expr_new_function(expr_new_symbol(SYM_Log),
        (Expr*[]){mod2_can}, 1);
    Expr* term1 = internal_times((Expr*[]){expr_copy(u_root), logmod}, 2);
    term1 = eval_and_free(term1);

    /* Term 2: v_root * LogToAtan[A, B, x]. */
    Expr* atan_part = intrat_log_to_atan(A, B, x);
    if (!atan_part) atan_part = expr_new_integer(0);
    Expr* term2 = internal_times((Expr*[]){expr_copy(v_root), atan_part}, 2);
    term2 = eval_and_free(term2);

    expr_free(A); expr_free(B); expr_free(u_root); expr_free(v_root);

    Expr* sum = internal_plus((Expr*[]){term1, term2}, 2);
    return eval_and_free(sum);
}

/* Handle a sparse nth-root depressed polynomial c_n t^n + c_0 (all
 * intermediate coefficients zero) for n >= 3.  Factor over R via the
 * standard cyclotomic decomposition:
 *
 *   roots are r e^{i θ_k} with r = |c_0 / c_n|^(1/n) and angles
 *     q := -c_0/c_n > 0:  θ_k = 2 k π / n   (k = 0..n-1)
 *     q < 0:              θ_k = (2 k + 1) π / n
 *
 * Real roots correspond to angles 0 and π.  Other roots come in
 * complex-conjugate pairs (k, n−k) for q > 0 or (k, n−1−k) for q < 0.
 * Each pair contributes a quadratic factor
 *     t^2 − 2 r cos(θ_k) t + r^2
 * whose discriminant −4 r^2 sin^2 θ_k is strictly negative, so
 * logtoreal_quadratic always lands in its ArcTan branch.
 *
 * Generalises the deg-3 and deg-4-Sophie-Germain branches that used
 * to live in logtoreal_dispatch.  Catches 1/(b ± a x^n) for all
 * n ≥ 3 once the LRT pipeline reduces to its sparse Q-in-t form.
 *
 * Returns the accumulated real form, or NULL when:
 *   • the polynomial isn't sparse (some intermediate c_k ≠ 0),
 *   • c_0 = 0 (caller should have factored t out upstream), or
 *   • the sign of −c_0/c_n cannot be decided symbolically or numerically.
 */
static Expr* logtoreal_nthroot_sparse(Expr* base, int deg,
                                      Expr* s, Expr* x, Expr* t) {
    if (deg < 3) return NULL;

    for (int k = 1; k < deg; k++) {
        Expr* ck = get_coeff(base, t, k);
        bool zero = (ck->type == EXPR_INTEGER && ck->data.integer == 0);
        expr_free(ck);
        if (!zero) return NULL;
    }

    Expr* c0 = get_coeff(base, t, 0);
    Expr* cn = get_coeff(base, t, deg);
    if (c0->type == EXPR_INTEGER && c0->data.integer == 0) {
        expr_free(c0); expr_free(cn);
        return NULL;
    }

    /* q := −c_0 / c_n. */
    Expr* neg_c0 = expand_and_free(internal_times(
        (Expr*[]){expr_new_integer(-1), expr_copy(c0)}, 2));
    Expr* q = eval_and_free(internal_divide(
        (Expr*[]){neg_c0, expr_copy(cn)}, 2));
    expr_free(c0); expr_free(cn);

    int q_sign = intsimp_sign_pos_assumption(q);
    if (q_sign == 0) q_sign = intsimp_numeric_sign(q);
    if (q_sign == 0) { expr_free(q); return NULL; }

    bool q_pos = (q_sign > 0);

    /* r = |q|^(1/n) using the principal real n-th root. */
    Expr* abs_q = q_pos
        ? q
        : eval_and_free(internal_times(
              (Expr*[]){expr_new_integer(-1), q}, 2));
    Expr* inv_n = expr_new_function(
        expr_new_symbol(SYM_Rational),
        (Expr*[]){expr_new_integer(1), expr_new_integer(deg)}, 2);
    inv_n = eval_and_free(inv_n);
    Expr* r = eval_and_free(internal_power(
        (Expr*[]){abs_q, inv_n}, 2));

    size_t cap = (size_t)deg / 2 + 4;
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * cap);
    size_t nterms = 0;
    bool fail = false;

    /* --- Linear factors for real roots. -------------------------- */
    if (q_pos) {
        /* (t − r) — always present when q > 0. */
        Expr* lin = build_log_term(r, s, t);
        if (!lin) { fail = true; goto done; }
        terms[nterms++] = lin;
        if (deg % 2 == 0) {
            /* (t + r) when n is even. */
            Expr* neg_r = eval_and_free(internal_times(
                (Expr*[]){expr_new_integer(-1), expr_copy(r)}, 2));
            Expr* lin2 = build_log_term(neg_r, s, t);
            expr_free(neg_r);
            if (!lin2) { fail = true; goto done; }
            terms[nterms++] = lin2;
        }
    } else {
        /* q < 0: real root only when n is odd, at t = −r. */
        if (deg % 2 == 1) {
            Expr* neg_r = eval_and_free(internal_times(
                (Expr*[]){expr_new_integer(-1), expr_copy(r)}, 2));
            Expr* lin = build_log_term(neg_r, s, t);
            expr_free(neg_r);
            if (!lin) { fail = true; goto done; }
            terms[nterms++] = lin;
        }
    }

    /* --- Quadratic factors for complex-conjugate pairs. ---------- */
    int k_start, k_end;
    if (q_pos) {
        k_start = 1;
        k_end = (deg - 1) / 2;          /* floor((n-1)/2) */
    } else {
        k_start = 0;
        k_end = (deg % 2 == 1)
            ? (deg - 3) / 2              /* odd n: skip k=(n-1)/2 (θ=π real root) */
            : (deg - 2) / 2;             /* even n: floor((n-2)/2) */
    }

    Expr* r_sq = expand_and_free(internal_power(
        (Expr*[]){expr_copy(r), expr_new_integer(2)}, 2));

    for (int k = k_start; k <= k_end; k++) {
        int angle_num_k = q_pos ? (2 * k) : (2 * k + 1);
        Expr* angle_num = eval_and_free(internal_times(
            (Expr*[]){expr_new_integer(angle_num_k),
                      expr_new_symbol(SYM_Pi)}, 2));
        Expr* angle = eval_and_free(internal_divide(
            (Expr*[]){angle_num, expr_new_integer(deg)}, 2));
        Expr* cos_v = expr_new_function(expr_new_symbol(SYM_Cos),
            (Expr*[]){angle}, 1);
        cos_v = eval_and_free(cos_v);

        Expr* b_quad = eval_and_free(internal_times(
            (Expr*[]){expr_new_integer(-2),
                      expr_copy(r),
                      cos_v}, 3));

        Expr* one_a = expr_new_integer(1);
        Expr* r_sq_copy = expr_copy(r_sq);
        Expr* quad = logtoreal_quadratic(
            one_a,
            b_quad,
            r_sq_copy,
            s, x, t);
        expr_free(one_a);
        expr_free(b_quad);
        expr_free(r_sq_copy);
        if (!quad) { fail = true; break; }
        terms[nterms++] = quad;
    }

    expr_free(r_sq);
    expr_free(r);

done:
    if (fail || nterms == 0) {
        for (size_t k = 0; k < nterms; k++) expr_free(terms[k]);
        free(terms);
        return NULL;
    }
    if (nterms == 1) {
        Expr* result = terms[0];
        free(terms);
        return result;
    }
    Expr* sum = internal_plus(terms, nterms);
    free(terms);
    return eval_and_free(sum);
}

/* Walk a Factor[Q]'s output and dispatch every factor.  Returns the
 * accumulated real form, or NULL if any factor exceeds the bounded-
 * Solve scope (degree > 2 in t). */
static Expr* logtoreal_dispatch(Expr* factored, Expr* s, Expr* x, Expr* t) {
    Expr** factors;
    size_t nfactors;
    if (factored->type == EXPR_FUNCTION
        && factored->data.function.head->type == EXPR_SYMBOL
        && factored->data.function.head->data.symbol.name == SYM_Times) {
        factors = factored->data.function.args;
        nfactors = factored->data.function.arg_count;
    } else {
        factors = &factored;
        nfactors = 1;
    }
    Expr** acc = (Expr**)malloc(sizeof(Expr*) * (nfactors * 2));
    size_t acc_n = 0;

    for (size_t i = 0; i < nfactors; i++) {
        Expr* fac = factors[i];
        Expr* base = fac;
        int multiplicity = 1;
        if (fac->type == EXPR_FUNCTION
            && fac->data.function.head->type == EXPR_SYMBOL
            && fac->data.function.head->data.symbol.name == SYM_Power
            && fac->data.function.arg_count == 2
            && fac->data.function.args[1]->type == EXPR_INTEGER
            && fac->data.function.args[1]->data.integer >= 1) {
            base = fac->data.function.args[0];
            multiplicity = (int)fac->data.function.args[1]->data.integer;
        }
        int deg = get_degree_poly(base, t);
        if (deg == 0) continue;  /* numeric coefficient */

        Expr* contribution = NULL;
        if (deg == 1) {
            Expr* alpha = extract_linear_root(base, t);
            if (!alpha) goto fail;
            contribution = build_log_term(alpha, s, t);
            expr_free(alpha);
        } else if (deg == 2) {
            Expr* a = get_coeff(base, t, 2);
            Expr* b = get_coeff(base, t, 1);
            Expr* c = get_coeff(base, t, 0);
            contribution = logtoreal_quadratic(a, b, c, s, x, t);
            expr_free(a); expr_free(b); expr_free(c);
            if (!contribution) goto fail;
        } else if (deg >= 3) {
            /* General sparse nth-root path: c_n t^n + c_0 (all
             * intermediate coefficients zero), dispatched via the
             * cyclotomic factorisation of t^n − q.  Subsumes the
             * previous deg-3 nth-root branch and the deg-4
             * Sophie-Germain shortcut; catches 1/(b ± a x^n) for
             * every n ≥ 3.  When sparse fails for deg == 4 we fall
             * through to the deg-4 biquadratic / palindromic
             * specialisations below. */
            contribution = logtoreal_nthroot_sparse(base, deg, s, x, t);
            if (!contribution && deg != 4) goto fail;
        }
        if (!contribution && deg == 4) {
            /* Non-sparse quartic: try biquadratic (c1 = c3 = 0,
             * c2 ≠ 0) or palindromic specialisations.  The fully
             * sparse case (c1 = c2 = c3 = 0) was already handled by
             * logtoreal_nthroot_sparse above; if it didn't fire, the
             * sign of c0/c4 was undecidable and there's nothing more
             * to do at this layer. */
            Expr* c0 = get_coeff(base, t, 0);
            Expr* c1 = get_coeff(base, t, 1);
            Expr* c2 = get_coeff(base, t, 2);
            Expr* c3 = get_coeff(base, t, 3);
            Expr* c4 = get_coeff(base, t, 4);

            bool c1z = (c1->type == EXPR_INTEGER && c1->data.integer == 0);
            bool c2z = (c2->type == EXPR_INTEGER && c2->data.integer == 0);
            bool c3z = (c3->type == EXPR_INTEGER && c3->data.integer == 0);

            if (c1z && c2z && c3z) {
                /* Fully sparse but sparse handler bailed (sign of
                 * c0/c4 undecidable both symbolically and numerically).
                 * No further branch can do better. */
                expr_free(c0); expr_free(c1); expr_free(c2);
                expr_free(c3); expr_free(c4);
                goto fail;
            } else if (c1z && c3z) {
                /* Genuine biquadratic in t: c4 t^4 + c2 t^2 + c0
                 * (c2 != 0).  Substitute u = t^2 to land in
                 *   c4 u^2 + c2 u + c0 = 0
                 * with roots u_± = (-c2 ± Sqrt[c2^2 - 4 c0 c4]) / (2 c4).
                 * Q(t) factors as c4 (t^2 - u_+)(t^2 - u_-).
                 *
                 * Soundness gate: only fire when c2^2 - 4 c0 c4 is
                 * provably positive under the positive-symbol
                 * assumption.  Otherwise u_± would be a complex
                 * conjugate pair, the (t^2 - u_±) factors would have
                 * complex coefficients, and logtoreal_quadratic
                 * (which expects a real-coefficient quadratic in t)
                 * would silently miscompute -- typically collapsing
                 * to 0 because the imaginary parts cancel between
                 * the two sub-factor calls.  Catches the parametric
                 * family 1/(a x^4 + b x^2 + c) under the
                 * b^2 > 4 a c assumption (which is always true in
                 * the positive-symbol model) while leaving cases
                 * with negative or sign-unknown inner discriminant
                 * to fall through to NaiveLogPart's held RootSum
                 * (where they were before this branch existed). */
                Expr* c2sq = expand_and_free(internal_power(
                    (Expr*[]){expr_copy(c2), expr_new_integer(2)}, 2));
                Expr* fourc0c4 = expand_and_free(internal_times(
                    (Expr*[]){expr_new_integer(4), expr_copy(c0), expr_copy(c4)}, 3));
                Expr* neg_fourc0c4 = expand_and_free(internal_times(
                    (Expr*[]){expr_new_integer(-1), fourc0c4}, 2));
                Expr* discsq = eval_and_free(internal_plus(
                    (Expr*[]){c2sq, neg_fourc0c4}, 2));
                int inner_disc_sign = intsimp_sign_pos_assumption(discsq);
                if (inner_disc_sign <= 0) {
                    /* Inner disc not provably positive.  Before giving
                     * up, try the Sophie-Germain-with-c2 factorisation
                     * over R:
                     *   c4 t^4 + c2 t^2 + c0
                     *      = c4 (t^2 + α t + β)(t^2 − α t + β)
                     * with β = Sqrt[c0/c4], α^2 = 2 β − c2/c4.
                     * Real when (a) c0/c4 > 0 (so β is real) and
                     * (b) 2 β − c2/c4 > 0 (so α is real).  Both gates
                     * fire under the positive-symbol assumption, and
                     * algebraically: c2^2 − 4 c0 c4 < 0 ⇒
                     * |c2| < 2 Sqrt[c0 c4] ⇒ 2 Sqrt[c0/c4] − c2/c4 =
                     * (2 Sqrt[c0 c4] − c2) / c4 > 0 when c4 > 0,
                     * so the gate is automatic in the negative-inner-
                     * disc regime we land in here.
                     *
                     * Each resulting real quadratic has its own
                     * discriminant α^2 − 4 β < 0, which routes through
                     * logtoreal_quadratic's ArcTan branch.  This is
                     * the path Mathematica's BronsteinRational takes
                     * for inputs like (-1+x^2)/(1-2x^2+2x^4) whose
                     * Q-in-t is 2 t^4 − 2 t^2 + 1 (inner disc −4). */
                    expr_free(discsq);

                    Expr* ratio = eval_and_free(internal_divide(
                        (Expr*[]){expr_copy(c0), expr_copy(c4)}, 2));
                    bool ratio_pos = (intsimp_sign_pos_assumption(ratio) > 0);

                    if (!ratio_pos) {
                        expr_free(ratio);
                        expr_free(c0); expr_free(c1); expr_free(c2);
                        expr_free(c3); expr_free(c4);
                        goto fail;
                    }

                    Expr* beta = intsimp_pos_sqrt(ratio); /* Sqrt[c0/c4] */
                    expr_free(ratio);

                    Expr* c2_over_c4 = eval_and_free(internal_divide(
                        (Expr*[]){expr_copy(c2), expr_copy(c4)}, 2));
                    Expr* neg_c2_over_c4 = expand_and_free(internal_times(
                        (Expr*[]){expr_new_integer(-1), c2_over_c4}, 2));

                    /* α^2 = 2 β − c2/c4. */
                    Expr* two_beta = eval_and_free(internal_times(
                        (Expr*[]){expr_new_integer(2), expr_copy(beta)}, 2));
                    Expr* alpha_sq = eval_and_free(internal_plus(
                        (Expr*[]){two_beta, neg_c2_over_c4}, 2));
                    Expr* alpha_sq_expanded = expr_expand(alpha_sq);
                    expr_free(alpha_sq);

                    if (intsimp_sign_pos_assumption(alpha_sq_expanded) <= 0) {
                        /* α^2 not provably positive — bail. */
                        expr_free(alpha_sq_expanded);
                        expr_free(beta);
                        expr_free(c0); expr_free(c1); expr_free(c2);
                        expr_free(c3); expr_free(c4);
                        goto fail;
                    }

                    Expr* alpha = intsimp_pos_sqrt(alpha_sq_expanded);
                    expr_free(alpha_sq_expanded);
                    Expr* neg_alpha = eval_and_free(internal_times(
                        (Expr*[]){expr_new_integer(-1),
                                  expr_copy(alpha)}, 2));

                    Expr* one_p1 = expr_new_integer(1);
                    Expr* beta_copy = expr_copy(beta);
                    Expr* part1 = logtoreal_quadratic(
                        one_p1,
                        alpha,
                        beta_copy,
                        s, x, t);
                    expr_free(one_p1); expr_free(beta_copy);

                    Expr* one_p2 = expr_new_integer(1);
                    Expr* part2 = logtoreal_quadratic(
                        one_p2,
                        neg_alpha,
                        beta,
                        s, x, t);
                    expr_free(one_p2);
                    expr_free(alpha); expr_free(neg_alpha); expr_free(beta);

                    expr_free(c0); expr_free(c1); expr_free(c2);
                    expr_free(c3); expr_free(c4);

                    if (part1 && part2) {
                        contribution = eval_and_free(internal_plus(
                            (Expr*[]){part1, part2}, 2));
                        for (int m = 0; m < multiplicity; m++) {
                            acc[acc_n++] = (m == 0) ? contribution : expr_copy(contribution);
                        }
                        continue; /* next factor */
                    }
                    if (part1) expr_free(part1);
                    if (part2) expr_free(part2);
                    goto fail;
                }
                Expr* sqrt_disc = intsimp_pos_sqrt(discsq);
                expr_free(discsq);

                Expr* neg_c2 = expand_and_free(internal_times(
                    (Expr*[]){expr_new_integer(-1), expr_copy(c2)}, 2));
                Expr* two_c4 = eval_and_free(internal_times(
                    (Expr*[]){expr_new_integer(2), expr_copy(c4)}, 2));

                Expr* num_plus = eval_and_free(internal_plus(
                    (Expr*[]){expr_copy(neg_c2), expr_copy(sqrt_disc)}, 2));
                Expr* num_minus = eval_and_free(internal_subtract(
                    (Expr*[]){neg_c2, sqrt_disc}, 2));
                Expr* u_plus = eval_and_free(internal_divide(
                    (Expr*[]){num_plus, expr_copy(two_c4)}, 2));
                Expr* u_minus = eval_and_free(internal_divide(
                    (Expr*[]){num_minus, two_c4}, 2));

                Expr* neg_u_plus = expand_and_free(internal_times(
                    (Expr*[]){expr_new_integer(-1), u_plus}, 2));
                Expr* neg_u_minus = expand_and_free(internal_times(
                    (Expr*[]){expr_new_integer(-1), u_minus}, 2));

                Expr* one_a1 = expr_new_integer(1);
                Expr* zero_b1 = expr_new_integer(0);
                Expr* part1 = logtoreal_quadratic(
                    one_a1,
                    zero_b1,
                    neg_u_plus,
                    s, x, t);
                expr_free(one_a1); expr_free(zero_b1); expr_free(neg_u_plus);

                Expr* one_a2 = expr_new_integer(1);
                Expr* zero_b2 = expr_new_integer(0);
                Expr* part2 = logtoreal_quadratic(
                    one_a2,
                    zero_b2,
                    neg_u_minus,
                    s, x, t);
                expr_free(one_a2); expr_free(zero_b2); expr_free(neg_u_minus);

                expr_free(c0); expr_free(c1); expr_free(c2);
                expr_free(c3); expr_free(c4);

                if (part1 && part2) {
                    contribution = eval_and_free(internal_plus(
                        (Expr*[]){part1, part2}, 2));
                } else {
                    if (part1) expr_free(part1);
                    if (part2) expr_free(part2);
                    goto fail;
                }
            } else {
                /* Scaled-palindromic quartic:
                 *   Q(t) = c4 t^4 + c3 t^3 + c2 t^2 + c1 t + c0
                 * with c4 c1^2 = c3^2 c0 (the symmetry condition that
                 * generalises the c0=c4, c1=c3 palindromic family).
                 * Under u = r t with r^2 = c3/c1 = Sqrt[c4/c0] the
                 * polynomial r^n Q(u/r) is palindromic in u:
                 *   c4 u^4 + (c3 r) u^3 + (c2 r^2) u^2 + (c1 r^3) u
                 *   + (c0 r^4)   with c0 r^4 = c4 and c1 r^3 = c3 r.
                 * Substitute v = u + 1/u (so u^2 + 1/u^2 = v^2 − 2):
                 *   c4 v^2 + (c3 r) v + (c2 r^2 − 2 c4) = 0.
                 * Solve for v_±; each (u^2 − v_± u + 1) translates to
                 *   r^2 t^2 − r v_± t + 1
                 * which is dispatched through logtoreal_quadratic.
                 *
                 * Catches Q-in-t for 1/(x^5+1) (Q = 625 t^4 + 125 t^3
                 * + 25 t^2 + 5 t + 1 — palindromic in u = 5t) and the
                 * scaled cyclotomic-quartic family.  When the symmetry
                 * condition fails or any sign-gate is undecidable we
                 * bail to RootSum. */

                /* Symmetry check: c4 c1^2 == c3^2 c0  AND  c0 != 0.
                 * General scaled-palindromic detection.  We further
                 * restrict to r = 1 (pure palindromic) below, because
                 * the scaled cases produce nested-radical S(t/r, x)
                 * polynomials that overwhelm the downstream LogToAtan
                 * GCD machinery.  Scaled-palindromic-with-clean-r is a
                 * future extension. */
                Expr* sym_lhs = expand_and_free(internal_times(
                    (Expr*[]){expr_copy(c4),
                              expand_and_free(internal_power(
                                  (Expr*[]){expr_copy(c1),
                                            expr_new_integer(2)}, 2))}, 2));
                Expr* sym_rhs = expand_and_free(internal_times(
                    (Expr*[]){expand_and_free(internal_power(
                                  (Expr*[]){expr_copy(c3),
                                            expr_new_integer(2)}, 2)),
                              expr_copy(c0)}, 2));
                Expr* sym_diff = internal_subtract(
                    (Expr*[]){sym_lhs, sym_rhs}, 2);
                bool symmetric = intsimp_zero_q(sym_diff);
                expr_free(sym_diff);

                bool c0_nonzero = !(c0->type == EXPR_INTEGER && c0->data.integer == 0);
                bool c1_nonzero = !(c1->type == EXPR_INTEGER && c1->data.integer == 0);

                if (!symmetric || !c0_nonzero || !c1_nonzero) {
                    expr_free(c0); expr_free(c1); expr_free(c2);
                    expr_free(c3); expr_free(c4);
                    goto fail;
                }

                /* r^2 = c3 / c1.  Restrict to r = 1 (c1 == c3 AND
                 * c0 == c4) so the resulting LogToReal substitution
                 * stays radical-clean. */
                Expr* c1mc3 = internal_subtract(
                    (Expr*[]){expr_copy(c1), expr_copy(c3)}, 2);
                Expr* c0mc4 = internal_subtract(
                    (Expr*[]){expr_copy(c0), expr_copy(c4)}, 2);
                bool r_is_one = intsimp_zero_q(c1mc3) && intsimp_zero_q(c0mc4);
                expr_free(c1mc3); expr_free(c0mc4);
                if (!r_is_one) {
                    /* Scaled-palindromic with r != 1: bail (the
                     * substituted S(t/r, x) blows up downstream). */
                    expr_free(c0); expr_free(c1); expr_free(c2);
                    expr_free(c3); expr_free(c4);
                    goto fail;
                }
                Expr* r_sq = expr_new_integer(1);
                Expr* r    = expr_new_integer(1);

                /* v-quadratic: c4 v^2 + (c3 r) v + (c2 r^2 − 2 c4) = 0.
                 * disc_v = (c3 r)^2 − 4 c4 (c2 r^2 − 2 c4)
                 *        = c3^2 r^2 − 4 c2 c4 r^2 + 8 c4^2.
                 * Soundness gate: disc_v >= 0 so v_± are real. */
                Expr* c3r = eval_and_free(internal_times(
                    (Expr*[]){expr_copy(c3), expr_copy(r)}, 2));
                Expr* c3r_sq = expand_and_free(internal_power(
                    (Expr*[]){expr_copy(c3r), expr_new_integer(2)}, 2));
                Expr* two_c4 = eval_and_free(internal_times(
                    (Expr*[]){expr_new_integer(2), expr_copy(c4)}, 2));
                Expr* c2_rsq = expand_and_free(internal_times(
                    (Expr*[]){expr_copy(c2), expr_copy(r_sq)}, 2));
                Expr* v_const = eval_and_free(internal_subtract(
                    (Expr*[]){c2_rsq, two_c4}, 2));
                Expr* four_c4 = expand_and_free(internal_times(
                    (Expr*[]){expr_new_integer(4), expr_copy(c4),
                              expr_copy(v_const)}, 3));
                Expr* neg_four_c4 = expand_and_free(internal_times(
                    (Expr*[]){expr_new_integer(-1), four_c4}, 2));
                Expr* disc_v = eval_and_free(internal_plus(
                    (Expr*[]){c3r_sq, neg_four_c4}, 2));
                Expr* disc_v_exp = expr_expand(disc_v); expr_free(disc_v);
                int disc_v_sign = intsimp_sign_pos_assumption(disc_v_exp);
                if (disc_v_sign == 0) disc_v_sign = intsimp_numeric_sign(disc_v_exp);

                if (disc_v_sign < 0) {
                    /* v complex — bail. */
                    expr_free(v_const); expr_free(disc_v_exp);
                    expr_free(c3r); expr_free(r); expr_free(r_sq);
                    expr_free(c0); expr_free(c1); expr_free(c2);
                    expr_free(c3); expr_free(c4);
                    goto fail;
                }
                Expr* sqrt_disc_v = intsimp_pos_sqrt(disc_v_exp);
                expr_free(disc_v_exp);

                /* v_± = (−c3 r ± Sqrt[disc_v]) / (2 c4). */
                Expr* neg_c3r = eval_and_free(internal_times(
                    (Expr*[]){expr_new_integer(-1), expr_copy(c3r)}, 2));
                Expr* v_num_p = eval_and_free(internal_plus(
                    (Expr*[]){expr_copy(neg_c3r), expr_copy(sqrt_disc_v)}, 2));
                Expr* v_num_m = eval_and_free(internal_subtract(
                    (Expr*[]){neg_c3r, sqrt_disc_v}, 2));
                Expr* two_c4_b = eval_and_free(internal_times(
                    (Expr*[]){expr_new_integer(2), expr_copy(c4)}, 2));
                Expr* v_plus = eval_and_free(internal_divide(
                    (Expr*[]){v_num_p, expr_copy(two_c4_b)}, 2));
                Expr* v_minus = eval_and_free(internal_divide(
                    (Expr*[]){v_num_m, two_c4_b}, 2));

                expr_free(v_const); expr_free(c3r);
                expr_free(c0); expr_free(c1); expr_free(c2);
                expr_free(c3); expr_free(c4);

                /* Each v_± yields the quadratic in t:
                 *   r^2 t^2 − (r v_±) t + 1
                 * Discriminant Δ_t = r^2 v^2 − 4 r^2 = r^2 (v^2 − 4).
                 * Verify Δ_t < 0 so logtoreal_quadratic's ArcTan
                 * branch fires.  When the symbolic walker is
                 * inconclusive we fall back to numerical N[]. */
                Expr* v_plus_sq = expand_and_free(internal_power(
                    (Expr*[]){expr_copy(v_plus), expr_new_integer(2)}, 2));
                Expr* v_plus_sq_m4 = eval_and_free(internal_plus(
                    (Expr*[]){v_plus_sq, expr_new_integer(-4)}, 2));
                Expr* v_plus_check_arg = expand_and_free(internal_times(
                    (Expr*[]){expr_copy(r_sq), expr_copy(v_plus_sq_m4)}, 2));
                expr_free(v_plus_sq_m4);
                int v_plus_check = intsimp_sign_pos_assumption(v_plus_check_arg);
                if (v_plus_check == 0) v_plus_check = intsimp_numeric_sign(v_plus_check_arg);
                expr_free(v_plus_check_arg);

                Expr* v_minus_sq = expand_and_free(internal_power(
                    (Expr*[]){expr_copy(v_minus), expr_new_integer(2)}, 2));
                Expr* v_minus_sq_m4 = eval_and_free(internal_plus(
                    (Expr*[]){v_minus_sq, expr_new_integer(-4)}, 2));
                Expr* v_minus_check_arg = expand_and_free(internal_times(
                    (Expr*[]){expr_copy(r_sq), expr_copy(v_minus_sq_m4)}, 2));
                expr_free(v_minus_sq_m4);
                int v_minus_check = intsimp_sign_pos_assumption(v_minus_check_arg);
                if (v_minus_check == 0) v_minus_check = intsimp_numeric_sign(v_minus_check_arg);
                expr_free(v_minus_check_arg);

                if (v_plus_check >= 0 || v_minus_check >= 0) {
                    /* One of the quadratic factors would have a non-
                     * negative discriminant.  Bail. */
                    expr_free(v_plus); expr_free(v_minus);
                    expr_free(r); expr_free(r_sq);
                    goto fail;
                }

                /* Build (r^2) t^2 + (−r v) t + 1 for each v_±. */
                Expr* rv_plus = expand_and_free(internal_times(
                    (Expr*[]){expr_copy(r), expr_copy(v_plus)}, 2));
                Expr* neg_rv_plus = expand_and_free(internal_times(
                    (Expr*[]){expr_new_integer(-1), rv_plus}, 2));
                Expr* rv_minus = expand_and_free(internal_times(
                    (Expr*[]){expr_copy(r), expr_copy(v_minus)}, 2));
                Expr* neg_rv_minus = expand_and_free(internal_times(
                    (Expr*[]){expr_new_integer(-1), rv_minus}, 2));

                Expr* r_sq_c1 = expr_copy(r_sq);
                Expr* one_c1 = expr_new_integer(1);
                Expr* part1 = logtoreal_quadratic(
                    r_sq_c1,
                    neg_rv_plus,
                    one_c1,
                    s, x, t);
                expr_free(r_sq_c1); expr_free(neg_rv_plus); expr_free(one_c1);

                Expr* r_sq_c2 = expr_copy(r_sq);
                Expr* one_c2 = expr_new_integer(1);
                Expr* part2 = logtoreal_quadratic(
                    r_sq_c2,
                    neg_rv_minus,
                    one_c2,
                    s, x, t);
                expr_free(r_sq_c2); expr_free(neg_rv_minus); expr_free(one_c2);

                expr_free(v_plus); expr_free(v_minus);
                expr_free(r); expr_free(r_sq);

                if (part1 && part2) {
                    contribution = eval_and_free(internal_plus(
                        (Expr*[]){part1, part2}, 2));
                } else {
                    if (part1) expr_free(part1);
                    if (part2) expr_free(part2);
                    goto fail;
                }
            }
        }
        if (!contribution) goto fail;
        for (int m = 0; m < multiplicity; m++) {
            acc[acc_n++] = (m == 0) ? contribution : expr_copy(contribution);
        }
    }
    if (acc_n == 0) { free(acc); return expr_new_integer(0); }
    if (acc_n == 1) { Expr* r = acc[0]; free(acc); return r; }
    Expr* sum = internal_plus(acc, acc_n);
    free(acc);
    return eval_and_free(sum);

fail:
    for (size_t k = 0; k < acc_n; k++) expr_free(acc[k]);
    free(acc);
    return NULL;
}

/* LogToReal[r, s, x, t]: factor r over Q, then dispatch each linear /
 * quadratic factor through logtoreal_quadratic.  Returns NULL when r
 * has any irreducible factor of degree >= 3 over Q — the caller is
 * expected to fall back to NaiveLogPart in that case.  We deliberately
 * do NOT guess algebraic extensions here; algebraic-extension closure
 * is the responsibility of Solve (in solve.c, forthcoming). */
static Expr* intrat_log_to_real(Expr* r, Expr* s, Expr* x, Expr* t) {
    intrat_trace("LogToReal", "IN", r);

    Expr* factored = internal_factor((Expr*[]){expr_copy(r)}, 1);
    Expr* factored_e = evaluate(factored);
    expr_free(factored);
    Expr* result = logtoreal_dispatch(factored_e, s, x, t);
    expr_free(factored_e);
    if (result) { intrat_trace("LogToReal", "OUT", result); return result; }

    intrat_trace("LogToReal", "OUT (failed)", r);
    return NULL;
}

/* ----- Linear-Q closer (Phase 2's contribution to closed forms) -----
 *
 * Try to evaluate the {Q_i, S_i} pair list to a closed-form Log sum
 * by factoring each Q_i over Q[t] (or Q[t, alpha]) and contributing
 * α Log[S_i(α, x)] for every linear factor (t - α).
 *
 * Phase 2 closes the case where every linear factor lies over Q —
 * which covers `1/((x-a)(x-b))`, `1/(x^n - c^n)` for rational c, and
 * the rational-roots fragment of any LRT problem.  Quadratic factors
 * with negative or symbolic discriminants stay unhandled here; Phase
 * 4's LogToReal will pick them up.
 */

/* Walk a Factor[poly]'s Times-headed result and append every linear
 * (t - α) factor's contribution to *terms.  Returns false if any
 * factor is non-linear in t with non-trivial t-degree (caller falls
 * back to the symbolic form). */
static bool collect_linear_factors(Expr* factored, Expr* S, Expr* t,
                                    Expr*** terms, size_t* nterms,
                                    size_t* tcap) {
    Expr** factors;
    size_t nfactors;
    if (factored->type == EXPR_FUNCTION
        && factored->data.function.head->type == EXPR_SYMBOL
        && factored->data.function.head->data.symbol.name == SYM_Times) {
        factors = factored->data.function.args;
        nfactors = factored->data.function.arg_count;
    } else {
        factors = &factored;
        nfactors = 1;
    }

    for (size_t i = 0; i < nfactors; i++) {
        Expr* fac = factors[i];

        /* (t - α)^n: bind α once, contribute n times.  Power[base, k]
         * with base linear in t. */
        Expr* base = fac;
        int multiplicity = 1;
        if (fac->type == EXPR_FUNCTION
            && fac->data.function.head->type == EXPR_SYMBOL
            && fac->data.function.head->data.symbol.name == SYM_Power
            && fac->data.function.arg_count == 2
            && fac->data.function.args[1]->type == EXPR_INTEGER
            && fac->data.function.args[1]->data.integer >= 1) {
            base = fac->data.function.args[0];
            multiplicity = (int)fac->data.function.args[1]->data.integer;
        }

        int deg = get_degree_poly(base, t);
        if (deg == 0) continue;  /* constant factor */
        if (deg != 1) return false;

        Expr* alpha = extract_linear_root(base, t);
        if (!alpha) return false;
        Expr* contribution = build_log_term(alpha, S, t);
        expr_free(alpha);
        for (int m = 0; m < multiplicity; m++) {
            if (*nterms >= *tcap) {
                *tcap *= 2;
                *terms = (Expr**)realloc(*terms, sizeof(Expr*) * (*tcap));
            }
            (*terms)[(*nterms)++] = (m == 0) ? contribution : expr_copy(contribution);
        }
    }
    return true;
}

/* Try to evaluate the {Q, S} pair list to a closed-form Log sum
 * by factoring each Q_i over Q.  Returns NULL when any Q_i has
 * factors of degree > 1 (Phase 4 territory). */
static Expr* intrat_linear_q_closer(Expr* pair_list, Expr* x, Expr* t) {
    (void)x;
    if (!pair_list || pair_list->type != EXPR_FUNCTION
        || pair_list->data.function.head->type != EXPR_SYMBOL
        || pair_list->data.function.head->data.symbol.name != SYM_List) return NULL;
    size_t n = pair_list->data.function.arg_count;
    /* Empty pair list = LRT producer couldn't decompose.  See the
     * matching guard in intrat_log_to_real_pairs above for why returning
     * Integer[0] would be a silent wrong-answer. */
    if (n == 0) return NULL;

    size_t tcap = 8, nterms = 0;
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * tcap);

    for (size_t i = 0; i < n; i++) {
        Expr* pair = pair_list->data.function.args[i];
        Expr* Qi = list_get(pair, 1);
        Expr* Si = list_get(pair, 2);
        if (!Qi || !Si) {
            for (size_t k = 0; k < nterms; k++) expr_free(terms[k]);
            free(terms);
            if (Qi) expr_free(Qi);
            if (Si) expr_free(Si);
            return NULL;
        }

        Expr* factored = internal_factor((Expr*[]){expr_copy(Qi)}, 1);
        Expr* factored_e = evaluate(factored);
        expr_free(factored);

        bool ok = collect_linear_factors(factored_e, Si, t, &terms, &nterms, &tcap);
        expr_free(factored_e);
        expr_free(Qi); expr_free(Si);
        if (!ok) {
            for (size_t k = 0; k < nterms; k++) expr_free(terms[k]);
            free(terms);
            return NULL;
        }
    }

    Expr* sum;
    if (nterms == 0) { sum = expr_new_integer(0); free(terms); }
    else if (nterms == 1) { sum = terms[0]; free(terms); }
    else { sum = internal_plus(terms, nterms); free(terms); sum = eval_and_free(sum); }
    return sum;
}

/* ====================================================================
 * Phase 6 + 7 — output simplification of the resulting integral
 * ====================================================================
 *
 * Log pairing (c Log[A] ± c Log[B] -> Log[A·B] / Log[A/B] / ArcTanh),
 * Log[c·p] -> Log[p] when c is free of x, Plus distribution, and
 * ArcTan / ArcTanh sign normalisation all live in intsimp.c — see
 * intsimp.h for the public surface.
 */


/* ====================================================================
 * Phase 5 — IntegrateRealRationalFunction-style per-summand loop.
 * ====================================================================
 *
 * Apart-split a squarefree-denominator residual h, ExtractConstants
 * each summand to keep coefficients out of the subresultant chain,
 * and process each piece independently.  If *every* piece closes via
 * the derivative-recognition / LRT / LogToReal pipeline we return
 * the sum; otherwise NULL so the caller can attempt the whole-h LRT
 * path (or surface the integral as unevaluated).
 *
 * Direct port of IntegrateRational.m:534-587, minus the `NaiveLogPart`
 * fallback (we prefer to surface the integral as unevaluated rather
 * than emit a symbolic RootSum that Mathilda can't yet differentiate
 * back to the integrand).
 */

/* LRT consumer: apply Rioboo's intrat_log_to_real to every {Q_i, S_i}
 * pair and sum the results.  Returns the closed real-elementary form
 * on success, or NULL if any pair's Q_i has an irreducible factor of
 * degree >= 5 (intrat_log_to_real's failure mode) or any list entry
 * is malformed.  Caller retains ownership of `pair_list`. */
static Expr* intrat_log_to_real_pairs(Expr* pair_list, Expr* x, Expr* t) {
    if (!pair_list || pair_list->type != EXPR_FUNCTION
        || pair_list->data.function.head->type != EXPR_SYMBOL
        || pair_list->data.function.head->data.symbol.name != SYM_List) {
        return NULL;
    }
    size_t npairs = pair_list->data.function.arg_count;
    /* Empty pair list means the LRT producer couldn't decompose the
     * input -- typically the denominator was irreducible over the
     * working ring (e.g. x^16 - x^8 - 1, an irreducible biquadratic-
     * in-x^8 with sqrt-5 roots).  Returning Integer[0] here would
     * make the caller treat the integral as identically zero --
     * silent wrong-answer.  Return NULL so try_lrt_close falls
     * through to NaiveLogPart's universal RootSum form. */
    if (npairs == 0) return NULL;
    Expr** lr_terms = (Expr**)malloc(sizeof(Expr*) * npairs);
    size_t lr_count = 0;

    for (size_t i = 0; i < npairs; i++) {
        Expr* pair = pair_list->data.function.args[i];
        Expr* Qi = list_get(pair, 1);
        Expr* Si = list_get(pair, 2);
        if (!Qi || !Si) {
            if (Qi) expr_free(Qi);
            if (Si) expr_free(Si);
            for (size_t k = 0; k < lr_count; k++) expr_free(lr_terms[k]);
            free(lr_terms);
            return NULL;
        }
        Expr* term = intrat_log_to_real(Qi, Si, x, t);
        expr_free(Qi); expr_free(Si);
        if (!term) {
            for (size_t k = 0; k < lr_count; k++) expr_free(lr_terms[k]);
            free(lr_terms);
            return NULL;
        }
        lr_terms[lr_count++] = term;
    }

    Expr* result;
    if (lr_count == 0) {
        result = expr_new_integer(0);
    } else if (lr_count == 1) {
        result = lr_terms[0];
        lr_terms[0] = NULL;  /* prevent double-free */
    } else {
        Expr* sum = internal_plus(lr_terms, lr_count);
        memset(lr_terms, 0, sizeof(Expr*) * lr_count);
        result = eval_and_free(sum);
    }
    free(lr_terms);
    return result;
}

/* Try to close `num/den` via the LRT-then-real-form pipeline.
 *
 * Architecture (mirrors AXIOM's INTRAT/IRRF separation in
 * intrf.spad.pamphlet:776-778, plus irexpand):
 *   - Producer:  intrat_int_rational_log_part -- single LRT call per
 *                try_lrt_close, emits the {Q_i, S_i} pair list.
 *   - LRT consumers (uniform `(pairs, x, t) -> Expr*|NULL` shape):
 *       * intrat_log_to_real_pairs  -- Rioboo real-elementary form
 *       * intrat_linear_q_closer    -- every Q_i splits over Q
 *   - Non-LRT fallback:  intrat_naive_log_part(num/den) -- Lagrange
 *     RootSum form, universal so try_lrt_close never returns NULL
 *     for a well-formed proper rational integrand.
 *
 * Closure preference matches IntegrateRational.m:560-572. */
static Expr* try_lrt_close(Expr* num, Expr* den, Expr* x) {
    if (is_zero_poly(num)) return expr_new_integer(0);

    Expr* h = internal_divide(
        (Expr*[]){expr_copy(num), expr_copy(den)}, 2);
    Expr* h_eval = eval_and_free(h);
    Expr* t_sym = expr_new_symbol("Integrate`Private`tt$");

    /* Single LRT producer call. */
    Expr* pairs = intrat_int_rational_log_part(h_eval, x, t_sym, false);

    Expr* result = NULL;
    if (pairs) {
        /* LRT consumer 1: Rioboo LogToReal on every pair. */
        result = intrat_log_to_real_pairs(pairs, x, t_sym);
        /* LRT consumer 2: linear-Q closer (factor each Q_i over Q). */
        if (!result) result = intrat_linear_q_closer(pairs, x, t_sym);
        expr_free(pairs);
    }
    expr_free(t_sym);

    /* Non-LRT universal fallback (Lagrange RootSum on full a/d). */
    if (!result) result = intrat_naive_log_part(h_eval, x);

    expr_free(h_eval);
    return result;
}

/* Process h piece-by-piece via Apart, ExtractConstants, and the
 * derivative-recognition / LRT pipeline.  Returns the closed sum or
 * NULL if any piece can't be closed. */
static Expr* intrat_integrate_summands(Expr* h, Expr* x) {
    Expr* pieces = intrat_apart_list(h, x, NULL);
    if (!pieces || pieces->type != EXPR_FUNCTION
        || pieces->data.function.head->type != EXPR_SYMBOL
        || pieces->data.function.head->data.symbol.name != SYM_List) {
        if (pieces) expr_free(pieces);
        return NULL;
    }
    size_t n = pieces->data.function.arg_count;
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    size_t nterms = 0;

    for (size_t i = 0; i < n; i++) {
        Expr* piece = pieces->data.function.args[i];

        /* Skip pure polynomial summands — they belong to the Q part. */
        if (intrat_polyq_test(piece, x)) {
            Expr* poly_int = intrat_integrate_polynomial(piece, x);
            if (!poly_int) goto fail;
            terms[nterms++] = poly_int;
            continue;
        }

        /* ExtractConstants — pull scalar prefactors out of num/den. */
        Expr* ec = intrat_extract_constants(piece, x);
        Expr* c_const = list_get(ec, 1);
        Expr* residue = list_get(ec, 2);
        expr_free(ec);
        if (!c_const || !residue) {
            if (c_const) expr_free(c_const);
            if (residue) expr_free(residue);
            goto fail;
        }

        Expr* num = intrat_numerator(residue);
        Expr* den = intrat_denominator(residue);
        Expr* num_e = eval_and_free(num);
        Expr* den_e = eval_and_free(den);
        expr_free(residue);

        /* Derivative-recognition fast path. */
        Expr* fast = intrat_derivative_recognition(num_e, den_e, x);
        Expr* piece_int = NULL;
        if (fast) piece_int = fast;
        else      piece_int = try_lrt_close(num_e, den_e, x);
        expr_free(num_e); expr_free(den_e);

        if (!piece_int) {
            expr_free(c_const);
            goto fail;
        }
        Expr* scaled = internal_times(
            (Expr*[]){c_const, piece_int}, 2);
        terms[nterms++] = eval_and_free(scaled);
    }

    expr_free(pieces);
    if (nterms == 0) { free(terms); return expr_new_integer(0); }
    if (nterms == 1) { Expr* r = terms[0]; free(terms); return r; }
    Expr* sum = internal_plus(terms, nterms);
    free(terms);
    return eval_and_free(sum);

fail:
    expr_free(pieces);
    for (size_t k = 0; k < nterms; k++) expr_free(terms[k]);
    free(terms);
    return NULL;
}

/* ------------------------------------------------------------------ */
/* BronsteinRational top-level (Phase 1 skeleton + Phase 2 hookup).    */
/* ------------------------------------------------------------------ */

static Expr* intrat_integrate_rational(Expr* f, Expr* x) {
    intrat_trace("BronsteinRational", "IN", f);

    /* Pure-polynomial fast path. */
    if (intrat_polyq_test(f, x)) {
        return intrat_integrate_polynomial(f, x);
    }

    /* Combine into a single fraction num/den. */
    Expr* tg = internal_together((Expr*[]){expr_copy(f)}, 1);
    Expr* combined = evaluate(tg);
    expr_free(tg);
    if (!combined) return NULL;

    Expr* num = intrat_numerator(combined);
    Expr* num_eval = eval_and_free(num);
    Expr* den = intrat_denominator(combined);
    Expr* den_eval = eval_and_free(den);
    expr_free(combined);

    if (!intrat_polyq_test(num_eval, x) || !intrat_polyq_test(den_eval, x)) {
        expr_free(num_eval); expr_free(den_eval);
        return NULL;
    }

    /* Polynomial part: Q + R/den. */
    Expr *Q = NULL, *R = NULL;
    if (!intrat_polyqr(num_eval, den_eval, x, &Q, &R)) {
        expr_free(num_eval); expr_free(den_eval);
        return NULL;
    }
    expr_free(num_eval);

    Expr* poly_int = expr_new_integer(0);
    if (!is_zero_poly(Q)) {
        Expr* tmp = intrat_integrate_polynomial(Q, x);
        if (tmp) { expr_free(poly_int); poly_int = tmp; }
    }
    expr_free(Q);

    if (is_zero_poly(R)) {
        expr_free(R); expr_free(den_eval);
        Expr* result = eval_and_free(poly_int);
        intrat_trace("BronsteinRational", "OUT", result);
        return result;
    }

    /* Try the derivative-recognition fast path on R/den. */
    Expr* fast = intrat_derivative_recognition(R, den_eval, x);
    if (fast) {
        Expr* sum = internal_plus(
            (Expr*[]){poly_int, fast}, 2);
        Expr* result = eval_and_free(sum);
        expr_free(R); expr_free(den_eval);
        intrat_trace("BronsteinRational", "OUT", result);
        return result;
    }

    /* HermiteReduce R/den. */
    Expr* proper = internal_divide(
        (Expr*[]){R, den_eval}, 2);
    Expr* hr = intrat_hermite_reduce(proper, x);
    expr_free(proper);
    if (!hr || hr->type != EXPR_FUNCTION || hr->data.function.arg_count != 2) {
        if (hr) expr_free(hr);
        expr_free(poly_int);
        return NULL;
    }
    Expr* g = expr_copy(hr->data.function.args[0]);
    Expr* h = expr_copy(hr->data.function.args[1]);
    expr_free(hr);

    /* If h is zero, we have a closed form. */
    if (is_zero_poly(h)) {
        expr_free(h);
        Expr* sum = internal_plus(
            (Expr*[]){poly_int, g}, 2);
        Expr* result = eval_and_free(sum);
        intrat_trace("BronsteinRational", "OUT", result);
        return result;
    }

    /* Phase 5 — IntegrateRealRationalFunction-style per-summand loop.
     *
     * We process the squarefree-denominator residual h piece-by-piece:
     * Apart-split it, ExtractConstants to keep coefficients out of the
     * subresultant chain, then run derivative-recognition / LRT /
     * LogToReal on each piece.  If every piece closes we sum them up
     * with the Hermite g and the polynomial part; otherwise we fall
     * back to running LRT on h whole (Phase 4 path).
     */
    Expr* per_summand = intrat_integrate_summands(h, x);
    if (per_summand) {
        expr_free(h);
        Expr* sum = internal_plus(
            (Expr*[]){poly_int, g, per_summand}, 3);
        Expr* sum_eval = eval_and_free(sum);
        /* Phase A-1: distribute scalar Times-over-Plus so per-summand
         * `c_k * piece_int_k` accumulators flatten into a single
         * Plus[...] of `c_k · term` summands.  Without this, the
         * subsequent log-pairing rule (intsimp_log_to_arctanh) cannot
         * see across the Times boundary.
         * Phase A-2: strip Log[c · p] -> Log[p] when FreeQ[c, x] so
         * downstream Log[A] + Log[B] merging sees the cleaned
         * arguments. */
        Expr* distributed = intsimp_distribute_plus(sum_eval);
        expr_free(sum_eval);
        Expr* stripped = intsimp_strip_log_constants(distributed, x);
        expr_free(distributed);
        /* Phase 6 final pass: combine c Log[A] ± c Log[B] terms.
         * Phase 7 final pass: normalise ArcTan/ArcTanh argument signs.*/
        Expr* simplified = intsimp_log_to_arctanh(stripped, x);
        expr_free(stripped);
        /* Re-distribute and re-strip in case the pairing rule
         * introduces new Times[scalar, Plus[…]] or Log[c · p] shapes. */
        Expr* simplified2 = intsimp_distribute_plus(simplified);
        expr_free(simplified);
        Expr* simplified3 = intsimp_strip_log_constants(simplified2, x);
        expr_free(simplified2);
        Expr* result = intsimp_normalize_inverse_trig_signs(simplified3);
        expr_free(simplified3);
        intrat_trace("BronsteinRational", "OUT", result);
        return result;
    }

    /* Beyond Phase 5's combined scope — leave unevaluated. */
    expr_free(h); expr_free(g); expr_free(poly_int);
    intrat_trace("BronsteinRational", "OUT (unresolved)", f);
    return NULL;
}

/* ------------------------------------------------------------------ */
/* BuiltinFunc wrappers.  Each follows the Mathilda convention: caller  */
/* owns res; return NULL on failure.                                   */
/* ------------------------------------------------------------------ */

/* Recognise a trailing option (Rule[lhs, rhs] / RuleDelayed[...]).
 * Returns true if `opt` is a Rule with one of the BronsteinRational
 * option names on the lhs. */
static bool is_intrat_option(Expr* opt) {
    if (!opt || opt->type != EXPR_FUNCTION
        || opt->data.function.head->type != EXPR_SYMBOL
        || (opt->data.function.head->data.symbol.name != SYM_Rule
            && opt->data.function.head->data.symbol.name != SYM_RuleDelayed)
        || opt->data.function.arg_count != 2) return false;
    Expr* lhs = opt->data.function.args[0];
    if (lhs->type == EXPR_SYMBOL) {
        return lhs->data.symbol.name == SYM_Extension;
    }
    if (lhs->type == EXPR_STRING) {
        const char* s = lhs->data.string;
        return strcmp(s, "PFD") == 0 || strcmp(s, "LogToArcTan") == 0
            || strcmp(s, "Radicals") == 0;
    }
    return false;
}

/* Recursively scan for any EXPR_REAL leaf.  Mirrors the same guard in
 * src/integrate.c; here it is silent (no Integrate::inexact message)
 * so the user-visible diagnostic comes from the public Integrate head. */
static bool intrat_contains_real(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
    if (e->type == EXPR_FUNCTION) {
        if (intrat_contains_real(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (intrat_contains_real(e->data.function.args[i])) return true;
        }
    }
    return false;
}

Expr* builtin_intrat_integraterational(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    /* Strip recognised trailing options.  Phase 7 stores option
     * values via extract_bool_option / extract_extension_option but
     * keeps them advisory for now — the algorithmic path uses the
     * defaults. */
    size_t argc = res->data.function.arg_count;
    while (argc > 0 && is_intrat_option(res->data.function.args[argc - 1])) {
        argc--;
    }
    if (argc != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    /* Try to coerce inexact integrands to exact rationals.  Silent
     * here — the public Integrate head owns the user-visible
     * diagnostic.  If Rationalize cannot eliminate every EXPR_REAL,
     * bubble back unevaluated. */
    Expr* coerced = NULL;
    if (intrat_contains_real(f)) {
        coerced = internal_rationalize_expr(f, 0.0, RATIONALIZE_DEFAULT);
        if (!coerced || intrat_contains_real(coerced)) {
            if (coerced) expr_free(coerced);
            return NULL;
        }
        f = coerced;
    }
    Expr* result = intrat_integrate_rational(f, x);
    if (coerced) expr_free(coerced);
    return result;
}

Expr* builtin_intrat_integratepolynomial(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    if (!intrat_polyq_test(f, x)) return NULL;
    return intrat_integrate_polynomial(f, x);
}

Expr* builtin_intrat_hermitereduce(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    return intrat_hermite_reduce(f, x);
}

Expr* builtin_intrat_helpers_content(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* p = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    if (is_zero_poly(p)) return expr_new_integer(0);
    return intrat_content(p, x);
}

Expr* builtin_intrat_helpers_primitive(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* p = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    return intrat_primitive(p, x);
}

Expr* builtin_intrat_helpers_monic(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* p = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    return intrat_monic_internal(p, x);
}

Expr* builtin_intrat_helpers_leadingcoefficient(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* p = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    return intrat_lc(p, x);
}

/* Phase 2 builtins. */

Expr* builtin_intrat_intrationallogpart(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    /* Strip a trailing RootSum -> True/False option. */
    size_t argc = res->data.function.arg_count;
    bool root_sum = false;
    if (argc >= 1) {
        Expr* last = res->data.function.args[argc - 1];
        if (last->type == EXPR_FUNCTION
            && last->data.function.head->type == EXPR_SYMBOL
            && (last->data.function.head->data.symbol.name == SYM_Rule
                || last->data.function.head->data.symbol.name == SYM_RuleDelayed)
            && last->data.function.arg_count == 2
            && last->data.function.args[0]->type == EXPR_SYMBOL
            && last->data.function.args[0]->data.symbol.name == SYM_RootSum) {
            Expr* val = last->data.function.args[1];
            root_sum = (val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_True);
            argc--;
        }
    }
    if (argc != 3) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    Expr* t = res->data.function.args[2];
    if (x->type != EXPR_SYMBOL || t->type != EXPR_SYMBOL) return NULL;
    return intrat_int_rational_log_part(f, x, t, root_sum);
}

Expr* builtin_intrat_helpers_squarefree(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    return intrat_squarefree_list(res->data.function.args[0]);
}

Expr* builtin_intrat_helpers_extractconstants(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    return intrat_extract_constants(f, x);
}

Expr* builtin_intrat_helpers_apartlist(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    const Expr* alpha = extract_extension_option(res, &argc);
    if (argc != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    return intrat_apart_list(f, x, alpha);
}

Expr* builtin_intrat_logtoatan(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    Expr* a = res->data.function.args[0];
    Expr* b = res->data.function.args[1];
    Expr* x = res->data.function.args[2];
    if (x->type != EXPR_SYMBOL) return NULL;
    return intrat_log_to_atan(a, b, x);
}

Expr* builtin_intrat_logtoreal(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 4) return NULL;
    Expr* r = res->data.function.args[0];
    Expr* s = res->data.function.args[1];
    Expr* x = res->data.function.args[2];
    Expr* t = res->data.function.args[3];
    if (x->type != EXPR_SYMBOL || t->type != EXPR_SYMBOL) return NULL;
    return intrat_log_to_real(r, s, x, t);
}

Expr* builtin_intrat_transcendental_log_part(Expr* res) {
    if (res->type != EXPR_FUNCTION
        || (res->data.function.arg_count != 6 && res->data.function.arg_count != 7))
        return NULL;
    Expr* a    = res->data.function.args[0];
    Expr* d    = res->data.function.args[1];
    Expr* tau  = res->data.function.args[2];
    Expr* z    = res->data.function.args[3];
    Expr* Dd   = res->data.function.args[4];
    Expr* xreal = res->data.function.args[5];
    /* Optional 7th argument is the "decide" marker (any expression): when present
     * the residue criterion, on an all-non-constant-residue resultant, returns
     * Integrate`PartialLogPart[0, a/d] instead of declining — so a caller seeking a
     * non-elementarity VERDICT (never re-integrating the remainder) learns the whole
     * simple part is non-elementary (Thm 5.6.1(ii)).  Absent for the integrator. */
    bool decide = (res->data.function.arg_count == 7);
    if (tau->type != EXPR_SYMBOL || z->type != EXPR_SYMBOL)
        return NULL;
    /* `xreal` is the residue-constant gate: a single symbol (single-kernel
     * extension) or a List of symbols (tower proper part — free of x and every
     * lower tower variable). */
    bool gate_ok = (xreal->type == EXPR_SYMBOL)
        || (xreal->type == EXPR_FUNCTION
            && xreal->data.function.head->type == EXPR_SYMBOL
            && xreal->data.function.head->data.symbol.name == SYM_List);
    if (!gate_ok) return NULL;
    Expr* rem = NULL;
    Expr* logs = intrat_transcendental_log_part(a, d, tau, z, Dd, xreal, &rem, decide);
    if (!logs) { if (rem) expr_free(rem); return NULL; }
    if (rem) {
        /* Partial log part: wrap the elementary logs and the unintegrated
         * non-constant-residue remainder in an inert head the recursive-Risch
         * caller (rt_field_lrt_logpart) decodes.  The head has no definition, so
         * it survives evaluation unchanged. */
        return expr_new_function(expr_new_symbol("Integrate`PartialLogPart"),
                                 (Expr*[]){ logs, rem }, 2);
    }
    return logs;
}

Expr* builtin_intrat_logtoarctanh(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* e = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    return intsimp_log_to_arctanh(e, x);
}

Expr* builtin_intrat_naive_log_part(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    return intrat_naive_log_part(f, x);
}

/* Integrate`RothsteinTragerResultant[num, den, z, x]
 *   = Resultant[num - z D[den, x], den, x],
 * the Rothstein-Trager resultant whose roots in z are the residues of num/den.
 * This is the named parametric-resultant primitive the LRT log part uses and that
 * Cherry's special-function argument generators reuse (CHERRY_DESIGN.md §3.3): the
 * log / exponential-integral (Ei) arguments are recovered from its roots. Built on
 * the existing `Resultant` builtin — no new resultant machinery. */
static Expr* builtin_intrat_rt_resultant(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 4) return NULL;
    Expr* num = res->data.function.args[0];
    Expr* den = res->data.function.args[1];
    Expr* z   = res->data.function.args[2];
    Expr* x   = res->data.function.args[3];
    Expr* dd = eval_and_free(expr_new_function(expr_new_symbol(SYM_D),
        (Expr*[]){ expr_copy(den), expr_copy(x) }, 2));            /* D[den, x] */
    if (!dd) return NULL;
    Expr* arg = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy(num),
            expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ expr_new_integer(-1), expr_copy(z), dd }, 3) }, 2)); /* num - z D[den] */
    Expr* r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Resultant),
        (Expr*[]){ arg, expr_copy(den), expr_copy(x) }, 3));       /* Resultant[.,den,x] */
    return r;
}

/* Integrate`ExpIntegralEiResultant[g1, p, q, alpha, x]
 *   = Resultant[g1, p + alpha q, x],
 * the parametric resultant whose roots in the parameter alpha are Cherry's
 * ei-argument constants (the u_i = f + alpha_i of the P1 generator, resolving
 * g~1/g1 = Sum c_i u_i'/u_i; Cherry 1989 §4.2 / CHERRY_PLAN.md §3.2).  Distinct
 * from RothsteinTragerResultant (which resultants against num - z D[den]); here
 * the parameter enters linearly through p + alpha q with q the exponent
 * denominator.  A thin, testable sibling built on the existing Resultant. */
static Expr* builtin_intrat_ei_resultant(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 5) return NULL;
    Expr* g1    = res->data.function.args[0];
    Expr* p     = res->data.function.args[1];
    Expr* q     = res->data.function.args[2];
    Expr* alpha = res->data.function.args[3];
    Expr* x     = res->data.function.args[4];
    Expr* arg = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy(p),
            expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ expr_copy(alpha), expr_copy(q) }, 2) }, 2));  /* p + alpha q */
    Expr* r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Resultant),
        (Expr*[]){ expr_copy(g1), arg, expr_copy(x) }, 3));             /* Resultant[g1,.,x] */
    return r;
}

/* ------------------------------------------------------------------ */
/* Init.                                                              */
/* ------------------------------------------------------------------ */

static void install(const char* name, Expr* (*fn)(Expr*), const char* docstring) {
    symtab_add_builtin(name, fn);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    if (docstring) symtab_set_docstring(name, docstring);
}

void intrat_init(void) {
    /* Package symbols live under `Integrate``. Names must be the fully
     * qualified strings — the symtab is keyed by interned name and the
     * parser routes `Integrate`Foo` directly to that key. */
    install("Integrate`BronsteinRational",
            builtin_intrat_integraterational,
            "Integrate`BronsteinRational[f, x] is the explicit form of the\n"
            "rational-function integrator. Equivalent to Integrate[f, x] when\n"
            "f is a polynomial or rational function in x. Phase 1 closes the\n"
            "polynomial / Hermite / derivative-recognition cases; non-rational\n"
            "or LRT-required inputs return unevaluated.");

    install("Integrate`IntegratePolynomial",
            builtin_intrat_integratepolynomial,
            "Integrate`IntegratePolynomial[p, x] integrates a polynomial p in x\n"
            "term-by-term using the rule a x^n -> a x^(n+1)/(n+1). Returns\n"
            "unevaluated for non-polynomial inputs.");

    install("Integrate`HermiteReduce",
            builtin_intrat_hermitereduce,
            "Integrate`HermiteReduce[f, x] returns {g, h} such that\n"
            "f == D[g, x] + h with Denominator[h] squarefree (Mack's linear\n"
            "version of Hermite reduction; Bronstein, Symbolic Integration I,\n"
            "p. 44).");

    install("Integrate`Helpers`Content",
            builtin_intrat_helpers_content,
            "Integrate`Helpers`Content[p, x] gives PolynomialGCD applied to\n"
            "the coefficient list of p as a polynomial in x.");

    install("Integrate`Helpers`Primitive",
            builtin_intrat_helpers_primitive,
            "Integrate`Helpers`Primitive[p, x] is p divided by Content[p, x],\n"
            "the primitive-part operation in the polynomial ring K[x].");

    install("Integrate`Helpers`Monic",
            builtin_intrat_helpers_monic,
            "Integrate`Helpers`Monic[p, x] divides p by its leading coefficient\n"
            "in x, returning a monic polynomial. For p free of x returns 1.");

    install("Integrate`Helpers`LeadingCoefficient",
            builtin_intrat_helpers_leadingcoefficient,
            "Integrate`Helpers`LeadingCoefficient[p, x] gives the coefficient\n"
            "of the highest power of x in p; for p free of x returns p.");

    /* Phase 2 — Lazard-Rioboo-Trager log part. */
    install("Integrate`IntRationalLogPart",
            builtin_intrat_intrationallogpart,
            "Integrate`IntRationalLogPart[A/D, x, t] computes the logarithmic\n"
            "part of an integral of a rational function with squarefree\n"
            "denominator D, returning a list of {Q_i(t), S_i(t,x)} pairs\n"
            "(default) or a sum of RootSum heads (option RootSum -> True).\n"
            "Bronstein, Symbolic Integration I, p. 51 (Lazard-Rioboo-Trager).");

    install("Integrate`Helpers`SquareFree",
            builtin_intrat_helpers_squarefree,
            "Integrate`Helpers`SquareFree[p] returns a list of\n"
            "{factor, multiplicity} pairs indexed densely by multiplicity\n"
            "1..max, with the i-th entry holding the product of all\n"
            "squarefree factors of multiplicity i (or {1, i} if none).");

    install("Integrate`Helpers`ExtractConstants",
            builtin_intrat_helpers_extractconstants,
            "Integrate`Helpers`ExtractConstants[f, x] returns {const, residue}\n"
            "with const free of x and residue holding the remaining\n"
            "polynomial-in-x parts of Numerator[f] and Denominator[f].");

    install("Integrate`Helpers`ApartList",
            builtin_intrat_helpers_apartlist,
            "Integrate`Helpers`ApartList[f, x] returns Apart[f, x]'s output\n"
            "as a List of summands. Accepts a trailing Extension -> alpha\n"
            "rule that is forwarded to the underlying Apart call.");

    /* Phase 3 — Rioboo's recursive LogToAtan. */
    install("Integrate`LogToAtan",
            builtin_intrat_logtoatan,
            "Integrate`LogToAtan[A, B, x] returns a sum of ArcTan of polynomials\n"
            "in K[x] whose derivative equals D[I Log[(A + I B)/(A - I B)], x].\n"
            "Used by LogToReal (Phase 4) to convert complex log pairs into real\n"
            "ArcTan summands. Bronstein, Symbolic Integration I, p. 63 (Rioboo).");

    /* Phase 6 — Log post-processing. */
    install("Integrate`LogToArcTanh",
            builtin_intrat_logtoarctanh,
            "Integrate`LogToArcTanh[expr, x] combines pairs of Log terms in\n"
            "expr's top-level Plus head: c Log[A] + c Log[B] -> c Log[A*B];\n"
            "c Log[A] - c Log[B] -> 2 c ArcTanh[(B-A)/(B+A)] when the\n"
            "ArcTanh argument is rational in x, else c Log[A/B] when A/B is\n"
            "rational in x. Beautifies the rational integrator's output\n"
            "without affecting differentiation; cf. IntegrateRational.m\n"
            "lines 1722-1761.");

    /* Phase 4 — Rioboo's LogToReal. */
    install("Integrate`LogToReal",
            builtin_intrat_logtoreal,
            "Integrate`LogToReal[r, s, x, t] converts the symbolic\n"
            "RootSum[r(α)==0, α Log[s(α, x)]] into a real Log + ArcTan form\n"
            "by factoring r over Q and dispatching each linear / quadratic\n"
            "factor (positive discriminant -> two real Logs; negative\n"
            "discriminant -> Log + LogToAtan-derived ArcTan). Returns the\n"
            "call unevaluated when r has factors of degree > 2 in t (the\n"
            "bounded-Solve scope of Phase 4).");

    /* Recursive-Risch transcendental LRT log part (consumed by
     * Integrate`RischTranscendental; not a user-facing surface form). */
    install("Integrate`RothsteinTragerResultant",
            builtin_intrat_rt_resultant,
            "Integrate`RothsteinTragerResultant[num, den, z, x] = Resultant[num - z D[den,x],\n"
            "den, x], the Rothstein-Trager resultant whose roots in z are the residues of\n"
            "num/den. Argument generator for the log / exponential-integral (Ei) parts of\n"
            "Cherry's special-function integration.");

    install("Integrate`ExpIntegralEiResultant",
            builtin_intrat_ei_resultant,
            "Integrate`ExpIntegralEiResultant[g1, p, q, a, x] = Resultant[g1, p + a q, x],\n"
            "the parametric resultant whose roots in the parameter a are the ExpIntegralEi\n"
            "argument constants alpha_i of Cherry's P1 generator (u_i = f + alpha_i, f = p/q).\n"
            "Built on Resultant; the ei sibling of RothsteinTragerResultant.");

    install("Integrate`TranscendentalLogPart",
            builtin_intrat_transcendental_log_part,
            "Integrate`TranscendentalLogPart[a, d, tau, z, Dd, g] computes the\n"
            "logarithmic part of Integrate[a/d] for a transcendental monomial\n"
            "tau (a/d proper with squarefree d in tau), where Dd is the monomial\n"
            "derivation D(d), z is the residue variable, and g is the residue-\n"
            "constant gate: a symbol (single kernel) or a List of symbols (tower\n"
            "proper part, gating x and every lower tower variable). Returns the\n"
            "real Log + ArcTan form as a function of tau (Rothstein-Trager\n"
            "resultant + Rioboo LogToReal), or unevaluated when the integral is\n"
            "not elementary in this form.");

    /* Phase 8b — NaiveLogPart RootSum fallback. */
    install("Integrate`NaiveLogPart",
            builtin_intrat_naive_log_part,
            "Integrate`NaiveLogPart[f, x] returns the held-symbolic\n"
            "RootSum form of the logarithmic part of Integrate[f, x]\n"
            "for a rational function f = a/d:\n"
            "  RootSum[Function[t, d(t)],\n"
            "          Function[t, a(t) Log(x - t) / d'(t)]].\n"
            "Used as the universal fallback by the rational integrator\n"
            "whenever LogToReal cannot close the log part to a real\n"
            "elementary expression.  The result differentiates back to f\n"
            "via the D[RootSum, x] rule wired in src/deriv.c.\n"
            "Direct port of IntegrateRational.m:1116-1124.");

    /* Trace flag: Integrate`$Verbose. Default False. */
    {
        Expr* pat = expr_new_symbol("Integrate`$Verbose");
        Expr* val = expr_new_symbol(SYM_False);
        symtab_add_own_value("Integrate`$Verbose", pat, val);
        expr_free(pat); expr_free(val);
    }

    /* RootSum is registered in src/root.c (Phase 8b-prereq) — the
     * canonical home of held algebraic-root machinery.  Nothing to
     * do here. */
}
