/* intrat.c — rational-function integration package (`Integrate``).
 *
 * Phase 1 of INTEGRATE_PLAN.md.  Implements the scaffolding,
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
 * Memory contract: every public BuiltinFunc follows the picocas
 * convention — caller (evaluator) owns `res`; the function returns a
 * freshly-allocated Expr* on success or NULL on failure.
 */

#include "intrat.h"
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

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
            && r->replacement->data.symbol == SYM_True) return true;
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
/* Small wrappers around picocas internal_* / evaluate idioms.        */
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

/* canonic[expr] = Cancel[Together[expr]].  Phase 1 keeps this
 * stripped down — the full Mathematica `canonic` adds a RootReduce /
 * ToRadicals tail that we'll bring online in later phases. */
static Expr* intrat_canonic(Expr* e) {
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

static Expr* intrat_denominator(Expr* e) {
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
        expr_new_symbol("PolynomialQuotientRemainder"), args, 3);
    Expr* result = evaluate(call);
    expr_free(call);
    if (!result || result->type != EXPR_FUNCTION
        || result->data.function.head->type != EXPR_SYMBOL
        || result->data.function.head->data.symbol != SYM_List
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
    bool ok = (val && val->type == EXPR_SYMBOL && val->data.symbol == SYM_True);
    if (val) expr_free(val);
    return ok;
}

static bool intrat_freeq_test(Expr* expr, Expr* var) {
    Expr* args[2] = {expr_copy(expr), expr_copy(var)};
    Expr* call = internal_freeq(args, 2);
    Expr* val  = evaluate(call);
    expr_free(call);
    bool ok = (val && val->type == EXPR_SYMBOL && val->data.symbol == SYM_True);
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
        || cl->data.function.head->data.symbol != SYM_List) {
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
        || coefs->data.function.head->data.symbol != SYM_List
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

    Expr* result = expr_new_function(expr_new_symbol("List"),
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
        || cl->data.function.head->data.symbol != SYM_List) {
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
        && fs_expanded->data.function.head->data.symbol == SYM_Times) {
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
                && arg->data.function.head->data.symbol == SYM_Power
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
        && fs_expanded->data.function.head->data.symbol == SYM_Power
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
        Expr* log_pol = expr_new_function(expr_new_symbol("Log"),
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
 * and the IntegrateRational test suite.
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
 * {1, i}.  Mirrors the IntegrateRational baseline's `SquareFree[]`
 * helper at IntegrateRational.m:1474-1487.
 *
 * Implementation note: picocas's FactorSquareFree returns the
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
        && fs_eval->data.function.head->data.symbol == SYM_Times) {
        for (size_t i = 0; i < fs_eval->data.function.arg_count; i++) {
            Expr* arg = fs_eval->data.function.args[i];
            if (arg->type == EXPR_FUNCTION
                && arg->data.function.head->type == EXPR_SYMBOL
                && arg->data.function.head->data.symbol == SYM_Power
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
        && fs_eval->data.function.head->data.symbol == SYM_Power
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
        Expr* pair = expr_new_function(expr_new_symbol("List"),
            (Expr*[]){expr_copy(p), expr_new_integer(1)}, 2);
        Expr* list = expr_new_function(expr_new_symbol("List"),
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
        out_args[m - 1] = expr_new_function(expr_new_symbol("List"),
            (Expr*[]){poly_for_m, expr_new_integer(m)}, 2);
    }
    Expr* list = expr_new_function(expr_new_symbol("List"), out_args, max_mult);
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
        || lst->data.function.head->data.symbol != SYM_List
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
        || p->data.function.head->data.symbol != SYM_Times) {
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

    return expr_new_function(expr_new_symbol("List"),
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
        args[2] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){expr_new_symbol("Extension"), expr_copy((Expr*)alpha)}, 2);
    }
    Expr* call = expr_new_function(expr_new_symbol("Apart"), args, argc);
    Expr* ap = evaluate(call);
    expr_free(call);
    free(args);

    if (!ap) return NULL;

    if (ap->type == EXPR_FUNCTION
        && ap->data.function.head->type == EXPR_SYMBOL
        && ap->data.function.head->data.symbol == SYM_Plus) {
        size_t n = ap->data.function.arg_count;
        Expr** terms = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) terms[i] = expr_copy(ap->data.function.args[i]);
        expr_free(ap);
        return expr_new_function(expr_new_symbol("List"), terms, n);
    }
    return expr_new_function(expr_new_symbol("List"), (Expr*[]){ap}, 1);
}

/* ----- Cyclotomic / nth-root pre-detection ----- */

/* nthrootQ[p, x]: True iff p has the form `a x^n + b` with a, b free
 * of x.  Caught by Mathematica's MatchQ pattern at
 * IntegrateRational.m:1093.  The integer n is returned via *n_out.
 * Reserved for Phase 4's LogToReal closer (n-th-root-of-unity case);
 * unused in Phase 2 but kept here so the helper does not need to be
 * resurrected when Phase 4 lands. */
__attribute__((unused))
static bool intrat_nthroot_q(Expr* p, Expr* x, int* n_out, Expr** a_out, Expr** b_out) {
    /* Run Collect[p, x] to gather the polynomial in canonical form. */
    Expr* coll = expr_new_function(expr_new_symbol("Collect"),
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
        && c->data.function.head->data.symbol == SYM_Plus) {
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
            && term->data.function.head->data.symbol == SYM_Power
            && term->data.function.arg_count == 2
            && expr_eq(term->data.function.args[0], x)
            && term->data.function.args[1]->type == EXPR_INTEGER) {
            deg = (int)term->data.function.args[1]->data.integer;
            coeff = expr_new_integer(1);
        } else if (expr_eq(term, x)) {
            deg = 1; coeff = expr_new_integer(1);
        } else if (term->type == EXPR_FUNCTION
            && term->data.function.head->type == EXPR_SYMBOL
            && term->data.function.head->data.symbol == SYM_Times) {
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
                    && f->data.function.head->data.symbol == SYM_Power
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
                if (a) expr_free(a); if (b) expr_free(b); expr_free(c); return false;
            }
            coeff = (ccn == 0) ? expr_new_integer(1)
                  : (ccn == 1) ? ccoll[0]
                  : eval_and_free(internal_times(ccoll, ccn));
            free(ccoll);
        } else {
            if (a) expr_free(a); if (b) expr_free(b); expr_free(c); return false;
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
        || prs->data.function.head->data.symbol != SYM_List) return NULL;
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

/* IntRationalLogPart core. Returns either the {Q, S} pair list
 * (root_sum=false) or the symbolic RootSum sum (root_sum=true).
 * On algorithmic failure returns NULL.
 *
 * Direct port of IntegrateRational.m:751-801. */
static Expr* intrat_int_rational_log_part(Expr* f, Expr* x, Expr* t,
                                          bool root_sum) {
    intrat_trace("IntRationalLogPart", "IN", f);

    /* a = Numerator[f], d = Denominator[f]. */
    Expr* a_raw = intrat_numerator(f);
    Expr* a = eval_and_free(a_raw);
    Expr* d_raw = intrat_denominator(f);
    Expr* d = expr_expand(d_raw); expr_free(d_raw);

    /* prs = SubresultantPolynomialRemainders[d, a - t*D[d,x], x] */
    Expr* dprime = intrat_d(d, x);
    Expr* dprime_e = expr_expand(dprime); expr_free(dprime);
    Expr* tdprime = internal_times(
        (Expr*[]){expr_copy(t), dprime_e}, 2);
    Expr* atdp_raw = internal_subtract(
        (Expr*[]){expr_copy(a), tdprime}, 2);
    Expr* atdp = expr_expand(atdp_raw); expr_free(atdp_raw);

    Expr* prs_args[3] = { expr_copy(d), expr_copy(atdp), expr_copy(x) };
    Expr* prs_call = expr_new_function(
        expr_new_symbol("SubresultantPolynomialRemainders"), prs_args, 3);
    Expr* prs = evaluate(prs_call);
    expr_free(prs_call);
    if (!prs) {
        expr_free(a); expr_free(d); expr_free(atdp); return NULL;
    }

    /* resultant = primitive[Resultant[d, a - t*D[d,x], x], t] */
    Expr* res_args[3] = { expr_copy(d), atdp, expr_copy(x) };
    Expr* res_call = expr_new_function(
        expr_new_symbol("Resultant"), res_args, 3);
    Expr* resultant_raw = evaluate(res_call);
    expr_free(res_call);
    if (!resultant_raw) {
        expr_free(a); expr_free(d); expr_free(prs); return NULL;
    }
    Expr* resultant = intrat_primitive(resultant_raw, t);
    expr_free(resultant_raw);

    /* Q = SquareFree[resultant]: list of {poly_i, mult_i}. */
    Expr* Q = intrat_squarefree_list(resultant);
    expr_free(resultant);
    if (!Q) {
        expr_free(a); expr_free(d); expr_free(prs); return NULL;
    }

    int deg_d = get_degree_poly(d, x);
    size_t nQ = Q->data.function.arg_count;

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
             * IntegrateRational.m:765-768): use lc[d, x] in the
             * Diophantine. */
            Expr* lc_d = intrat_lc(d, x);
            Expr *r_eg = NULL, *_unused = NULL;
            if (!intrat_extended_euclidean(lc_d, Qi, expr_new_integer(1), t,
                                          &r_eg, &_unused)) {
                expr_free(lc_d); expr_free(Qi); continue;
            }
            expr_free(lc_d);
            if (_unused) expr_free(_unused);

            /* d := primitive[PolynomialRemainder[r_eg * d, Qi, t], x] */
            Expr* rd_raw = internal_times(
                (Expr*[]){r_eg, expr_copy(d)}, 2);
            Expr* rd = expr_expand(rd_raw); expr_free(rd_raw);
            Expr* rem = intrat_polyrem_t(rd, Qi, t);  /* defined just below */
            expr_free(rd);
            if (!rem) { expr_free(Qi); continue; }
            S[i_idx] = intrat_primitive(rem, x);
            expr_free(rem);
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
                Expr* gj = expr_expand(internal_power(
                    (Expr*[]){g, expr_new_integer((int64_t)jpow)}, 2));
                Expr* sn = intrat_exquo(s, gj, t);
                expr_free(gj);
                if (sn) { expr_free(s); s = sn; }
            }
            expr_free(A);

            /* r = ExtendedEuclidean[lc[s, x], Qi, 1, t][[1]]. */
            Expr* lc_s2 = intrat_lc(s, x);
            Expr *r_eg = NULL, *_unused = NULL;
            if (!intrat_extended_euclidean(lc_s2, Qi, expr_new_integer(1), t,
                                          &r_eg, &_unused)) {
                expr_free(lc_s2); expr_free(s); expr_free(Qi); continue;
            }
            expr_free(lc_s2);
            if (_unused) expr_free(_unused);

            /* s := primitive[PolynomialRemainder[r * s, Qi, t], x]. */
            Expr* rs_raw = internal_times(
                (Expr*[]){r_eg, s}, 2);
            Expr* rs = expr_expand(rs_raw); expr_free(rs_raw);
            Expr* rem = intrat_polyrem_t(rs, Qi, t);
            expr_free(rs);
            if (!rem) { expr_free(Qi); continue; }
            S[i_idx] = intrat_primitive(rem, x);
            expr_free(rem);
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
            Expr* Qi_pair = Q->data.function.args[i_idx];
            Expr* Qi = list_get(Qi_pair, 1);
            if (!Qi) continue;
            if (get_degree_poly(Qi, t) <= 0) { expr_free(Qi); continue; }

            /* RootSum[Function[t, Qi], Function[t, t Log[S[i]]]]. */
            Expr* func1 = expr_new_function(expr_new_symbol("Function"),
                (Expr*[]){expr_copy(t), Qi}, 2);
            Expr* logS = expr_new_function(expr_new_symbol("Log"),
                (Expr*[]){expr_copy(S[i_idx])}, 1);
            Expr* tlog = internal_times(
                (Expr*[]){expr_copy(t), logS}, 2);
            Expr* func2 = expr_new_function(expr_new_symbol("Function"),
                (Expr*[]){expr_copy(t), tlog}, 2);
            terms[nterms++] = expr_new_function(expr_new_symbol("RootSum"),
                (Expr*[]){func1, func2}, 2);
        }
        for (size_t i_idx = 0; i_idx < nQ; i_idx++) if (S[i_idx]) expr_free(S[i_idx]);
        free(S);
        expr_free(Q);
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
        Expr* Qi_pair = Q->data.function.args[i_idx];
        Expr* Qi = list_get(Qi_pair, 1);
        if (!Qi) continue;
        if (get_degree_poly(Qi, t) <= 0) { expr_free(Qi); continue; }
        pairs[npairs++] = expr_new_function(expr_new_symbol("List"),
            (Expr*[]){Qi, expr_copy(S[i_idx])}, 2);
    }
    for (size_t i_idx = 0; i_idx < nQ; i_idx++) if (S[i_idx]) expr_free(S[i_idx]);
    free(S);
    expr_free(Q);

    Expr* result = expr_new_function(expr_new_symbol("List"), pairs, npairs);
    free(pairs);
    intrat_trace("IntRationalLogPart", "OUT", result);
    return result;
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
        Expr* arctan = expr_new_function(expr_new_symbol("ArcTan"),
            (Expr*[]){canon}, 1);
        Expr* res = internal_times((Expr*[]){expr_new_integer(2), arctan}, 2);
        return eval_and_free(res);
    }

    /* A := Collect[a, x] (we use intrat_canonic for the closest
     * picocas equivalent — it folds rational coefficients to canonical
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
            Expr* arctan = expr_new_function(expr_new_symbol("ArcTan"),
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
        Expr* negB = expr_expand(internal_times(
            (Expr*[]){expr_new_integer(-1), expr_copy(B)}, 2));
        Expr* recurse = intrat_log_to_atan(negB, A, x);
        expr_free(negB); expr_free(A); expr_free(B);
        return recurse;
    }

    /* {g, {d, c}} = PolynomialExtendedGCD[B, -A, x]   so B*d - A*c = g. */
    Expr* negA = expr_expand(internal_times(
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
    Expr* arctan = expr_new_function(expr_new_symbol("ArcTan"),
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
    Expr* rule = expr_new_function(expr_new_symbol("Rule"),
        (Expr*[]){expr_copy(t), expr_copy(alpha)}, 2);
    Expr* substituted = expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){expr_copy(S), rule}, 2);
    Expr* Sat = evaluate(substituted);
    expr_free(substituted);
    Expr* logS = expr_new_function(expr_new_symbol("Log"),
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
        && term->data.function.head->data.symbol == SYM_Complex
        && term->data.function.arg_count == 2) {
        *re_out = expr_copy(term->data.function.args[0]);
        *im_out = expr_copy(term->data.function.args[1]);
        return;
    }
    if (term->type == EXPR_FUNCTION
        && term->data.function.head->type == EXPR_SYMBOL
        && term->data.function.head->data.symbol == SYM_Times) {
        /* Look for a Complex factor among the Times args. */
        Expr* complex_factor = NULL;
        size_t complex_idx = 0;
        size_t n = term->data.function.arg_count;
        for (size_t i = 0; i < n; i++) {
            Expr* arg = term->data.function.args[i];
            if (arg->type == EXPR_FUNCTION
                && arg->data.function.head->type == EXPR_SYMBOL
                && arg->data.function.head->data.symbol == SYM_Complex) {
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
        && expanded->data.function.head->data.symbol == SYM_Plus)) {
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
    Expr* rule = expr_new_function(expr_new_symbol("Rule"),
        (Expr*[]){expr_copy(t), expr_copy(sub_expr)}, 2);
    Expr* call = expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){expr_copy(e), rule}, 2);
    Expr* result = evaluate(call);
    expr_free(call);
    return result;
}

/* Try to dispatch a single quadratic factor `Q(t) = a t^2 + b t + c`
 * to its real-form contribution.  Returns NULL on failure (caller
 * falls back to symbolic RootSum). */
static Expr* logtoreal_quadratic(Expr* a, Expr* b, Expr* c,
                                 Expr* s, Expr* x, Expr* t) {
    /* discriminant = b^2 - 4 a c. */
    Expr* b2 = expr_expand(internal_power((Expr*[]){expr_copy(b), expr_new_integer(2)}, 2));
    Expr* fourac = expr_expand(
        internal_times((Expr*[]){expr_new_integer(4), expr_copy(a), expr_copy(c)}, 3));
    Expr* neg_fourac = internal_times(
        (Expr*[]){expr_new_integer(-1), fourac}, 2);
    Expr* disc = eval_and_free(internal_plus((Expr*[]){b2, neg_fourac}, 2));

    /* Decide the sign of disc.  We commit only when it's an integer or
     * rational; otherwise we punt to the symbolic form. */
    int disc_sign = 0;
    if (disc->type == EXPR_INTEGER) {
        disc_sign = (disc->data.integer > 0) ? 1 : (disc->data.integer < 0) ? -1 : 0;
    } else if (disc->type == EXPR_FUNCTION
        && disc->data.function.head->type == EXPR_SYMBOL
        && disc->data.function.head->data.symbol == SYM_Rational
        && disc->data.function.arg_count == 2
        && disc->data.function.args[0]->type == EXPR_INTEGER) {
        int64_t n = disc->data.function.args[0]->data.integer;
        disc_sign = (n > 0) ? 1 : (n < 0) ? -1 : 0;
    } else {
        expr_free(disc); return NULL;
    }

    Expr* two_a = eval_and_free(internal_times(
        (Expr*[]){expr_new_integer(2), expr_copy(a)}, 2));

    if (disc_sign >= 0) {
        /* Two real roots.  α± = (-b ± Sqrt[disc])/(2 a). */
        Expr* sqrt_d = expr_new_function(expr_new_symbol("Sqrt"),
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
        Expr* log_plus = expr_new_function(expr_new_symbol("Log"),
            (Expr*[]){s_at_plus}, 1);
        Expr* log_minus = expr_new_function(expr_new_symbol("Log"),
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
    Expr* neg_disc = eval_and_free(internal_times(
        (Expr*[]){expr_new_integer(-1), disc}, 2));
    Expr* sqrt_neg_disc = expr_new_function(expr_new_symbol("Sqrt"),
        (Expr*[]){neg_disc}, 1);
    sqrt_neg_disc = eval_and_free(sqrt_neg_disc);
    Expr* v_root = eval_and_free(internal_divide(
        (Expr*[]){sqrt_neg_disc, two_a}, 2));

    /* Substitute t -> u_root + I*v_root in s and split into A + I B. */
    Expr* iv = internal_times((Expr*[]){
        expr_new_function(expr_new_symbol("Complex"),
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
    Expr* A2 = expr_expand(internal_power(
        (Expr*[]){expr_copy(A), expr_new_integer(2)}, 2));
    Expr* B2 = expr_expand(internal_power(
        (Expr*[]){expr_copy(B), expr_new_integer(2)}, 2));
    Expr* mod2 = eval_and_free(internal_plus((Expr*[]){A2, B2}, 2));
    Expr* mod2_can = intrat_canonic(mod2); expr_free(mod2);
    Expr* logmod = expr_new_function(expr_new_symbol("Log"),
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

/* Walk a Factor[Q]'s output and dispatch every factor.  Returns the
 * accumulated real form, or NULL if any factor exceeds the bounded-
 * Solve scope (degree > 2 in t). */
static Expr* logtoreal_dispatch(Expr* factored, Expr* s, Expr* x, Expr* t) {
    Expr** factors;
    size_t nfactors;
    if (factored->type == EXPR_FUNCTION
        && factored->data.function.head->type == EXPR_SYMBOL
        && factored->data.function.head->data.symbol == SYM_Times) {
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
            && fac->data.function.head->data.symbol == SYM_Power
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
        } else {
            /* Higher-degree factor: out of bounded-Solve scope. */
            goto fail;
        }
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
        && factored->data.function.head->data.symbol == SYM_Times) {
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
            && fac->data.function.head->data.symbol == SYM_Power
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
        || pair_list->data.function.head->data.symbol != SYM_List) return NULL;
    size_t n = pair_list->data.function.arg_count;
    if (n == 0) return expr_new_integer(0);

    size_t tcap = 8, nterms = 0;
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * tcap);

    for (size_t i = 0; i < n; i++) {
        Expr* pair = pair_list->data.function.args[i];
        Expr* Qi = list_get(pair, 1);
        Expr* Si = list_get(pair, 2);
        if (!Qi || !Si) {
            for (size_t k = 0; k < nterms; k++) expr_free(terms[k]);
            free(terms);
            if (Qi) expr_free(Qi); if (Si) expr_free(Si);
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
 * Phase 6 — LogToArcTan / LogToArcTanh post-processing.
 * ====================================================================
 *
 * Combines `c Log[A] + c Log[B] -> c Log[A B]` and
 * `c Log[A] - c Log[B] -> c Log[A/B]` into single logarithms, then
 * recognises the `c Log[A] - c Log[B] -> 2 c ArcTanh[(B-A)/(B+A)]`
 * pattern when the ArcTanh argument simplifies to a rational
 * function in x.  Direct port of IntegrateRational.m:1722-1761
 * (LogToArcTanh) and :1896-1958 (LogToArcTan).
 *
 * Implemented as direct C transformations on Plus[...] of Log[...]
 * terms (per the plan note) rather than pattern-rewriting at the
 * picocas rule-engine layer.  The transformations are
 * differentiation-equivalent — they only beautify the output — so
 * the universal correctness check stays green either way.
 */

/* Decompose a term of a Plus head into (coeff, log_arg) when it
 * matches `coeff * Log[log_arg]`.  Returns true on success. */
static bool decompose_log_term(Expr* term, Expr** coeff_out, Expr** log_arg_out) {
    /* Plain Log[x]. */
    if (term->type == EXPR_FUNCTION
        && term->data.function.head->type == EXPR_SYMBOL
        && term->data.function.head->data.symbol == SYM_Log
        && term->data.function.arg_count == 1) {
        *coeff_out = expr_new_integer(1);
        *log_arg_out = expr_copy(term->data.function.args[0]);
        return true;
    }
    /* Times[..., Log[arg], ...].  Find the Log factor; the rest is
     * the coefficient. */
    if (term->type != EXPR_FUNCTION
        || term->data.function.head->type != EXPR_SYMBOL
        || term->data.function.head->data.symbol != SYM_Times) return false;

    size_t n = term->data.function.arg_count;
    Expr* log_arg = NULL;
    size_t log_idx = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* arg = term->data.function.args[i];
        if (arg->type == EXPR_FUNCTION
            && arg->data.function.head->type == EXPR_SYMBOL
            && arg->data.function.head->data.symbol == SYM_Log
            && arg->data.function.arg_count == 1) {
            if (log_arg) return false;  /* Multiple Log factors. */
            log_arg = arg->data.function.args[0];
            log_idx = i;
        }
    }
    if (!log_arg) return false;

    /* Build coeff = product of all other factors. */
    if (n == 1) {
        *coeff_out = expr_new_integer(1);
        *log_arg_out = expr_copy(log_arg);
        return true;
    }
    Expr** rest = (Expr**)malloc(sizeof(Expr*) * (n - 1));
    size_t k = 0;
    for (size_t i = 0; i < n; i++) {
        if (i != log_idx) rest[k++] = expr_copy(term->data.function.args[i]);
    }
    Expr* coeff;
    if (k == 1) coeff = rest[0];
    else coeff = eval_and_free(internal_times(rest, k));
    free(rest);
    *coeff_out = coeff;
    *log_arg_out = expr_copy(log_arg);
    return true;
}

/* zeroQ[e] = TrueQ[Cancel[Together[e]] === 0] for a single expression. */
static bool intrat_zero_q(Expr* e) {
    Expr* canon = intrat_canonic(e);
    bool ok = is_zero_poly(canon);
    expr_free(canon);
    return ok;
}

/* Simplify a Plus of Log terms by:
 *  1) merging `c Log[A] + c Log[B] -> c Log[Expand[A B]]`,
 *  2) merging `c Log[A] - c Log[B] -> c Log[A/B]` when
 *     Denominator[Cancel[Together[A/B]]] is FreeQ[..., x],
 *  3) rewriting `c Log[A] - c Log[B] -> 2 c ArcTanh[(B-A)/(B+A)]`
 *     when Denominator[Cancel[Together[(B-A)/(B+A)]]] is FreeQ[..., x].
 *
 * Always preserves D[result, x] up to canonic-zero. */
static Expr* intrat_log_to_arctanh(Expr* e, Expr* x) {
    if (!e) return NULL;
    if (!(e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Plus)) {
        return expr_copy(e);
    }

    size_t n = e->data.function.arg_count;
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * n);
    bool* used = (bool*)calloc(n, sizeof(bool));
    Expr** coeffs = (Expr**)calloc(n, sizeof(Expr*));
    Expr** logargs = (Expr**)calloc(n, sizeof(Expr*));
    bool* is_log = (bool*)calloc(n, sizeof(bool));

    for (size_t i = 0; i < n; i++) {
        terms[i] = e->data.function.args[i];
        is_log[i] = decompose_log_term(terms[i], &coeffs[i], &logargs[i]);
    }

    Expr** out = (Expr**)malloc(sizeof(Expr*) * (n * 2));
    size_t out_n = 0;

    for (size_t i = 0; i < n; i++) {
        if (used[i]) continue;
        if (!is_log[i]) {
            out[out_n++] = expr_copy(terms[i]);
            used[i] = true;
            continue;
        }
        /* Try to pair with another Log term j > i. */
        bool merged = false;
        for (size_t j = i + 1; j < n; j++) {
            if (used[j] || !is_log[j]) continue;
            if (!intrat_freeq_test(coeffs[i], x) || !intrat_freeq_test(coeffs[j], x)) continue;

            /* Same coefficient -> Log[A*B]. */
            Expr* delta = internal_subtract(
                (Expr*[]){expr_copy(coeffs[i]), expr_copy(coeffs[j])}, 2);
            if (intrat_zero_q(delta)) {
                expr_free(delta);
                Expr* prod = expr_expand(internal_times(
                    (Expr*[]){expr_copy(logargs[i]), expr_copy(logargs[j])}, 2));
                Expr* logp = expr_new_function(expr_new_symbol("Log"),
                    (Expr*[]){prod}, 1);
                Expr* combined = internal_times(
                    (Expr*[]){expr_copy(coeffs[i]), logp}, 2);
                out[out_n++] = eval_and_free(combined);
                used[i] = used[j] = true; merged = true; break;
            }
            expr_free(delta);

            /* Opposite coefficient: try ArcTanh first, then Log[A/B]. */
            Expr* sumcoef = internal_plus(
                (Expr*[]){expr_copy(coeffs[i]), expr_copy(coeffs[j])}, 2);
            if (intrat_zero_q(sumcoef)) {
                expr_free(sumcoef);
                /* Mathematica's log2ArcTanhRule emits
                 *   (c2 - c1) ArcTanh[(A + B) / (B - A)]
                 * for c1 Log[A] + c2 Log[B] with c1 + c2 == 0.
                 * Equivalent to 2 c2 ArcTanh[…].  We require the
                 * ArcTanh argument's Denominator to be free of x; if
                 * not, fall back to the simpler Log[A/B] rewrite. */
                Expr* sumAB_raw = internal_plus(
                    (Expr*[]){expr_copy(logargs[i]), expr_copy(logargs[j])}, 2);
                Expr* sumAB = expr_expand(sumAB_raw); expr_free(sumAB_raw);
                Expr* diffBA_raw = internal_subtract(
                    (Expr*[]){expr_copy(logargs[j]), expr_copy(logargs[i])}, 2);
                Expr* diffBA = expr_expand(diffBA_raw); expr_free(diffBA_raw);

                Expr* arg_raw = internal_divide(
                    (Expr*[]){expr_copy(sumAB), expr_copy(diffBA)}, 2);
                Expr* arg_can = intrat_canonic(arg_raw); expr_free(arg_raw);
                Expr* arg_den = intrat_denominator(arg_can);
                Expr* arg_den_eval = eval_and_free(arg_den);
                bool atanh_ok = intrat_freeq_test(arg_den_eval, x)
                              && !is_zero_poly(arg_can);
                expr_free(arg_den_eval);
                if (atanh_ok) {
                    /* (c2 - c1) ArcTanh[(A + B) / (B - A)] */
                    Expr* coef_diff = internal_subtract(
                        (Expr*[]){expr_copy(coeffs[j]), expr_copy(coeffs[i])}, 2);
                    coef_diff = eval_and_free(coef_diff);
                    Expr* atanh = expr_new_function(expr_new_symbol("ArcTanh"),
                        (Expr*[]){arg_can}, 1);
                    Expr* term = internal_times(
                        (Expr*[]){coef_diff, atanh}, 2);
                    out[out_n++] = eval_and_free(term);
                    expr_free(sumAB); expr_free(diffBA);
                    used[i] = used[j] = true; merged = true; break;
                }
                expr_free(arg_can);

                /* Log[A/B] with denominator free of x. */
                Expr* AoverB_raw = internal_divide(
                    (Expr*[]){expr_copy(logargs[i]), expr_copy(logargs[j])}, 2);
                Expr* AoverB = intrat_canonic(AoverB_raw); expr_free(AoverB_raw);
                Expr* AoverB_den = intrat_denominator(AoverB);
                Expr* AoverB_den_eval = eval_and_free(AoverB_den);
                bool divlog_ok = intrat_freeq_test(AoverB_den_eval, x);
                expr_free(AoverB_den_eval);
                if (divlog_ok) {
                    Expr* logp = expr_new_function(expr_new_symbol("Log"),
                        (Expr*[]){AoverB}, 1);
                    Expr* term = internal_times(
                        (Expr*[]){expr_copy(coeffs[i]), logp}, 2);
                    out[out_n++] = eval_and_free(term);
                    expr_free(sumAB); expr_free(diffBA);
                    used[i] = used[j] = true; merged = true; break;
                }
                expr_free(AoverB);
                expr_free(sumAB); expr_free(diffBA);
                continue;
            }
            expr_free(sumcoef);
        }
        if (!merged) {
            out[out_n++] = expr_copy(terms[i]);
            used[i] = true;
        }
    }

    for (size_t i = 0; i < n; i++) {
        if (coeffs[i])  expr_free(coeffs[i]);
        if (logargs[i]) expr_free(logargs[i]);
    }
    free(coeffs); free(logargs); free(is_log); free(used); free(terms);

    Expr* result;
    if (out_n == 0) { free(out); result = expr_new_integer(0); }
    else if (out_n == 1) { result = out[0]; free(out); }
    else { result = eval_and_free(internal_plus(out, out_n)); free(out); }
    return result;
}

/* ====================================================================
 * Phase 7 — Top-level options and output cleanup.
 * ====================================================================
 *
 * Adds option parsing to `Integrate`IntegrateRational[f, x, ...]` and
 * the final output-canonicalisation pass described in the
 * "final-pass output cleanup" section of INTEGRATE_PLAN.md:
 *
 *   1. ArcTan / ArcTanh sign normalisation — `ArcTan[-arg] -> -ArcTan[arg]`
 *      (and similarly for ArcTanh) so the printed form is independent of
 *      sign conventions deep inside the resultant chain.
 *   2. Plus-head sort happens automatically via picocas's ATTR_ORDERLESS,
 *      which fires when we eval_and_free the assembled sum.
 *
 * Options:
 *   "PFD"        -> True    Apart-split per-summand loop
 *   "LogToArcTan" -> True   Phase 6 Log -> ArcTanh post-processing
 *   Extension    -> Automatic   forwarded to Apart / Factor when set
 *
 * The defaults reproduce Mathematica's IntegrateRational.m behaviour.
 */

/* ArcTan / ArcTanh sign normalisation: pull a leading minus out of
 * the argument so the printed form is canonical.  Walks the input
 * expression top-down, rewriting just the relevant heads. */
static Expr* normalize_inverse_trig_signs(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy(e);

    /* ArcTan[-(...)] -> -ArcTan[...], same for ArcTanh.  We look for
     * an argument that starts with a unary minus or a Times whose
     * first factor is a literal -1. */
    if (e->data.function.head->type == EXPR_SYMBOL
        && (e->data.function.head->data.symbol == SYM_ArcTan
            || e->data.function.head->data.symbol == SYM_ArcTanh)
        && e->data.function.arg_count == 1) {
        Expr* arg = e->data.function.args[0];
        bool negative = false;
        Expr* stripped = NULL;
        if (arg->type == EXPR_INTEGER && arg->data.integer < 0) {
            negative = true;
            stripped = expr_new_integer(-arg->data.integer);
        } else if (arg->type == EXPR_FUNCTION
            && arg->data.function.head->type == EXPR_SYMBOL
            && arg->data.function.head->data.symbol == SYM_Times
            && arg->data.function.arg_count >= 1
            && arg->data.function.args[0]->type == EXPR_INTEGER
            && arg->data.function.args[0]->data.integer == -1) {
            negative = true;
            size_t n = arg->data.function.arg_count;
            if (n == 2) stripped = expr_copy(arg->data.function.args[1]);
            else {
                Expr** rest = (Expr**)malloc(sizeof(Expr*) * (n - 1));
                for (size_t i = 1; i < n; i++) rest[i - 1] = expr_copy(arg->data.function.args[i]);
                stripped = eval_and_free(internal_times(rest, n - 1));
                free(rest);
            }
        }
        if (negative && stripped) {
            Expr* inner = expr_new_function(
                expr_copy(e->data.function.head),
                (Expr*[]){stripped}, 1);
            Expr* neg = internal_times(
                (Expr*[]){expr_new_integer(-1), inner}, 2);
            return eval_and_free(neg);
        }
    }

    /* Recurse: rebuild the function with normalised children. */
    size_t n = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        new_args[i] = normalize_inverse_trig_signs(e->data.function.args[i]);
    }
    Expr* head = expr_copy(e->data.function.head);
    Expr* result = expr_new_function(head, new_args, n);
    /* Run a single evaluator pass so Plus / Times canonicalise. */
    return eval_and_free(result);
}

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
 * than emit a symbolic RootSum that picocas can't yet differentiate
 * back to the integrand).
 */

/* Try to close `r/d` via the LRT-then-LogToReal pipeline.  Returns
 * NULL on failure. */
static Expr* try_lrt_close(Expr* num, Expr* den, Expr* x) {
    if (is_zero_poly(num)) return expr_new_integer(0);
    Expr* h = internal_divide(
        (Expr*[]){expr_copy(num), expr_copy(den)}, 2);
    Expr* h_eval = eval_and_free(h);

    Expr* t_sym = expr_new_symbol("Integrate`Private`tt$");
    Expr* lrt_pairs = intrat_int_rational_log_part(h_eval, x, t_sym, false);
    expr_free(h_eval);
    if (!lrt_pairs) { expr_free(t_sym); return NULL; }

    size_t npairs = lrt_pairs->data.function.arg_count;
    Expr** lr_terms = (Expr**)malloc(sizeof(Expr*) * (npairs ? npairs : 1));
    size_t lr_count = 0;
    bool all_ok = true;
    for (size_t i = 0; i < npairs; i++) {
        Expr* pair = lrt_pairs->data.function.args[i];
        Expr* Qi = list_get(pair, 1);
        Expr* Si = list_get(pair, 2);
        if (!Qi || !Si) { all_ok = false; if (Qi) expr_free(Qi); if (Si) expr_free(Si); break; }
        Expr* term = intrat_log_to_real(Qi, Si, x, t_sym);
        expr_free(Qi); expr_free(Si);
        if (!term) { all_ok = false; break; }
        lr_terms[lr_count++] = term;
    }
    if (!all_ok) {
        for (size_t k = 0; k < lr_count; k++) expr_free(lr_terms[k]);
        free(lr_terms);
        /* Linear-Q fallback: every Q_i factors completely over Q. */
        Expr* closed = intrat_linear_q_closer(lrt_pairs, x, t_sym);
        expr_free(lrt_pairs); expr_free(t_sym);
        return closed;
    }
    expr_free(lrt_pairs); expr_free(t_sym);
    if (lr_count == 0) { free(lr_terms); return expr_new_integer(0); }
    if (lr_count == 1) { Expr* r = lr_terms[0]; free(lr_terms); return r; }
    Expr* sum = internal_plus(lr_terms, lr_count);
    free(lr_terms);
    return eval_and_free(sum);
}

/* Process h piece-by-piece via Apart, ExtractConstants, and the
 * derivative-recognition / LRT pipeline.  Returns the closed sum or
 * NULL if any piece can't be closed. */
static Expr* intrat_integrate_summands(Expr* h, Expr* x) {
    Expr* pieces = intrat_apart_list(h, x, NULL);
    if (!pieces || pieces->type != EXPR_FUNCTION
        || pieces->data.function.head->type != EXPR_SYMBOL
        || pieces->data.function.head->data.symbol != SYM_List) {
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
/* IntegrateRational top-level (Phase 1 skeleton + Phase 2 hookup).    */
/* ------------------------------------------------------------------ */

static Expr* intrat_integrate_rational(Expr* f, Expr* x) {
    intrat_trace("IntegrateRational", "IN", f);

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
        intrat_trace("IntegrateRational", "OUT", result);
        return result;
    }

    /* Try the derivative-recognition fast path on R/den. */
    Expr* fast = intrat_derivative_recognition(R, den_eval, x);
    if (fast) {
        Expr* sum = internal_plus(
            (Expr*[]){poly_int, fast}, 2);
        Expr* result = eval_and_free(sum);
        expr_free(R); expr_free(den_eval);
        intrat_trace("IntegrateRational", "OUT", result);
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
        intrat_trace("IntegrateRational", "OUT", result);
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
        /* Phase 6 final pass: combine c Log[A] ± c Log[B] terms.
         * Phase 7 final pass: normalise ArcTan/ArcTanh argument signs.*/
        Expr* simplified = intrat_log_to_arctanh(sum_eval, x);
        expr_free(sum_eval);
        Expr* result = normalize_inverse_trig_signs(simplified);
        expr_free(simplified);
        intrat_trace("IntegrateRational", "OUT", result);
        return result;
    }

    /* Beyond Phase 5's combined scope — leave unevaluated. */
    expr_free(h); expr_free(g); expr_free(poly_int);
    intrat_trace("IntegrateRational", "OUT (unresolved)", f);
    return NULL;
}

/* ------------------------------------------------------------------ */
/* BuiltinFunc wrappers.  Each follows the picocas convention: caller  */
/* owns res; return NULL on failure.                                   */
/* ------------------------------------------------------------------ */

/* Recognise a trailing option (Rule[lhs, rhs] / RuleDelayed[...]).
 * Returns true if `opt` is a Rule with one of the IntegrateRational
 * option names on the lhs. */
static bool is_intrat_option(Expr* opt) {
    if (!opt || opt->type != EXPR_FUNCTION
        || opt->data.function.head->type != EXPR_SYMBOL
        || (opt->data.function.head->data.symbol != SYM_Rule
            && opt->data.function.head->data.symbol != SYM_RuleDelayed)
        || opt->data.function.arg_count != 2) return false;
    Expr* lhs = opt->data.function.args[0];
    if (lhs->type == EXPR_SYMBOL) {
        return lhs->data.symbol == SYM_Extension;
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
            && (last->data.function.head->data.symbol == SYM_Rule
                || last->data.function.head->data.symbol == SYM_RuleDelayed)
            && last->data.function.arg_count == 2
            && last->data.function.args[0]->type == EXPR_SYMBOL
            && last->data.function.args[0]->data.symbol == SYM_RootSum) {
            Expr* val = last->data.function.args[1];
            root_sum = (val->type == EXPR_SYMBOL && val->data.symbol == SYM_True);
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

Expr* builtin_intrat_logtoarctanh(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* e = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    return intrat_log_to_arctanh(e, x);
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
    install("Integrate`IntegrateRational",
            builtin_intrat_integraterational,
            "Integrate`IntegrateRational[f, x] is the explicit form of the\n"
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

    /* Trace flag: Integrate`$Verbose. Default False. */
    {
        Expr* pat = expr_new_symbol("Integrate`$Verbose");
        Expr* val = expr_new_symbol("False");
        symtab_add_own_value("Integrate`$Verbose", pat, val);
        expr_free(pat); expr_free(val);
    }

    /* RootSum is registered in src/root.c (Phase 8b-prereq) — the
     * canonical home of held algebraic-root machinery.  Nothing to
     * do here. */
}
