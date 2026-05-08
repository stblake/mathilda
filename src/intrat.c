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

/* ------------------------------------------------------------------ */
/* IntegrateRational top-level (Phase 1 skeleton).                     */
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

    /* Otherwise try derivative recognition on h. */
    Expr* h_num = intrat_numerator(h);
    Expr* h_den = intrat_denominator(h);
    h_num = eval_and_free(h_num);
    h_den = eval_and_free(h_den);
    expr_free(h);

    Expr* h_int = intrat_derivative_recognition(h_num, h_den, x);
    expr_free(h_num); expr_free(h_den);
    if (!h_int) {
        /* Phase 1 cannot close this — leave the call unevaluated.    */
        expr_free(g); expr_free(poly_int);
        intrat_trace("IntegrateRational", "OUT (unresolved)", f);
        return NULL;
    }

    Expr* sum = internal_plus(
        (Expr*[]){poly_int, g, h_int}, 3);
    Expr* result = eval_and_free(sum);
    intrat_trace("IntegrateRational", "OUT", result);
    return result;
}

/* ------------------------------------------------------------------ */
/* BuiltinFunc wrappers.  Each follows the picocas convention: caller  */
/* owns res; return NULL on failure.                                   */
/* ------------------------------------------------------------------ */

Expr* builtin_intrat_integraterational(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    return intrat_integrate_rational(f, x);
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

    /* Trace flag: Integrate`$Verbose. Default False. */
    {
        Expr* pat = expr_new_symbol("Integrate`$Verbose");
        Expr* val = expr_new_symbol("False");
        symtab_add_own_value("Integrate`$Verbose", pat, val);
        expr_free(pat); expr_free(val);
    }

    /* RootSum: passive head, registered with no builtin so it appears
     * verbatim in output (Phase 2 will upgrade this). */
    symtab_get_def("RootSum")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("RootSum",
        "RootSum[fn, summand] is a symbolic placeholder used by the rational\n"
        "integrator to denote a sum over the roots of fn. Phase 2 upgrades\n"
        "this to a fully-evaluating head once LogToReal is implemented.");
}
