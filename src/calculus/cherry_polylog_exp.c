/* cherry_polylog_exp.c — exponential-tower polylogarithm-ladder integration.
 *
 * General-weight, algebraic-root generalisation of the exponential-tower
 * dilogarithm engine (cherry_dilog_exp.c).  For an integrand P(x)/Q(E^(c x))
 * — rational in the single exponential kernel theta = E^(c x), polynomial in x —
 * partial-fraction over the roots rho of Q (rational OR algebraic) reduces to
 * simple poles x^n/(theta - rho), each closed by the EXACT Cherry polylogarithm
 * ladder built from d/dx PolyLog[k, rho/theta] = -c PolyLog[k-1, rho/theta]:
 *
 *   INT x^n/(theta - rho) dx
 *     = Sum_{k=0}^{n} -(1/rho) (n!/(n-k)!) / c^(k+1) x^(n-k) PolyLog[k+1, rho/theta]
 *
 * with PolyLog[1, z] = -Log[1-z].  Weight rises to n+1, so this is the ladder
 * above the weight-2 dilogarithm.  The construction is exact (a telescoping
 * integration by parts), re-checked by a PowerExpand diff-back over the (possibly
 * algebraic) constant field; a mis-generation declines (NULL), never a wrong form.
 *
 * Scope of THIS increment: a single commensurable exponential kernel; Q a proper
 * fraction in theta with SIMPLE nonzero roots (rational or algebraic).  Declines
 * cleanly: an improper fraction / polynomial-in-theta part, a theta=0 (Laurent)
 * pole, repeated roots (higher-order-pole ladders), and outer-log integrands
 * (Log[P(theta)], left to cherry_dilog_exp).  Registered ahead of
 * rt_cherry_dilog_exp so it supplies the general (and cleaner) rational forms.
 */

#include "cherry_polylog_exp.h"
#include "cherry_dilog_exp.h"     /* shares the top-level depth-gate rationale */
#include "risch_util.h"
#include "risch_singleext.h"      /* rt_exp_kernelize */
#include "integrate.h"            /* g_integrate_depth (top-level-only gating) */

#include "expr.h"
#include "eval.h"
#include "print.h"
#include "sym_intern.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static Expr* mk_sym(const char* s)  { return expr_new_symbol(s); }
static Expr* mk_int(long n)         { return expr_new_integer(n); }
static Expr* mk_pow(Expr* b, Expr* e) {
    return expr_new_function(mk_sym("Power"), (Expr*[]){ b, e }, 2);
}
static Expr* mk_neg(Expr* a) {
    return expr_new_function(mk_sym("Times"), (Expr*[]){ mk_int(-1), a }, 2);
}
static Expr* mk_plus2(Expr* a, Expr* b) {
    return expr_new_function(mk_sym("Plus"), (Expr*[]){ a, b }, 2);
}
static Expr* mk_times2(Expr* a, Expr* b) {
    return expr_new_function(mk_sym("Times"), (Expr*[]){ a, b }, 2);
}
static Expr* mk_log(Expr* a) { return expr_new_function(mk_sym("Log"), (Expr*[]){ a }, 1); }
static Expr* mk_fn2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(mk_sym(h), (Expr*[]){ a, b }, 2);
}

/* Falling factorial n!/(n-k)! as an owned integer Expr (bignum-safe: a product of
 * integers evaluated by Mathilda, so no int64 overflow at large n). */
static Expr* falling_factorial(long n, long k) {
    if (k <= 0) return mk_int(1);
    Expr** fac = malloc((size_t)k * sizeof(Expr*));
    for (long i = 0; i < k; i++) fac[i] = mk_int(n - i);   /* n (n-1) ... (n-k+1) */
    Expr* prod = expr_new_function(mk_sym("Times"), fac, (size_t)k);
    free(fac);
    return rt_eval_own(prod);
}

Expr* rt_cherry_polylog_exp(Expr* f, Expr* x) {
    /* Top-level only: a DerivativeDivides u=Log[x] substitution can turn a
     * LOG-tower integrand into this exp-tower shape; handling it in nested
     * recursion would pre-empt the direct log-tower engines.  Same gate as
     * cherry_dilog_exp. */
    if (g_integrate_depth > 1) return NULL;
    if (!rt_find_exp_of_x(f, x)) return NULL;

    /* 1. kernelize E^(c x) -> theta = rmT; derive the rate c from the primitive. */
    Expr* u = NULL;
    Expr* F = rt_exp_kernelize(f, x, &u);
    if (!F || !u) { if (F) expr_free(F); if (u) expr_free(u); return NULL; }
    Expr* th = mk_sym("rmT");
    if (!rt_is_poly(u, x) || rt_degree(u, x) != 1) {
        expr_free(F); expr_free(u); expr_free(th); return NULL;
    }
    Expr* u0 = rt_coeff(u, x, 0);
    bool linok = u0 && rt_is_zero(u0);
    if (u0) expr_free(u0);
    Expr* c = rt_coeff(u, x, 1);
    expr_free(u);
    if (!linok || !c || rt_is_zero(c)) {
        if (c) expr_free(c);
        expr_free(F); expr_free(th);
        return NULL;
    }

    /* 2. F = num(x,theta)/den(theta), a PROPER fraction in theta, den free of x,
     * no residual exp/log of x (outer-log integrands go to cherry_dilog_exp). */
    Expr* Ft  = rt_eval1("Together", expr_copy(F));
    Expr* num = Ft ? rt_eval1("Numerator", expr_copy(Ft)) : NULL;
    Expr* den = Ft ? rt_eval1("Denominator", expr_copy(Ft)) : NULL;
    if (Ft) expr_free(Ft);
    bool ok = num && den
        && rt_free_of_x(F, th) == false        /* actually involves theta */
        && rt_free_of_head(F, "Log")
        && rt_find_exp_of_x(F, x) == NULL
        && rt_is_poly(den, th) && rt_free_of_x(den, x)
        && rt_is_poly(num, th)
        && rt_degree(num, th) < rt_degree(den, th)   /* proper */
        && rt_degree(den, th) >= 1;
    if (!ok) {
        if (num) expr_free(num);
        if (den) expr_free(den);
        expr_free(F); expr_free(th); expr_free(c);
        return NULL;
    }
    long ddeg = rt_degree(den, th);

    /* 3. simple nonzero roots rho of den (rational or algebraic) via Solve. */
    Expr* sols = rt_eval2("Solve", mk_fn2("Equal", expr_copy(den), mk_int(0)),
                          expr_copy(th));
    Expr* roots[16]; size_t nroot = 0;
    bool bad = false;
    if (sols && sols->type == EXPR_FUNCTION && rt_head_is(sols, "List")) {
        for (size_t i = 0; i < sols->data.function.arg_count; i++) {
            Expr* rule = sols->data.function.args[i];
            if (rule->type != EXPR_FUNCTION || !rt_head_is(rule, "List")
                || rule->data.function.arg_count != 1) { bad = true; break; }
            Expr* rr = rule->data.function.args[0];
            if (rr->type != EXPR_FUNCTION || !rt_head_is(rr, "Rule")) { bad = true; break; }
            Expr* val = rr->data.function.args[1];
            if (rt_is_zero(val)) { bad = true; break; }        /* theta=0 Laurent pole */
            bool dup = false;
            for (size_t j = 0; j < nroot; j++) if (expr_eq(roots[j], val)) { dup = true; break; }
            if (dup) { bad = true; break; }                    /* repeated root */
            if (nroot < 16) roots[nroot++] = expr_copy(val); else bad = true;
        }
    } else bad = true;
    if (sols) expr_free(sols);
    /* every root must be simple: distinct-root count == denominator degree. */
    if (bad || (long)nroot != ddeg) {
        for (size_t i = 0; i < nroot; i++) expr_free(roots[i]);
        expr_free(num); expr_free(den); expr_free(F); expr_free(th); expr_free(c);
        return NULL;
    }

    /* 4. per-pole polylog ladder.  For each rho: residue A(x) = Cancel[num(th-rho)
     * /den] /. th->rho (poly in x), then add the ladder terms. */
    Expr** terms = NULL; size_t nterms = 0, tcap = 0;
    bool fail = false;
    for (size_t r = 0; r < nroot && !fail; r++) {
        Expr* rho = roots[r];
        /* A = Expand[ Cancel[ num (theta - rho) / den ] /. theta -> rho ] */
        Expr* piece = mk_times2(mk_times2(expr_copy(num),
                          mk_plus2(expr_copy(th), mk_neg(expr_copy(rho)))),
                          mk_pow(expr_copy(den), mk_int(-1)));
        Expr* canc = rt_eval1("Cancel", piece);
        Expr* Asub = canc ? rt_eval_own(mk_fn2("ReplaceAll", canc,
                          mk_fn2("Rule", expr_copy(th), expr_copy(rho)))) : NULL;
        Expr* A = Asub ? rt_eval1("Expand", Asub) : NULL;
        if (!A || !rt_is_poly(A, x)) { if (A) expr_free(A); fail = true; break; }
        long n = rt_degree(A, x);
        if (n < 0) { expr_free(A); continue; }           /* zero residue */
        /* rho/theta, reused per PolyLog argument (kept as rmT here). */
        for (long j = 0; j <= n && !fail; j++) {
            Expr* aj = rt_coeff(A, x, j);
            if (!aj || rt_is_zero(aj)) { if (aj) expr_free(aj); continue; }
            for (long k = 0; k <= j; k++) {
                /* coeff = Simplify[ -(aj/rho) (j!/(j-k)!) / c^(k+1) ] */
                Expr* ff = falling_factorial(j, k);
                Expr* coeff = rt_eval1("Simplify", mk_times2(
                    mk_times2(mk_int(-1), mk_times2(expr_copy(aj),
                        mk_pow(expr_copy(rho), mk_int(-1)))),
                    mk_times2(ff, mk_pow(expr_copy(c), mk_int(-(k + 1))))));
                Expr* arg = mk_times2(expr_copy(rho), mk_pow(expr_copy(th), mk_int(-1)));
                Expr* plog = expr_new_function(mk_sym("PolyLog"),
                    (Expr*[]){ mk_int(k + 1), arg }, 2);
                Expr* term = mk_times2(coeff,
                    mk_times2(mk_pow(expr_copy(x), mk_int(j - k)), plog));
                if (nterms == tcap) {
                    tcap = tcap ? tcap * 2 : 16;
                    terms = realloc(terms, tcap * sizeof(Expr*));
                }
                terms[nterms++] = term;
            }
            expr_free(aj);
        }
        expr_free(A);
    }

    /* 5. assemble, back-substitute theta -> E^(c x), reduce PolyLog[1] -> -Log. */
    Expr* result = NULL;
    if (!fail && nterms > 0) {
        Expr** tc = malloc(nterms * sizeof(Expr*));
        for (size_t i = 0; i < nterms; i++) tc[i] = expr_copy(terms[i]);
        Expr* Q = expr_new_function(mk_sym("Plus"), tc, nterms);
        free(tc);
        Q = rt_eval_own(mk_fn2("ReplaceAll", Q,
                mk_fn2("Rule", expr_copy(th),
                       mk_pow(mk_sym("E"), mk_times2(expr_copy(c), expr_copy(x))))));
        /* PolyLog[1, z_] :> -Log[1 - z] (branch-safe weight-1 reduction). */
        if (Q) {
            Expr* zz = mk_sym("chpl$z");
            Expr* p1rule = expr_new_function(mk_sym("RuleDelayed"),
                (Expr*[]){ expr_new_function(mk_sym("PolyLog"),
                              (Expr*[]){ mk_int(1),
                                  expr_new_function(mk_sym("Pattern"),
                                     (Expr*[]){ expr_copy(zz),
                                         expr_new_function(mk_sym("Blank"), NULL, 0) }, 2) }, 2),
                           mk_neg(mk_log(mk_plus2(mk_int(1), mk_neg(expr_copy(zz))))) }, 2);
            expr_free(zz);
            Q = rt_eval_own(mk_fn2("ReplaceAll", Q, p1rule));
        }
        /* 6. PowerExpand diff-back verify (exact over the algebraic field). */
        if (Q && rt_free_of_head(Q, "Integrate")) {
            Expr* diff = mk_plus2(rt_eval2("D", expr_copy(Q), expr_copy(x)),
                                  mk_neg(expr_copy(f)));
            Expr* aa = mk_sym("chpl$a");
            Expr* lfr = expr_new_function(mk_sym("RuleDelayed"),
                (Expr*[]){ mk_log(expr_new_function(mk_sym("Pattern"),
                              (Expr*[]){ expr_copy(aa),
                                         expr_new_function(mk_sym("Blank"), NULL, 0) }, 2)),
                           mk_log(expr_new_function(mk_sym("Factor"),
                              (Expr*[]){ expr_new_function(mk_sym("Together"),
                                           (Expr*[]){ expr_copy(aa) }, 1) }, 1)) }, 2);
            expr_free(aa);
            Expr* diff2 = rt_eval_own(mk_fn2("ReplaceAll", diff, lfr));
            Expr* chk = diff2 ? rt_eval1("Simplify", rt_eval1("PowerExpand", diff2)) : NULL;
            if (chk && chk->type == EXPR_INTEGER && chk->data.integer == 0) result = Q;
            else expr_free(Q);
            if (chk) expr_free(chk);
        } else if (Q) expr_free(Q);
    }

    for (size_t i = 0; i < nterms; i++) expr_free(terms[i]);
    free(terms);
    for (size_t i = 0; i < nroot; i++) expr_free(roots[i]);
    expr_free(num); expr_free(den); expr_free(F); expr_free(th); expr_free(c);
    return result;
}

/* Integrate`Cherry`PolyLogExp[f, x] — direct debuggable surface, registered by
 * cherry_builtins_init (cherry_driver.c). */
Expr* builtin_cherry_polylog_exp(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    return rt_cherry_polylog_exp(res->data.function.args[0],
                                 res->data.function.args[1]);
}
