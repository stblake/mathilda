/*
 * reduce_algfiber.c
 *
 * Real-algebraic-coefficient fibre isolation for `Reduce`'s CAD (Phase 6b).
 * See reduce_algfiber.h for the method and the soundness contract.
 *
 * All polynomial work is delegated to the evaluator (ReplaceAll / Resultant /
 * Exponent) and to the exact primitives it already ships: rru_collect_roots for
 * integer-coefficient real-root isolation and flint_qqbar_equal for the exact
 * zero test that discards the conjugate-spurious roots the resultant introduces.
 */
#include "reduce_algfiber.h"
#include "reduce_real_util.h"

#include "eval.h"
#include "sym_names.h"
#include "flint_qqbar.h"

#include <stdlib.h>

/* ------------------------------------------------------------------ *
 *  Resource ceiling (a sound decline, never a wrong answer)           *
 * ------------------------------------------------------------------ *
 * Iterated resultants multiply the var-degree per algebraic tower level, so a
 * deep/high-degree tower can blow up.  A budget overrun bails (Reduce stays
 * unevaluated) rather than running unbounded -- matching the rest of the engine's
 * budget-exhausted declines.  The ceilings are generous: the sphere/ball towers
 * that motivate 6b stay far below them. */
#define ALGFIBER_MAX_VAR_DEGREE 128
#define ALGFIBER_MAX_NODES      50000

/* ------------------------------------------------------------------ *
 *  Small node builders (each CONSUMES its Expr* arguments)            *
 * ------------------------------------------------------------------ */

static Expr* mkfun2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2);
}
static Expr* mkfun3(const char* h, Expr* a, Expr* b, Expr* c) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b, c }, 3);
}

static bool is_zero_expr(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 0;
}

static size_t node_count(const Expr* e, size_t budget) {
    if (!e || budget == 0) return 0;
    size_t c = 1;
    if (e->type == EXPR_FUNCTION) {
        c += node_count(e->data.function.head, budget);
        for (size_t i = 0; i < e->data.function.arg_count && c < budget; i++)
            c += node_count(e->data.function.args[i], budget - c);
    }
    return c;
}

/* Exponent[p, v] as an int; -1 when the degree is not a plain integer. */
static int degree_in(const Expr* p, const Expr* v) {
    Expr* e = eval_and_free(mkfun2(SYM_Exponent, expr_copy((Expr*)p), expr_copy((Expr*)v)));
    int d = (e->type == EXPR_INTEGER) ? (int)e->data.integer : -1;
    expr_free(e);
    return d;
}

/* ReplaceAll[p, v -> s], CONSUMING the owned `p`. */
static Expr* repl(Expr* p, const Expr* v, const Expr* s) {
    Expr* rule = mkfun2(SYM_Rule, expr_copy((Expr*)v), expr_copy((Expr*)s));
    return eval_and_free(mkfun2(SYM_ReplaceAll, p, rule));
}

/* factor with vv[0..nlev-1] -> vals[0..nlev-1] and var -> pt, evaluated to a
 * constant.  Owns the result. */
static Expr* subst_point(const Expr* factor, Expr** vv, Expr** vals, int nlev,
                         const Expr* var, const Expr* pt) {
    Expr* e = expr_copy((Expr*)factor);
    for (int i = 0; i < nlev; i++) e = repl(e, vv[i], vals[i]);
    e = repl(e, var, pt);
    return e;
}

static void prov_push(Expr*** arr, int* n, int* cap, int** prov, int factor_id, Expr* v) {
    if (*n == *cap) {
        *cap = *cap ? *cap * 2 : 8;
        *arr = realloc(*arr, (size_t)*cap * sizeof(Expr*));
        if (prov) *prov = realloc(*prov, (size_t)*cap * sizeof(int));
    }
    if (prov) (*prov)[*n] = factor_id;
    (*arr)[(*n)++] = v;
}

bool rru_algebraic_fiber_roots(const Expr* factor, const Expr* var,
                               Expr** vv, Expr** vals, Expr** defs, int nlev,
                               Expr*** arr, int* n, int* cap,
                               int** prov, int factor_id,
                               const ReduceOpts* opts) {
#ifndef USE_FLINT
    (void)factor; (void)var; (void)vv; (void)vals; (void)defs; (void)nlev;
    (void)arr; (void)n; (void)cap; (void)prov; (void)factor_id; (void)opts;
    return false;                         /* no qqbar oracle -> decline */
#else
    /* 1. Project the fibre back down to Q: substitute the rational levels
     *    directly, then eliminate each algebraic tower variable by resultant. */
    Expr* R = expr_copy((Expr*)factor);
    for (int i = 0; i < nlev; i++)
        if (!defs[i]) R = repl(R, vv[i], vals[i]);      /* rational coordinate */

    for (int i = nlev - 1; i >= 0; i--) {
        if (!defs[i]) continue;                          /* rational: already done */
        /* defs[i] is a factor over Q in vv[0..i]; substitute its own rational
         * lower levels, keep the remaining algebraic tower variables symbolic. */
        Expr* D = expr_copy(defs[i]);
        for (int t = 0; t < i; t++)
            if (!defs[t]) D = repl(D, vv[t], vals[t]);
        R = eval_and_free(mkfun3(SYM_Resultant, R, D, expr_copy(vv[i])));
        if (is_zero_expr(R)) { expr_free(R); return false; }        /* degenerate */
        if (degree_in(R, var) > ALGFIBER_MAX_VAR_DEGREE
            || node_count(R, ALGFIBER_MAX_NODES + 1) > ALGFIBER_MAX_NODES) {
            expr_free(R); return false;                             /* budget overrun */
        }
    }

    /* 2. Isolate the real roots of the rational-coefficient univariate R. */
    Expr** cand = NULL; int nc = 0, cc = 0;
    bool ok = rru_collect_roots(R, var, &cand, &nc, &cc, NULL, 0, opts);
    expr_free(R);
    if (!ok) { for (int j = 0; j < nc; j++) expr_free(cand[j]); free(cand); return false; }

    /* 3. Keep exactly the candidates at which the ORIGINAL fibre vanishes at the
     *    true assignment -- discarding the conjugate-spurious roots.  An
     *    undecidable qqbar equality bails (soundness). */
    Expr* zero = expr_new_integer(0);
    bool bail = false;
    for (int j = 0; j < nc; j++) {
        Expr* fval = subst_point(factor, vv, vals, nlev, var, cand[j]);
        int eq = flint_qqbar_equal(fval, zero);
        expr_free(fval);
        if (eq == -1) { bail = true; break; }             /* undecidable -> bail */
        if (eq == 1) prov_push(arr, n, cap, prov, factor_id, cand[j]);
        else expr_free(cand[j]);                          /* spurious -> drop */
        cand[j] = NULL;
    }
    if (bail) for (int j = 0; j < nc; j++) if (cand[j]) expr_free(cand[j]);
    expr_free(zero);
    free(cand);
    return !bail;
#endif
}
