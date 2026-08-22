/*
 * solveint_exp.c
 *
 * Part of the Solve[..., Integers] engine; split out of solveint.c.
 * See solveint_internal.h for the shared SICtx/SearchState substrate.
 */
#include "solveint.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "internal.h"
#include "sym_names.h"
#include "symtab.h"
#include "checked_int.h"
#include "poly/mpoly.h"
#include "numbertheory/numbertheory_internal.h"
#include "linalg/hnf.h"
#include "solvethue.h"
#include "solveint_internal.h"


/* ------------------------------------------------------------------ *
 *  Exponential Diophantine  (variable exponents).                     *
 *                                                                     *
 *  Equations like  x^a - y^b == 1  cannot be represented as polynomials
 *  (the exponent is a variable), so they are handled here, before the  *
 *  MPoly stage.  A fully bounded box is enumerated exactly; the        *
 *  unbounded Catalan shape  x^a - y^b == +/-1  (bases and exponents    *
 *  >= 2) is settled by Mihailescu's theorem: the only solution is      *
 *  3^2 - 2^3 = 1.                                                      *
 * ------------------------------------------------------------------ */

/* SIExpTerm is defined in solveint_internal.h (shared with solveint_rn.c). */

/* Collect additive power terms of `e` (multiplied by outer `sign`) into t[],
 * folding integer constants into K.  Returns false on an unrecognised shape. */
bool si_exp_collect(Expr* e, int sign, Expr** var, int n,
                    SIExpTerm* t, int* nt, mpz_t K) {
    if (!e) return false;
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Plus) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!si_exp_collect(e->data.function.args[i], sign, var, n, t, nt, K)) return false;
        return true;
    }
    int coef = sign; Expr* powfac = e;
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Times
        && e->data.function.arg_count == 2
        && e->data.function.args[0]->type == EXPR_INTEGER) {
        coef = sign * (int)e->data.function.args[0]->data.integer;
        powfac = e->data.function.args[1];
    }
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) {   /* constant */
        mpz_t v; expr_to_mpz(e, v);
        if (sign < 0) mpz_sub(K, K, v); else mpz_add(K, K, v);
        mpz_clear(v);
        return true;
    }
    /* powfac must be base^exp, a bare variable (exp 1), or Power[base, exp]. */
    Expr *base = NULL, *exq = NULL;
    if (powfac->type == EXPR_SYMBOL) { base = powfac; exq = NULL; }
    else if (powfac->type == EXPR_FUNCTION && powfac->data.function.head->type == EXPR_SYMBOL
             && powfac->data.function.head->data.symbol.name == SYM_Power
             && powfac->data.function.arg_count == 2) {
        base = powfac->data.function.args[0];
        exq  = powfac->data.function.args[1];
    } else return false;

    if (*nt >= 8) return false;
    SIExpTerm* pt = &t[*nt];
    pt->coef = coef; pt->bvar = -1; pt->bconst = 0; pt->evar = -1; pt->econst = 1;
    if (base->type == EXPR_SYMBOL) {
        int bi = -1; for (int i = 0; i < n; i++) if (var[i]->data.symbol.name == base->data.symbol.name) bi = i;
        if (bi < 0) return false;
        pt->bvar = bi;
    } else if (base->type == EXPR_INTEGER) pt->bconst = base->data.integer;
    else return false;
    if (exq == NULL) pt->econst = 1;
    else if (exq->type == EXPR_SYMBOL) {
        int ei = -1; for (int i = 0; i < n; i++) if (var[i]->data.symbol.name == exq->data.symbol.name) ei = i;
        if (ei < 0) return false;
        pt->evar = ei;
    } else if (exq->type == EXPR_INTEGER) pt->econst = exq->data.integer;
    else return false;
    (*nt)++;
    return true;
}


/* Build {{v -> val, ...}} for the single tuple `vals`. */
Expr* si_exp_one_tuple(Expr** var, const int64_t* vals, int n) {
    Expr** rules = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++) rules[i] = mk_rule(expr_copy(var[i]), mk_int(vals[i]));
    Expr* tup = mk_list(rules, (size_t)n); free(rules);
    return tup;
}


/* Is  target == base^exp  for some exp >= 1?  (base >= 2, target >= 1.) */
static bool si_int_ppow(int64_t base, int64_t target, int64_t* exp) {
    if (base < 2 || target < 2) return false;        /* base^exp >= 2 for exp >= 1 */
    int64_t v = base, e = 1;
    while (v < target) {
        if (v > INT64_MAX / base) return false;
        v *= base; e++;
    }
    if (v == target) { *exp = e; return true; }
    return false;
}


/* Verify a full integer assignment `vals` against the original conjunction
 * `expr` symbolically (used by paths whose equation is not an MPoly). */
bool si_verify_symbolic(Expr* expr, Expr** var, const int64_t* vals, int n) {
    Expr** rules = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++) rules[i] = mk_rule(expr_copy(var[i]), mk_int(vals[i]));
    Expr* rl = mk_list(rules, (size_t)n); free(rules);
    Expr* out = eval_and_free(internal_replace_all((Expr*[]){ expr_copy(expr), rl }, 2));
    bool ok = is_sym(out, SYM_True);
    expr_free(out);
    return ok;
}


/* Try to solve an exponential Diophantine equation.  Returns the owned result
 * List, or NULL to decline (not this shape, or not decidable here). */
Expr* si_solve_exponential(Expr* expr, Expr** var, int n) {
    Expr** conj; int ncj;
    flatten_conjuncts(expr, &conj, &ncj);

    /* Find the (first) equation whose residual carries a variable exponent. */
    Expr* eqn = NULL;
    for (int i = 0; i < ncj; i++)
        if (is_fun(conj[i], SYM_Equal, 2)) { eqn = conj[i]; break; }
    if (!eqn) return NULL;

    SIExpTerm t[8]; int nt = 0;
    mpz_t K; mpz_init_set_ui(K, 0);
    bool okp = si_exp_collect(eqn->data.function.args[0], +1, var, n, t, &nt, K)
            && si_exp_collect(eqn->data.function.args[1], -1, var, n, t, &nt, K);
    mpz_neg(K, K);                                     /* sum of power terms == K */
    bool has_var_exp = false;
    for (int i = 0; i < nt; i++) if (t[i].evar >= 0) has_var_exp = true;
    if (!okp || !has_var_exp) { mpz_clear(K); return NULL; }

    /* Parse the constraints for per-variable bounds (skip the equation). */
    SICtx c; memset(&c, 0, sizeof(c)); c.var = var; c.n = n; c.original = expr; c.all_captured = true;
    for (int i = 0; i < ncj; i++) if (conj[i] != eqn) classify_conjunct(&c, conj[i]);

    Expr* result = NULL;

    /* Catalan / Mihailescu shape: p^m - q^n == +/-1, bases and exponents all
     * variables constrained to >= 2.  Unique solution 3^2 - 2^3 = 1. */
    if (nt == 2 && (mpz_cmp_si(K, 1) == 0 || mpz_cmp_si(K, -1) == 0)) {
        int ip = (t[0].coef > 0) ? 0 : 1, in = 1 - ip;   /* + and - terms */
        bool shape = t[ip].coef == 1 && t[in].coef == -1
            && t[ip].bvar >= 0 && t[ip].evar >= 0 && t[in].bvar >= 0 && t[in].evar >= 0;
        int b9 = -1, e9 = -1, b8 = -1, e8 = -1;          /* the ^ that is 9 (3^2) / 8 (2^3) */
        if (shape) {
            if (mpz_cmp_si(K, 1) == 0) { b9 = t[ip].bvar; e9 = t[ip].evar; b8 = t[in].bvar; e8 = t[in].evar; }
            else                        { b9 = t[in].bvar; e9 = t[in].evar; b8 = t[ip].bvar; e8 = t[ip].evar; }
            /* every base/exponent must admit the required 3,2 / 2,3 and be >= 2. */
            int need_lo[4] = {2, 2, 2, 2};
            int idx[4] = {b9, e9, b8, e8}; int val[4] = {3, 2, 2, 3};
            bool feas = true;
            for (int q = 0; q < 4 && feas; q++) {
                int v = idx[q];
                if (!(c.has_lo[v] && c.lo[v] >= need_lo[q])) feas = false;   /* bases,exps >= 2 */
                if (c.has_lo[v] && val[q] < c.lo[v]) feas = false;
                if (c.has_hi[v] && val[q] > c.hi[v]) feas = false;
            }
            /* the four variables must be distinct (x,y,a,b). */
            if (b9==e9||b9==b8||b9==e8||e9==b8||e9==e8||b8==e8) feas = false;
            if (feas) {
                int64_t vals[SI_MAX_VARS];
                for (int i = 0; i < n; i++) vals[i] = 0;
                vals[b9] = 3; vals[e9] = 2; vals[b8] = 2; vals[e8] = 3;
                Expr* tup = si_exp_one_tuple(var, vals, n);
                result = mk_list((Expr*[]){ tup }, 1);
            } else {
                result = mk_list(NULL, 0);                /* provably no solution */
            }
        }
    }

    /* Fixed-base Pillai:  P^m - Q^n == +/-1  with CONSTANT bases P,Q >= 2 and
     * VARIABLE exponents m,n.  Sound & complete via Mihailescu: solutions with
     * both exponents >= 2 are consecutive perfect powers, i.e. only 3^2 - 2^3 = 1;
     * the exponent-1 cases are  Q^n = P-1  and  P^m = Q+1.  Settles 3^m - 2^n == 1
     * -> {(1,1),(2,3)} even though m,n are unbounded. */
    if (!result && n == 2 && (mpz_cmp_si(K, 1) == 0 || mpz_cmp_si(K, -1) == 0)) {
        int ip = (t[0].coef > 0) ? 0 : 1, in = 1 - ip;
        bool shape = t[ip].coef == 1 && t[in].coef == -1
            && t[ip].bvar < 0 && t[in].bvar < 0          /* constant bases */
            && t[ip].evar >= 0 && t[in].evar >= 0         /* variable exponents */
            && t[ip].bconst >= 2 && t[in].bconst >= 2;
        if (shape) {
            /* Normalise to  P^a - Q^b == 1  (av,bv = exponent variable indices). */
            int64_t P, Q; int av, bv;
            if (mpz_cmp_si(K, 1) == 0) { P = t[ip].bconst; av = t[ip].evar; Q = t[in].bconst; bv = t[in].evar; }
            else                       { P = t[in].bconst; av = t[in].evar; Q = t[ip].bconst; bv = t[ip].evar; }
            if (av != bv) {
                int64_t sa[8], sb[8]; int ns = 0; int64_t e;
                if (P == 3 && Q == 2) { sa[ns] = 2; sb[ns] = 3; ns++; }       /* both exps >= 2 */
                if (si_int_ppow(Q, P - 1, &e)) { sa[ns] = 1; sb[ns] = e; ns++; } /* a == 1 */
                if (si_int_ppow(P, Q + 1, &e)) { sa[ns] = e; sb[ns] = 1; ns++; } /* b == 1 */
                for (int u = 0; u + 1 < ns; u++)         /* sort ascending by (a, b) */
                    for (int w = u + 1; w < ns; w++)
                        if (sa[w] < sa[u] || (sa[w] == sa[u] && sb[w] < sb[u])) {
                            int64_t ta = sa[u], tb = sb[u]; sa[u] = sa[w]; sb[u] = sb[w]; sa[w] = ta; sb[w] = tb;
                        }
                Expr** sols = NULL; int nsol = 0, cap = 0;
                for (int s = 0; s < ns; s++) {
                    int64_t vals[SI_MAX_VARS];
                    for (int i = 0; i < n; i++) vals[i] = 0;
                    vals[av] = sa[s]; vals[bv] = sb[s];
                    bool dup = false;                     /* de-dup (a=b=1 arises twice) */
                    for (int s2 = 0; s2 < s && !dup; s2++) if (sa[s2] == sa[s] && sb[s2] == sb[s]) dup = true;
                    if (dup) continue;
                    if (si_verify_symbolic(expr, var, vals, n)) {
                        if (nsol == cap) { cap = cap ? cap * 2 : 4; sols = realloc(sols, sizeof(Expr*) * (size_t)cap); }
                        sols[nsol++] = si_exp_one_tuple(var, vals, n);
                    }
                }
                result = mk_list(sols, (size_t)nsol);
                free(sols);
            }
        }
    }

    /* Otherwise: fully bounded enumeration.  Every base and exponent variable
     * must have a finite box; enumerate and test exactly. */
    if (!result) {
        int uvar[SI_MAX_VARS], nu = 0;
        bool bounded = true;
        for (int i = 0; i < n; i++) {
            bool used = false;
            for (int k = 0; k < nt; k++) if (t[k].bvar == i || t[k].evar == i) used = true;
            if (!used) continue;
            if (!(c.has_lo[i] && c.has_hi[i])) { bounded = false; break; }
            uvar[nu++] = i;
        }
        long double space = 1.0L;
        for (int q = 0; q < nu; q++) space *= (long double)(c.hi[uvar[q]] - c.lo[uvar[q]] + 1);
        if (bounded && space <= 5.0e7L) {
            Expr** sols = NULL; int nsol = 0, cap = 0;
            int64_t val[SI_MAX_VARS];
            for (int i = 0; i < n; i++) val[i] = 0;
            for (int q = 0; q < nu; q++) val[uvar[q]] = c.lo[uvar[q]];
            mpz_t acc, term, be; mpz_init(acc); mpz_init(term); mpz_init(be);
            for (;;) {
                mpz_set_ui(acc, 0);
                for (int k = 0; k < nt; k++) {
                    int64_t bb = (t[k].bvar >= 0) ? val[t[k].bvar] : t[k].bconst;
                    int64_t ee = (t[k].evar >= 0) ? val[t[k].evar] : t[k].econst;
                    if (ee < 0) { mpz_set_ui(be, 0); }     /* negative exponent -> skip term as 0 (nonintegral) */
                    else { mpz_set_si(be, bb); mpz_pow_ui(be, be, (unsigned long)ee); }
                    mpz_mul_si(term, be, t[k].coef); mpz_add(acc, acc, term);
                }
                if (mpz_cmp(acc, K) == 0) {
                    if (nsol == cap) { cap = cap ? cap * 2 : 8; sols = realloc(sols, sizeof(Expr*) * (size_t)cap); }
                    sols[nsol++] = si_exp_one_tuple(var, val, n);
                }
                int q = 0;
                for (; q < nu; q++) { if (++val[uvar[q]] <= c.hi[uvar[q]]) break; val[uvar[q]] = c.lo[uvar[q]]; }
                if (q == nu) break;
            }
            mpz_clear(acc); mpz_clear(term); mpz_clear(be);
            result = mk_list(sols, (size_t)nsol);
            free(sols);
        }
    }

    for (int i = 0; i < c.nbc; i++) mpoly_free(c.bc[i].Q);
    for (int i = 0; i < c.neq; i++) mpoly_free(c.eq[i]);
    mpz_clear(K);
    return result;
}
