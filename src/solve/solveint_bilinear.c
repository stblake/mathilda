/*
 * solveint_bilinear.c
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


/* --- Linear elimination + bilinear divisor solver (Pythagorean). --- */

/* Solve the reduced bilinear equation P (only vars u, w free) by factoring
 * M = b*c - a*d over its divisors.  Eliminated variables are reconstructed
 * numerically, in reverse elimination order, from their stored integer-
 * polynomial formulas.  Returns false if P is not a genuine hyperbola. */
static bool si_bilinear_divisor_solve(const MPoly* P, int u, int w, SICtx* c,
                                      SearchState* st, MPoly* const* formula_mp,
                                      const int* elim_order, int n_elim) {
    mpz_t a, b, cc, d, M, tmp; mpz_init(a); mpz_init(b); mpz_init(cc);
    mpz_init(d); mpz_init(M); mpz_init(tmp);
    int* ex = (int*)calloc((size_t)c->n, sizeof(int));
    ex[u] = 1; ex[w] = 1; { const mpz_t* p = mpoly_get_coef(P, ex); if (p) mpz_set(a, *p); }
    ex[u] = 1; ex[w] = 0; { const mpz_t* p = mpoly_get_coef(P, ex); if (p) mpz_set(b, *p); }
    ex[u] = 0; ex[w] = 1; { const mpz_t* p = mpoly_get_coef(P, ex); if (p) mpz_set(cc, *p); }
    ex[u] = 0; ex[w] = 0; { const mpz_t* p = mpoly_get_coef(P, ex); if (p) mpz_set(d, *p); }
    free(ex);
    mpz_mul(M, b, cc); mpz_mul(tmp, a, d); mpz_sub(M, M, tmp);  /* M = bc - ad */
    bool handled = (mpz_sgn(a) != 0 && mpz_sgn(M) != 0);
    if (handled) {
        st->max_visits = SI_MAX_NODES;
        mpz_t absM; mpz_init(absM); mpz_abs(absM, M);
        Expr* dl = divisors_ordinary(absM);
        mpz_clear(absM);
        if (dl && dl->type == EXPR_FUNCTION) {
            mpz_t P1, Q1, uu, ww; mpz_init(P1); mpz_init(Q1); mpz_init(uu); mpz_init(ww);
            for (size_t i = 0; i < dl->data.function.arg_count && !st->overflow; i++) {
                if (!expr_is_integer_like(dl->data.function.args[i])) continue;
                mpz_t dv; mpz_init(dv);
                expr_to_mpz(dl->data.function.args[i], dv);
                for (int sgn = 1; sgn >= -1; sgn -= 2) {
                    if (++st->visits > st->max_visits) { st->overflow = true; break; }
                    mpz_set(P1, dv); if (sgn < 0) mpz_neg(P1, P1);   /* P1 | M */
                    mpz_divexact(Q1, M, P1);                          /* Q1 = M/P1 */
                    mpz_sub(uu, P1, cc); mpz_sub(ww, Q1, b);          /* u=(P1-c)/a, w=(Q1-b)/a */
                    if (!mpz_divisible_p(uu, a) || !mpz_divisible_p(ww, a)) continue;
                    mpz_divexact(uu, uu, a); mpz_divexact(ww, ww, a);
                    if (!mpz_fits_slong_p(uu) || !mpz_fits_slong_p(ww)) continue;
                    int64_t vals[SI_MAX_VARS];
                    for (int k = 0; k < c->n; k++) vals[k] = 0;
                    vals[u] = mpz_get_si(uu);
                    vals[w] = mpz_get_si(ww);
                    /* Reconstruct eliminated variables in reverse order: a
                     * later-eliminated variable's formula never references an
                     * earlier-eliminated one, so all references are resolved. */
                    bool okrec = true;
                    mpz_t rv; mpz_init(rv);
                    for (int e = n_elim - 1; e >= 0 && okrec; e--) {
                        int v = elim_order[e];
                        si_eval_mpoly(formula_mp[v], vals, rv);
                        if (!mpz_fits_slong_p(rv)) okrec = false;
                        else vals[v] = mpz_get_si(rv);
                    }
                    mpz_clear(rv);
                    if (okrec && si_verify(c, vals)) emit_full(st, vals);
                }
                mpz_clear(dv);
            }
            mpz_clear(P1); mpz_clear(Q1); mpz_clear(uu); mpz_clear(ww);
        }
        if (dl) expr_free(dl);
    }
    mpz_clear(a); mpz_clear(b); mpz_clear(cc); mpz_clear(d); mpz_clear(M); mpz_clear(tmp);
    return handled;
}


/* Try to reduce the equation system to a single bilinear equation in the
 * variables {keepU, keepW} by eliminating every other variable with a
 * unit-coefficient linear equation, then divisor-solve it. */
static bool si_attempt_pair(SICtx* c, SearchState* st, int keepU, int keepW) {
    Expr** conj; int ncj;
    flatten_conjuncts(c->original, &conj, &ncj);
    Expr* eqs[SI_MAX_VARS * 2]; int neq = 0;
    for (int i = 0; i < ncj && neq < (int)(sizeof(eqs)/sizeof(eqs[0])); i++)
        if (is_fun(conj[i], SYM_Equal, 2))
            eqs[neq++] = mk_fn2("Plus", expr_copy(conj[i]->data.function.args[0]),
                mk_fn2("Times", mk_int(-1), expr_copy(conj[i]->data.function.args[1])));

    bool elim_done[SI_MAX_VARS]; for (int i = 0; i < c->n; i++) elim_done[i] = false;
    MPoly* formula_mp[SI_MAX_VARS]; for (int i = 0; i < c->n; i++) formula_mp[i] = NULL;
    int elim_order[SI_MAX_VARS], n_elim = 0;
    bool alive[SI_MAX_VARS * 2]; for (int i = 0; i < neq; i++) alive[i] = true;

    /* Eliminate every non-kept variable that is unit-linear in some equation.
     * The reconstruction formula is kept as an MPoly (evaluated numerically per
     * candidate); only the equation reduction goes through the symbolic
     * substitution, which runs a handful of times, not per candidate. */
    for (bool progress = true; progress; ) {
        progress = false;
        for (int e = 0; e < neq && !progress; e++) {
            if (!alive[e]) continue;
            MPoly* P = si_resid_to_mpoly(eqs[e], c->var, c->n);
            if (!P) continue;
            for (int v = 0; v < c->n && !progress; v++) {
                if (v == keepU || v == keepW || elim_done[v] || mpoly_deg_var(P, v) != 1) continue;
                MPoly* lc = mpoly_coef_of_var(P, v, 1);
                bool unit = (mpoly_total_deg(lc) == 0 && lc->n_terms == 1
                             && (mpz_cmp_si(lc->coefs[0], 1) == 0
                                 || mpz_cmp_si(lc->coefs[0], -1) == 0));
                long coef = (lc->n_terms == 1) ? mpz_get_si(lc->coefs[0]) : 0;
                mpoly_free(lc);
                if (!unit) continue;
                MPoly* rest = mpoly_subst_var_int(P, v, 0);
                MPoly* fpoly = mpoly_scale_si(rest, -coef);   /* v = -(rest)/coef, coef = +/-1 */
                mpoly_free(rest);
                Expr* fexpr = mpoly_to_expr(fpoly, c->var);
                for (int e2 = 0; e2 < neq; e2++) {
                    if (!alive[e2] || e2 == e) continue;
                    Expr* rl = mk_list((Expr*[]){ mk_rule(expr_copy(c->var[v]), expr_copy(fexpr)) }, 1);
                    eqs[e2] = eval_and_free(internal_replace_all((Expr*[]){ eqs[e2], rl }, 2));
                }
                expr_free(fexpr);
                formula_mp[v] = fpoly; elim_order[n_elim++] = v;
                elim_done[v] = true; alive[e] = false; progress = true;
            }
            mpoly_free(P);
        }
    }

    /* Success needs every non-kept variable eliminated and exactly one live
     * equation, bilinear in {keepU, keepW}. */
    bool all_elim = true;
    for (int i = 0; i < c->n; i++)
        if (i != keepU && i != keepW && !elim_done[i]) { all_elim = false; break; }
    int live = -1, nlive = 0;
    for (int e = 0; e < neq; e++) if (alive[e]) { live = e; nlive++; }

    bool handled = false;
    if (all_elim && nlive == 1) {
        MPoly* P = si_resid_to_mpoly(eqs[live], c->var, c->n);
        if (P && mpoly_deg_var(P, keepU) <= 1 && mpoly_deg_var(P, keepW) <= 1
            && mpoly_total_deg(P) <= 2)
            handled = si_bilinear_divisor_solve(P, keepU, keepW, c, st,
                                                formula_mp, elim_order, n_elim);
        if (P) mpoly_free(P);
    }

    for (int i = 0; i < neq; i++) if (eqs[i]) expr_free(eqs[i]);
    for (int i = 0; i < c->n; i++) if (formula_mp[i]) mpoly_free(formula_mp[i]);
    return handled;
}


/* Try every pair of variables to keep; the first that reduces to a genuine
 * bilinear hyperbola wins. */
bool si_solve_linelim_bilinear(SICtx* c, SearchState* st) {
    for (int u = 0; u < c->n; u++)
        for (int w = u + 1; w < c->n; w++)
            if (si_attempt_pair(c, st, u, w)) return true;
    return false;
}
