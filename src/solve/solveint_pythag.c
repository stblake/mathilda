/*
 * solveint_pythag.c
 *
 * Part of the Solve[..., Integers] engine; split out of solveint.c.
 * See solveint_internal.h for the shared SICtx/SearchState substrate.
 *
 * General ("extended") Pythagorean equation
 *
 *     x_1^2 + x_2^2 + ... + x_k^2 == y^2      (k >= 3 summands, so n >= 4 vars)
 *
 * over an unbounded domain.  The homogeneous quadric  Sum x_i^2 = y^2  is a cone
 * with the rational point (1,0,...,0,1), so it is rational (genus 0) and admits a
 * complete parametric solution.  With parameters C[1..k]:
 *
 *     x_i -> 2 C[i] C[k]                       (i = 1 .. k-1)
 *     x_k -> C[1]^2 + ... + C[k-1]^2 - C[k]^2
 *     y   -> C[1]^2 + ... + C[k-1]^2 + C[k]^2
 *
 * Identity:  Sum_{i<k}(2 C[i] C[k])^2 + (T - C[k]^2)^2 = (T + C[k]^2)^2  with
 * T = Sum_{i<k} C[i]^2, since 4 C[k]^2 T + (T - C[k]^2)^2 = (T + C[k]^2)^2.  This
 * is the standard stereographic parametrisation (the k=2 case is the classic
 * (2mn, m^2-n^2, m^2+n^2) Pythagorean triple).  A single family is emitted,
 * matching Mathematica / sympy's representation of the general solution.
 *
 * The k <= 2 cases are the binary / ternary quadratic and are handled by their
 * dedicated solvers (si_solve_pell_parametric, si_solve_ternary_quadratic and
 * the general ternary solver); this path requires k >= 3.  Weighted coefficients
 * (a_i != a_j on the summands) decline -- a documented follow-up.
 */
#include "solveint.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "eval.h"
#include "expr.h"
#include "poly/mpoly.h"
#include "solveint_internal.h"


/* C[i]. */
static Expr* py_Ck(int i) { return mk_fn1("C", mk_int(i)); }
/* C[i]^2. */
static Expr* py_Ck2(int i) { return mk_fn2("Power", py_Ck(i), mk_int(2)); }


/* Detect the homogeneous  Sum_{i} x_i^2 == y^2  shape (after dividing content):
 * every term is a pure square of a single variable, exactly one variable has a
 * coefficient of the opposite sign to the others, and ALL coefficients have equal
 * magnitude.  On success fills summ[0..k-1] (ascending summand indices), *k, and
 * *hyp (the odd-sign variable).  Returns false on any cross / linear / constant
 * term, unequal magnitudes, or fewer than 3 summands. */
static bool py_detect(const MPoly* eq, int n, int* summ, int* k, int* hyp) {
    /* coefficient of x_v^2 for each variable (0 if absent). */
    mpz_t coef[SI_MAX_VARS]; bool present[SI_MAX_VARS];
    for (int v = 0; v < n; v++) { mpz_init_set_ui(coef[v], 0); present[v] = false; }
    bool ok = true;
    for (size_t t = 0; t < eq->n_terms && ok; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        int sqv = -1, others = 0;
        for (int v = 0; v < n; v++) {
            if (ex[v] == 0) continue;
            if (ex[v] == 2 && sqv < 0) sqv = v; else others++;
        }
        if (sqv < 0 || others != 0) { ok = false; break; }   /* cross/linear/const/quartic */
        mpz_set(coef[sqv], eq->coefs[t]); present[sqv] = true;
    }

    int hyp_idx = -1, k_summ = 0;
    mpz_t mag; mpz_init(mag); bool have_mag = false;
    if (ok) {
        /* All present coefficients must share |coef|; exactly one differs in sign. */
        int npos = 0, nneg = 0;
        for (int v = 0; v < n && ok; v++) {
            if (!present[v]) continue;
            mpz_t a; mpz_init(a); mpz_abs(a, coef[v]);
            if (!have_mag) { mpz_set(mag, a); have_mag = true; }
            else if (mpz_cmp(mag, a) != 0) ok = false;
            mpz_clear(a);
            if (mpz_sgn(coef[v]) > 0) npos++; else nneg++;
        }
        /* The odd one out (opposite sign) is the hypotenuse; the rest are summands.
         * Whichever sign is in the minority (and unique) is the hypotenuse. */
        if (ok) {
            int hyp_sign;
            if (npos == 1 && nneg >= 3) hyp_sign = 1;
            else if (nneg == 1 && npos >= 3) hyp_sign = -1;
            else ok = false;
            if (ok) {
                for (int v = 0; v < n; v++) {
                    if (!present[v]) continue;
                    if (mpz_sgn(coef[v]) == hyp_sign) hyp_idx = v;
                    else summ[k_summ++] = v;                  /* ascending by construction */
                }
            }
        }
    }
    mpz_clear(mag);
    for (int v = 0; v < n; v++) mpz_clear(coef[v]);
    if (!ok || hyp_idx < 0 || k_summ < 3) return false;
    *k = k_summ; *hyp = hyp_idx;
    return true;
}


/* General Pythagorean solver.  Returns the owned one-element family List, or NULL
 * to decline (not this shape, bounded / constrained, or k < 3). */
Expr* si_solve_general_pythagorean(SICtx* c) {
    if (c->neq != 1 || c->n < 4) return NULL;
    /* Fully unbounded and unconstrained only: a bounded / ordered instance is left
     * to the finite leaf search, which already enumerates it. */
    for (int i = 0; i < c->n; i++) if (c->has_lo[i] || c->has_hi[i]) return NULL;
    if (c->n_ord != 0 || c->n_neq != 0 || c->n_abs_ord != 0 || !c->all_captured)
        return NULL;

    int summ[SI_MAX_VARS], k = 0, hyp = -1;
    if (!py_detect(c->eq[0], c->n, summ, &k, &hyp)) return NULL;

    /* T = C[1]^2 + ... + C[k-1]^2. */
    Expr* tterms[SI_MAX_VARS]; int nt = 0;
    for (int i = 1; i <= k - 1; i++) tterms[nt++] = py_Ck2(i);
    Expr* T = (nt == 1) ? tterms[0] : expr_new_function(mk_sym("Plus"), tterms, nt);

    /* Rules, placed by variable index for an ascending canonical tuple. */
    Expr* rules[SI_MAX_VARS];
    for (int i = 0; i < c->n; i++) rules[i] = NULL;

    /* Summands 1..k-1 -> 2 C[i] C[k];  summand k -> T - C[k]^2. */
    for (int i = 1; i <= k - 1; i++) {
        Expr* val = expr_new_function(mk_sym("Times"),
            (Expr*[]){ mk_int(2), py_Ck(i), py_Ck(k) }, 3);
        rules[summ[i - 1]] = mk_rule(expr_copy(c->var[summ[i - 1]]), val);
    }
    rules[summ[k - 1]] = mk_rule(expr_copy(c->var[summ[k - 1]]),
        mk_fn2("Plus", expr_copy(T), mk_fn2("Times", mk_int(-1), py_Ck2(k))));

    /* Hypotenuse -> T + C[k]^2. */
    rules[hyp] = mk_rule(expr_copy(c->var[hyp]), mk_fn2("Plus", T, py_Ck2(k)));

    Expr* rlist[SI_MAX_VARS]; int nr = 0;
    for (int i = 0; i < c->n; i++) if (rules[i]) rlist[nr++] = rules[i];
    Expr* tuple = mk_list(rlist, (size_t)nr);
    return eval_and_free(mk_list((Expr*[]){ tuple }, 1));
}
