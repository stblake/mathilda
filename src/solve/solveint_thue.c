/*
 * solveint_thue.c
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
 *  Thue equations  F(x, y) == m  (F irreducible homogeneous, deg >=3). *
 *                                                                     *
 *  Detect a single homogeneous binary form of degree n >= 3 equal to a *
 *  constant, extract the form coefficients a_j (of x^(n-j) y^j) and m, *
 *  and hand them to the Tzanakis-de Weger engine (src/solvethue.c).    *
 *  Returns the complete finite solution set (sorted), the empty list   *
 *  (a proof), or NULL to decline -- the engine itself declines unless  *
 *  it can certify completeness.                                        */
Expr* si_solve_thue(SICtx* c) {
    if (c->neq != 1 || c->n != 2) return NULL;
    /* pure equation only: extra bound / ordering / disequation constraints
     * would need a bignum-aware post-filter (a later extension). */
    if (c->n_ord != 0 || c->n_neq != 0 || c->n_abs_ord != 0) return NULL;
    for (int i = 0; i < c->n; i++) if (c->has_lo[i] || c->has_hi[i]) return NULL;

    const MPoly* eq = c->eq[0];
    /* active variables */
    int act[SI_MAX_VARS], na = 0;
    for (int v = 0; v < c->n; v++) if (mpoly_deg_var(eq, v) >= 1) act[na++] = v;
    if (na != 2) return NULL;
    int xv = act[0], yv = act[1];

    /* determine the form degree (max total degree of a non-constant term) */
    int deg = 0;
    for (size_t t = 0; t < eq->n_terms; t++) {
        const int* ex = eq->exps + t * (size_t)eq->n_vars;
        for (int v = 0; v < c->n; v++) if (v != xv && v != yv && ex[v] != 0) return NULL;
        int d = ex[xv] + ex[yv];
        if (d > deg) deg = d;
    }
    if (deg < 3) return NULL;

    /* form[j] = coeff of x^(n-j) y^j (j = 0..deg); constant term = -m */
    mpz_t* form = malloc(sizeof(mpz_t) * (size_t)(deg + 1));
    for (int j = 0; j <= deg; j++) mpz_init(form[j]);
    mpz_t m; mpz_init(m);
    bool ok = true;
    for (size_t t = 0; t < eq->n_terms && ok; t++) {
        const int* ex = eq->exps + t * (size_t)eq->n_vars;
        int dx = ex[xv], dy = ex[yv], d = dx + dy;
        if (d == 0) { mpz_neg(m, eq->coefs[t]); }         /* residual is F - m */
        else if (d == deg) { mpz_set(form[dy], eq->coefs[t]); }
        else ok = false;                                   /* inhomogeneous -> not Thue */
    }
    if (ok && (mpz_sgn(form[0]) == 0 || mpz_sgn(form[deg]) == 0)) ok = false;
    if (ok && mpz_sgn(m) == 0) ok = false;                 /* F == 0 is not a Thue eq */

    Expr* result = NULL;
    if (ok) {
        ThueSol* sols = NULL;
        int cnt = thue_solve_binary_form((const mpz_t*)form, deg, m, &sols);
        if (cnt >= 0) {
            /* sort by (x, y) ascending for a canonical result */
            for (int i = 1; i < cnt; i++) {
                ThueSol key = sols[i]; int j = i - 1;
                while (j >= 0 && (mpz_cmp(sols[j].x, key.x) > 0 ||
                       (mpz_cmp(sols[j].x, key.x) == 0 && mpz_cmp(sols[j].y, key.y) > 0)))
                    { sols[j + 1] = sols[j]; j--; }
                sols[j + 1] = key;
            }
            Expr** rows = (cnt > 0) ? malloc(sizeof(Expr*) * (size_t)cnt) : NULL;
            for (int i = 0; i < cnt; i++) {
                Expr* rx = mk_rule(expr_copy(c->var[xv]), mk_mpz(sols[i].x));
                Expr* ry = mk_rule(expr_copy(c->var[yv]), mk_mpz(sols[i].y));
                rows[i] = mk_list((Expr*[]){ rx, ry }, 2);
            }
            result = mk_list(rows, (size_t)cnt);
            free(rows);
            thue_sols_free(sols, cnt);
        }
    }

    for (int j = 0; j <= deg; j++) mpz_clear(form[j]);
    free(form); mpz_clear(m);
    return result;
}
