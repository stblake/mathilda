/*
 * solveint_powerleaf.c
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
 *  A4: non-polynomial bounded power-leaf.                             *
 *                                                                     *
 *  An equation like  n! + 1 == m^2  cannot become an MPoly (Factorial *
 *  is not a monomial), so Stage A declines it.  When one side is a     *
 *  pure power  m^e  of a leaf variable m that appears nowhere else,    *
 *  and every OTHER variable is bounded to a small box, we enumerate    *
 *  those variables, evaluate the other side numerically through the    *
 *  interpreter (which knows Factorial, Binomial, ...), and solve m by  *
 *  an exact integer e-th root.  m itself may be unbounded (it is       *
 *  determined, not enumerated), so its value may be a bignum.          *
 * ------------------------------------------------------------------ */
static bool si_expr_mentions(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == name;
    if (e->type != EXPR_FUNCTION) return false;
    if (si_expr_mentions(e->data.function.head, name)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (si_expr_mentions(e->data.function.args[i], name)) return true;
    return false;
}


/* Is `side` a bare  Power[m, e]  with m a solve variable (index -> *leaf),
 * e a positive integer literal (-> *e), and m absent from `other`? */
static bool si_power_leaf_side(Expr* side, Expr* other, Expr** var, int n,
                               int* leaf, int* e) {
    if (!is_fun(side, SYM_Power, 2)) return false;
    Expr* base = side->data.function.args[0];
    Expr* exq  = side->data.function.args[1];
    if (base->type != EXPR_SYMBOL || exq->type != EXPR_INTEGER || exq->data.integer < 1)
        return false;
    int mi = -1;
    for (int i = 0; i < n; i++)
        if (var[i]->data.symbol.name == base->data.symbol.name) mi = i;
    if (mi < 0 || si_expr_mentions(other, base->data.symbol.name)) return false;
    *leaf = mi; *e = (int)exq->data.integer;
    return true;
}

Expr* si_solve_bounded_powerleaf(Expr* expr, Expr** var, int n) {
    Expr** conj; int ncj;
    flatten_conjuncts(expr, &conj, &ncj);
    Expr* eqn = NULL;
    for (int i = 0; i < ncj; i++)
        if (is_fun(conj[i], SYM_Equal, 2)) { eqn = conj[i]; break; }
    if (!eqn) return NULL;

    SICtx c; memset(&c, 0, sizeof(c)); c.var = var; c.n = n; c.original = expr; c.all_captured = true;

    /* Gate: engage only for a NON-polynomial equation (the polynomial leaf
     * search is better and already handles the polynomial power-leaf case). */
    MPoly* testQ = relation_to_mpoly(&c, eqn->data.function.args[0], eqn->data.function.args[1]);
    if (testQ) { mpoly_free(testQ); for (int i = 0; i < c.nbc; i++) mpoly_free(c.bc[i].Q); return NULL; }

    /* Identify the power-leaf side and the evaluable target side. */
    Expr *A = eqn->data.function.args[0], *B = eqn->data.function.args[1];
    int leaf = -1, e = 0; Expr* target = NULL;
    if (si_power_leaf_side(A, B, var, n, &leaf, &e)) target = B;
    else if (si_power_leaf_side(B, A, var, n, &leaf, &e)) target = A;
    if (leaf < 0) { for (int i = 0; i < c.nbc; i++) mpoly_free(c.bc[i].Q); return NULL; }

    /* Parse constraints for the OTHER variables' bounds (skip the equation). */
    for (int i = 0; i < ncj; i++) if (conj[i] != eqn) classify_conjunct(&c, conj[i]);

    /* Every non-leaf variable must be finitely bounded and enumerable. */
    int ov[SI_MAX_VARS], nov = 0; long double box = 1.0L;
    bool ok = true;
    for (int i = 0; i < n; i++) {
        if (i == leaf) continue;
        if (!(c.has_lo[i] && c.has_hi[i]) || c.hi[i] < c.lo[i]) { ok = false; break; }
        ov[nov++] = i;
        box *= (long double)(c.hi[i] - c.lo[i] + 1);
    }
    if (!ok || box > 5.0e6L) { for (int i = 0; i < c.nbc; i++) mpoly_free(c.bc[i].Q); return NULL; }

    /* Enumerate the other variables; solve m by an exact e-th root. */
    Expr** sols = NULL; int nsol = 0, cap = 0;
    int64_t val[SI_MAX_VARS];
    for (int i = 0; i < n; i++) val[i] = 0;
    for (int q = 0; q < nov; q++) val[ov[q]] = c.lo[ov[q]];
    mpz_t T, root, chk; mpz_init(root); mpz_init(chk);  /* T is init-set by expr_to_mpz */
    for (;;) {
        /* target evaluated at the current other-variable assignment. */
        Expr** rules = (Expr**)malloc(sizeof(Expr*) * (size_t)(nov ? nov : 1));
        for (int q = 0; q < nov; q++) rules[q] = mk_rule(expr_copy(var[ov[q]]), mk_int(val[ov[q]]));
        Expr* rl = mk_list(rules, (size_t)nov); free(rules);
        Expr* tv = eval_and_free(internal_replace_all((Expr*[]){ expr_copy(target), rl }, 2));
        bool got = expr_is_integer_like(tv);
        if (got) expr_to_mpz(tv, T);    /* initialises T */
        expr_free(tv);

        if (got) {
            /* m^e == T.  Even e needs T >= 0 and yields +/- root; odd e keeps sign. */
            bool have_m = false; int msign = 1;
            if (e % 2 == 0) {
                if (mpz_sgn(T) >= 0 && mpz_root(root, T, (unsigned long)e)) have_m = true;
            } else {
                mpz_abs(chk, T);
                if (mpz_root(root, chk, (unsigned long)e)) { have_m = true; msign = mpz_sgn(T) < 0 ? -1 : 1; }
            }
            if (have_m) {
                for (int sg = 1; sg >= -1; sg -= 2) {
                    if (e % 2 == 1 && sg != msign) continue;      /* odd: single sign */
                    if (e % 2 == 0 && sg < 0 && mpz_sgn(root) == 0) continue; /* avoid -0 dup */
                    mpz_t mval; mpz_init(mval); mpz_set(mval, root); if (sg < 0) mpz_neg(mval, mval);
                    /* Verify the full conjunction symbolically (handles Factorial etc.). */
                    Expr** vr = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
                    for (int i = 0; i < n; i++)
                        vr[i] = (i == leaf) ? mk_rule(expr_copy(var[i]), mk_mpz(mval))
                                            : mk_rule(expr_copy(var[i]), mk_int(val[i]));
                    Expr* vrl = mk_list(vr, (size_t)n); free(vr);
                    Expr* chk2 = eval_and_free(internal_replace_all((Expr*[]){ expr_copy(expr), vrl }, 2));
                    bool good = is_sym(chk2, SYM_True);
                    expr_free(chk2);
                    if (good) {
                        Expr** tr = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
                        for (int i = 0; i < n; i++)
                            tr[i] = (i == leaf) ? mk_rule(expr_copy(var[i]), mk_mpz(mval))
                                                : mk_rule(expr_copy(var[i]), mk_int(val[i]));
                        if (nsol == cap) { cap = cap ? cap * 2 : 8; sols = realloc(sols, sizeof(Expr*) * (size_t)cap); }
                        sols[nsol++] = mk_list(tr, (size_t)n); free(tr);
                    }
                    mpz_clear(mval);
                }
            }
        }
        if (got) mpz_clear(T);          /* T was init-set this iteration */
        int q = 0;
        for (; q < nov; q++) { if (++val[ov[q]] <= c.hi[ov[q]]) break; val[ov[q]] = c.lo[ov[q]]; }
        if (q == nov) break;
    }
    mpz_clear(root); mpz_clear(chk);
    for (int i = 0; i < c.nbc; i++) mpoly_free(c.bc[i].Q);
    for (int i = 0; i < c.neq; i++) mpoly_free(c.eq[i]);

    Expr* result = mk_list(sols, (size_t)nsol);
    free(sols);
    return result;
}
