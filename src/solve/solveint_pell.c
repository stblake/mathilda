/*
 * solveint_pell.c
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


/* --- Pell equation  x^2 - D y^2 == N  (N = +/-1) via continued fractions. --- */

/* Recognise a single equation whose only terms are xside^2, yside^2 and a
 * constant, with opposite-sign squares and the +square coefficient equal to 1.
 * Fills xside, yside, D (>0, non-square) and N (+/-1). */
static bool si_pell_detect(const MPoly* eq, int n, int* xside, int* yside,
                           mpz_t D, mpz_t N) {
    int active[SI_MAX_VARS], ka = 0;
    for (int v = 0; v < n; v++) if (mpoly_deg_var(eq, v) >= 1) active[ka++] = v;
    if (ka != 2) return false;
    int A = active[0], B = active[1];
    mpz_t cA, cB, e; mpz_init_set_ui(cA, 0); mpz_init_set_ui(cB, 0); mpz_init_set_ui(e, 0);
    bool ok = true;
    for (size_t t = 0; t < eq->n_terms && ok; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        for (int v = 0; v < n; v++) if (v != A && v != B && ex[v] != 0) ok = false;
        if (!ok) break;
        int dA = ex[A], dB = ex[B];
        if (dA == 2 && dB == 0) mpz_set(cA, eq->coefs[t]);
        else if (dA == 0 && dB == 2) mpz_set(cB, eq->coefs[t]);
        else if (dA == 0 && dB == 0) mpz_set(e, eq->coefs[t]);
        else ok = false;                                     /* linear or cross term */
    }
    if (ok && (mpz_sgn(cA) == 0 || mpz_sgn(cB) == 0 || mpz_sgn(cA) == mpz_sgn(cB)))
        ok = false;
    if (ok) {
        if (mpz_sgn(cA) > 0) {
            if (mpz_cmp_si(cA, 1) != 0) ok = false;
            else { *xside = A; *yside = B; mpz_neg(D, cB); mpz_neg(N, e); }
        } else {
            if (mpz_cmp_si(cB, 1) != 0) ok = false;
            else { *xside = B; *yside = A; mpz_neg(D, cA); mpz_neg(N, e); }
        }
    }
    mpz_clear(cA); mpz_clear(cB); mpz_clear(e);
    if (!ok) return false;
    if (mpz_sgn(D) <= 0 || mpz_perfect_square_p(D)) return false;
    return (mpz_cmp_si(N, 1) == 0 || mpz_cmp_si(N, -1) == 0);
}


/* Continued fraction of sqrt(D): sets the fundamental unit (u,v) of
 * x^2 - D y^2 = 1, and the first convergent (bx,by) reaching value N.
 * Returns true if the N-base was found. */
bool si_pell_cf(const mpz_t D, const mpz_t N, mpz_t u, mpz_t v,
                       mpz_t bx, mpz_t by) {
    mpz_t a0, m, d, a, hprev, h, kprev, k, val, t2;
    mpz_init(a0); mpz_sqrt(a0, D);
    mpz_init_set_ui(m, 0); mpz_init_set_ui(d, 1); mpz_init_set(a, a0);
    mpz_init_set_ui(hprev, 1); mpz_init_set(h, a0);
    mpz_init_set_ui(kprev, 0); mpz_init_set_ui(k, 1);
    mpz_init(val); mpz_init(t2);
    bool haveU = false, haveB = false;
    for (int iter = 0; iter < 200000 && !haveU; iter++) {
        mpz_mul(val, h, h); mpz_mul(t2, k, k); mpz_mul(t2, t2, D);
        mpz_sub(val, val, t2);                               /* h^2 - D k^2 */
        if (!haveB && mpz_cmp(val, N) == 0) { mpz_set(bx, h); mpz_set(by, k); haveB = true; }
        if (mpz_cmp_si(val, 1) == 0) { mpz_set(u, h); mpz_set(v, k); haveU = true; break; }
        /* next term: m=d*a-m; d=(D-m^2)/d; a=floor((a0+m)/d) */
        mpz_mul(t2, d, a); mpz_sub(m, t2, m);
        mpz_mul(t2, m, m); mpz_sub(t2, D, t2); mpz_divexact(d, t2, d);
        mpz_add(t2, a0, m); mpz_fdiv_q(a, t2, d);
        mpz_mul(t2, a, h); mpz_add(t2, t2, hprev); mpz_set(hprev, h); mpz_set(h, t2);
        mpz_mul(t2, a, k); mpz_add(t2, t2, kprev); mpz_set(kprev, k); mpz_set(k, t2);
    }
    mpz_clear(a0); mpz_clear(m); mpz_clear(d); mpz_clear(a);
    mpz_clear(hprev); mpz_clear(h); mpz_clear(kprev); mpz_clear(k);
    mpz_clear(val); mpz_clear(t2);
    return haveB;
}


/* Emit (+/-X, +/-Y) at (xside, yside), each sign combo verified. */
static void si_pell_emit(SICtx* c, SearchState* st, int xs, int ys,
                         const mpz_t X, const mpz_t Y) {
    if (!mpz_fits_slong_p(X) || !mpz_fits_slong_p(Y)) return;
    int64_t xa = mpz_get_si(X), ya = mpz_get_si(Y);
    for (int sx = 1; sx >= -1; sx -= 2)
        for (int sy = 1; sy >= -1; sy -= 2) {
            int64_t vals[SI_MAX_VARS];
            for (int i = 0; i < c->n; i++) vals[i] = 0;
            vals[xs] = sx * xa; vals[ys] = sy * ya;
            if (si_verify(c, vals)) emit_full(st, vals);
        }
}

bool si_solve_pell(SICtx* c, SearchState* st) {
    if (c->neq != 1) return false;
    int xs, ys; mpz_t D, N; mpz_init(D); mpz_init(N);
    if (!si_pell_detect(c->eq[0], c->n, &xs, &ys, D, N)) { mpz_clear(D); mpz_clear(N); return false; }
    /* A finite bound is required to terminate the (otherwise infinite) orbit. */
    if (!c->has_hi[xs] && !c->has_hi[ys]) { mpz_clear(D); mpz_clear(N); return false; }

    mpz_t u, v, bx, by; mpz_init(u); mpz_init(v); mpz_init(bx); mpz_init(by);
    bool haveB = si_pell_cf(D, N, u, v, bx, by);
    st->max_visits = SI_MAX_NODES;

    /* Base class of solutions to x^2 - D y^2 = N. */
    mpz_t cx, cy, nx, ny, t1, t2; mpz_init(cx); mpz_init(cy);
    mpz_init(nx); mpz_init(ny); mpz_init(t1); mpz_init(t2);
    bool go = true;
    if (mpz_cmp_si(N, 1) == 0) {
        mpz_set_ui(cx, 1); mpz_set_ui(cy, 0);                /* trivial (1,0) */
        si_pell_emit(c, st, xs, ys, cx, cy);
        mpz_set(cx, u); mpz_set(cy, v);                      /* first nontrivial = U */
    } else if (haveB) {
        mpz_set(cx, bx); mpz_set(cy, by);                    /* N = -1 base */
    } else {
        go = false;                                          /* x^2 - D y^2 = -1 unsolvable */
    }

    for (int guard = 0; go && guard < 100000 && !st->overflow; guard++) {
        bool past = (c->has_hi[xs] && mpz_cmp_si(cx, c->hi[xs]) > 0)
                 || (c->has_hi[ys] && mpz_cmp_si(cy, c->hi[ys]) > 0);
        if (past) break;
        si_pell_emit(c, st, xs, ys, cx, cy);
        /* compose with U: (x,y) -> (x u + D y v, x v + y u) */
        mpz_mul(nx, cx, u); mpz_mul(t1, cy, v); mpz_mul(t1, t1, D); mpz_add(nx, nx, t1);
        mpz_mul(ny, cx, v); mpz_mul(t1, cy, u); mpz_add(ny, ny, t1);
        mpz_set(cx, nx); mpz_set(cy, ny);
    }

    mpz_clear(cx); mpz_clear(cy); mpz_clear(nx); mpz_clear(ny); mpz_clear(t1); mpz_clear(t2);
    mpz_clear(u); mpz_clear(v); mpz_clear(bx); mpz_clear(by);
    mpz_clear(D); mpz_clear(N);
    return true;
}


/* --- Unbounded Pell  x^2 - D y^2 == 1,  x > 0 && y > 0  -> parametric family. ---
 *
 * With no finite bound the positive-orthant solution set is infinite: the k-th
 * solution (k >= 1) is  x_k + y_k sqrt(D) = (x1 + y1 sqrt(D))^k  for the
 * fundamental unit (x1, y1).  Emit the closed form
 *   x -> ((x1+y1 Sqrt[D])^C[1] + (x1-y1 Sqrt[D])^C[1]) / 2,
 *   y -> ((x1+y1 Sqrt[D])^C[1] - (x1-y1 Sqrt[D])^C[1]) / (2 Sqrt[D]),
 * each a ConditionalExpression valid for the integer parameter C[1] >= 1 --
 * mirroring Mathematica's representation.  Returns the owned family, or NULL to
 * decline (bounded, wrong sign pattern, or N = -1, all handled elsewhere). */
Expr* si_solve_pell_parametric(SICtx* c) {
    if (c->neq != 1) return NULL;
    int xs, ys; mpz_t D, N; mpz_init(D); mpz_init(N);
    if (!si_pell_detect(c->eq[0], c->n, &xs, &ys, D, N)) { mpz_clear(D); mpz_clear(N); return NULL; }
    bool unbounded = !c->has_hi[xs] && !c->has_hi[ys];
    bool posx = c->has_lo[xs] && c->lo[xs] >= 1;
    bool posy = c->has_lo[ys] && c->lo[ys] >= 1;
    if (!(unbounded && posx && posy && mpz_cmp_si(N, 1) == 0)
        || c->n_ord != 0 || c->n_neq != 0 || c->n_abs_ord != 0 || !c->all_captured) {
        mpz_clear(D); mpz_clear(N); return NULL;
    }
    mpz_t u, v, bx, by; mpz_init(u); mpz_init(v); mpz_init(bx); mpz_init(by);
    si_pell_cf(D, N, u, v, bx, by);                  /* fundamental unit (x1,y1) = (u,v) */

    Expr* sqrtD = mk_fn1("Sqrt", mk_mpz(D));
    mpz_t ny; mpz_init(ny); mpz_neg(ny, v);
    /* uu = x1 + y1 Sqrt[D],  vv = x1 - y1 Sqrt[D]. */
    Expr* uu = mk_fn2("Plus", mk_mpz(u), mk_fn2("Times", mk_mpz(v), expr_copy(sqrtD)));
    Expr* vv = mk_fn2("Plus", mk_mpz(u), mk_fn2("Times", mk_mpz(ny), expr_copy(sqrtD)));
    mpz_clear(ny);
    Expr* uk = mk_fn2("Power", uu, mk_fn1("C", mk_int(1)));
    Expr* vk = mk_fn2("Power", vv, mk_fn1("C", mk_int(1)));
    /* xval = (uk + vk)/2 */
    Expr* xval = mk_fn2("Times",
        mk_fn2("Plus", expr_copy(uk), expr_copy(vk)),
        mk_fn2("Power", mk_int(2), mk_int(-1)));
    /* yval = (uk - vk)/(2 Sqrt[D]) */
    Expr* yval = mk_fn2("Times",
        mk_fn2("Plus", uk, mk_fn2("Times", mk_int(-1), vk)),
        mk_fn2("Power", mk_fn2("Times", mk_int(2), expr_copy(sqrtD)), mk_int(-1)));
    expr_free(sqrtD);

    Expr* cond = mk_fn2("GreaterEqual", mk_fn1("C", mk_int(1)), mk_int(1));
    Expr* cex = mk_fn2("ConditionalExpression", xval, expr_copy(cond));
    Expr* cey = mk_fn2("ConditionalExpression", yval, cond);

    Expr* rules[SI_MAX_VARS];
    for (int i = 0; i < c->n; i++) rules[i] = NULL;
    rules[xs] = mk_rule(expr_copy(c->var[xs]), cex);
    rules[ys] = mk_rule(expr_copy(c->var[ys]), cey);
    /* Emit in variable order (xs, ys are the only two active). */
    Expr* rlist[SI_MAX_VARS]; int nr = 0;
    for (int i = 0; i < c->n; i++) if (rules[i]) rlist[nr++] = rules[i];
    Expr* tuple = mk_list(rlist, (size_t)nr);
    Expr* result = eval_and_free(mk_list((Expr*[]){ tuple }, 1));

    mpz_clear(u); mpz_clear(v); mpz_clear(bx); mpz_clear(by);
    mpz_clear(D); mpz_clear(N);
    return result;
}
