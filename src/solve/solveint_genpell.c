/*
 * solveint_genpell.c
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


/* --- Generalised Pell  x^2 - D y^2 == N  (any nonzero N), unbounded positive
 *     orthant -> a parametric family PER solution class. ---
 *
 * Detect x^2 - D y^2 == N with a +1 square coefficient, D > 0 non-square, and
 * ANY nonzero N (si_pell_detect handles only N = +/-1). */
static bool si_genpell_detect(const MPoly* eq, int n, int* xside, int* yside,
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
    return (mpz_sgn(N) != 0);
}


/* Largest y a Nagell fundamental solution can have.  Ceiling over-estimate of
 * u*sqrt(|N|/(2(t +/- 1))) (+ for N>0, - for N<0), a safe bound validated
 * against brute force.  Sets *ok false if it would be impractically large. */
static void si_genpell_ybound(const mpz_t N, const mpz_t t, const mpz_t u,
                              mpz_t out, bool* ok) {
    mpz_t num, den, absN; mpz_init(num); mpz_init(den); mpz_init(absN);
    mpz_abs(absN, N);
    mpz_mul(num, u, u); mpz_mul(num, num, absN);             /* u^2 |N| */
    if (mpz_sgn(N) > 0) { mpz_add_ui(den, t, 1); }           /* 2(t+1) */
    else { mpz_sub_ui(den, t, 1); }                          /* 2(t-1) */
    mpz_mul_ui(den, den, 2);
    if (mpz_sgn(den) <= 0) { mpz_set_ui(out, 0); *ok = true; }
    else {
        mpz_fdiv_q(num, num, den);
        mpz_sqrt(out, num);
        mpz_add_ui(out, out, 2);                             /* safety margin */
    }
    /* Guard against a pathological search (huge fundamental unit and |N|). */
    *ok = (mpz_cmp_ui(out, 50000000UL) <= 0);
    mpz_clear(num); mpz_clear(den); mpz_clear(absN);
}


/* Compose (x, y) with the fundamental unit (t, u): (x,y) <- (x t + D y u,
 * x u + y t). */
static void si_pell_step(mpz_t x, mpz_t y, const mpz_t t, const mpz_t u,
                         const mpz_t D, mpz_t s1, mpz_t s2) {
    mpz_mul(s1, x, t); mpz_mul(s2, D, y); mpz_mul(s2, s2, u); mpz_add(s1, s1, s2); /* x t + D y u */
    mpz_mul(s2, x, u); mpz_mul(x, y, t); mpz_add(s2, s2, x);                       /* x u + y t */
    mpz_set(x, s1); mpz_set(y, s2);
}


/* Fundamental unit (t, u) of x^2 - D y^2 = 1, and the minimal positive-orthant
 * representative (x > 0, y > 0) of every solution class of x^2 - D y^2 == N.
 * The bases are written (deduplicated) into bx[]/by[]; *nb is their count.
 * Returns false to decline (no unit found, or the Nagell box too large). */
static bool si_genpell_bases(const mpz_t D, const mpz_t N, mpz_t t, mpz_t u,
                             mpz_t* bx, mpz_t* by, int* nb, int cap) {
    *nb = 0;
    mpz_t one, dummy1, dummy2;
    mpz_init_set_ui(one, 1); mpz_init(dummy1); mpz_init(dummy2);
    if (!si_pell_cf(D, one, t, u, dummy1, dummy2)) {  /* (t,u): x^2 - D y^2 = 1 */
        mpz_clear(one); mpz_clear(dummy1); mpz_clear(dummy2); return false;
    }
    mpz_clear(one); mpz_clear(dummy1); mpz_clear(dummy2);

    mpz_t Ymax; mpz_init(Ymax); bool okb = true;
    si_genpell_ybound(N, t, u, Ymax, &okb);
    if (!okb || !mpz_fits_slong_p(Ymax)) { mpz_clear(Ymax); return false; }
    int64_t ymax = mpz_get_si(Ymax);
    mpz_clear(Ymax);

    mpz_t val, root, x, y, s1, s2, px, py;
    mpz_inits(val, root, x, y, s1, s2, px, py, NULL);
    for (int64_t yy = 0; yy <= ymax && *nb < cap; yy++) {
        mpz_set_si(y, yy); mpz_mul(val, y, y); mpz_mul(val, val, D);
        mpz_add(val, val, N);                                /* N + D y^2 */
        if (mpz_sgn(val) <= 0) continue;
        if (!mpz_perfect_square_p(val)) continue;
        mpz_sqrt(root, val);                                 /* x0 > 0 */
        /* All four sign combos (+/-x0, +/-yy) are candidate fundamentals of
         * (possibly distinct) classes; advance each into the positive orthant. */
        for (int sx = 1; sx >= -1; sx -= 2)
            for (int sy = 1; sy >= -1; sy -= 2) {
                mpz_set(x, root); if (sx < 0) mpz_neg(x, x);
                mpz_set_si(s1, sy * yy); mpz_set(y, s1);
                int guard = 0;
                while (!(mpz_sgn(x) > 0 && mpz_sgn(y) > 0) && guard < 1000) {
                    si_pell_step(x, y, t, u, D, s1, s2); guard++;
                }
                if (!(mpz_sgn(x) > 0 && mpz_sgn(y) > 0)) continue;
                /* Reduce to the MINIMAL positive-orthant representative: apply
                 * eps^{-1} = (t, -u) while the result stays in the orthant, so
                 * two solutions of the same class (e.g. (1,1) and (5,3)=(1,1)eps
                 * for D=3, N=-2) collapse to one base rather than two overlapping
                 * families.  eps^{-1}: x' = x t - D y u,  y' = -x u + y t. */
                for (int g2 = 0; g2 < 1000; g2++) {
                    mpz_mul(px, x, t); mpz_mul(s2, D, y); mpz_mul(s2, s2, u);
                    mpz_sub(px, px, s2);                       /* x' = x t - D y u */
                    mpz_mul(py, y, t); mpz_mul(s2, x, u);
                    mpz_sub(py, py, s2);                       /* y' = y t - x u */
                    if (mpz_sgn(px) > 0 && mpz_sgn(py) > 0) { mpz_set(x, px); mpz_set(y, py); }
                    else break;
                }
                bool dup = false;
                for (int i = 0; i < *nb; i++)
                    if (mpz_cmp(bx[i], x) == 0 && mpz_cmp(by[i], y) == 0) { dup = true; break; }
                if (dup || *nb >= cap) continue;
                mpz_set(bx[*nb], x); mpz_set(by[*nb], y); (*nb)++;
            }
    }
    mpz_clears(val, root, x, y, s1, s2, px, py, NULL);
    return true;
}


/* Order bases ascending by (x, then y) for a deterministic emission order. */
static void si_genpell_sort_bases(mpz_t* bx, mpz_t* by, int nb) {
    for (int i = 0; i < nb; i++)
        for (int j = i + 1; j < nb; j++) {
            int c = mpz_cmp(bx[i], bx[j]);
            if (c > 0 || (c == 0 && mpz_cmp(by[i], by[j]) > 0)) {
                mpz_swap(bx[i], bx[j]); mpz_swap(by[i], by[j]);
            }
        }
}


/* Build the class family tuple for base (a, b), fundamental unit (t, u):
 *   x -> (P eps^C + Pbar epsbar^C) / 2,
 *   y -> (P eps^C - Pbar epsbar^C) / (2 Sqrt[D]),
 * where P = a + b Sqrt[D], eps = t + u Sqrt[D] (bars negate the Sqrt[D] part),
 * each a ConditionalExpression on C[1] >= 0.  Rules placed at xs/ys. */
static Expr* si_genpell_family(SICtx* c, int xs, int ys, const mpz_t a,
                               const mpz_t b, const mpz_t t, const mpz_t u,
                               const mpz_t D) {
    Expr* sqrtD = mk_fn1("Sqrt", mk_mpz(D));
    mpz_t nb, nu; mpz_init(nb); mpz_init(nu); mpz_neg(nb, b); mpz_neg(nu, u);
    Expr* P    = mk_fn2("Plus", mk_mpz(a), mk_fn2("Times", mk_mpz(b),  expr_copy(sqrtD)));
    Expr* Pbar = mk_fn2("Plus", mk_mpz(a), mk_fn2("Times", mk_mpz(nb), expr_copy(sqrtD)));
    Expr* eps    = mk_fn2("Plus", mk_mpz(t), mk_fn2("Times", mk_mpz(u),  expr_copy(sqrtD)));
    Expr* epsbar = mk_fn2("Plus", mk_mpz(t), mk_fn2("Times", mk_mpz(nu), expr_copy(sqrtD)));
    mpz_clear(nb); mpz_clear(nu);
    Expr* epsC    = mk_fn2("Power", eps,    mk_fn1("C", mk_int(1)));
    Expr* epsbarC = mk_fn2("Power", epsbar, mk_fn1("C", mk_int(1)));
    Expr* Pe    = mk_fn2("Times", P, epsC);
    Expr* Pbare = mk_fn2("Times", Pbar, epsbarC);
    /* xval = (Pe + Pbare)/2 */
    Expr* xval = mk_fn2("Times",
        mk_fn2("Plus", expr_copy(Pe), expr_copy(Pbare)),
        mk_fn2("Power", mk_int(2), mk_int(-1)));
    /* yval = (Pe - Pbare)/(2 Sqrt[D]) */
    Expr* yval = mk_fn2("Times",
        mk_fn2("Plus", Pe, mk_fn2("Times", mk_int(-1), Pbare)),
        mk_fn2("Power", mk_fn2("Times", mk_int(2), expr_copy(sqrtD)), mk_int(-1)));
    expr_free(sqrtD);

    Expr* cond = mk_fn2("GreaterEqual", mk_fn1("C", mk_int(1)), mk_int(0));
    Expr* cex = mk_fn2("ConditionalExpression", xval, expr_copy(cond));
    Expr* cey = mk_fn2("ConditionalExpression", yval, cond);

    Expr* rules[SI_MAX_VARS];
    for (int i = 0; i < c->n; i++) rules[i] = NULL;
    rules[xs] = mk_rule(expr_copy(c->var[xs]), cex);
    rules[ys] = mk_rule(expr_copy(c->var[ys]), cey);
    Expr* rlist[SI_MAX_VARS]; int nr = 0;
    for (int i = 0; i < c->n; i++) if (rules[i]) rlist[nr++] = rules[i];
    return mk_list(rlist, (size_t)nr);
}


/* Unbounded x^2 - D y^2 == N (N != +1: any N < 0 or N >= 2) with x > 0 && y > 0
 * -> one parametric family per solution class.  Returns the owned family List
 * (empty {} when the equation is provably unsolvable, including negative Pell
 * over a D whose sqrt has even CF period), or NULL to decline (not this shape,
 * bounded, N = +1, or the Nagell search too large -- handled elsewhere). */
Expr* si_solve_genpell_parametric(SICtx* c) {
    if (c->neq != 1) return NULL;
    int xs, ys; mpz_t D, N; mpz_init(D); mpz_init(N);
    if (!si_genpell_detect(c->eq[0], c->n, &xs, &ys, D, N)
        || mpz_sgn(N) == 0 || mpz_cmp_si(N, 1) == 0) {  /* N = +1 -> pell_parametric */
        mpz_clear(D); mpz_clear(N); return NULL;
    }
    bool unbounded = !c->has_hi[xs] && !c->has_hi[ys];
    bool posx = c->has_lo[xs] && c->lo[xs] >= 1;
    bool posy = c->has_lo[ys] && c->lo[ys] >= 1;
    if (!(unbounded && posx && posy) || c->n_ord != 0 || c->n_neq != 0 || c->n_abs_ord != 0 || !c->all_captured) {
        mpz_clear(D); mpz_clear(N); return NULL;
    }

    enum { GP_CAP = 4096 };
    mpz_t t, u; mpz_init(t); mpz_init(u);
    mpz_t *bx = (mpz_t*)malloc(sizeof(mpz_t) * GP_CAP);
    mpz_t *by = (mpz_t*)malloc(sizeof(mpz_t) * GP_CAP);
    for (int i = 0; i < GP_CAP; i++) { mpz_init(bx[i]); mpz_init(by[i]); }
    int nb = 0;
    bool solved = si_genpell_bases(D, N, t, u, bx, by, &nb, GP_CAP);

    Expr* result = NULL;
    if (solved) {
        si_genpell_sort_bases(bx, by, nb);
        if (nb == 0) {
            result = mk_list(NULL, 0);            /* provably unsolvable -> {} */
        } else {
            Expr** tuples = (Expr**)malloc(sizeof(Expr*) * (size_t)nb);
            for (int i = 0; i < nb; i++)
                tuples[i] = si_genpell_family(c, xs, ys, bx[i], by[i], t, u, D);
            result = eval_and_free(mk_list(tuples, (size_t)nb));
            free(tuples);
        }
    }

    for (int i = 0; i < GP_CAP; i++) { mpz_clear(bx[i]); mpz_clear(by[i]); }
    free(bx); free(by);
    mpz_clear(t); mpz_clear(u); mpz_clear(D); mpz_clear(N);
    return result;                               /* NULL -> decline (unevaluated) */
}
