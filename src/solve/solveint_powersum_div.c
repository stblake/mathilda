/*
 * solveint_powersum_div.c
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


/* --- Power-sum divisor method for separable additive equations. ---
 *
 * For an equation  coef*(x^e + y^e) + (terms in the other variables) == 0
 * with e ODD, enumerate the "outer" variables over their box; for each, the
 * inner pair satisfies  x^e + y^e == m.  Because e is odd, s = x + y divides m,
 * and the power sum p_e(s, p) (p = x*y) is a degree-e/2 polynomial in p, so for
 * each divisor s of m the integer roots p give (x, y) via t^2 - s t + p = 0.
 * This turns the O(N) inner loop into O(#divisors(m)) -- the difference between
 * O(N^2) brute force and O(N * factoring) for e.g. sums of three cubes. */

/* Solve x^e + y^e == m (e odd) for the inner pair (ip, jp), outer variables
 * already fixed in vals.  Emits verified full assignments. */
void si_two_power_solve(SICtx* c, SearchState* st, int ip, int jp,
                               int e, const mpz_t m, int64_t* vals) {
    int64_t loi = c->lo[ip], hii = c->hi[ip], loj = c->lo[jp], hij = c->hi[jp];
    /* Tighten the pair windows by any abs-ordering that bounds a pair variable
     * below an already-assigned outer variable (|pair| < |vals[outer]|).  Sound
     * (prune-only; si_verify still backstops), and for the ordered three-cubes
     * query it shrinks the s = x+y range to ~2|outer| as the outer sweeps. */
    for (int k = 0; k < c->n_abs_ord; k++) {
        int a = c->abs_ord_a[k], b = c->abs_ord_b[k]; int s = c->abs_ord_strict[k] ? 1 : 0;
        if (b == ip || b == jp) continue;            /* larger side must be an outer */
        if (a != ip && a != jp) continue;            /* smaller side must be a pair var */
        int64_t vb = vals[b] < 0 ? -vals[b] : vals[b];
        int64_t bnd = vb - s; if (bnd < 0) bnd = -1; /* |a| <= bnd */
        if (a == ip) { if ( bnd < hii) hii =  bnd; if (-bnd > loi) loi = -bnd; }
        else         { if ( bnd < hij) hij =  bnd; if (-bnd > loj) loj = -bnd; }
    }
    if (mpz_sgn(m) == 0) {                        /* x^e + y^e = 0  ->  y = -x */
        int64_t t0 = (loi > -hij) ? loi : -hij;
        int64_t t1 = (hii < -loj) ? hii : -loj;
        for (int64_t t = t0; t <= t1 && !st->overflow; t++) {
            if (++st->visits > st->max_visits) { st->overflow = true; return; }
            vals[ip] = t; vals[jp] = -t; if (si_verify(c, vals)) emit_full(st, vals);
        }
        return;
    }
    int deg = e / 2;                              /* degree in p = x*y */
    if (deg > SI_LEAF_MAXDEG) return;
    mpz_t absm; mpz_init(absm); mpz_abs(absm, m);
    Expr* dl = divisors_ordinary(absm);
    mpz_clear(absm);
    if (!dl || dl->type != EXPR_FUNCTION) { if (dl) expr_free(dl); return; }

    /* Newton-recurrence scratch for p_k(s, p) as a polynomial in p:
     * p_0 = 2, p_1 = s, p_k = s*p_{k-1} - p*p_{k-2}. */
    mpz_t sz, pm2[SI_LEAF_MAXDEG + 1], pm1[SI_LEAF_MAXDEG + 1], cur[SI_LEAF_MAXDEG + 1];
    mpz_t a[SI_LEAF_MAXDEG + 1], sq, r, num, tmp;
    mpz_init(sz); mpz_init(sq); mpz_init(r); mpz_init(num); mpz_init(tmp);
    for (int k = 0; k <= SI_LEAF_MAXDEG; k++) { mpz_init(pm2[k]); mpz_init(pm1[k]); mpz_init(cur[k]); mpz_init(a[k]); }

    /* |s| = |x+y| is at most the box span. */
    int64_t smax = (hii - loj);
    int64_t smin = (loi - hij);
    int64_t sbound = (smax > -smin) ? smax : -smin;

    for (size_t di = 0; di < dl->data.function.arg_count && !st->overflow; di++) {
        int64_t dv;
        if (!expr_as_i64(dl->data.function.args[di], &dv)) continue;
        for (int sgn = 1; sgn >= -1; sgn -= 2) {
            int64_t sval = sgn * dv;
            if (sval < smin || sval > smax) continue;      /* s = x+y out of range */
            if (++st->visits > st->max_visits) { st->overflow = true; break; }
            mpz_set_si(sz, sval);
            /* Build p_e(s, p) via the recurrence (polynomials in p). */
            for (int k = 0; k <= deg; k++) { mpz_set_ui(pm2[k], 0); mpz_set_ui(pm1[k], 0); }
            mpz_set_ui(pm2[0], 2);                          /* p_0 = 2 */
            mpz_set(pm1[0], sz);                            /* p_1 = s */
            int dm2 = 0, dm1 = 0;
            for (int k = 2; k <= e; k++) {
                for (int t = 0; t <= deg; t++) mpz_set_ui(cur[t], 0);
                for (int t = 0; t <= dm1; t++) mpz_mul(cur[t], sz, pm1[t]);        /* s*p_{k-1} */
                for (int t = 0; t <= dm2; t++) mpz_sub(cur[t + 1], cur[t + 1], pm2[t]); /* - p*p_{k-2} */
                int dc = dm1; if (dm2 + 1 > dc) dc = dm2 + 1;
                for (int t = 0; t <= deg; t++) { mpz_set(pm2[t], pm1[t]); mpz_set(pm1[t], cur[t]); }
                dm2 = dm1; dm1 = dc;
            }
            for (int t = 0; t <= deg; t++) mpz_set(a[t], pm1[t]);
            mpz_sub(a[0], a[0], m);                         /* p_e(s,p) - m == 0 */

            int64_t proots[SI_LEAF_MAXDEG * 2];
            int npr = univariate_roots(a, deg, INT64_MIN / 2, INT64_MAX / 2, proots);
            for (int pi = 0; pi < npr; pi++) {
                /* x, y from t^2 - s t + p = 0: disc = s^2 - 4p. */
                mpz_mul(sq, sz, sz);
                mpz_set_si(tmp, proots[pi]); mpz_mul_ui(tmp, tmp, 4); mpz_sub(sq, sq, tmp);
                if (mpz_sgn(sq) < 0 || !mpz_perfect_square_p(sq)) continue;
                mpz_sqrt(r, sq);
                for (int rs = 1; rs >= -1; rs -= 2) {
                    mpz_set(num, sz); if (rs < 0) mpz_sub(num, num, r); else mpz_add(num, num, r);
                    if (!mpz_divisible_ui_p(num, 2)) continue;
                    mpz_divexact_ui(num, num, 2);           /* x = (s +/- r)/2 */
                    if (!mpz_fits_slong_p(num)) continue;
                    int64_t xv = mpz_get_si(num);
                    int64_t yv = sval - xv;
                    if (xv < loi || xv > hii || yv < loj || yv > hij) continue;
                    vals[ip] = xv; vals[jp] = yv;
                    if (si_verify(c, vals)) emit_full(st, vals);
                }
            }
            (void)sbound;
        }
    }
    for (int k = 0; k <= SI_LEAF_MAXDEG; k++) { mpz_clear(pm2[k]); mpz_clear(pm1[k]); mpz_clear(cur[k]); mpz_clear(a[k]); }
    mpz_clear(sz); mpz_clear(sq); mpz_clear(r); mpz_clear(num); mpz_clear(tmp);
    expr_free(dl);
}


/* Detect a separable additive equation with two variables sharing an odd
 * exponent, and solve it by the divisor method (outer variables enumerated,
 * inner pair solved per fixed outer assignment).  Returns true if handled. */
bool si_solve_powersum_divisor(SICtx* c, SearchState* st) {
    if (c->neq != 1) return false;
    const MPoly* eq = c->eq[0];
    int n = c->n;
    for (int i = 0; i < n; i++) if (!(c->has_lo[i] && c->has_hi[i])) return false;

    /* Separable: every term touches at most one variable.  Record each
     * variable's single power term (exponent, coefficient), if it is a lone
     * power c*v^e. */
    int vexp[SI_MAX_VARS]; mpz_t vcoef[SI_MAX_VARS];
    int nterm[SI_MAX_VARS];
    for (int i = 0; i < n; i++) { vexp[i] = -1; nterm[i] = 0; mpz_init_set_ui(vcoef[i], 0); }
    bool ok = true;
    for (size_t t = 0; t < eq->n_terms && ok; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        int nz = -1;
        for (int v = 0; v < n; v++) if (ex[v] > 0) { if (nz >= 0) ok = false; nz = v; }
        if (!ok) break;
        if (nz < 0) continue;                        /* constant term */
        nterm[nz]++; vexp[nz] = ex[nz]; mpz_set(vcoef[nz], eq->coefs[t]);
    }
    int ip = -1, jp = -1;
    if (ok) {
        for (int i = 0; i < n && ip < 0; i++) {
            if (nterm[i] != 1 || vexp[i] < 3 || (vexp[i] % 2) == 0) continue;
            for (int j = i + 1; j < n; j++) {
                if (nterm[j] == 1 && vexp[j] == vexp[i] && mpz_cmp(vcoef[j], vcoef[i]) == 0)
                    { ip = i; jp = j; break; }
            }
        }
    }

    bool handled = false;
    if (ip >= 0) {
        int e = vexp[ip];
        /* Outer variables enumeration size. */
        long double outer = 1.0L;
        for (int i = 0; i < n; i++) if (i != ip && i != jp)
            outer *= (long double)(c->hi[i] - c->lo[i] + 1);
        /* |m| = |x^e + y^e| ~ B^e is the number the method factors once per outer
         * assignment; keep it in the fast-factoring range so the O(N*factoring)
         * budget stays reasonable (this admits cubes over |.|<10^6, and higher
         * powers only over correspondingly smaller boxes). */
        int64_t Bin = 0;
        for (int q = 0; q < 2; q++) {
            int vv = (q == 0) ? ip : jp;
            int64_t b = (c->hi[vv] > -c->lo[vv]) ? c->hi[vv] : -c->lo[vv];
            if (b > Bin) Bin = b;
        }
        long double mmag = 1.0L;
        for (int k = 0; k < e; k++) mmag *= (long double)Bin;
        if (outer <= 3.0e5L && mmag <= 1.0e18L) {    /* factoring budget bounded */
            handled = true;
            st->max_visits = SI_MAX_NODES;
            /* Odometer over the outer variables. */
            int outv[SI_MAX_VARS], nouter = 0;
            for (int i = 0; i < n; i++) if (i != ip && i != jp) outv[nouter++] = i;
            int64_t vals[SI_MAX_VARS];
            for (int i = 0; i < n; i++) vals[i] = 0;
            for (int q = 0; q < nouter; q++) vals[outv[q]] = c->lo[outv[q]];
            mpz_t rest, mm, coef; mpz_init(rest); mpz_init(mm); mpz_init_set(coef, vcoef[ip]);
            for (;;) {
                vals[ip] = 0; vals[jp] = 0;
                si_eval_mpoly(eq, vals, rest);       /* = sum of other terms + const */
                mpz_neg(mm, rest);
                if (mpz_divisible_p(mm, coef)) {
                    mpz_divexact(mm, mm, coef);      /* x^e + y^e == mm */
                    si_two_power_solve(c, st, ip, jp, e, mm, vals);
                }
                int q = 0;
                for (; q < nouter; q++) {
                    if (++vals[outv[q]] <= c->hi[outv[q]]) break;
                    vals[outv[q]] = c->lo[outv[q]];
                }
                if (q == nouter || st->overflow) break;
            }
            mpz_clear(rest); mpz_clear(mm); mpz_clear(coef);
        }
    }
    for (int i = 0; i < n; i++) mpz_clear(vcoef[i]);
    return handled;
}
