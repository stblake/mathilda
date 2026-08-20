/*
 * solvethue.c -- see solvethue.h.
 *
 * Scope of the current implementation: monic forms (|a_0| == 1) with
 * |m| == 1, over a monogenic field with a real embedding (torsion {+-1}).
 * Then beta = x - theta*y is a UNIT, so every solution is
 *     x - theta*y = +- prod_k eps_k^{b_k}
 * for the certified fundamental units eps_k.  Enumerate the exponent
 * vectors (b_k) up to a bound, reconstruct (x, y) from the coordinates of
 * the resulting unit (it must lie in Z + Z*theta), and verify F(x,y)==m
 * exactly.  The bound is the one piece that makes this COMPLETE: the
 * rigorous Baker/de-Weger bound (thue_exponent_bound) currently declines,
 * so the Solve-facing entry declines until it lands; the *_bounded entry
 * validates the reconstruction mechanism with an explicit bound.
 */
#include "solvethue.h"

#include "symtab.h"
#include "eval.h"
#include "sym_names.h"

#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef USE_FLINT
#include "numberfield.h"
#include "numberfield_internal.h"
#include "nfunits.h"
#include "linalg.h"                 /* lll_reduce_q */
#include <math.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
#include <flint/fmpz_vec.h>
#include <flint/nf_elem.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arf.h>

/* ---- solution list ---- */
void thue_sols_free(ThueSol* sols, int count) {
    if (!sols) return;
    for (int i = 0; i < count; i++) { mpz_clear(sols[i].x); mpz_clear(sols[i].y); }
    free(sols);
}

/* ---- nf_elem <-> integer coordinates ---- */

static void nfelem_from_coords(nf_t nf, const mpz_t* coords, int deg, nf_elem_t out) {
    fmpq_poly_t g; fmpq_poly_init(g);
    fmpz_t z; fmpz_init(z);
    for (int i = 0; i < deg; i++) { fmpz_set_mpz(z, coords[i]); fmpq_poly_set_coeff_fmpz(g, i, z); }
    fmpz_clear(z);
    nf_elem_set_fmpq_poly(out, g, nf);
    fmpq_poly_clear(g);
}

/* out = (1/D) * sum coords[i] theta^i  -- an O_K element from its theta-power
 * NUMERATOR and the shared denominator D (D = 1 recovers nfelem_from_coords). */
static void nfelem_from_coords_den(nf_t nf, const mpz_t* coords, const mpz_t D, int deg, nf_elem_t out) {
    fmpq_poly_t g; fmpq_poly_init(g);
    fmpz_t z, d; fmpz_init(z); fmpz_init(d); fmpz_set_mpz(d, D);
    for (int i = 0; i < deg; i++) { fmpz_set_mpz(z, coords[i]); fmpq_poly_set_coeff_fmpz(g, i, z); }
    if (!fmpz_is_one(d)) fmpq_poly_scalar_div_fmpz(g, g, d);
    fmpz_clear(z); fmpz_clear(d);
    nf_elem_set_fmpq_poly(out, g, nf);
    fmpq_poly_clear(g);
}

/* Extract integer coords of e into out[0..deg-1]; false if any non-integer. */
static bool coords_from_nfelem(nf_t nf, nf_elem_t e, int deg, mpz_t* out) {
    fmpq_poly_t g; fmpq_poly_init(g);
    nf_elem_get_fmpq_poly(g, e, nf);
    bool ok = true;
    fmpq_t c; fmpq_init(c);
    for (int i = 0; i < deg; i++) {
        fmpq_poly_get_coeff_fmpq(c, g, i);
        if (!fmpz_is_one(fmpq_denref(c))) { ok = false; break; }
        fmpz_get_mpz(out[i], fmpq_numref(c));
    }
    fmpq_clear(c);
    fmpq_poly_clear(g);
    return ok;
}

/* F(x,y) = sum_j form[j] x^(n-j) y^j, into `out`. */
static void eval_form(const mpz_t* form, int n, const mpz_t x, const mpz_t y, mpz_t out) {
    mpz_set_ui(out, 0);
    mpz_t xp, yp, term; mpz_init(xp); mpz_init(yp); mpz_init(term);
    for (int j = 0; j <= n; j++) {
        mpz_pow_ui(xp, x, (unsigned long)(n - j));
        mpz_pow_ui(yp, y, (unsigned long)j);
        mpz_mul(term, form[j], xp);
        mpz_mul(term, term, yp);
        mpz_add(out, out, term);
    }
    mpz_clear(xp); mpz_clear(yp); mpz_clear(term);
}

/* ================================================================== *
 *  Reducible binary forms  F(x,y) == m  (F reducible over Q).          *
 *                                                                      *
 *  Not a Thue equation (the Tzanakis-de Weger method needs F           *
 *  irreducible), but finite when F factors into >= 2 COPRIME           *
 *  homogeneous forms: at an integer point each factor value            *
 *  G_i(x,y) is an integer and prod G_i(x,y)^{e_i} = m / content.  We   *
 *  enumerate the (finitely many) signed factorisations and solve each  *
 *  system { G_i(x,y) = d_i }: parametrise the line of one LINEAR       *
 *  factor, substitute into a second factor to get a univariate whose   *
 *  integer roots are the candidate points, and verify every factor     *
 *  equation and F == m exactly.  A pure power of a single irreducible  *
 *  factor (e.g. (x-y)^3 == 1) has infinitely many or PARI-refused      *
 *  solutions -> DECLINE.  Returns the complete finite count, or -1.    */

/* True iff g factors over Q into more than one irreducible (or a repeated one):
 * FLINT 3 has no fmpz_poly_is_irreducible, so read it off the factorisation. */
static bool solvethue_poly_reducible(const fmpz_poly_t g) {
    fmpz_poly_factor_t fac; fmpz_poly_factor_init(fac); fmpz_poly_factor(fac, g);
    bool red = !(fac->num == 1 && fac->exp[0] == 1
                 && fmpz_poly_degree(&fac->p[0]) == fmpz_poly_degree(g));
    fmpz_poly_factor_clear(fac);
    return red;
}

/* value of a homogeneous binary factor G (coeffs G[0..dG], G[p]=coeff of
 * x^(dG-p) y^p) at integer (x,y). */
static void eval_binform(const mpz_t* G, int dG, const mpz_t x, const mpz_t y, mpz_t out) {
    mpz_set_ui(out, 0);
    mpz_t xp, yp, term; mpz_init(xp); mpz_init(yp); mpz_init(term);
    for (int p = 0; p <= dG; p++) {
        mpz_pow_ui(xp, x, (unsigned long)(dG - p));
        mpz_pow_ui(yp, y, (unsigned long)p);
        mpz_mul(term, G[p], xp); mpz_mul(term, term, yp); mpz_add(out, out, term);
    }
    mpz_clear(xp); mpz_clear(yp); mpz_clear(term);
}

typedef struct { mpz_t* G; int deg; int mult; } RedFactor;   /* homogeneous factor */

/* push (x,y) into a dynamic mpz pair buffer with dedup */
static void redpush(mpz_t** xs, mpz_t** ys, int* n, int* cap, const mpz_t x, const mpz_t y) {
    for (int i = 0; i < *n; i++) if (mpz_cmp((*xs)[i], x) == 0 && mpz_cmp((*ys)[i], y) == 0) return;
    if (*n == *cap) { *cap = *cap ? *cap * 2 : 16; *xs = realloc(*xs, sizeof(mpz_t) * (size_t)*cap); *ys = realloc(*ys, sizeof(mpz_t) * (size_t)*cap); }
    mpz_init_set((*xs)[*n], x); mpz_init_set((*ys)[*n], y); (*n)++;
}

/* integer roots of an fmpz_poly P (via its factorisation): each linear factor
 * (a t + b) with a | b contributes t = -b/a.  Appends to roots[] (mpz). */
static void fmpz_poly_int_roots(const fmpz_poly_t P, mpz_t** roots, int* nr, int* cap) {
    if (fmpz_poly_degree(P) < 1) return;
    fmpz_poly_factor_t fac; fmpz_poly_factor_init(fac); fmpz_poly_factor(fac, P);
    for (slong i = 0; i < fac->num; i++) {
        if (fmpz_poly_degree(&fac->p[i]) != 1) continue;
        fmpz_t a, b; fmpz_init(a); fmpz_init(b);
        fmpz_poly_get_coeff_fmpz(a, &fac->p[i], 1);
        fmpz_poly_get_coeff_fmpz(b, &fac->p[i], 0);
        if (!fmpz_is_zero(a) && fmpz_divisible(b, a)) {
            fmpz_t r; fmpz_init(r); fmpz_neg(r, b); fmpz_divexact(r, r, a);
            if (*nr == *cap) { *cap = *cap ? *cap * 2 : 8; *roots = realloc(*roots, sizeof(mpz_t) * (size_t)*cap); }
            mpz_init(( *roots)[*nr]); fmpz_get_mpz((*roots)[*nr], r); (*nr)++;
            fmpz_clear(r);
        }
        fmpz_clear(a); fmpz_clear(b);
    }
    fmpz_poly_factor_clear(fac);
}

/* Build the univariate (in x) poly  F(x, y) - m  for a fixed integer y:
 * F(x,y) = sum_j form[j] x^(n-j) y^j, so coeff of x^(n-j) is form[j]*y^j
 * (the powers n-j are all distinct), minus m at x^0.  Its integer roots are
 * exactly the integer x with F(x, y) == m -- used to close the small-|Y| gap
 * exactly, in place of an O(Xmax) scan over x. */
static void thue_form_xpoly(const mpz_t* form, int n, const mpz_t y, const mpz_t m, fmpz_poly_t out) {
    fmpz_poly_zero(out);
    fmpz_t z, yp; fmpz_init(z); fmpz_init(yp);
    for (int j = 0; j <= n; j++) {
        fmpz_set_mpz(yp, y); fmpz_pow_ui(yp, yp, (unsigned long)j);
        fmpz_set_mpz(z, form[j]); fmpz_mul(z, z, yp);
        fmpz_poly_set_coeff_fmpz(out, n - j, z);
    }
    fmpz_t c0; fmpz_init(c0); fmpz_poly_get_coeff_fmpz(c0, out, 0);
    fmpz_set_mpz(z, m); fmpz_sub(c0, c0, z); fmpz_poly_set_coeff_fmpz(out, 0, c0);
    fmpz_clear(c0); fmpz_clear(z); fmpz_clear(yp);
}

/* Build the univariate (in y) poly  G(X, y) - target  for a fixed integer X:
 * coeff of y^p is G[p]*X^(dG-p), minus target at y^0. */
static void red_yspec(fmpz_poly_t out, const mpz_t* G, int dG, const mpz_t X, const mpz_t target) {
    fmpz_poly_zero(out);
    fmpz_t z, xp; fmpz_init(z); fmpz_init(xp);
    for (int p = 0; p <= dG; p++) {
        fmpz_set_mpz(xp, X); fmpz_pow_ui(xp, xp, (unsigned long)(dG - p));
        fmpz_set_mpz(z, G[p]); fmpz_mul(z, z, xp);
        fmpz_poly_set_coeff_fmpz(out, p, z);
    }
    fmpz_t c0; fmpz_init(c0); fmpz_poly_get_coeff_fmpz(c0, out, 0);
    fmpz_set_mpz(z, target); fmpz_sub(c0, c0, z); fmpz_poly_set_coeff_fmpz(out, 0, c0);
    fmpz_clear(c0); fmpz_clear(z); fmpz_clear(xp);
}

/* No linear factor: eliminate y between factors 0 and 1 via the resultant
 * R(x) = Res_y(G_0(x,y)-d_0, G_1(x,y)-d_1) (obtained by sampling + integer
 * interpolation), whose integer roots are the candidate x; solve G_0(x,y)=d_0
 * for y and verify every factor and F == m. */
static void red_solve_noline(const RedFactor* fac, int k, const mpz_t* d,
                             const mpz_t* form, int n, const mpz_t m,
                             mpz_t** xs, mpz_t** ys, int* np, int* cap) {
    int da = fac[0].deg, db = fac[1].deg;
    if (da > 8 || db > 8) return;                        /* keep the sampling small */
    int nb = 2 * da * db + 4;                            /* > deg R(x): over-sample is safe */
    fmpz *Xs = _fmpz_vec_init(nb), *Ys = _fmpz_vec_init(nb);
    fmpz_poly_t Pa, Pb, R; fmpz_poly_init(Pa); fmpz_poly_init(Pb); fmpz_poly_init(R);
    fmpz_t rr; fmpz_init(rr);
    mpz_t Xi; mpz_init(Xi);
    for (int s = 0; s < nb; s++) {
        fmpz_set_si(Xs + s, s); mpz_set_si(Xi, s);
        red_yspec(Pa, fac[0].G, da, Xi, d[0]);
        red_yspec(Pb, fac[1].G, db, Xi, d[1]);
        fmpz_poly_resultant(rr, Pa, Pb);
        fmpz_set(Ys + s, rr);
    }
    fmpz_poly_interpolate_fmpz_vec(R, Xs, Ys, nb);

    mpz_t* xroots = NULL; int nxr = 0, xcap = 0;
    fmpz_poly_int_roots(R, &xroots, &nxr, &xcap);
    mpz_t X, Y, val; mpz_init(X); mpz_init(Y); mpz_init(val);
    for (int xi = 0; xi < nxr; xi++) {
        mpz_set(X, xroots[xi]);
        red_yspec(Pa, fac[0].G, da, X, d[0]);            /* G_0(X,y) = d_0 in y */
        mpz_t* yroots = NULL; int nyr = 0, ycap = 0;
        fmpz_poly_int_roots(Pa, &yroots, &nyr, &ycap);
        for (int yi = 0; yi < nyr; yi++) {
            mpz_set(Y, yroots[yi]);
            bool good = true;
            for (int i = 0; i < k && good; i++) { eval_binform(fac[i].G, fac[i].deg, X, Y, val); if (mpz_cmp(val, d[i]) != 0) good = false; }
            if (good) { eval_form(form, n, X, Y, val); if (mpz_cmp(val, m) == 0) redpush(xs, ys, np, cap, X, Y); }
        }
        for (int yi = 0; yi < nyr; yi++) mpz_clear(yroots[yi]);
        free(yroots);
    }
    for (int xi = 0; xi < nxr; xi++) mpz_clear(xroots[xi]);
    free(xroots);
    mpz_clear(X); mpz_clear(Y); mpz_clear(val); mpz_clear(Xi); fmpz_clear(rr);
    fmpz_poly_clear(Pa); fmpz_poly_clear(Pb); fmpz_poly_clear(R);
    _fmpz_vec_clear(Xs, nb); _fmpz_vec_clear(Ys, nb);
}

/* Solve the system { fac[i].G(x,y) = d[i] } given signed targets d[], collecting
 * verified points (that also satisfy F == m) into (xs,ys). */
static void red_solve_system(const RedFactor* fac, int k, const mpz_t* d,
                             const mpz_t* form, int n, const mpz_t m,
                             mpz_t** xs, mpz_t** ys, int* np, int* cap) {
    /* pick a LINEAR factor as the line, and any other factor as the pivot. */
    int li = -1;
    for (int i = 0; i < k; i++) if (fac[i].deg == 1) { li = i; break; }
    if (li < 0) { red_solve_noline(fac, k, d, form, n, m, xs, ys, np, cap); return; }
    int pv = -1;
    for (int i = 0; i < k; i++) if (i != li) { pv = i; break; }
    if (pv < 0) return;

    /* line  a x + b y = d[li],  a = G[0], b = G[1] (degree-1 form a x + b y). */
    mpz_t a, b, g, x0, y0, dx, dy, s, t2;
    mpz_init_set(a, fac[li].G[0]); mpz_init_set(b, fac[li].G[1]);
    mpz_init(g); mpz_init(x0); mpz_init(y0); mpz_init(dx); mpz_init(dy); mpz_init(s); mpz_init(t2);
    mpz_gcdext(g, x0, y0, a, b);                          /* a*x0 + b*y0 = g */
    if (mpz_sgn(g) != 0 && mpz_divisible_p(d[li], g)) {
        mpz_divexact(s, d[li], g);
        mpz_mul(x0, x0, s); mpz_mul(y0, y0, s);          /* particular a x0 + b y0 = d[li] */
        mpz_divexact(dx, b, g); mpz_divexact(dy, a, g); mpz_neg(dy, dy);   /* direction (b/g, -a/g) */

        /* pivot P(t) = G_pv(x0+dx t, y0+dy t) - d[pv], an fmpz_poly in t. */
        fmpz_poly_t px, py, xt, yt, acc, Pt; fmpz_poly_init(px); fmpz_poly_init(py);
        fmpz_poly_init(xt); fmpz_poly_init(yt); fmpz_poly_init(acc); fmpz_poly_init(Pt);
        fmpz_t z; fmpz_init(z);
        fmpz_set_mpz(z, x0); fmpz_poly_set_coeff_fmpz(px, 0, z); fmpz_set_mpz(z, dx); fmpz_poly_set_coeff_fmpz(px, 1, z);
        fmpz_set_mpz(z, y0); fmpz_poly_set_coeff_fmpz(py, 0, z); fmpz_set_mpz(z, dy); fmpz_poly_set_coeff_fmpz(py, 1, z);
        int dG = fac[pv].deg;
        for (int p = 0; p <= dG; p++) {                  /* acc += G[p] * px^(dG-p) * py^p */
            fmpz_poly_one(xt); for (int e = 0; e < dG - p; e++) fmpz_poly_mul(xt, xt, px);
            fmpz_poly_one(yt); for (int e = 0; e < p; e++) fmpz_poly_mul(yt, yt, py);
            fmpz_poly_mul(xt, xt, yt);
            fmpz_set_mpz(z, fac[pv].G[p]); fmpz_poly_scalar_mul_fmpz(xt, xt, z);
            fmpz_poly_add(acc, acc, xt);
        }
        fmpz_set_mpz(z, d[pv]); fmpz_poly_set(Pt, acc);
        { fmpz_t c0; fmpz_init(c0); fmpz_poly_get_coeff_fmpz(c0, Pt, 0); fmpz_sub(c0, c0, z); fmpz_poly_set_coeff_fmpz(Pt, 0, c0); fmpz_clear(c0); }

        mpz_t* roots = NULL; int nr = 0, rcap = 0;
        fmpz_poly_int_roots(Pt, &roots, &nr, &rcap);
        mpz_t X, Y, val; mpz_init(X); mpz_init(Y); mpz_init(val);
        for (int ri = 0; ri < nr; ri++) {
            mpz_mul(X, dx, roots[ri]); mpz_add(X, X, x0);
            mpz_mul(Y, dy, roots[ri]); mpz_add(Y, Y, y0);
            bool good = true;
            for (int i = 0; i < k && good; i++) { eval_binform(fac[i].G, fac[i].deg, X, Y, val); if (mpz_cmp(val, d[i]) != 0) good = false; }
            if (good) { eval_form(form, n, X, Y, val); if (mpz_cmp(val, m) == 0) redpush(xs, ys, np, cap, X, Y); }
        }
        for (int ri = 0; ri < nr; ri++) mpz_clear(roots[ri]);
        free(roots);
        mpz_clear(X); mpz_clear(Y); mpz_clear(val);
        fmpz_poly_clear(px); fmpz_poly_clear(py); fmpz_poly_clear(xt); fmpz_poly_clear(yt);
        fmpz_poly_clear(acc); fmpz_poly_clear(Pt); fmpz_clear(z);
    }
    mpz_clear(a); mpz_clear(b); mpz_clear(g); mpz_clear(x0); mpz_clear(y0);
    mpz_clear(dx); mpz_clear(dy); mpz_clear(s); mpz_clear(t2);
}

/* recursively assign signed targets d[i] (i>=idx) with prod d[i]^{mult_i} =
 * remaining, then solve the resulting system. */
static void red_enumerate(const RedFactor* fac, int k, int idx, const mpz_t remaining,
                          mpz_t* d, const mpz_t* form, int n, const mpz_t m,
                          mpz_t** xs, mpz_t** ys, int* np, int* cap) {
    if (idx == k - 1) {
        /* last factor: d^{mult} must equal `remaining` exactly */
        int e = fac[idx].mult;
        if (e == 1) { mpz_set(d[idx], remaining); red_solve_system(fac, k, (const mpz_t*)d, form, n, m, xs, ys, np, cap); return; }
        /* e-th root of |remaining| (mpz_root on a negative with even e is
         * undefined -> SIGFPE; take |.| and check both signs against remaining). */
        mpz_t r, ar; mpz_init(r); mpz_init(ar); mpz_abs(ar, remaining);
        if (mpz_root(r, ar, (unsigned long)e)) {
            for (int sgn = 0; sgn < 2; sgn++) {
                mpz_set(d[idx], r); if (sgn) mpz_neg(d[idx], d[idx]);
                mpz_t chk; mpz_init(chk); mpz_pow_ui(chk, d[idx], (unsigned long)e);
                bool hit = (mpz_cmp(chk, remaining) == 0);
                mpz_clear(chk);
                if (hit) red_solve_system(fac, k, (const mpz_t*)d, form, n, m, xs, ys, np, cap);
                if (mpz_sgn(r) == 0) break;                     /* r == 0: +0 == -0 */
            }
        }
        mpz_clear(r); mpz_clear(ar);
        return;
    }
    /* factor idx: d ranges over signed divisors whose e-th power divides `remaining` */
    int e = fac[idx].mult;
    mpz_t absrem; mpz_init(absrem); mpz_abs(absrem, remaining);
    unsigned long lim = mpz_fits_ulong_p(absrem) ? mpz_get_ui(absrem) : 0;
    if (lim == 0 || lim > 2000000UL) { mpz_clear(absrem); return; }         /* too many divisors -> give up branch */
    mpz_t dd, de, rem2; mpz_init(dd); mpz_init(de); mpz_init(rem2);
    for (unsigned long v = 1; v <= lim; v++) {
        if (lim % v != 0) continue;
        mpz_set_ui(dd, v); mpz_pow_ui(de, dd, (unsigned long)e);
        if (!mpz_divisible_p(remaining, de)) continue;
        for (int sgn = 0; sgn < 2; sgn++) {
            mpz_set_ui(dd, v); if (sgn) mpz_neg(dd, dd);
            mpz_pow_ui(de, dd, (unsigned long)e);
            mpz_divexact(rem2, remaining, de);
            mpz_set(d[idx], dd);
            red_enumerate(fac, k, idx + 1, rem2, d, form, n, m, xs, ys, np, cap);
            if (e % 2 == 1 && v == 0) break;
        }
    }
    mpz_clear(dd); mpz_clear(de); mpz_clear(rem2); mpz_clear(absrem);
}

static int thue_solve_reducible_form(const mpz_t* form, int n, const mpz_t m, ThueSol** out) {
    /* F(x,1) = g(x): g_coeff[i] = form[n-i]. */
    fmpz_poly_t g; fmpz_poly_init(g);
    { fmpz_t z; fmpz_init(z); for (int i = 0; i <= n; i++) { fmpz_set_mpz(z, form[n - i]); fmpz_poly_set_coeff_fmpz(g, i, z); } fmpz_clear(z); }

    fmpz_poly_factor_t fac; fmpz_poly_factor_init(fac); fmpz_poly_factor(fac, g);
    int rc = -1;

    /* distinct factors -> homogeneous binary forms */
    int k = fac->num;
    if (k < 2) { fmpz_poly_factor_clear(fac); fmpz_poly_clear(g); return -1; }  /* single irreducible -> Thue/decline */

    RedFactor* F = malloc(sizeof(RedFactor) * (size_t)k);
    bool have_linear = false;
    for (int i = 0; i < k; i++) {
        int dG = (int)fmpz_poly_degree(&fac->p[i]);
        F[i].deg = dG; F[i].mult = (int)fac->exp[i];
        F[i].G = malloc(sizeof(mpz_t) * (size_t)(dG + 1));
        /* homogenise: G[p] = coeff of x^(dG-p) y^p = g_i coeff of x^(dG-p). */
        for (int p = 0; p <= dG; p++) { mpz_init(F[i].G[p]); fmpz_t z; fmpz_init(z); fmpz_poly_get_coeff_fmpz(z, &fac->p[i], dG - p); fmpz_get_mpz(F[i].G[p], z); fmpz_clear(z); }
        if (dG == 1) have_linear = true;
    }

    /* content c: F(x,1) = c * prod p_i^{e_i} (p_i primitive), stored in fac->c. */
    mpz_t c; mpz_init(c); fmpz_get_mpz(c, &fac->c);

    mpz_t M; mpz_init(M);
    (void)have_linear;   /* linear-factor path preferred; resultant handles the rest */
    bool ok = mpz_sgn(c) != 0 && mpz_divisible_p(m, c);
    mpz_t* xs = NULL; mpz_t* ys = NULL; int np = 0, cap = 0;
    if (ok) {
        mpz_divexact(M, m, c);
        mpz_t* d = malloc(sizeof(mpz_t) * (size_t)k);
        for (int i = 0; i < k; i++) mpz_init(d[i]);
        red_enumerate((const RedFactor*)F, k, 0, M, d, form, n, m, &xs, &ys, &np, &cap);
        for (int i = 0; i < k; i++) mpz_clear(d[i]);
        free(d);
        /* package (np may be 0 -> proven empty set) */
        ThueSol* res = malloc(sizeof(ThueSol) * (size_t)(np > 0 ? np : 1));
        for (int i = 0; i < np; i++) { mpz_init_set(res[i].x, xs[i]); mpz_init_set(res[i].y, ys[i]); }
        *out = res; rc = np;
    }

    for (int i = 0; i < np; i++) { mpz_clear(xs[i]); mpz_clear(ys[i]); }
    free(xs); free(ys);
    for (int i = 0; i < k; i++) { for (int p = 0; p <= F[i].deg; p++) mpz_clear(F[i].G[p]); free(F[i].G); }
    free(F);
    mpz_clear(c); mpz_clear(M);
    fmpz_poly_factor_clear(fac); fmpz_poly_clear(g);
    return rc;
}

/* Sentinel bound: compute the rigorous Baker/de-Weger bound internally. */
#define THUE_BOUND_RIGOROUS (-2)

/* ================================================================== *
 *  Rigorous exponent bound: Tzanakis-de Weger (Baker + LLL).          *
 *  Implements J. Number Theory 31 (1989) 99-132 (see                  *
 *  docs/references/thue/ALGORITHM_NOTES.md).                          *
 * ================================================================== */

typedef struct { bool ok; long B; long Y2p; long Xmax; } ThueBound;

/* Evaluate the algebraic integer with integer coords c[0..deg-1] at the
 * complex root z (acb Horner). */
static void embed_acb(const mpz_t* c, int deg, const acb_t z, acb_t out, slong prec) {
    acb_zero(out);
    fmpz_t ci; fmpz_init(ci);
    for (int i = deg - 1; i >= 0; i--) {
        acb_mul(out, out, z, prec);
        fmpz_set_mpz(ci, c[i]);
        acb_add_fmpz(out, out, ci, prec);
    }
    fmpz_clear(ci);
}

static double arb_to_d(const arb_t x)      { return arf_get_d(arb_midref(x), ARF_RND_NEAR); }
static double arb_ubound_d(const arb_t x, slong prec) {
    arf_t u; arf_init(u); arb_get_ubound_arf(u, x, prec);
    double d = arf_get_d(u, ARF_RND_UP); arf_clear(u); return d;
}
static double arb_lbound_d(const arb_t x, slong prec) {
    arf_t l; arf_init(l); arb_get_lbound_arf(l, x, prec);
    double d = arf_get_d(l, ARF_RND_DOWN); arf_clear(l); return d;
}

/* Exact q x q rational solve  A s = b  (A,b integer, mpz).  Returns false
 * if singular; else s[] (mpq) holds the solution. */
static bool solve_rational(const mpz_t* A, const mpz_t* b, int q, mpq_t* s) {
    mpq_t* M = malloc(sizeof(mpq_t) * (size_t)q * (size_t)(q + 1));
    for (int i = 0; i < q; i++) {
        for (int j = 0; j < q; j++) { mpq_init(M[i*(q+1)+j]); mpq_set_z(M[i*(q+1)+j], A[i*q+j]); }
        mpq_init(M[i*(q+1)+q]); mpq_set_z(M[i*(q+1)+q], b[i]);
    }
    bool ok = true;
    mpq_t f, t; mpq_init(f); mpq_init(t);
    for (int col = 0; col < q && ok; col++) {
        int piv = -1;
        for (int row = col; row < q; row++) if (mpq_sgn(M[row*(q+1)+col]) != 0) { piv = row; break; }
        if (piv < 0) { ok = false; break; }
        if (piv != col) for (int j = 0; j <= q; j++) mpq_swap(M[piv*(q+1)+j], M[col*(q+1)+j]);
        for (int row = 0; row < q; row++) {
            if (row == col || mpq_sgn(M[row*(q+1)+col]) == 0) continue;
            mpq_div(f, M[row*(q+1)+col], M[col*(q+1)+col]);
            for (int j = col; j <= q; j++) { mpq_mul(t, f, M[col*(q+1)+j]); mpq_sub(M[row*(q+1)+j], M[row*(q+1)+j], t); }
        }
    }
    if (ok) for (int i = 0; i < q; i++) mpq_div(s[i], M[i*(q+1)+q], M[i*(q+1)+i]);
    mpq_clear(f); mpq_clear(t);
    for (int i = 0; i < q * (q + 1); i++) mpq_clear(M[i]);
    free(M);
    return ok;
}

/* One de Weger reduction (Prop 3.2, inhomogeneous).  delta, mu[0..q-1] are
 * arb at `prec` bits; returns the new bound (>=0) if the hypothesis holds,
 * or -1 to signal "grow c0 and retry". */
static double thue_reduce_once(int q, const arb_t delta, arb_srcptr mu,
                               double K1, double K2, double K3,
                               const mpz_t c0, slong prec) {
    /* lattice basis vectors as rows (integer): row j = e_j + round(c0*mu_j)*e_{q-1}
     * for j<q-1, and row q-1 = round(c0*mu_{q-1})*e_{q-1}. */
    mpz_t* basis = malloc(sizeof(mpz_t) * (size_t)q * (size_t)q);
    for (int i = 0; i < q * q; i++) mpz_init(basis[i]);
    mpq_t* rows = malloc(sizeof(mpq_t) * (size_t)q * (size_t)q);
    for (int i = 0; i < q * q; i++) mpq_init(rows[i]);

    arb_t prod; arb_init(prod);
    fmpz_t c0f, rc; fmpz_init(c0f); fmpz_init(rc);
    fmpz_set_mpz(c0f, c0);
    for (int j = 0; j < q; j++) {
        if (j < q - 1) mpz_set_ui(basis[j*q+j], 1);
        arb_mul_fmpz(prod, &mu[j], c0f, prec);
        arf_get_fmpz(rc, arb_midref(prod), ARF_RND_NEAR);
        fmpz_get_mpz(basis[j*q + (q-1)], rc);
        for (int i = 0; i < q; i++) mpq_set_z(rows[j*q + i], basis[j*q + i]);
    }
    arb_clear(prod); fmpz_clear(c0f); fmpz_clear(rc);

    mpq_t mingso; mpq_init(mingso);
    int dep = lll_reduce_q(rows, q, q, &mingso);
    mpq_clear(mingso);
    double result = -1.0;
    if (!dep) {
        /* reduced basis rows r_i (integer-valued mpq; take the numerator) */
        for (int i = 0; i < q * q; i++) mpz_set(basis[i], mpq_numref(rows[i]));
        /* |b_1|^2 = sum r_0[i]^2 (first reduced row) */
        mpz_t b1sq, tmp; mpz_init_set_ui(b1sq, 0); mpz_init(tmp);
        for (int i = 0; i < q; i++) { mpz_mul(tmp, basis[0*q+i], basis[0*q+i]); mpz_add(b1sq, b1sq, tmp); }

        /* solve M^T s = x, M^T[a][b] = r_b[a] (columns = reduced rows), x = (0..0,-round(c0 delta)) */
        mpz_t* MT = malloc(sizeof(mpz_t) * (size_t)q * (size_t)q);
        for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) { mpz_init_set(MT[a*q+b], basis[b*q+a]); }
        mpz_t* xv = malloc(sizeof(mpz_t) * (size_t)q);
        for (int i = 0; i < q; i++) mpz_init_set_ui(xv[i], 0);
        { arb_t pr; arb_init(pr); fmpz_t c0f2, rc2; fmpz_init(c0f2); fmpz_init(rc2);
          fmpz_set_mpz(c0f2, c0); arb_mul_fmpz(pr, delta, c0f2, prec);
          arf_get_fmpz(rc2, arb_midref(pr), ARF_RND_NEAR); fmpz_get_mpz(xv[q-1], rc2);
          mpz_neg(xv[q-1], xv[q-1]); arb_clear(pr); fmpz_clear(c0f2); fmpz_clear(rc2); }

        mpq_t* s = malloc(sizeof(mpq_t) * (size_t)q);
        for (int i = 0; i < q; i++) mpq_init(s[i]);
        if (solve_rational(MT, xv, q, s)) {
            /* i* = largest i with s_i not integer; ||s_i*|| = dist to nearest int */
            int istar = -1;
            mpq_t dist; mpq_init(dist);
            double sdist = 0.0;
            for (int i = q - 1; i >= 0; i--) {
                if (mpz_cmp_ui(mpq_denref(s[i]), 1) != 0) {
                    istar = i;
                    /* distance to nearest integer */
                    mpz_t nq, rem2; mpz_init(nq); mpz_init(rem2);
                    mpz_t num2, den2; mpz_init_set(num2, mpq_numref(s[i])); mpz_init_set(den2, mpq_denref(s[i]));
                    mpz_fdiv_qr(nq, rem2, num2, den2);           /* 0 <= rem2 < den2 */
                    mpq_set_z(dist, rem2); mpq_t dq; mpq_init(dq); mpq_set_z(dq, den2);
                    mpq_div(dist, dist, dq);                      /* frac part in [0,1) */
                    double fr = mpq_get_d(dist);
                    sdist = fr < 0.5 ? fr : 1.0 - fr;
                    mpq_clear(dq); mpz_clear(nq); mpz_clear(rem2); mpz_clear(num2); mpz_clear(den2);
                    break;
                }
            }
            mpq_clear(dist);
            /* hypothesis: 2^{-(q-1)} * sdist^2 * |b_1|^2 >= (4q^2+3q-3/4) * K3^2 */
            double b1sq_d = mpz_get_d(b1sq);
            double lhs = ldexp(sdist * sdist * b1sq_d, -(q - 1));
            double rhs = (4.0*q*q + 3.0*q - 0.75) * K3 * K3;
            if (getenv("THUE_DEBUG"))
                fprintf(stderr, "[thue]     reduce: istar=%d sdist=%.4g |b1|^2=%.4g lhs=%.4g rhs=%.4g\n",
                        istar, sdist, b1sq_d, lhs, rhs);
            if (istar >= 0) {
                if (lhs >= rhs) {
                    /* A < (1/K2) * log(c0 * K1 / (q * K3)) */
                    double c0d = mpz_get_d(c0);
                    result = (1.0 / K2) * log(c0d * K1 / (q * K3));
                }
            }
        }
        for (int i = 0; i < q; i++) mpq_clear(s[i]);
        free(s);
        for (int i = 0; i < q; i++) mpz_clear(xv[i]);
        free(xv);
        for (int i = 0; i < q * q; i++) mpz_clear(MT[i]);
        free(MT);
        mpz_clear(b1sq); mpz_clear(tmp);
    }

    for (int i = 0; i < q * q; i++) { mpz_clear(basis[i]); mpq_clear(rows[i]); }
    free(basis); free(rows);
    return result;
}

/* Maximal Q-independent subset of {mu_0..mu_{q-1}}, greedily.  mu_i is dropped
 * when a small-integer relation c_i*mu_i + sum_{j in nz} c_j*mu_{nz[j]} = 0
 * (|c| <= CMAX, c_i != 0) is found -- detected rigorously as an arb enclosure
 * of the combination whose whole ball lies below 1e-20 (a genuine vanishing;
 * a non-relation is O(1)).  Writes the independent indices to nz[], returns
 * their count.  (A relation with coefficients above CMAX is simply not found,
 * so those mu stay "independent", the later reduction degenerates, and the
 * caller DECLINEs -- never a wrong bound.) */
static int find_independent_subset(arb_srcptr mu, int q, slong prec, int* nz) {
    enum { CMAX = 8 };
    int p = 0;
    arb_t val, term; arb_init(val); arb_init(term);
    for (int i = 0; i < q; i++) {
        bool dependent = false;
        long range = 2 * CMAX + 1;
        long total = 1; for (int j = 0; j < p; j++) total *= range;
        for (int ci = 1; ci <= CMAX && !dependent; ci++) {
            for (long idx = 0; idx < total && !dependent; idx++) {
                long tt = idx; int cj[16];
                for (int j = 0; j < p; j++) { cj[j] = (int)(tt % range) - CMAX; tt /= range; }
                arb_mul_si(val, mu + i, ci, prec);
                for (int j = 0; j < p; j++) { arb_mul_si(term, mu + nz[j], cj[j], prec); arb_add(val, val, term, prec); }
                arb_abs(val, val);
                if (arb_ubound_d(val, prec) < 1e-20) dependent = true;   /* genuine relation */
            }
        }
        if (!dependent) nz[p++] = i;
    }
    arb_clear(val); arb_clear(term);
    return p;
}

/* Iterate thue_reduce_once (growing c0 until the hypothesis holds) to shrink
 * the bound on max|a_i| for the linear form (delta, mu[0..q-1]).  Returns the
 * reduced bound, or -1 if the form is degenerate (LLL keeps finding a short
 * vector: mu's are Q-dependent) and cannot be reduced. */
static double reduce_form(int q, const arb_t delta, arb_srcptr mu,
                          double K1, double K2, double K3init, slong prec, bool dbg) {
    double curK3 = K3init;
    bool improved_ever = false;
    /* c0 digits capped so log2(c0) <= prec - 128 (128 guard bits): floor(c0*mu)
     * stays EXACT at the working precision.  If the reduction hypothesis cannot
     * be met within this cap, reduce_form returns -1 and the caller escalates. */
    double lc_cap = ((double)prec - 128.0) / 3.3219280949;
    for (int round = 0; round < 30; round++) {
        double newB = -1.0, margin = 12.0;
        while (margin < 220.0) {
            double lc = q * log10(curK3) + margin + 3.0;
            if (lc > lc_cap) break;
            mpz_t c0; mpz_init(c0);
            mpz_ui_pow_ui(c0, 10, (unsigned long)ceil(lc));
            newB = thue_reduce_once(q, delta, mu, K1, K2, curK3, c0, prec);
            if (dbg) fprintf(stderr, "[thue]     reduce q=%d round=%d c0~10^%ld %.6g -> %.6g\n",
                             q, round, (long)mpz_sizeinbase(c0,10)-1, curK3, newB);
            mpz_clear(c0);
            if (newB >= 0) break;
            margin += 24.0;
        }
        if (newB < 0) break;
        improved_ever = true;
        if (newB >= curK3 - 0.5) { curK3 = newB; break; }
        curK3 = newB;
        if (curK3 < 5) break;
    }
    return improved_ever ? curK3 : -1.0;
}

/* Solve  A < a + b*log(A)  for the largest A (downward fixed point). */
static double solve_log_bound(double a, double b) {
    double A = 1e300;
    for (int it = 0; it < 2000; it++) {
        double next = a + b * log(A > 2.718 ? A : 2.718);
        if (next >= A) break;               /* converged / would grow */
        if (next < 2.0) { A = 2.0; break; }
        if (A - next < 1e-6 * A) { A = next; break; }
        A = next;
    }
    return A;
}

/* ------------------------------------------------------------------ *
 *  Bounded-norm representatives (general |m| != 1).                    *
 *                                                                      *
 *  For a monic form, N(x - theta*y) = F(x,y) = m, so beta = x - theta*y *
 *  is an algebraic integer of norm m.  Every such beta is mu * unit,   *
 *  where mu ranges over representatives of {N(.)=m}/units.  This        *
 *  enumerates a set of mu that COVERS every unit orbit (over-coverage   *
 *  is safe: solbuf dedups the final (x,y)).                             *
 *                                                                      *
 *  RANK-1 complex cubic (r1=1, r2=1).  A canonical orbit rep can be     *
 *  reduced by the fundamental unit eps to |sigma1(mu)| in [1, e^R)      *
 *  (R = regulator); the norm then forces |sigma2(mu)| <= sqrt|m|.       *
 *  Those embedding bounds map, through the inverse Vandermonde, to a    *
 *  finite integer coordinate box; keep the points with N == m.         *
 *                                                                      *
 *  Returns the count (>= 0) and mallocs *reps (each a mpz_t[3]); the    *
 *  caller frees via thue_free_reps.  Returns -1 to DECLINE (box too     *
 *  large, or not a rank-1 complex cubic).                               */
static int thue_norm_reps_cubic11(NumberField* K, NFUnits* U, const mpz_t m,
                                  mpz_t*** reps_out) {
    /* Cubic; rank-1 complex (r1=1,r2=1) or rank-2 totally real (r1=3,r2=0). */
    if (K->deg != 3) return -1;
    bool rank1 = (nf_r1(K) == 1 && nf_r2(K) == 1 && nf_units_rank(U) == 1);
    bool rank2 = (nf_r1(K) == 3 && nf_r2(K) == 0 && nf_units_rank(U) == 2);
    if (!rank1 && !rank2) return -1;
    if (nf_ok_index(K) != 1) return -1;   /* non-monogenic: mu-enumeration over Z[theta]
                                           * would miss O_K\Z[theta] reps -> DECLINE (M2 is
                                           * monogenic-only; general m over O_K is a follow-on). */
    const int n = 3;
    slong prec = K->prec + 16;
    acb_srcptr rt = K->roots;                 /* rt[0] real, rt[1] Im>0, rt[2] conj */
    double sqm = sqrt(fabs(mpz_get_d(m)));

    /* per-embedding upper bounds bnd[i] on |sigma_i(mu)| for a mu reduced into a
     * fundamental domain of the unit log-lattice. */
    double bnd[3];
    acb_t tmpc; acb_init(tmpc); arb_t ab; arb_init(ab);
    if (rank1) {
        /* |sigma1(mu)| in [1, e^R); norm forces |sigma2| = |sigma3| <= sqrt|m|. */
        const mpz_t* uc = nf_units_coords(U, 0);
        embed_acb(uc, n, rt + 0, tmpc, prec); acb_abs(ab, tmpc, prec);
        double R = fabs(log(arb_to_d(ab)));
        bnd[0] = exp(R); bnd[1] = sqm; bnd[2] = sqm;
    } else {
        /* rank-2 totally real: a rep reduced into the fundamental parallelogram
         * of <L(eps1),L(eps2)> has, PER embedding i, |L_i(mu)| <= |L_i(eps1)| +
         * |L_i(eps2)| + |log|m||/3 (over-cover is safe). */
        double Li[2][3];
        mpz_t Dd; mpz_init(Dd); nf_units_denom(U, Dd); double logD = log(mpz_get_d(Dd)); mpz_clear(Dd);
        for (int k = 0; k < 2; k++) {
            const mpz_t* uc = nf_units_coords(U, k);
            for (int i = 0; i < 3; i++) { embed_acb(uc, n, rt + i, tmpc, prec); acb_abs(ab, tmpc, prec);
                Li[k][i] = fabs(log(arb_to_d(ab)) - logD); }
        }
        double shift = fabs(log(fabs(mpz_get_d(m)))) / 3.0;
        for (int i = 0; i < 3; i++) bnd[i] = exp(Li[0][i] + Li[1][i] + shift + 1.0);
    }

    /* inverse Vandermonde V^{-1}, V[i][j] = rt[i]^j; coords = V^{-1}(s1,s2,s3) */
    acb_mat_t V, Vi; acb_mat_init(V, 3, 3); acb_mat_init(Vi, 3, 3);
    for (int i = 0; i < 3; i++) {
        acb_one(acb_mat_entry(V, i, 0));
        acb_set(acb_mat_entry(V, i, 1), rt + i);
        acb_mul(acb_mat_entry(V, i, 2), rt + i, rt + i, prec);
    }
    int inv_ok = acb_mat_inv(Vi, V, prec);
    long C[3]; bool feasible = (inv_ok != 0); double boxprod = 1.0;
    for (int i = 0; i < 3 && feasible; i++) {
        double s = 0;
        for (int j = 0; j < 3; j++) { acb_abs(ab, acb_mat_entry(Vi, i, j), prec); s += arb_ubound_d(ab, prec) * bnd[j]; }
        long Ci = (long)ceil(s * 1.05) + 2;    /* +5% + 2: over-cover (safe) */
        C[i] = Ci; boxprod *= (2.0 * Ci + 1);
        if (boxprod > 5e6) feasible = false;   /* too large -> decline */
    }
    acb_mat_clear(V); acb_mat_clear(Vi); acb_clear(tmpc); arb_clear(ab);
    if (!feasible) return -1;

    /* enumerate the coordinate box; keep N(mu) == m exactly */
    mpz_t** reps = NULL; int nrep = 0, cap = 0;
    mpz_t cc[3], nrm; for (int i = 0; i < 3; i++) mpz_init(cc[i]); mpz_init(nrm);
    long W0 = 2*C[0]+1, W1 = 2*C[1]+1, W2 = 2*C[2]+1;
    for (long a0 = 0; a0 < W0; a0++) { mpz_set_si(cc[0], a0 - C[0]);
      for (long a1 = 0; a1 < W1; a1++) { mpz_set_si(cc[1], a1 - C[1]);
        for (long a2 = 0; a2 < W2; a2++) { mpz_set_si(cc[2], a2 - C[2]);
          if (mpz_sgn(cc[0]) == 0 && mpz_sgn(cc[1]) == 0 && mpz_sgn(cc[2]) == 0) continue;
          if (!nf_norm_int(K, (const mpz_t*)cc, nrm)) continue;
          if (mpz_cmp(nrm, m) != 0) continue;
          if (nrep == cap) { cap = cap ? cap*2 : 8; reps = realloc(reps, sizeof(mpz_t*) * (size_t)cap); }
          reps[nrep] = malloc(sizeof(mpz_t) * 3);
          for (int k = 0; k < 3; k++) mpz_init_set(reps[nrep][k], cc[k]);
          nrep++;
        } } }
    for (int i = 0; i < 3; i++) mpz_clear(cc[i]);
    mpz_clear(nrm);
    *reps_out = reps;
    return nrep;
}

static void thue_free_reps(mpz_t** reps, int nrep) {
    if (!reps) return;
    for (int i = 0; i < nrep; i++) { for (int k = 0; k < 3; k++) mpz_clear(reps[i][k]); free(reps[i]); }
    free(reps);
}

/* Rigorous exponent bound.  Mreps/nM give the bounded-norm representatives mu
 * (each a coord vector, |m| != 1 case); pass NULL/0 for the |m| == 1 case, where
 * mu is implicitly {+-1} (mu_ = mu_+ = 1, no mu-ratio in the linear form) and the
 * function behaves exactly as before.  Every mu-dependent constant is an
 * OVER-estimate: a too-large bound is safe (the reduction shrinks it, or we
 * decline), a too-small bound would miss solutions. */
static ThueBound thue_exponent_bound(NumberField* K, NFUnits* U,
                                     const mpz_t* form, int n, const mpz_t m,
                                     mpz_t* const* Mreps, int nM) {
    ThueBound R = { false, 0, 0, 0 };
    int s = nf_r1(K), t = nf_r2(K), r = nf_units_rank(U);
    if (s < 1 || r < 1) return R;               /* need a real root */

    /* Adaptive precision: the constants and K3 only need a modest base; the
     * reduction's precision is sized from K3 and escalated on demand (below).
     * The old blanket 1600 bits paid the worst case on every equation. */
    const slong BASE = 320;
    if (!nf_ensure_prec(K, BASE)) return R;
    slong prec = K->prec;
    acb_srcptr rt = K->roots;                   /* n roots, ordered */

    /* --- M-set constants (mu_lo, mu_hi, log-height) --- over-estimates --- *
     * mu_lo_eff = min over M of |mu^(i)|, clamped <= 1 so C4 >= the |m|=1
     * value; mu_hi_eff = max, clamped >= 1; logH_mu bounds the extra height the
     * mu-ratio adds to the linear form's leading algebraic number (-> C7/K3). */
    double mu_lo_eff = 1.0, mu_hi_eff = 1.0, logH_mu = 0.0;
    if (nM > 0) {
        double lo = 1e300, hi = 0.0;
        acb_t z; acb_init(z); arb_t ab; arb_init(ab);
        for (int mr = 0; mr < nM; mr++)
            for (int h = 0; h < n; h++) {
                embed_acb((const mpz_t*)Mreps[mr], n, rt + h, z, prec);
                acb_abs(ab, z, prec);
                double v = arb_to_d(ab);
                if (v < lo) lo = v;
                if (v > hi) hi = v;
            }
        arb_clear(ab); acb_clear(z);
        if (lo > 0 && lo < 1e300) {
            mu_lo_eff = fmin(lo, 1.0);
            mu_hi_eff = fmax(hi, 1.0);
            logH_mu = (log(hi) - log(lo)) + fabs(log(fabs(mpz_get_d(m)))) + 1.0;
        }
    }

    /* g(x) = F(x,1): g_coeff[i] = form[n-i]; g'(x) coeffs */
    /* --- root-difference constants --- */
    double C1, C2, C3, C4, maxsep = 0.0;
    {
        /* min |g'(root_i)| over real i, and 2^{n-1}|m| */
        acb_t gp, zp, term; acb_init(gp); acb_init(zp); acb_init(term);
        arb_t ab; arb_init(ab);
        double mingp = 1e300;
        fmpz_t cf; fmpz_init(cf);
        for (int i = 0; i < s; i++) {
            /* g'(z) = sum_{i=1}^{n} i*g_coeff[i] z^{i-1}, g_coeff[i]=form[n-i] */
            acb_zero(gp);
            for (int d = n; d >= 1; d--) {
                acb_mul(gp, gp, rt + i, prec);
                fmpz_set_mpz(cf, form[n - d]);           /* g_coeff[d] */
                fmpz_mul_si(cf, cf, d);
                acb_add_fmpz(gp, gp, cf, prec);
            }
            acb_abs(ab, gp, prec);
            double v = arb_lbound_d(ab, prec);           /* lower bound -> C1 upper */
            if (v < mingp) mingp = v;
        }
        fmpz_clear(cf);
        double twopow = ldexp(1.0, n - 1) * fabs(mpz_get_d(m));
        C1 = (mingp > 0) ? twopow / mingp : 1e300;

        /* C2 = 0.5 min sep ; maxsep = max sep ; C3 = max ratio */
        double minsep = 1e300;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
            acb_sub(term, rt + i, rt + j, prec); acb_abs(ab, term, prec);
            double lo = arb_lbound_d(ab, prec), hi = arb_ubound_d(ab, prec);
            if (lo < minsep) minsep = lo;
            if (hi > maxsep) maxsep = hi;
        }
        C2 = 0.5 * minsep;
        C3 = 0.0;
        for (int i0 = 0; i0 < s; i0++)
            for (int j = 0; j < n; j++) { if (j == i0) continue;
                for (int k = 0; k < n; k++) { if (k == i0 || k == j) continue;
                    acb_sub(zp, rt + i0, rt + j, prec);
                    acb_sub(term, rt + i0, rt + k, prec);
                    acb_div(zp, zp, term, prec); acb_abs(ab, zp, prec);
                    double v = arb_ubound_d(ab, prec);
                    if (v > C3) C3 = v;
                }
            }
        C4 = (0.5 + maxsep) / mu_lo_eff;                 /* mu_ = min_M |mu^(i)| */
        acb_clear(gp); acb_clear(zp); acb_clear(term); arb_clear(ab);
    }
    if (!(C1 > 0 && C2 > 0 && C3 > 0)) return R;

    /* --- unit embeddings: eps_acb[l][h], logabs[l][h] --- */
    acb_t* eabs = malloc(sizeof(acb_t) * (size_t)r * (size_t)n);
    double* logabs = malloc(sizeof(double) * (size_t)r * (size_t)n);
    double* Hunit = malloc(sizeof(double) * (size_t)r);          /* log H_l */
    {
        arb_t ab; arb_init(ab);
        for (int l = 0; l < r; l++) {
            const mpz_t* uc = nf_units_coords(U, l);
            double mx = -1e300, mn = 1e300;
            for (int h = 0; h < n; h++) {
                acb_init(eabs[l*n + h]);
                embed_acb(uc, n, rt + h, eabs[l*n + h], prec);
                acb_abs(ab, eabs[l*n + h], prec);
                double la = log(arb_to_d(ab));
                logabs[l*n + h] = la;
                if (la > mx) mx = la;
                if (la < mn) mn = la;
            }
            Hunit[l] = mx - mn;                                  /* log H_l */
        }
        arb_clear(ab);
    }

    /* --- C5 = min((n-1) min_I N[U_I^-1], max_I N[U_I^-1]) over r-subsets I --- */
    double C5;
    {
        int idx[8]; for (int i = 0; i < r; i++) idx[i] = i;
        double minN = 1e300, maxN = 0.0;
        /* iterate all r-subsets of {0..n-1} */
        bool done = false;
        while (!done) {
            /* build r x r matrix Ul[i][l] = logabs[l][idx[i]]; invert (double) */
            double Mmat[9], Inv[9];
            for (int i = 0; i < r; i++) for (int l = 0; l < r; l++) Mmat[i*r+l] = logabs[l*n + idx[i]];
            /* Gauss-Jordan inverse */
            double Aug[9*2]; for (int i=0;i<r;i++){ for(int j=0;j<r;j++) Aug[i*(2*r)+j]=Mmat[i*r+j]; for(int j=0;j<r;j++) Aug[i*(2*r)+r+j]=(i==j)?1.0:0.0; }
            bool sing=false;
            for (int col=0; col<r && !sing; col++){
                int piv=-1; double best=1e-12; for(int row=col;row<r;row++){ double a=fabs(Aug[row*(2*r)+col]); if(a>best){best=a;piv=row;} }
                if(piv<0){sing=true;break;}
                if(piv!=col) for(int j=0;j<2*r;j++){ double tt=Aug[piv*(2*r)+j]; Aug[piv*(2*r)+j]=Aug[col*(2*r)+j]; Aug[col*(2*r)+j]=tt; }
                double d=Aug[col*(2*r)+col]; for(int j=0;j<2*r;j++) Aug[col*(2*r)+j]/=d;
                for(int row=0;row<r;row++){ if(row==col) continue; double f=Aug[row*(2*r)+col]; if(f!=0) for(int j=0;j<2*r;j++) Aug[row*(2*r)+j]-=f*Aug[col*(2*r)+j]; }
            }
            if(!sing){
                for(int i=0;i<r;i++) for(int j=0;j<r;j++) Inv[i*r+j]=Aug[i*(2*r)+r+j];
                double N=0.0; for(int i=0;i<r;i++){ double rs=0; for(int j=0;j<r;j++) rs+=fabs(Inv[i*r+j]); if(rs>N) N=rs; }
                N *= 1.0001;                                     /* tiny safety -> upper bound */
                if(N<minN) minN=N;
                if(N>maxN) maxN=N;
            }
            /* next r-subset (idx ascending) */
            int p=r-1; while(p>=0 && idx[p]==n-r+p) p--;
            if(p<0) done=true; else { idx[p]++; for(int i=p+1;i<r;i++) idx[i]=idx[i-1]+1; }
        }
        if (minN > 1e299) {
            for (int i = 0; i < r * n; i++) acb_clear(eabs[i]);
            free(eabs); free(logabs); free(Hunit); return R;
        }
        double a1 = (n-1)*minN, a2 = maxN;
        C5 = a1 < a2 ? a1 : a2;
    }

    double C6 = 1.39 * C1 * C3 * pow(C4, n) / C2;
    double K1 = C6, K2 = (double)n / C5;
    bool dbg = getenv("THUE_DEBUG") != NULL;
    if (dbg) fprintf(stderr, "[thue] n=%d s=%d t=%d r=%d C1=%.4g C2=%.4g C3=%.4g C4=%.4g C5=%.4g C6=%.4g K2=%.4g\n",
                     n, s, t, r, C1, C2, C3, C4, C5, C6, K2);

    /* --- Y-thresholds --- */
    double Y1d, Y2p;
    {
        double y0 = 1.0;
        if (t >= 1) {
            /* min|g'(complex root)| and min|Im| over complex reps */
            acb_t gp; acb_init(gp); arb_t ab; arb_init(ab); fmpz_t cf; fmpz_init(cf);
            double mingp = 1e300, minim = 1e300;
            for (int i = s; i < s + t; i++) {
                acb_zero(gp);
                for (int d = n; d >= 1; d--) { acb_mul(gp, gp, rt+i, prec); fmpz_set_mpz(cf, form[n-d]); fmpz_mul_si(cf,cf,d); acb_add_fmpz(gp,gp,cf,prec); }
                acb_abs(ab, gp, prec); double v = arb_lbound_d(ab, prec); if (v<mingp) mingp=v;
                arb_abs(ab, acb_imagref(rt+i)); double im = arb_lbound_d(ab, prec); if (im<minim) minim=im;
            }
            if (mingp>0 && minim>0) y0 = ceil(pow(ldexp(1.0,n-1)*fabs(mpz_get_d(m))/(mingp*minim), 1.0/n));
            acb_clear(gp); arb_clear(ab); fmpz_clear(cf);
        }
        Y1d = fmax(y0, ceil(pow(4*C1, 1.0/(n-2))));
        double y2s = fmax(Y1d, ceil(pow(2*C1*C3/C2, 1.0/n)));
        Y2p = fmax(y2s, fmax(2*pow(fabs(mpz_get_d(m)),1.0/n), mu_hi_eff/C2));  /* mu_+ = max_M |mu^(i)| */
    }

    /* --- Baker (Waldschmidt) initial bound K3 --- */
    double K3;
    {
        int N = r + 1;
        mpz_t discabs; mpz_init(discabs); nf_disc(K, discabs); mpz_abs(discabs, discabs);
        double V0 = log(mpz_get_d(discabs)) + log(C3) + logH_mu;   /* +mu-ratio height */
        mpz_clear(discabs);
        double* V = malloc(sizeof(double) * (size_t)N);
        V[0] = V0; for (int l = 0; l < r; l++) V[l+1] = Hunit[l];
        /* sort ascending */
        for (int i=1;i<N;i++){ double key=V[i]; int j=i-1; while(j>=0&&V[j]>key){V[j+1]=V[j];j--;} V[j+1]=key; }
        double D = 1.0; for (int i = 2; i <= n; i++) D *= i;        /* n! over-estimate */
        double e_n = fmin(fmin(8.0*N+51, 10.0*N+33), 9.0*N+39);
        double prodV = 1.0; for (int i=0;i<N;i++) prodV *= V[i];
        double Vn1p = fmax(V[N-2 >= 0 ? N-2 : 0], 1.0);
        double Vnp  = fmax(V[N-1], 1.0);
        double C7 = pow(2.0, e_n) * pow((double)N, 2.0*N) * pow(D, N+2) * prodV * log(exp(1.0)*D*Vn1p);
        double C8 = log(exp(1.0)*D*Vnp);
        double realcase_min = (s >= 3) ? 0.0 : 1.0;   /* placeholder; C8' below */
        (void)realcase_min;
        double C8p = C8;                              /* real case; complex adds log r (handled per-i0 below via +log r) */
        /* use complex-safe C8' = C8 + log(max(r,1)) as an upper bound (valid for both) */
        C8p = C8 + log((double)(r > 1 ? r : 1));
        double a = (log(C6) + C7 * C8p) / K2;
        double b = C7 / K2;
        K3 = solve_log_bound(a, b);
        if (dbg) fprintf(stderr, "[thue] N=%d C7=%.4g C8'=%.4g K3=%.6g Y2p=%.4g\n", N, C7, C8p, K3, Y2p);
        free(V);
    }
    if (!(K3 > 0) || K3 > 1e300) {
        for (int i = 0; i < r * n; i++) acb_clear(eabs[i]);
        free(eabs); free(logabs); free(Hunit); return R;
    }

    /* mu representative embeddings muemb[mr*n + h] (|m| != 1 case): recomputed
     * on every precision escalation, alongside the units. */
    int nMr = (nM > 0) ? nM : 0;
    acb_t* muemb = (nMr > 0) ? malloc(sizeof(acb_t) * (size_t)nMr * (size_t)n) : NULL;
    for (int i = 0; i < nMr * n; i++) acb_init(muemb[i]);
    for (int mr = 0; mr < nMr; mr++)
        for (int h = 0; h < n; h++)
            embed_acb((const mpz_t*)Mreps[mr], n, rt + h, muemb[mr*n + h], prec);

    /* --- adaptive-precision reduction ---
     * Size the working precision from K3 (the reduction's c0 ~ K3^q * small
     * margin), then ESCALATE if the reduction cannot certify at that precision.
     * Easy fields (small K3) finish at a few hundred bits instead of 1600. */
    slong red_prec = (slong)ceil(((double)(r + 1) * log10(K3 > 10 ? K3 : 10) + 30.0) * 3.3219280949) + 224;
    if (red_prec < BASE) red_prec = BASE;
    const slong PREC_CAP = 3200;
    double Bmax = 0.0;
    bool bound_ok = false;
    for (int esc = 0; esc < 6 && !bound_ok; esc++) {
        if (red_prec > prec) {
            if (!nf_ensure_prec(K, red_prec)) { red_prec *= 2; if (red_prec > PREC_CAP) break; continue; }
            prec = K->prec; rt = K->roots;               /* roots re-isolated at higher prec */
            for (int l = 0; l < r; l++)                  /* re-embed the units at the new prec */
                for (int h = 0; h < n; h++)
                    embed_acb(nf_units_coords(U, l), n, rt + h, eabs[l*n + h], prec);
            for (int mr = 0; mr < nMr; mr++)             /* re-embed the mu reps too */
                for (int h = 0; h < n; h++)
                    embed_acb((const mpz_t*)Mreps[mr], n, rt + h, muemb[mr*n + h], prec);
        }
        if (dbg) fprintf(stderr, "[thue] reduction at prec=%ld bits (esc=%d)\n", (long)prec, esc);
        Bmax = 0.0;
        /* Loop over the bounded-norm representatives (once, mu=1, for |m|=1). */
        for (int mrep = 0; mrep < (nMr > 0 ? nMr : 1) && K3 < 1e301; mrep++) {
        for (int i0 = 0; i0 < s && K3 < 1e301; i0++) {
        /* choose j,k and case */
        bool realcase = (s - 1) >= 2;
        int j, k, q;
        if (realcase) { q = r; j = (i0==0)?1:0; k = (j+1==i0)?(j+2):(j+1); if(k==i0)k++; }
        else          { q = r + 1; j = s; k = s + t; }   /* conjugate pair */
        /* mu[] and delta as arb at prec */
        arb_ptr mu = _arb_vec_init(q);
        arb_t delta; arb_init(delta);
        acb_t z0, zj, zk, ratio; acb_init(z0); acb_init(zj); acb_init(zk); acb_init(ratio);
        acb_sub(zj, rt+i0, rt+j, prec); acb_sub(zk, rt+i0, rt+k, prec);
        acb_div(ratio, zj, zk, prec);
        if (realcase) {
            arb_t ab, lj, lk; arb_init(ab); arb_init(lj); arb_init(lk);
            acb_abs(ab, ratio, prec); arb_log(delta, ab, prec);
            if (nMr > 0) {                                /* delta += log|mu^(k)/mu^(j)| */
                acb_abs(lk, muemb[mrep*n + k], prec); arb_log(lk, lk, prec);
                acb_abs(lj, muemb[mrep*n + j], prec); arb_log(lj, lj, prec);
                arb_add(delta, delta, lk, prec); arb_sub(delta, delta, lj, prec);
            }
            for (int l = 0; l < r; l++) {
                /* mu_l = log|eps_l^(k)| - log|eps_l^(j)|, at FULL precision
                 * (double would lose every digit past ~16, but c0 ~ 10^60). */
                acb_abs(lk, eabs[l*n + k], prec); arb_log(lk, lk, prec);
                acb_abs(lj, eabs[l*n + j], prec); arb_log(lj, lj, prec);
                arb_sub(mu + l, lk, lj, prec);
            }
            arb_clear(ab); arb_clear(lj); arb_clear(lk);
        } else {
            acb_arg(delta, ratio, prec);
            if (nMr > 0) {                                /* delta += arg(mu^(k)/mu^(j)) */
                acb_t rr; arb_t a2; acb_init(rr); arb_init(a2);
                acb_div(rr, muemb[mrep*n + k], muemb[mrep*n + j], prec);
                acb_arg(a2, rr, prec); arb_add(delta, delta, a2, prec);
                acb_clear(rr); arb_clear(a2);
            }
            for (int l = 0; l < r; l++) {
                acb_t rr; acb_init(rr);
                acb_div(rr, eabs[l*n + k], eabs[l*n + j], prec);
                acb_arg(mu + l, rr, prec);
                acb_clear(rr);
            }
            arb_const_pi(mu + r, prec); arb_mul_2exp_si(mu + r, mu + r, 1);   /* 2*pi */
        }
        acb_clear(z0); acb_clear(zj); acb_clear(zk); acb_clear(ratio);

        /* Q-dependence handling (Prop 3.2 case iii): a fundamental unit in a
         * proper subfield has mu_i = 0 (real at the paired embedding), so its
         * exponent does not appear in Lambda and the lattice degenerates.
         * Reduce the independent subform, then bound the FULL exponent max via
         * the "L-trick": A < (1/K2) log(K1 / L), L = min|Lambda'| in range. */
        int nz[16];
        int p = find_independent_subset(mu, q, prec, nz);

        double curbound;
        if (p == q) {
            curbound = reduce_form(q, delta, mu, K1, K2, K3, prec, dbg);
            if (curbound < 0) curbound = 1e18;                  /* unexpected degeneracy -> cap */
        } else {
            double Ared = 0.0;
            if (p > 0) {
                arb_ptr mu2 = _arb_vec_init(p);
                for (int i = 0; i < p; i++) arb_set(mu2 + i, mu + nz[i]);
                Ared = reduce_form(p, delta, mu2, K1, K2, K3, prec, dbg);
                _arb_vec_clear(mu2, p);
            }
            if (Ared < 0) { curbound = 1e18; }                  /* subform still degenerate */
            else {
                long AR = (long)ceil(Ared); if (AR < 0) AR = 0;
                double combos = 1; for (int i = 0; i < p; i++) combos *= (2.0*AR + 1);
                double L = 1e300; bool okL = (combos <= 2e6);
                if (okL) {
                    arb_t val, term; arb_init(val); arb_init(term);
                    long total = 1; for (int i = 0; i < p; i++) total *= (2*AR + 1);
                    int* av = malloc(sizeof(int) * (size_t)(p > 0 ? p : 1));
                    for (long idx = 0; idx < total && okL; idx++) {
                        long tt = idx; for (int i = 0; i < p; i++) { av[i] = (int)(tt % (2*AR+1)) - AR; tt /= (2*AR+1); }
                        arb_set(val, delta);
                        for (int i = 0; i < p; i++) { arb_mul_si(term, mu + nz[i], av[i], prec); arb_add(val, val, term, prec); }
                        arb_abs(term, val);
                        double lo = arb_lbound_d(term, prec);   /* rigorous lower bound on |Lambda'| */
                        if (lo < L) L = lo;
                    }
                    free(av); arb_clear(val); arb_clear(term);
                }
                if (okL && L > 1e-250) curbound = (1.0/K2) * log(K1 / L);
                else curbound = 1e18;                            /* L uncertain/tiny -> decline */
                if (dbg) fprintf(stderr, "[thue]   i0=%d degenerate p=%d/%d Ared=%.4g L=%.4g -> A_full=%.6g\n",
                                 i0, p, q, Ared, L, curbound);
            }
        }
        if (curbound > Bmax) Bmax = curbound;
        _arb_vec_clear(mu, q); arb_clear(delta);
        }  /* i0 */
        }  /* mrep */
        if (Bmax > 0 && Bmax <= 1e6) bound_ok = true;       /* certified at this precision */
        else { red_prec *= 2; if (red_prec > PREC_CAP) break; }   /* escalate precision */
    }
    for (int i = 0; i < nMr * n; i++) acb_clear(muemb[i]);
    free(muemb);

    for (int i = 0; i < r * n; i++) acb_clear(eabs[i]);
    free(eabs); free(logabs); free(Hunit);

    if (!bound_ok) return R;                                 /* could not certify -> DECLINE */

    /* Xmax for the small-|Y| brute force: |X| <= max|root|*Y2p + C1 */
    double maxrootabs = 0.0;
    { arb_t ab; arb_init(ab);
      for (int i = 0; i < n; i++) { acb_abs(ab, rt + i, prec); double v = arb_ubound_d(ab, prec); if (v > maxrootabs) maxrootabs = v; }
      arb_clear(ab); }

    R.ok = true;
    R.B = (long)ceil(Bmax);
    R.Y2p = (long)ceil(Y2p);
    /* For |Y|<=Y2p, min_i|beta^(i)| <= |m|^(1/n) = 1, so some |X-Re(xi^(i))Y|<=1,
     * giving |X| <= max|root|*|Y| + 1. */
    R.Xmax = (long)ceil(maxrootabs * (double)R.Y2p) + 2;
    return R;
}

/* Is the rigorous automatic bound available (needs FLINT). */
static bool thue_bound_available(void) { return true; }

/* ---- core reconstruction/enumeration ---- */

typedef struct { mpz_t* x; mpz_t* y; int n; int cap; } SolBuf;
static void solbuf_push(SolBuf* s, const mpz_t x, const mpz_t y) {
    for (int i = 0; i < s->n; i++)
        if (mpz_cmp(s->x[i], x) == 0 && mpz_cmp(s->y[i], y) == 0) return;  /* dedupe */
    if (s->n == s->cap) {
        s->cap = s->cap ? s->cap * 2 : 16;
        s->x = realloc(s->x, sizeof(mpz_t) * (size_t)s->cap);
        s->y = realloc(s->y, sizeof(mpz_t) * (size_t)s->cap);
    }
    mpz_init_set(s->x[s->n], x); mpz_init_set(s->y[s->n], y); s->n++;
}

/* ================================================================== *
 *  Totally complex field (r1 = 0): the elementary imaginary-part      *
 *  bound.                                                             *
 *                                                                     *
 *  Every root theta_i is non-real, so for real integers x, y          *
 *      |x - theta_i y| >= |Im(theta_i)| * |y|                         *
 *  (the imaginary part alone).  Unlike the real case there is no real *
 *  root for x/y to approach, so NO factor can be small, and           *
 *      |m| = |a0| * prod_i |x - theta_i y|                            *
 *          >= |a0| * |y|^n * prod_i |Im theta_i|,                     *
 *  giving the elementary, rigorous bound                              *
 *      |y| <= ( |m| / (|a0| * prod_i |Im theta_i|) )^{1/n}.           *
 *  No Baker bound, no fundamental units, no torsion enumeration are   *
 *  needed -- exactly the pieces the real-i0 engine lacks here.  For   *
 *  each y in that (small) range the integer x with F(x,y)==m are the  *
 *  integer roots of the univariate F(x,y)-m, found exactly.  a0 = 1   *
 *  (the caller normalised form[0] = +1).  Returns the complete count  *
 *  (>= 0), or -1 to DECLINE (bound not certifiable / infeasible).     */
static int thue_solve_totally_complex(NumberField* K, const mpz_t* form, int n,
                                      const mpz_t m, ThueSol** out) {
    if (nf_r1(K) != 0) return -1;                       /* totally complex only */
    bool edbg = getenv("THUE_DEBUG") != NULL;

    /* rigorous LOWER bound on prod_i |Im(theta_i)| via arb (so |m|/prod is an
     * UPPER bound); escalate precision if any |Im| lower bound underflows. */
    slong prec = K->prec;
    double sum_log_im = 0.0;                            /* sum of log lower-bounds */
    bool ok = false;
    for (int esc = 0; esc < 4 && !ok; esc++) {
        if (!nf_ensure_prec(K, prec)) return -1;
        prec = K->prec;
        acb_srcptr rt = K->roots;
        arb_t im; arb_init(im);
        double s = 0.0; bool good = true;
        for (int i = 0; i < n; i++) {
            arb_abs(im, acb_imagref(rt + i));
            double lo = arb_lbound_d(im, prec);        /* certified |Im| >= lo */
            if (!(lo > 0.0)) { good = false; break; }
            s += log(lo);
        }
        arb_clear(im);
        if (good) { sum_log_im = s; ok = true; } else prec *= 2;
    }
    if (!ok) { if (edbg) fprintf(stderr, "[thue] DECLINE: totally-complex |Im| bound underflow\n"); return -1; }

    /* logY = (log|m| - sum log|Im|) / n is an UPPER bound on log(true |y| max). */
    double logY = (log(fabs(mpz_get_d(m))) - sum_log_im) / (double)n;
    double Yd = exp(logY);
    if (!(Yd >= 0.0) || Yd > 4e6) { if (edbg) fprintf(stderr, "[thue] DECLINE: totally-complex Y too large (%.4g)\n", Yd); return -1; }
    long Y = (long)floor(Yd) + 1;                       /* +1: safe over-estimate */
    if (Y < 0) Y = 0;
    if (edbg) fprintf(stderr, "[thue-prof] totally-complex (r1=0): |y| <= %ld\n", Y);

    /* For each y, the integer x with F(x,y)==m are the integer roots of the
     * univariate F(x,y)-m (exact); verify F==m and collect (solbuf dedups). */
    SolBuf sb; memset(&sb, 0, sizeof(sb));
    mpz_t Yv, fval; mpz_init(Yv); mpz_init(fval);
    fmpz_poly_t Px; fmpz_poly_init(Px);
    mpz_t* xroots = NULL; int nxr = 0, xcap = 0;
    for (long yy = -Y; yy <= Y; yy++) {
        mpz_set_si(Yv, yy);
        thue_form_xpoly(form, n, Yv, m, Px);
        nxr = 0;
        fmpz_poly_int_roots(Px, &xroots, &nxr, &xcap);
        for (int i = 0; i < nxr; i++) {
            eval_form(form, n, xroots[i], Yv, fval);
            if (mpz_cmp(fval, m) == 0) solbuf_push(&sb, xroots[i], Yv);
            mpz_clear(xroots[i]);
        }
    }
    free(xroots); fmpz_poly_clear(Px);
    mpz_clear(Yv); mpz_clear(fval);

    ThueSol* res = malloc(sizeof(ThueSol) * (size_t)(sb.n > 0 ? sb.n : 1));
    for (int i = 0; i < sb.n; i++) { mpz_init_set(res[i].x, sb.x[i]); mpz_init_set(res[i].y, sb.y[i]); }
    *out = res; int rc = sb.n;
    for (int i = 0; i < sb.n; i++) { mpz_clear(sb.x[i]); mpz_clear(sb.y[i]); }
    free(sb.x); free(sb.y);
    return rc;
}

/* Shared engine: monic |a0|=1, |m|=1, monogenic field, torsion {+-1}.
 * Enumerate |b_k| <= bound.  Returns count or -1 (setup decline). */
static int thue_enumerate(const mpz_t* form_in, int n, const mpz_t m_in,
                          int bound, ThueSol** out) {
    if (n < 3) return -1;

    /* normalise leading coeff to +1 (a0 = +-1 required for the monic path) */
    mpz_t* form = malloc(sizeof(mpz_t) * (size_t)(n + 1));
    for (int i = 0; i <= n; i++) mpz_init_set(form[i], form_in[i]);
    mpz_t m; mpz_init_set(m, m_in);
    int rc = -1;

    if (mpz_cmpabs_ui(form[0], 1) != 0) goto done;           /* need |a0| == 1 */
    if (mpz_sgn(form[0]) < 0) { for (int i = 0; i <= n; i++) mpz_neg(form[i], form[i]); mpz_neg(m, m); }
    if (mpz_sgn(form[n]) == 0) goto done;                    /* not a genuine binary form */
    /* |m| == 1: beta = x - theta*y is a UNIT (mu implicitly {+-1}).  |m| != 1:
     * beta is a norm-m integer = mu * unit, mu over bounded-norm reps (M2). */
    bool general_m = (mpz_cmpabs_ui(m, 1) != 0);

    /* monic defining polynomial f(t) = F(t,1): coeff of t^(n-j) is a_j */
    mpz_t* fc = malloc(sizeof(mpz_t) * (size_t)(n + 1));
    for (int i = 0; i <= n; i++) mpz_init_set(fc[i], form[n - i]);   /* fc[i] = a_{n-i} */

    bool edbg = getenv("THUE_DEBUG") != NULL;
    clock_t tprof = clock();
    NumberField* K = nf_field_create((const mpz_t*)fc, n);
    for (int i = 0; i <= n; i++) mpz_clear(fc[i]);
    free(fc);
    if (!K) { if (edbg) fprintf(stderr, "[thue] DECLINE: Gate 1 (non-monogenic/reducible/non-monic)\n"); goto done; }
    if (nf_r1(K) < 1) {
        /* Totally complex (no real embedding): the Baker/unit machinery needs a
         * real type-index i0, but the elementary imaginary-part bound solves it
         * directly and rigorously (any m).  Solves e.g. the cyclotomic quartic
         * x^4+x^3y+x^2y^2+xy^3+y^4 == 1 over Q(zeta_5). */
        int tc = thue_solve_totally_complex(K, (const mpz_t*)form, n, m, out);
        nf_field_free(K);
        if (tc >= 0) { rc = tc; goto done; }
        if (edbg) fprintf(stderr, "[thue] DECLINE: totally-complex bound not certified\n");
        goto done;
    }
    if (edbg) { fprintf(stderr, "[thue-prof] field_create: %.1fms\n", (clock()-tprof)*1000.0/CLOCKS_PER_SEC); tprof = clock(); }

    NFUnits* U = nf_fundamental_units(K);
    if (!U) { if (edbg) fprintf(stderr, "[thue] DECLINE: Gate 2 (fundamental units not certified)\n"); nf_field_free(K); goto done; }
    int r = nf_units_rank(U);
    if (edbg) { fprintf(stderr, "[thue-prof] units (rank %d): %.1fms\n", r, (clock()-tprof)*1000.0/CLOCKS_PER_SEC); tprof = clock(); }

    /* --- bounded-norm representatives mu (general |m| != 1) --- */
    mpz_t** Mreps = NULL; int nM = 0;
    if (general_m) {
        nM = thue_norm_reps_cubic11(K, U, m, &Mreps);   /* rank-1 complex cubic only (M2) */
        if (nM < 0) { if (edbg) fprintf(stderr, "[thue] DECLINE: |m|!=1 unsupported field / mu-box too large\n"); nf_units_free(U); nf_field_free(K); goto done; }
        if (edbg) fprintf(stderr, "[thue-prof] mu-reps (%d, N=m=%ld): %.1fms\n", nM, mpz_get_si(m), (clock()-tprof)*1000.0/CLOCKS_PER_SEC);
        if (nM == 0) {                                  /* no norm-m integer => no solutions */
            *out = malloc(sizeof(ThueSol)); rc = 0;
            nf_units_free(U); nf_field_free(K);
            goto done;
        }
        tprof = clock();
    }

    /* effective enumeration bound: rigorous (Baker/de-Weger) when the caller
     * passes the sentinel, else the explicit test bound. */
    long Y2p = -1, Xmax = -1;
    int B;
    if (bound == THUE_BOUND_RIGOROUS) {
        ThueBound tb = thue_exponent_bound(K, U, (const mpz_t*)form, n, m,
                                           general_m ? Mreps : NULL, general_m ? nM : 0);
        if (!tb.ok) { if (edbg) fprintf(stderr, "[thue] DECLINE: exponent bound not certified\n"); thue_free_reps(Mreps, nM); nf_units_free(U); nf_field_free(K); goto done; }
        B = (int)tb.B; Y2p = tb.Y2p; Xmax = tb.Xmax;
        /* both enumeration boxes must be feasible, else DECLINE (never truncate).
         * Small-|Y| gap: a NARROW x-window (2*Xmax+1 <= 256) is scanned exactly
         * as before (cheapest when Xmax is tiny, e.g. Xmax=0 with a huge Y2p); a
         * WIDE window uses exact univariate root-finding (one degree-n
         * factorisation per y), whose cost is O(Y2p) independent of Xmax.  The
         * gate matches the method, so no previously-solved case regresses and
         * wide-window cases (old O(Y2p*Xmax) blowup) now become feasible. */
        double ebox = (double)(nM > 0 ? nM : 1); for (int k = 0; k < r; k++) ebox *= (2.0 * B + 1);
        bool narrowX = (2.0 * Xmax + 1.0) <= 256.0;
        double ybox = narrowX ? (2.0 * Y2p + 1.0) * (2.0 * Xmax + 1.0)   /* scan area */
                              : (2.0 * Y2p + 1.0);                        /* one factor / row */
        double ybox_cap = narrowX ? 2e8 : 4e6;
        if (ebox > 2e8 || ybox > ybox_cap) { if (edbg) fprintf(stderr, "[thue] DECLINE: enumeration box too large (ebox=%.3g ybox=%.3g narrowX=%d B=%d Y2p=%ld Xmax=%ld)\n", ebox, ybox, (int)narrowX, B, Y2p, Xmax); thue_free_reps(Mreps, nM); nf_units_free(U); nf_field_free(K); goto done; }
        if (edbg) { fprintf(stderr, "[thue-prof] exponent_bound (B=%d Y2p=%ld): %.1fms\n", B, Y2p, (clock()-tprof)*1000.0/CLOCKS_PER_SEC); tprof = clock(); }
    } else {
        B = bound;
    }
    if (B < 0) { thue_free_reps(Mreps, nM); nf_units_free(U); nf_field_free(K); goto done; }

    /* precompute the fundamental units AND their powers eps_k^j, j in [-B,B]:
     * then each of the (2B+1)^r lattice points is just r-1 nf_elem_mul instead
     * of r nf_elem_pow -- the enumeration is the hot loop. */
    mpz_t Udn; mpz_init(Udn); nf_units_denom(U, Udn);   /* O_K basis denominator (1 if monogenic) */
    nf_elem_t* E = malloc(sizeof(nf_elem_t) * (size_t)(r > 0 ? r : 1));
    for (int k = 0; k < r; k++) { nf_elem_init(E[k], K->nf); nfelem_from_coords_den(K->nf, nf_units_coords(U, k), Udn, n, E[k]); }
    mpz_clear(Udn);
    long W = 2L * B + 1;
    nf_elem_t** Epow = malloc(sizeof(nf_elem_t*) * (size_t)(r > 0 ? r : 1));
    for (int k = 0; k < r; k++) {
        Epow[k] = malloc(sizeof(nf_elem_t) * (size_t)W);
        for (long j = 0; j < W; j++) nf_elem_init(Epow[k][j], K->nf);
        nf_elem_one(Epow[k][B], K->nf);                                             /* exp 0 */
        for (long j = B + 1; j < W; j++) nf_elem_mul(Epow[k][j], Epow[k][j-1], E[k], K->nf);   /* exp > 0 */
        nf_elem_t einv; nf_elem_init(einv, K->nf); nf_elem_inv(einv, E[k], K->nf);
        for (long j = B - 1; j >= 0; j--) nf_elem_mul(Epow[k][j], Epow[k][j+1], einv, K->nf);  /* exp < 0 */
        nf_elem_clear(einv, K->nf);
    }

    /* mu factors: the bounded-norm reps (|m| != 1) or the single identity 1
     * (|m| == 1, where beta = +-unit).  beta = mu * prod eps_k^{b_k}. */
    int nmu = general_m ? nM : 1;
    nf_elem_t* Mel = malloc(sizeof(nf_elem_t) * (size_t)nmu);
    for (int i = 0; i < nmu; i++) {
        nf_elem_init(Mel[i], K->nf);
        if (general_m) nfelem_from_coords(K->nf, (const mpz_t*)Mreps[i], n, Mel[i]);
        else           nf_elem_one(Mel[i], K->nf);
    }

    SolBuf sb; memset(&sb, 0, sizeof(sb));
    nf_elem_t u, ubeta, usign;
    nf_elem_init(u, K->nf); nf_elem_init(ubeta, K->nf); nf_elem_init(usign, K->nf);
    mpz_t* coord = malloc(sizeof(mpz_t) * (size_t)n);
    for (int i = 0; i < n; i++) mpz_init(coord[i]);
    mpz_t yneg, fval; mpz_init(yneg); mpz_init(fval);

    /* mixed-radix enumeration of the exponent-index vector in [0, 2B]^r
     * (index j maps to exponent j-B), times each mu, plus the +-1 torsion. */
    int* b = malloc(sizeof(int) * (size_t)(r > 0 ? r : 1));
    long total = 1; for (int k = 0; k < r; k++) total *= W;

    for (int mi = 0; mi < nmu; mi++) {
      for (long idx = 0; idx < total; idx++) {
        long t = idx;
        for (int k = 0; k < r; k++) { b[k] = (int)(t % W); t /= W; }        /* index into Epow */
        nf_elem_set(u, Epow[0][b[0]], K->nf);                              /* r >= 1 always */
        for (int k = 1; k < r; k++) nf_elem_mul(u, u, Epow[k][b[k]], K->nf);
        nf_elem_mul(ubeta, Mel[mi], u, K->nf);                             /* beta_core = mu * unit */
        for (int sgn = 0; sgn < 2; sgn++) {
            if (sgn == 0) nf_elem_set(usign, ubeta, K->nf);
            else          nf_elem_neg(usign, ubeta, K->nf);

            if (!coords_from_nfelem(K->nf, usign, n, coord)) continue;   /* not an algebraic integer? skip */
            bool inplane = true;
            for (int i = 2; i < n; i++) if (mpz_sgn(coord[i]) != 0) { inplane = false; break; }
            if (!inplane) continue;                                /* not of the form x - theta y */

            /* x = coord[0], y = -coord[1] */
            mpz_neg(yneg, coord[1]);
            eval_form(form, n, coord[0], yneg, fval);
            if (mpz_cmp(fval, m) == 0) solbuf_push(&sb, coord[0], yneg);
        }
      }
    }

    /* Small-|Y| gap (|Y| <= Y2p): the exponent enumeration above proves the
     * bound only for |Y| > Y2p, so this closes the gap.  Two exact methods,
     * chosen per the x-window width (see the feasibility gate):
     *  - NARROW window: scan x in [-Xmax, Xmax] (as before -- cheapest when the
     *    window is tiny, and Xmax is a proven x-bound for |y| <= Y2p);
     *  - WIDE window: for each fixed y the integer x with F(x,y)==m are EXACTLY
     *    the integer roots of the univariate F(x,y)-m, found by factorisation
     *    (magnitude-independent, so no O(Xmax) blowup).
     * Both are exact; union with the exponent orbit = the complete set (solbuf
     * dedups). */
    if (Y2p >= 0) {
        bool narrowX = (2.0 * Xmax + 1.0) <= 256.0;
        mpz_t X, Y; mpz_init(X); mpz_init(Y);
        if (narrowX) {
            for (long yy = -Y2p; yy <= Y2p; yy++) {
                mpz_set_si(Y, yy);
                for (long xx = -Xmax; xx <= Xmax; xx++) {
                    mpz_set_si(X, xx);
                    eval_form(form, n, X, Y, fval);
                    if (mpz_cmp(fval, m) == 0) solbuf_push(&sb, X, Y);
                }
            }
        } else {
            fmpz_poly_t Px; fmpz_poly_init(Px);
            mpz_t* xroots = NULL; int nxr = 0, xcap = 0;
            for (long yy = -Y2p; yy <= Y2p; yy++) {
                mpz_set_si(Y, yy);
                thue_form_xpoly(form, n, Y, m, Px);
                nxr = 0;
                fmpz_poly_int_roots(Px, &xroots, &nxr, &xcap);
                for (int i = 0; i < nxr; i++) {
                    eval_form(form, n, xroots[i], Y, fval);      /* keep the exact F==m contract */
                    if (mpz_cmp(fval, m) == 0) solbuf_push(&sb, xroots[i], Y);
                    mpz_clear(xroots[i]);
                }
            }
            free(xroots);
            fmpz_poly_clear(Px);
        }
        mpz_clear(X); mpz_clear(Y);
    }

    if (edbg) fprintf(stderr, "[thue-prof] enumerate+brute (B=%d, %d sols): %.1fms\n", B, sb.n, (clock()-tprof)*1000.0/CLOCKS_PER_SEC);

    /* package results */
    ThueSol* res = malloc(sizeof(ThueSol) * (size_t)(sb.n > 0 ? sb.n : 1));
    for (int i = 0; i < sb.n; i++) { mpz_init_set(res[i].x, sb.x[i]); mpz_init_set(res[i].y, sb.y[i]); }
    *out = res; rc = sb.n;

    /* cleanup */
    for (int i = 0; i < sb.n; i++) { mpz_clear(sb.x[i]); mpz_clear(sb.y[i]); }
    free(sb.x); free(sb.y); free(b);
    for (int i = 0; i < n; i++) mpz_clear(coord[i]);
    free(coord);
    mpz_clear(yneg); mpz_clear(fval);
    nf_elem_clear(u, K->nf); nf_elem_clear(ubeta, K->nf); nf_elem_clear(usign, K->nf);
    for (int i = 0; i < nmu; i++) nf_elem_clear(Mel[i], K->nf);
    free(Mel);
    thue_free_reps(Mreps, nM);
    for (int k = 0; k < r; k++) {
        for (long j = 0; j < W; j++) nf_elem_clear(Epow[k][j], K->nf);
        free(Epow[k]);
        nf_elem_clear(E[k], K->nf);
    }
    free(Epow); free(E);
    nf_units_free(U);
    nf_field_free(K);

done:
    for (int i = 0; i <= n; i++) mpz_clear(form[i]);
    free(form);
    mpz_clear(m);
    return rc;
}

int thue_solve_binary_form_bounded(const mpz_t* form, int n, const mpz_t m,
                                   int bound, ThueSol** out) {
    return thue_enumerate(form, n, m, bound, out);
}

int thue_solve_binary_form(const mpz_t* form, int n, const mpz_t m, ThueSol** out) {
    /* Reducible F(x,1) is not a Thue equation (the TdW method needs F
     * irreducible); when F factors into >= 2 coprime forms the solution set is
     * finite and found by the factorisation solver. */
    { fmpz_poly_t g; fmpz_poly_init(g);
      fmpz_t z; fmpz_init(z);
      for (int i = 0; i <= n; i++) { fmpz_set_mpz(z, form[n - i]); fmpz_poly_set_coeff_fmpz(g, i, z); }
      fmpz_clear(z);
      bool reducible = solvethue_poly_reducible(g);
      fmpz_poly_clear(g);
      if (reducible) return thue_solve_reducible_form(form, n, m, out);
    }
    /* Rigorous path: enumerate with the internally-computed Baker/de-Weger
     * bound.  Until that bound lands, DECLINE cheaply -- never an unproven
     * (incomplete) answer, and no number field is built per call. */
    if (!thue_bound_available()) return -1;
    return thue_enumerate(form, n, m, THUE_BOUND_RIGOROUS, out);
}

/* ---- test builtin: Solve`ThueSolveForm[{a0,..,an}, m, bound] ---- */

static Expr* builtin_thue_solve_form(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    Expr* coeffs = res->data.function.args[0];
    Expr* mE     = res->data.function.args[1];
    Expr* bE     = res->data.function.args[2];
    if (coeffs->type != EXPR_FUNCTION) return NULL;
    if (coeffs->data.function.head->type != EXPR_SYMBOL ||
        coeffs->data.function.head->data.symbol.name != SYM_List) return NULL;
    if (mE->type != EXPR_INTEGER || bE->type != EXPR_INTEGER) return NULL;

    int n = (int)coeffs->data.function.arg_count - 1;
    if (n < 3) return NULL;
    mpz_t* form = malloc(sizeof(mpz_t) * (size_t)(n + 1));
    for (int i = 0; i <= n; i++) {
        Expr* ci = coeffs->data.function.args[i];
        mpz_init(form[i]);
        if (ci->type == EXPR_INTEGER) mpz_set_si(form[i], (long)ci->data.integer);
        else if (ci->type == EXPR_BIGINT) mpz_set(form[i], ci->data.bigint);
        else { for (int k = 0; k <= i; k++) mpz_clear(form[k]); free(form); return NULL; }
    }
    mpz_t m; mpz_init_set_si(m, (long)mE->data.integer);
    int bound = (int)bE->data.integer;

    ThueSol* sols = NULL;
    int cnt = thue_solve_binary_form_bounded((const mpz_t*)form, n, m, bound, &sols);

    for (int i = 0; i <= n; i++) mpz_clear(form[i]);
    free(form); mpz_clear(m);

    if (cnt < 0) return NULL;   /* decline */

    Expr** items = (cnt > 0) ? malloc(sizeof(Expr*) * (size_t)cnt) : NULL;
    for (int i = 0; i < cnt; i++) {
        Expr* pr[2] = { expr_new_bigint_from_mpz(sols[i].x), expr_new_bigint_from_mpz(sols[i].y) };
        pr[0] = expr_bigint_normalize(pr[0]); pr[1] = expr_bigint_normalize(pr[1]);
        items[i] = expr_new_function(expr_new_symbol("List"), pr, 2);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), items, (size_t)cnt);
    free(items);
    thue_sols_free(sols, cnt);
    return out;
}

void solvethue_init(void) {
    symtab_add_builtin("Solve`ThueSolveForm", builtin_thue_solve_form);
    symtab_set_docstring("Solve`ThueSolveForm",
        "Solve`ThueSolveForm[{a0,...,an}, m, bound]\n"
        "\tInternal/testing: integer solutions {x,y} of the binary-form Thue\n"
        "\tequation Sum a_j x^(n-j) y^j == m, enumerating unit exponents up to\n"
        "\t|b_k| <= bound.  Complete only when `bound` provably covers all\n"
        "\tsolutions; used to validate the reconstruction engine.");
}

#else  /* !USE_FLINT */

void thue_sols_free(ThueSol* sols, int count) { (void)sols; (void)count; }
int thue_solve_binary_form(const mpz_t* form, int n, const mpz_t m, ThueSol** out) {
    (void)form; (void)n; (void)m; (void)out; return -1;
}
int thue_solve_binary_form_bounded(const mpz_t* form, int n, const mpz_t m,
                                   int bound, ThueSol** out) {
    (void)form; (void)n; (void)m; (void)bound; (void)out; return -1;
}
void solvethue_init(void) {}

#endif /* USE_FLINT */
