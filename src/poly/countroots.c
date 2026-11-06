/* Mathilda -- CountRoots[poly, x] and CountRoots[poly, {x, a, b}]: the number
 * of roots of a univariate polynomial, counted with multiplicity.
 *
 *   CountRoots[p, x]            real roots
 *   CountRoots[p, {x, a, b}]    real a, b (or +-Infinity): roots in [a, b];
 *                               complex a, b: roots in the closed rectangle
 *                               with opposite corners a and b
 *
 * The answer is exact. The coefficients must be rational numbers (an
 * approximate coefficient stands for its exact binary value); a Gaussian
 * rational coefficient is allowed for a real count. Anything else -- a second
 * variable (CountRoots::nuniv), an irrational or non-polynomial term
 * (CountRoots::nmet) -- stays unevaluated.
 *
 * METHOD. p is split into its squarefree decomposition p = c prod f_i^e_i
 * (FLINT), so the f_i are squarefree and pairwise coprime; each contributes
 * e_i times its number of distinct roots.
 *
 *   real interval   Sturm's theorem on f_i: the number of distinct roots in
 *                   (a, b] is V(a) - V(b), V the sign changes of the Sturm
 *                   sequence f, f', -rem(f, f'), ... (computed as primitive
 *                   pseudo-remainders with the sign of lc^delta restored);
 *                   one more when f(a) = 0. Endpoints are exact rationals or
 *                   infinite (the sign of the leading terms).
 *   complex box     the roots on the boundary are counted EXACTLY: on the edge
 *                   Im z = c (c rational), f(t + I c) = U(t) + I V(t) with U,
 *                   V in Q[t], and the roots there are the real roots of
 *                   gcd(U, V), counted on the edge by Sturm (likewise for the
 *                   vertical edges; a root at a corner is on two edges and is
 *                   counted once). The interior count comes from certified
 *                   isolating balls of every complex root (arb_fmpz_poly_
 *                   complex_roots): a ball strictly inside is inside; a ball
 *                   meeting the boundary is either a boundary root or not
 *                   yet resolved, and the precision is raised until the
 *                   number of boundary-meeting balls equals the exact number
 *                   of boundary roots.
 *
 * Attributes: Protected. Surfaces: symbolic (a polynomial in, an Integer out),
 * so no NDArray / packed / Compile[] kernel; the audits' EXEMPT lists record
 * this.
 *
 * Memory contract: takes ownership of res but never frees it; returns NULL
 * (unevaluated) or an owned result.
 */
#include "countroots.h"

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "message.h"
#include "print.h"
#include "sym_names.h"
#include "symtab.h"

#if defined(USE_FLINT) && defined(USE_MPFR)
#include <gmp.h>
#include <mpfr.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
#include <flint/fmpq_poly.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/arb_fmpz_poly.h>

/* ---- Expr -> exact numbers ------------------------------------------------ */

/* An exact rational from a real number leaf (an approximate one by its exact
 * binary value). */
static bool cr_real_to_fmpq(const Expr* e, fmpq_t out) {
    switch (e->type) {
        case EXPR_INTEGER: fmpq_set_si(out, e->data.integer, 1); return true;
        case EXPR_BIGINT:
            fmpz_set_mpz(fmpq_numref(out), e->data.bigint);
            fmpz_one(fmpq_denref(out));
            return true;
        case EXPR_REAL: {
            if (!isfinite(e->data.real)) return false;
            mpq_t q; mpq_init(q);
            mpq_set_d(q, e->data.real);
            fmpq_set_mpq(out, q);
            mpq_clear(q);
            return true;
        }
        case EXPR_MPFR: {
            if (!mpfr_number_p(e->data.mpfr)) return false;
            mpq_t q; mpq_init(q);
            mpfr_get_q(q, e->data.mpfr);
            fmpq_set_mpq(out, q);
            mpq_clear(q);
            return true;
        }
        case EXPR_FUNCTION: {
            const Expr* h = e->data.function.head;
            if (h->type != EXPR_SYMBOL || h->data.symbol.name != SYM_Rational ||
                e->data.function.arg_count != 2) return false;
            fmpq_t a, b; fmpq_init(a); fmpq_init(b);
            bool ok = cr_real_to_fmpq(e->data.function.args[0], a) &&
                      cr_real_to_fmpq(e->data.function.args[1], b) && !fmpq_is_zero(b);
            if (ok) fmpq_div(out, a, b);
            fmpq_clear(a); fmpq_clear(b);
            return ok;
        }
        default:
            return false;
    }
}

/* A Gaussian rational (re, im) from a number leaf or Complex[re, im]. */
static bool cr_to_gauss(const Expr* e, fmpq_t re, fmpq_t im) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Complex &&
        e->data.function.arg_count == 2)
        return cr_real_to_fmpq(e->data.function.args[0], re) &&
               cr_real_to_fmpq(e->data.function.args[1], im);
    fmpq_zero(im);
    return cr_real_to_fmpq(e, re);
}

/* An interval endpoint: a finite Gaussian rational, or +-Infinity (inf = +1
 * or -1). */
static bool cr_endpoint(const Expr* e, fmpq_t re, fmpq_t im, int* inf) {
    *inf = 0;
    if (e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_Infinity) { *inf = 1; return true; }
    /* -Infinity is Times[-1, Infinity] */
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Times &&
        e->data.function.arg_count == 2 &&
        e->data.function.args[0]->type == EXPR_INTEGER &&
        e->data.function.args[0]->data.integer == -1 &&
        e->data.function.args[1]->type == EXPR_SYMBOL &&
        e->data.function.args[1]->data.symbol.name == SYM_Infinity) {
        *inf = -1;
        return true;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_DirectedInfinity &&
        e->data.function.arg_count == 1) {
        const Expr* d = e->data.function.args[0];
        if (d->type == EXPR_INTEGER && (d->data.integer == 1 || d->data.integer == -1)) {
            *inf = (int)d->data.integer;
            return true;
        }
        return false;
    }
    return cr_to_gauss(e, re, im);
}

/* ---- Sturm counting --------------------------------------------------------- */

/* Sign of f at a point: x finite (inf == 0) or +-infinity. */
static int cr_sign_at(const fmpz_poly_t f, const fmpq_t x, int inf) {
    slong d = fmpz_poly_degree(f);
    if (d < 0) return 0;
    if (inf != 0) {
        int s = fmpz_sgn(fmpz_poly_lead(f));
        return (inf < 0 && (d & 1)) ? -s : s;
    }
    fmpq_t v; fmpq_init(v);
    fmpz_poly_evaluate_fmpq(v, f, x);
    int s = fmpq_sgn(v);
    fmpq_clear(v);
    return s;
}

typedef struct { fmpz_poly_struct* p; slong n; } CrSturm;

/* Divide by the (positive) content, keeping the sign -- FLINT's
 * fmpz_poly_primitive_part normalises the leading coefficient positive, which
 * would break the sign pattern a Sturm sequence is made of. */
static void cr_div_content(fmpz_poly_t r) {
    fmpz_t c;
    fmpz_init(c);
    fmpz_poly_content(c, r);
    if (!fmpz_is_zero(c) && !fmpz_is_one(c)) fmpz_poly_scalar_divexact_fmpz(r, r, c);
    fmpz_clear(c);
}

/* The Sturm sequence of a squarefree f of degree >= 1. */
static void cr_sturm_init(CrSturm* S, const fmpz_poly_t f) {
    slong cap = fmpz_poly_degree(f) + 2;
    S->p = (fmpz_poly_struct*)flint_malloc((size_t)cap * sizeof(fmpz_poly_struct));
    for (slong i = 0; i < cap; i++) fmpz_poly_init(S->p + i);
    fmpz_poly_set(S->p + 0, f);
    fmpz_poly_derivative(S->p + 1, f);
    cr_div_content(S->p + 1);
    S->n = 2;
    while (fmpz_poly_degree(S->p + S->n - 1) > 0 && S->n < cap) {
        fmpz_poly_struct* a = S->p + S->n - 2;
        fmpz_poly_struct* b = S->p + S->n - 1;
        fmpz_poly_struct* r = S->p + S->n;
        ulong dexp;
        fmpz_poly_pseudo_rem(r, &dexp, a, b);
        if (fmpz_poly_is_zero(r)) break;
        /* lc(b)^d a = Q b + R, so -rem = -R / lc(b)^d: keep the sign. */
        bool neg_lc = fmpz_sgn(fmpz_poly_lead(b)) < 0 && (dexp & 1);
        if (!neg_lc) fmpz_poly_neg(r, r);
        cr_div_content(r);
        S->n++;
    }
}

static void cr_sturm_clear(CrSturm* S) {
    slong cap = fmpz_poly_degree(S->p + 0) + 2;
    for (slong i = 0; i < cap; i++) fmpz_poly_clear(S->p + i);
    flint_free(S->p);
}

static slong cr_sturm_V(const CrSturm* S, const fmpq_t x, int inf) {
    slong v = 0;
    int last = 0;
    for (slong i = 0; i < S->n; i++) {
        int s = cr_sign_at(S->p + i, x, inf);
        if (s == 0) continue;
        if (last != 0 && s != last) v++;
        last = s;
    }
    return v;
}

/* Distinct real roots of a squarefree f in [lo, hi] (each end finite or
 * infinite; a finite end is included). */
static slong cr_count_real(const fmpz_poly_t f, const fmpq_t lo, int lo_inf,
                           const fmpq_t hi, int hi_inf) {
    if (fmpz_poly_degree(f) < 1) return 0;
    CrSturm S;
    cr_sturm_init(&S, f);
    slong c = cr_sturm_V(&S, lo, lo_inf) - cr_sturm_V(&S, hi, hi_inf);
    if (lo_inf == 0 && cr_sign_at(f, lo, 0) == 0) c++;
    cr_sturm_clear(&S);
    return c;
}

/* The squarefree part of f (primitive). */
static void cr_squarefree(fmpz_poly_t out, const fmpz_poly_t f) {
    fmpz_poly_t g, d;
    fmpz_poly_init(g); fmpz_poly_init(d);
    fmpz_poly_derivative(d, f);
    fmpz_poly_gcd(g, f, d);
    fmpz_poly_div(out, f, g);
    fmpz_poly_primitive_part(out, out);
    fmpz_poly_clear(g); fmpz_poly_clear(d);
}

/* ---- the complex box ----------------------------------------------------------- */

/* f(t + I c) (horizontal, vert = 0) or f(c + I t) (vertical, vert = 1) as
 * U + I V with U, V in Q[t]. */
static void cr_restrict_line(fmpq_poly_t U, fmpq_poly_t V, const fmpz_poly_t f,
                             const fmpq_t c, int vert) {
    fmpq_poly_t t1, t2, a, b;
    fmpq_poly_init(t1); fmpq_poly_init(t2); fmpq_poly_init(a); fmpq_poly_init(b);
    fmpq_poly_zero(U); fmpq_poly_zero(V);
    /* z = a + I b with a, b in Q[t]. */
    if (!vert) { fmpq_poly_set_coeff_si(a, 1, 1); fmpq_poly_set_coeff_fmpq(b, 0, c); }
    else       { fmpq_poly_set_coeff_fmpq(a, 0, c); fmpq_poly_set_coeff_si(b, 1, 1); }
    fmpz_t co; fmpz_init(co);
    for (slong i = fmpz_poly_degree(f); i >= 0; i--) {
        /* (U + I V)(a + I b) = (U a - V b) + I (U b + V a), then + f_i. */
        fmpq_poly_mul(t1, U, a);
        fmpq_poly_mul(t2, V, b);
        fmpq_poly_sub(t1, t1, t2);
        fmpq_poly_mul(t2, U, b);
        fmpq_poly_mul(V, V, a);
        fmpq_poly_add(V, V, t2);
        fmpq_poly_set(U, t1);
        fmpz_poly_get_coeff_fmpz(co, f, i);
        fmpq_poly_zero(t1);
        fmpq_poly_set_coeff_fmpz(t1, 0, co);
        fmpq_poly_add(U, U, t1);
    }
    fmpz_clear(co);
    fmpq_poly_clear(t1); fmpq_poly_clear(t2); fmpq_poly_clear(a); fmpq_poly_clear(b);
}

/* Distinct roots of f on the closed edge {t + I c : lo <= t <= hi} (vert 0)
 * or {c + I t : lo <= t <= hi} (vert 1). */
static slong cr_count_edge(const fmpz_poly_t f, const fmpq_t c, int vert,
                           const fmpq_t lo, const fmpq_t hi) {
    fmpq_poly_t U, V, G;
    fmpq_poly_init(U); fmpq_poly_init(V); fmpq_poly_init(G);
    cr_restrict_line(U, V, f, c, vert);
    fmpq_poly_gcd(G, U, V);
    slong n = 0;
    if (fmpq_poly_degree(G) >= 1) {
        fmpz_poly_t g, gs;
        fmpz_poly_init(g); fmpz_poly_init(gs);
        fmpq_poly_get_numerator(g, G);
        cr_squarefree(gs, g);
        n = cr_count_real(gs, lo, 0, hi, 0);
        fmpz_poly_clear(g); fmpz_poly_clear(gs);
    }
    fmpq_poly_clear(U); fmpq_poly_clear(V); fmpq_poly_clear(G);
    return n;
}

/* f(x + I y) == 0 exactly. */
static bool cr_zero_at(const fmpz_poly_t f, const fmpq_t x, const fmpq_t y) {
    fmpq_poly_t U, V;
    fmpq_poly_init(U); fmpq_poly_init(V);
    cr_restrict_line(U, V, f, y, 0);
    fmpq_t u, v; fmpq_init(u); fmpq_init(v);
    fmpq_poly_evaluate_fmpq(u, U, x);
    fmpq_poly_evaluate_fmpq(v, V, x);
    bool z = fmpq_is_zero(u) && fmpq_is_zero(v);
    fmpq_clear(u); fmpq_clear(v);
    fmpq_poly_clear(U); fmpq_poly_clear(V);
    return z;
}

/* Distinct roots of the squarefree f in the closed box [xl, xh] x [yl, yh]
 * (xl <= xh, yl <= yh, not both degenerate). -1 when the isolation does not
 * resolve within the precision cap. */
static slong cr_count_box(const fmpz_poly_t f, const fmpq_t xl, const fmpq_t xh,
                          const fmpq_t yl, const fmpq_t yh) {
    slong deg = fmpz_poly_degree(f);
    if (deg < 1) return 0;
    bool flat_x = fmpq_equal(xl, xh), flat_y = fmpq_equal(yl, yh);
    if (flat_y) return cr_count_edge(f, yl, 0, xl, xh);
    if (flat_x) return cr_count_edge(f, xl, 1, yl, yh);
    /* Exact boundary count: four closed edges, minus the corners counted
     * twice. */
    slong B = cr_count_edge(f, yl, 0, xl, xh) + cr_count_edge(f, yh, 0, xl, xh) +
              cr_count_edge(f, xl, 1, yl, yh) + cr_count_edge(f, xh, 1, yl, yh);
    const fmpq* cx[2] = { xl, xh };
    const fmpq* cy[2] = { yl, yh };
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            if (cr_zero_at(f, cx[i], cy[j])) B--;

    acb_ptr roots = _acb_vec_init(deg);
    slong result = -1;
    for (slong prec = 64; prec <= 16384; prec *= 2) {
        arb_fmpz_poly_complex_roots(roots, f, 0, prec);
        arb_t bxl, bxh, byl, byh;
        arb_init(bxl); arb_init(bxh); arb_init(byl); arb_init(byh);
        arb_set_fmpq(bxl, xl, prec + 32); arb_set_fmpq(bxh, xh, prec + 32);
        arb_set_fmpq(byl, yl, prec + 32); arb_set_fmpq(byh, yh, prec + 32);
        slong inside = 0, amb = 0;
        for (slong k = 0; k < deg; k++) {
            const arb_struct* re = acb_realref(roots + k);
            const arb_struct* im = acb_imagref(roots + k);
            bool in = arb_gt(re, bxl) && arb_lt(re, bxh) && arb_gt(im, byl) && arb_lt(im, byh);
            bool out = arb_lt(re, bxl) || arb_gt(re, bxh) || arb_lt(im, byl) || arb_gt(im, byh);
            if (in) inside++;
            else if (!out) amb++;
        }
        arb_clear(bxl); arb_clear(bxh); arb_clear(byl); arb_clear(byh);
        if (amb == B) { result = inside + B; break; }
    }
    _acb_vec_clear(roots, deg);
    return result;
}

/* ---- the builtin ----------------------------------------------------------------- */

static void cr_msg_nuniv(const Expr* p, const Expr* x) {
    char* ps = expr_to_string((Expr*)p);
    char* xs = expr_to_string((Expr*)x);
    mth_message("CountRoots", "nuniv", "%s is not a univariate function of %s.",
                ps ? ps : "", xs ? xs : "");
    free(ps); free(xs);
}

static void cr_msg_nmet(void) {
    mth_message("CountRoots", "nmet",
                "Unable to compute the number of roots with the available methods.");
}

/* The coefficient list of p in x (evaluated), or NULL when p is not a
 * polynomial in x. */
static Expr* cr_coefficients(const Expr* p, const Expr* x) {
    Expr* a1[2] = { expr_copy((Expr*)p), expr_copy((Expr*)x) };
    Expr* pq = eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ), a1, 2));
    bool poly = pq && pq->type == EXPR_SYMBOL && pq->data.symbol.name == SYM_True;
    if (pq) expr_free(pq);
    if (!poly) return NULL;
    Expr* a2[2] = { expr_copy((Expr*)p), expr_copy((Expr*)x) };
    Expr* cl = eval_and_free(expr_new_function(expr_new_symbol(SYM_CoefficientList), a2, 2));
    if (cl && cl->type == EXPR_FUNCTION && cl->data.function.head->type == EXPR_SYMBOL &&
        cl->data.function.head->data.symbol.name == SYM_List)
        return cl;
    if (cl) expr_free(cl);
    return NULL;
}

/* Is e free of every symbol other than numeric constants? (A coefficient
 * that mentions another variable makes p multivariate.) */
static bool cr_mentions_symbol(const Expr* e) {
    if (e->type == EXPR_SYMBOL) {
        const char* n = e->data.symbol.name;
        return !(n == SYM_Pi || n == SYM_E || n == SYM_I || n == SYM_EulerGamma ||
                 n == SYM_GoldenRatio || n == SYM_Catalan || n == SYM_Degree);
    }
    if (e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (cr_mentions_symbol(e->data.function.args[i])) return true;
    return false;
}

Expr* builtin_countroots(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    const Expr* p = res->data.function.args[0];
    const Expr* spec = res->data.function.args[1];
    const Expr* x = spec;
    const Expr *ea = NULL, *eb = NULL;
    if (spec->type == EXPR_FUNCTION && spec->data.function.head->type == EXPR_SYMBOL &&
        spec->data.function.head->data.symbol.name == SYM_List) {
        if (spec->data.function.arg_count != 3) return NULL;
        x = spec->data.function.args[0];
        ea = spec->data.function.args[1];
        eb = spec->data.function.args[2];
    }
    if (x->type != EXPR_SYMBOL) return NULL;

    /* The interval. */
    fmpq_t are, aim, bre, bim;
    fmpq_init(are); fmpq_init(aim); fmpq_init(bre); fmpq_init(bim);
    int ainf = -1, binf = 1;
    Expr* out = NULL;
    Expr* cl = NULL;
    fmpq_poly_t P, Q;            /* Re and Im parts of p */
    fmpq_poly_init(P); fmpq_poly_init(Q);
    if (ea) {
        if (!cr_endpoint(ea, are, aim, &ainf) || !cr_endpoint(eb, bre, bim, &binf)) goto done;
    }
    bool complex_box = ea && ainf == 0 && binf == 0 &&
                       (!fmpq_is_zero(aim) || !fmpq_is_zero(bim));
    if (ea && (ainf != 0 || binf != 0) && (!fmpq_is_zero(aim) || !fmpq_is_zero(bim)))
        goto done;   /* an infinite end with a complex other end */

    /* The polynomial. */
    cl = cr_coefficients(p, x);
    if (!cl) {
        cr_msg_nmet();   /* not a polynomial in x: Sin[x], 1/x, Sqrt[x] */
        goto done;
    }
    {
        size_t n = cl->data.function.arg_count;
        fmpq_t cre, cim;
        fmpq_init(cre); fmpq_init(cim);
        bool ok = true;
        for (size_t i = 0; i < n && ok; i++) {
            const Expr* c = cl->data.function.args[i];
            if (!cr_to_gauss(c, cre, cim)) {
                if (cr_mentions_symbol(c)) cr_msg_nuniv(p, x);
                else cr_msg_nmet();
                ok = false;
                break;
            }
            fmpq_poly_set_coeff_fmpq(P, (slong)i, cre);
            fmpq_poly_set_coeff_fmpq(Q, (slong)i, cim);
        }
        fmpq_clear(cre); fmpq_clear(cim);
        if (!ok) goto done;
    }
    if (!fmpq_poly_is_zero(Q)) {
        /* Gaussian coefficients: the real roots are those of gcd(Re p, Im p);
         * a complex box would need the complex polynomial itself. */
        if (complex_box) { cr_msg_nmet(); goto done; }
        fmpq_poly_gcd(P, P, Q);
    }
    if (fmpq_poly_degree(P) < 1) { out = expr_new_integer(0); goto done; }

    {
        /* Normalise the interval: ordered ends. */
        int lo_inf = ea ? ainf : -1, hi_inf = ea ? binf : 1;
        if (ea && !complex_box) {
            bool swap = (lo_inf > hi_inf) ||
                        (lo_inf == 0 && hi_inf == 0 && fmpq_cmp(are, bre) > 0);
            if (swap) {
                fmpq_swap(are, bre);
                int t = lo_inf; lo_inf = hi_inf; hi_inf = t;
            }
            if (lo_inf == hi_inf && lo_inf != 0) { out = expr_new_integer(0); goto done; }
        }
        fmpq_t xl, xh, yl, yh;
        fmpq_init(xl); fmpq_init(xh); fmpq_init(yl); fmpq_init(yh);
        if (complex_box) {
            if (fmpq_cmp(are, bre) <= 0) { fmpq_set(xl, are); fmpq_set(xh, bre); }
            else { fmpq_set(xl, bre); fmpq_set(xh, are); }
            if (fmpq_cmp(aim, bim) <= 0) { fmpq_set(yl, aim); fmpq_set(yh, bim); }
            else { fmpq_set(yl, bim); fmpq_set(yh, aim); }
        }
        fmpz_poly_t F;
        fmpz_poly_init(F);
        fmpq_poly_get_numerator(F, P);
        fmpz_poly_factor_t fac;
        fmpz_poly_factor_init(fac);
        fmpz_poly_factor_squarefree(fac, F);
        slong total = 0;
        bool ok = true;
        for (slong i = 0; i < fac->num && ok; i++) {
            slong c;
            if (complex_box) c = cr_count_box(fac->p + i, xl, xh, yl, yh);
            else c = cr_count_real(fac->p + i, are, lo_inf, bre, hi_inf);
            if (c < 0) ok = false;
            else total += c * fac->exp[i];
        }
        fmpz_poly_factor_clear(fac);
        fmpz_poly_clear(F);
        fmpq_clear(xl); fmpq_clear(xh); fmpq_clear(yl); fmpq_clear(yh);
        if (ok) out = expr_new_integer(total);
        else cr_msg_nmet();
    }

done:
    if (cl) expr_free(cl);
    fmpq_poly_clear(P); fmpq_poly_clear(Q);
    fmpq_clear(are); fmpq_clear(aim); fmpq_clear(bre); fmpq_clear(bim);
    return out;
}

#else  /* !(USE_FLINT && USE_MPFR) */

Expr* builtin_countroots(Expr* res) {
    (void)res;
    return NULL;
}

#endif

void countroots_init(void) {
    symtab_add_builtin("CountRoots", builtin_countroots);
    symtab_get_def("CountRoots")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("CountRoots",
        "CountRoots[poly, x]\n"
        "\tThe number of real roots of the univariate polynomial poly, counted\n"
        "\twith multiplicity.\n"
        "CountRoots[poly, {x, a, b}]\n"
        "\tThe number of roots in the closed interval [a, b] (a, b real or\n"
        "\t+-Infinity), or in the closed rectangle with opposite corners a and b\n"
        "\twhen they are complex. Exact: Sturm sequences on the squarefree\n"
        "\tfactors; exact boundary counts and certified isolation for a\n"
        "\trectangle. Rational (or approximate) coefficients.");
}
