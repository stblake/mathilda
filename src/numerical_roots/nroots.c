/*
 * nroots.c — NRoots[lhs == rhs, var, opts]   (see nroots.h)
 *
 * Numerically finds every root of a univariate polynomial equation and returns
 * a disjunction of equations  var==r1 || var==r2 || ...  (a bare equation when
 * there is a single root), repeating identical equations for roots of
 * multiplicity > 1.
 *
 * Pipeline:
 *   1. parse  Equal[lhs, rhs]  and the variable; form poly = lhs - rhs.
 *   2. Expand, validate polynomial-in-var, extract coefficients.
 *   3. numericalise each coefficient to a complex MPFR value (working prec).
 *   4. strip a trailing x^m factor (exact zero roots), then dispatch the
 *      reduced polynomial to the selected engine (Aberth / CompanionMatrix /
 *      JenkinsTraub), which returns all roots with multiplicity.
 *   5. chop noise, sort canonically, cluster multiple roots to identical
 *      values, round to the target precision, assemble the disjunction.
 *
 * All numeric work is MPFR; without USE_MPFR, NRoots returns unevaluated.
 *
 * Memory contract: never frees `res`; returns a fresh Expr* or NULL.
 */

#include "nroots.h"
#include "nroots_internal.h"

#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "expr.h"
#include "eval.h"
#include "numeric.h"
#include "sym_names.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"     /* is_complex, make_complex, is_rational */
#include "expand.h"         /* expr_expand */
#include "poly/poly.h"      /* is_polynomial, get_degree_poly, get_all_coeffs_expanded */
#include "poly/zupoly.h"    /* exact integer polynomials for squarefree decomposition */
#include "linalg/eigen.h"   /* eigen_all_eigenvalues_real_mpfr */

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* ------------------------------------------------------------------ *
 *  Diagnostics                                                        *
 * ------------------------------------------------------------------ */
static void nr_warn(const char* tag, const char* fmt, ...) {
    if (arith_warnings_muted()) return;
    va_list ap;
    fprintf(stderr, "NRoots::%s: ", tag);
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
}

/* Method selector. */
typedef enum { NR_AUTO = 0, NR_ABERTH, NR_COMPANION, NR_JT } NrMethod;

#ifdef USE_MPFR

/* ================================================================== *
 *  Shared NrPoly utilities (declared in nroots_internal.h)
 * ================================================================== */

void nr_poly_eval(const NrPoly* p, const ncpx* z, ncpx* val, ncpx* der,
                  mpfr_prec_t wp) {
    int n = p->deg;
    ncpx b, d;
    ncpx_init(&b, wp);
    ncpx_init(&d, wp);
    ncpx_set(&b, &p->c[n]);
    ncpx_set_ui(&d, 0);
    for (int i = n - 1; i >= 0; i--) {
        /* d = d*z + b ; b = b*z + c[i]   (Horner fused with derivative). */
        ncpx_mul(&d, &d, z, wp);
        ncpx_add(&d, &d, &b);
        ncpx_mul(&b, &b, z, wp);
        ncpx_add(&b, &b, &p->c[i]);
    }
    ncpx_set(val, &b);
    if (der) ncpx_set(der, &d);
    ncpx_clear(&b);
    ncpx_clear(&d);
}

void nr_newton_polish(const NrPoly* p, ncpx* z, mpfr_prec_t wp, int max_iter) {
    ncpx val, der, step;
    mpfr_t an, az, bound, tol, mag;
    ncpx_init(&val, wp); ncpx_init(&der, wp); ncpx_init(&step, wp);
    mpfr_init2(an, wp); mpfr_init2(az, wp); mpfr_init2(bound, wp);
    mpfr_init2(tol, wp); mpfr_init2(mag, wp);

    /* tol = 2^-(wp-6) relative step floor. */
    mpfr_set_ui(tol, 1, MPFR_RNDN);
    mpfr_div_2si(tol, tol, (long)wp - 6, MPFR_RNDN);

    for (int it = 0; it < max_iter; it++) {
        nr_poly_eval(p, z, &val, &der, wp);
        ncpx_abs(mag, &der);
        if (mpfr_zero_p(mag)) break;
        ncpx_div(&step, &val, &der, wp);
        ncpx_sub(z, z, &step);

        ncpx_abs(an, &step);
        ncpx_abs(az, z);
        mpfr_set_ui(bound, 1, MPFR_RNDN);
        if (mpfr_cmp(az, bound) > 0) mpfr_set(bound, az, MPFR_RNDN);
        mpfr_mul(bound, bound, tol, MPFR_RNDN);
        if (mpfr_cmp(an, bound) <= 0) break;
    }

    ncpx_clear(&val); ncpx_clear(&der); ncpx_clear(&step);
    mpfr_clear(an); mpfr_clear(az); mpfr_clear(bound);
    mpfr_clear(tol); mpfr_clear(mag);
}

/* ------------------------------------------------------------------ *
 *  Numeric leaf -> ncpx
 * ------------------------------------------------------------------ */
static bool nr_real_to_mpfr(Expr* e, mpfr_t out) {
    if (!e) return false;
    int64_t n, d;
    switch (e->type) {
        case EXPR_INTEGER: mpfr_set_si(out, e->data.integer, MPFR_RNDN); return true;
        case EXPR_REAL:    mpfr_set_d (out, e->data.real,    MPFR_RNDN); return true;
        case EXPR_BIGINT:  mpfr_set_z (out, e->data.bigint,  MPFR_RNDN); return true;
        case EXPR_MPFR:    mpfr_set   (out, e->data.mpfr,    MPFR_RNDN); return true;
        case EXPR_FUNCTION:
            if (is_rational(e, &n, &d) && d != 0) {
                mpfr_set_si(out, n, MPFR_RNDN);
                mpfr_t dd; mpfr_init2(dd, mpfr_get_prec(out));
                mpfr_set_si(dd, d, MPFR_RNDN);
                mpfr_div(out, out, dd, MPFR_RNDN);
                mpfr_clear(dd);
                return true;
            }
            return false;
        default: return false;
    }
}

/* Convert an (already numericalised) Expr to a complex MPFR value. */
static bool nr_expr_to_ncpx(Expr* e, ncpx* out) {
    mpfr_set_zero(out->re, 1);
    mpfr_set_zero(out->im, 1);
    if (nr_real_to_mpfr(e, out->re)) return true;
    Expr *re, *im;
    if (is_complex(e, &re, &im)) {
        if (nr_real_to_mpfr(re, out->re) && nr_real_to_mpfr(im, out->im))
            return true;
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Coefficient extraction:  equation + var -> NrPoly at working prec
 *
 *  Returns the polynomial on success.  On a constant polynomial sets
 *  *constant_zero (true if the constant is numerically zero) and returns NULL
 *  with *was_constant = true.  On any error returns NULL with *was_constant
 *  = false (a diagnostic has been emitted).
 * ------------------------------------------------------------------ */
static NrPoly* nr_build_poly(Expr* equation, Expr* var, mpfr_prec_t wp,
                             bool* was_constant, bool* constant_zero,
                             ZUPoly** out_zpoly) {
    *was_constant = false;
    *constant_zero = false;
    if (out_zpoly) *out_zpoly = NULL;

    /* The equation is evaluated before NRoots sees it, so a numeric (in)equality
     * may already have collapsed:  1 == 0 -> False,  1 == 1 -> True. */
    if (equation && equation->type == EXPR_SYMBOL) {
        if (equation->data.symbol == SYM_False) { *was_constant = true; *constant_zero = false; return NULL; }
        if (equation->data.symbol == SYM_True)  { *was_constant = true; *constant_zero = true;  return NULL; }
    }

    if (!equation || equation->type != EXPR_FUNCTION
        || equation->data.function.head->type != EXPR_SYMBOL
        || equation->data.function.head->data.symbol != SYM_Equal
        || equation->data.function.arg_count != 2) {
        nr_warn("neqn", "first argument is not an equation lhs == rhs.");
        return NULL;
    }

    Expr* lhs = equation->data.function.args[0];
    Expr* rhs = equation->data.function.args[1];

    /* poly = lhs - rhs = Plus[lhs, Times[-1, rhs]]. */
    Expr* neg = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(rhs) }, 2);
    Expr* sub = expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy(lhs), neg }, 2);
    Expr* poly = eval_and_free(sub);
    Expr* expanded = expr_expand(poly);
    expr_free(poly);

    if (!is_polynomial(expanded, &var, 1)) {
        nr_warn("npoly", "the equation is not polynomial in the given variable.");
        expr_free(expanded);
        return NULL;
    }

    int d = get_degree_poly(expanded, var);
    NumericSpec spec; spec.mode = NUMERIC_MODE_MPFR; spec.bits = (long)wp;
    spec.preserve_inexact = false;

    if (d < 1) {
        /* Constant: decide zero vs non-zero numerically. */
        Expr* c0 = get_coeff_expanded(expanded, var, 0);
        Expr* num = eval_and_free(numericalize(c0, spec));
        expr_free(c0);
        ncpx v; ncpx_init(&v, wp);
        bool ok = nr_expr_to_ncpx(num, &v);
        expr_free(num);
        mpfr_t mag; mpfr_init2(mag, wp); ncpx_abs(mag, &v);
        *was_constant = true;
        *constant_zero = ok ? (mpfr_zero_p(mag) != 0) : false;
        mpfr_clear(mag); ncpx_clear(&v);
        expr_free(expanded);
        return NULL;
    }

    Expr** coeffs = NULL;
    if (!get_all_coeffs_expanded(expanded, var, d, &coeffs)) {
        nr_warn("npoly", "could not extract coefficients (variable not atomic).");
        expr_free(expanded);
        return NULL;
    }

    NrPoly* p = (NrPoly*)malloc(sizeof(NrPoly));
    p->deg = d;
    p->prec = wp;
    p->is_real = 1;
    p->c = (ncpx*)malloc(sizeof(ncpx) * (size_t)(d + 1));

    bool fail = false;
    int built = 0;
    for (int k = 0; k <= d; k++) {
        ncpx_init(&p->c[k], wp);
        built = k + 1;
        Expr* num = eval_and_free(numericalize(coeffs[k], spec));
        bool ok = nr_expr_to_ncpx(num, &p->c[k]);
        expr_free(num);
        if (!ok) { fail = true; k = d; /* stop loop */ }
        else if (!mpfr_zero_p(p->c[k].im)) p->is_real = 0;
    }
    for (int k = 0; k <= d; k++) expr_free(coeffs[k]);
    free(coeffs);
    /* Exact integer polynomial (if any) for the squarefree multiplicity path. */
    if (out_zpoly && !fail) *out_zpoly = expr_to_zupoly(expanded, var);
    expr_free(expanded);

    if (fail) {
        nr_warn("nnum", "polynomial has a non-numeric coefficient.");
        for (int k = 0; k < built; k++) ncpx_clear(&p->c[k]);
        free(p->c); free(p);
        return NULL;
    }
    return p;
}

static void nr_poly_free(NrPoly* p) {
    if (!p) return;
    for (int k = 0; k <= p->deg; k++) ncpx_clear(&p->c[k]);
    free(p->c);
    free(p);
}

/* ------------------------------------------------------------------ *
 *  Engine dispatch with trailing-zero (x^m) deflation.
 *  `roots` is caller-owned, length p->deg, every cell ncpx_init'd at wp.
 * ------------------------------------------------------------------ */
static int nr_solve(const NrPoly* p, NrMethod method, int max_iter,
                    mpfr_prec_t wp, ncpx* roots) {
    int d = p->deg;

    /* m = number of leading zero coefficients = multiplicity of the root 0. */
    int m = 0;
    {
        mpfr_t mag; mpfr_init2(mag, wp);
        for (; m < d; m++) {
            ncpx_abs(mag, &p->c[m]);
            if (!mpfr_zero_p(mag)) break;
        }
        mpfr_clear(mag);
    }

    /* Reduced polynomial shares the contiguous coefficient tail (no copy). */
    NrPoly r = *p;
    r.c = p->c + m;
    r.deg = d - m;

    int rc = 0;
    if (r.deg >= 1) {
        switch (method) {
            case NR_AUTO:
            case NR_ABERTH:
                rc = nr_aberth(&r, roots, max_iter, wp);
                break;
            case NR_COMPANION:
                rc = nr_companion(&r, roots, wp);
                break;
            case NR_JT:
                rc = nr_jenkinstraub(&r, roots, max_iter, wp);
                break;
        }
    }
    if (rc != 0) return rc;

    /* Append m exact-zero roots. */
    for (int i = r.deg; i < d; i++) ncpx_set_ui(&roots[i], 0);
    return 0;
}

/* ------------------------------------------------------------------ *
 *  Exact (integer-coefficient) path: squarefree decomposition so that
 *  multiple roots are solved as well-conditioned squarefree factors with
 *  known multiplicity, instead of as a catastrophically ill-conditioned
 *  high-degree polynomial (e.g. (x^2-2)^30).
 * ------------------------------------------------------------------ */

/* p'(x) as a ZUPoly (zupoly_derivative is not public). */
static ZUPoly* nr_zupoly_derivative(const ZUPoly* p) {
    if (p->deg <= 0) return zupoly_new(1);
    ZUPoly* q = zupoly_new(p->deg);
    mpz_t c; mpz_init(c);
    for (int i = 1; i <= p->deg; i++) {
        mpz_mul_ui(c, *zupoly_getcoef(p, i), (unsigned long)i);
        zupoly_setcoef(q, i - 1, c);
    }
    mpz_clear(c);
    return q;
}

/* Yun squarefree decomposition: p = prod fac[k]^mult[k] with each fac[k]
 * squarefree and pairwise coprime.  Returns 0 and fills the facs, mults and nf
 * outputs on success (caller frees each fac with zupoly_free, then the arrays);
 * returns -1 on any inexact division (caller falls back to the numeric path). */
static int nr_yun(const ZUPoly* p, ZUPoly*** facs_out, int** mults_out, int* nf_out) {
    ZUPoly* fp = nr_zupoly_derivative(p);
    ZUPoly* g  = zupoly_gcd(p, fp);
    zupoly_free(fp);
    ZUPoly* w = zupoly_divexact(p, g);
    if (!w) { zupoly_free(g); return -1; }
    ZUPoly* y = g;   /* ownership moves to y */

    int cap = 8, nf = 0;
    ZUPoly** facs = (ZUPoly**)malloc(sizeof(ZUPoly*) * (size_t)cap);
    int* mults = (int*)malloc(sizeof(int) * (size_t)cap);
    int rc = 0, i = 1;
    while (w->deg > 0) {
        ZUPoly* z  = zupoly_gcd(w, y);
        ZUPoly* fi = zupoly_divexact(w, z);
        ZUPoly* ny = zupoly_divexact(y, z);
        if (!fi || !ny) { zupoly_free(z); zupoly_free(fi); zupoly_free(ny); rc = -1; break; }
        if (fi->deg > 0) {
            if (nf == cap) { cap *= 2; facs = realloc(facs, sizeof(ZUPoly*) * (size_t)cap); mults = realloc(mults, sizeof(int) * (size_t)cap); }
            facs[nf] = zupoly_copy(fi); mults[nf] = i; nf++;
        }
        zupoly_free(fi);
        zupoly_free(w); w = z;
        zupoly_free(y); y = ny;
        i++;
        if (i > p->deg + 2) break;   /* safety */
    }
    zupoly_free(w); zupoly_free(y);
    if (rc != 0) { for (int k = 0; k < nf; k++) zupoly_free(facs[k]); free(facs); free(mults); return -1; }
    *facs_out = facs; *mults_out = mults; *nf_out = nf;
    return 0;
}

/* Convert an integer ZUPoly to a (real) NrPoly at working precision. */
static NrPoly* nr_zupoly_to_nrpoly(const ZUPoly* z, mpfr_prec_t wp) {
    NrPoly* p = (NrPoly*)malloc(sizeof(NrPoly));
    p->deg = z->deg; p->prec = wp; p->is_real = 1;
    p->c = (ncpx*)malloc(sizeof(ncpx) * (size_t)(z->deg + 1));
    for (int i = 0; i <= z->deg; i++) {
        ncpx_init(&p->c[i], wp);
        mpfr_set_z(p->c[i].re, *zupoly_getcoef(z, i), MPFR_RNDN);
        mpfr_set_zero(p->c[i].im, 1);
    }
    return p;
}

/* Solve an exact polynomial via squarefree decomposition.  Fills roots[0..d-1]
 * (d = zp->deg) with each squarefree factor's roots repeated by multiplicity.
 * Returns 0 on success, -1 to request the numeric fallback. */
static int nr_solve_exact(const ZUPoly* zp, NrMethod method, int max_iter,
                          mpfr_prec_t wp, ncpx* roots) {
    int d = zp->deg;
    ZUPoly** facs = NULL; int* mults = NULL; int nf = 0;
    if (nr_yun(zp, &facs, &mults, &nf) != 0) return -1;

    int fill = 0, rc = 0;
    for (int k = 0; k < nf && rc == 0; k++) {
        NrPoly* fp = nr_zupoly_to_nrpoly(facs[k], wp);
        int di = fp->deg;
        ncpx* fr = (ncpx*)malloc(sizeof(ncpx) * (size_t)di);
        for (int j = 0; j < di; j++) ncpx_init(&fr[j], wp);
        if (nr_solve(fp, method, max_iter, wp, fr) == 0) {
            for (int j = 0; j < di; j++)
                for (int m = 0; m < mults[k]; m++)
                    if (fill < d) ncpx_set(&roots[fill++], &fr[j]);
        } else rc = -1;
        for (int j = 0; j < di; j++) ncpx_clear(&fr[j]);
        free(fr);
        nr_poly_free(fp);
    }
    for (int k = 0; k < nf; k++) zupoly_free(facs[k]);
    free(facs); free(mults);
    if (rc == 0 && fill != d) rc = -1;
    return rc;
}

/* ------------------------------------------------------------------ *
 *  Post-processing: chop, canonical sort, multiplicity clustering
 * ------------------------------------------------------------------ */

/* Canonical ordering, mirroring root_numeric.c:
 *  reals (ascending), then complex by Re ascending, |Im| ascending,
 *  more-negative Im first. */
static int nr_canon_cmp(const ncpx* a, int areal, const ncpx* b, int breal,
                        mpfr_prec_t wp) {
    if (areal && !breal) return -1;
    if (!areal && breal) return  1;
    int c = mpfr_cmp(a->re, b->re);
    if (c != 0) return c < 0 ? -1 : 1;
    if (areal && breal) return 0;
    mpfr_t aa, bb; mpfr_init2(aa, wp); mpfr_init2(bb, wp);
    mpfr_abs(aa, a->im, MPFR_RNDN);
    mpfr_abs(bb, b->im, MPFR_RNDN);
    c = mpfr_cmp(aa, bb);
    mpfr_clear(aa); mpfr_clear(bb);
    if (c != 0) return c < 0 ? -1 : 1;
    int sa = mpfr_sgn(a->im), sb = mpfr_sgn(b->im);
    if (sa < sb) return -1;
    if (sa > sb) return  1;
    return 0;
}

/* |a - b| <= tol * max(1, |a|) ? */
static bool nr_close(const ncpx* a, const ncpx* b, const mpfr_t tol,
                     mpfr_prec_t wp) {
    ncpx diff; ncpx_init(&diff, wp);
    mpfr_t dm, am, bound;
    mpfr_init2(dm, wp); mpfr_init2(am, wp); mpfr_init2(bound, wp);
    ncpx_sub(&diff, a, b);
    ncpx_abs(dm, &diff);
    ncpx_abs(am, a);
    mpfr_set_ui(bound, 1, MPFR_RNDN);
    if (mpfr_cmp(am, bound) > 0) mpfr_set(bound, am, MPFR_RNDN);
    mpfr_mul(bound, bound, tol, MPFR_RNDN);
    bool ok = mpfr_cmp(dm, bound) <= 0;
    ncpx_clear(&diff);
    mpfr_clear(dm); mpfr_clear(am); mpfr_clear(bound);
    return ok;
}

/* Shallow swap of two ncpx cells (moves the mpfr limb handles). */
static void nr_swap_ncpx(ncpx* a, ncpx* b) { ncpx t = *a; *a = *b; *b = t; }

static void nr_postprocess(ncpx* roots, int* realf, int n, int poly_is_real,
                           mpfr_prec_t wp, long target_bits) {
    /* 1. Chop a component that is pure numerical noise relative to the other
     *    (for all polynomials): a near-real root prints as  a, a near-imaginary
     *    root as  b I, instead of carrying a 1e-30 partner.  realf records the
     *    real/complex split after chopping. */
    mpfr_t chop, scale, aim, are, bound;
    mpfr_init2(chop, wp); mpfr_init2(scale, wp);
    mpfr_init2(aim, wp); mpfr_init2(are, wp); mpfr_init2(bound, wp);
    mpfr_set_ui(chop, 1, MPFR_RNDN);
    mpfr_div_2si(chop, chop, target_bits / 2, MPFR_RNDN);  /* 2^-(target/2) */

    for (int i = 0; i < n; i++) {
        mpfr_abs(aim, roots[i].im, MPFR_RNDN);
        mpfr_abs(are, roots[i].re, MPFR_RNDN);
        /* chop Im if negligible vs Re. */
        mpfr_set_ui(bound, 1, MPFR_RNDN);
        if (mpfr_cmp(are, bound) > 0) mpfr_set(bound, are, MPFR_RNDN);
        mpfr_mul(bound, bound, chop, MPFR_RNDN);
        if (mpfr_cmp(aim, bound) <= 0) mpfr_set_zero(roots[i].im, 1);
        /* else chop Re if negligible vs Im. */
        else {
            mpfr_set_ui(bound, 1, MPFR_RNDN);
            if (mpfr_cmp(aim, bound) > 0) mpfr_set(bound, aim, MPFR_RNDN);
            mpfr_mul(bound, bound, chop, MPFR_RNDN);
            if (mpfr_cmp(are, bound) <= 0) mpfr_set_zero(roots[i].re, 1);
        }
        realf[i] = mpfr_zero_p(roots[i].im) ? 1 : 0;
    }

    /* 2. General pairwise clustering: group all roots within a relative
     *    tolerance (regardless of position) and snap each group to its
     *    centroid, so roots of multiplicity > 1 become numerically identical
     *    and print as identical equations.  Single-linkage anchored at each
     *    group's first member; n is small. */
    mpfr_t tolc, sre, sim;
    mpfr_init2(tolc, wp); mpfr_init2(sre, wp); mpfr_init2(sim, wp);
    {
        double td = numeric_bits_to_digits(target_bits);
        double t = pow(10.0, -0.5 * td);
        mpfr_set_d(tolc, t, MPFR_RNDN);
    }
    int* grp = (int*)malloc(sizeof(int) * (size_t)n);
    for (int i = 0; i < n; i++) grp[i] = -1;
    int ng = 0;
    for (int i = 0; i < n; i++) {
        if (grp[i] >= 0) continue;
        grp[i] = ng;
        for (int j = i + 1; j < n; j++)
            if (grp[j] < 0 && nr_close(&roots[j], &roots[i], tolc, wp)) grp[j] = ng;
        ng++;
    }
    for (int g = 0; g < ng; g++) {
        int cnt = 0, allreal = 1;
        mpfr_set_zero(sre, 1); mpfr_set_zero(sim, 1);
        for (int i = 0; i < n; i++) if (grp[i] == g) {
            mpfr_add(sre, sre, roots[i].re, MPFR_RNDN);
            mpfr_add(sim, sim, roots[i].im, MPFR_RNDN);
            if (!realf[i]) allreal = 0;
            cnt++;
        }
        if (cnt <= 1) continue;
        mpfr_div_ui(sre, sre, (unsigned long)cnt, MPFR_RNDN);
        mpfr_div_ui(sim, sim, (unsigned long)cnt, MPFR_RNDN);
        if (allreal) mpfr_set_zero(sim, 1);
        for (int i = 0; i < n; i++) if (grp[i] == g) {
            mpfr_set(roots[i].re, sre, MPFR_RNDN);
            mpfr_set(roots[i].im, sim, MPFR_RNDN);
            if (allreal) realf[i] = 1;
        }
    }
    free(grp);

    /* 3. Conjugate-pair symmetrization (real-coefficient polynomials only):
     *     complex roots come in exact conjugate pairs.  Average each pair's Re
     *     and |Im| so the two share an identical real part and ±Im, which both
     *     improves accuracy and makes ordering deterministic (the negative-Im
     *     root then sorts first, matching the Wolfram Language). */
    if (poly_is_real) {
        mpfr_t t1, t2, scl, tolp;
        mpfr_init2(t1, wp); mpfr_init2(t2, wp); mpfr_init2(scl, wp);
        mpfr_init2(tolp, wp);
        mpfr_set_ui(tolp, 1, MPFR_RNDN);
        mpfr_div_2si(tolp, tolp, target_bits / 2, MPFR_RNDN);   /* 2^-(target/2) */
        int* paired = (int*)calloc((size_t)n, sizeof(int));
        for (int a = 0; a < n; a++) {
            if (realf[a] || paired[a]) continue;
            int best = -1;
            for (int b = a + 1; b < n; b++) {
                if (realf[b] || paired[b]) continue;
                mpfr_sub(t1, roots[a].re, roots[b].re, MPFR_RNDN); mpfr_abs(t1, t1, MPFR_RNDN);
                mpfr_add(t2, roots[a].im, roots[b].im, MPFR_RNDN); mpfr_abs(t2, t2, MPFR_RNDN);
                mpfr_abs(scl, roots[a].re, MPFR_RNDN);
                mpfr_abs(are, roots[a].im, MPFR_RNDN);
                if (mpfr_cmp(are, scl) > 0) mpfr_set(scl, are, MPFR_RNDN);
                mpfr_set_ui(bound, 1, MPFR_RNDN);
                if (mpfr_cmp(scl, bound) > 0) mpfr_set(bound, scl, MPFR_RNDN);
                mpfr_mul(bound, bound, tolp, MPFR_RNDN);
                if (mpfr_cmp(t1, bound) <= 0 && mpfr_cmp(t2, bound) <= 0) { best = b; break; }
            }
            if (best >= 0) {
                int b = best;
                mpfr_add(t1, roots[a].re, roots[b].re, MPFR_RNDN); mpfr_div_ui(t1, t1, 2, MPFR_RNDN);
                mpfr_abs(t2, roots[a].im, MPFR_RNDN);
                mpfr_abs(scl, roots[b].im, MPFR_RNDN);
                mpfr_add(t2, t2, scl, MPFR_RNDN); mpfr_div_ui(t2, t2, 2, MPFR_RNDN);
                mpfr_set(roots[a].re, t1, MPFR_RNDN); mpfr_neg(roots[a].im, t2, MPFR_RNDN);
                mpfr_set(roots[b].re, t1, MPFR_RNDN); mpfr_set(roots[b].im, t2, MPFR_RNDN);
                paired[a] = paired[b] = 1;
            }
        }
        free(paired);
        mpfr_clear(t1); mpfr_clear(t2); mpfr_clear(scl); mpfr_clear(tolp);
    }

    /* 3b. Chop a noise real part on complex roots of real polynomials
     *     (a pure-imaginary root prints as  b I, not  1e-16 + b I). */
    if (poly_is_real) {
        for (int i = 0; i < n; i++) {
            if (realf[i]) continue;
            mpfr_abs(are, roots[i].re, MPFR_RNDN);
            mpfr_abs(aim, roots[i].im, MPFR_RNDN);
            mpfr_set_ui(bound, 1, MPFR_RNDN);
            if (mpfr_cmp(aim, bound) > 0) mpfr_set(bound, aim, MPFR_RNDN);
            mpfr_mul(bound, bound, chop, MPFR_RNDN);
            if (mpfr_cmp(are, bound) <= 0) mpfr_set_zero(roots[i].re, 1);
        }
    }

    /* 4. Canonical insertion sort (n small). */
    for (int i = 1; i < n; i++) {
        for (int j = i; j > 0; j--) {
            if (nr_canon_cmp(&roots[j], realf[j], &roots[j-1], realf[j-1], wp) >= 0)
                break;
            nr_swap_ncpx(&roots[j], &roots[j-1]);
            int t = realf[j]; realf[j] = realf[j-1]; realf[j-1] = t;
        }
    }

    mpfr_clear(chop); mpfr_clear(scale);
    mpfr_clear(aim); mpfr_clear(are); mpfr_clear(bound);
    mpfr_clear(tolc); mpfr_clear(sre); mpfr_clear(sim);
}

/* Build a numeric Expr from one root at the target precision. */
static Expr* nr_root_to_expr(const ncpx* z, int isreal, long target_bits,
                             bool want_machine) {
    if (isreal) {
        if (want_machine) return expr_new_real(mpfr_get_d(z->re, MPFR_RNDN));
        mpfr_t r; mpfr_init2(r, target_bits);
        mpfr_set(r, z->re, MPFR_RNDN);
        Expr* e = expr_new_mpfr_copy(r);
        mpfr_clear(r);
        return e;
    }
    if (want_machine) {
        return make_complex(expr_new_real(mpfr_get_d(z->re, MPFR_RNDN)),
                            expr_new_real(mpfr_get_d(z->im, MPFR_RNDN)));
    }
    mpfr_t r, im; mpfr_init2(r, target_bits); mpfr_init2(im, target_bits);
    mpfr_set(r, z->re, MPFR_RNDN);
    mpfr_set(im, z->im, MPFR_RNDN);
    Expr* e = make_complex(expr_new_mpfr_copy(r), expr_new_mpfr_copy(im));
    mpfr_clear(r); mpfr_clear(im);
    return e;
}

#endif /* USE_MPFR */

/* ================================================================== *
 *  Option parsing
 * ================================================================== */
typedef struct {
    NrMethod   method;
    const char* method_name;  /* requested name, for diagnostics; NULL = auto */
    double     prec_goal;     /* digits; -1 => Automatic                       */
    int        max_iter;      /* -1 => Automatic                               */
} NrOpts;

static bool nr_is_known_option(const char* s) {
    return s == SYM_Method || s == SYM_MaxIterations
        || s == SYM_PrecisionGoal || s == SYM_StepMonitor;
}

static bool nr_is_option_arg(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* lhs = e->data.function.args[0];
    return lhs->type == EXPR_SYMBOL && nr_is_known_option(lhs->data.symbol);
}

static NrMethod nr_method_from_name(const char* s) {
    if (!strcmp(s, "Aberth")) return NR_ABERTH;
    if (!strcmp(s, "CompanionMatrix")) return NR_COMPANION;
    if (!strcmp(s, "JenkinsTraub")) return NR_JT;
    return NR_AUTO;  /* unknown handled by caller via method_name */
}

static bool nr_to_double_real(Expr* e, double* out) {
    int64_t n, d;
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer; return true;
        case EXPR_REAL:    *out = e->data.real;            return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint); return true;
        case EXPR_FUNCTION:
            if (is_rational(e, &n, &d) && d != 0) { *out = (double)n/(double)d; return true; }
            return false;
        default: return false;
    }
}

/* Returns false (with a diagnostic) on an invalid option value. */
static bool nr_apply_option(Expr* rule, NrOpts* o) {
    Expr* lhs = rule->data.function.args[0];
    Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol;

    if (name == SYM_Method) {
        Expr* m = rhs;
        if (m->type == EXPR_FUNCTION && m->data.function.head->type == EXPR_SYMBOL
            && m->data.function.head->data.symbol == SYM_List
            && m->data.function.arg_count >= 1)
            m = m->data.function.args[0];
        if (m->type == EXPR_SYMBOL && m->data.symbol == SYM_Automatic) {
            o->method = NR_AUTO; o->method_name = NULL; return true;
        }
        if (m->type == EXPR_STRING) {
            o->method_name = m->data.string;
            o->method = nr_method_from_name(m->data.string);
            if (o->method == NR_AUTO) { nr_warn("bdmtd", "unknown Method \"%s\".", m->data.string); return false; }
            return true;
        }
        if (m->type == EXPR_SYMBOL) {
            o->method_name = m->data.symbol;
            o->method = nr_method_from_name(m->data.symbol);
            if (o->method == NR_AUTO) { nr_warn("bdmtd", "unknown Method %s.", m->data.symbol); return false; }
            return true;
        }
        nr_warn("bdmtd", "invalid Method value."); return false;
    }
    if (name == SYM_PrecisionGoal) {
        if (rhs->type == EXPR_SYMBOL
            && (rhs->data.symbol == SYM_Automatic || rhs->data.symbol == SYM_Infinity)) {
            o->prec_goal = -1.0; return true;
        }
        double v;
        if (nr_to_double_real(rhs, &v) && v > 0.0) { o->prec_goal = v; return true; }
        nr_warn("badopt", "invalid PrecisionGoal value."); return false;
    }
    if (name == SYM_MaxIterations) {
        if (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_Automatic) { o->max_iter = -1; return true; }
        if (rhs->type == EXPR_INTEGER && rhs->data.integer > 0) { o->max_iter = (int)rhs->data.integer; return true; }
        nr_warn("badopt", "MaxIterations must be a positive integer or Automatic."); return false;
    }
    if (name == SYM_StepMonitor) return true;   /* accepted; best-effort */
    return false;
}

/* ================================================================== *
 *  Builtin entry
 * ================================================================== */
Expr* builtin_nroots(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    /* Peel trailing options; exactly two positional args must remain. */
    size_t pos_end = argc;
    while (pos_end > 0 && nr_is_option_arg(res->data.function.args[pos_end - 1]))
        pos_end--;
    for (size_t i = pos_end; i < argc; i++) {
        if (!nr_is_option_arg(res->data.function.args[i])) {
            nr_warn("nopt", "argument %zu is not a valid option.", i + 1);
            return NULL;
        }
    }
    if (pos_end != 2) return NULL;

    Expr* equation = res->data.function.args[0];
    Expr* var      = res->data.function.args[1];
    if (!var || var->type != EXPR_SYMBOL) {
        nr_warn("ivar", "second argument is not a variable.");
        return NULL;
    }

    NrOpts opts = { NR_AUTO, NULL, -1.0, -1 };
    for (size_t i = pos_end; i < argc; i++) {
        if (!nr_apply_option(res->data.function.args[i], &opts)) return NULL;
    }

#ifndef USE_MPFR
    (void)equation;
    nr_warn("nompfr", "NRoots requires the MPFR backend.");
    return NULL;
#else
    /* Target / working precision (driven by PrecisionGoal; Automatic = machine). */
    long target_bits;
    bool want_machine;
    if (opts.prec_goal < 0.0) {
        target_bits = 53; want_machine = true;
    } else {
        target_bits = numeric_digits_to_bits(opts.prec_goal);
        if (target_bits <= 53) { target_bits = 53; want_machine = true; }
        else want_machine = false;
    }
    long guard = target_bits / 2;
    if (guard < 48) guard = 48;
    mpfr_prec_t wp = (mpfr_prec_t)(target_bits + guard);

    bool was_constant = false, constant_zero = false;
    ZUPoly* zp = NULL;
    NrPoly* p = nr_build_poly(equation, var, wp, &was_constant, &constant_zero, &zp);
    if (!p) {
        if (zp) zupoly_free(zp);
        if (was_constant)
            return expr_new_symbol(constant_zero ? SYM_True : SYM_False);
        return NULL;  /* diagnostic already emitted */
    }

    int d = p->deg;
    int max_iter = (opts.max_iter > 0) ? opts.max_iter : (100 + 20 * d);

    ncpx* roots = (ncpx*)malloc(sizeof(ncpx) * (size_t)d);
    for (int i = 0; i < d; i++) ncpx_init(&roots[i], wp);

    /* Exact integer polynomials of degree >= 2 take the squarefree path so that
     * multiple roots stay well-conditioned; everything else (and any Yun
     * failure) uses the direct numeric engine. */
    int rc;
    int poly_is_real;
    if (zp && zp->deg == d && d >= 2
        && nr_solve_exact(zp, opts.method, max_iter, wp, roots) == 0) {
        rc = 0; poly_is_real = 1;
    } else {
        rc = nr_solve(p, opts.method, max_iter, wp, roots);
        poly_is_real = p->is_real;
    }
    if (zp) zupoly_free(zp);
    nr_poly_free(p);

    if (rc != 0) {
        nr_warn("conv", "root-finding did not converge.");
        for (int i = 0; i < d; i++) ncpx_clear(&roots[i]);
        free(roots);
        return NULL;
    }

    int* realf = (int*)malloc(sizeof(int) * (size_t)d);
    nr_postprocess(roots, realf, d, poly_is_real, wp, target_bits);

    /* Assemble  var==r1 || ... || var==rd  (bare equation for d == 1). */
    Expr* out;
    if (d == 1) {
        Expr* val = nr_root_to_expr(&roots[0], realf[0], target_bits, want_machine);
        out = expr_new_function(expr_new_symbol(SYM_Equal),
                  (Expr*[]){ expr_copy(var), val }, 2);
    } else {
        Expr** eqs = (Expr**)malloc(sizeof(Expr*) * (size_t)d);
        for (int i = 0; i < d; i++) {
            Expr* val = nr_root_to_expr(&roots[i], realf[i], target_bits, want_machine);
            eqs[i] = expr_new_function(expr_new_symbol(SYM_Equal),
                         (Expr*[]){ expr_copy(var), val }, 2);
        }
        out = expr_new_function(expr_new_symbol(SYM_Or), eqs, (size_t)d);
        free(eqs);
    }

    for (int i = 0; i < d; i++) ncpx_clear(&roots[i]);
    free(roots);
    free(realf);
    return out;
#endif /* USE_MPFR */
}

#ifdef USE_MPFR
/* ================================================================== *
 *  CompanionMatrix engine
 *
 *  Real coefficients: the n*n Frobenius companion matrix is real, so its
 *  eigenvalues come straight from the existing real MPFR QR kernel.
 *
 *  Complex coefficients: the companion matrix is complex; we use the standard
 *  C^{n*n} -> R^{2n*2n} isomorphism  M = [[Re, -Im], [Im, Re]],  whose spectrum
 *  is {lambda_j} U {conj lambda_j}.  We run the real QR on the 2n matrix,
 *  classify which eigenvalues are genuine roots of p by their residual |p(z)|,
 *  Newton-polish those, deduplicate, and read each root's multiplicity directly
 *  from p by repeated synthetic division (Taylor coefficients at the root).
 *  This reuses the proven real kernel and is exact except for the measure-zero
 *  case of a complex-coefficient polynomial having a conjugate pair of roots
 *  with *unequal* multiplicities, which the embedding cannot distinguish.
 * ================================================================== */

/* Multiplicity of root r in p via successive synthetic division: the k-th
 * remainder is the Taylor coefficient a_k = p^(k)(r)/k!; m is the index of the
 * first a_k that is non-negligible. */
static int nr_multiplicity(const NrPoly* p, const ncpx* r, mpfr_prec_t wp) {
    int d = p->deg;
    ncpx* w = (ncpx*)malloc(sizeof(ncpx) * (size_t)(d + 1));
    ncpx* q = (ncpx*)malloc(sizeof(ncpx) * (size_t)(d + 1));
    for (int i = 0; i <= d; i++) { ncpx_init(&w[i], wp); ncpx_init(&q[i], wp); ncpx_set(&w[i], &p->c[i]); }

    mpfr_t maxc, ak, tol, t;
    mpfr_init2(maxc, wp); mpfr_init2(ak, wp); mpfr_init2(tol, wp); mpfr_init2(t, wp);
    mpfr_set_zero(maxc, 1);
    for (int i = 0; i <= d; i++) { ncpx_abs(t, &p->c[i]); if (mpfr_cmp(t, maxc) > 0) mpfr_set(maxc, t, MPFR_RNDN); }
    /* tol = maxc * 2^-(wp/3): comfortably above the O(eps^(1/m)) residual a
     * genuine multiple root leaves, below an O(maxc) non-root coefficient. */
    mpfr_set(tol, maxc, MPFR_RNDN);
    mpfr_div_2si(tol, tol, (long)wp / 3, MPFR_RNDN);

    int m = 0, curdeg = d;
    for (int k = 0; k <= d; k++) {
        /* remainder of w (degree curdeg) divided by (x - r); quotient -> q. */
        ncpx rem;
        ncpx_init(&rem, wp);
        if (curdeg == 0) {
            ncpx_set(&rem, &w[0]);
        } else {
            ncpx_set(&q[curdeg - 1], &w[curdeg]);
            for (int i = curdeg - 1; i >= 1; i--) {
                ncpx_mul(&q[i - 1], &q[i], r, wp);
                ncpx_add(&q[i - 1], &q[i - 1], &w[i]);
            }
            ncpx_mul(&rem, &q[0], r, wp);
            ncpx_add(&rem, &rem, &w[0]);
        }
        ncpx_abs(ak, &rem);
        ncpx_clear(&rem);
        if (mpfr_cmp(ak, tol) > 0) break;   /* a_k != 0  ->  multiplicity is k */
        m = k + 1;
        if (curdeg == 0) break;
        for (int i = 0; i < curdeg; i++) ncpx_set(&w[i], &q[i]);
        curdeg--;
    }

    for (int i = 0; i <= d; i++) { ncpx_clear(&w[i]); ncpx_clear(&q[i]); }
    free(w); free(q);
    mpfr_clear(maxc); mpfr_clear(ak); mpfr_clear(tol); mpfr_clear(t);
    return m;
}

int nr_companion(const NrPoly* p, ncpx* roots, mpfr_prec_t wp) {
    int n = p->deg;
    if (n == 1) {
        ncpx_div(&roots[0], &p->c[0], &p->c[1], wp);
        mpfr_neg(roots[0].re, roots[0].re, MPFR_RNDN);
        mpfr_neg(roots[0].im, roots[0].im, MPFR_RNDN);
        return 0;
    }

    /* Normalized last column  -c_i/c_n  (complex). */
    ncpx* lastcol = (ncpx*)malloc(sizeof(ncpx) * (size_t)n);
    for (int i = 0; i < n; i++) {
        ncpx_init(&lastcol[i], wp);
        ncpx_div(&lastcol[i], &p->c[i], &p->c[n], wp);
        mpfr_neg(lastcol[i].re, lastcol[i].re, MPFR_RNDN);
        mpfr_neg(lastcol[i].im, lastcol[i].im, MPFR_RNDN);
    }

    int rc = 0;
    if (p->is_real) {
        /* ---- Real n*n Frobenius companion ---- */
        mpfr_t* M = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)(n * n));
        for (int i = 0; i < n * n; i++) { mpfr_init2(M[i], wp); mpfr_set_zero(M[i], 1); }
        for (int i = 0; i + 1 < n; i++) mpfr_set_ui(M[(i + 1) * n + i], 1, MPFR_RNDN);
        for (int i = 0; i < n; i++) mpfr_set(M[i * n + (n - 1)], lastcol[i].re, MPFR_RNDN);

        mpfr_t* er = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)n);
        mpfr_t* ei = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)n);
        for (int i = 0; i < n; i++) { mpfr_init2(er[i], wp); mpfr_init2(ei[i], wp); }

        rc = eigen_all_eigenvalues_real_mpfr(M, (size_t)n, wp, er, ei);
        if (rc == 0) {
            for (int i = 0; i < n; i++) {
                mpfr_set(roots[i].re, er[i], MPFR_RNDN);
                mpfr_set(roots[i].im, ei[i], MPFR_RNDN);
                nr_newton_polish(p, &roots[i], wp, 12);
            }
        }
        for (int i = 0; i < n * n; i++) mpfr_clear(M[i]);
        for (int i = 0; i < n; i++) { mpfr_clear(er[i]); mpfr_clear(ei[i]); }
        free(M); free(er); free(ei);
    } else {
        /* ---- Complex companion via 2n real embedding ---- */
        int N = 2 * n;
        mpfr_t* M = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)(N * N));
        for (int i = 0; i < N * N; i++) { mpfr_init2(M[i], wp); mpfr_set_zero(M[i], 1); }
        /* Companion real/imag parts cre, cim (n*n). Sub-diagonal ones are real;
         * the last column carries lastcol's re/im. */
        for (int rr = 0; rr < n; rr++) {
            for (int cc = 0; cc < n; cc++) {
                /* sub-diagonal: companion entry (rr, cc) == 1 when rr == cc + 1. */
                int sub = (rr == cc + 1) ? 1 : 0;
                int last = (cc == n - 1) ? 1 : 0;
                /* Re part. */
                if (sub) mpfr_set_ui(M[rr * N + cc], 1, MPFR_RNDN);
                if (last) mpfr_add(M[rr * N + cc], M[rr * N + cc], lastcol[rr].re, MPFR_RNDN);
                /* top-right block = -cim. */
                if (last) mpfr_neg(M[rr * N + (cc + n)], lastcol[rr].im, MPFR_RNDN);
                /* bottom-left block = cim. */
                if (last) mpfr_set(M[(rr + n) * N + cc], lastcol[rr].im, MPFR_RNDN);
                /* bottom-right block = cre. */
                if (sub) mpfr_set_ui(M[(rr + n) * N + (cc + n)], 1, MPFR_RNDN);
                if (last) mpfr_add(M[(rr + n) * N + (cc + n)], M[(rr + n) * N + (cc + n)], lastcol[rr].re, MPFR_RNDN);
            }
        }

        mpfr_t* er = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)N);
        mpfr_t* ei = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)N);
        for (int i = 0; i < N; i++) { mpfr_init2(er[i], wp); mpfr_init2(ei[i], wp); }
        rc = eigen_all_eigenvalues_real_mpfr(M, (size_t)N, wp, er, ei);

        if (rc == 0) {
            /* Classify genuine roots by raw residual |p(eig)|; collect them. */
            ncpx cand; ncpx_init(&cand, wp);
            ncpx val; ncpx_init(&val, wp);
            mpfr_t resid, maxc, rtol, t;
            mpfr_init2(resid, wp); mpfr_init2(maxc, wp); mpfr_init2(rtol, wp); mpfr_init2(t, wp);
            mpfr_set_zero(maxc, 1);
            for (int i = 0; i <= n; i++) { ncpx_abs(t, &p->c[i]); if (mpfr_cmp(t, maxc) > 0) mpfr_set(maxc, t, MPFR_RNDN); }
            mpfr_set(rtol, maxc, MPFR_RNDN);
            mpfr_div_2si(rtol, rtol, (long)wp / 4, MPFR_RNDN);   /* root residual gate */

            ncpx* kept = (ncpx*)malloc(sizeof(ncpx) * (size_t)N);
            int nk = 0;
            for (int i = 0; i < N; i++) {
                mpfr_set(cand.re, er[i], MPFR_RNDN);
                mpfr_set(cand.im, ei[i], MPFR_RNDN);
                nr_poly_eval(p, &cand, &val, NULL, wp);
                ncpx_abs(resid, &val);
                if (mpfr_cmp(resid, rtol) <= 0) {
                    ncpx_init(&kept[nk], wp);
                    nr_newton_polish(p, &cand, wp, 12);
                    ncpx_set(&kept[nk], &cand);
                    nk++;
                }
            }

            /* Deduplicate kept roots; multiplicity read from p. */
            mpfr_t ctol; mpfr_init2(ctol, wp);
            mpfr_set(ctol, maxc, MPFR_RNDN);
            mpfr_div_2si(ctol, ctol, (long)wp / 3, MPFR_RNDN);
            int* used = (int*)calloc((size_t)nk, sizeof(int));
            int out = 0;
            for (int i = 0; i < nk && rc == 0; i++) {
                if (used[i]) continue;
                used[i] = 1;
                for (int j = i + 1; j < nk; j++) {
                    if (used[j]) continue;
                    if (nr_close(&kept[j], &kept[i], ctol, wp)) used[j] = 1;
                }
                int m = nr_multiplicity(p, &kept[i], wp);
                if (m < 1) m = 1;
                for (int t2 = 0; t2 < m; t2++) {
                    if (out >= n) { rc = -1; break; }
                    ncpx_set(&roots[out], &kept[i]);
                    out++;
                }
            }
            if (rc == 0 && out != n) rc = -1;   /* multiplicity bookkeeping failed */

            free(used);
            mpfr_clear(ctol);
            for (int i = 0; i < nk; i++) ncpx_clear(&kept[i]);
            free(kept);
            ncpx_clear(&cand); ncpx_clear(&val);
            mpfr_clear(resid); mpfr_clear(maxc); mpfr_clear(rtol); mpfr_clear(t);
        }

        for (int i = 0; i < N * N; i++) mpfr_clear(M[i]);
        for (int i = 0; i < N; i++) { mpfr_clear(er[i]); mpfr_clear(ei[i]); }
        free(M); free(er); free(ei);
    }

    for (int i = 0; i < n; i++) ncpx_clear(&lastcol[i]);
    free(lastcol);
    return rc;
}

/* JenkinsTraub engine lives in nroots_jt.c. */
#endif

/* ================================================================== *
 *  Registration
 * ================================================================== */
void nroots_init(void) {
    symtab_add_builtin("NRoots", builtin_nroots);
    /* Protected only.  The equation is evaluated normally (constants fold, the
     * variable stays symbolic); coefficients are numericalised after symbolic
     * extraction, so no argument holding is required.  Not Listable. */
    symtab_get_def("NRoots")->attributes |= ATTR_PROTECTED;
}
