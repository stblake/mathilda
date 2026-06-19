/*
 * nsolve_system.c — numerical solver for zero-dimensional polynomial systems.
 *
 * Implements nsolve_polynomial_system() (see nsolve_system.h) for square /
 * zero-dimensional systems  f_1 == ... == f_m == 0  in variables x_1..x_n.
 *
 * Two methods:
 *
 *   NSYS_ENDOMORPHISM (default) — the eigenvalue / multiplication-matrix
 *      (Möller–Stetter) method.  A Gröbner basis over Q gives the quotient
 *      ring  A = Q[x]/I.  When I is zero-dimensional, A is a finite-
 *      dimensional Q-vector space with the standard-monomial basis B (the
 *      monomials not divisible by any leading term).  For a generic linear
 *      form  l = sum c_i x_i  the multiplication map  M_l : A -> A,
 *      g |-> l*g mod I  has eigenvalues l(p) over the solutions p, and the
 *      shared eigenvectors v_p satisfy  M_{x_i} v_p = x_i(p) v_p, so each
 *      coordinate x_i(p) is read off as (M_{x_i} v_p)[j] / v_p[j].  The
 *      eigenproblem is solved at MPFR precision via the real-matrix
 *      eigenvector backend (eigen_all_eigenvectors_real_mpfr).
 *
 *   NSYS_HOMOTOPY — currently routed to the same eigenvalue engine (a true
 *      homotopy continuation tracker is a documented non-goal for now).
 *
 * The caller (nsolve.c) selects the elimination / triangular fallback
 * separately (Method -> "Symbolic"); see nsolve_system_eliminate().
 *
 * Memory: returns a freshly owned List of rule-lists  {{x->r, ...}, ...},
 * or NULL when the system is outside the supported envelope (not a numeric
 * polynomial system, positive-dimensional, or too large), so the caller can
 * fall back or leave NSolve unevaluated.
 */

#include "nsolve_system.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <gmp.h>
#include <mpfr.h>

#include "expr.h"
#include "eval.h"
#include "sym_names.h"
#include "numeric.h"
#include "numeric_complex.h"
#include "arithmetic.h"
#include "poly/groebner.h"
#include "poly/poly.h"        /* is_polynomial, get_degree_poly */
#include "linalg/eigen.h"
#include "nroots.h"           /* builtin_nroots */

/* Safety caps: above these the system path bails (returns NULL) rather than
 * risk a combinatorial / memory explosion. */
#define NSYS_MAX_DIM      256   /* quotient-ring dimension (matrix size)     */
#define NSYS_MAX_BOX   200000   /* product of per-variable degree bounds     */
#define NSYS_MAX_TDEG     60    /* per-generator total degree (anti-hang gate) */

/* Total degree (max over terms of the exponent sum) of a GBPoly. */
static int gbpoly_total_degree(const GBPoly* p) {
    int best = 0;
    for (size_t t = 0; t < p->n_terms; t++) {
        const int* e = p->exps + t * (size_t)p->n_vars;
        int s = 0;
        for (int k = 0; k < p->n_vars; k++) s += e[k];
        if (s > best) best = s;
    }
    return best;
}

/* True iff every generator is within the total-degree gate (else the
 * Gröbner computation could blow up; the caller bails to NULL). */
static bool gb_set_degree_ok(GBPoly* const* F, int nF) {
    for (int i = 0; i < nF; i++)
        if (gbpoly_total_degree(F[i]) > NSYS_MAX_TDEG) return false;
    return true;
}

/* ================================================================== *
 *  Monomial helpers (exponent vectors are int[nvar]).
 * ================================================================== */

static bool mono_divides(const int* a, const int* b, int n) {
    /* a | b ? */
    for (int k = 0; k < n; k++) if (a[k] > b[k]) return false;
    return true;
}

/* Index of exponent vector `e` within basis `B` (dB rows of n ints), or -1. */
static int mono_index(const int* B, int dB, int n, const int* e) {
    for (int j = 0; j < dB; j++) {
        const int* b = B + (size_t)j * n;
        int k = 0;
        for (; k < n; k++) if (b[k] != e[k]) break;
        if (k == n) return j;
    }
    return -1;
}

/* ================================================================== *
 *  Verification: max |f_i(point)| over the residual polynomials.
 *  `vals` holds the candidate coordinates (ncpx, length nvar).  Uses the
 *  evaluation engine (ReplaceAll + numericalize), returning the residual
 *  magnitude as a double (verification only — precision is not critical).
 * ================================================================== */

static double value_magnitude(const Expr* e) {
    if (!e) return INFINITY;
    if (e->type == EXPR_INTEGER) return fabs((double)e->data.integer);
    if (e->type == EXPR_REAL)    return fabs(e->data.real);
    if (e->type == EXPR_BIGINT)  return fabs(mpz_get_d(e->data.bigint));
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    return fabs(mpfr_get_d(e->data.mpfr, MPFR_RNDN));
#endif
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Complex
        && e->data.function.arg_count == 2) {
        double re = value_magnitude(e->data.function.args[0]);
        double im = value_magnitude(e->data.function.args[1]);
        return hypot(re, im);
    }
    return INFINITY; /* unevaluated / symbolic -> treat as non-zero */
}

/* Build a Complex[Real,Real] (or Real) Expr from an ncpx at machine precision,
 * for substitution into residuals. */
static Expr* ncpx_to_expr_machine(const ncpx* z) {
    double re = mpfr_get_d(z->re, MPFR_RNDN);
    double im = mpfr_get_d(z->im, MPFR_RNDN);
    if (im == 0.0) return expr_new_real(re);
    return expr_new_function(expr_new_symbol(SYM_Complex),
               (Expr*[]){ expr_new_real(re), expr_new_real(im) }, 2);
}

/* Max |f_i| after applying the rule-list `rulelist` (borrowed) to each poly. */
static double residual_with_rules(Expr** polys, int npoly, Expr* rulelist) {
    double worst = 0.0;
    NumericSpec spec = numeric_machine_spec();
    for (int p = 0; p < npoly; p++) {
        Expr* ra = expr_new_function(expr_new_symbol(SYM_ReplaceAll),
                       (Expr*[]){ expr_copy(polys[p]), expr_copy(rulelist) }, 2);
        Expr* sub = eval_and_free(ra);
        Expr* num = eval_and_free(numericalize(sub, spec));
        expr_free(sub);
        double m = value_magnitude(num);
        expr_free(num);
        if (m > worst) worst = m;
    }
    return worst;
}

static double system_residual(Expr** polys, int npoly, Expr** vars, int nvar,
                              const ncpx* vals) {
    Expr** rules = (Expr**)malloc(sizeof(Expr*) * (size_t)nvar);
    for (int i = 0; i < nvar; i++) {
        rules[i] = expr_new_function(expr_new_symbol(SYM_Rule),
                       (Expr*[]){ expr_copy(vars[i]), ncpx_to_expr_machine(&vals[i]) }, 2);
    }
    Expr* rulelist = expr_new_function(expr_new_symbol(SYM_List), rules, (size_t)nvar);
    free(rules);
    double worst = residual_with_rules(polys, npoly, rulelist);
    expr_free(rulelist);
    return worst;
}

/* ================================================================== *
 *  Eigenvalue / multiplication-matrix method.
 * ================================================================== */

/* Build the d*d rational multiplication-by-x_i matrix on basis B, stored
 * column-major-by-meaning as M[row*d + col]: column j (basis monomial B_j)
 * maps to the normal form of x_i * B_j expanded in B.  Returns false on any
 * structural surprise (a reduced term outside B). */
static bool build_mult_matrix(int var, const int* B, int d, int nvar,
                              GBPoly* const* basis, size_t nbasis,
                              GBOrder order, mpq_t* M) {
    bool ok = true;
    int* e = (int*)malloc(sizeof(int) * (size_t)nvar);
    for (int j = 0; j < d && ok; j++) {
        /* monomial x_var * B_j */
        const int* bj = B + (size_t)j * nvar;
        for (int k = 0; k < nvar; k++) e[k] = bj[k];
        e[var] += 1;

        GBPoly* mono = gb_poly_new(nvar, order, 0);
        gb_poly_push_term_si(mono, e, 1, 1);
        gb_poly_normalize(mono);

        GBPoly* rem = gb_reduce(mono, basis, nbasis);
        gb_poly_free(mono);

        for (size_t t = 0; t < rem->n_terms; t++) {
            const int* re = rem->exps + t * (size_t)nvar;
            int row = mono_index(B, d, nvar, re);
            if (row < 0) { ok = false; break; }
            mpq_set(M[(size_t)row * d + j], rem->coefs[t]);
        }
        gb_poly_free(rem);
    }
    free(e);
    return ok;
}

/* Deterministic small pseudo-random coefficient stream (seeded). */
static long nsys_rand(unsigned long* state) {
    *state = *state * 6364136223846793005UL + 1442695040888963407UL;
    /* small nonzero coefficient in [1, 9] with alternating sign */
    long v = (long)((*state >> 33) % 9UL) + 1;
    return ((*state >> 17) & 1UL) ? v : -v;
}

Expr* nsolve_polynomial_system(Expr** polys, int npoly, Expr** vars, int nvar,
                               NSysMethod method, bool reals_only,
                               bool want_machine, long bits, long max_roots,
                               int verify, unsigned long seed) {
    (void)method; /* HOMOTOPY and ENDOMORPHISM both use the eigenvalue engine */
    if (npoly < 1 || nvar < 1) return NULL;

    GBOrder order = GB_ORDER_GREVLEX;

    /* 1. Convert each polynomial to a GBPoly over Q. */
    GBPoly** F = (GBPoly**)malloc(sizeof(GBPoly*) * (size_t)npoly);
    int nF = 0;
    bool conv_ok = true;
    for (int i = 0; i < npoly; i++) {
        GBPoly* g = gb_from_expr(polys[i], vars, nvar, order, 0, NULL);
        if (!g) { conv_ok = false; break; }
        if (!gb_poly_is_zero(g)) F[nF++] = g;
        else gb_poly_free(g);
    }
    if (!conv_ok || nF == 0 || !gb_set_degree_ok(F, nF)) {
        for (int i = 0; i < nF; i++) gb_poly_free(F[i]);
        free(F);
        return NULL;
    }

    /* 2. Gröbner basis. */
    size_t nG = 0;
    GBPoly** G = gb_buchberger(F, (size_t)nF, &nG);
    for (int i = 0; i < nF; i++) gb_poly_free(F[i]);
    free(F);
    if (!G || nG == 0) { if (G) gb_basis_free(G, nG); return NULL; }

    /* Inconsistent system: GB contains a nonzero constant -> no solutions. */
    for (size_t i = 0; i < nG; i++) {
        if (gb_poly_is_constant(G[i]) && !gb_poly_is_zero(G[i])) {
            gb_basis_free(G, nG);
            return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
        }
    }

    /* 3. Zero-dimensionality + per-variable degree bound: for each variable
     *    some leading monomial must be a pure power of that variable. */
    int* bound = (int*)malloc(sizeof(int) * (size_t)nvar);
    for (int k = 0; k < nvar; k++) bound[k] = 0;
    for (size_t i = 0; i < nG; i++) {
        const int* lm = gb_poly_lm(G[i]);
        if (!lm) continue;
        int nz = -1, cnt = 0;
        for (int k = 0; k < nvar; k++) if (lm[k] > 0) { nz = k; cnt++; }
        if (cnt == 1) {
            if (bound[nz] == 0 || lm[nz] < bound[nz]) bound[nz] = lm[nz];
        }
    }
    bool zero_dim = true;
    double box = 1.0;
    for (int k = 0; k < nvar; k++) {
        if (bound[k] == 0) { zero_dim = false; break; }
        box *= (double)bound[k];
    }
    if (!zero_dim || box > (double)NSYS_MAX_BOX) {
        free(bound); gb_basis_free(G, nG);
        return NULL;
    }

    /* 4. Enumerate the standard-monomial basis B (box minus LM multiples). */
    int* B = NULL; int d = 0; int Bcap = 0;
    {
        int* e = (int*)calloc((size_t)nvar, sizeof(int));
        bool overflow = false;
        for (;;) {
            /* test e: keep iff not divisible by any leading monomial. */
            bool standard = true;
            for (size_t i = 0; i < nG; i++) {
                const int* lm = gb_poly_lm(G[i]);
                if (lm && mono_divides(lm, e, nvar)) { standard = false; break; }
            }
            if (standard) {
                if (d == Bcap) {
                    Bcap = Bcap ? Bcap * 2 : 16;
                    B = (int*)realloc(B, sizeof(int) * (size_t)Bcap * (size_t)nvar);
                }
                memcpy(B + (size_t)d * nvar, e, sizeof(int) * (size_t)nvar);
                d++;
                if (d > NSYS_MAX_DIM) { overflow = true; break; }
            }
            /* increment e as a mixed-radix counter over the box. */
            int k = 0;
            for (; k < nvar; k++) {
                e[k]++;
                if (e[k] < bound[k]) break;
                e[k] = 0;
            }
            if (k == nvar) break; /* wrapped around -> done */
        }
        free(e);
        if (overflow || d == 0) { free(B); free(bound); gb_basis_free(G, nG); return NULL; }
    }
    free(bound);

    /* 5. Multiplication matrices M_{x_i} (rational), and a generic M_l. */
    mpq_t* Mx = (mpq_t*)malloc(sizeof(mpq_t) * (size_t)nvar * (size_t)d * (size_t)d);
    for (size_t t = 0; t < (size_t)nvar * d * d; t++) mpq_init(Mx[t]);
    bool mok = true;
    for (int i = 0; i < nvar && mok; i++) {
        mok = build_mult_matrix(i, B, d, nvar, G, nG, order,
                                Mx + (size_t)i * d * d);
    }
    gb_basis_free(G, nG);
    if (!mok) {
        for (size_t t = 0; t < (size_t)nvar * d * d; t++) mpq_clear(Mx[t]);
        free(Mx); free(B);
        return NULL;
    }

    /* Working precision: target + guard. */
    long target_bits = want_machine ? 53 : bits;
    long guard = target_bits / 2; if (guard < 32) guard = 32;
    mpfr_prec_t wp = (mpfr_prec_t)(target_bits + guard);

    /* Generic linear form l = sum c_i x_i ; M_l = sum c_i M_{x_i} (real). */
    unsigned long st = seed ? seed : 1234UL;
    long* cform = (long*)malloc(sizeof(long) * (size_t)nvar);
    for (int i = 0; i < nvar; i++) cform[i] = nsys_rand(&st);

    mpfr_t* Ml = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)d * (size_t)d);
    for (size_t t = 0; t < (size_t)d * d; t++) mpfr_init2(Ml[t], wp);
    {
        mpq_t acc, tmp; mpq_init(acc); mpq_init(tmp);
        for (size_t cell = 0; cell < (size_t)d * d; cell++) {
            mpq_set_ui(acc, 0, 1);
            for (int i = 0; i < nvar; i++) {
                mpq_set_si(tmp, cform[i], 1);
                mpq_mul(tmp, tmp, Mx[(size_t)i * d * d + cell]);
                mpq_add(acc, acc, tmp);
            }
            mpfr_set_q(Ml[cell], acc, MPFR_RNDN);
        }
        mpq_clear(acc); mpq_clear(tmp);
    }
    free(cform);

    /* 6. Eigenvalues + eigenvectors of M_l at MPFR precision. */
    mpfr_t* ev_re = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)d);
    mpfr_t* ev_im = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)d);
    mpfr_t* V_re  = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)d * (size_t)d);
    mpfr_t* V_im  = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)d * (size_t)d);
    for (int i = 0; i < d; i++) { mpfr_init2(ev_re[i], wp); mpfr_init2(ev_im[i], wp); }
    for (size_t t = 0; t < (size_t)d * d; t++) { mpfr_init2(V_re[t], wp); mpfr_init2(V_im[t], wp); }

    int est = eigen_all_eigenvectors_real_mpfr(Ml, (size_t)d, wp,
                                               ev_re, ev_im, V_re, V_im);
    for (size_t t = 0; t < (size_t)d * d; t++) mpfr_clear(Ml[t]);
    free(Ml);
    if (est != 0) {
        for (int i = 0; i < d; i++) { mpfr_clear(ev_re[i]); mpfr_clear(ev_im[i]); }
        for (size_t t = 0; t < (size_t)d * d; t++) { mpfr_clear(V_re[t]); mpfr_clear(V_im[t]); }
        free(ev_re); free(ev_im); free(V_re); free(V_im);
        for (size_t t = 0; t < (size_t)nvar * d * d; t++) mpq_clear(Mx[t]);
        free(Mx); free(B);
        return NULL;
    }

    /* 7. Coordinate recovery: for eigenvector row k (v_k), pick its largest
     *    component j*, then x_i(p_k) = (M_{x_i} v_k)[j*] / v_k[j*]. */
    double verify_tol = pow(10.0, -0.5 * numeric_bits_to_digits(target_bits));
    if (verify_tol < 1e-6) verify_tol = 1e-6;

    /* Real multiplication matrices as MPFR (reused for every eigenvector). */
    mpfr_t* MxF = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)nvar * (size_t)d * (size_t)d);
    for (size_t t = 0; t < (size_t)nvar * d * d; t++) {
        mpfr_init2(MxF[t], wp);
        mpfr_set_q(MxF[t], Mx[t], MPFR_RNDN);
    }
    for (size_t t = 0; t < (size_t)nvar * d * d; t++) mpq_clear(Mx[t]);
    free(Mx);

    Expr** sols = (Expr**)malloc(sizeof(Expr*) * (size_t)d);
    int nsol = 0;

    ncpx *vk = (ncpx*)malloc(sizeof(ncpx) * (size_t)d);
    for (int t = 0; t < d; t++) ncpx_init(&vk[t], wp);
    ncpx *coord = (ncpx*)malloc(sizeof(ncpx) * (size_t)nvar);
    for (int i = 0; i < nvar; i++) ncpx_init(&coord[i], wp);
    ncpx acc, prod, pivot, num;
    ncpx_init(&acc, wp); ncpx_init(&prod, wp); ncpx_init(&pivot, wp); ncpx_init(&num, wp);

    for (int k = 0; k < d && (max_roots < 0 || nsol < max_roots); k++) {
        /* load v_k and find pivot component. */
        int jstar = 0; mpfr_t best, mag; mpfr_init2(best, wp); mpfr_init2(mag, wp);
        mpfr_set_si(best, -1, MPFR_RNDN);
        for (int j = 0; j < d; j++) {
            mpfr_set(vk[j].re, V_re[(size_t)k * d + j], MPFR_RNDN);
            mpfr_set(vk[j].im, V_im[(size_t)k * d + j], MPFR_RNDN);
            mpfr_hypot(mag, vk[j].re, vk[j].im, MPFR_RNDN);
            if (mpfr_cmp(mag, best) > 0) { mpfr_set(best, mag, MPFR_RNDN); jstar = j; }
        }
        mpfr_clear(best); mpfr_clear(mag);
        mpfr_set(pivot.re, vk[jstar].re, MPFR_RNDN);
        mpfr_set(pivot.im, vk[jstar].im, MPFR_RNDN);

        for (int i = 0; i < nvar; i++) {
            /* num = (M_{x_i} v_k)[jstar] = sum_m MxF[i][jstar][m] * v_k[m] */
            mpfr_set_zero(num.re, 1); mpfr_set_zero(num.im, 1);
            const mpfr_t* Mi = MxF + (size_t)i * d * d;
            for (int m = 0; m < d; m++) {
                /* prod = Mi[jstar*d+m] * vk[m]  (real * complex) */
                mpfr_mul(prod.re, Mi[(size_t)jstar * d + m], vk[m].re, MPFR_RNDN);
                mpfr_mul(prod.im, Mi[(size_t)jstar * d + m], vk[m].im, MPFR_RNDN);
                mpfr_add(num.re, num.re, prod.re, MPFR_RNDN);
                mpfr_add(num.im, num.im, prod.im, MPFR_RNDN);
            }
            ncpx_div(&coord[i], &num, &pivot, wp);
        }

        /* Verify residual; drop if too large. */
        double resid = verify ? system_residual(polys, npoly, vars, nvar, coord) : 0.0;
        if (verify && (resid > verify_tol || !(resid == resid))) continue;

        /* Reals filter. */
        if (reals_only) {
            bool real_ok = true;
            mpfr_t a; mpfr_init2(a, wp);
            for (int i = 0; i < nvar; i++) {
                mpfr_abs(a, coord[i].im, MPFR_RNDN);
                if (mpfr_get_d(a, MPFR_RNDN) > verify_tol) { real_ok = false; break; }
            }
            mpfr_clear(a);
            if (!real_ok) continue;
        }

        /* Build {x_i -> value} rule list at the requested precision. */
        Expr** rules = (Expr**)malloc(sizeof(Expr*) * (size_t)nvar);
        for (int i = 0; i < nvar; i++) {
            Expr* val;
            if (want_machine) {
                double re = mpfr_get_d(coord[i].re, MPFR_RNDN);
                double im = mpfr_get_d(coord[i].im, MPFR_RNDN);
                if (reals_only || fabs(im) <= verify_tol * (1.0 + fabs(re)))
                    val = expr_new_real(re);
                else
                    val = expr_new_function(expr_new_symbol(SYM_Complex),
                              (Expr*[]){ expr_new_real(re), expr_new_real(im) }, 2);
            } else {
                mpfr_t rr, ii; mpfr_init2(rr, target_bits); mpfr_init2(ii, target_bits);
                mpfr_set(rr, coord[i].re, MPFR_RNDN);
                if (reals_only) mpfr_set_zero(ii, 1);
                else            mpfr_set(ii, coord[i].im, MPFR_RNDN);
                val = numeric_mpfr_make_complex(rr, ii);
                mpfr_clear(rr); mpfr_clear(ii);
            }
            rules[i] = expr_new_function(expr_new_symbol(SYM_Rule),
                           (Expr*[]){ expr_copy(vars[i]), val }, 2);
        }
        sols[nsol++] = expr_new_function(expr_new_symbol(SYM_List), rules, (size_t)nvar);
        free(rules);
    }

    /* cleanup */
    for (int t = 0; t < d; t++) ncpx_clear(&vk[t]);
    free(vk);
    for (int i = 0; i < nvar; i++) ncpx_clear(&coord[i]);
    free(coord);
    ncpx_clear(&acc); ncpx_clear(&prod); ncpx_clear(&pivot); ncpx_clear(&num);
    for (size_t t = 0; t < (size_t)nvar * d * d; t++) mpfr_clear(MxF[t]);
    free(MxF);
    for (int i = 0; i < d; i++) { mpfr_clear(ev_re[i]); mpfr_clear(ev_im[i]); }
    for (size_t t = 0; t < (size_t)d * d; t++) { mpfr_clear(V_re[t]); mpfr_clear(V_im[t]); }
    free(ev_re); free(ev_im); free(V_re); free(V_im);
    free(B);

    /* No accepted solution though the quotient ring is non-trivial: signal
     * failure (NULL) so the caller can fall back to the elimination method
     * (rather than reporting a false empty set). */
    if (nsol == 0 && d > 0) { free(sols); return NULL; }

    Expr* out = expr_new_function(expr_new_symbol(SYM_List), sols, (size_t)nsol);
    free(sols);
    return out;
}

/* ================================================================== *
 *  Elimination / triangular method (Method -> "Symbolic"; also the
 *  eigenvalue method's fallback).  A lexicographic Gröbner basis is
 *  triangular for a shape-position ideal: solve the univariate generator
 *  in the last variable with NRoots, back-substitute, and recurse,
 *  verifying every completed tuple against the original residuals.
 * ================================================================== */

/* Numeric roots of a univariate residual `uni` in `var` via NRoots.  Returns
 * a freshly malloc'd array of value Exprs (owned); count in *out_n.  NULL/0
 * when NRoots cannot reduce `uni` (e.g. it is constant). */
static Expr** elim_nroots(Expr* uni, Expr* var, bool want_machine,
                          long target_bits, int* out_n) {
    *out_n = 0;
    Expr* eqn = expr_new_function(expr_new_symbol(SYM_Equal),
                    (Expr*[]){ expr_copy(uni), expr_new_integer(0) }, 2);
    Expr* call;
    if (!want_machine) {
        Expr* pg = expr_new_function(expr_new_symbol(SYM_Rule),
                       (Expr*[]){ expr_new_symbol(SYM_PrecisionGoal),
                                  expr_new_real(numeric_bits_to_digits(target_bits)) }, 2);
        call = expr_new_function(expr_new_symbol(SYM_NRoots),
                   (Expr*[]){ eqn, expr_copy(var), pg }, 3);
    } else {
        call = expr_new_function(expr_new_symbol(SYM_NRoots),
                   (Expr*[]){ eqn, expr_copy(var) }, 2);
    }
    Expr* r = builtin_nroots(call);
    expr_free(call);
    if (!r) return NULL;
    if (r->type == EXPR_SYMBOL) { expr_free(r); return NULL; } /* True / False */

    Expr** vals = NULL; int n = 0, cap = 0;
    Expr** eqs; size_t neq; Expr* single[1];
    if (r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL
        && r->data.function.head->data.symbol == SYM_Or) {
        eqs = r->data.function.args; neq = r->data.function.arg_count;
    } else { single[0] = r; eqs = single; neq = 1; }
    for (size_t i = 0; i < neq; i++) {
        Expr* e = eqs[i];
        if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
            && e->data.function.head->data.symbol == SYM_Equal
            && e->data.function.arg_count == 2) {
            if (n == cap) { cap = cap ? cap * 2 : 4; vals = realloc(vals, sizeof(Expr*) * (size_t)cap); }
            vals[n++] = expr_copy(e->data.function.args[1]);
        }
    }
    expr_free(r);
    *out_n = n;
    return vals;
}

typedef struct {
    Expr** orig; int norig;
    Expr** gb;   int ngb;
    Expr** vars; int nvar;
    bool   reals_only, want_machine;
    long   target_bits, max_roots;
    double tol;
    Expr** sols; int nsol, solcap;
} ElimCtx;

static bool value_is_complex(const Expr* v) {
    return v->type == EXPR_FUNCTION
        && v->data.function.head->type == EXPR_SYMBOL
        && v->data.function.head->data.symbol == SYM_Complex;
}

/* `assigned` is a List of rules (x->val) for the last `depth` variables. */
static void elim_rec(ElimCtx* c, Expr* assigned, int depth) {
    if (c->max_roots >= 0 && c->nsol >= c->max_roots) return;

    if (depth == c->nvar) {
        if (residual_with_rules(c->orig, c->norig, assigned) > c->tol) return;
        Expr** rules = (Expr**)malloc(sizeof(Expr*) * (size_t)c->nvar);
        for (int i = 0; i < c->nvar; i++) {
            Expr* val = NULL;
            for (size_t a = 0; a < assigned->data.function.arg_count; a++) {
                Expr* rule = assigned->data.function.args[a];
                if (expr_eq(rule->data.function.args[0], c->vars[i])) {
                    val = rule->data.function.args[1]; break;
                }
            }
            rules[i] = expr_new_function(expr_new_symbol(SYM_Rule),
                          (Expr*[]){ expr_copy(c->vars[i]),
                                     val ? expr_copy(val) : expr_new_integer(0) }, 2);
        }
        if (c->nsol == c->solcap) {
            c->solcap = c->solcap ? c->solcap * 2 : 8;
            c->sols = realloc(c->sols, sizeof(Expr*) * (size_t)c->solcap);
        }
        c->sols[c->nsol++] = expr_new_function(expr_new_symbol(SYM_List),
                                               rules, (size_t)c->nvar);
        free(rules);
        return;
    }

    Expr* tvar = c->vars[c->nvar - 1 - depth];
    /* Lowest-degree generator that, after substituting the assigned (higher)
     * variables, is a non-constant polynomial in tvar ALONE — i.e. free of
     * the other still-unassigned variables vars[0 .. nvar-2-depth].  (Without
     * the freeness check, is_polynomial treats a remaining variable as a
     * symbolic coefficient, and NRoots then chokes on it.) */
    Expr* best = NULL; int bestdeg = 0;
    for (int g = 0; g < c->ngb; g++) {
        Expr* ra = expr_new_function(expr_new_symbol(SYM_ReplaceAll),
                       (Expr*[]){ expr_copy(c->gb[g]), expr_copy(assigned) }, 2);
        Expr* sub = eval_and_free(ra);
        bool usable = is_polynomial(sub, &tvar, 1)
                   && get_degree_poly(sub, tvar) >= 1;
        for (int j = 0; usable && j <= c->nvar - 2 - depth; j++)
            if (contains_any_symbol_from(sub, c->vars[j])) usable = false;
        if (usable) {
            int dg = get_degree_poly(sub, tvar);
            if (best == NULL || dg < bestdeg) {
                if (best) expr_free(best);
                best = sub; bestdeg = dg; continue;
            }
        }
        expr_free(sub);
    }
    if (!best) return;  /* no triangular generator on this branch */

    int nv = 0;
    Expr** vals = elim_nroots(best, tvar, c->want_machine, c->target_bits, &nv);
    expr_free(best);
    for (int i = 0; i < nv; i++) {
        if (!(c->reals_only && value_is_complex(vals[i]))) {
            size_t na = assigned->data.function.arg_count;
            Expr** na2 = (Expr**)malloc(sizeof(Expr*) * (na + 1));
            for (size_t a = 0; a < na; a++) na2[a] = expr_copy(assigned->data.function.args[a]);
            na2[na] = expr_new_function(expr_new_symbol(SYM_Rule),
                          (Expr*[]){ expr_copy(tvar), expr_copy(vals[i]) }, 2);
            Expr* na_list = expr_new_function(expr_new_symbol(SYM_List), na2, na + 1);
            free(na2);
            elim_rec(c, na_list, depth + 1);
            expr_free(na_list);
        }
    }
    for (int i = 0; i < nv; i++) expr_free(vals[i]);
    free(vals);
}

Expr* nsolve_system_eliminate(Expr** polys, int npoly, Expr** vars, int nvar,
                              bool reals_only, bool want_machine, long bits,
                              long max_roots, int verify) {
    (void)verify;
    if (npoly < 1 || nvar < 1) return NULL;

    GBPoly** F = (GBPoly**)malloc(sizeof(GBPoly*) * (size_t)npoly);
    int nF = 0; bool ok = true;
    for (int i = 0; i < npoly; i++) {
        GBPoly* g = gb_from_expr(polys[i], vars, nvar, GB_ORDER_LEX, 0, NULL);
        if (!g) { ok = false; break; }
        if (!gb_poly_is_zero(g)) F[nF++] = g; else gb_poly_free(g);
    }
    if (!ok || nF == 0 || !gb_set_degree_ok(F, nF)) {
        for (int i = 0; i < nF; i++) gb_poly_free(F[i]);
        free(F);
        return NULL;
    }

    size_t nG = 0;
    GBPoly** G = gb_groebner_walk(F, (size_t)nF, GB_ORDER_LEX, NULL, &nG);
    for (int i = 0; i < nF; i++) gb_poly_free(F[i]);
    free(F);
    if (!G || nG == 0) { if (G) gb_basis_free(G, nG); return NULL; }

    for (size_t i = 0; i < nG; i++)
        if (gb_poly_is_constant(G[i]) && !gb_poly_is_zero(G[i])) {
            gb_basis_free(G, nG);
            return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
        }

    /* Zero-dimensional?  Each variable needs a pure-power leading monomial.
     * A positive-dimensional (e.g. underdetermined) ideal is not handled
     * here -> NULL so the caller can fall back to symbolic Solve. */
    {
        int* haspp = (int*)calloc((size_t)nvar, sizeof(int));
        for (size_t i = 0; i < nG; i++) {
            const int* lm = gb_poly_lm(G[i]);
            if (!lm) continue;
            int nz = -1, cnt = 0;
            for (int k = 0; k < nvar; k++) if (lm[k] > 0) { nz = k; cnt++; }
            if (cnt == 1) haspp[nz] = 1;
        }
        bool zerodim = true;
        for (int k = 0; k < nvar; k++) if (!haspp[k]) { zerodim = false; break; }
        free(haspp);
        if (!zerodim) { gb_basis_free(G, nG); return NULL; }
    }

    Expr** gbE = (Expr**)malloc(sizeof(Expr*) * nG);
    for (size_t i = 0; i < nG; i++) gbE[i] = gb_to_expr(G[i], vars);
    gb_basis_free(G, nG);

    long target_bits = want_machine ? 53 : bits;
    double tol = pow(10.0, -0.5 * numeric_bits_to_digits(target_bits));
    if (tol < 1e-6) tol = 1e-6;

    ElimCtx c;
    c.orig = polys; c.norig = npoly; c.gb = gbE; c.ngb = (int)nG;
    c.vars = vars; c.nvar = nvar; c.reals_only = reals_only;
    c.want_machine = want_machine; c.target_bits = target_bits;
    c.max_roots = max_roots; c.tol = tol;
    c.sols = NULL; c.nsol = 0; c.solcap = 0;

    Expr* empty = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    elim_rec(&c, empty, 0);
    expr_free(empty);

    for (size_t i = 0; i < nG; i++) expr_free(gbE[i]);
    free(gbE);

    Expr* out = expr_new_function(expr_new_symbol(SYM_List), c.sols, (size_t)c.nsol);
    free(c.sols);
    return out;
}
