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
#ifdef USE_MPFR
#include <mpfr.h>
#endif

#include "expr.h"
#include "eval.h"
#include "sym_names.h"
#include "numeric.h"
#include "numeric_complex.h"
#include "arithmetic.h"
#include "poly/groebner.h"
#include "poly/gb_flint.h"  /* gb_flint_quotient */
#include "poly/poly.h"        /* is_polynomial, get_degree_poly */
#include "linalg/eigen.h"
#include "linalg/lapack.h"  /* mat_lapack_dgeev: double eigenvector seeds */
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

/* These monomial helpers serve only the MPFR polynomial-system solver below, so
 * they are compiled with MPFR only; otherwise they are unused
 * (-Werror=unused-function). */
#ifdef USE_MPFR
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
#endif

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
        && e->data.function.head->data.symbol.name == SYM_Complex
        && e->data.function.arg_count == 2) {
        double re = value_magnitude(e->data.function.args[0]);
        double im = value_magnitude(e->data.function.args[1]);
        return hypot(re, im);
    }
    return INFINITY; /* unevaluated / symbolic -> treat as non-zero */
}

#ifdef USE_MPFR
/* Build a Complex[Real,Real] (or Real) Expr from an ncpx at machine precision,
 * for substitution into residuals. */
static Expr* ncpx_to_expr_machine(const ncpx* z) {
    double re = mpfr_get_d(z->re, MPFR_RNDN);
    double im = mpfr_get_d(z->im, MPFR_RNDN);
    if (im == 0.0) return expr_new_real(re);
    return expr_new_function(expr_new_symbol(SYM_Complex),
               (Expr*[]){ expr_new_real(re), expr_new_real(im) }, 2);
}
#endif

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

#ifdef USE_MPFR
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

        /* x_var * B_j itself standard: its normal form is itself. */
        int self = mono_index(B, d, nvar, e);
        if (self >= 0) { mpq_set_ui(M[(size_t)self * d + j], 1, 1); continue; }

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

/* Deterministic small pseudo-random coefficient stream (seeded). A 64-bit
 * LCG, so the state is uint64_t: `unsigned long` is only 32 bits on LLP64
 * Windows, where the `>> 33` would overflow (issue #92). */
static long nsys_rand(uint64_t* state) {
    *state = *state * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);
    /* small nonzero coefficient in [1, 9] with alternating sign */
    long v = (long)((*state >> 33) % 9u) + 1;
    return ((*state >> 17) & 1u) ? v : -v;
}

/* ------------------------------------------------------------------ *
 *  The original system in evaluation form, for the Newton polish and the
 *  verification of each eigenvector candidate.
 * ------------------------------------------------------------------ */
/* cd: |c_t| as a double (the residual scale); cv: c_t rounded to a double
 * (the double-precision polish, sys_polish_d). */
typedef struct { int nt; int* ex; double* cd; double* cv; mpq_t* c; } SysP;

static SysP* sys_build(GBPoly* const* F, int nF, int nvar) {
    SysP* S = (SysP*)malloc(sizeof(SysP) * (size_t)(nF ? nF : 1));
    for (int i = 0; i < nF; i++) {
        size_t nt = F[i]->n_terms;
        S[i].nt = (int)nt;
        S[i].ex = (int*)malloc(sizeof(int) * (nt ? nt : 1) * (size_t)nvar);
        S[i].c = (mpq_t*)malloc(sizeof(mpq_t) * (nt ? nt : 1));
        S[i].cd = (double*)malloc(sizeof(double) * (nt ? nt : 1));
        S[i].cv = (double*)malloc(sizeof(double) * (nt ? nt : 1));
        memcpy(S[i].ex, F[i]->exps, sizeof(int) * nt * (size_t)nvar);
        for (size_t t = 0; t < nt; t++) {
            mpq_init(S[i].c[t]);
            mpq_set(S[i].c[t], F[i]->coefs[t]);
            S[i].cv[t] = mpq_get_d(F[i]->coefs[t]);
            S[i].cd[t] = fabs(S[i].cv[t]);
        }
    }
    return S;
}

static void sys_free(SysP* S, int nF) {
    for (int i = 0; i < nF; i++) {
        for (int t = 0; t < S[i].nt; t++) mpq_clear(S[i].c[t]);
        free(S[i].c); free(S[i].ex); free(S[i].cd); free(S[i].cv);
    }
    free(S);
}

/* f_i(x) and (J != NULL) the Jacobian J[i*nv + k] at x, in MPFR complex at
 * wp; scale[i] = sum |c_t| |x^e_t| (double), the size of the terms whose
 * cancellation f_i is. */
static void sys_eval(const SysP* S, int np, int nv, const ncpx* x, mpfr_prec_t wp,
                     ncpx* f, ncpx* J, double* scale) {
    int maxd = 0;
    for (int i = 0; i < np; i++)
        for (int t = 0; t < S[i].nt; t++)
            for (int k = 0; k < nv; k++)
                if (S[i].ex[(size_t)t * nv + k] > maxd) maxd = S[i].ex[(size_t)t * nv + k];
    int W = maxd + 1;
    ncpx* pw = (ncpx*)malloc(sizeof(ncpx) * (size_t)nv * (size_t)W);
    double* ax = (double*)malloc(sizeof(double) * (size_t)nv);
    for (int k = 0; k < nv; k++) {
        ncpx_init(&pw[k * W], wp);
        ncpx_set_ui(&pw[k * W], 1);
        for (int e = 1; e < W; e++) {
            ncpx_init(&pw[k * W + e], wp);
            ncpx_mul(&pw[k * W + e], &pw[k * W + e - 1], &x[k], wp);
        }
        /* |x_k| floored at 1: a coordinate that vanishes at the solution makes
         * every term of some equation vanish (Katsura: u1 = ... = u5 = 0), and
         * |f| / sum |terms| is then ~1 however accurate the point is.  The
         * floor measures such coordinates absolutely, as the duplicate test
         * does. */
        ax[k] = hypot(mpfr_get_d(x[k].re, MPFR_RNDN), mpfr_get_d(x[k].im, MPFR_RNDN));
        if (ax[k] < 1.0) ax[k] = 1.0;
    }
    ncpx m, cm; ncpx_init(&m, wp); ncpx_init(&cm, wp);
    mpfr_t cq; mpfr_init2(cq, wp);
    for (int i = 0; i < np; i++) {
        ncpx_set_ui(&f[i], 0);
        if (J) for (int k = 0; k < nv; k++) ncpx_set_ui(&J[(size_t)i * nv + k], 0);
        double sc = 0.0;
        for (int t = 0; t < S[i].nt; t++) {
            const int* e = S[i].ex + (size_t)t * nv;
            mpfr_set_q(cq, S[i].c[t], MPFR_RNDN);
            ncpx_set_ui(&m, 1);
            double am = S[i].cd[t];
            for (int k = 0; k < nv; k++) {
                if (e[k]) ncpx_mul(&m, &m, &pw[k * W + e[k]], wp);
                am *= pow(ax[k], (double)e[k]);
            }
            sc += am;
            mpfr_mul(cm.re, m.re, cq, MPFR_RNDN);
            mpfr_mul(cm.im, m.im, cq, MPFR_RNDN);
            ncpx_add(&f[i], &f[i], &cm);
            if (!J) continue;
            for (int k = 0; k < nv; k++) {
                if (!e[k]) continue;
                /* d/dx_k: e_k c x_k^(e_k-1) prod_{j != k} x_j^e_j */
                ncpx_set_ui(&m, 1);
                for (int j = 0; j < nv; j++) {
                    int ej = (j == k) ? e[j] - 1 : e[j];
                    if (ej) ncpx_mul(&m, &m, &pw[j * W + ej], wp);
                }
                mpfr_mul(cm.re, m.re, cq, MPFR_RNDN);
                mpfr_mul(cm.im, m.im, cq, MPFR_RNDN);
                mpfr_mul_si(cm.re, cm.re, e[k], MPFR_RNDN);
                mpfr_mul_si(cm.im, cm.im, e[k], MPFR_RNDN);
                ncpx_add(&J[(size_t)i * nv + k], &J[(size_t)i * nv + k], &cm);
            }
        }
        scale[i] = sc;
    }
    ncpx_clear(&m); ncpx_clear(&cm); mpfr_clear(cq);
    for (int k = 0; k < nv * W; k++) ncpx_clear(&pw[k]);
    free(pw); free(ax);
}

/* max_i |f_i| / scale_i from an evaluation. */
static double sys_relres(const ncpx* f, const double* scale, int np) {
    double worst = 0.0;
    for (int i = 0; i < np; i++) {
        double a = hypot(mpfr_get_d(f[i].re, MPFR_RNDN), mpfr_get_d(f[i].im, MPFR_RNDN));
        double r = scale[i] > 0.0 ? a / scale[i] : a;
        if (!(r == r)) return INFINITY;
        if (r > worst) worst = r;
    }
    return worst;
}

/* Solve J d = f (n x n complex) in place by partial-pivot elimination; the
 * solution overwrites f.  false when J is singular. */
static bool cpx_solve(ncpx* J, ncpx* f, int n, mpfr_prec_t wp) {
    ncpx t, q; ncpx_init(&t, wp); ncpx_init(&q, wp);
    mpfr_t a, b; mpfr_init2(a, wp); mpfr_init2(b, wp);
    bool ok = true;
    for (int c = 0; c < n && ok; c++) {
        int pr = c;
        ncpx_abs(a, &J[(size_t)c * n + c]);
        for (int r = c + 1; r < n; r++) {
            ncpx_abs(b, &J[(size_t)r * n + c]);
            if (mpfr_cmp(b, a) > 0) { mpfr_set(a, b, MPFR_RNDN); pr = r; }
        }
        if (mpfr_zero_p(a)) { ok = false; break; }
        if (pr != c) {
            for (int k = 0; k < n; k++) {
                ncpx_set(&t, &J[(size_t)pr * n + k]);
                ncpx_set(&J[(size_t)pr * n + k], &J[(size_t)c * n + k]);
                ncpx_set(&J[(size_t)c * n + k], &t);
            }
            ncpx_set(&t, &f[pr]); ncpx_set(&f[pr], &f[c]); ncpx_set(&f[c], &t);
        }
        for (int r = c + 1; r < n; r++) {
            ncpx_div(&q, &J[(size_t)r * n + c], &J[(size_t)c * n + c], wp);
            for (int k = c; k < n; k++) {
                ncpx_mul(&t, &q, &J[(size_t)c * n + k], wp);
                ncpx_sub(&J[(size_t)r * n + k], &J[(size_t)r * n + k], &t);
            }
            ncpx_mul(&t, &q, &f[c], wp);
            ncpx_sub(&f[r], &f[r], &t);
        }
    }
    for (int c = n - 1; c >= 0 && ok; c--) {
        for (int k = c + 1; k < n; k++) {
            ncpx_mul(&t, &J[(size_t)c * n + k], &f[k], wp);
            ncpx_sub(&f[c], &f[c], &t);
        }
        ncpx_div(&f[c], &f[c], &J[(size_t)c * n + c], wp);
    }
    ncpx_clear(&t); ncpx_clear(&q); mpfr_clear(a); mpfr_clear(b);
    return ok;
}

/* Newton on the original (square) system from the eigenvector estimate x,
 * in place, at wp.  Returns the final relative residual.  At a multiple
 * solution Newton only creeps; the residual then decides. */
static double sys_polish(const SysP* S, int n, ncpx* x, mpfr_prec_t wp) {
    ncpx* f = (ncpx*)malloc(sizeof(ncpx) * (size_t)n);
    ncpx* J = (ncpx*)malloc(sizeof(ncpx) * (size_t)n * (size_t)n);
    double* sc = (double*)malloc(sizeof(double) * (size_t)n);
    for (int i = 0; i < n; i++) ncpx_init(&f[i], wp);
    for (int i = 0; i < n * n; i++) ncpx_init(&J[i], wp);
    double floor_rel = pow(2.0, -(double)wp + 8.0);
    double rr = INFINITY;
    for (int it = 0; it < 12; it++) {
        sys_eval(S, n, n, x, wp, f, J, sc);
        rr = sys_relres(f, sc, n);
        if (rr <= floor_rel) break;
        if (!cpx_solve(J, f, n, wp)) break;
        double step = 0.0, xs = 1.0;
        for (int k = 0; k < n; k++) {
            ncpx_sub(&x[k], &x[k], &f[k]);
            double dk = hypot(mpfr_get_d(f[k].re, MPFR_RNDN), mpfr_get_d(f[k].im, MPFR_RNDN));
            double xk = hypot(mpfr_get_d(x[k].re, MPFR_RNDN), mpfr_get_d(x[k].im, MPFR_RNDN));
            if (dk > step) step = dk;
            if (xk > xs) xs = xk;
        }
        if (!(step == step)) { rr = INFINITY; break; }
        if (step <= floor_rel * xs) {
            sys_eval(S, n, n, x, wp, f, NULL, sc);
            rr = sys_relres(f, sc, n);
            break;
        }
    }
    for (int i = 0; i < n; i++) ncpx_clear(&f[i]);
    for (int i = 0; i < n * n; i++) ncpx_clear(&J[i]);
    free(f); free(J); free(sc);
    return rr;
}

/* ------------------------------------------------------------------ *
 *  Double-precision polish of a double seed
 * ------------------------------------------------------------------ *
 *
 * The LAPACK eigenvectors are double-precision seeds; polishing them in MPFR
 * from the start pays every Newton step at the working precision, with an
 * allocation per complex product, although only the last step needs it (a
 * quadratically convergent step from a double-accurate point lands at the
 * working precision).  So the seed is polished in double first: for a
 * machine-precision result that is the whole polish (the answer is rounded to
 * a double anyway, and the verification threshold 2^-26.5 is far above the
 * double residual floor); for a higher target the MPFR polish then starts
 * one step from its floor.  Complex numbers are re/im pairs. */

/* f_i(x) and J (n x n, row-major) at x; scale[i] as in sys_eval. */
static void sys_eval_d(const SysP* S, int n, const double* xr, const double* xi,
                       double* fr, double* fi, double* Jr, double* Ji,
                       double* scale, double* pwr, double* pwi, int W) {
    for (int k = 0; k < n; k++) {
        pwr[k * W] = 1.0; pwi[k * W] = 0.0;
        for (int e = 1; e < W; e++) {
            double ar = pwr[k * W + e - 1], ai = pwi[k * W + e - 1];
            pwr[k * W + e] = ar * xr[k] - ai * xi[k];
            pwi[k * W + e] = ar * xi[k] + ai * xr[k];
        }
    }
    for (int i = 0; i < n; i++) {
        double sr = 0.0, si = 0.0, sc = 0.0;
        for (int k = 0; k < n; k++) { Jr[i * n + k] = 0.0; Ji[i * n + k] = 0.0; }
        for (int t = 0; t < S[i].nt; t++) {
            const int* e = S[i].ex + (size_t)t * n;
            double mr = 1.0, mi = 0.0, am = S[i].cd[t];
            for (int k = 0; k < n; k++) {
                if (!e[k]) continue;
                double pr = pwr[k * W + e[k]], pi = pwi[k * W + e[k]];
                double tr = mr * pr - mi * pi;
                mi = mr * pi + mi * pr; mr = tr;
                double ax = hypot(xr[k], xi[k]);
                am *= pow(ax < 1.0 ? 1.0 : ax, (double)e[k]);
            }
            double c = S[i].cv[t];
            sr += c * mr; si += c * mi; sc += am;
            for (int k = 0; k < n; k++) {
                if (!e[k]) continue;
                double dr = 1.0, di = 0.0;
                for (int j = 0; j < n; j++) {
                    int ej = (j == k) ? e[j] - 1 : e[j];
                    if (!ej) continue;
                    double pr = pwr[j * W + ej], pi = pwi[j * W + ej];
                    double tr = dr * pr - di * pi;
                    di = dr * pi + di * pr; dr = tr;
                }
                Jr[i * n + k] += c * e[k] * dr;
                Ji[i * n + k] += c * e[k] * di;
            }
        }
        fr[i] = sr; fi[i] = si; scale[i] = sc;
    }
}

/* Solve J d = f in place (partial pivoting); false when J is singular. */
static bool cpx_solve_d(double* Jr, double* Ji, double* fr, double* fi, int n) {
    for (int c = 0; c < n; c++) {
        int pr = c;
        double best = hypot(Jr[c * n + c], Ji[c * n + c]);
        for (int r = c + 1; r < n; r++) {
            double a = hypot(Jr[r * n + c], Ji[r * n + c]);
            if (a > best) { best = a; pr = r; }
        }
        if (!(best > 0.0)) return false;
        if (pr != c) {
            for (int k = 0; k < n; k++) {
                double t = Jr[pr * n + k]; Jr[pr * n + k] = Jr[c * n + k]; Jr[c * n + k] = t;
                t = Ji[pr * n + k]; Ji[pr * n + k] = Ji[c * n + k]; Ji[c * n + k] = t;
            }
            double t = fr[pr]; fr[pr] = fr[c]; fr[c] = t;
            t = fi[pr]; fi[pr] = fi[c]; fi[c] = t;
        }
        double pr_ = Jr[c * n + c], pi_ = Ji[c * n + c], den = pr_ * pr_ + pi_ * pi_;
        for (int r = c + 1; r < n; r++) {
            double ar = Jr[r * n + c], ai = Ji[r * n + c];
            double qr = (ar * pr_ + ai * pi_) / den, qi = (ai * pr_ - ar * pi_) / den;
            for (int k = c; k < n; k++) {
                double br = Jr[c * n + k], bi = Ji[c * n + k];
                Jr[r * n + k] -= qr * br - qi * bi;
                Ji[r * n + k] -= qr * bi + qi * br;
            }
            fr[r] -= qr * fr[c] - qi * fi[c];
            fi[r] -= qr * fi[c] + qi * fr[c];
        }
    }
    for (int c = n - 1; c >= 0; c--) {
        double sr = fr[c], si = fi[c];
        for (int k = c + 1; k < n; k++) {
            sr -= Jr[c * n + k] * fr[k] - Ji[c * n + k] * fi[k];
            si -= Jr[c * n + k] * fi[k] + Ji[c * n + k] * fr[k];
        }
        double pr_ = Jr[c * n + c], pi_ = Ji[c * n + c], den = pr_ * pr_ + pi_ * pi_;
        fr[c] = (sr * pr_ + si * pi_) / den;
        fi[c] = (si * pr_ - sr * pi_) / den;
    }
    return true;
}

/* Newton in double from (xr, xi), in place; returns the final relative
 * residual (INFINITY on breakdown).  Stops at the double floor: a residual or
 * a step of a few ulps, or a step that no longer shrinks. */
static double sys_polish_d(const SysP* S, int n, double* xr, double* xi) {
    int maxd = 0;
    for (int i = 0; i < n; i++)
        for (int t = 0; t < S[i].nt; t++)
            for (int k = 0; k < n; k++)
                if (S[i].ex[(size_t)t * n + k] > maxd) maxd = S[i].ex[(size_t)t * n + k];
    int W = maxd + 1;
    size_t nn = (size_t)n * (size_t)n;
    double* buf = (double*)malloc(sizeof(double) * (5 * (size_t)n + 2 * nn + 2 * (size_t)n * W));
    double *fr = buf, *fi = fr + n, *sc = fi + n, *Jr = sc + n, *Ji = Jr + nn;
    double *pwr = Ji + nn, *pwi = pwr + (size_t)n * W;
    double floor_rel = 0x1p-50;
    double rr = INFINITY, prev_step = INFINITY;
    for (int it = 0; it < 16; it++) {
        sys_eval_d(S, n, xr, xi, fr, fi, Jr, Ji, sc, pwr, pwi, W);
        rr = 0.0;
        for (int i = 0; i < n; i++) {
            double a = hypot(fr[i], fi[i]);
            double r = sc[i] > 0.0 ? a / sc[i] : a;
            if (!(r == r)) { rr = INFINITY; break; }
            if (r > rr) rr = r;
        }
        if (rr <= floor_rel || rr == INFINITY) break;
        if (!cpx_solve_d(Jr, Ji, fr, fi, n)) break;
        double step = 0.0, xs = 1.0;
        for (int k = 0; k < n; k++) {
            xr[k] -= fr[k]; xi[k] -= fi[k];
            double dk = hypot(fr[k], fi[k]), xk = hypot(xr[k], xi[k]);
            if (dk > step) step = dk;
            if (xk > xs) xs = xk;
        }
        if (!(step == step)) { rr = INFINITY; break; }
        /* at the floor, or stagnating (a multiple root creeps linearly):
         * one last evaluation gives the residual of the final point */
        if (step <= floor_rel * xs || (it >= 4 && step >= 0.5 * prev_step)) {
            sys_eval_d(S, n, xr, xi, fr, fi, Jr, Ji, sc, pwr, pwi, W);
            rr = 0.0;
            for (int i = 0; i < n; i++) {
                double a = hypot(fr[i], fi[i]);
                double r = sc[i] > 0.0 ? a / sc[i] : a;
                if (!(r == r)) { rr = INFINITY; break; }
                if (r > rr) rr = r;
            }
            break;
        }
        prev_step = step;
    }
    free(buf);
    return rr;
}

/* Accepted solutions of one attempt: coordinates (n x nvar ncpx), the
 * multiplicity of each, the eigenvalue it came from, and the bookkeeping
 * that decides whether the attempt is complete. */
typedef struct {
    int n, cap, nvar, total, lost, rejected;
    ncpx* x; int* mult; double* lre; double* lim;
} NsysSols;

static void nsys_sols_clear(NsysSols* s) {
    for (int i = 0; i < s->n * s->nvar; i++) ncpx_clear(&s->x[i]);
    free(s->x); free(s->mult); free(s->lre); free(s->lim);
    memset(s, 0, sizeof *s);
}

/* One eigenvalue attempt at working precision wp with a fresh random linear
 * form.  Returns 0 when the eigen solver fails. */
static int nsys_attempt(mpq_t* Mx, const int* B, int d, int nvar, mpfr_prec_t wp,
                        uint64_t* st, const SysP* S, int nS,
                        Expr** polys, int npoly, Expr** vars, int verify,
                        long target_bits, int use_lapack, NsysSols* out) {
    (void)B; (void)nS;
    memset(out, 0, sizeof *out);
    out->nvar = nvar;
    long* cform = (long*)malloc(sizeof(long) * (size_t)nvar);
    /* Coefficients with ~30 random bits: small integers make l(p) = l(q) for
     * distinct solutions of a structured system (Katsura, cyclic) far too
     * likely, and a repeated eigenvalue of M_l has no reliable eigenvectors. */
    for (int i = 0; i < nvar; i++) {
        long a = nsys_rand(st), b = nsys_rand(st), c = nsys_rand(st);
        cform[i] = a * 100000000L + b * 10000L + c;   /* |cform| < 2^30 */
    }
    /* The LAPACK attempt (double seeds) works in double throughout up to the
     * polish: M_l, its eigenvectors and the coordinates read off them are
     * only double-accurate anyway, so building M_l in rationals and carrying
     * the eigenvectors and M_{x_i} at wp (d^2 nvar MPFR values) bought no
     * accuracy.  The MPFR attempts keep the exact M_l rounded once to wp. */
    double* MxD = NULL;         /* lapack: M_{x_i} in double, row-major     */
    double* Vd = NULL;          /* lapack: eigenvector k = Vd[2(k d + j)..] */
    double* wr = NULL; double* wi = NULL;
    mpfr_t *ev_re = NULL, *ev_im = NULL, *V_re = NULL, *V_im = NULL, *MxF = NULL;
    int est = -1;
    if (use_lapack) {
        /* Double-precision eigenvectors (LAPACK dgeev) as Newton seeds: the
         * polish on the original equations supplies the accuracy, so the
         * eigenproblem need not be solved at wp.  A seed that misses (an
         * ill-conditioned eigenvalue) shows up as a rejected or lost
         * candidate and the caller redoes the attempt in MPFR. */
        size_t dd = (size_t)d * (size_t)d;
        MxD = (double*)malloc(sizeof(double) * (size_t)nvar * dd);
        for (size_t t = 0; t < (size_t)nvar * dd; t++) MxD[t] = mpq_get_d(Mx[t]);
        double* A = (double*)malloc(sizeof(double) * dd);
        wr = (double*)malloc(sizeof(double) * (size_t)d);
        wi = (double*)malloc(sizeof(double) * (size_t)d);
        double* VR = (double*)malloc(sizeof(double) * dd);
        for (int r = 0; r < d; r++)
            for (int c = 0; c < d; c++) {
                double acc = 0.0;
                for (int v = 0; v < nvar; v++)
                    acc += ((double)cform[v] * 1e-9) * MxD[(size_t)v * dd + (size_t)r * d + c];
                A[(size_t)r + (size_t)c * d] = acc;      /* column-major */
            }
        if (mat_lapack_dgeev(d, A, d, wr, wi, VR, d) == 0) {
            est = 0;
            Vd = (double*)malloc(sizeof(double) * 2 * dd);
            for (int k = 0; k < d; k++) {
                if (wi[k] == 0.0) {
                    for (int j = 0; j < d; j++) {
                        Vd[2 * ((size_t)k * d + j)] = VR[(size_t)j + (size_t)k * d];
                        Vd[2 * ((size_t)k * d + j) + 1] = 0.0;
                    }
                } else if (k + 1 < d && wi[k] > 0.0) {
                    /* a conjugate pair: columns k (real) and k+1 (imaginary) */
                    for (int j = 0; j < d; j++) {
                        double re = VR[(size_t)j + (size_t)k * d];
                        double im = VR[(size_t)j + (size_t)(k + 1) * d];
                        Vd[2 * ((size_t)k * d + j)] = re;
                        Vd[2 * ((size_t)k * d + j) + 1] = im;
                        Vd[2 * ((size_t)(k + 1) * d + j)] = re;
                        Vd[2 * ((size_t)(k + 1) * d + j) + 1] = -im;
                    }
                    k++;
                } else {
                    est = -1; break;                     /* unpaired: give up */
                }
            }
        }
        free(A); free(VR);
        free(cform);
        if (est != 0) { free(MxD); free(Vd); free(wr); free(wi); return 0; }
    } else {
        mpfr_t* Ml = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)d * (size_t)d);
        for (size_t t = 0; t < (size_t)d * d; t++) mpfr_init2(Ml[t], wp);
        {
            mpq_t acc, tmp; mpq_init(acc); mpq_init(tmp);
            for (size_t cell = 0; cell < (size_t)d * d; cell++) {
                mpq_set_ui(acc, 0, 1);
                for (int i = 0; i < nvar; i++) {
                    mpq_set_si(tmp, cform[i], 1000000000UL);
                    mpq_canonicalize(tmp);
                    mpq_mul(tmp, tmp, Mx[(size_t)i * d * d + cell]);
                    mpq_add(acc, acc, tmp);
                }
                mpfr_set_q(Ml[cell], acc, MPFR_RNDN);
            }
            mpq_clear(acc); mpq_clear(tmp);
        }
        free(cform);
        ev_re = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)d);
        ev_im = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)d);
        V_re  = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)d * (size_t)d);
        V_im  = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)d * (size_t)d);
        for (int i = 0; i < d; i++) { mpfr_init2(ev_re[i], wp); mpfr_init2(ev_im[i], wp); }
        for (size_t t = 0; t < (size_t)d * d; t++) { mpfr_init2(V_re[t], wp); mpfr_init2(V_im[t], wp); }
        est = eigen_all_eigenvectors_real_mpfr(Ml, (size_t)d, wp, ev_re, ev_im, V_re, V_im);
        for (size_t t = 0; t < (size_t)d * d; t++) mpfr_clear(Ml[t]);
        free(Ml);
        if (est != 0) {
            for (int i = 0; i < d; i++) { mpfr_clear(ev_re[i]); mpfr_clear(ev_im[i]); }
            for (size_t t = 0; t < (size_t)d * d; t++) { mpfr_clear(V_re[t]); mpfr_clear(V_im[t]); }
            free(ev_re); free(ev_im); free(V_re); free(V_im);
            return 0;
        }
        MxF = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)nvar * (size_t)d * (size_t)d);
        for (size_t t = 0; t < (size_t)nvar * d * d; t++) {
            mpfr_init2(MxF[t], wp);
            mpfr_set_q(MxF[t], Mx[t], MPFR_RNDN);
        }
    }
    double tol = pow(2.0, -0.5 * (double)target_bits);   /* relative residual */
    double old_tol = pow(10.0, -0.5 * numeric_bits_to_digits(target_bits));
    if (old_tol < 1e-6) old_tol = 1e-6;
    double same_tol = 4.0 * tol;                          /* duplicate point  */

    out->cap = d;
    out->x = (ncpx*)malloc(sizeof(ncpx) * (size_t)d * (size_t)nvar);
    out->mult = (int*)malloc(sizeof(int) * (size_t)d);
    out->lre = (double*)malloc(sizeof(double) * (size_t)d);
    out->lim = (double*)malloc(sizeof(double) * (size_t)d);

    ncpx* vk = (ncpx*)malloc(sizeof(ncpx) * (size_t)d);
    for (int t = 0; t < d; t++) ncpx_init(&vk[t], wp);
    ncpx* coord = (ncpx*)malloc(sizeof(ncpx) * (size_t)nvar);
    for (int i = 0; i < nvar; i++) ncpx_init(&coord[i], wp);
    ncpx prod, pivot, num;
    ncpx_init(&prod, wp); ncpx_init(&pivot, wp); ncpx_init(&num, wp);
    mpfr_t bestm, mag; mpfr_init2(bestm, wp); mpfr_init2(mag, wp);
    double* xdr = (double*)malloc(sizeof(double) * 2 * (size_t)nvar);
    double* xdi = xdr + nvar;

    for (int k = 0; k < d; k++) {
        /* x_i(p_k) = (M_{x_i} v_k)[j*] / v_k[j*], j* the largest component */
        if (use_lapack) {
            const double* v = Vd + 2 * (size_t)k * d;
            int jstar = 0; double best = -1.0;
            for (int j = 0; j < d; j++) {
                double m = hypot(v[2 * j], v[2 * j + 1]);
                if (m > best) { best = m; jstar = j; }
            }
            double pr = v[2 * jstar], pi = v[2 * jstar + 1], den = pr * pr + pi * pi;
            for (int i = 0; i < nvar; i++) {
                const double* row = MxD + (size_t)i * d * d + (size_t)jstar * d;
                double nr = 0.0, ni = 0.0;
                for (int m = 0; m < d; m++) { nr += row[m] * v[2 * m]; ni += row[m] * v[2 * m + 1]; }
                xdr[i] = (nr * pr + ni * pi) / den;
                xdi[i] = (ni * pr - nr * pi) / den;
                mpfr_set_d(coord[i].re, xdr[i], MPFR_RNDN);
                mpfr_set_d(coord[i].im, xdi[i], MPFR_RNDN);
            }
        } else {
            int jstar = 0;
            mpfr_set_si(bestm, -1, MPFR_RNDN);
            for (int j = 0; j < d; j++) {
                mpfr_set(vk[j].re, V_re[(size_t)k * d + j], MPFR_RNDN);
                mpfr_set(vk[j].im, V_im[(size_t)k * d + j], MPFR_RNDN);
                mpfr_hypot(mag, vk[j].re, vk[j].im, MPFR_RNDN);
                if (mpfr_cmp(mag, bestm) > 0) { mpfr_set(bestm, mag, MPFR_RNDN); jstar = j; }
            }
            mpfr_set(pivot.re, vk[jstar].re, MPFR_RNDN);
            mpfr_set(pivot.im, vk[jstar].im, MPFR_RNDN);
            for (int i = 0; i < nvar; i++) {
                mpfr_set_zero(num.re, 1); mpfr_set_zero(num.im, 1);
                const mpfr_t* Mi = MxF + (size_t)i * d * d;
                for (int m = 0; m < d; m++) {
                    mpfr_mul(prod.re, Mi[(size_t)jstar * d + m], vk[m].re, MPFR_RNDN);
                    mpfr_mul(prod.im, Mi[(size_t)jstar * d + m], vk[m].im, MPFR_RNDN);
                    mpfr_add(num.re, num.re, prod.re, MPFR_RNDN);
                    mpfr_add(num.im, num.im, prod.im, MPFR_RNDN);
                }
                ncpx_div(&coord[i], &num, &pivot, wp);
            }
        }

        /* Polish on the original equations, then verify: relative residual
         * (the size of f against the terms it cancels) for a square system;
         * the absolute residual for an overdetermined one. */
        bool pass;
        if (S) {
            double rr = INFINITY;
            bool done = false;
            if (use_lapack) {
                /* a double seed: polish it in double first (sys_polish_d) */
                double rd = sys_polish_d(S, nvar, xdr, xdi);
                if (rd <= 0x1p-40) {
                    for (int i = 0; i < nvar; i++) {
                        mpfr_set_d(coord[i].re, xdr[i], MPFR_RNDN);
                        mpfr_set_d(coord[i].im, xdi[i], MPFR_RNDN);
                    }
                    /* a machine result: the double polish is the polish */
                    if (target_bits <= 53) { rr = rd; done = true; }
                }
            }
            if (!done) rr = sys_polish(S, nvar, coord, wp);
            pass = !verify || rr <= tol;
        } else {
            double resid = verify ? system_residual(polys, npoly, vars, nvar, coord) : 0.0;
            pass = !verify || (resid <= old_tol && resid == resid);
        }
        if (!pass) { out->rejected++; continue; }

        double lr = use_lapack ? wr[k] : mpfr_get_d(ev_re[k], MPFR_RNDN);
        double li = use_lapack ? wi[k] : mpfr_get_d(ev_im[k], MPFR_RNDN);
        int dup = -1;
        for (int s = 0; s < out->n && dup < 0; s++) {
            double worst = 0.0;
            for (int i = 0; i < nvar; i++) {
                const ncpx* y = out->x + (size_t)s * nvar + i;
                double dr = mpfr_get_d(coord[i].re, MPFR_RNDN) - mpfr_get_d(y->re, MPFR_RNDN);
                double di = mpfr_get_d(coord[i].im, MPFR_RNDN) - mpfr_get_d(y->im, MPFR_RNDN);
                double sy = hypot(mpfr_get_d(y->re, MPFR_RNDN), mpfr_get_d(y->im, MPFR_RNDN));
                double rel = hypot(dr, di) / (sy > 1.0 ? sy : 1.0);
                if (rel > worst) worst = rel;
            }
            if (worst <= same_tol) dup = s;
        }
        if (dup >= 0) {
            /* A multiple solution: its eigenvalues cluster.  Two
             * well-separated eigenvalues on one point: a solution lost. */
            double dl = hypot(lr - out->lre[dup], li - out->lim[dup]);
            double sl = hypot(lr, li);
            if (dl <= 0x1p-16 * (1.0 + sl)) { out->mult[dup]++; out->total++; }
            else out->lost++;
            continue;
        }
        int s = out->n++;
        for (int i = 0; i < nvar; i++) {
            ncpx_init(&out->x[(size_t)s * nvar + i], wp);
            ncpx_set(&out->x[(size_t)s * nvar + i], &coord[i]);
        }
        out->mult[s] = 1; out->lre[s] = lr; out->lim[s] = li;
        out->total++;
    }

    mpfr_clear(bestm); mpfr_clear(mag);
    free(xdr);
    for (int t = 0; t < d; t++) ncpx_clear(&vk[t]);
    free(vk);
    for (int i = 0; i < nvar; i++) ncpx_clear(&coord[i]);
    free(coord);
    ncpx_clear(&prod); ncpx_clear(&pivot); ncpx_clear(&num);
    if (use_lapack) {
        free(MxD); free(Vd); free(wr); free(wi);
    } else {
        for (size_t t = 0; t < (size_t)nvar * d * d; t++) mpfr_clear(MxF[t]);
        free(MxF);
        for (int i = 0; i < d; i++) { mpfr_clear(ev_re[i]); mpfr_clear(ev_im[i]); }
        for (size_t t = 0; t < (size_t)d * d; t++) { mpfr_clear(V_re[t]); mpfr_clear(V_im[t]); }
        free(ev_re); free(ev_im); free(V_re); free(V_im);
    }
    return 1;
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

    /* 2-5, fast path: the quotient ring (standard monomials and exact
     * multiplication matrices) straight from the FLINT DRL engine, instead of
     * a GBPoly basis and one rational reduction per border monomial (that
     * reduction was a third of a Katsura-5 solve). */
    int* B = NULL; int d = 0;
    mpq_t* Mx = NULL;
    GBPoly** Fsys = F; int nFsys = nF;
    if (gb_flint_quotient(F, (size_t)nF, NSYS_MAX_DIM, &B, &d, &Mx)) {
        if (d == 0) {                                 /* <1>: no solutions */
            for (int i = 0; i < nFsys; i++) gb_poly_free(Fsys[i]);
            free(Fsys);
            return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
        }
        goto eigen_stage;
    }

    /* 2. Gröbner basis. */
    size_t nG = 0;
    GBPoly** G = gb_buchberger(F, (size_t)nF, &nG);
    /* The input polynomials stay alive: the candidates are polished and
     * verified against them (sys_build). */
#define NSYS_FREE_F() do { for (int i_ = 0; i_ < nFsys; i_++) gb_poly_free(Fsys[i_]); free(Fsys); } while (0)
    if (!G || nG == 0) { if (G) gb_basis_free(G, nG); NSYS_FREE_F(); return NULL; }

    /* Inconsistent system: GB contains a nonzero constant -> no solutions. */
    for (size_t i = 0; i < nG; i++) {
        if (gb_poly_is_constant(G[i]) && !gb_poly_is_zero(G[i])) {
            gb_basis_free(G, nG);
            NSYS_FREE_F();
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
        NSYS_FREE_F();
        return NULL;
    }

    /* 4. Enumerate the standard-monomial basis B (box minus LM multiples). */
    int Bcap = 0;
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
        if (overflow || d == 0) { free(B); free(bound); gb_basis_free(G, nG); NSYS_FREE_F(); return NULL; }
    }
    free(bound);

    /* 5. Multiplication matrices M_{x_i} (rational), and a generic M_l. */
    Mx = (mpq_t*)malloc(sizeof(mpq_t) * (size_t)nvar * (size_t)d * (size_t)d);
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
        NSYS_FREE_F();
        return NULL;
    }
#undef NSYS_FREE_F

eigen_stage:
    /* 6-7. Eigenvectors of a generic M_l, coordinates, Newton polish on the
     * original system, verification and multiplicity accounting
     * (nsys_attempt).  An attempt that rejects a candidate, or in which two
     * candidates from well-separated eigenvalues land on one solution (a
     * solution lost: the eigenvector of an ill-conditioned eigenvalue
     * pointed at a neighbour), is redone with another linear form at twice
     * the precision; the best attempt is kept.  The quotient dimension d is
     * the solution count with multiplicity, so the target is known. */
    long target_bits = want_machine ? 53 : bits;
    long guard = target_bits / 2; if (guard < 32) guard = 32;
    mpfr_prec_t wp = (mpfr_prec_t)(target_bits + guard);
    uint64_t st = seed ? (uint64_t)seed : 1234u;
    SysP* S = (npoly == nvar) ? sys_build(Fsys, nFsys, nvar) : NULL;

    NsysSols best; memset(&best, 0, sizeof best);
    for (int attempt = 0; attempt < 5; attempt++) {
        NsysSols cur;
        /* First attempt: LAPACK seeds (square systems only: they need the
         * polish); then MPFR eigenvectors at wp, 2 wp, 4 wp, 8 wp. */
        int lap = (attempt == 0 && S != NULL);
        if (attempt == 0 && !lap) continue;
        int ok = nsys_attempt(Mx, B, d, nvar, wp, &st, S, nFsys, polys, npoly, vars,
                              verify, target_bits, lap, &cur);
        if (ok && cur.total > best.total) { nsys_sols_clear(&best); best = cur; }
        else if (ok) nsys_sols_clear(&cur);
        if (best.total == d && best.lost == 0 && best.rejected == 0) break;
        if (attempt > 0) wp *= 2;
    }
    if (S) sys_free(S, nFsys);
    for (int i = 0; i < nFsys; i++) gb_poly_free(Fsys[i]);
    free(Fsys);
    for (size_t t = 0; t < (size_t)nvar * d * d; t++) mpq_clear(Mx[t]);
    free(Mx);
    free(B);

    /* No accepted solution though the quotient ring is non-trivial: signal
     * failure (NULL) so the caller can fall back to the elimination method
     * (rather than reporting a false empty set). */
    if (best.n == 0 && d > 0) { nsys_sols_clear(&best); return NULL; }

    double verify_tol = pow(2.0, -0.5 * (double)target_bits);
    Expr** sols = (Expr**)malloc(sizeof(Expr*) * (size_t)(best.total ? best.total : 1));
    int nsol = 0;
    for (int s = 0; s < best.n && (max_roots < 0 || nsol < max_roots); s++) {
        const ncpx* coord = best.x + (size_t)s * nvar;
        if (reals_only) {
            bool real_ok = true;
            for (int i = 0; i < nvar; i++) {
                double im = fabs(mpfr_get_d(coord[i].im, MPFR_RNDN));
                double re = fabs(mpfr_get_d(coord[i].re, MPFR_RNDN));
                if (im > verify_tol * (1.0 + re)) { real_ok = false; break; }
            }
            if (!real_ok) continue;
        }
        for (int rep = 0; rep < best.mult[s] && (max_roots < 0 || nsol < max_roots); rep++) {
            Expr** rules = (Expr**)malloc(sizeof(Expr*) * (size_t)nvar);
            for (int i = 0; i < nvar; i++) {
                Expr* val;
                double re = mpfr_get_d(coord[i].re, MPFR_RNDN);
                double im = mpfr_get_d(coord[i].im, MPFR_RNDN);
                bool is_real = reals_only || fabs(im) <= verify_tol * (1.0 + fabs(re));
                if (want_machine) {
                    if (is_real) val = expr_new_real(re);
                    else val = expr_new_function(expr_new_symbol(SYM_Complex),
                                  (Expr*[]){ expr_new_real(re), expr_new_real(im) }, 2);
                } else {
                    mpfr_t rr, ii; mpfr_init2(rr, target_bits); mpfr_init2(ii, target_bits);
                    mpfr_set(rr, coord[i].re, MPFR_RNDN);
                    if (is_real) mpfr_set_zero(ii, 1);
                    else         mpfr_set(ii, coord[i].im, MPFR_RNDN);
                    val = numeric_mpfr_make_complex(rr, ii);
                    mpfr_clear(rr); mpfr_clear(ii);
                }
                rules[i] = expr_new_function(expr_new_symbol(SYM_Rule),
                               (Expr*[]){ expr_copy(vars[i]), val }, 2);
            }
            sols[nsol++] = expr_new_function(expr_new_symbol(SYM_List), rules, (size_t)nvar);
            free(rules);
        }
    }
    nsys_sols_clear(&best);

    Expr* out = expr_new_function(expr_new_symbol(SYM_List), sols, (size_t)nsol);
    free(sols);
    return out;
}
#else /* !USE_MPFR */
/* Without MPFR the ncpx (arbitrary-precision complex) eigenvalue engine is
 * unavailable; return NULL so the caller falls back to the elimination method. */
Expr* nsolve_polynomial_system(Expr** polys, int npoly, Expr** vars, int nvar,
                               NSysMethod method, bool reals_only,
                               bool want_machine, long bits, long max_roots,
                               int verify, unsigned long seed) {
    (void)polys; (void)npoly; (void)vars; (void)nvar; (void)method;
    (void)reals_only; (void)want_machine; (void)bits; (void)max_roots;
    (void)verify; (void)seed;
    return NULL;
}
#endif /* USE_MPFR */

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
        && r->data.function.head->data.symbol.name == SYM_Or) {
        eqs = r->data.function.args; neq = r->data.function.arg_count;
    } else { single[0] = r; eqs = single; neq = 1; }
    for (size_t i = 0; i < neq; i++) {
        Expr* e = eqs[i];
        if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
            && e->data.function.head->data.symbol.name == SYM_Equal
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
        && v->data.function.head->data.symbol.name == SYM_Complex;
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
