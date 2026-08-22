/*
 * solveint_linear.c
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


static bool si_linear_detect(const MPoly* eq, int n, mpz_t* a, mpz_t b);

/* Fraction-free (Bareiss) determinant of an m x m integer matrix.  M is
 * destroyed; `out` (pre-init'd) receives the exact determinant. */
static void si_det_bareiss(mpz_t** M, int m, mpz_t out) {
    if (m == 0) { mpz_set_ui(out, 1); return; }
    mpz_t prev, t1, t2; mpz_init_set_ui(prev, 1); mpz_init(t1); mpz_init(t2);
    int sign = 1;
    bool zero = false;
    for (int k = 0; k < m; k++) {
        if (mpz_sgn(M[k][k]) == 0) {                 /* pivot: swap in a nonzero row */
            int r = -1;
            for (int i = k + 1; i < m; i++) if (mpz_sgn(M[i][k]) != 0) { r = i; break; }
            if (r < 0) { zero = true; break; }
            for (int j = 0; j < m; j++) mpz_swap(M[k][j], M[r][j]);
            sign = -sign;
        }
        for (int i = k + 1; i < m; i++)
            for (int j = k + 1; j < m; j++) {
                mpz_mul(t1, M[i][j], M[k][k]);
                mpz_mul(t2, M[i][k], M[k][j]);
                mpz_sub(t1, t1, t2);
                mpz_divexact(M[i][j], t1, prev);     /* Bareiss: exact division */
            }
        mpz_set(prev, M[k][k]);
    }
    if (zero) mpz_set_ui(out, 0);
    else { mpz_set(out, M[m - 1][m - 1]); if (sign < 0) mpz_neg(out, out); }
    mpz_clear(prev); mpz_clear(t1); mpz_clear(t2);
}


/* --- Homogeneous linear SYSTEM with positivity -> parametric ray. ---
 *
 * A system of  m = n-1  homogeneous linear equations in n positive unknowns has
 * (generically) a 1-dimensional integer kernel spanned by a primitive vector v
 * (the generalised cross product: v_j = (-1)^j * det(A with column j removed)).
 * If v (or -v) is entirely positive the positive solutions are exactly
 * { k v : k = C[1] >= 1 }; if v has mixed signs the positive orthant meets the
 * kernel only at 0, so there are no positive solutions.  This settles rational
 * linear systems like  w == 5/6 x + y && x == 9/20 y + z && y == 13/42 z + w
 * (denominators cleared by the MPoly conversion). */
Expr* si_solve_linear_system_ray(SICtx* c) {
    int n = c->n;
    if (c->neq != n - 1 || n < 2 || c->n_ord != 0 || c->n_neq != 0 || c->n_abs_ord != 0 || !c->all_captured) return NULL;
    for (int i = 0; i < n; i++)                       /* every variable strictly positive */
        if (!(c->has_lo[i] && c->lo[i] >= 1)) return NULL;

    int m = c->neq;
    /* Build the m x n integer coefficient matrix; require every equation to be
     * linear and HOMOGENEOUS (constant term 0). */
    mpz_t* a = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n);
    mpz_t bb; for (int i = 0; i < n; i++) mpz_init(a[i]);
    mpz_init(bb);
    mpz_t** A = (mpz_t**)malloc(sizeof(mpz_t*) * (size_t)m);
    for (int q = 0; q < m; q++) { A[q] = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n); for (int i = 0; i < n; i++) mpz_init(A[q][i]); }
    bool ok = true;
    for (int q = 0; q < m && ok; q++) {
        if (!si_linear_detect(c->eq[q], n, a, bb) || mpz_sgn(bb) != 0) ok = false;
        else for (int i = 0; i < n; i++) mpz_set(A[q][i], a[i]);
    }

    Expr* result = NULL;
    if (ok) {
        /* v_j = (-1)^j det(A without column j). */
        mpz_t* v = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n);
        for (int i = 0; i < n; i++) mpz_init(v[i]);
        mpz_t** Msub = (mpz_t**)malloc(sizeof(mpz_t*) * (size_t)m);
        for (int q = 0; q < m; q++) { Msub[q] = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)m); for (int i = 0; i < m; i++) mpz_init(Msub[q][i]); }
        mpz_t det; mpz_init(det);
        for (int j = 0; j < n; j++) {
            int cc = 0;
            for (int col = 0; col < n; col++) { if (col == j) continue;
                for (int q = 0; q < m; q++) mpz_set(Msub[q][cc], A[q][col]);
                cc++;
            }
            si_det_bareiss(Msub, m, det);
            mpz_set(v[j], det);
            if (j & 1) mpz_neg(v[j], v[j]);
        }
        /* Primitive + sign-normalise (make the first nonzero positive). */
        mpz_t g; mpz_init_set_ui(g, 0);
        for (int i = 0; i < n; i++) mpz_gcd(g, g, v[i]);
        bool nonzero = mpz_sgn(g) != 0;
        if (nonzero) {
            for (int i = 0; i < n; i++) mpz_divexact(v[i], v[i], g);
            int lead = 0; while (lead < n && mpz_sgn(v[lead]) == 0) lead++;
            if (lead < n && mpz_sgn(v[lead]) < 0) for (int i = 0; i < n; i++) mpz_neg(v[i], v[i]);
        }
        bool all_pos = nonzero;
        for (int i = 0; i < n && all_pos; i++) if (mpz_sgn(v[i]) <= 0) all_pos = false;

        if (!nonzero) {
            result = NULL;                            /* degenerate kernel: decline */
        } else if (!all_pos) {
            result = mk_list(NULL, 0);                /* no positive solution: {} */
        } else {
            Expr** rules = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
            for (int i = 0; i < n; i++) {
                Expr* term = mpz_cmp_ui(v[i], 1) == 0
                    ? mk_fn1("C", mk_int(1))
                    : mk_fn2("Times", mk_mpz(v[i]), mk_fn1("C", mk_int(1)));
                Expr* cond = mk_fn2("GreaterEqual", mk_fn1("C", mk_int(1)), mk_int(1));
                rules[i] = mk_rule(expr_copy(c->var[i]),
                                   mk_fn2("ConditionalExpression", term, cond));
            }
            Expr* tuple = mk_list(rules, (size_t)n); free(rules);
            result = eval_and_free(mk_list((Expr*[]){ tuple }, 1));
        }
        mpz_clear(det); mpz_clear(g);
        for (int q = 0; q < m; q++) { for (int i = 0; i < m; i++) mpz_clear(Msub[q][i]); free(Msub[q]); }
        free(Msub);
        for (int i = 0; i < n; i++) mpz_clear(v[i]);
        free(v);
    }

    for (int q = 0; q < m; q++) { for (int i = 0; i < n; i++) mpz_clear(A[q][i]); free(A[q]); }
    free(A);
    for (int i = 0; i < n; i++) mpz_clear(a[i]);
    free(a); mpz_clear(bb);
    return result;
}


/* --- Linear Diophantine: parametric (unbounded) family. --- */

/* Detect a single linear equation  sum a_i x_i == b.  Fills a[i] (coef of
 * x_i) and b (= -constant term); false if any term has degree > 1. */
static bool si_linear_detect(const MPoly* eq, int n, mpz_t* a, mpz_t b) {
    for (int i = 0; i < n; i++) mpz_set_ui(a[i], 0);
    mpz_set_ui(b, 0);
    for (size_t t = 0; t < eq->n_terms; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        int nz = -1, deg = 0;
        for (int v = 0; v < n; v++) { deg += ex[v]; if (ex[v] > 0) { if (nz >= 0) return false; nz = v; } }
        if (deg > 1) return false;
        if (nz < 0) mpz_neg(b, eq->coefs[t]);           /* constant -> b = -const */
        else mpz_set(a[nz], eq->coefs[t]);
    }
    return true;
}


/* Particular solution x0[n] and homogeneous-lattice basis (n-1 vectors,
 * basis[j][i] = component i of vector j) of  sum a_i x_i == b, via the gcd
 * staircase.  Returns false (no solution) when gcd(a) does not divide b.  All
 * mpz arrays are pre-init'd by the caller. */
static bool si_linear_lattice(const mpz_t* a, int n, const mpz_t b,
                              mpz_t* x0, mpz_t** basis) {
    mpz_t g, s, t, ng, q1, q2, tmp;
    mpz_init(g); mpz_init(s); mpz_init(t); mpz_init(ng);
    mpz_init(q1); mpz_init(q2); mpz_init(tmp);
    mpz_t* r = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n);
    for (int i = 0; i < n; i++) mpz_init_set_ui(r[i], 0);

    mpz_set(g, a[0]); mpz_set_ui(r[0], 1);              /* g_0 = a_0, r_0 = e_0 */
    for (int i = 1; i < n; i++) {
        mpz_gcdext(ng, s, t, g, a[i]);                  /* ng = s*g + t*a_i */
        mpz_t* d = basis[i - 1];                        /* d = (a_i/ng) r - (g/ng) e_i */
        if (mpz_sgn(ng) != 0) { mpz_divexact(q1, a[i], ng); mpz_divexact(q2, g, ng); }
        else { mpz_set_ui(q1, 0); mpz_set_ui(q2, 0); }
        for (int k = 0; k < n; k++) mpz_mul(d[k], q1, r[k]);
        mpz_sub(d[i], d[i], q2);
        for (int k = 0; k < n; k++) mpz_mul(r[k], s, r[k]);   /* r = s*r + t*e_i */
        mpz_add(r[i], r[i], t);
        mpz_set(g, ng);
    }

    bool solvable;
    if (mpz_sgn(g) == 0) {                               /* all coefficients zero */
        solvable = (mpz_sgn(b) == 0);
        for (int i = 0; i < n; i++) mpz_set_ui(x0[i], 0);
    } else if (mpz_divisible_p(b, g)) {
        solvable = true;
        mpz_divexact(tmp, b, g);
        for (int i = 0; i < n; i++) mpz_mul(x0[i], tmp, r[i]);   /* x0 = (b/g) r */
    } else solvable = false;

    for (int i = 0; i < n; i++) mpz_clear(r[i]);
    free(r);
    mpz_clear(g); mpz_clear(s); mpz_clear(t); mpz_clear(ng);
    mpz_clear(q1); mpz_clear(q2); mpz_clear(tmp);
    return solvable;
}


/* A single linear equation with NO other constraints -> the full parametric
 * family  {{x_i -> x0_i + sum_j basis[j][i] C[j+1]}}, C[k] integer parameters.
 * Returns the owned result Expr, {} for an unsolvable equation, or NULL to
 * decline (not this shape). */
Expr* si_solve_linear_parametric(SICtx* c) {
    if (c->neq != 1 || c->n_ord != 0 || c->n_neq != 0 || c->n_abs_ord != 0 || !c->all_captured) return NULL;
    for (int i = 0; i < c->n; i++) if (c->has_lo[i] || c->has_hi[i]) return NULL;

    int n = c->n;
    mpz_t* a = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n);
    mpz_t* x0 = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n);
    for (int i = 0; i < n; i++) { mpz_init(a[i]); mpz_init(x0[i]); }
    mpz_t b; mpz_init(b);
    mpz_t** basis = (mpz_t**)malloc(sizeof(mpz_t*) * (size_t)(n > 1 ? n - 1 : 1));
    for (int j = 0; j < n - 1; j++) {
        basis[j] = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n);
        for (int i = 0; i < n; i++) mpz_init_set_ui(basis[j][i], 0);
    }

    Expr* result = NULL;
    if (!si_linear_detect(c->eq[0], n, a, b)) {
        result = NULL;
    } else if (!si_linear_lattice(a, n, b, x0, basis)) {
        result = mk_list(NULL, 0);                       /* provably no solution */
    } else {
        Expr** rules = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
        for (int i = 0; i < n; i++) {
            Expr** terms = (Expr**)malloc(sizeof(Expr*) * (size_t)(n + 1));
            int nt = 0;
            if (mpz_sgn(x0[i]) != 0) terms[nt++] = mk_mpz(x0[i]);
            for (int j = 0; j < n - 1; j++) {
                if (mpz_sgn(basis[j][i]) == 0) continue;
                Expr* Ck = expr_new_function(mk_sym("C"), (Expr*[]){ mk_int(j + 1) }, 1);
                terms[nt++] = mk_fn2("Times", mk_mpz(basis[j][i]), Ck);
            }
            Expr* val = (nt == 0) ? mk_int(0)
                      : (nt == 1) ? terms[0]
                      : expr_new_function(mk_sym("Plus"), terms, (size_t)nt);
            free(terms);
            val = eval_and_free(val);
            rules[i] = mk_rule(expr_copy(c->var[i]), val);
        }
        Expr* tuple = mk_list(rules, (size_t)n);
        free(rules);
        result = mk_list((Expr*[]){ tuple }, 1);
    }

    for (int j = 0; j < n - 1; j++) { for (int i = 0; i < n; i++) mpz_clear(basis[j][i]); free(basis[j]); }
    free(basis);
    for (int i = 0; i < n; i++) { mpz_clear(a[i]); mpz_clear(x0[i]); }
    free(a); free(x0); mpz_clear(b);
    return result;
}


/* --- General linear Diophantine SYSTEM (m >= 2 equations) via HNF. --- */

/* Solve the unconstrained integer linear system  A x = b  (m = c->neq linear
 * equations in n = c->n unknowns) completely, using the Hermite normal form.
 *
 * With  P A^T = R  (row HNF of A^T, P unimodular n x n, R n x m, rank rr), set
 * U = P^T and H = R^T so that  A U = H  and H is in column-echelon form.  The
 * substitution x = U y turns A x = b into H y = b; forward-substitution over
 * the rr pivot columns yields the pivot coordinates y_0..y_{rr-1} (each step an
 * exact-division / divisibility test -- a failure proves NO integer solution),
 * the free coordinates y_{rr}..y_{n-1} parameterise the kernel, and a check on
 * the non-pivot rows detects inconsistency.  The particular solution is
 * x0 = U y|_{free = 0}; the kernel basis is the free columns of U (= the free
 * rows of P).  Emits the family  {{x_i -> x0_i + sum_k Ker[k][i] C[k+1]}}, the
 * empty List for a provably unsolvable system, or NULL to decline (not the
 * shape: fewer than two equations, any bound / ordering / disequation present,
 * or a non-linear equation).  Single equations keep the dedicated
 * si_solve_linear_parametric path. */
Expr* si_solve_linear_system_hnf(SICtx* c) {
    int n = c->n, m = c->neq;
    if (m < 2 || c->n_ord != 0 || c->n_neq != 0 || c->n_abs_ord != 0 || !c->all_captured) return NULL;
    for (int i = 0; i < n; i++) if (c->has_lo[i] || c->has_hi[i]) return NULL;

    /* Build A (m x n) and b (m) from the linear equations; decline on any
     * nonlinear equation (si_linear_detect returns false). */
    mpz_t* A  = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)m * (size_t)n);
    mpz_t* bv = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)m);
    for (int i = 0; i < m * n; i++) mpz_init(A[i]);
    for (int i = 0; i < m; i++) mpz_init(bv[i]);
    mpz_t* arow = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n);
    for (int i = 0; i < n; i++) mpz_init(arow[i]);
    mpz_t brow; mpz_init(brow);

    bool shape_ok = true;
    for (int q = 0; q < m && shape_ok; q++) {
        if (!si_linear_detect(c->eq[q], n, arow, brow)) { shape_ok = false; break; }
        for (int i = 0; i < n; i++) mpz_set(A[(size_t)q * n + i], arow[i]);
        mpz_set(bv[q], brow);
    }
    for (int i = 0; i < n; i++) mpz_clear(arow[i]);
    free(arow); mpz_clear(brow);

    if (!shape_ok) {
        for (int i = 0; i < m * n; i++) mpz_clear(A[i]);
        for (int i = 0; i < m; i++) mpz_clear(bv[i]);
        free(A); free(bv);
        return NULL;
    }

    /* AT = A^T  (n x m). */
    mpz_t* AT = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n * (size_t)m);
    for (int i = 0; i < n * m; i++) mpz_init(AT[i]);
    for (int q = 0; q < m; q++)
        for (int i = 0; i < n; i++)
            mpz_set(AT[(size_t)i * m + q], A[(size_t)q * n + i]);

    mpz_t* R = NULL; mpz_t* P = NULL;
    int rr = linalg_hnf(AT, n, m, &R, &P);            /* P (n x n) * AT = R (n x m) */
    for (int i = 0; i < n * m; i++) mpz_clear(AT[i]);
    free(AT);
    if (rr < 0) {
        for (int i = 0; i < m * n; i++) mpz_clear(A[i]);
        for (int i = 0; i < m; i++) mpz_clear(bv[i]);
        free(A); free(bv);
        linalg_hnf_free(R, n * m); linalg_hnf_free(P, n * n);
        return NULL;
    }

    /* Pivot column of each HNF row k (k = 0..rr-1): the H-row where pivot
     * column k of H = R^T leads.  R row k's first nonzero column. */
    int* pivrow = (int*)malloc(sizeof(int) * (size_t)(rr > 0 ? rr : 1));
    for (int k = 0; k < rr; k++) {
        pivrow[k] = -1;
        for (int j = 0; j < m; j++)
            if (mpz_sgn(R[(size_t)k * m + j]) != 0) { pivrow[k] = j; break; }
    }

    /* Forward substitution:  y_k = ( b[pivrow[k]] - sum_{l<k} H[pivrow[k]][l] y_l )
     *                               / H[pivrow[k]][k],   H[i][l] = R[l][i]. */
    mpz_t* y = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)(rr > 0 ? rr : 1));
    for (int k = 0; k < rr; k++) mpz_init(y[k]);
    mpz_t acc, prod; mpz_init(acc); mpz_init(prod);
    bool solvable = true;
    for (int k = 0; k < rr && solvable; k++) {
        int pr = pivrow[k];
        mpz_set(acc, bv[pr]);
        for (int l = 0; l < k; l++) {
            mpz_mul(prod, R[(size_t)l * m + pr], y[l]);   /* H[pr][l] = R[l][pr] */
            mpz_sub(acc, acc, prod);
        }
        const mpz_t* piv = &R[(size_t)k * m + pr];        /* H[pr][k] = R[k][pr] > 0 */
        if (!mpz_divisible_p(acc, *piv)) { solvable = false; break; }
        mpz_divexact(y[k], acc, *piv);
    }

    /* Consistency over every row i: sum_{l<rr} H[i][l] y_l == b[i]. */
    for (int i = 0; i < m && solvable; i++) {
        mpz_set_ui(acc, 0);
        for (int l = 0; l < rr; l++) {
            mpz_mul(prod, R[(size_t)l * m + i], y[l]);
            mpz_add(acc, acc, prod);
        }
        if (mpz_cmp(acc, bv[i]) != 0) solvable = false;
    }

    Expr* result;
    if (!solvable) {
        result = mk_list(NULL, 0);                        /* provably no solution */
    } else {
        /* x0_i = sum_{k<rr} U[i][k] y_k = sum_k P[k][i] y_k.
         * Kernel vector for free coordinate f (rr <= f < n): comp i = P[f][i]. */
        int nfree = n - rr;
        Expr** rules = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
        for (int i = 0; i < n; i++) {
            mpz_t x0i; mpz_init(x0i);
            for (int k = 0; k < rr; k++) {
                mpz_mul(prod, P[(size_t)k * n + i], y[k]);
                mpz_add(x0i, x0i, prod);
            }
            Expr** terms = (Expr**)malloc(sizeof(Expr*) * (size_t)(nfree + 1));
            int nt = 0;
            if (mpz_sgn(x0i) != 0) terms[nt++] = mk_mpz(x0i);
            for (int f = 0; f < nfree; f++) {
                const mpz_t* comp = &P[(size_t)(rr + f) * n + i];   /* P[rr+f][i] */
                if (mpz_sgn(*comp) == 0) continue;
                Expr* Ck = expr_new_function(mk_sym("C"),
                                             (Expr*[]){ mk_int(f + 1) }, 1);
                terms[nt++] = mk_fn2("Times", mk_mpz(*comp), Ck);
            }
            Expr* val = (nt == 0) ? mk_int(0)
                      : (nt == 1) ? terms[0]
                      : expr_new_function(mk_sym("Plus"), terms, (size_t)nt);
            free(terms);
            val = eval_and_free(val);
            rules[i] = mk_rule(expr_copy(c->var[i]), val);
            mpz_clear(x0i);
        }
        Expr* tuple = mk_list(rules, (size_t)n);
        free(rules);
        result = mk_list((Expr*[]){ tuple }, 1);
    }

    for (int k = 0; k < rr; k++) mpz_clear(y[k]);
    free(y); free(pivrow);
    mpz_clear(acc); mpz_clear(prod);
    for (int i = 0; i < m * n; i++) mpz_clear(A[i]);
    for (int i = 0; i < m; i++) mpz_clear(bv[i]);
    free(A); free(bv);
    linalg_hnf_free(R, n * m); linalg_hnf_free(P, n * n);
    return result;
}


/* Invert an m x m matrix G (double) into Ginv via Gauss-Jordan.  Returns false
 * if singular. */
static bool si_mat_inverse(double G[SI_MAX_VARS][SI_MAX_VARS], int m,
                           double Ginv[SI_MAX_VARS][SI_MAX_VARS]) {
    double A[SI_MAX_VARS][2 * SI_MAX_VARS];
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++) { A[i][j] = G[i][j]; A[i][m + j] = (i == j) ? 1.0 : 0.0; }
    for (int col = 0; col < m; col++) {
        int piv = col; double best = 0;
        for (int r = col; r < m; r++) {
            double av = A[r][col] < 0 ? -A[r][col] : A[r][col];
            if (av > best) { best = av; piv = r; }
        }
        if (best < 1e-9) return false;
        for (int j = 0; j < 2 * m; j++) { double tmp = A[col][j]; A[col][j] = A[piv][j]; A[piv][j] = tmp; }
        double d = A[col][col];
        for (int j = 0; j < 2 * m; j++) A[col][j] /= d;
        for (int r = 0; r < m; r++) if (r != col) {
            double f = A[r][col];
            for (int j = 0; j < 2 * m; j++) A[r][j] -= f * A[col][j];
        }
    }
    for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) Ginv[i][j] = A[i][m + j];
    return true;
}


/* A single linear equation over a finite box.  gcd(a) not dividing b -> {}.
 * Otherwise the solution set is  x0 + Z-span of an (n-1)-vector lattice; the
 * lattice is LLL-reduced (LatticeReduce) and the coefficient box is the exact
 * projection of the value box under the pseudoinverse, so the search box is
 * small when the coefficients are large (few solutions) and provably too large
 * (declined) when the lattice is dense.  Every candidate is checked exactly.
 * Returns true when it settled the answer (emitting nothing means {}). */
bool si_solve_linear_bounded(SICtx* c, SearchState* st) {
    if (c->neq != 1) return false;
    int n = c->n, m = n - 1;
    if (m < 1 || n > SI_MAX_VARS) return false;
    for (int i = 0; i < n; i++) if (!(c->has_lo[i] && c->has_hi[i])) return false;

    mpz_t* a  = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n);
    mpz_t* x0 = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n);
    for (int i = 0; i < n; i++) { mpz_init(a[i]); mpz_init(x0[i]); }
    mpz_t b, g; mpz_init(b); mpz_init_set_ui(g, 0);
    mpz_t** basis = (mpz_t**)malloc(sizeof(mpz_t*) * (size_t)m);
    for (int j = 0; j < m; j++) { basis[j] = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n);
        for (int i = 0; i < n; i++) mpz_init_set_ui(basis[j][i], 0); }

    bool handled = false;
    if (si_linear_detect(c->eq[0], n, a, b)) {
        for (int i = 0; i < n; i++) mpz_gcd(g, g, a[i]);
        bool solvable = (mpz_sgn(g) == 0) ? (mpz_sgn(b) == 0)
                                          : mpz_divisible_p(b, g);
        if (!solvable) {
            handled = true;                              /* provably {} */
        } else if (si_linear_lattice(a, n, b, x0, basis)) {
            /* LLL-reduce the homogeneous basis via LatticeReduce[...]. */
            Expr** rowsE = (Expr**)malloc(sizeof(Expr*) * (size_t)m);
            for (int j = 0; j < m; j++) {
                Expr** comp = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
                for (int i = 0; i < n; i++) comp[i] = mk_mpz(basis[j][i]);
                rowsE[j] = mk_list(comp, (size_t)n); free(comp);
            }
            Expr* red = eval_and_free(expr_new_function(mk_sym("LatticeReduce"),
                                      (Expr*[]){ mk_list(rowsE, (size_t)m) }, 1));
            free(rowsE);

            mpz_t Bz[SI_MAX_VARS][SI_MAX_VARS];
            double Bd[SI_MAX_VARS][SI_MAX_VARS];
            for (int j = 0; j < SI_MAX_VARS; j++) for (int i = 0; i < SI_MAX_VARS; i++) mpz_init(Bz[j][i]);
            int mm = 0;
            if (red && red->type == EXPR_FUNCTION
                && red->data.function.head->type == EXPR_SYMBOL
                && red->data.function.head->data.symbol.name == SYM_List
                && (int)red->data.function.arg_count == m) {
                mm = m;
                for (int j = 0; j < m; j++) {
                    Expr* row = red->data.function.args[j];
                    for (int i = 0; i < n; i++) {
                        mpz_set(Bz[j][i], basis[j][i]);      /* fallback: unreduced */
                        if (row->type == EXPR_FUNCTION && (int)row->data.function.arg_count == n) {
                            int64_t vv;
                            if (expr_as_i64(row->data.function.args[i], &vv)) mpz_set_si(Bz[j][i], vv);
                        }
                        Bd[j][i] = mpz_get_d(Bz[j][i]);
                    }
                }
            }
            if (red) expr_free(red);

            /* P = (B B^T)^{-1} B  (m x n). */
            double Gm[SI_MAX_VARS][SI_MAX_VARS], Ginv[SI_MAX_VARS][SI_MAX_VARS];
            for (int j = 0; j < mm; j++) for (int k = 0; k < mm; k++) {
                double s = 0; for (int i = 0; i < n; i++) s += Bd[j][i] * Bd[k][i]; Gm[j][k] = s;
            }
            if (mm == m && si_mat_inverse(Gm, m, Ginv)) {
                /* Coefficient box = exact projection of the value box. */
                double Lb[SI_MAX_VARS], Ub[SI_MAX_VARS];
                for (int i = 0; i < n; i++) {
                    double x0d = mpz_get_d(x0[i]);
                    Lb[i] = (double)c->lo[i] - x0d;
                    Ub[i] = (double)c->hi[i] - x0d;
                }
                int64_t Clo[SI_MAX_VARS], Chi[SI_MAX_VARS];
                long double space = 1.0L; bool ok = true;
                for (int j = 0; j < m && ok; j++) {
                    double lo = 0, hi = 0;
                    for (int i = 0; i < n; i++) {
                        double P = 0; for (int k = 0; k < m; k++) P += Ginv[j][k] * Bd[k][i];
                        double t1 = P * Lb[i], t2 = P * Ub[i];
                        lo += (t1 < t2) ? t1 : t2;
                        hi += (t1 < t2) ? t2 : t1;
                    }
                    Clo[j] = (int64_t)floor(lo) - 1;
                    Chi[j] = (int64_t)ceil(hi) + 1;
                    space *= (long double)(Chi[j] - Clo[j] + 1);
                    if (space > (long double)SI_MAX_NODES) ok = false;
                }
                if (ok) {
                    st->max_visits = SI_MAX_NODES;
                    int64_t cc[SI_MAX_VARS];
                    for (int j = 0; j < m; j++) cc[j] = Clo[j];
                    mpz_t xi, term; mpz_init(xi); mpz_init(term);
                    for (;;) {
                        int64_t vals[SI_MAX_VARS]; bool inbox = true;
                        for (int i = 0; i < n && inbox; i++) {
                            mpz_set(xi, x0[i]);
                            for (int j = 0; j < m; j++) { mpz_mul_si(term, Bz[j][i], (long)cc[j]); mpz_add(xi, xi, term); }
                            if (!mpz_fits_slong_p(xi)) { inbox = false; break; }
                            int64_t xv = mpz_get_si(xi);
                            if (xv < c->lo[i] || xv > c->hi[i]) inbox = false;
                            vals[i] = xv;
                        }
                        if (inbox && si_verify(c, vals)) emit_full(st, vals);
                        int j = 0;
                        for (; j < m; j++) { if (++cc[j] <= Chi[j]) break; cc[j] = Clo[j]; }
                        if (j == m) break;
                    }
                    mpz_clear(xi); mpz_clear(term);
                    handled = true;
                }
            }
            for (int j = 0; j < SI_MAX_VARS; j++) for (int i = 0; i < SI_MAX_VARS; i++) mpz_clear(Bz[j][i]);
        }
    }

    for (int j = 0; j < m; j++) { for (int i = 0; i < n; i++) mpz_clear(basis[j][i]); free(basis[j]); }
    free(basis);
    for (int i = 0; i < n; i++) { mpz_clear(a[i]); mpz_clear(x0[i]); }
    free(a); free(x0); mpz_clear(b); mpz_clear(g);
    return handled;
}
