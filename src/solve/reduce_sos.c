/*
 * reduce_sos.c -- Positivstellensatz / sum-of-squares emptiness certificate.
 * See reduce_sos.h for the contract.  D5 of tasks/reduce_deficiencies.md.
 *
 * Pipeline (all sound -- a wrong certificate can never be produced because the
 * final check is exact):
 *   1. Extract {g_i >= 0}, {h_j == 0} and a strict target q < 0 from the DNF
 *      conjunction.  Prove q >= 0 on {g_i >= 0, h_j == 0} (so the region is
 *      empty) -- try each strict atom as the target.
 *   2. Eliminate LINEAR equalities by substitution (decline nonlinear ones).
 *   3. Putinar template  q == s0 + sum_i s_i g_i,  each s SOS = basis^T Q basis,
 *      Q >= 0.  Matching coefficients is a linear map from the symmetric Gram
 *      entries to the polynomial's coefficients: an affine system  A vec = b.
 *   4. Find PSD Gram matrices satisfying A vec = b NUMERICALLY, by alternating
 *      projection onto the PSD cone (LAPACK dsyev) and the affine subspace.
 *   5. If a boundary/interior zero x* is detected (near-kernel of the numeric
 *      Gram), recognise it as a rational point, verify q(x*)=0 exactly, and do
 *      EXACT facial reduction: shift every basis monomial by its value at x* so
 *      each basis polynomial vanishes at x* (the zero then lies in every Gram
 *      kernel by construction), which makes the reduced SDP strictly feasible
 *      and thus roundable.
 *   6. Round the numeric solution to an EXACT rational point of the affine
 *      subspace (round the free coordinates of A vec = b; the pivots follow
 *      exactly), rebuild each rational Gram matrix, and test it PSD by an exact
 *      rational LDL^T.  All PSD  => certificate verified  => the region is empty.
 *
 * Needs LAPACK for step 4; without it the stage declines.
 */

#include "reduce_sos.h"

#include "sym_names.h"
#include "eval.h"
#include "groebner.h"
#include "expr.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <gmp.h>

/* The whole certifier needs the LAPACK numeric guide; without it reduce_sos is
 * a stub that declines (the #else at the bottom of the file). */
#ifdef USE_LAPACK
#include "lapack.h"

/* Opt-in stderr tracing (never Head::tag form, so check-messages ignores it). */
static int sos_dbg(void) {
    static int v = -1;
    if (v < 0) { const char* e = getenv("MATHILDA_SOS_DEBUG"); v = (e && *e) ? 1 : 0; }
    return v;
}
#define DBG(...) do { if (sos_dbg()) fprintf(stderr, __VA_ARGS__); } while (0)

/* ================================================================== *
 *  Exponent-vector polynomial over mpq (a thin wrapper on GBPoly)     *
 * ================================================================== */

/* Full product of two GBPolys (GBPoly ships only mul_by_monomial + add). */
static GBPoly* sos_poly_mul(const GBPoly* a, const GBPoly* b) {
    if (!a || !b) return NULL;
    GBPoly* acc = gb_poly_new(a->n_vars, a->order, a->elim_pivot);
    if (!acc) return NULL;
    for (size_t t = 0; t < b->n_terms; t++) {
        GBPoly* term = gb_poly_mul_by_monomial(a, &b->exps[t * (size_t)b->n_vars], b->coefs[t]);
        if (!term) { gb_poly_free(acc); return NULL; }
        GBPoly* sum = gb_poly_add(acc, term);
        gb_poly_free(term);
        gb_poly_free(acc);
        if (!sum) return NULL;
        acc = sum;
    }
    return acc;
}

/* A single-term GBPoly: coefficient `c` times the monomial `exp` (length nv).
 * (gb_poly_new is the ZERO polynomial, and mul_by_monomial of zero stays zero,
 * so a monomial has to be built term-first.) */
static GBPoly* sos_monomial(int nv, const int* exp, const mpq_t c) {
    GBPoly* p = gb_poly_new(nv, GB_ORDER_GREVLEX, 0);
    if (!p) return NULL;
    if (mpq_sgn(c) == 0) return p;            /* the zero polynomial */
    p->exps  = malloc((size_t)nv * sizeof(int));
    p->coefs = malloc(sizeof(mpq_t));
    memcpy(p->exps, exp, (size_t)nv * sizeof(int));
    mpq_init(p->coefs[0]); mpq_set(p->coefs[0], c);
    p->n_terms = 1; p->cap = 1;               /* gb_poly_free clears `cap` coefs */
    return p;
}

/* Coefficient of the monomial `exp` (length n_vars) in p, into out (0 if absent). */
static void sos_coeff_of(const GBPoly* p, const int* exp, mpq_t out) {
    mpq_set_ui(out, 0, 1);
    for (size_t t = 0; t < p->n_terms; t++) {
        if (memcmp(&p->exps[t * (size_t)p->n_vars], exp, sizeof(int) * (size_t)p->n_vars) == 0) {
            mpq_set(out, p->coefs[t]);
            return;
        }
    }
}

/* Evaluate p at the rational point pt[n_vars], into out. */
static void sos_eval(const GBPoly* p, mpq_t* pt, mpq_t out) {
    mpq_set_ui(out, 0, 1);
    mpq_t term, pw;
    mpq_init(term); mpq_init(pw);
    for (size_t t = 0; t < p->n_terms; t++) {
        mpq_set(term, p->coefs[t]);
        for (int v = 0; v < p->n_vars; v++) {
            int e = p->exps[t * (size_t)p->n_vars + v];
            for (int k = 0; k < e; k++) mpq_mul(term, term, pt[v]);
        }
        mpq_add(out, out, term);
    }
    mpq_clear(term); mpq_clear(pw);
    mpq_canonicalize(out);
}

/* Evaluate p at a double point pt[n_vars]. */
static double sos_eval_d(const GBPoly* p, const double* pt) {
    double acc = 0;
    for (size_t t = 0; t < p->n_terms; t++) {
        double term = mpq_get_d(p->coefs[t]);
        for (int v = 0; v < p->n_vars; v++) {
            int e = p->exps[t * (size_t)p->n_vars + v];
            for (int k = 0; k < e; k++) term *= pt[v];
        }
        acc += term;
    }
    return acc;
}

/* Partial derivative d p / d x_var as a new GBPoly. */
static GBPoly* gb_partial(const GBPoly* p, int var) {
    GBPoly* d = gb_poly_new(p->n_vars, p->order, p->elim_pivot);
    for (size_t t = 0; t < p->n_terms; t++) {
        int e = p->exps[t * (size_t)p->n_vars + var];
        if (e <= 0) continue;
        int* ex = malloc((size_t)p->n_vars * sizeof(int));
        memcpy(ex, &p->exps[t * (size_t)p->n_vars], (size_t)p->n_vars * sizeof(int));
        ex[var] = e - 1;
        mpq_t c; mpq_init(c); mpq_set_si(c, e, 1); mpq_mul(c, c, p->coefs[t]);
        GBPoly* term = sos_monomial(p->n_vars, ex, c);
        GBPoly* s = gb_poly_add(d, term);
        gb_poly_free(d); gb_poly_free(term); free(ex); mpq_clear(c);
        d = s;
    }
    return d;
}

static int gb_total_degree(const GBPoly* p) {
    int d = 0;
    for (size_t t = 0; t < p->n_terms; t++) {
        int s = 0;
        for (int v = 0; v < p->n_vars; v++) s += p->exps[t * (size_t)p->n_vars + v];
        if (s > d) d = s;
    }
    return d;
}

/* ================================================================== *
 *  Monomial enumeration                                               *
 * ================================================================== */

/* All exponent vectors in `nv` variables with total degree <= dmax, as a flat
 * int array (count * nv), count returned in *cnt.  Caller frees. */
static int* enum_monomials(int nv, int dmax, int* cnt) {
    /* count first */
    int cap = 16, n = 0;
    int* out = malloc((size_t)cap * (size_t)nv * sizeof(int));
    int* e = calloc((size_t)nv, sizeof(int));
    for (;;) {
        int s = 0;
        for (int i = 0; i < nv; i++) s += e[i];
        if (s <= dmax) {
            if (n == cap) { cap *= 2; out = realloc(out, (size_t)cap * (size_t)nv * sizeof(int)); }
            memcpy(&out[n * nv], e, sizeof(int) * (size_t)nv);
            n++;
        }
        /* increment odometer e[] over 0..dmax each (prune by total degree) */
        int i = 0;
        for (; i < nv; i++) {
            e[i]++;
            if (e[i] <= dmax) break;
            e[i] = 0;
        }
        if (i == nv) break;
    }
    free(e);
    *cnt = n;
    return out;
}

/* ================================================================== *
 *  Exact rational symmetric-matrix PSD test (LDL^T, kernel-aware)     *
 * ================================================================== */

/* True iff the symmetric n*n rational matrix M (row-major, mpq) is positive
 * semidefinite.  Symmetric LDL with symmetric pivoting on a zero diagonal: a
 * zero pivot is admissible only if its whole (reduced) row/column is zero,
 * otherwise the matrix is indefinite. */
static bool mpq_matrix_psd(mpq_t* M, int n) {
    bool psd = true;
    mpq_t* A = malloc((size_t)n * (size_t)n * sizeof(mpq_t));
    for (int i = 0; i < n * n; i++) { mpq_init(A[i]); mpq_set(A[i], M[i]); }
    mpq_t factor, tmp;
    mpq_init(factor); mpq_init(tmp);
    int* perm = malloc((size_t)n * sizeof(int));
    for (int i = 0; i < n; i++) perm[i] = i;
    for (int k = 0; k < n && psd; k++) {
        /* find a nonzero diagonal pivot at k..n-1 (symmetric permutation) */
        int piv = -1;
        for (int i = k; i < n; i++) if (mpq_sgn(A[i * n + i]) != 0) { piv = i; break; }
        if (piv < 0) {
            /* all remaining diagonals zero: PSD only if the whole trailing
             * block is zero (a nonzero off-diagonal => indefinite). */
            for (int i = k; i < n && psd; i++)
                for (int j = k; j < n; j++)
                    if (mpq_sgn(A[i * n + j]) != 0) { psd = false; break; }
            break;
        }
        if (piv != k) {
            for (int j = 0; j < n; j++) { mpq_swap(A[k * n + j], A[piv * n + j]); }
            for (int i = 0; i < n; i++) { mpq_swap(A[i * n + k], A[i * n + piv]); }
        }
        if (mpq_sgn(A[k * n + k]) < 0) { psd = false; break; }   /* negative pivot */
        /* eliminate column k below */
        for (int i = k + 1; i < n; i++) {
            if (mpq_sgn(A[i * n + k]) == 0) continue;
            mpq_div(factor, A[i * n + k], A[k * n + k]);
            for (int j = k; j < n; j++) {
                mpq_mul(tmp, factor, A[k * n + j]);
                mpq_sub(A[i * n + j], A[i * n + j], tmp);
            }
        }
    }
    for (int i = 0; i < n * n; i++) mpq_clear(A[i]);
    free(A);
    mpq_clear(factor); mpq_clear(tmp); free(perm);
    return psd;
}

/* ================================================================== *
 *  double -> nearby rational (continued fractions, bounded denom)     *
 * ================================================================== */
static void double_to_mpq(double x, mpq_t out) {
    if (!isfinite(x)) { mpq_set_ui(out, 0, 1); return; }
    int neg = (x < 0); double ax = neg ? -x : x;
    /* continued fraction expansion with denominator cap */
    long p0 = 0, q0 = 1, p1 = 1, q1 = 0;
    double v = ax;
    for (int it = 0; it < 40; it++) {
        long a = (long)floor(v);
        long p2 = a * p1 + p0, q2 = a * q1 + q0;
        if (q2 > 1000000L) break;
        p0 = p1; q0 = q1; p1 = p2; q1 = q2;
        double frac = v - (double)a;
        if (frac < 1e-12) break;
        v = 1.0 / frac;
    }
    if (q1 == 0) { mpq_set_si(out, neg ? -(long)llround(ax) : (long)llround(ax), 1); return; }
    mpq_set_si(out, neg ? -p1 : p1, (unsigned long)q1);
    mpq_canonicalize(out);
}

/* ================================================================== *
 *  SOS block: a Gram block over a basis of rational "basis polys"     *
 * ================================================================== */
typedef struct {
    int         nb;          /* basis size */
    GBPoly**    basis;       /* nb basis polynomials (each owns) */
    const GBPoly* mult;      /* borrowed multiplier g_i (NULL => constant 1) */
    /* precomputed products P[a*nb+b] = basis[a]*basis[b]*mult (a<=b) */
    GBPoly**    prod;        /* nb*nb, only a<=b filled */
} SosBlock;

static void sosblock_free(SosBlock* B) {
    if (!B) return;
    if (B->basis) { for (int i = 0; i < B->nb; i++) gb_poly_free(B->basis[i]); free(B->basis); }
    if (B->prod)  { for (int i = 0; i < B->nb * B->nb; i++) if (B->prod[i]) gb_poly_free(B->prod[i]); free(B->prod); }
    free(B);
}

/* ================================================================== *
 *  The certifier: prove q >= 0 on { g_e[i] >= 0 }  (no equalities)    *
 * ================================================================== */
/* Minimise q on K = {g_i >= 0} numerically (quadratic-penalty gradient descent,
 * multistart), round the best point to rationals, and verify q(x*)==0 and every
 * g_i(x*)>=0 in EXACT arithmetic.  Fills zero[nrv] (caller inits/clears) and
 * returns true on a verified rational zero in K; else false.  Only rational
 * zeros are usable for the exact facial reduction, which is exactly the
 * interior-zero class these certificates need. */
static bool sos_find_zero(const GBPoly* q, GBPoly** g, int ng, int nrv, mpq_t* zero) {
    double qscale = 1.0;
    for (size_t t = 0; t < q->n_terms; t++) { double v = fabs(mpq_get_d(q->coefs[t])); if (v > qscale) qscale = v; }
    GBPoly** dq = malloc((size_t)nrv * sizeof(GBPoly*));
    for (int i = 0; i < nrv; i++) dq[i] = gb_partial(q, i);
    GBPoly*** dg = malloc((size_t)ng * sizeof(GBPoly**));
    for (int k = 0; k < ng; k++) { dg[k] = malloc((size_t)nrv * sizeof(GBPoly*)); for (int i = 0; i < nrv; i++) dg[k][i] = gb_partial(g[k], i); }

    double* best = malloc((size_t)nrv * sizeof(double));
    double bestq = 1e300; bool have = false;
    double* x = malloc((size_t)nrv * sizeof(double));
    double* grad = malloc((size_t)nrv * sizeof(double));
    unsigned seed = 12345;
    for (int start = 0; start < 60; start++) {
        for (int i = 0; i < nrv; i++) {
            if (start == 0) x[i] = 1.0 / (nrv + 1);              /* simplex centroid guess */
            else { seed = seed * 1103515245u + 12345u; x[i] = ((seed >> 8) & 0xffff) / 65535.0 * 1.5 - 0.25; }
        }
        double mu = 4.0;
        for (int phase = 0; phase < 3; phase++, mu *= 20.0) {
            double lr = 0.05;
            for (int step = 0; step < 600; step++) {
                for (int i = 0; i < nrv; i++) grad[i] = sos_eval_d(dq[i], x) / qscale;
                for (int k = 0; k < ng; k++) {
                    double gv = sos_eval_d(g[k], x);
                    if (gv < 0) for (int i = 0; i < nrv; i++) grad[i] += mu * 2.0 * (-gv) * (-sos_eval_d(dg[k][i], x)) / qscale;
                }
                double gn = 0; for (int i = 0; i < nrv; i++) gn += grad[i] * grad[i];
                gn = sqrt(gn); if (gn < 1e-14) break;
                for (int i = 0; i < nrv; i++) x[i] -= lr * grad[i] / (gn > 1 ? gn : 1);
                if (step % 200 == 199) lr *= 0.5;
            }
        }
        /* feasible? */
        bool feas = true;
        for (int k = 0; k < ng; k++) if (sos_eval_d(g[k], x) < -1e-6) feas = false;
        double qv = fabs(sos_eval_d(q, x)) / qscale;
        if (feas && qv < bestq) { bestq = qv; for (int i = 0; i < nrv; i++) best[i] = x[i]; have = true; }
    }

    bool ok = false;
    if (have && bestq < 1e-4) {
        for (int i = 0; i < nrv; i++) double_to_mpq(best[i], zero[i]);
        mpq_t val; mpq_init(val);
        sos_eval(q, zero, val);
        ok = (mpq_sgn(val) == 0);
        for (int k = 0; k < ng && ok; k++) { sos_eval(g[k], zero, val); if (mpq_sgn(val) < 0) ok = false; }
        mpq_clear(val);
        if (sos_dbg()) {
            DBG("[sos] find_zero: bestq=%.3e verified=%d  x*=", bestq, ok);
            for (int i = 0; i < nrv; i++) { char* s = mpq_get_str(NULL, 10, zero[i]); DBG("%s ", s); free(s); }
            DBG("\n");
        }
    }
    for (int i = 0; i < nrv; i++) gb_poly_free(dq[i]);
    free(dq);
    for (int k = 0; k < ng; k++) { for (int i = 0; i < nrv; i++) gb_poly_free(dg[k][i]); free(dg[k]); }
    free(dg);
    free(best); free(x); free(grad);
    return ok;
}

static bool sos_certify(Expr* q_e, Expr** g_e, int ng, Expr** rvars, int nrv) {
    bool proven = false;
    int D = 0;
    GBPoly* q = gb_from_expr(q_e, rvars, nrv, GB_ORDER_GREVLEX, 0, NULL);
    GBPoly** g = calloc((size_t)ng, sizeof(GBPoly*));
    if (!q || !g) { gb_poly_free(q); free(g); return false; }
    for (int i = 0; i < ng; i++) {
        g[i] = gb_from_expr(g_e[i], rvars, nrv, GB_ORDER_GREVLEX, 0, NULL);
        if (!g[i]) goto cleanup;
    }
    D = gb_total_degree(q);
    if (D < 2) D = 2;
    if (D % 2) D++;
    if (D > 8) goto cleanup;   /* degree ladder cap (bounded; declines above) */

    /* The coefficient coordinates: all monomials up to degree D. */
    ;
    {
    int ncoef = 0;
    int* coef_exp = enum_monomials(nrv, D, &ncoef);
    /* index lookup for a monomial exponent */
    /* (linear scan is fine at these sizes) */

    /* Build blocks: block 0 = s0 (mult 1, basis deg<=D/2); block i = s_i g_i
     * (mult g[i-1], basis deg <= floor((D-deg g_i)/2)). */
    int nblocks = 1 + ng;
    SosBlock** blk = calloc((size_t)nblocks, sizeof(SosBlock*));

    for (int bi = 0; bi < nblocks; bi++) {
        const GBPoly* mult = (bi == 0) ? NULL : g[bi - 1];
        int dmult = (bi == 0) ? 0 : gb_total_degree(g[bi - 1]);
        int dbasis = (D - dmult) / 2;
        if (dbasis < 0) dbasis = 0;
        int nmono = 0;
        int* mex = enum_monomials(nrv, dbasis, &nmono);
        SosBlock* B = calloc(1, sizeof(SosBlock));
        B->nb = nmono; B->mult = mult;
        B->basis = calloc((size_t)nmono, sizeof(GBPoly*));
        mpq_t one; mpq_init(one); mpq_set_ui(one, 1, 1);
        for (int a = 0; a < nmono; a++)
            B->basis[a] = sos_monomial(nrv, &mex[a * nrv], one);
        mpq_clear(one);
        free(mex);
        blk[bi] = B;
    }

    /* Facial reduction: if q has a rational zero x* on K then every Gram matrix
     * must annihilate the basis evaluated at x* (an interior/boundary zero makes
     * the raw SDP rank-deficient and so unroundable).  Find x*, then shift every
     * basis monomial by its value at x* so each basis polynomial vanishes at x*;
     * the reduced SDP is then strictly feasible on its face and roundable. */
    mpq_t* zero_pt = malloc((size_t)nrv * sizeof(mpq_t));
    for (int i = 0; i < nrv; i++) mpq_init(zero_pt[i]);
    bool have_zero = sos_find_zero(q, g, ng, nrv, zero_pt);
    if (have_zero) {
        for (int bi = 0; bi < nblocks; bi++) {
            SosBlock* B = blk[bi];
            GBPoly** nb2 = calloc((size_t)B->nb, sizeof(GBPoly*));
            int cnt = 0;
            for (int a = 0; a < B->nb; a++) {
                if (B->basis[a]->n_terms == 0) continue;
                int tot = 0; for (int vv = 0; vv < nrv; vv++) tot += B->basis[a]->exps[vv];
                if (tot == 0) continue;                   /* drop the constant monomial */
                mpq_t valx; mpq_init(valx);
                sos_eval(B->basis[a], zero_pt, valx);
                GBPoly* shifted = gb_poly_copy(B->basis[a]);
                if (mpq_sgn(valx) != 0) {
                    int* z = calloc((size_t)nrv, sizeof(int));
                    GBPoly* cstp = sos_monomial(nrv, z, valx);  /* valx * 1 */
                    free(z);
                    GBPoly* s = gb_poly_sub(shifted, cstp);
                    gb_poly_free(shifted); gb_poly_free(cstp);
                    shifted = s;
                }
                mpq_clear(valx);
                nb2[cnt++] = shifted;
            }
            for (int a = 0; a < B->nb; a++) gb_poly_free(B->basis[a]);
            free(B->basis); B->basis = nb2; B->nb = cnt;
        }
    }

    /* Precompute block products P[a,b] = basis_a * basis_b * mult. */
    for (int bi = 0; bi < nblocks; bi++) {
        SosBlock* B = blk[bi];
        B->prod = calloc((size_t)B->nb * (size_t)B->nb, sizeof(GBPoly*));
        for (int a = 0; a < B->nb; a++)
            for (int b = a; b < B->nb; b++) {
                GBPoly* ab = sos_poly_mul(B->basis[a], B->basis[b]);
                GBPoly* abm = B->mult ? sos_poly_mul(ab, B->mult) : gb_poly_copy(ab);
                gb_poly_free(ab);
                B->prod[a * B->nb + b] = abm;
            }
    }

    /* Enumerate Gram variables: (block, a<=b) -> var index. */
    int V = 0;
    for (int bi = 0; bi < nblocks; bi++) V += blk[bi]->nb * (blk[bi]->nb + 1) / 2;

    /* Affine system A (ncoef x V) and rhs b: coeff of each monomial in
     * s0 + sum s_i g_i must equal coeff in q. */
    double* A = calloc((size_t)ncoef * (size_t)V, sizeof(double));
    double* rhs = calloc((size_t)ncoef, sizeof(double));
    mpq_t cq; mpq_init(cq);
    for (int c = 0; c < ncoef; c++) {
        sos_coeff_of(q, &coef_exp[c * nrv], cq);
        rhs[c] = mpq_get_d(cq);
    }
    {
        int var = 0;
        mpq_t cf; mpq_init(cf);
        for (int bi = 0; bi < nblocks; bi++) {
            SosBlock* B = blk[bi];
            for (int a = 0; a < B->nb; a++)
                for (int b = a; b < B->nb; b++) {
                    GBPoly* P = B->prod[a * B->nb + b];
                    double mult2 = (a == b) ? 1.0 : 2.0;
                    for (int c = 0; c < ncoef; c++) {
                        sos_coeff_of(P, &coef_exp[c * nrv], cf);
                        double val = mpq_get_d(cf);
                        if (val != 0.0) A[c * V + var] += mult2 * val;
                    }
                    var++;
                }
        }
        mpq_clear(cf);
    }
    mpq_clear(cq);

    if (sos_dbg()) {
        double bn = 0; for (int c = 0; c < ncoef; c++) bn += rhs[c]*rhs[c];
        double an = 0; for (int i = 0; i < ncoef*V; i++) an += A[i]*A[i];
        double qmax = 0; for (size_t t = 0; t < q->n_terms; t++) { double v = fabs(mpq_get_d(q->coefs[t])); if (v > qmax) qmax = v; }
        DBG("[sos] facial=%d: nrv=%d D=%d q_terms=%zu q_deg=%d maxcoef=%.3g ncoef=%d V=%d ||b||^2=%.3e ||A||^2=%.3e\n",
            have_zero, nrv, D, q->n_terms, gb_total_degree(q), qmax, ncoef, V, bn, an);
    }

    /* ---- numeric SDP by alternating projection (Dykstra-free POCS) ---- */
    /* Variable vector x (length V) <-> per-block symmetric matrices. */
    double* x = calloc((size_t)V, sizeof(double));
    /* init at a small positive-definite guess: identity on each block */
    {
        int var = 0;
        for (int bi = 0; bi < nblocks; bi++) {
            SosBlock* B = blk[bi];
            for (int a = 0; a < B->nb; a++)
                for (int b = a; b < B->nb; b++) { x[var] = (a == b) ? 1.0 : 0.0; var++; }
        }
    }

    /* Precompute A A^T (ncoef x ncoef) for the affine projection solve. */
    double* AAt = calloc((size_t)ncoef * (size_t)ncoef, sizeof(double));
    for (int i = 0; i < ncoef; i++)
        for (int j = 0; j < ncoef; j++) {
            double s = 0; for (int k = 0; k < V; k++) s += A[i * V + k] * A[j * V + k];
            AAt[i * ncoef + j] = s;
        }

    bool near_feasible = false;
    /* Tikhonov regularisation for the (A A^T) solve -- A can be rank-deficient
     * (dependent coefficient rows), which makes the raw normal-equation solve
     * blow up; a relative ridge bounds the condition number so the projection is
     * a stable (lightly damped) contraction and POCS stays convergent. */
    double aat_tr = 0; for (int i = 0; i < ncoef; i++) aat_tr += AAt[i * ncoef + i];
    double reg = 1e-7 * (aat_tr / (ncoef ? ncoef : 1)) + 1e-12;
    /* Strict-interior eigenvalue floor for the PSD projection (see the clamp
     * below).  Scaled to the data so it is meaningful against b's magnitude. */
    double bmax = 0; for (int c = 0; c < ncoef; c++) { double v = fabs(rhs[c]); if (v > bmax) bmax = v; }
    double ef = 1e-4 * (bmax > 1 ? bmax : 1);
    double* work = calloc((size_t)ncoef, sizeof(double));
    int* ipiv = calloc((size_t)ncoef, sizeof(int));
    double* AAt_f = malloc((size_t)ncoef * (size_t)ncoef * sizeof(double));
    double best_res = 1e300;
    for (int iter = 0; iter < 1500; iter++) {
        /* affine projection: x <- x - A^T (A A^T + reg I)^{-1}(A x - b) */
        for (int i = 0; i < ncoef; i++) {
            double s = 0; for (int k = 0; k < V; k++) s += A[i * V + k] * x[k];
            work[i] = s - rhs[i];
        }
        /* solve (A A^T + reg I) y = work  (column-major copy for LAPACK) */
        for (int i = 0; i < ncoef; i++)
            for (int j = 0; j < ncoef; j++)
                AAt_f[j * ncoef + i] = AAt[i * ncoef + j] + (i == j ? reg : 0.0);
        int info = mat_lapack_dgesv(ncoef, 1, AAt_f, ncoef, ipiv, work, ncoef);
        if (info != 0) continue;   /* skip this projection step; ridge makes it rare */
        for (int k = 0; k < V; k++) {
            double s = 0; for (int i = 0; i < ncoef; i++) s += A[i * V + k] * work[i];
            x[k] -= s;
        }
        /* PSD projection per block */
        double maxneg = 0.0;
        int var = 0;
        for (int bi = 0; bi < nblocks; bi++) {
            SosBlock* B = blk[bi];
            int n = B->nb;
            if (n == 0) continue;
            double* Q = malloc((size_t)n * (size_t)n * sizeof(double));
            int v2 = var;
            for (int a = 0; a < n; a++)
                for (int b = a; b < n; b++) { double val = x[v2++]; Q[a * n + b] = val; Q[b * n + a] = val; }
            double* w = malloc((size_t)n * sizeof(double));
            double* Qc = malloc((size_t)n * (size_t)n * sizeof(double));
            for (int i = 0; i < n * n; i++) Qc[i] = Q[i];
            int einfo = mat_lapack_dsyev(n, Qc, n, w);   /* Qc columns = eigenvectors */
            if (einfo == 0) {
                if (w[0] < maxneg) maxneg = w[0];
                /* Project onto { Q >= ef I }: clamp every eigenvalue UP to the
                 * positive floor ef.  This steers POCS to a STRICTLY interior
                 * (full-rank, lambda_min >= ef) feasible point when one exists on
                 * the facial-reduced face -- the boundary point plain clamping to
                 * 0 lands on is rank-deficient and so unroundable. */
                for (int i = 0; i < n * n; i++) Q[i] = 0.0;
                for (int e = 0; e < n; e++) {
                    double lam = w[e] > ef ? w[e] : ef;
                    if (lam == 0.0) continue;
                    for (int a = 0; a < n; a++)
                        for (int b = 0; b < n; b++)
                            Q[a * n + b] += lam * Qc[a + e * n] * Qc[b + e * n];  /* col-major evec */
                }
                v2 = var;
                for (int a = 0; a < n; a++)
                    for (int b = a; b < n; b++) x[v2++] = Q[a * n + b];
            }
            free(Q); free(w); free(Qc);
            var += n * (n + 1) / 2;
        }
        /* convergence: affine residual after the PSD step */
        double res = 0; for (int i = 0; i < ncoef; i++) { double s = 0; for (int k = 0; k < V; k++) s += A[i*V+k]*x[k]; s -= rhs[i]; res += s*s; }
        if (res < best_res) best_res = res;
        if (res < 1e-9 && maxneg > -1e-5) near_feasible = true;
        if (sos_dbg() && (iter % 300 == 0 || iter == 1499))
            DBG("[sos] iter %d affine-res=%.3e maxneg=%.3e\n", iter, res, maxneg);
    }
    free(work); free(ipiv); free(AAt_f); free(AAt);

    /* ---- exact rounding + verification ----
     * Always attempt: the exact LDL^T + consistency check is the sound arbiter.
     * A diverged/garbage numeric solution simply fails to verify (declines). */
    (void)near_feasible;
    {
        DBG("[sos] exact-verifying (have_zero=%d)\n", have_zero);
        /* Build exact affine system A_q vec = b_q (rational) and round x onto it.
         * Parametrise by rref: pick pivot (basic) vars, free vars take rounded x. */
        /* Exact A_q and b_q */
        mpq_t* Aq = malloc((size_t)ncoef * (size_t)V * sizeof(mpq_t));
        mpq_t* bq = malloc((size_t)ncoef * sizeof(mpq_t));
        for (int i = 0; i < ncoef * V; i++) mpq_init(Aq[i]);
        for (int i = 0; i < ncoef; i++) mpq_init(bq[i]);
        {
            mpq_t cf, cq2; mpq_init(cf); mpq_init(cq2);
            for (int c = 0; c < ncoef; c++) { sos_coeff_of(q, &coef_exp[c*nrv], cq2); mpq_set(bq[c], cq2); }
            int var = 0;
            for (int bi = 0; bi < nblocks; bi++) {
                SosBlock* B = blk[bi];
                for (int a = 0; a < B->nb; a++)
                    for (int b = a; b < B->nb; b++) {
                        GBPoly* P = B->prod[a * B->nb + b];
                        for (int c = 0; c < ncoef; c++) {
                            sos_coeff_of(P, &coef_exp[c*nrv], cf);
                            if (mpq_sgn(cf) != 0) {
                                if (a != b) { mpq_t t; mpq_init(t); mpq_set_ui(t,2,1); mpq_mul(t,t,cf); mpq_add(Aq[c*V+var], Aq[c*V+var], t); mpq_clear(t); }
                                else mpq_add(Aq[c*V+var], Aq[c*V+var], cf);
                            }
                        }
                        var++;
                    }
            }
            mpq_clear(cf); mpq_clear(cq2);
        }

        /* Gaussian elimination (rref) on [Aq | bq] to identify pivot columns. */
        int* pivcol = malloc((size_t)ncoef * sizeof(int));
        int rank = 0;
        {
            mpq_t f, t; mpq_init(f); mpq_init(t);
            int row = 0;
            char* is_pivot_col = calloc((size_t)V, 1);
            for (int col = 0; col < V && row < ncoef; col++) {
                int sel = -1;
                for (int r = row; r < ncoef; r++) if (mpq_sgn(Aq[r*V+col]) != 0) { sel = r; break; }
                if (sel < 0) continue;
                if (sel != row) {
                    for (int cc = 0; cc < V; cc++) mpq_swap(Aq[row*V+cc], Aq[sel*V+cc]);
                    mpq_swap(bq[row], bq[sel]);
                }
                /* normalise row */
                mpq_set(f, Aq[row*V+col]);
                for (int cc = 0; cc < V; cc++) mpq_div(Aq[row*V+cc], Aq[row*V+cc], f);
                mpq_div(bq[row], bq[row], f);
                for (int r = 0; r < ncoef; r++) {
                    if (r == row) continue;
                    if (mpq_sgn(Aq[r*V+col]) == 0) continue;
                    mpq_set(f, Aq[r*V+col]);
                    for (int cc = 0; cc < V; cc++) { mpq_mul(t, f, Aq[row*V+cc]); mpq_sub(Aq[r*V+cc], Aq[r*V+cc], t); }
                    mpq_mul(t, f, bq[row]); mpq_sub(bq[r], bq[r], t);
                }
                pivcol[row] = col; is_pivot_col[col] = 1; row++;
            }
            rank = row;
            /* consistency: a zero row with nonzero rhs => the numeric guide was
             * infeasible exactly (shouldn't happen if near_feasible, but guard). */
            bool inconsistent = false;
            for (int r = rank; r < ncoef; r++) if (mpq_sgn(bq[r]) != 0) inconsistent = true;
            if (!inconsistent) {
                /* exact solution vec: free vars = rounded x[col]; basic solved. */
                mpq_t* vec = malloc((size_t)V * sizeof(mpq_t));
                for (int i = 0; i < V; i++) { mpq_init(vec[i]); }
                for (int col = 0; col < V; col++)
                    if (!is_pivot_col[col]) double_to_mpq(x[col], vec[col]);
                for (int r = rank - 1; r >= 0; r--) {
                    int pc = pivcol[r];
                    mpq_t acc; mpq_init(acc); mpq_set(acc, bq[r]);
                    for (int col = 0; col < V; col++) {
                        if (col == pc) continue;
                        if (mpq_sgn(Aq[r*V+col]) == 0) continue;
                        mpq_mul(t, Aq[r*V+col], vec[col]); mpq_sub(acc, acc, t);
                    }
                    mpq_set(vec[pc], acc); mpq_clear(acc);
                }
                /* rebuild each block's rational Gram and test PSD */
                bool all_psd = true;
                int var = 0;
                for (int bi = 0; bi < nblocks && all_psd; bi++) {
                    int n = blk[bi]->nb;
                    if (n == 0) continue;
                    mpq_t* Q = malloc((size_t)n * (size_t)n * sizeof(mpq_t));
                    for (int i = 0; i < n*n; i++) mpq_init(Q[i]);
                    int v2 = var;
                    for (int a = 0; a < n; a++)
                        for (int b = a; b < n; b++) { mpq_set(Q[a*n+b], vec[v2]); mpq_set(Q[b*n+a], vec[v2]); v2++; }
                    if (!mpq_matrix_psd(Q, n)) all_psd = false;
                    for (int i = 0; i < n*n; i++) mpq_clear(Q[i]);
                    free(Q);
                    var += n*(n+1)/2;
                }
                DBG("[sos] exact verify: rank=%d all_psd=%d\n", rank, all_psd);
                if (all_psd) proven = true;
                for (int i = 0; i < V; i++) mpq_clear(vec[i]);
                free(vec);
            } else DBG("[sos] exact system inconsistent (guide infeasible)\n");
            free(is_pivot_col);
            mpq_clear(f); mpq_clear(t);
        }
        for (int i = 0; i < ncoef * V; i++) mpq_clear(Aq[i]);
        for (int i = 0; i < ncoef; i++) mpq_clear(bq[i]);
        free(Aq); free(bq); free(pivcol);
    }

    free(A); free(rhs); free(x); free(coef_exp);
    for (int bi = 0; bi < nblocks; bi++) sosblock_free(blk[bi]);
    free(blk);
    if (zero_pt) { for (int vv = 0; vv < nrv; vv++) mpq_clear(zero_pt[vv]); free(zero_pt); }
    }

cleanup:
    gb_poly_free(q);
    for (int i = 0; i < ng; i++) if (g[i]) gb_poly_free(g[i]);
    free(g);
    return proven;
}

/* ================================================================== *
 *  Front end: extract the refutation shape from the RForm             *
 * ================================================================== */

static bool is_sym_name(const Expr* e, const char* n) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == n;
}

/* Is `p` a polynomial in rvars? (reuse the evaluator's PolynomialQ) */
static bool poly_in_vars(Expr* p, Expr** rvars, int nrv) {
    Expr** va = malloc((size_t)nrv * sizeof(Expr*));
    for (int i = 0; i < nrv; i++) va[i] = expr_copy(rvars[i]);
    Expr* vl = expr_new_function(expr_new_symbol(SYM_List), va, (size_t)nrv);
    free(va);
    Expr* call = expr_new_function(expr_new_symbol("PolynomialQ"),
                                   (Expr*[]){ expr_copy(p), vl }, 2);
    Expr* r = eval_and_free(call);
    bool ok = is_sym_name(r, SYM_True);
    expr_free(r);
    return ok;
}

Expr* reduce_sos(const RForm* F, Expr** vlist, int nv) {
    if (!F || F->is_true || F->n != 1) return NULL;   /* single conjunction only */
    RConj* cj = F->c[0];
    if (cj->is_false) return NULL;

    /* classify atoms */
    int ng = 0, nh = 0, nstrict = 0;
    for (int k = 0; k < cj->n; k++) {
        RAtom* a = &cj->a[k];
        if (a->rel == R_ELEM || a->nonconst_denom) return NULL;   /* out of scope */
        if (a->rel == R_EQ) nh++;
        else if (a->rel == R_LE) ng++;
        else if (a->rel == R_LT) { ng++; nstrict++; }
        /* R_NE: dropped (sound; enlarges the region) */
    }
    if (nstrict == 0) return NULL;   /* no strict atom to refute */

    /* Only LINEAR equalities handled (eliminate by substitution). */
    /* Gather equality polys, inequality polys g = -poly (>=0), strict targets. */
    Expr** heq = nh ? calloc((size_t)nh, sizeof(Expr*)) : NULL;
    int nhe = 0;
    for (int k = 0; k < cj->n; k++)
        if (cj->a[k].rel == R_EQ) heq[nhe++] = cj->a[k].poly;   /* borrowed */

    Expr* result = NULL;

    /* try each strict atom as the refutation target */
    for (int tk = 0; tk < cj->n && !result; tk++) {
        if (cj->a[tk].rel != R_LT) continue;
        /* q = poly of the strict target (prove q >= 0) */
        /* g list = -(poly) for every R_LE and every OTHER R_LT */
        Expr* q = expr_copy(cj->a[tk].poly);
        int gc = 0;
        Expr** gl = calloc((size_t)(ng), sizeof(Expr*));
        for (int k = 0; k < cj->n; k++) {
            if (k == tk) continue;
            RAtom* a = &cj->a[k];
            if (a->rel == R_LE || a->rel == R_LT) {
                /* g = -poly >= 0 */
                gl[gc++] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                             (Expr*[]){ expr_new_integer(-1), expr_copy(a->poly) }, 2));
            }
        }

        /* reduced variable set + linear-equality elimination */
        Expr** rvars = malloc((size_t)nv * sizeof(Expr*));
        int nrv = 0;
        for (int i = 0; i < nv; i++) rvars[nrv++] = expr_copy(vlist[i]);

        bool ok = true;
        for (int j = 0; j < nhe && ok; j++) {
            /* solve heq[j] == 0 for some remaining variable with constant coeff */
            int chosen = -1; Expr* rhs = NULL;
            for (int vi = 0; vi < nrv && chosen < 0; vi++) {
                /* coeff of rvars[vi]^1 in heq[j] */
                Expr* co = eval_and_free(expr_new_function(expr_new_symbol(SYM_Coefficient),
                             (Expr*[]){ expr_copy(heq[j]), expr_copy(rvars[vi]), expr_new_integer(1) }, 3));
                /* constant (free of all rvars)? and nonzero */
                bool is_const = true;
                for (int w = 0; w < nrv; w++) if (!poly_in_vars(co, &rvars[w], 1) ? false : false) {}
                /* check co does not contain any rvar */
                for (int w = 0; w < nrv; w++) {
                    Expr* d = eval_and_free(expr_new_function(expr_new_symbol("Coefficient"),
                                (Expr*[]){ expr_copy(co), expr_copy(rvars[w]), expr_new_integer(1) }, 3));
                    bool zero = (d->type == EXPR_INTEGER && d->data.integer == 0);
                    expr_free(d);
                    if (!zero) { is_const = false; break; }
                }
                bool nonzero = !(co->type == EXPR_INTEGER && co->data.integer == 0);
                if (is_const && nonzero) {
                    /* heq = co*v + rest ; v = -rest/co ; rest = heq - co*v */
                    Expr* rest = eval_and_free(expr_new_function(expr_new_symbol("Expand"),
                                   (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Plus),
                                       (Expr*[]){ expr_copy(heq[j]),
                                           expr_new_function(expr_new_symbol(SYM_Times),
                                               (Expr*[]){ expr_new_integer(-1), expr_copy(co), expr_copy(rvars[vi]) }, 3) }, 2) }, 1));
                    rhs = eval_and_free(expr_new_function(expr_new_symbol("Expand"),
                            (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Times),
                                (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Power),
                                    (Expr*[]){ expr_copy(co), expr_new_integer(-1) }, 2),
                                    expr_new_function(expr_new_symbol(SYM_Times),
                                        (Expr*[]){ expr_new_integer(-1), rest }, 2) }, 2) }, 1));
                    chosen = vi;
                }
                expr_free(co);
            }
            if (chosen < 0) { ok = false; break; }
            /* substitute rvars[chosen] -> rhs in q and all g */
            Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
                           (Expr*[]){ expr_copy(rvars[chosen]), expr_copy(rhs) }, 2);
            Expr* nq = eval_and_free(expr_new_function(expr_new_symbol("Expand"),
                         (Expr*[]){ expr_new_function(expr_new_symbol(SYM_ReplaceAll),
                             (Expr*[]){ q, expr_copy(rule) }, 2) }, 1));
            q = nq;
            for (int gi = 0; gi < gc; gi++)
                gl[gi] = eval_and_free(expr_new_function(expr_new_symbol("Expand"),
                           (Expr*[]){ expr_new_function(expr_new_symbol(SYM_ReplaceAll),
                               (Expr*[]){ gl[gi], expr_copy(rule) }, 2) }, 1));
            /* also substitute into remaining equalities */
            for (int j2 = j + 1; j2 < nhe; j2++) {
                /* heq[j2] is borrowed from the atom; we need an owned working copy.
                 * Keep it simple: copy-substitute into a temp and swap via a side array. */
            }
            expr_free(rule); expr_free(rhs);
            /* drop rvars[chosen] */
            expr_free(rvars[chosen]);
            for (int w = chosen; w < nrv - 1; w++) rvars[w] = rvars[w + 1];
            nrv--;
        }

        if (ok && nhe <= 1 && nrv >= 1) {
            /* verify q and all g are polynomials in the reduced vars */
            bool allpoly = poly_in_vars(q, rvars, nrv);
            for (int gi = 0; gi < gc && allpoly; gi++) allpoly = poly_in_vars(gl[gi], rvars, nrv);
            if (allpoly && sos_certify(q, gl, gc, rvars, nrv))
                result = expr_new_symbol(SYM_False);
        }

        expr_free(q);
        for (int gi = 0; gi < gc; gi++) expr_free(gl[gi]);
        free(gl);
        for (int i = 0; i < nrv; i++) expr_free(rvars[i]);
        free(rvars);
    }

    free(heq);
    return result;
}

#else /* !USE_LAPACK: the numeric guide is unavailable, so decline everything. */

Expr* reduce_sos(const RForm* F, Expr** vlist, int nv) {
    (void)F; (void)vlist; (void)nv;
    return NULL;
}

#endif /* USE_LAPACK */
