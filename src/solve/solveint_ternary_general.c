/*
 * solveint_ternary_general.c
 *
 * Part of the Solve[..., Integers] engine; split out of solveint.c.
 * See solveint_internal.h for the shared SICtx/SearchState substrate.
 *
 * The GENERAL homogeneous ternary quadratic
 *
 *     Q(x,y,z) = a x^2 + b y^2 + c z^2 + d x y + e y z + f z x == 0
 *
 * over an unbounded domain (cross terms and/or non-symmetric diagonal
 * coefficients -- everything the symmetric x^2+y^2==k z^2 solver in
 * solveint_ternary.c declines).  A non-degenerate ternary quadratic that has one
 * rational point is a genus-0 curve, hence rational: it has either only the
 * trivial point (definite / Legendre fails) or a complete 2-parameter integer
 * family obtained by the chord construction.
 *
 * Method:
 *   1. Extract the integer symmetric matrix M (2x convention: G(v)=v^T M v=2Q(v)).
 *   2. Congruently diagonalise M over Q (mpq), tracking the transform T with
 *      T^T M T = diag(D_0,D_1,D_2), so v = T w maps a diagonal solution back.
 *   3. Clear denominators -> an integer diagonal form  C_0 u^2+C_1 v^2+C_2 w^2=0,
 *      proportional to sum D_i w_i^2.  Reduce (content / squarefree / coprime) to
 *      decide solvability by LEGENDRE's theorem (real indefiniteness + the three
 *      quadratic-residue conditions).  Not solvable -> trivial-only {{0,0,0}}
 *      (Mathematica's answer, a proof of no nontrivial solution).
 *   4. Solvable -> find a witness of the integer diagonal form within Holzer's
 *      bound (guaranteed to exist), map it back through T, clear denominators and
 *      divide by content -> a primitive integer witness P0 with G(P0)==0.
 *   5. Chord-parametrise IN THE ORIGINAL COORDINATES from P0: the second
 *      intersection of the line P0 + (s e_i + t e_j) with the conic is
 *      G(V)*P0 - 2 (P0^T M V) V, a quadratic in (s,t)=(C[1],C[2]) with integer
 *      coefficients; scale by C[3]; plus the tangent family C[1]*P0.  Only ONE
 *      point (the witness) is ever mapped back, so the family never inherits the
 *      diagonalisation's denominators.
 *
 * Scope: non-degenerate (rank-3) forms; a rank-deficient form (a product of
 * linear factors) declines.  Homogeneous only (any linear / constant term
 * declines).  Coefficients small enough for the Holzer witness box and the
 * squarefree extraction; otherwise declines (never a wrong answer).
 */
#include "solveint.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "eval.h"
#include "expr.h"
#include "poly/mpoly.h"
#include "numbertheory/numbertheory_internal.h"
#include "solveint_internal.h"


/* Extract the homogeneous ternary quadratic coefficients over the three active
 * variables vs[0..2].  a,b,c are the x^2,y^2,z^2 coefficients; d,e,f the xy,yz,zx
 * cross coefficients.  Returns false on any linear / constant / higher term or a
 * term touching a variable outside vs (i.e. not a pure homogeneous degree-2 form
 * in exactly these three variables). */
static bool tg_extract(const MPoly* eq, int n, const int* vs,
                       mpz_t a, mpz_t b, mpz_t cc, mpz_t d, mpz_t e, mpz_t f) {
    mpz_set_ui(a, 0); mpz_set_ui(b, 0); mpz_set_ui(cc, 0);
    mpz_set_ui(d, 0); mpz_set_ui(e, 0); mpz_set_ui(f, 0);
    for (size_t t = 0; t < eq->n_terms; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        for (int v = 0; v < n; v++)
            if (v != vs[0] && v != vs[1] && v != vs[2] && ex[v] != 0) return false;
        int d0 = ex[vs[0]], d1 = ex[vs[1]], d2 = ex[vs[2]];
        if (d0 + d1 + d2 != 2) return false;            /* not homogeneous degree 2 */
        if      (d0 == 2)              mpz_set(a, eq->coefs[t]);
        else if (d1 == 2)              mpz_set(b, eq->coefs[t]);
        else if (d2 == 2)              mpz_set(cc, eq->coefs[t]);
        else if (d0 == 1 && d1 == 1)   mpz_set(d, eq->coefs[t]);   /* x y */
        else if (d1 == 1 && d2 == 1)   mpz_set(e, eq->coefs[t]);   /* y z */
        else if (d0 == 1 && d2 == 1)   mpz_set(f, eq->coefs[t]);   /* z x */
        else return false;
    }
    return true;
}


/* Remove square factors from n (n != 0): n <- squarefree part, sign preserved.
 * Trial division; declines (returns false) for |n| too large to factor quickly. */
static bool tg_squarefree(mpz_t n) {
    if (mpz_sgn(n) == 0) return false;
    mpz_t an; mpz_init(an); mpz_abs(an, n);
    if (mpz_cmp_ui(an, 1) == 0) { mpz_clear(an); return true; }
    /* Cap: the diagonal coefficients of a real form are small; decline giants. */
    if (mpz_sizeinbase(an, 2) > 80) { mpz_clear(an); return false; }
    int sgn = mpz_sgn(n);
    mpz_t p, p2, q; mpz_init_set_ui(p, 2); mpz_init(p2); mpz_init(q);
    mpz_set(q, an);
    while (mpz_mul(p2, p, p), mpz_cmp(p2, q) <= 0) {
        if (mpz_divisible_p(q, p2)) { mpz_divexact(q, q, p2); }  /* drop one p^2 */
        else mpz_add_ui(p, p, 1);
    }
    mpz_set(n, q); if (sgn < 0) mpz_neg(n, n);
    mpz_clear(an); mpz_clear(p); mpz_clear(p2); mpz_clear(q);
    return true;
}


/* One Legendre condition:  -o1*o2  is a quadratic residue modulo every odd prime
 * factor of |self|.  (Vacuous for the prime 2.) */
static bool tg_legendre_one(const mpz_t self, const mpz_t o1, const mpz_t o2) {
    mpz_t prod, m; mpz_init(prod); mpz_init(m);
    mpz_mul(prod, o1, o2); mpz_neg(prod, prod);
    mpz_abs(m, self);
    mpz_t primes[PR_MAX_DISTINCT_PRIMES]; size_t np = 0;
    for (size_t k = 0; k < PR_MAX_DISTINCT_PRIMES; k++) mpz_init(primes[k]);
    bool ok = pr_collect_distinct_primes(m, primes, &np);
    for (size_t k = 0; k < np && ok; k++) {
        if (mpz_cmp_ui(primes[k], 2) == 0) continue;          /* mod 2 vacuous */
        if (mpz_legendre(prod, primes[k]) == -1) ok = false;
    }
    for (size_t k = 0; k < PR_MAX_DISTINCT_PRIMES; k++) mpz_clear(primes[k]);
    mpz_clear(prod); mpz_clear(m);
    return ok;
}

/* Legendre solvability of  a X^2 + b Y^2 + c Z^2 == 0  with a,b,c squarefree,
 * pairwise coprime and NOT all the same sign: nontrivial iff -b c is a QR mod
 * |a|, -a c mod |b|, -a b mod |c| (odd primes only). */
static bool tg_legendre_ok(const mpz_t a, const mpz_t b, const mpz_t cc) {
    return tg_legendre_one(a, b, cc)
        && tg_legendre_one(b, a, cc)
        && tg_legendre_one(cc, a, b);
}


/* ---- mpq 3x3 helpers for the congruent diagonalisation. ---- */
typedef struct { mpq_t m[3][3]; } QMat;
static void qm_init(QMat* A) { for (int i=0;i<3;i++) for (int j=0;j<3;j++) mpq_init(A->m[i][j]); }
static void qm_clear(QMat* A){ for (int i=0;i<3;i++) for (int j=0;j<3;j++) mpq_clear(A->m[i][j]); }
static void qm_identity(QMat* A){ for (int i=0;i<3;i++) for (int j=0;j<3;j++) mpq_set_ui(A->m[i][j], i==j?1:0, 1); }
/* col_j <- col_j + lam*col_k (in place, over cols). */
static void qm_col_add(QMat* A, int j, int k, const mpq_t lam) {
    mpq_t t; mpq_init(t);
    for (int r = 0; r < 3; r++) { mpq_mul(t, lam, A->m[r][k]); mpq_add(A->m[r][j], A->m[r][j], t); }
    mpq_clear(t);
}
/* row_j <- row_j + lam*row_k. */
static void qm_row_add(QMat* A, int j, int k, const mpq_t lam) {
    mpq_t t; mpq_init(t);
    for (int cc = 0; cc < 3; cc++) { mpq_mul(t, lam, A->m[k][cc]); mpq_add(A->m[j][cc], A->m[j][cc], t); }
    mpq_clear(t);
}
static void qm_col_swap(QMat* A, int j, int k) { for (int r=0;r<3;r++) mpq_swap(A->m[r][j], A->m[r][k]); }
static void qm_row_swap(QMat* A, int j, int k) { for (int cc=0;cc<3;cc++) mpq_swap(A->m[j][cc], A->m[k][cc]); }

/* Congruently diagonalise the integer symmetric M (given as int64-free mpz):
 * find rational T with T^T M T = diag(D0,D1,D2).  Mw is worked in place; T
 * accumulates the column operations.  Returns false if the form is degenerate
 * (some diagonal entry stays 0 -> rank < 3). */
static bool tg_diagonalise(const mpz_t Mint[3][3], QMat* T, mpq_t D[3]) {
    QMat Mw; qm_init(&Mw);
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++)
        mpq_set_z(Mw.m[i][j], Mint[i][j]);
    qm_identity(T);
    mpq_t lam; mpq_init(lam);
    bool ok = true;
    for (int i = 0; i < 3 && ok; i++) {
        if (mpq_sgn(Mw.m[i][i]) == 0) {
            int piv = -1;
            for (int j = i + 1; j < 3; j++) if (mpq_sgn(Mw.m[j][j]) != 0) { piv = j; break; }
            if (piv >= 0) { qm_col_swap(&Mw, i, piv); qm_row_swap(&Mw, i, piv);
                            qm_col_swap(T, i, piv); }
            else {
                /* all diagonal (>= i) zero: use an off-diagonal to build one. */
                int fj = -1, fk = -1;
                for (int j = i; j < 3 && fj < 0; j++)
                    for (int k = j + 1; k < 3; k++)
                        if (mpq_sgn(Mw.m[j][k]) != 0) { fj = j; fk = k; break; }
                if (fj < 0) { ok = false; break; }        /* zero block -> degenerate */
                mpq_set_ui(lam, 1, 1);
                qm_col_add(&Mw, fj, fk, lam); qm_row_add(&Mw, fj, fk, lam);
                qm_col_add(T, fj, fk, lam);
                if (fj != i) { qm_col_swap(&Mw, i, fj); qm_row_swap(&Mw, i, fj);
                               qm_col_swap(T, i, fj); }
            }
        }
        if (!ok || mpq_sgn(Mw.m[i][i]) == 0) { ok = false; break; }
        for (int r = i + 1; r < 3; r++) {
            if (mpq_sgn(Mw.m[i][r]) == 0) continue;
            mpq_div(lam, Mw.m[i][r], Mw.m[i][i]); mpq_neg(lam, lam);   /* -M[i][r]/M[i][i] */
            qm_col_add(&Mw, r, i, lam); qm_row_add(&Mw, r, i, lam);
            qm_col_add(T, r, i, lam);
        }
        mpq_set(D[i], Mw.m[i][i]);
    }
    mpq_clear(lam); qm_clear(&Mw);
    if (ok) for (int i = 0; i < 3; i++) if (mpq_sgn(D[i]) == 0) ok = false;
    return ok;
}


/* Integer square root floor (>= 0). */
static void tg_isqrt(const mpz_t n, mpz_t r) { if (mpz_sgn(n) <= 0) mpz_set_ui(r, 0); else mpz_sqrt(r, n); }


/* Search for a witness of  Ca1*u_a1^2 + Ca2*u_a2^2 + Csi*u_si^2 == 0  over the
 * box  0 <= u_a1 <= b1, 0 <= u_a2 <= b2, solving u_si exactly (perfect k-th root).
 * On success writes the witness components to o_a1/o_a2/o_si and returns true. */
static bool tg_witness_search(const mpz_t Ca1, const mpz_t Ca2, const mpz_t Csi,
                              int64_t b1, int64_t b2,
                              mpz_t o_a1, mpz_t o_a2, mpz_t o_si) {
    mpz_t val, root, x, y; mpz_inits(val, root, x, y, NULL);
    bool found = false;
    for (int64_t x1 = 0; x1 <= b1 && !found; x1++)
        for (int64_t x2 = 0; x2 <= b2 && !found; x2++) {
            if (x1 == 0 && x2 == 0) continue;
            mpz_set_si(x, x1); mpz_set_si(y, x2);
            mpz_mul(val, Ca1, x); mpz_mul(val, val, x);
            mpz_mul(root, Ca2, y); mpz_mul(root, root, y);
            mpz_add(val, val, root); mpz_neg(val, val);       /* -(Ca1 x1^2 + Ca2 x2^2) */
            if (!mpz_divisible_p(val, Csi)) continue;
            mpz_divexact(val, val, Csi);                      /* u_si^2 */
            if (mpz_sgn(val) < 0 || !mpz_perfect_square_p(val)) continue;
            mpz_sqrt(root, val);
            mpz_set_si(o_a1, x1); mpz_set_si(o_a2, x2); mpz_set(o_si, root);
            found = true;
        }
    mpz_clears(val, root, x, y, NULL);
    return found;
}


/* Homogeneous ternary quadratic general solver.  Returns the owned family List,
 * the trivial-only list, or NULL to decline. */
Expr* si_solve_ternary_general(SICtx* c) {
    if (c->neq != 1 || c->n != 3) return NULL;
    for (int i = 0; i < 3; i++) if (c->has_lo[i] || c->has_hi[i]) return NULL;
    if (c->n_ord != 0 || c->n_neq != 0 || c->n_abs_ord != 0 || !c->all_captured)
        return NULL;

    int vs[3] = {0, 1, 2};
    for (int i = 0; i < 3; i++) if (mpoly_deg_var(c->eq[0], vs[i]) < 1) return NULL;

    mpz_t a, b, cc, d, e, f;
    mpz_inits(a, b, cc, d, e, f, NULL);
    Expr* result = NULL; bool declined = true;   /* pessimistic: set false on success */

    if (!tg_extract(c->eq[0], c->n, vs, a, b, cc, d, e, f)) goto done;

    /* Integer symmetric 2x matrix  M (G(v)=v^T M v = 2 Q(v)). */
    mpz_t M[3][3];
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) mpz_init(M[i][j]);
    mpz_mul_ui(M[0][0], a, 2); mpz_mul_ui(M[1][1], b, 2); mpz_mul_ui(M[2][2], cc, 2);
    mpz_set(M[0][1], d); mpz_set(M[1][0], d);
    mpz_set(M[1][2], e); mpz_set(M[2][1], e);
    mpz_set(M[0][2], f); mpz_set(M[2][0], f);

    QMat T; qm_init(&T);
    mpq_t D[3]; for (int i = 0; i < 3; i++) mpq_init(D[i]);
    bool diag_ok = tg_diagonalise(M, &T, D);

    if (diag_ok) {
        /* Integer diagonal form  C_i = num(D_i) * (L/den(D_i)),  L = lcm(dens),
         * proportional to sum D_i w_i^2 (so a witness u solves both). */
        mpz_t L, tmp; mpz_init_set_ui(L, 1); mpz_init(tmp);
        for (int i = 0; i < 3; i++) mpz_lcm(L, L, mpq_denref(D[i]));
        mpz_t C[3]; for (int i = 0; i < 3; i++) {
            mpz_init(C[i]); mpz_divexact(tmp, L, mpq_denref(D[i]));
            mpz_mul(C[i], mpq_numref(D[i]), tmp);
        }
        mpz_clear(L); mpz_clear(tmp);

        /* Reduced (squarefree, pairwise-coprime, content-free) copy for Legendre. */
        mpz_t r0, r1, r2, g; mpz_inits(r0, r1, r2, g, NULL);
        mpz_set(r0, C[0]); mpz_set(r1, C[1]); mpz_set(r2, C[2]);
        mpz_gcd(g, r0, r1); mpz_gcd(g, g, r2);
        if (mpz_sgn(g) != 0 && mpz_cmp_ui(g, 1) != 0) {
            mpz_divexact(r0, r0, g); mpz_divexact(r1, r1, g); mpz_divexact(r2, r2, g);
        }
        bool red_ok = tg_squarefree(r0) && tg_squarefree(r1) && tg_squarefree(r2);
        /* Make pairwise coprime: move a shared prime p|r_i,r_j to r_k (keeps them
         * squarefree since p divides neither r_k, by content-freeness). */
        for (int pass = 0; pass < 64 && red_ok; pass++) {
            bool changed = false;
            int pr[3][2] = {{0,1},{0,2},{1,2}};
            for (int t = 0; t < 3 && !changed; t++) {
                mpz_t* Ri = (pr[t][0]==0)?&r0:(pr[t][0]==1?&r1:&r2);
                mpz_t* Rj = (pr[t][1]==0)?&r0:(pr[t][1]==1?&r1:&r2);
                mpz_t* Rk = (pr[t][0]+pr[t][1]==1)?&r2:((pr[t][0]+pr[t][1]==2)?&r1:&r0);
                mpz_gcd(g, *Ri, *Rj);
                if (mpz_cmp_ui(g, 1) > 0) {
                    mpz_divexact(*Ri, *Ri, g); mpz_divexact(*Rj, *Rj, g); mpz_mul(*Rk, *Rk, g);
                    changed = true;
                }
            }
            if (!changed) break;
        }

        bool indefinite = !(mpz_sgn(r0) == mpz_sgn(r1) && mpz_sgn(r1) == mpz_sgn(r2));
        bool solvable = red_ok && indefinite && tg_legendre_ok(r0, r1, r2);

        if (red_ok && !solvable) {
            /* Definite or Legendre fails -> only the trivial point (a proof). */
            Expr* rules[3];
            for (int i = 0; i < 3; i++) rules[i] = mk_rule(expr_copy(c->var[vs[i]]), mk_int(0));
            Expr* tuple = mk_list(rules, 3);
            result = mk_list(&tuple, 1); declined = false;
        } else if (solvable) {
            /* Witness of  C_0 u0^2 + C_1 u1^2 + C_2 u2^2 = 0  within Holzer's box.
             * Pick the coordinate to solve for as the one whose C has sign opposite
             * to the other two's dominant sign; search the other two. */
            int si = 0;                                   /* solve-for index */
            for (int i = 0; i < 3; i++) {
                int other1 = (i+1)%3, other2 = (i+2)%3;
                if (mpz_sgn(C[i]) != mpz_sgn(C[other1]) && mpz_sgn(C[i]) != mpz_sgn(C[other2])) { si = i; break; }
                if (i == 2) si = 0;
            }
            int a1 = (si+1)%3, a2 = (si+2)%3;
            mpz_t B1, B2, prod, u[3];
            mpz_inits(B1, B2, prod, u[0], u[1], u[2], NULL);
            /* Holzer bounds: a solution exists with |u_a1| <= sqrt|C_a2 * C_si|
             * and |u_a2| <= sqrt|C_a1 * C_si| (each searched variable's bound
             * involves the solved-for coefficient, not the other searched one). */
            mpz_mul(prod, C[a2], C[si]); mpz_abs(prod, prod); tg_isqrt(prod, B1);
            mpz_mul(prod, C[a1], C[si]); mpz_abs(prod, prod); tg_isqrt(prod, B2);
            bool found = false;
            /* Small pre-scan first: a witness is often tiny even when the
             * coefficients (and hence the full Holzer box) are huge -- e.g.
             * 99991 x^2 - 99989 y^2 - 2 z^2 == 0 has the solution (1,1,1).  Cap
             * each dimension at SMALL so this never costs more than SMALL^2. */
            const int64_t SMALL = 3000;
            int64_t sb1 = mpz_cmp_ui(B1, SMALL) <= 0 ? mpz_get_si(B1) : SMALL;
            int64_t sb2 = mpz_cmp_ui(B2, SMALL) <= 0 ? mpz_get_si(B2) : SMALL;
            found = tg_witness_search(C[a1], C[a2], C[si], sb1, sb2, u[a1], u[a2], u[si]);
            /* Otherwise the full Holzer box, if it fits the node cap (Holzer
             * guarantees a witness inside it whenever the form is isotropic). */
            if (!found) {
                bool box_ok = mpz_cmp_ui(B1, 200000) <= 0 && mpz_cmp_ui(B2, 200000) <= 0;
                long double box = box_ok ? (2.0L*mpz_get_d(B1)+1)*(2.0L*mpz_get_d(B2)+1) : 1e30L;
                if (box_ok && box <= (long double)SI_MAX_NODES)
                    found = tg_witness_search(C[a1], C[a2], C[si],
                                              mpz_get_si(B1), mpz_get_si(B2), u[a1], u[a2], u[si]);
            }
            if (found) {
                /* P0 = T * u (rational); clear denominators, divide by content. */
                mpq_t pq[3]; for (int i = 0; i < 3; i++) mpq_init(pq[i]);
                mpq_t acc, uq; mpq_init(acc); mpq_init(uq);
                for (int i = 0; i < 3; i++) {
                    mpq_set_ui(pq[i], 0, 1);
                    for (int k = 0; k < 3; k++) {
                        mpq_set_z(uq, u[k]); mpq_mul(acc, T.m[i][k], uq); mpq_add(pq[i], pq[i], acc);
                    }
                }
                mpq_clear(acc); mpq_clear(uq);
                mpz_t P0[3], den, cont;
                mpz_inits(P0[0], P0[1], P0[2], den, cont, NULL);
                mpz_set_ui(den, 1);
                for (int i = 0; i < 3; i++) mpz_lcm(den, den, mpq_denref(pq[i]));
                for (int i = 0; i < 3; i++) {
                    mpz_divexact(cont, den, mpq_denref(pq[i]));
                    mpz_mul(P0[i], mpq_numref(pq[i]), cont);
                }
                mpz_set(cont, P0[0]); mpz_gcd(cont, cont, P0[1]); mpz_gcd(cont, cont, P0[2]);
                if (mpz_sgn(cont) != 0 && mpz_cmp_ui(cont, 1) != 0)
                    for (int i = 0; i < 3; i++) mpz_divexact(P0[i], P0[i], cont);
                for (int i = 0; i < 3; i++) mpq_clear(pq[i]);

                /* Verify G(P0) = 0 exactly (P0^T M P0). */
                mpz_t gcheck, tt; mpz_init_set_ui(gcheck, 0); mpz_init(tt);
                for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
                    mpz_mul(tt, M[i][j], P0[j]); mpz_mul(tt, tt, P0[i]); mpz_add(gcheck, gcheck, tt);
                }
                bool witness_ok = (mpz_sgn(gcheck) == 0)
                    && (mpz_sgn(P0[0]) || mpz_sgn(P0[1]) || mpz_sgn(P0[2]));
                mpz_clear(gcheck); mpz_clear(tt);

                if (witness_ok) {
                    /* Dominant axis L (largest |P0|); parametrise over the other two. */
                    int lp = 0;
                    for (int i = 1; i < 3; i++) if (mpz_cmpabs(P0[i], P0[lp]) > 0) lp = i;
                    int ii = (lp+1)%3, jj = (lp+2)%3;
                    /* gii=M[ii][ii], gij=M[ii][jj], gjj=M[jj][jj];
                     * mi = row_ii(M).P0,  mj = row_jj(M).P0. */
                    mpz_t gii, gij, gjj, mi, mj, s2;
                    mpz_inits(gii, gij, gjj, mi, mj, s2, NULL);
                    mpz_set(gii, M[ii][ii]); mpz_set(gij, M[ii][jj]); mpz_set(gjj, M[jj][jj]);
                    for (int k = 0; k < 3; k++) { mpz_mul(s2, M[ii][k], P0[k]); mpz_add(mi, mi, s2); }
                    for (int k = 0; k < 3; k++) { mpz_mul(s2, M[jj][k], P0[k]); mpz_add(mj, mj, s2); }

                    /* Component quadratic-form coefficients (s^2, s t, t^2). */
                    mpz_t co[3][3]; for (int r=0;r<3;r++) for (int cix=0;cix<3;cix++) mpz_init(co[r][cix]);
                    /* value_ii = (P0ii*gii - 2mi) s^2 + (2 P0ii gij - 2mj) st + (P0ii gjj) t^2 */
                    mpz_mul(co[ii][0], P0[ii], gii); mpz_submul_ui(co[ii][0], mi, 2);
                    mpz_mul(co[ii][1], P0[ii], gij); mpz_mul_ui(co[ii][1], co[ii][1], 2); mpz_submul_ui(co[ii][1], mj, 2);
                    mpz_mul(co[ii][2], P0[ii], gjj);
                    /* value_jj = (P0jj gii) s^2 + (2 P0jj gij - 2mi) st + (P0jj gjj - 2mj) t^2 */
                    mpz_mul(co[jj][0], P0[jj], gii);
                    mpz_mul(co[jj][1], P0[jj], gij); mpz_mul_ui(co[jj][1], co[jj][1], 2); mpz_submul_ui(co[jj][1], mi, 2);
                    mpz_mul(co[jj][2], P0[jj], gjj); mpz_submul_ui(co[jj][2], mj, 2);
                    /* value_lp = (P0lp gii) s^2 + (2 P0lp gij) st + (P0lp gjj) t^2 */
                    mpz_mul(co[lp][0], P0[lp], gii);
                    mpz_mul(co[lp][1], P0[lp], gij); mpz_mul_ui(co[lp][1], co[lp][1], 2);
                    mpz_mul(co[lp][2], P0[lp], gjj);

                    /* Content of the whole map (constant) -> divide out for completeness. */
                    mpz_t cont2; mpz_init_set_ui(cont2, 0);
                    for (int r=0;r<3;r++) for (int cix=0;cix<3;cix++) mpz_gcd(cont2, cont2, co[r][cix]);
                    if (mpz_sgn(cont2) != 0 && mpz_cmp_ui(cont2, 1) != 0)
                        for (int r=0;r<3;r++) for (int cix=0;cix<3;cix++) mpz_divexact(co[r][cix], co[r][cix], cont2);
                    mpz_clear(cont2);

                    /* Build the chord family (scale C[3]) and the tangent family
                     * (C[1]*P0), then dedup after evaluation. */
                    Expr* raw[2]; int nraw = 0;
                    {   Expr* rules[3];
                        for (int r = 0; r < 3; r++) {
                            Expr* terms[3]; int nt = 0;
                            if (mpz_sgn(co[r][0])) terms[nt++] = mk_fn2("Times", mk_mpz(co[r][0]),
                                mk_fn2("Power", mk_fn1("C", mk_int(1)), mk_int(2)));
                            if (mpz_sgn(co[r][1])) terms[nt++] = expr_new_function(mk_sym("Times"),
                                (Expr*[]){ mk_mpz(co[r][1]), mk_fn1("C", mk_int(1)), mk_fn1("C", mk_int(2)) }, 3);
                            if (mpz_sgn(co[r][2])) terms[nt++] = mk_fn2("Times", mk_mpz(co[r][2]),
                                mk_fn2("Power", mk_fn1("C", mk_int(2)), mk_int(2)));
                            Expr* q = (nt == 0) ? mk_int(0)
                                    : (nt == 1) ? terms[0]
                                    : expr_new_function(mk_sym("Plus"), terms, nt);
                            Expr* val = mk_fn2("Times", mk_fn1("C", mk_int(3)), q);
                            rules[r] = mk_rule(expr_copy(c->var[vs[r]]), val);
                        }
                        raw[nraw++] = mk_list(rules, 3);
                    }
                    {   Expr* rules[3];
                        for (int r = 0; r < 3; r++)
                            rules[r] = mk_rule(expr_copy(c->var[vs[r]]),
                                mk_fn2("Times", mk_mpz(P0[r]), mk_fn1("C", mk_int(1))));
                        raw[nraw++] = mk_list(rules, 3);
                    }
                    Expr* uniq[2]; int nu = 0;
                    for (int i = 0; i < nraw; i++) {
                        Expr* ev = eval_and_free(raw[i]);
                        bool dup = false;
                        for (int j = 0; j < nu; j++) if (expr_eq(ev, uniq[j])) { dup = true; break; }
                        if (dup) expr_free(ev); else uniq[nu++] = ev;
                    }
                    result = mk_list(uniq, (size_t)nu); declined = false;

                    for (int r=0;r<3;r++) for (int cix=0;cix<3;cix++) mpz_clear(co[r][cix]);
                    mpz_clears(gii, gij, gjj, mi, mj, s2, NULL);
                }
                mpz_clears(P0[0], P0[1], P0[2], den, cont, NULL);
            }
            mpz_clears(B1, B2, prod, u[0], u[1], u[2], NULL);
        }
        /* else: red_ok false (couldn't reduce) -> decline. */

        mpz_clears(r0, r1, r2, g, NULL);
        for (int i = 0; i < 3; i++) mpz_clear(C[i]);
    }

    for (int i = 0; i < 3; i++) mpq_clear(D[i]);
    qm_clear(&T);
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) mpz_clear(M[i][j]);

done:
    mpz_clears(a, b, cc, d, e, f, NULL);
    if (declined) { if (result) { expr_free(result); result = NULL; } }
    return result;
}
