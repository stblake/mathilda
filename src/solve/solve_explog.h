/*
 * solve_explog.h -- exact log-linear algebra for exponential equations.
 *
 * An equation whose variable occurs only in exponents of variable-free
 * numeric bases,
 *
 *      c1 * prod b_i^g_i(x)  ==  c2 * prod b_j^g_j(x)     (or  == K),
 *
 * is, after taking logarithms, LINEAR in the logarithms of its bases:
 *
 *      sum_i s_i g_i(x) Log[b_i] - Log[K]  ==  2 Pi I C,   C in Integers.
 *
 * This module makes that linear form EXACT.  Every positive rational that
 * occurs (bases, the constant K) is written over a multiplicative basis --
 * pairwise-coprime integers p_k that are not perfect powers, computed by gcd
 * refinement without any integer factoring -- so each Log[q] becomes an
 * integer/rational vector over the Q-linearly independent Log[p_k].  (Coprime
 * integers > 1 have disjoint prime supports, so their logarithms, together
 * with 1 = Log[E], are linearly independent over Q by unique factorisation
 * and Lindemann's theorem.)  The consequences are decision procedures, not
 * heuristics:
 *
 *   - Log[2^201]/Log[8] is recognised as 67: the coefficient matrix has rank 1
 *     and the equation collapses to a RATIONAL polynomial in x,
 *     x - 67 == 2 Pi I C/Log[8].
 *   - a rational x solves sum_k E_k(x) Log[p_k] == 0 iff it is a common root
 *     of every row polynomial E_k, so Solve over Integers/Rationals is exact.
 *   - for real x and positive real bases the logarithm is a bijection, so over
 *     the Reals only the C = 0 branch exists (and K <= 0 has no solution).
 *
 * The kernel helper handles sums of exponentials over one common generator
 * (4^x - 3 2^x + 2 == 0) by the substitution u = B^(d x).
 *
 * The functions here are pure algebra; the Solve glue (parameter minting,
 * branch wrapping, domain dispatch) lives in solveinv.c.
 */
#ifndef SOLVE_EXPLOG_H
#define SOLVE_EXPLOG_H

#include <gmp.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "expr.h"

/* ------------------------------------------------------------------ *
 *  Multiplicative basis of positive integers.                         *
 * ------------------------------------------------------------------ */

typedef struct {
    size_t n;     /* number of basis elements                          */
    mpz_t* p;     /* ascending, pairwise coprime, each > 1, primitive  */
} ExplogBasis;

/* Build the basis of the positive integers vals[0..nv) (values <= 1 are
 * ignored).  Every input is then an exact product of basis powers.  Uses gcd
 * refinement (mpz_remove), never factoring, so a 10^6-bit power is cheap. */
void explog_basis_build(const mpz_t* vals, size_t nv, ExplogBasis* B);

/* Exponent vector of the positive integer v over B: v = prod p_k^e[k].
 * `e` has B->n entries.  Returns false when v is not representable (a prime
 * factor outside the basis). */
bool explog_basis_exponents(const ExplogBasis* B, const mpz_t v, int64_t* e);

void explog_basis_clear(ExplogBasis* B);

/* ------------------------------------------------------------------ *
 *  Linearisation of an exponential equation.                          *
 * ------------------------------------------------------------------ */

typedef enum {
    EXPLOG_NO_SHAPE = 0,   /* not an exponential equation of this form     */
    EXPLOG_NO_SOLUTION,    /* provably no solution (prod of exps == 0, ...) */
    EXPLOG_LINEAR          /* lhs/period below are filled                   */
} ExplogKind;

typedef enum {
    EXPLOG_RHS_POSITIVE = 0,   /* the constant side is a positive real      */
    EXPLOG_RHS_NONE,           /* provably <= 0 or non-real: no real branch */
    EXPLOG_RHS_COND            /* positive exactly when real_cond holds     */
} ExplogRhs;

typedef struct {
    ExplogKind kind;
    /* Every complex solution satisfies  lhs == period * C,  C in Integers. */
    Expr* lhs;
    Expr* period;
    /* Real branch.  real_ok: all bases are positive reals and every exponent
     * is real-valued on the reals, so for real x the solutions are exactly
     * those of lhs == 0 (when rhs_real != NONE).  real_cond is owned and set
     * only for EXPLOG_RHS_COND. */
    bool      real_ok;
    ExplogRhs rhs_real;
    Expr*     real_cond;
    /* Rational decision polynomial: when non-NULL, the rational solutions of
     * lhs == 0 are exactly the rational roots of this polynomial (a nonzero
     * constant means there are none).  NULL when the logarithms involved are
     * not known to be independent (opaque bases such as Pi). */
    Expr* rational_poly;
} ExplogLinear;

/* Linearise `residual` (an evaluated lhs - rhs) in `var`.  Returns false and
 * sets kind = EXPLOG_NO_SHAPE when the residual is not of the supported
 * shape; otherwise fills *out (release with explog_linear_free). */
bool explog_linearize(const Expr* residual, const Expr* var, ExplogLinear* out);
void explog_linear_free(ExplogLinear* L);

/* ------------------------------------------------------------------ *
 *  Polynomial-in-one-exponential kernel.                              *
 * ------------------------------------------------------------------ */

typedef struct {
    Expr* poly;     /* polynomial in u (owned)                         */
    Expr* u;        /* the kernel symbol (owned)                       */
    Expr* kernel;   /* B^(d var): substitute back u == kernel (owned)  */
} ExplogKernel;

/* A sum of at least two exponential terms over ONE multiplicative generator,
 * with exponents linear in var with rational slopes, plus an optional free
 * term:  sum_t c_t B^(n_t d x) + f  ->  sum_t c_t u^n_t + f  (u = B^(d x)).
 * Returns false when the residual is not of that shape. */
bool explog_kernel(const Expr* residual, const Expr* var, ExplogKernel* out);
void explog_kernel_free(ExplogKernel* K);

#endif /* SOLVE_EXPLOG_H */
