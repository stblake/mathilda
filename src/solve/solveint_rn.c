/*
 * solveint_rn.c
 *
 * Part of the Solve[..., Integers] engine; split out of solveint.c.
 * See solveint_internal.h for the shared SICtx/SearchState substrate.
 *
 * Tier 3 F (route A) of SOLVE_INTEGERS.md: the unbounded Ramanujan-Nagell-type
 * exponential Diophantine  x^2 + D == 2^n  (b = 2).  For the class-number-1
 * imaginary-quadratic case (D squarefree, D == 7 (mod 8) so 2 splits in the
 * half-integer ring O_K = Z[(1+sqrt(-D))/2], and h(Q(sqrt(-D))) == 1), the
 * factorisation
 *     (x + sqrt(-D))/2 = +- alpha^(n-2),   alpha = (1+sqrt(-D))/2
 * forces  U_{n-2} = +-1  for the Lucas sequence  U_0=0, U_1=1,
 * U_m = U_{m-1} - Q U_{m-2},  Q = (1+D)/4.  By the Bilu-Hanrot-Voutier
 * primitive-divisor theorem every term with m > 30 has a primitive divisor,
 * so |U_m| = 1 forces m <= 30, i.e. n <= 32.  We then simply scan n in that
 * finite window for  2^n - D  a perfect square -- an exact, complete search --
 * and (defensively) cross-check the found set against the Lucas condition.
 *
 * For b = 2 the gate singles out exactly D = 7 (the classical Ramanujan-Nagell
 * equation, 2^n - 7 = x^2 with n in {3,4,5,7,15}); the checks are kept general
 * (not hard-coded to 7) so the method documents its own correctness and extends
 * cleanly.  Everything outside the gate DECLINES (route B -- linear forms in
 * logarithms -- is a documented future extension, never a guess).
 */
#include "solveint.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

#include "eval.h"
#include "expr.h"
#include "sym_names.h"
#include "poly/mpoly.h"
#include "solveint_internal.h"


static long rn_gcd(long a, long b) { if (a < 0) a = -a; if (b < 0) b = -b;
    while (b) { long t = a % b; a = b; b = t; } return a; }

/* Squarefree test for a positive long. */
static bool rn_squarefree(long m) {
    if (m <= 0) return false;
    for (long p = 2; p * p <= m; p++)
        if (m % p == 0) { m /= p; if (m % p == 0) return false; }
    return true;
}

/* Class number of the imaginary quadratic field of fundamental discriminant
 * D < 0 (here D = -Ddisc with Ddisc == 3 (mod 4), so D == 1 (mod 4) is
 * fundamental), by counting primitive reduced positive-definite forms (A,B,C)
 * with B^2 - 4AC = D: reduced iff |B| <= A <= C, and B >= 0 on the boundary
 * (|B| == A or A == C).  Mirrors si_class_number_neg4c in solveint_mordell.c
 * but for the odd fundamental discriminant. */
static long rn_class_number_disc(long D) {          /* D < 0 */
    long h = 0;
    long Amax = (long)floor(sqrt((double)(-D) / 3.0)) + 1;
    for (long A = 1; A <= Amax; A++)
        for (long B = -A + 1; B <= A; B++) {
            long num = B * B - D;                   /* = 4AC */
            if (num % (4 * A) != 0) continue;
            long C = num / (4 * A);
            if (C < A) continue;
            if (rn_gcd(rn_gcd(A, B), C) != 1) continue;
            if ((labs(B) == A || A == C) && B < 0) continue;
            h++;
        }
    return h;
}


/* Ramanujan-Nagell-type solver.  Returns the owned finite solution List (a
 * proven {} if empty within the gate), or NULL to decline. */
Expr* si_solve_ramanujan_nagell(Expr* expr, Expr** var, int n) {
    if (n != 2) return NULL;                          /* exactly {n, x} */

    Expr** conj; int ncj;
    flatten_conjuncts(expr, &conj, &ncj);
    Expr* eqn = NULL;
    for (int i = 0; i < ncj; i++)
        if (is_fun(conj[i], SYM_Equal, 2)) { eqn = conj[i]; break; }
    if (!eqn) return NULL;

    /* Parse both sides into additive power terms. */
    SIExpTerm t[8]; int nt = 0;
    mpz_t K; mpz_init_set_ui(K, 0);
    bool okp = si_exp_collect(eqn->data.function.args[0], +1, var, n, t, &nt, K)
            && si_exp_collect(eqn->data.function.args[1], -1, var, n, t, &nt, K);
    mpz_neg(K, K);                                    /* sum of power terms == K */
    if (!okp || nt != 2) { mpz_clear(K); return NULL; }

    /* Identify  ib = the b^n term (constant base >= 2, variable exponent) and
     * ix = the x^2 term (variable base, constant exponent 2). */
    int ib = -1, ix = -1;
    for (int i = 0; i < 2; i++) {
        if (t[i].bvar < 0 && t[i].bconst >= 2 && t[i].evar >= 0) ib = i;
        else if (t[i].bvar >= 0 && t[i].evar < 0 && t[i].econst == 2) ix = i;
    }
    if (ib < 0 || ix < 0 || ib == ix) { mpz_clear(K); return NULL; }
    if (!((t[ib].coef == 1 && t[ix].coef == -1) ||
          (t[ib].coef == -1 && t[ix].coef == 1))) { mpz_clear(K); return NULL; }

    int64_t b = t[ib].bconst;
    int nidx = t[ib].evar, xidx = t[ix].bvar;
    if (nidx == xidx) { mpz_clear(K); return NULL; }

    /* Normalise to  x^2 + D == b^n :  D = coef(b^n) * K. */
    mpz_t D; mpz_init(D);
    if (t[ib].coef == 1) mpz_set(D, K); else mpz_neg(D, K);
    mpz_clear(K);

    /* --- Route A gate --- */
    if (b != 2 || mpz_sgn(D) <= 0 || !mpz_fits_slong_p(D)) { mpz_clear(D); return NULL; }
    long Dl = mpz_get_si(D);
    mpz_clear(D);
    if (!rn_squarefree(Dl)) return NULL;
    if (((Dl % 8) + 8) % 8 != 7) return NULL;         /* 2 splits, half-integer ring */
    if (rn_class_number_disc(-Dl) != 1) return NULL;  /* h(Q(sqrt(-D))) == 1 */

    /* --- BHV window scan --- */
    const int N_MAX = 40;                             /* > 32 (BHV), harmless superset */
    int64_t sn[64], sx[64]; int ns = 0;               /* verified (n, x) pairs */
    int scan_ns[64]; int n_scan = 0;                  /* n with 2^n - D a square */
    mpz_t T, rt; mpz_init(T); mpz_init(rt);
    for (int nn = 0; nn <= N_MAX; nn++) {
        mpz_ui_pow_ui(T, 2, (unsigned)nn);
        mpz_sub_ui(T, T, (unsigned long)Dl);          /* 2^nn - D */
        if (mpz_sgn(T) <= 0) continue;
        if (!mpz_perfect_square_p(T)) continue;
        mpz_sqrt(rt, T);
        if (!mpz_fits_slong_p(rt)) continue;          /* unreachable for nn <= 40 */
        int64_t xval = mpz_get_si(rt);
        scan_ns[n_scan++] = nn;
        for (int sgn = 1; sgn >= -1; sgn -= 2) {      /* +x and -x (x^2 is even) */
            int64_t vals[SI_MAX_VARS];
            for (int i = 0; i < n; i++) vals[i] = 0;
            vals[nidx] = nn; vals[xidx] = sgn * xval;
            if (si_verify_symbolic(expr, var, vals, n) && ns < 64) {
                bool dup = false;
                for (int j = 0; j < ns; j++)
                    if (sn[j] == nn && sx[j] == sgn * xval) { dup = true; break; }
                if (!dup) { sn[ns] = nn; sx[ns] = sgn * xval; ns++; }
            }
            if (xval == 0) break;                     /* +0 == -0 */
        }
    }
    mpz_clear(T); mpz_clear(rt);

    /* --- Lucas cross-check (defensive): {m+2 : |U_m| = 1, 1 <= m <= N_MAX-2}
     * must equal the perfect-square scan set (n >= 3).  A mismatch means the
     * number theory / gate is wrong -> DECLINE rather than emit. */
    {
        mpz_t Um2, Um1, Um, tmp; mpz_init_set_ui(Um2, 0); mpz_init_set_ui(Um1, 1);
        mpz_init(Um); mpz_init(tmp);
        long Q = (Dl + 1) / 4;                        /* integer since D == 3 (mod 4) */
        int lucas_ns[64]; int nl = 0;
        /* m = 1 -> U_1 = 1 -> n = 3 */
        if (mpz_cmpabs_ui(Um1, 1) == 0) lucas_ns[nl++] = 1 + 2;
        for (int m = 2; m <= N_MAX - 2; m++) {
            mpz_mul_si(tmp, Um2, Q);                  /* Q U_{m-2} */
            mpz_sub(Um, Um1, tmp);                    /* U_m = U_{m-1} - Q U_{m-2} */
            if (mpz_cmpabs_ui(Um, 1) == 0 && nl < 64) lucas_ns[nl++] = m + 2;
            mpz_set(Um2, Um1); mpz_set(Um1, Um);
        }
        mpz_clear(Um2); mpz_clear(Um1); mpz_clear(Um); mpz_clear(tmp);

        /* Compare {scan n >= 3} with lucas_ns as sets. */
        int cs = 0; for (int i = 0; i < n_scan; i++) if (scan_ns[i] >= 3) cs++;
        bool match = (cs == nl);
        for (int i = 0; i < n_scan && match; i++) {
            if (scan_ns[i] < 3) continue;
            bool found = false;
            for (int j = 0; j < nl; j++) if (lucas_ns[j] == scan_ns[i]) { found = true; break; }
            if (!found) match = false;
        }
        if (!match) return NULL;                      /* cross-check failed: decline */
    }

    /* Sort verified pairs ascending by (n, x). */
    for (int i = 0; i < ns; i++)
        for (int j = i + 1; j < ns; j++)
            if (sn[j] < sn[i] || (sn[j] == sn[i] && sx[j] < sx[i])) {
                int64_t a = sn[i]; sn[i] = sn[j]; sn[j] = a;
                a = sx[i]; sx[i] = sx[j]; sx[j] = a;
            }

    Expr** sols = (Expr**)malloc(sizeof(Expr*) * (size_t)(ns ? ns : 1));
    for (int i = 0; i < ns; i++) {
        int64_t vals[SI_MAX_VARS];
        for (int k = 0; k < n; k++) vals[k] = 0;
        vals[nidx] = sn[i]; vals[xidx] = sx[i];
        sols[i] = si_exp_one_tuple(var, vals, n);
    }
    Expr* result = mk_list(sols, (size_t)ns);         /* {} is a proof within the gate */
    free(sols);
    return result;
}
