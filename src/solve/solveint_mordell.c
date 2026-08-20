/*
 * solveint_mordell.c
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

static long si_lgcd(long a, long b) { a = labs(a); b = labs(b); while (b) { long t = a % b; a = b; b = t; } return a; }


/* Is m (|m|) squarefree? */
static bool si_is_squarefree_l(long m) {
    m = labs(m);
    if (m == 0) return false;
    for (long p = 2; p * p <= m; p++)
        if (m % p == 0) { m /= p; if (m % p == 0) return false; }
    return true;
}


/* Class number of the imaginary quadratic field of fundamental discriminant
 * D = -4c (c >= 1 squarefree, c = 1,2 mod 4), by counting primitive reduced
 * positive-definite binary quadratic forms (A,B,C) with B^2 - 4AC = D:
 * reduced iff |B| <= A <= C and B >= 0 when |B| == A or A == C. */
static long si_class_number_neg4c(long c) {
    long D = -4 * c, h = 0;
    long Amax = (long)floor(sqrt((double)(4 * c) / 3.0)) + 1;
    for (long A = 1; A <= Amax; A++)
        for (long B = -A + 1; B <= A; B++) {
            long num = B * B - D;                     /* = 4AC */
            if (num % (4 * A) != 0) continue;
            long C = num / (4 * A);
            if (C < A) continue;
            if (si_lgcd(si_lgcd(A, B), C) != 1) continue;    /* primitive */
            if ((labs(B) == A || A == C) && B < 0) continue; /* boundary: B >= 0 */
            h++;
        }
    return h;
}


/* The Z[sqrt k] cube descent below is sound AND complete exactly when:
 *   - k < 0 and |k| is squarefree (so D = -4|k| is a fundamental discriminant
 *     and the field is imaginary),
 *   - k = 2,3 (mod 4)  (ring of integers is Z[sqrt k]; unit group {+/-1}, and a
 *     mod-8 argument then forces x ODD and coprime to k, so the two factors of
 *     x^3 = (y - sqrt k)(y + sqrt k) are ALWAYS coprime -- no non-coprime case
 *     analysis is needed),
 *   - 3 does not divide the class number h  (so an ideal whose cube is principal
 *     is itself principal, and every point comes from y + sqrt k = (a+b sqrt k)^3).
 * Outside this set (k = 1 mod 4 -> half-integer ring; 3 | h; k > 0 -> real field
 * with infinite units) the descent is declined, never guessed. */
static bool si_mordell_descent_sound(long k) {
    if (k >= 0 || !si_is_squarefree_l(k)) return false;
    int r = (int)(((k % 4) + 4) % 4);
    if (r != 2 && r != 3) return false;
    return si_class_number_neg4c(-k) % 3 != 0;
}


/* --- Unbounded Mordell  y^2 == x^3 + k  via factorisation in Z[sqrt k]. ---
 *
 * x^3 = y^2 - k = (y - sqrt k)(y + sqrt k).  For the sound set of k (see
 * si_mordell_descent_sound) the factors are coprime, the class number is prime
 * to 3, and the units are {+/-1}, so every point satisfies
 *   y + sqrt k = (a + b sqrt k)^3,   1 = 3 a^2 b + k b^3 = b (3 a^2 + k b^2),
 * forcing b = +/-1 and a fixed, giving the COMPLETE integer-point set
 *   x = a^2 - k b^2,   y = +/- (a^3 + 3 a k b^2).
 * So y^2 = x^3 - 2 -> (3, +/-5), y^2 = x^3 - 1 -> (1, 0), and y^2 = x^3 - 5 ->
 * {} (proved, not merely unfound).  Only engaged when x is unbounded (a bounded
 * box is handled by the ordinary leaf search). */
Expr* si_solve_mordell(SICtx* c) {
    if (c->neq != 1) return NULL;
    const MPoly* eq = c->eq[0]; int n = c->n;
    int act[SI_MAX_VARS], ka = 0;
    for (int v = 0; v < n; v++) if (mpoly_deg_var(eq, v) >= 1) act[ka++] = v;
    if (ka != 2) return NULL;

    /* Identify Y (only Y^2) and X (only X^3), monic and opposite sign, plus a
     * constant.  Normalise to  Y^2 - X^3 - k == 0  (so k = -const). */
    for (int choice = 0; choice < 2; choice++) {
        int Yv = act[choice], Xv = act[1 - choice];
        mpz_t cY, cX, c0; mpz_init_set_ui(cY, 0); mpz_init_set_ui(cX, 0); mpz_init_set_ui(c0, 0);
        bool shape = true;
        for (size_t t = 0; t < eq->n_terms && shape; t++) {
            const int* ex = eq->exps + t * (size_t)n;
            for (int v = 0; v < n; v++) if (v != Yv && v != Xv && ex[v] != 0) shape = false;
            if (!shape) break;
            int ey = ex[Yv], exx = ex[Xv];
            if (ey == 2 && exx == 0) mpz_add(cY, cY, eq->coefs[t]);
            else if (ey == 0 && exx == 3) mpz_add(cX, cX, eq->coefs[t]);
            else if (ey == 0 && exx == 0) mpz_add(c0, c0, eq->coefs[t]);
            else shape = false;
        }
        bool ok = shape && (mpz_cmpabs_ui(cY, 1) == 0) && (mpz_cmpabs_ui(cX, 1) == 0)
               && (mpz_sgn(cY) == -mpz_sgn(cX));
        long k = 0;
        if (ok) {
            if (mpz_sgn(cY) < 0) { mpz_neg(cY, cY); mpz_neg(cX, cX); mpz_neg(c0, c0); }
            if (mpz_fits_slong_p(c0)) k = -mpz_get_si(c0);      /* y^2 = x^3 + k */
            else ok = false;
        }
        /* Engage on the sound descent set (any imaginary k with k = 2,3 mod 4,
         * |k| squarefree, 3 not dividing the class number), and only when X is
         * unbounded (a bounded box is handled by the ordinary leaf search). */
        bool unbounded = !(c->has_lo[Xv] && c->has_hi[Xv]);
        if (ok && si_mordell_descent_sound(k) && unbounded) {
            Expr** sols = NULL; int nsol = 0, cap = 0;
            for (int b = 1; b >= -1; b -= 2) {
                /* b (3 a^2 + k b^2) = 1  ->  3 a^2 = (1/b) - k b^2 = b - k b^2. */
                long rhs = b - k * b * b;                /* = 3 a^2 */
                if (rhs < 0 || rhs % 3 != 0) continue;
                long a2 = rhs / 3;
                long a0 = (long)llround(sqrt((double)a2));
                for (long a = a0 - 1; a <= a0 + 1; a++) {  /* pin the exact root */
                    if (a < 0 || a * a != a2) continue;
                    for (int as = (a == 0 ? 1 : 1); as >= (a == 0 ? 1 : -1); as -= 2) {
                        long av = as * a;
                        long xv = av * av - k * b * b;             /* norm */
                        long yv = av * av * av + 3 * av * k * b * b;
                        for (int ys = 1; ys >= -1; ys -= 2) {      /* unit +/-1 -> +/- y */
                            int64_t vals[SI_MAX_VARS];
                            for (int i = 0; i < n; i++) vals[i] = 0;
                            vals[Xv] = xv; vals[Yv] = ys * yv;
                            bool dup = false;
                            if (si_verify_symbolic(c->original, c->var, vals, n)) {
                                Expr* cand = si_exp_one_tuple(c->var, vals, n);
                                for (int u = 0; u < nsol && !dup; u++) if (expr_eq(sols[u], cand)) dup = true;
                                if (dup) { expr_free(cand); }
                                else { if (nsol == cap) { cap = cap ? cap * 2 : 4; sols = realloc(sols, sizeof(Expr*) * (size_t)cap); } sols[nsol++] = cand; }
                            }
                        }
                    }
                }
            }
            mpz_clear(cY); mpz_clear(cX); mpz_clear(c0);
            /* sort ascending by (x, y) for a canonical result */
            for (int i = 0; i + 1 < nsol; i++) for (int j = i + 1; j < nsol; j++)
                if (expr_compare(sols[j], sols[i]) < 0) { Expr* tt = sols[i]; sols[i] = sols[j]; sols[j] = tt; }
            Expr* result = mk_list(sols, (size_t)nsol);
            free(sols);
            return result;
        }
        mpz_clear(cY); mpz_clear(cX); mpz_clear(c0);
    }
    return NULL;
}
