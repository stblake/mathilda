/* chineseremainder.c -- ChineseRemainder[].
 * Split from numbertheory.c; see numbertheory.h and
 * numbertheory_internal.h for the subsystem layout. */

#include "numbertheory.h"
#include "numbertheory_internal.h"
#include "arithmetic.h"
#include "eval.h"
#include "sym_names.h"
#include "internal.h"
#include "print.h"
#include "symtab.h"
#include "attr.h"
#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

/* ----------------------------------------------------------------------
 * ChineseRemainder[{r1, r2, ...}, {m1, m2, ...}]
 *     the smallest x >= 0 satisfying x mod mi == ri mod mi for all i,
 *     lying in 0 <= x < lcm(m1, m2, ...).
 * ChineseRemainder[{r1, r2, ...}, {m1, m2, ...}, d]
 *     the smallest x >= d with the same congruences, in d <= x < d + lcm.
 *
 * The moduli need NOT be pairwise coprime.  A solution exists iff every
 * pair of congruences agrees modulo gcd(mi, mj); when the system is
 * inconsistent ChineseRemainder returns unevaluated (this builtin returns
 * NULL), matching Mathematica -- e.g. ChineseRemainder[{1, 2}, {6, 10}].
 *
 * Integer-only: machine ints and GMP bigints share one path (expr_to_mpz
 * auto-promotes; the result demotes back to EXPR_INTEGER via
 * expr_bigint_normalize when it fits), so a residue-number-system recovery
 * over large coprime moduli yields the full bignum.
 *
 * Algorithm: a streaming pairwise CRT fold.  The accumulator carries a
 * single congruence  x == x_acc (mod m_acc),  starting from (0, 1) (the
 * always-true congruence).  Each step merges in (ri, |mi|):
 *
 *     g = gcd(m_acc, mi),   m_acc*u + mi*v = g          (mpz_gcdext)
 *     solvable iff  g | (ri - x_acc)
 *     x_acc <- x_acc + ((ri - x_acc) / g) * u * m_acc
 *     m_acc <- lcm(m_acc, mi) = m_acc / g * mi
 *     x_acc <- x_acc mod m_acc                          (reduce to [0, L))
 *
 * The new x_acc is == x_acc (mod m_acc) since the added term is a multiple
 * of m_acc, and == ri (mod mi) since m_acc*u == g (mod mi) makes the added
 * term == (ri - x_acc) (mod mi).  Folding keeps the accumulator reduced, so
 * the intermediates stay bounded by the running lcm.
 *
 * Finally the smallest x >= d (with d = 0 for the two-argument form) is
 *     x = d + ((x_acc - d) mod L),
 * which lands in [d, d + L) and preserves x == x_acc (mod L); this also
 * handles a negative d and the empty-list case (L = 1) correctly.
 * ------------------------------------------------------------------- */

Expr* builtin_chineseremainder(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3)
        return builtin_arg_error("ChineseRemainder", argc, 2, 3);

    Expr* rlist = res->data.function.args[0];
    Expr* mlist = res->data.function.args[1];

    /* The residues and moduli must be Lists of equal length; anything else
     * (symbolic arguments, a visible NDArray, mismatched lengths) leaves the
     * call unevaluated so symbolic flow and rules still apply. */
    if (!head_is(rlist, SYM_List) || !head_is(mlist, SYM_List)) return NULL;
    size_t n = rlist->data.function.arg_count;
    if (mlist->data.function.arg_count != n) return NULL;

    /* Optional offset d must be an exact integer. */
    Expr* d_expr = (argc == 3) ? res->data.function.args[2] : NULL;
    if (d_expr && !expr_is_integer_like(d_expr)) return NULL;

    /* Every residue and modulus must be an exact integer. */
    for (size_t i = 0; i < n; i++) {
        if (!expr_is_integer_like(rlist->data.function.args[i]) ||
            !expr_is_integer_like(mlist->data.function.args[i]))
            return NULL;
    }

    /* --- Streaming pairwise CRT fold (all GMP). --- */
    mpz_t x_acc, m_acc, g, u, v, diff, t, lcm;
    mpz_inits(x_acc, m_acc, g, u, v, diff, t, lcm, NULL);
    mpz_set_ui(x_acc, 0);
    mpz_set_ui(m_acc, 1);

    bool ok = true;
    for (size_t i = 0; i < n; i++) {
        mpz_t ri, mi;
        expr_to_mpz(rlist->data.function.args[i], ri); /* inits ri */
        expr_to_mpz(mlist->data.function.args[i], mi); /* inits mi */

        /* A zero modulus has no residues to constrain against -- reject it
         * (it would also make the lcm zero and break the [0, L) contract). */
        if (mpz_sgn(mi) == 0) { ok = false; mpz_clears(ri, mi, NULL); break; }
        mpz_abs(mi, mi);

        mpz_gcdext(g, u, v, m_acc, mi);      /* g = gcd(m_acc, mi) */
        mpz_sub(diff, ri, x_acc);
        if (!mpz_divisible_p(diff, g)) {     /* congruences disagree mod g */
            ok = false;
            mpz_clears(ri, mi, NULL);
            break;
        }
        mpz_divexact(diff, diff, g);         /* (ri - x_acc) / g */
        mpz_mul(t, diff, u);                 /* * u */
        mpz_mul(t, t, m_acc);                /* * m_acc */
        mpz_add(x_acc, x_acc, t);            /* x_acc + ((ri - x_acc)/g) u m_acc */

        mpz_divexact(lcm, m_acc, g);         /* lcm = m_acc / g * mi */
        mpz_mul(lcm, lcm, mi);
        mpz_set(m_acc, lcm);
        mpz_mod(x_acc, x_acc, m_acc);        /* reduce into [0, lcm) */

        mpz_clears(ri, mi, NULL);
    }

    if (!ok) {
        mpz_clears(x_acc, m_acc, g, u, v, diff, t, lcm, NULL);
        return NULL;                         /* no solution -> unevaluated */
    }

    /* Smallest x >= d congruent to x_acc modulo the total lcm. */
    if (d_expr) {
        mpz_t d;
        expr_to_mpz(d_expr, d);              /* inits d */
        mpz_sub(t, x_acc, d);
        mpz_mod(t, t, m_acc);                /* [0, lcm) */
        mpz_add(x_acc, d, t);
        mpz_clear(d);
    } else {
        mpz_mod(x_acc, x_acc, m_acc);        /* [0, lcm); == 0 for empty lists */
    }

    Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(x_acc));
    mpz_clears(x_acc, m_acc, g, u, v, diff, t, lcm, NULL);
    return out;
}
