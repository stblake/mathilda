/* chineseremainder.c -- ChineseRemainder[].
 * Split from numbertheory.c; see numbertheory.h and
 * numbertheory_internal.h for the subsystem layout. */

#include "numbertheory.h"
#include "linkedlist.h"
#include "numbertheory_internal.h"
#include "arithmetic.h"
#include "eval.h"
#include "sym_names.h"
#include "internal.h"
#include "print.h"
#include "symtab.h"
#include "attr.h"
#include "message.h"   /* mth_message_gated: Quiet/Check funnel */
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

/* A number of any kind: what separates ::pilist from a symbolic call. */
static bool cr_is_number(Expr* e) {
    Expr *re, *im;
    return is_rational_like(e) || e->type == EXPR_REAL || is_complex(e, &re, &im)
#ifdef USE_MPFR
        || e->type == EXPR_MPFR
#endif
        ;
}

/* The CRT over residues rv[0..nr) and moduli mv[0..nm): the elements of two
 * Lists, or of LinkedLists read from their slots. */
/* Element i of a List (borrowed) or a LinkedList (its slot, held in *hold
 * until the next read): residues and moduli are read one at a time. */
static Expr* cr_at(const Expr* seq, size_t i, Expr** hold) {
    if (seq->type != EXPR_LINKEDLIST) return seq->data.function.args[i];
    expr_free(*hold);
    *hold = ll_get(seq->data.llist, (int64_t)i);
    return *hold;
}

static Expr* cr_core_impl(const Expr* rs, size_t nr, const Expr* ms, size_t nm, Expr* d_expr,
                          Expr** hr, Expr** hm) {
    /* Symbolic arguments leave the call unevaluated silently, so symbolic flow
     * and rules still apply. */
    size_t n = nr;
    if (n == 0 && nm == 0) return NULL;
    for (size_t k = 0; k < 2; k++) {
        const Expr* l = k ? ms : rs;
        size_t ln = k ? nm : nr;
        for (size_t i = 0; i < ln; i++)
            if (!cr_is_number(cr_at(l, i, k ? hm : hr))) return NULL;
    }

    /* Numbers, but not two equal-length lists of integers with positive moduli
     * -- ChineseRemainder::pilist. */
    bool valid = nm == n;
    for (size_t i = 0; i < n && valid; i++) {
        Expr* r = cr_at(rs, i, hr);
        Expr* m = cr_at(ms, i, hm);
        if (!expr_is_integer_like(r) || !expr_is_integer_like(m)) valid = false;
        else if (m->type == EXPR_INTEGER ? m->data.integer <= 0
                                         : mpz_sgn(m->data.bigint) <= 0) valid = false;
    }
    if (!valid) {
        mth_message_gated(g_arith_warnings_muted, "ChineseRemainder", "pilist",
            "The arguments to ChineseRemainder must be two lists of integers of "
            "identical length, with the second list containing only positive integers.");
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
        expr_to_mpz(cr_at(rs, i, hr), ri); /* inits ri */
        expr_to_mpz(cr_at(ms, i, hm), mi); /* inits mi */

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

static Expr* cr_core(const Expr* rs, const Expr* ms, Expr* d_expr) {
    Expr *hr = NULL, *hm = NULL;
    Expr* r = cr_core_impl(rs, (size_t)ll_seq_len(rs), ms, (size_t)ll_seq_len(ms), d_expr, &hr, &hm);
    expr_free(hr);
    expr_free(hm);
    return r;
}

Expr* builtin_chineseremainder(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    /* The message says "2 arguments are expected" even though 3 are accepted:
     * argr for one argument, argrx otherwise. */
    if (argc < 2 || argc > 3) {
        if (argc == 1)
            mth_message_gated(g_arith_warnings_muted, "ChineseRemainder", "argr",
                "ChineseRemainder called with 1 argument; 2 arguments are expected.");
        else
            mth_message_gated(g_arith_warnings_muted, "ChineseRemainder", "argrx",
                "ChineseRemainder called with %zu arguments; 2 arguments are expected.", argc);
        return NULL;
    }

    Expr* rlist = res->data.function.args[0];
    Expr* mlist = res->data.function.args[1];

    /* Optional offset d must be an exact integer (ChineseRemainder::int). */
    Expr* d_expr = (argc == 3) ? res->data.function.args[2] : NULL;
    if (d_expr && !expr_is_integer_like(d_expr))
        return nt_emit_int("ChineseRemainder", 3, res);

    /* Symbolic arguments leave the call unevaluated silently, so symbolic flow
     * and rules still apply. A LinkedList is the List of its elements. */
    if (ll_seq_len(rlist) < 0 || ll_seq_len(mlist) < 0) return NULL;
    return cr_core(rlist, mlist, d_expr);
}
