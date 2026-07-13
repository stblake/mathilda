/* multiplicativeorder.c -- MultiplicativeOrder[].
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

/* =====================================================================
 * MultiplicativeOrder
 *
 *   MultiplicativeOrder[k, n]
 *       smallest positive m with k^m == 1 (mod n).
 *   MultiplicativeOrder[k, n, {r1, r2, ...}]
 *       smallest positive m with k^m == r_i (mod n) for some i.
 *
 * Returns unevaluated when:
 *   - n is 0 (no group);
 *   - gcd(k, n) != 1 (k is not a unit mod n, no finite order exists);
 *   - in the 3-arg form, no power of k lands in the residue set, or the
 *     order is too large to iterate over.
 * All arithmetic is on GMP mpz_t so bignum k and n share the same code
 * path with machine integers.  Negative n is treated as |n|; negative or
 * out-of-range k is reduced modulo n.
 * ===================================================================*/

/* phi(n) for n >= 1.  Uses the pr_collect_distinct_primes helper above to
 * gather the distinct prime factors of n, then applies the multiplicative
 * identity phi(n) = n * prod_{p | n} (1 - 1/p), iteratively realised as
 * (phi / p) * (p - 1).  This avoids needing the exact multiplicities. */
static bool mo_eulerphi_mpz(const mpz_t n, mpz_t phi_out) {
    if (mpz_sgn(n) <= 0) return false;
    if (mpz_cmp_ui(n, 1) == 0) { mpz_set_ui(phi_out, 1); return true; }
    mpz_t primes[PR_MAX_DISTINCT_PRIMES];
    size_t nprimes = 0;
    if (!pr_collect_distinct_primes(n, primes, &nprimes)) {
        for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
        return false;
    }
    mpz_set(phi_out, n);
    mpz_t pm1;
    mpz_init(pm1);
    for (size_t i = 0; i < nprimes; i++) {
        mpz_divexact(phi_out, phi_out, primes[i]);
        mpz_sub_ui(pm1, primes[i], 1);
        mpz_mul(phi_out, phi_out, pm1);
    }
    mpz_clear(pm1);
    for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
    return true;
}

/* Multiplicative order of k_in modulo n_in.  Pre: order_out is mpz_init'd
 * by the caller; n_in / k_in are read-only.  Returns true on success.
 * Returns false (without setting order_out) when n is zero, k is not a
 * unit modulo n, or phi(n) cannot be factored.
 *
 * Strategy: any element order divides phi(n).  Compute phi, then for each
 * distinct prime p | phi successively divide phi by p as long as the
 * residual exponent still maps k to 1.  After the loop, the residual is
 * exactly the order. */
static bool mo_order_mpz(const mpz_t k_in, const mpz_t n_in, mpz_t order_out) {
    if (mpz_sgn(n_in) == 0) return false;
    mpz_t n;
    mpz_init_set(n, n_in);
    mpz_abs(n, n);                               /* allow negative n */
    if (mpz_cmp_ui(n, 1) == 0) { mpz_set_ui(order_out, 1); mpz_clear(n); return true; }

    mpz_t k, g;
    mpz_inits(k, g, NULL);
    mpz_mod(k, k_in, n);                         /* k in [0, n-1] */
    mpz_gcd(g, k, n);
    if (mpz_cmp_ui(g, 1) != 0) {
        mpz_clears(k, g, n, NULL);
        return false;                            /* not coprime */
    }
    mpz_clear(g);

    mpz_t phi;
    mpz_init(phi);
    if (!mo_eulerphi_mpz(n, phi)) {
        mpz_clears(k, phi, n, NULL);
        return false;
    }

    /* phi == 1 means the unit group is trivial; the only unit is 1 and
     * its order is 1. */
    if (mpz_cmp_ui(phi, 1) == 0) {
        mpz_set_ui(order_out, 1);
        mpz_clears(k, phi, n, NULL);
        return true;
    }

    mpz_t primes[PR_MAX_DISTINCT_PRIMES];
    size_t nprimes = 0;
    if (!pr_collect_distinct_primes(phi, primes, &nprimes)) {
        for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
        mpz_clears(k, phi, n, NULL);
        return false;
    }

    mpz_t order, q, pow;
    mpz_inits(order, q, pow, NULL);
    mpz_set(order, phi);
    for (size_t i = 0; i < nprimes; i++) {
        while (mpz_divisible_p(order, primes[i])) {
            mpz_divexact(q, order, primes[i]);
            mpz_powm(pow, k, q, n);
            if (mpz_cmp_ui(pow, 1) == 0) {
                mpz_set(order, q);
            } else {
                break;
            }
        }
    }

    mpz_set(order_out, order);
    mpz_clears(order, q, pow, NULL);
    for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
    mpz_clears(k, phi, n, NULL);
    return true;
}

/* Search for the smallest m in [1, order] such that k^m mod n equals one
 * of the residues.  Pre: m_out is mpz_init'd by caller; residues are all
 * already reduced modulo n.  Returns true on success.
 *
 * Cap: bails out if order doesn't fit in unsigned long, or exceeds the
 * MO_SEARCH_CAP iteration budget — for huge groups the enumeration would
 * not terminate in a useful timeframe.  Callers see this as "no match". */
#define MO_SEARCH_CAP 100000000UL
static bool mo_search_residues(const mpz_t k_in, const mpz_t n,
                               const mpz_t order,
                               const mpz_t* residues, size_t nres,
                               mpz_t m_out) {
    if (nres == 0) return false;
    if (!mpz_fits_ulong_p(order)) return false;
    unsigned long d = mpz_get_ui(order);
    if (d > MO_SEARCH_CAP) return false;

    mpz_t cur, k;
    mpz_inits(cur, k, NULL);
    mpz_mod(k, k_in, n);
    mpz_set_ui(cur, 1);
    for (unsigned long m = 1; m <= d; m++) {
        mpz_mul(cur, cur, k);
        mpz_mod(cur, cur, n);
        for (size_t i = 0; i < nres; i++) {
            if (mpz_cmp(cur, residues[i]) == 0) {
                mpz_set_ui(m_out, m);
                mpz_clears(cur, k, NULL);
                return true;
            }
        }
    }
    mpz_clears(cur, k, NULL);
    return false;
}

static Expr* mo_emit_argt(size_t argc) {
    fprintf(stderr,
            "MultiplicativeOrder::argt: MultiplicativeOrder called with %zu "
            "argument%s; 2 or 3 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_multiplicativeorder(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return mo_emit_argt(argc);

    Expr* k_expr = res->data.function.args[0];
    Expr* n_expr = res->data.function.args[1];

    /* Integer-only contract: non-integer numerics (Real, Complex, Rational)
     * and symbolic args flow through unevaluated with no diagnostic, matching
     * MultiplicativeOrder[10., 21] -> MultiplicativeOrder[10., 21]. */
    if (!expr_is_integer_like(k_expr) || !expr_is_integer_like(n_expr)) return NULL;

    mpz_t k, n;
    expr_to_mpz(k_expr, k);
    expr_to_mpz(n_expr, n);
    if (mpz_sgn(n) == 0) { mpz_clears(k, n, NULL); return NULL; }

    mpz_t order;
    mpz_init(order);
    if (!mo_order_mpz(k, n, order)) {
        mpz_clears(k, n, order, NULL);
        return NULL;                             /* not coprime / failure */
    }

    if (argc == 2) {
        Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(order));
        mpz_clears(k, n, order, NULL);
        return out;
    }

    /* 3-arg form: third argument must be a List of integer residues.  A
     * non-List or any non-integer element leaves the call unevaluated. */
    Expr* list = res->data.function.args[2];
    if (list->type != EXPR_FUNCTION ||
        list->data.function.head->type != EXPR_SYMBOL ||
        list->data.function.head->data.symbol.name != SYM_List) {
        mpz_clears(k, n, order, NULL);
        return NULL;
    }
    size_t lcount = list->data.function.arg_count;
    if (lcount == 0) {
        mpz_clears(k, n, order, NULL);
        return NULL;                             /* empty target list */
    }

    mpz_t abs_n;
    mpz_init_set(abs_n, n);
    mpz_abs(abs_n, abs_n);

    mpz_t* residues = (mpz_t*)malloc(sizeof(mpz_t) * lcount);
    if (!residues) {
        mpz_clears(k, n, order, abs_n, NULL);
        return NULL;
    }
    size_t nres = 0;
    bool list_ok = true;
    for (size_t i = 0; i < lcount; i++) {
        Expr* r = list->data.function.args[i];
        if (!expr_is_integer_like(r)) { list_ok = false; break; }
        mpz_init(residues[nres]);
        expr_to_mpz(r, residues[nres]);
        mpz_mod(residues[nres], residues[nres], abs_n);
        nres++;
    }
    if (!list_ok) {
        for (size_t j = 0; j < nres; j++) mpz_clear(residues[j]);
        free(residues);
        mpz_clears(k, n, order, abs_n, NULL);
        return NULL;
    }

    mpz_t m;
    mpz_init(m);
    bool found = mo_search_residues(k, abs_n, order,
                                    (const mpz_t*)residues, nres, m);
    Expr* out = NULL;
    if (found) out = expr_bigint_normalize(expr_new_bigint_from_mpz(m));

    for (size_t i = 0; i < nres; i++) mpz_clear(residues[i]);
    free(residues);
    mpz_clears(k, n, order, abs_n, m, NULL);
    return out;
}
