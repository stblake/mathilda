/* primitiveroot.c -- PrimitiveRoot[] / PrimitiveRootList[] and the pr_* cluster.
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
 * PrimitiveRoot / PrimitiveRootList
 *
 *   PrimitiveRoot[n]        a primitive root of n
 *   PrimitiveRoot[n, k]     smallest primitive root of n >= k
 *   PrimitiveRootList[n]    sorted list of primitive roots of n in [1, n-1]
 *
 * A primitive root of n is a generator of the multiplicative group
 * (Z/nZ)^*, which is cyclic iff n in {1, 2, 4} or n = p^k or n = 2 p^k
 * for an odd prime p and k >= 1.  For non-cyclic n the call is left
 * unevaluated (PrimitiveRoot) or returns {} (PrimitiveRootList).
 *
 * All arithmetic is on GMP mpz_t so machine and big integers share the
 * same code path; int64 inputs are coerced via expr_to_mpz.
 * ===================================================================*/

/* Naive prime test for unsigned long; only used to skip composite trial
 * exponents in pr_decompose_odd_prime_power.  q <= log2(n) is tiny in
 * practice (a few hundred at most for any realistic input). */
static bool pr_ulong_is_prime(unsigned long q) {
    if (q < 2) return false;
    if (q < 4) return true;
    if ((q & 1UL) == 0) return false;
    for (unsigned long d = 3; d * d <= q; d += 2) {
        if (q % d == 0) return false;
    }
    return true;
}

/* Decompose an odd integer n >= 3 as p^k with p prime.  Returns true on
 * success; on success p_out is set to p and *k_out to k.  Returns false
 * if n is composite but not a prime power, or if the base turns out to
 * be 2 (n was supposed to be odd).
 *
 * Algorithm: iteratively strip prime exponents q from the current value
 * by computing the exact q-th root.  After O(omega(k)) reductions the
 * current value is the prime base p.  Each mpz_root is O(M(size)) so
 * the whole routine is dominated by the primality test on the final p. */
static bool pr_decompose_odd_prime_power(const mpz_t n, mpz_t p_out, uint64_t* k_out) {
    if (mpz_cmp_ui(n, 1) <= 0) return false;
    mpz_t current;
    mpz_init_set(current, n);
    uint64_t total_k = 1;

    while (mpz_probab_prime_p(current, 25) == 0) {
        unsigned long log2_cur = (unsigned long)mpz_sizeinbase(current, 2);
        bool factored = false;
        for (unsigned long q = 2; q <= log2_cur; q++) {
            if (!pr_ulong_is_prime(q)) continue;
            mpz_t cand;
            mpz_init(cand);
            int exact = mpz_root(cand, current, q);
            if (exact) {
                mpz_set(current, cand);
                total_k *= q;
                mpz_clear(cand);
                factored = true;
                break;
            }
            mpz_clear(cand);
        }
        if (!factored) {
            mpz_clear(current);
            return false;
        }
    }
    /* current is prime.  Reject p == 2 since we only recognise odd prime
     * power moduli here; the 2 / 4 cases are handled at the caller. */
    if (mpz_cmp_ui(current, 2) == 0) {
        mpz_clear(current);
        return false;
    }
    mpz_set(p_out, current);
    *k_out = total_k;
    mpz_clear(current);
    return true;
}

/* Classify n as a cyclic modulus.  On success sets:
 *   *is_two_pk_out  true iff n = 2 p^k (k >= 1, odd prime p)
 *   p_out, *k_out   the odd prime base and exponent (only meaningful for
 *                   the n = p^k and n = 2 p^k cases; for n in {2, 4}
 *                   they are not written).
 * Returns one of CYCLIC_NONE / CYCLIC_TWO / CYCLIC_FOUR / CYCLIC_PRIME_POWER
 * / CYCLIC_TWO_PRIME_POWER.  n must be a positive mpz_t. */
typedef enum {
    PR_CYC_NONE = 0,
    PR_CYC_TWO,            /* n == 2 */
    PR_CYC_FOUR,           /* n == 4 */
    PR_CYC_PRIME_POWER,    /* n == p^k, p odd prime, k >= 1 */
    PR_CYC_TWO_PRIME_POWER /* n == 2 p^k */
} pr_cyclic_kind;

static pr_cyclic_kind pr_classify(const mpz_t n, mpz_t p_out, uint64_t* k_out) {
    if (mpz_cmp_ui(n, 2) == 0) return PR_CYC_TWO;
    if (mpz_cmp_ui(n, 4) == 0) return PR_CYC_FOUR;
    if (mpz_cmp_ui(n, 1) <= 0) return PR_CYC_NONE;

    if (mpz_odd_p(n)) {
        if (pr_decompose_odd_prime_power(n, p_out, k_out)) return PR_CYC_PRIME_POWER;
        return PR_CYC_NONE;
    }
    /* even: n must be 2 * (odd prime power).  Strip exactly one factor of 2. */
    if (mpz_divisible_2exp_p(n, 2)) return PR_CYC_NONE; /* 4 | n => non-cyclic */
    mpz_t half;
    mpz_init(half);
    mpz_tdiv_q_2exp(half, n, 1);
    if (mpz_cmp_ui(half, 1) == 0) {
        /* n == 2, already handled above. */
        mpz_clear(half);
        return PR_CYC_NONE;
    }
    bool ok = pr_decompose_odd_prime_power(half, p_out, k_out);
    mpz_clear(half);
    return ok ? PR_CYC_TWO_PRIME_POWER : PR_CYC_NONE;
}

/* Maximum number of distinct prime factors of phi(n) we expect.  phi for
 * any cyclic modulus is bounded by n, and omega(m) <= ~15 for m <= 10^18,
 * <= ~50 for m <= 10^60, <= ~330 for m <= 10^900.  256 leaves headroom
 * for moderately huge p while staying stack-friendly. */
/* PR_MAX_DISTINCT_PRIMES is defined in numbertheory_internal.h */

/* Collect distinct prime divisors of m (m >= 1) into out[0..*count).
 * out has capacity PR_MAX_DISTINCT_PRIMES; *count must be 0 on entry and
 * each written out[i] is mpz_init_set by the routine.  Returns true on
 * success.  Caller must mpz_clear out[0..*count) on both true and false.
 *
 * Strategy:
 *   - Trial-divide by primes up to PR_TRIAL_LIMIT (covers p - 1 for any
 *     p up to PR_TRIAL_LIMIT^2 ~ 1e10 in one pass).
 *   - If a non-trivial cofactor remains, descend via Pollard rho until
 *     all factors are prime.
 * This avoids pulling in the full FactorInteger pipeline for what is
 * almost always a small number (p - 1). */
#define PR_TRIAL_LIMIT 100000UL

static void pr_pollard_rho(mpz_t f_out, const mpz_t n_in, unsigned long c) {
    mpz_t x, y, d, diff, n;
    mpz_inits(x, y, d, diff, n, NULL);
    mpz_set(n, n_in);
    mpz_set_ui(x, 2);
    mpz_set_ui(y, 2);
    mpz_set_ui(d, 1);
    unsigned long steps = 0;
    while (mpz_cmp_ui(d, 1) == 0 && steps < 1000000UL) {
        /* x = x^2 + c mod n */
        mpz_mul(x, x, x); mpz_add_ui(x, x, c); mpz_mod(x, x, n);
        /* y = (y^2 + c)^2 + c mod n  (advance two steps) */
        mpz_mul(y, y, y); mpz_add_ui(y, y, c); mpz_mod(y, y, n);
        mpz_mul(y, y, y); mpz_add_ui(y, y, c); mpz_mod(y, y, n);
        mpz_sub(diff, x, y); mpz_abs(diff, diff);
        mpz_gcd(d, diff, n);
        steps++;
    }
    if (mpz_cmp(d, n) == 0 || mpz_cmp_ui(d, 1) == 0) {
        mpz_set_ui(f_out, 0); /* failure signal */
    } else {
        mpz_set(f_out, d);
    }
    mpz_clears(x, y, d, diff, n, NULL);
}

static bool pr_add_unique(mpz_t* out, size_t* count, const mpz_t p) {
    for (size_t i = 0; i < *count; i++) {
        if (mpz_cmp(out[i], p) == 0) return true;
    }
    if (*count >= PR_MAX_DISTINCT_PRIMES) return false;
    mpz_init_set(out[*count], p);
    (*count)++;
    return true;
}

static bool pr_distinct_primes_recursive(mpz_t n, mpz_t* out, size_t* count) {
    if (mpz_cmp_ui(n, 1) <= 0) return true;
    if (mpz_probab_prime_p(n, 25)) return pr_add_unique(out, count, n);
    mpz_t f, q;
    mpz_inits(f, q, NULL);
    bool ok = false;
    for (unsigned long c = 1; c < 64 && !ok; c++) {
        pr_pollard_rho(f, n, c);
        if (mpz_sgn(f) != 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, n) < 0) {
            ok = true;
        }
    }
    if (!ok) {
        mpz_clears(f, q, NULL);
        return false;
    }
    mpz_divexact(q, n, f);
    bool result = pr_distinct_primes_recursive(f, out, count) &&
                  pr_distinct_primes_recursive(q, out, count);
    mpz_clears(f, q, NULL);
    return result;
}

/* Append distinct prime factors of m_in to out[*count..]; preserves any
 * entries already present (callers chain this with pr_add_unique). */
bool pr_collect_distinct_primes(const mpz_t m_in, mpz_t* out, size_t* count) {
    if (mpz_cmp_ui(m_in, 1) <= 0) return true;
    mpz_t m;
    mpz_init_set(m, m_in);
    for (unsigned long p = 2; p < PR_TRIAL_LIMIT && mpz_cmp_ui(m, 1) > 0; p++) {
        if (mpz_divisible_ui_p(m, p)) {
            mpz_t pz;
            mpz_init_set_ui(pz, p);
            if (!pr_add_unique(out, count, pz)) {
                mpz_clear(pz);
                mpz_clear(m);
                return false;
            }
            mpz_clear(pz);
            while (mpz_divisible_ui_p(m, p)) mpz_divexact_ui(m, m, p);
        }
    }
    bool ok = true;
    if (mpz_cmp_ui(m, 1) > 0) ok = pr_distinct_primes_recursive(m, out, count);
    mpz_clear(m);
    return ok;
}

/* Build phi(n) and the distinct prime factors of phi(n) for a cyclic
 * modulus described by `kind`, `p`, `k`.  Returns true on success.
 * phi_out is mpz_init'd by the caller; primes[i] are mpz_init_set by
 * this routine and must be cleared by the caller. */
static bool pr_phi_and_primes(pr_cyclic_kind kind, const mpz_t p, uint64_t k,
                              mpz_t phi_out, mpz_t* primes, size_t* nprimes) {
    *nprimes = 0;
    switch (kind) {
        case PR_CYC_TWO:
            mpz_set_ui(phi_out, 1);
            return true;
        case PR_CYC_FOUR: {
            mpz_set_ui(phi_out, 2);
            mpz_t two;
            mpz_init_set_ui(two, 2);
            bool ok = pr_add_unique(primes, nprimes, two);
            mpz_clear(two);
            return ok;
        }
        case PR_CYC_PRIME_POWER:
        case PR_CYC_TWO_PRIME_POWER: {
            /* phi(p^k)   = p^(k-1) (p - 1)
             * phi(2 p^k) = phi(2) phi(p^k) = p^(k-1) (p - 1)  (same value) */
            mpz_t pm1, pkm1;
            mpz_inits(pm1, pkm1, NULL);
            mpz_sub_ui(pm1, p, 1);
            if (k == 1) {
                mpz_set_ui(pkm1, 1);
            } else {
                mpz_pow_ui(pkm1, p, (unsigned long)(k - 1));
            }
            mpz_mul(phi_out, pkm1, pm1);

            /* Distinct prime factors of phi: {p if k >= 2} ∪ primes(p - 1). */
            bool ok = true;
            if (k >= 2) ok = pr_add_unique(primes, nprimes, p);
            if (ok) ok = pr_collect_distinct_primes(pm1, primes, nprimes);
            mpz_clears(pm1, pkm1, NULL);
            return ok;
        }
        default:
            return false;
    }
}

/* True iff g is a primitive root of n, given phi = phi(n) and the
 * distinct prime divisors of phi.  Pre: 1 <= g < n + something
 * (we always reduce modulo n).  Uses mpz_powm for the order test. */
static bool pr_is_primitive_root(const mpz_t g_in, const mpz_t n, const mpz_t phi,
                                 const mpz_t* primes, size_t nprimes) {
    mpz_t g, gcd_gn;
    mpz_inits(g, gcd_gn, NULL);
    mpz_mod(g, g_in, n);
    if (mpz_sgn(g) == 0) {
        mpz_clears(g, gcd_gn, NULL);
        return false;
    }
    mpz_gcd(gcd_gn, g, n);
    if (mpz_cmp_ui(gcd_gn, 1) != 0) {
        mpz_clears(g, gcd_gn, NULL);
        return false;
    }
    /* Phi == 1 (n in {1, 2}): only the residue 1 generates the trivial
     * group, so any unit is a "primitive root".  All units of Z/2Z is just
     * {1}, which we've already accepted by gcd above. */
    if (mpz_cmp_ui(phi, 1) == 0) {
        mpz_clears(g, gcd_gn, NULL);
        return true;
    }
    mpz_t exp, power;
    mpz_inits(exp, power, NULL);
    bool is_pr = true;
    for (size_t i = 0; i < nprimes; i++) {
        mpz_divexact(exp, phi, primes[i]);
        mpz_powm(power, g, exp, n);
        if (mpz_cmp_ui(power, 1) == 0) { is_pr = false; break; }
    }
    mpz_clears(g, gcd_gn, exp, power, NULL);
    return is_pr;
}

/* Find the smallest primitive root of n that is >= start.  On success
 * writes the root to g_out and returns true.  Caller pre-inits g_out.
 *
 * The density of primitive roots in [1, n-1] is phi(phi(n)) / phi(n),
 * which for typical p is on the order of 1 / ln(ln(p)), so the search
 * terminates after O(ln(ln(p))) candidates on average. */
static bool pr_smallest_primitive_root(mpz_t g_out, const mpz_t n, const mpz_t phi,
                                       const mpz_t* primes, size_t nprimes,
                                       const mpz_t start) {
    mpz_t g;
    mpz_init(g);
    if (mpz_sgn(start) < 1) mpz_set_ui(g, 1); else mpz_set(g, start);

    /* For tiny n (2, 4) the loop terminates immediately. */
    while (mpz_cmp(g, n) < 0 || mpz_cmp_ui(n, 1) <= 0) {
        if (pr_is_primitive_root(g, n, phi, primes, nprimes)) {
            mpz_set(g_out, g);
            mpz_clear(g);
            return true;
        }
        mpz_add_ui(g, g, 1);
        /* Defensive cap: a primitive root always exists in [1, n-1] for
         * cyclic n, so we should never iterate beyond n.  The break below
         * catches the n == 2 case (where the only PR is 1 and start >= 2
         * would otherwise loop forever). */
        if (mpz_cmp(g, n) >= 0 && mpz_cmp(start, n) >= 0) break;
    }
    mpz_clear(g);
    return false;
}

/* Emit `PrimitiveRoot::argt: PrimitiveRoot called with N arguments; 1 or
 * 2 arguments are expected.` to stderr and return NULL. */
static Expr* pr_emit_argt(size_t argc) {
    fprintf(stderr,
            "PrimitiveRoot::argt: PrimitiveRoot called with %zu argument%s; "
            "1 or 2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

/* Emit `PrimitiveRoot::intg: Integer greater than 1 expected at position
 * <pos> in <call>.` */
static Expr* pr_emit_intg(size_t pos, Expr* res) {
    char* call_str = expr_to_string(res);
    fprintf(stderr,
            "PrimitiveRoot::intg: Integer greater than 1 expected at position "
            "%zu in %s.\n",
            pos, call_str ? call_str : "?");
    free(call_str);
    return NULL;
}

/* Emit `PrimitiveRootList::argx: PrimitiveRootList called with N arguments;
 * 1 argument is expected.` */
static Expr* prl_emit_argx(size_t argc) {
    fprintf(stderr,
            "PrimitiveRootList::argx: PrimitiveRootList called with %zu argument%s; "
            "1 argument is expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_primitiveroot(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return pr_emit_argt(argc);

    Expr* n_expr = res->data.function.args[0];
    if (!expr_is_integer_like(n_expr)) {
        /* Numeric non-integer (Real, Rational, Complex, ...) triggers the
         * Mathematica-style ::intg diagnostic; symbolic args flow through
         * silently so user-supplied DownValues / pattern matching apply. */
        if (expr_is_numeric_like(n_expr)) return pr_emit_intg(1, res);
        return NULL;
    }

    mpz_t n;
    expr_to_mpz(n_expr, n);
    if (mpz_cmp_ui(n, 2) < 0) {
        mpz_clear(n);
        return pr_emit_intg(1, res);
    }

    mpz_t p; uint64_t k;
    mpz_init(p);
    pr_cyclic_kind kind = pr_classify(n, p, &k);
    if (kind == PR_CYC_NONE) {
        mpz_clears(n, p, NULL);
        return NULL;
    }

    mpz_t phi;
    mpz_init(phi);
    mpz_t primes[PR_MAX_DISTINCT_PRIMES];
    size_t nprimes = 0;
    if (!pr_phi_and_primes(kind, p, k, phi, primes, &nprimes)) {
        for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
        mpz_clears(n, p, phi, NULL);
        return NULL;
    }

    /* For 2-arg PrimitiveRoot[n, k]: forward search from k.
     * For 1-arg PrimitiveRoot[n]:
     *   - n = 2          -> 1
     *   - n = 4          -> 3
     *   - n = p^k        -> smallest primitive root of n
     *   - n = 2 p^k      -> g if odd else g + p^k, where g is the
     *                       smallest primitive root of p^k
     */
    Expr* out = NULL;
    if (argc == 2) {
        Expr* k_expr = res->data.function.args[1];
        if (!expr_is_integer_like(k_expr)) {
            if (expr_is_numeric_like(k_expr)) {
                for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
                mpz_clears(n, p, phi, NULL);
                return pr_emit_intg(2, res);
            }
            for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
            mpz_clears(n, p, phi, NULL);
            return NULL;
        }
        mpz_t start, g;
        mpz_inits(start, g, NULL);
        expr_to_mpz(k_expr, start);
        if (mpz_sgn(start) < 0) mpz_set_ui(start, 0);
        bool found = pr_smallest_primitive_root(g, n, phi, primes, nprimes, start);
        if (found) out = expr_bigint_normalize(expr_new_bigint_from_mpz(g));
        mpz_clears(start, g, NULL);
    } else {
        switch (kind) {
            case PR_CYC_TWO:
                out = expr_new_integer(1);
                break;
            case PR_CYC_FOUR:
                out = expr_new_integer(3);
                break;
            case PR_CYC_PRIME_POWER: {
                mpz_t g, two;
                mpz_inits(g, two, NULL);
                mpz_set_ui(two, 2);
                bool found = pr_smallest_primitive_root(g, n, phi, primes, nprimes, two);
                if (found) out = expr_bigint_normalize(expr_new_bigint_from_mpz(g));
                mpz_clears(g, two, NULL);
                break;
            }
            case PR_CYC_TWO_PRIME_POWER: {
                /* Smallest PR of p^k.  If odd, return as-is; else add p^k. */
                mpz_t pk, g, two;
                mpz_inits(pk, g, two, NULL);
                mpz_pow_ui(pk, p, (unsigned long)k);
                mpz_set_ui(two, 2);
                bool found = pr_smallest_primitive_root(g, pk, phi, primes, nprimes, two);
                if (found) {
                    if (mpz_even_p(g)) mpz_add(g, g, pk);
                    out = expr_bigint_normalize(expr_new_bigint_from_mpz(g));
                }
                mpz_clears(pk, g, two, NULL);
                break;
            }
            default:
                break;
        }
    }

    for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
    mpz_clears(n, p, phi, NULL);
    return out;
}

/* Comparator for sorting an mpz_t* array via qsort. */
static int pr_mpz_cmp_qsort(const void* a, const void* b) {
    const mpz_t* ma = (const mpz_t*)a;
    const mpz_t* mb = (const mpz_t*)b;
    return mpz_cmp(*ma, *mb);
}

Expr* builtin_primitiverootlist(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return prl_emit_argx(argc);

    Expr* n_expr = res->data.function.args[0];
    if (!expr_is_integer_like(n_expr)) {
        /* Per spec: PrimitiveRootList[non-integer] stays unevaluated with
         * no diagnostic (Real / Complex / Symbolic all flow through). */
        return NULL;
    }

    mpz_t n;
    expr_to_mpz(n_expr, n);
    if (mpz_cmp_ui(n, 1) <= 0) {
        mpz_clear(n);
        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    }

    mpz_t p; uint64_t k;
    mpz_init(p);
    pr_cyclic_kind kind = pr_classify(n, p, &k);
    if (kind == PR_CYC_NONE) {
        mpz_clears(n, p, NULL);
        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    }

    mpz_t phi;
    mpz_init(phi);
    mpz_t primes[PR_MAX_DISTINCT_PRIMES];
    size_t nprimes = 0;
    if (!pr_phi_and_primes(kind, p, k, phi, primes, &nprimes)) {
        for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
        mpz_clears(n, p, phi, NULL);
        return NULL;
    }

    /* The enumeration is g^i mod n for i in [1, phi] coprime to phi;
     * the number of such i is phi(phi(n)).  Guard against pathological
     * sizes: if phi doesn't fit in an unsigned long we cannot index the
     * loop, and even when it does the resulting list could be huge. */
    if (!mpz_fits_ulong_p(phi)) {
        for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
        mpz_clears(n, p, phi, NULL);
        return NULL;
    }
    unsigned long phi_ul = mpz_get_ui(phi);

    /* Smallest primitive root of n itself (not of p^k) — for the n = 2 p^k
     * case we want the residue in [1, n-1], which lets us walk g^i mod n
     * directly. */
    mpz_t g, two;
    mpz_inits(g, two, NULL);
    mpz_set_ui(two, 1);
    bool have_g = pr_smallest_primitive_root(g, n, phi, primes, nprimes, two);
    if (!have_g) {
        mpz_clears(g, two, NULL);
        for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
        mpz_clears(n, p, phi, NULL);
        return NULL;
    }

    /* Walk g^i mod n iteratively; collect i with gcd(i, phi) == 1. */
    mpz_t cur, i_mpz, g_local;
    mpz_inits(cur, i_mpz, g_local, NULL);
    mpz_set_ui(cur, 1);
    mpz_set(g_local, g);

    mpz_t* roots = (mpz_t*)malloc(sizeof(mpz_t) * (phi_ul + 1));
    if (!roots) {
        mpz_clears(g, two, cur, i_mpz, g_local, NULL);
        for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
        mpz_clears(n, p, phi, NULL);
        return NULL;
    }
    size_t root_count = 0;

    for (unsigned long i = 1; i <= phi_ul; i++) {
        mpz_mul(cur, cur, g_local);
        mpz_mod(cur, cur, n);
        unsigned long ggcd;
        if (phi_ul == 0) {
            ggcd = i;
        } else {
            /* gcd(i, phi_ul) in unsigned long. */
            unsigned long a = i, b = phi_ul;
            while (b) { unsigned long t = a % b; a = b; b = t; }
            ggcd = a;
        }
        if (ggcd == 1) {
            mpz_init_set(roots[root_count], cur);
            root_count++;
        }
    }

    qsort(roots, root_count, sizeof(mpz_t), pr_mpz_cmp_qsort);

    Expr** items = NULL;
    if (root_count > 0) {
        items = (Expr**)malloc(sizeof(Expr*) * root_count);
        for (size_t i = 0; i < root_count; i++) {
            items[i] = expr_bigint_normalize(expr_new_bigint_from_mpz(roots[i]));
        }
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, root_count);
    free(items);

    for (size_t i = 0; i < root_count; i++) mpz_clear(roots[i]);
    free(roots);

    mpz_clears(g, two, cur, i_mpz, g_local, NULL);
    for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
    mpz_clears(n, p, phi, NULL);
    return out;
}
