/*
 * fi_trial.c — batched trial division and the shared small-prime table.
 *
 * The naive loop divides the bignum once per prime.  Here the primes are
 * grouped into word-sized products M = p1 p2 ... pk < 2^64 and n is reduced
 * once per group (r = n mod M, a single mpn_mod_1-class pass), after which
 * each p only needs the word operation r mod p.  For primes below 2^16 that
 * is four primes per bignum pass, three up to 2^21, two up to 2^32.  A prime
 * found is divided out with mpz_divexact_ui until it no longer divides.  The
 * scan stops as soon as p^2 exceeds the remaining cofactor, which is then
 * prime.  Primes above the table (2^20) come from a segmented sieve.
 */
#include "fi_internal.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static uint32_t* g_primes = NULL;   /* odd primes < FI_SMALL_PRIME_LIMIT, plus 2 first */
static size_t g_nprimes = 0;

const uint32_t* fi_small_primes(size_t* count) {
    if (!g_primes) {
        uint32_t L = FI_SMALL_PRIME_LIMIT;
        unsigned char* comp = calloc(L / 2 + 1, 1);   /* comp[i]: 2i+1 composite */
        size_t cnt = 1;
        for (uint32_t i = 1; (2 * i + 1) * (2 * i + 1) < L; i++)
            if (!comp[i])
                for (uint32_t j = (2 * i + 1) * (2 * i + 1) / 2; j < L / 2; j += 2 * i + 1)
                    comp[j] = 1;
        for (uint32_t i = 1; i < L / 2; i++) cnt += !comp[i];
        uint32_t* P = malloc(cnt * sizeof(uint32_t));
        size_t k = 0;
        P[k++] = 2;
        for (uint32_t i = 1; i < L / 2; i++) if (!comp[i]) P[k++] = 2 * i + 1;
        free(comp);
        g_nprimes = k;
        g_primes = P;
    }
    if (count) *count = g_nprimes;
    return g_primes;
}

/* Divide p out of n completely; returns the exponent. */
static int64_t strip_prime(mpz_t n, uint64_t p) {
    int64_t e = 0;
    while (mth_mpz_divisible_u64_p(n, (uint64_t)p)) {
        mth_mpz_divexact_u64(n, n, (uint64_t)p);
        e++;
    }
    return e;
}

/* Process one batch of primes ps[0..k) whose product is M.  Returns 1 when
 * the scan can stop (n == 1 or the cofactor is known prime). */
static int batch(mpz_t n, const uint64_t* ps, int k, uint64_t M,
                 int64_t mult, FiFactors* out, size_t* found) {
    uint64_t r = mth_mpz_fdiv_u64(n, (uint64_t)M);
    for (int i = 0; i < k; i++) {
        if (r % ps[i] == 0) {
            int64_t e = strip_prime(n, ps[i]);
            if (e) { fi_factors_add_ui(out, ps[i], e * mult); (*found)++; }
        }
    }
    if (mpz_cmp_ui(n, 1) == 0) return 1;
    /* Cofactor below (largest prime tried)^2 has no factor left: prime. */
    uint64_t pl = ps[k - 1];
    if (mpz_sizeinbase(n, 2) <= 64) {
        uint64_t v = mth_mpz_get_u64(n);
        uint64_t h, l;
        fi_mul64(pl, pl, &h, &l);
        if (!h && v < l) {
            fi_factors_add(out, n, mult);
            (*found)++;
            mpz_set_ui(n, 1);
            return 1;
        }
    }
    return 0;
}

size_t fi_trial_divide(mpz_t n, uint64_t bound, int64_t mult, FiFactors* out) {
    size_t found = 0;
    if (mpz_cmp_ui(n, 1) <= 0) return 0;
    mp_bitcnt_t z = mpz_scan1(n, 0);
    if (z) {
        mpz_tdiv_q_2exp(n, n, z);
        fi_factors_add_ui(out, 2, (int64_t)z * mult);
        found++;
    }
    if (bound < 3 || mpz_cmp_ui(n, 1) == 0) return found;

    size_t np;
    const uint32_t* P = fi_small_primes(&np);
    uint64_t ps[8];
    int k = 0;
    uint64_t M = 1;
    size_t i;
    for (i = 1; i < np && P[i] <= bound; i++) {
        uint64_t p = P[i];
        uint64_t h, l;
        fi_mul64(M, p, &h, &l);
        if (h || k == 8) {
            if (batch(n, ps, k, M, mult, out, &found)) return found;
            k = 0; M = 1; l = p;
        }
        ps[k++] = p;
        M = l;
    }
    if (k && batch(n, ps, k, M, mult, out, &found)) return found;
    if (bound <= FI_SMALL_PRIME_LIMIT || mpz_cmp_ui(n, 1) == 0) return found;

    /* Beyond the table: segmented sieve over (2^20, bound], same batching.
     * The segment primes are sieved with the table primes. */
    const uint64_t SEG = 1u << 18;
    unsigned char* seg = malloc(SEG);
    k = 0; M = 1;
    for (uint64_t lo = FI_SMALL_PRIME_LIMIT; lo <= bound; lo += SEG) {
        uint64_t hi = lo + SEG - 1;
        if (hi > bound) hi = bound;
        memset(seg, 1, (size_t)(hi - lo + 1));
        for (size_t j = 0; j < np; j++) {
            uint64_t p = P[j];
            if (p * p > hi) break;
            uint64_t s = ((lo + p - 1) / p) * p;
            if (s < p * p) s = p * p;
            for (uint64_t x = s; x <= hi; x += p) seg[x - lo] = 0;
        }
        for (uint64_t x = lo | 1; x <= hi; x += 2) {
            if (!seg[x - lo]) continue;
            uint64_t h, l;
            fi_mul64(M, x, &h, &l);
            if (h || k == 8) {
                if (batch(n, ps, k, M, mult, out, &found)) { free(seg); return found; }
                k = 0; M = 1; l = x;
            }
            ps[k++] = x;
            M = l;
        }
    }
    if (k) batch(n, ps, k, M, mult, out, &found);
    free(seg);
    return found;
}
