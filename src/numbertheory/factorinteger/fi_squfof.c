/*
 * fi_squfof.c — Shanks' square forms factorisation beyond one word.
 *
 * SQUFOF walks the continued-fraction expansion of sqrt(kN) looking for a
 * square form; the form coefficients stay below 2 sqrt(kN).  For kN < 2^124
 * they fit a machine word (signed, < 2^63), so the main loop runs in plain
 * word arithmetic even though N itself is two words: only the start and the
 * reverse cycle need kN - P^2, done once each in GMP.  The expected work is
 * O(N^(1/4)) per multiplier, so beyond 2^124 the method is declined.
 *
 * Multipliers: the 16 squarefree products of 3, 5, 7, 11 (Gower–Wagstaff),
 * tried in turn, each with a budget of `iterations` forward steps (default
 * 4 sqrt(2 sqrt(kN)), the bound that makes a single multiplier succeed with
 * high probability).
 */
#include "fi_internal.h"

#include <math.h>

static int is_sq(uint64_t x, uint64_t* r) {
    static const uint64_t QR64 = UINT64_C(0x0202021202030213);   /* squares mod 64 */
    if (!((QR64 >> (x & 63)) & 1)) return 0;
    uint64_t s = fi_u64_isqrt(x);
    if (s * s != x) return 0;
    *r = s;
    return 1;
}

static int squfof_k(mpz_t d, const mpz_t n, unsigned long k, uint64_t budget) {
    mpz_t kn, t;
    mpz_inits(kn, t, NULL);
    mpz_mul_ui(kn, n, k);
    int ok = 0;
    if (mpz_sizeinbase(kn, 2) > 124) goto done;
    mpz_sqrtrem(t, d, kn);                      /* d used as scratch: remainder */
    if (mpz_sgn(d) == 0) {
        mpz_gcd(d, n, t);
        ok = mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, n) < 0;
        goto done;
    }
    uint64_t s = mth_mpz_get_u64(t);
    uint64_t P = s, Qp = 1, Q = mth_mpz_get_u64(d);  /* Q = kN - s^2 */
    if (!budget) budget = 4 * (uint64_t)sqrt(2.0 * sqrt(mpz_get_d(kn))) + 64;
    uint64_t r = 0, i;
    for (i = 1; i < budget; i++) {
        uint64_t b = (s + P) / Q;
        uint64_t Pn = b * Q - P;
        int64_t dq = (int64_t)Qp + (int64_t)b * ((int64_t)P - (int64_t)Pn);
        Qp = Q; Q = (uint64_t)dq; P = Pn;
        /* After step i, Q holds Q_{i+1}: test the even-indexed forms. */
        if ((i & 1) && is_sq(Q, &r)) break;
    }
    if (i >= budget) goto done;
    {
        /* Reverse cycle: Q0' = r, P0' = b r + P, Q1' = (kN - P0'^2) / r. */
        uint64_t b = (s - P) / r;
        P = b * r + P;
        Qp = r;
        mth_mpz_set_u64(t, (uint64_t)P);
        mpz_mul(t, t, t);
        mpz_sub(t, kn, t);
        mth_mpz_divexact_u64(t, t, (uint64_t)r);
        Q = mth_mpz_get_u64(t);
        for (uint64_t j = 0; j < budget; j++) {
            uint64_t bb = (s + P) / Q;
            uint64_t Pn = bb * Q - P;
            if (Pn == P) break;
            int64_t dq = (int64_t)Qp + (int64_t)bb * ((int64_t)P - (int64_t)Pn);
            Qp = Q; Q = (uint64_t)dq; P = Pn;
        }
        mth_mpz_set_u64(t, (uint64_t)P);
        mpz_gcd(d, n, t);
        ok = mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, n) < 0;
    }
done:
    mpz_clears(kn, t, NULL);
    return ok;
}

static const unsigned short MULT[] = {
    1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155
};

int fi_split_squfof(mpz_t d, const mpz_t n, FiCtx* c, uint64_t max_iters) {
    (void)c;
    for (size_t i = 0; i < sizeof MULT / sizeof MULT[0]; i++)
        if (squfof_k(d, n, MULT[i], max_iters)) return 1;
    return 0;
}
