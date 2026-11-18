/*
 * fi_ecm.c — the elliptic-curve method, Pollard p-1 and Williams p+1.
 *
 * Backends, best first:
 *   GMP-ECM   (default build)    ECM, p-1, p+1 — asymptotically fast stage 2,
 *                                batch stage 1, base-2 special forms;
 *   FLINT     (USE_ECM=0)        ECM (fmpz_factor_ecm), p+1 (fmpz_factor_pp1),
 *                                p-1 native;
 *   native    (neither)          all three on GMP mpz: Montgomery curves with
 *                                Suyama's parametrisation, prime-power ladder
 *                                stage 1 and baby-step/giant-step stage 2;
 *                                p-1 / p+1 with a prime-gap stage 2.
 *
 * Determinism: every random choice (ECM sigma, p+1 seed) is drawn from the
 * FiCtx generator, seeded from n, and GMP-ECM's own RNG is reseeded from it,
 * so a given input always runs the same curves.
 *
 * Effort levels (Automatic's ECM tiers and the ECM method's default) follow
 * GMP-ECM's table of optimal B1 and expected curve counts for finding a
 * factor of D digits with probability 1 - 1/e.
 */
#include "fi_internal.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef NO_ECM
#include "ecm.h"
#endif
#ifdef USE_FLINT
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>
#endif

int fi_have_gmp_ecm(void) {
#ifndef NO_ECM
    return 1;
#else
    return 0;
#endif
}

/* ------------------------------------------------------------------ */
/* Effort table                                                        */
/* ------------------------------------------------------------------ */

static const struct { int digits; double B1; uint64_t curves; } LEVELS[] = {
    { 15, 2000.0,     25 },
    { 20, 11000.0,    74 },
    { 25, 50000.0,    214 },
    { 30, 250000.0,   430 },
    { 35, 1000000.0,  904 },
    { 40, 3000000.0,  2350 },
    { 45, 11000000.0, 4480 },
    { 50, 43000000.0, 7553 },
    { 55, 110000000.0, 17769 },
    { 60, 260000000.0, 42017 },
    { 65, 850000000.0, 69408 },
};
#define NLEVELS (sizeof LEVELS / sizeof LEVELS[0])

void fi_ecm_level_params(int digits, double* B1, uint64_t* curves) {
    size_t i = 0;
    while (i + 1 < NLEVELS && LEVELS[i].digits < digits) i++;
    *B1 = LEVELS[i].B1;
    *curves = LEVELS[i].curves;
}

/* Curves for an explicit B1: the table's count for the nearest level. */
static uint64_t curves_for_B1(double B1) {
    size_t i = 0;
    while (i + 1 < NLEVELS && LEVELS[i].B1 < B1) i++;
    return LEVELS[i].curves;
}

static int nontrivial(const mpz_t d, const mpz_t n) {
    return mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, n) < 0;
}

/* Primes below 2^10 are split off before any curve or p-+1 run: GMP-ECM
 * cannot invert a starting value that shares a factor with n (it reports
 * "Failed to invert" and then crashes), and a curve run on a number with
 * several tiny primes catches them all at once (gcd = n) and fails. */
static int tiny_factor(mpz_t d, const mpz_t n) {
    size_t np;
    const uint32_t* P = fi_small_primes(&np);
    for (size_t i = 0; i < np && P[i] < 1024; i++)
        if (mpz_divisible_ui_p(n, P[i]) && mpz_cmp_ui(n, P[i]) != 0) { mpz_set_ui(d, P[i]); return 1; }
    return 0;
}

/* ------------------------------------------------------------------ */
/* Native ECM (fallback)                                               */
/* ------------------------------------------------------------------ */

#if defined(NO_ECM) && !defined(USE_FLINT)
#define NEED_NATIVE_ECM 1
#endif
#if defined(NO_ECM)
#define NEED_NATIVE_PM1 1
#endif
#if defined(NO_ECM) && !defined(USE_FLINT)
#define NEED_NATIVE_PP1 1
#endif

#ifdef NEED_NATIVE_ECM
typedef struct { mpz_t X, Z; } NPt;

typedef struct {
    mpz_srcptr n;
    mpz_t a24, t1, t2, t3, t4;
} NCurve;

static void nmod(mpz_t r, const NCurve* E) { mpz_mod(r, r, E->n); }

static void n_dbl(NPt* R, const NPt* P, NCurve* E) {
    mpz_add(E->t1, P->X, P->Z); mpz_mul(E->t1, E->t1, E->t1); nmod(E->t1, E);   /* (X+Z)^2 */
    mpz_sub(E->t2, P->X, P->Z); mpz_mul(E->t2, E->t2, E->t2); nmod(E->t2, E);   /* (X-Z)^2 */
    mpz_sub(E->t3, E->t1, E->t2);                                               /* 4XZ     */
    mpz_mul(R->X, E->t1, E->t2); nmod(R->X, E);
    mpz_mul(E->t4, E->a24, E->t3); mpz_add(E->t4, E->t4, E->t2);
    mpz_mul(R->Z, E->t3, E->t4); nmod(R->Z, E);
}

/* R = P + Q, D = P - Q.  R may alias P or Q. */
static void n_add(NPt* R, const NPt* P, const NPt* Q, const NPt* D, NCurve* E) {
    mpz_sub(E->t1, P->X, P->Z); mpz_add(E->t2, Q->X, Q->Z); mpz_mul(E->t1, E->t1, E->t2);
    mpz_add(E->t3, P->X, P->Z); mpz_sub(E->t4, Q->X, Q->Z); mpz_mul(E->t3, E->t3, E->t4);
    mpz_add(E->t2, E->t1, E->t3); mpz_mul(E->t2, E->t2, E->t2); nmod(E->t2, E);
    mpz_sub(E->t4, E->t1, E->t3); mpz_mul(E->t4, E->t4, E->t4); nmod(E->t4, E);
    mpz_mul(E->t2, E->t2, D->Z);
    mpz_mul(E->t4, E->t4, D->X);
    mpz_mod(R->X, E->t2, E->n);
    mpz_mod(R->Z, E->t4, E->n);
}

static void n_mul(NPt* P, uint64_t k, NCurve* E) {
    if (k <= 1) return;
    NPt R0, R1, B;
    mpz_init_set(R0.X, P->X); mpz_init_set(R0.Z, P->Z);
    mpz_init_set(B.X, P->X);  mpz_init_set(B.Z, P->Z);
    mpz_init(R1.X); mpz_init(R1.Z);
    n_dbl(&R1, P, E);
    int top = 63;
    while (!((k >> top) & 1)) top--;
    for (int b = top - 1; b >= 0; b--) {
        if ((k >> b) & 1) { n_add(&R0, &R1, &R0, &B, E); n_dbl(&R1, &R1, E); }
        else              { n_add(&R1, &R1, &R0, &B, E); n_dbl(&R0, &R0, E); }
    }
    mpz_set(P->X, R0.X); mpz_set(P->Z, R0.Z);
    mpz_clears(R0.X, R0.Z, R1.X, R1.Z, B.X, B.Z, NULL);
}

static int native_ecm_curve(mpz_t d, const mpz_t n, uint64_t sigma, double B1d) {
    uint64_t B1 = (uint64_t)B1d;
    uint64_t B2 = 50 * B1;
    /* The fallback walks the shared prime table, so both bounds stop at it. */
    if (B2 > FI_SMALL_PRIME_LIMIT - 512) B2 = FI_SMALL_PRIME_LIMIT - 512;
    if (B1 > B2) B1 = B2;
    NCurve E;
    E.n = n;
    mpz_inits(E.a24, E.t1, E.t2, E.t3, E.t4, NULL);
    mpz_t u, v, x, z, num, den, g;
    mpz_inits(u, v, x, z, num, den, g, NULL);
    /* Suyama. */
    mpz_set_ui(u, (unsigned long)sigma); mpz_mul(u, u, u); mpz_sub_ui(u, u, 5); mpz_mod(u, u, n);
    mpz_set_ui(v, (unsigned long)(4 * sigma)); mpz_mod(v, v, n);
    mpz_powm_ui(x, u, 3, n);
    mpz_powm_ui(z, v, 3, n);
    mpz_sub(num, v, u); mpz_powm_ui(num, num, 3, n);
    mpz_mul_ui(g, u, 3); mpz_add(g, g, v); mpz_mul(num, num, g); mpz_mod(num, num, n);
    mpz_mul_ui(den, x, 16); mpz_mul(den, den, v); mpz_mod(den, den, n);
    int found = 0;
    if (!mpz_invert(g, den, n)) {
        mpz_gcd(g, den, n);
        if (nontrivial(g, n)) { mpz_set(d, g); found = 1; }
        goto out;
    }
    mpz_mul(E.a24, num, g); mpz_mod(E.a24, E.a24, n);
    {
        NPt P;
        mpz_init_set(P.X, x); mpz_init_set(P.Z, z);
        size_t np;
        const uint32_t* PR = fi_small_primes(&np);
        uint64_t acc = 1;
        for (size_t i = 0; i < np && PR[i] <= B1; i++) {
            uint64_t q = PR[i];
            while (q * PR[i] <= B1) q *= PR[i];
            if (acc > UINT64_MAX / q) { n_mul(&P, acc, &E); acc = q; }
            else acc *= q;
        }
        n_mul(&P, acc, &E);
        mpz_gcd(g, P.Z, n);
        if (nontrivial(g, n)) { mpz_set(d, g); found = 1; }
        else if (mpz_cmp_ui(g, 1) == 0) {
            /* Stage 2, baby-step/giant-step: every prime q in (B1, B2] is
             * m D +- j with j < D/2 coprime to D; X_M Z_j - X_j Z_M vanishes
             * mod p exactly when q P is the identity mod p. */
            enum { D = 210, NJ = 24 };
            static const unsigned char J[NJ] = { 1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
                                                 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103 };
            NPt baby[NJ], P2, prev, cur, nx, G, Mc, Mp, Mn;
            mpz_init(P2.X); mpz_init(P2.Z); mpz_init(nx.X); mpz_init(nx.Z);
            mpz_init_set(prev.X, P.X); mpz_init_set(prev.Z, P.Z);
            mpz_init(cur.X); mpz_init(cur.Z);
            n_dbl(&P2, &P, &E);
            n_add(&cur, &P2, &P, &P, &E);                     /* 3P */
            size_t bi = 0;
            mpz_init_set(baby[0].X, P.X); mpz_init_set(baby[0].Z, P.Z);
            bi = 1;
            for (unsigned j = 3; j < D / 2 && bi < NJ; j += 2) {
                if (j == J[bi]) { mpz_init_set(baby[bi].X, cur.X); mpz_init_set(baby[bi].Z, cur.Z); bi++; }
                n_add(&nx, &cur, &P2, &prev, &E);              /* (j+2)P */
                mpz_set(prev.X, cur.X); mpz_set(prev.Z, cur.Z);
                mpz_set(cur.X, nx.X); mpz_set(cur.Z, nx.Z);
            }
            mpz_init_set(G.X, P.X); mpz_init_set(G.Z, P.Z);
            n_mul(&G, D, &E);
            uint64_t m0 = (B1 + D / 2) / D;
            if (m0 < 1) m0 = 1;
            mpz_init_set(Mc.X, G.X); mpz_init_set(Mc.Z, G.Z);
            n_mul(&Mc, m0, &E);
            mpz_init_set(Mp.X, G.X); mpz_init_set(Mp.Z, G.Z);
            if (m0 > 1) n_mul(&Mp, m0 - 1, &E);
            mpz_init(Mn.X); mpz_init(Mn.Z);
            /* Prime bitmap over (B1, B2]. */
            unsigned char* isp = calloc((size_t)(B2 + D + 1), 1);
            for (size_t i = 0; i < np && PR[i] <= B2 + D; i++) isp[PR[i]] = 1;
            mpz_t acc2, t;
            mpz_init_set_ui(acc2, 1);
            mpz_init(t);
            for (uint64_t m = m0; m * D <= B2 + D / 2; m++) {
                for (int k = 0; k < NJ; k++) {
                    uint64_t q1 = m * D - J[k], q2 = m * D + J[k];
                    int hit = (q1 > B1 && q1 <= B2 && isp[q1]) || (q2 > B1 && q2 <= B2 && isp[q2]);
                    if (!hit) continue;
                    mpz_mul(t, Mc.X, baby[k].Z);
                    mpz_submul(t, baby[k].X, Mc.Z);
                    mpz_mul(acc2, acc2, t);
                    mpz_mod(acc2, acc2, n);
                }
                if (m == 1) n_dbl(&Mn, &Mc, &E);
                else n_add(&Mn, &Mc, &G, &Mp, &E);
                mpz_set(Mp.X, Mc.X); mpz_set(Mp.Z, Mc.Z);
                mpz_set(Mc.X, Mn.X); mpz_set(Mc.Z, Mn.Z);
            }
            free(isp);
            mpz_gcd(g, acc2, n);
            if (nontrivial(g, n)) { mpz_set(d, g); found = 1; }
            mpz_clears(acc2, t, NULL);
            for (size_t k = 0; k < bi; k++) mpz_clears(baby[k].X, baby[k].Z, NULL);
            mpz_clears(P2.X, P2.Z, nx.X, nx.Z, prev.X, prev.Z, cur.X, cur.Z, G.X, G.Z,
                       Mc.X, Mc.Z, Mp.X, Mp.Z, Mn.X, Mn.Z, NULL);
        }
        mpz_clears(P.X, P.Z, NULL);
    }
out:
    mpz_clears(u, v, x, z, num, den, g, NULL);
    mpz_clears(E.a24, E.t1, E.t2, E.t3, E.t4, NULL);
    return found;
}
#endif /* NEED_NATIVE_ECM */

/* ------------------------------------------------------------------ */
/* ECM front end                                                       */
/* ------------------------------------------------------------------ */


#ifndef NO_ECM
/* One GMP-ECM curve, run on a worker thread (GMP-ECM is thread-safe since
 * 7.0).  On gcd = n (every prime caught at once) the curve is replayed with
 * B1 shrunk until the primes separate. */
typedef struct {
    mpz_srcptr    n;
    double        B1, B2;
    int           param;
    unsigned long sigma;
    mpz_t         f;
    int           ok;
} CurveJob;

static void curve_job(void* arg) {
    CurveJob* j = (CurveJob*)arg;
    j->ok = 0;
    for (double b = j->B1; b >= 50 && !fi_stop_requested(); b /= 8) {
        ecm_params p;
        ecm_init(p);
        p->stop_asap = fi_stop_requested;   /* TimeConstrained deadline */
        p->param = j->param;
        mpz_set_ui(p->sigma, j->sigma);
        if (j->B2 > 0 && b == j->B1) mpz_set_d(p->B2, j->B2);
        int r = ecm_factor(j->f, (mpz_ptr)j->n, b, p);
        ecm_clear(p);
        if (r > 0 && nontrivial(j->f, j->n)) { j->ok = 1; return; }
        if (!(r > 0 && mpz_cmp(j->f, j->n) == 0)) return;
    }
}

/* Curves below this size are too cheap to be worth a thread each. */
#define ECM_PAR_BITS 160
#endif

/* Run `curves` curves at B1.  With GMP-ECM and threads, the curves go out in
 * regions of one curve per thread; the sigmas are drawn in curve order and a
 * region's successes are taken in curve order, so the factor found depends
 * only on n (and the thread count), not on scheduling. */
static int run_curves(mpz_t d, const mpz_t n, FiCtx* c, double B1, uint64_t curves) {
    int found = 0;
#ifndef NO_ECM
    int nt = mpz_sizeinbase(n, 2) >= ECM_PAR_BITS ? fi_nthreads() : 1;
    CurveJob job[FI_MAXTHREADS];
    for (int t = 0; t < nt; t++) mpz_init(job[t].f);
    for (uint64_t k = 0; k < curves && !found; ) {
        int m = 0;
        for (; m < nt && k < curves; m++, k++) {
            CurveJob* j = &job[m];
            j->n = n;
            j->B1 = B1;
            j->B2 = c->B2;
            /* A user "Sigma" is a Suyama parameter (the classical
             * convention); our own curves use GMP-ECM's fastest
             * parametrisation (batch stage 1) with sigmas drawn from the
             * n-seeded generator (GMP-ECM's own draw is time-seeded). */
            if (c->sigma) { j->param = ECM_PARAM_SUYAMA; j->sigma = (unsigned long)(c->sigma + k); }
            else { j->param = ECM_PARAM_BATCH_SQUARE; j->sigma = (unsigned long)(2 + fi_rand(c) % (UINT64_C(1) << 31)); }
        }
        fi_run_parallel(curve_job, job, sizeof(CurveJob), m);
        for (int t = 0; t < m && !found; t++)
            if (job[t].ok) { mpz_set(d, job[t].f); found = 1; }
    }
    for (int t = 0; t < nt; t++) mpz_clear(job[t].f);
#elif defined(USE_FLINT)
    flint_rand_t st;
    flint_rand_init(st);
    flint_rand_set_seed(st, fi_rand(c), fi_rand(c));
    fmpz_t f, nn;
    fmpz_init(f);
    fmpz_init(nn);
    fmpz_set_mpz(nn, n);
    ulong b1 = (ulong)B1, b2 = c->B2 > 0 ? (ulong)c->B2 : 100 * (ulong)B1;
    if (fmpz_factor_ecm(f, (ulong)curves, b1, b2, st, nn)) {
        mpz_t g; mpz_init(g);
        fmpz_get_mpz(g, f);
        if (nontrivial(g, n)) { mpz_set(d, g); found = 1; }
        mpz_clear(g);
    }
    fmpz_clear(f);
    fmpz_clear(nn);
    flint_rand_clear(st);
#else
    for (uint64_t k = 0; k < curves && !found; k++) {
        uint64_t sigma = c->sigma ? c->sigma + k : 6 + fi_rand(c) % (UINT64_C(1) << 30);
        found = native_ecm_curve(d, n, sigma, B1);
    }
#endif
    return found;
}

int fi_ecm_level(mpz_t d, const mpz_t n, FiCtx* c, int digits) {
    if (tiny_factor(d, n)) return 1;
    double B1;
    uint64_t curves;
    fi_ecm_level_params(digits, &B1, &curves);
    return run_curves(d, n, c, B1, curves);
}

int fi_split_ecm(mpz_t d, const mpz_t n, FiCtx* c, double B1, uint64_t curves) {
    if (tiny_factor(d, n)) return 1;
    if (B1 > 0) return run_curves(d, n, c, B1, curves ? curves : curves_for_B1(B1));
    if (curves) return run_curves(d, n, c, 11000.0, curves);
    /* Default: climb the levels up to half the digit count. */
    int digits = (int)(mpz_sizeinbase(n, 10));
    for (size_t i = 0; i < NLEVELS && LEVELS[i].digits <= (digits + 1) / 2 + 5; i++)
        if (run_curves(d, n, c, LEVELS[i].B1, LEVELS[i].curves)) return 1;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Pollard p-1                                                         */
/* ------------------------------------------------------------------ */

#ifdef NEED_NATIVE_PM1
/* Stage 1: a = 3^E with E = prod of prime powers <= B1; stage 2: for primes
 * q in (B1, B2], accumulate a^q - 1 by multiplying a^(gap) from a table. */
static int native_pm1(mpz_t d, const mpz_t n, double B1d, double B2d) {
    uint64_t B1 = (uint64_t)B1d, B2 = (uint64_t)B2d;
    size_t np;
    const uint32_t* P = fi_small_primes(&np);
    mpz_t a, g, t, acc;
    mpz_inits(a, g, t, acc, NULL);
    mpz_set_ui(a, 3);
    int found = 0;
    uint64_t e = 1;
    for (size_t i = 0; i < np && P[i] <= B1; i++) {
        uint64_t q = P[i];
        while (q * P[i] <= B1) q *= P[i];
        if (e > UINT64_MAX / q) { mpz_powm_ui(a, a, (unsigned long)e, n); e = q; }
        else e *= q;
        if ((i & 1023) == 1023) {                    /* periodic gcd */
            mpz_powm_ui(a, a, (unsigned long)e, n); e = 1;
            mpz_sub_ui(t, a, 1); mpz_gcd(g, t, n);
            if (mpz_cmp_ui(g, 1) != 0) break;
        }
    }
    mpz_powm_ui(a, a, (unsigned long)e, n);
    mpz_sub_ui(t, a, 1); mpz_gcd(g, t, n);
    if (nontrivial(g, n)) { mpz_set(d, g); found = 1; goto out; }
    if (mpz_cmp_ui(g, 1) != 0) goto out;             /* g == n: give up */
    {
        enum { MAXGAP = 512 };
        mpz_t T[MAXGAP / 2 + 1];
        mpz_t a2; mpz_init(a2);
        mpz_powm_ui(a2, a, 2, n);
        mpz_init_set(T[1], a2);
        for (int i = 2; i <= MAXGAP / 2; i++) { mpz_init(T[i]); mpz_mul(T[i], T[i - 1], a2); mpz_mod(T[i], T[i], n); }
        mpz_init(T[0]);
        size_t i0 = 0;
        while (i0 < np && P[i0] <= B1) i0++;
        if (i0 < np) {
            mpz_t x; mpz_init(x);
            mpz_powm_ui(x, a, P[i0], n);
            mpz_set_ui(acc, 1);
            uint64_t qprev = P[i0];
            for (size_t k = i0; k < np && P[k] <= B2; k++) {
                if (k > i0) {
                    uint64_t gap = P[k] - qprev;
                    if (gap > MAXGAP) break;
                    mpz_mul(x, x, T[gap / 2]); mpz_mod(x, x, n);
                    qprev = P[k];
                }
                mpz_sub_ui(t, x, 1);
                mpz_mul(acc, acc, t); mpz_mod(acc, acc, n);
            }
            mpz_gcd(g, acc, n);
            if (nontrivial(g, n)) { mpz_set(d, g); found = 1; }
            mpz_clear(x);
        }
        for (int i = 0; i <= MAXGAP / 2; i++) mpz_clear(T[i]);
        mpz_clear(a2);
    }
out:
    mpz_clears(a, g, t, acc, NULL);
    return found;
}
#endif

/* Words: the u64 walk, which backtracks prime by prime on overshoot. */
static int word_smooth(mpz_t d, const mpz_t n, double B1, int plus, FiCtx* c) {
    uint64_t w = mth_mpz_get_u64(n), g = 0;
    if (plus) { for (int s = 0; s < 3 && !g; s++) g = fi_u64_pp1(w, (uint64_t)B1, fi_rand(c)); }
    else g = fi_u64_pm1(w, (uint64_t)B1);
    if (g > 1 && g < w) { mth_mpz_set_u64(d, (uint64_t)g); return 1; }
    return 0;
}

#ifndef NO_ECM
/* One GMP-ECM p-1 / p+1 run.  Returns 1 = split, 0 = nothing, -1 = every
 * prime of n caught at once (gcd = n). */
static int gmpecm_smooth(mpz_t d, const mpz_t n, unsigned long rngseed, int method, double B1,
                         double B2, int seed) {
    ecm_params p;
    ecm_init(p);
    p->stop_asap = fi_stop_requested;       /* TimeConstrained deadline */
    gmp_randseed_ui(p->rng, rngseed);
    p->method = method;
    /* Fixed starting values (GMP-ECM would otherwise draw one at random):
     * 3 for p-1; for p+1 the seeds 2/7 and 6/5 that GMP-ECM recommends (they
     * make 6 resp. 4 divide the group order), then 3, 5, ...  A p+1 seed x0
     * finds p | n when p - (x0^2-4 / p) is smooth, so each seed is a fresh
     * coin flip between "p + 1 smooth" and "p - 1 smooth". */
    if (method == ECM_PM1 || seed >= 2) {
        mpz_set_ui(p->x, method == ECM_PM1 ? 3 : (unsigned long)(2 * seed - 1));
    } else {
        mpz_t den;
        mpz_init_set_ui(den, seed == 0 ? 7 : 5);
        mpz_invert(den, den, n);
        mpz_mul_ui(p->x, den, seed == 0 ? 2 : 6);
        mpz_mod(p->x, p->x, n);
        mpz_clear(den);
    }
    if (B2 > 0) mpz_set_d(p->B2, B2);
    mpz_t f;
    mpz_init(f);
    int r = ecm_factor(f, (mpz_ptr)n, B1, p);
    int out = 0;
    if (r > 0 && nontrivial(f, n)) { mpz_set(d, f); out = 1; }
    else if (r > 0 && mpz_cmp(f, n) == 0) out = -1;
    mpz_clear(f);
    ecm_clear(p);
    return out;
}

/* Run, and on overshoot retry with B1 shrunk 16x at a time: the primes of n
 * separate once B1 drops below the largest prime of some p -/+ 1.  Touches
 * no shared state, so it runs on a worker thread. */
static int gmpecm_smooth_backoff(mpz_t d, const mpz_t n, unsigned long rngseed, int method, double B1,
                                 double B2, int seed) {
    for (double b = B1; b >= 20 && !fi_stop_requested(); b /= 16) {
        int r = gmpecm_smooth(d, n, rngseed, method, b, b == B1 ? B2 : 0, seed);
        if (r >= 0) return r;
    }
    return 0;
}
#endif

int fi_split_pm1(mpz_t d, const mpz_t n, FiCtx* c, double B1, double B2) {
    if (B1 <= 0) B1 = 1e6;
    if (tiny_factor(d, n)) return 1;
    if (mpz_sizeinbase(n, 2) <= 64) return word_smooth(d, n, B1, 0, c);
#ifndef NO_ECM
    return gmpecm_smooth_backoff(d, n, (unsigned long)fi_rand(c), ECM_PM1, B1, B2, 0);
#else
    (void)c;
    return native_pm1(d, n, B1, B2 > 0 ? B2 : 100 * B1);
#endif
}

/* ------------------------------------------------------------------ */
/* Williams p+1                                                        */
/* ------------------------------------------------------------------ */

#ifdef NEED_NATIVE_PP1
/* V_{2k} = V_k^2 - 2, V_{2k+1} = V_k V_{k+1} - V_1: binary ladder. */
static void lucas_v(mpz_t v, const mpz_t v1, uint64_t k, const mpz_t n) {
    mpz_t a, b;
    mpz_init_set(a, v1);                    /* V_1 */
    mpz_init(b);                            /* V_2 */
    mpz_mul(b, v1, v1); mpz_sub_ui(b, b, 2); mpz_mod(b, b, n);
    int top = 63;
    while (top > 0 && !((k >> top) & 1)) top--;
    for (int i = top - 1; i >= 0; i--) {
        if ((k >> i) & 1) {
            mpz_mul(a, a, b); mpz_sub(a, a, v1); mpz_mod(a, a, n);
            mpz_mul(b, b, b); mpz_sub_ui(b, b, 2); mpz_mod(b, b, n);
        } else {
            mpz_mul(b, a, b); mpz_sub(b, b, v1); mpz_mod(b, b, n);
            mpz_mul(a, a, a); mpz_sub_ui(a, a, 2); mpz_mod(a, a, n);
        }
    }
    if (k == 0) mpz_set_ui(v, 2); else mpz_set(v, a);
    mpz_clears(a, b, NULL);
}

static int native_pp1(mpz_t d, const mpz_t n, double B1d, uint64_t seed) {
    uint64_t B1 = (uint64_t)B1d;
    size_t np;
    const uint32_t* P = fi_small_primes(&np);
    mpz_t v, g, t;
    mpz_inits(v, g, t, NULL);
    mpz_set_ui(v, (unsigned long)(3 + seed % 1000));
    for (size_t i = 0; i < np && P[i] <= B1; i++) {
        uint64_t q = P[i];
        while (q * P[i] <= B1) q *= P[i];
        lucas_v(v, v, q, n);
    }
    mpz_sub_ui(t, v, 2);
    mpz_gcd(g, t, n);
    int ok = nontrivial(g, n);
    if (ok) mpz_set(d, g);
    mpz_clears(v, g, t, NULL);
    return ok;
}
#endif

/* One p+1 seed, run on a worker thread: x0 from the seed index (GMP-ECM),
 * or from the pre-drawn word (FLINT, native). */
typedef struct {
    mpz_srcptr n;
    double     B1, B2;
    int        seed;
    uint64_t   rng;
    mpz_t      f;
    int        ok;
} SeedJob;

static void pp1_seed_job(void* arg) {
    SeedJob* j = (SeedJob*)arg;
#ifndef NO_ECM
    j->ok = gmpecm_smooth_backoff(j->f, j->n, (unsigned long)j->rng, ECM_PP1, j->B1, j->B2, j->seed) == 1;
#elif defined(USE_FLINT)
    fmpz_t f, nn;
    fmpz_init(f); fmpz_init(nn);
    fmpz_set_mpz(nn, j->n);
    ulong b2s = (ulong)sqrt(100.0 * j->B1);
    j->ok = fmpz_factor_pp1(f, nn, (ulong)j->B1, b2s, (ulong)(3 + j->rng % 1000));
    if (j->ok) {
        fmpz_get_mpz(j->f, f);
        j->ok = nontrivial(j->f, j->n);
    }
    fmpz_clear(f); fmpz_clear(nn);
#else
    j->ok = native_pp1(j->f, j->n, j->B1, j->rng);
#endif
}

int fi_split_pp1(mpz_t d, const mpz_t n, FiCtx* c, double B1, double B2) {
    return fi_split_pp1_seeds(d, n, c, B1, B2, 3);
}

/* p+1 only works when the seed's discriminant is a non-residue mod p, so
 * the method runs several seeds (three: GMP-ECM's recommendation).  Seeds are
 * independent: with threads they go out one per thread, and a region's
 * successes are taken in seed order, so the factor found is the serial one. */
int fi_split_pp1_seeds(mpz_t d, const mpz_t n, FiCtx* c, double B1, double B2, int nseeds) {
    if (B1 <= 0) B1 = 1e6;
    if (tiny_factor(d, n)) return 1;
    if (mpz_sizeinbase(n, 2) <= 64) return word_smooth(d, n, B1, 1, c);
    int nt = fi_nthreads();
    if (nt > nseeds) nt = nseeds;
    if (nt < 1) nt = 1;
    SeedJob job[FI_MAXTHREADS];
    for (int t = 0; t < nt; t++) mpz_init(job[t].f);
    int found = 0;
    for (int s = 0; s < nseeds && !found; ) {
        int m = 0;
        for (; m < nt && s < nseeds; m++, s++) {
            SeedJob* j = &job[m];
            j->n = n;
            j->B1 = B1;
            j->B2 = B2;
            j->seed = s;
            j->rng = fi_rand(c);
            j->ok = 0;
        }
        fi_run_parallel(pp1_seed_job, job, sizeof(SeedJob), m);
        for (int t = 0; t < m && !found; t++)
            if (job[t].ok) { mpz_set(d, job[t].f); found = 1; }
    }
    for (int t = 0; t < nt; t++) mpz_clear(job[t].f);
    return found;
}
