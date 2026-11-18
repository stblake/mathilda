/*
 * fi_qs.c — the self-initialising quadratic sieve (SIQS).
 *
 * For kN (k a Knuth–Schroeppel multiplier) and a factor base of primes p
 * with (kN/p) = 1, polynomials g(x) = ((Ax + B)^2 - kN) / A are sieved over
 * x in [-M, M); (Ax + B)^2 == A g(x) (mod N), so every x with g(x) smooth
 * over the base (up to one large prime) is a relation, and a GF(2)
 * dependency among the relations splits N (fi_gf2.c).
 *
 *   A       product of s base primes q_l, chosen so A ~ sqrt(2kN)/M
 *           (keeps |g| ~ M sqrt(kN/2) over the whole interval);
 *   B       one of 2^(s-1) square roots of kN mod A, B = +-B_1 +- ... +- B_s;
 *           successive B differ in one sign (Gray code), so every sieve root
 *           moves by a precomputed 2 B_l A^-1 mod p: no modular inversion
 *           per polynomial ("self-initialisation", Alford–Pomerance / Contini);
 *   sieve   byte logarithms, 32 KB blocks (one L1 data cache): primes below
 *           the block size are sieved block by block; larger primes are
 *           bucket-sieved once per polynomial into per-block hit lists,
 *           which double as the trial-division lists for candidates;
 *           primes below 2^5 are not sieved (small-prime variation) and the
 *           threshold is lowered to compensate;
 *   scan    8 bytes at a time: the array is initialised to 128 - T so a
 *           candidate is any byte with its top bit set;
 *   LP      single large prime up to lpmult * pmax (partials pair up in
 *           fi_gf2.c).
 *
 * Parameters (factor-base size, interval, large-prime multiplier) come from
 * a table by bit size of kN, interpolated.
 *
 * Parallelism (MATHILDA_THREADS): polynomial families are independent, so
 * from ~45 digits each worker thread takes one A and sieves all 2^(s-1) of
 * its polynomials into a private relation log; the main thread then replays
 * the logs into the relation store in worker order (so the result is the
 * same for a given thread count) and decides whether to solve.  Workers run
 * with every signal blocked inside a tc_async_defer region, with the main
 * thread parked in pthread_join: a TimeConstrained timeout is taken on the
 * main thread once the region (one A per worker) ends.
 */
#include "fi_internal.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define QS_BLOCK      32768u
#define QS_BLOCK_BITS 15
#define QS_SMALL_P    64u        /* primes below this are not sieved */
#define QS_MED_P      32768u     /* primes from here up are bucket-sieved */
#define QS_DLP_BITS   240u       /* kN this large uses double large primes */
#define QS_PAR_BITS   150u       /* kN this large sieves in parallel */
#define QS_MAXS       20

/* ------------------------------------------------------------------ */
/* Parameters                                                          */
/* ------------------------------------------------------------------ */

/* Tuned on an i9-9880H (32 KB L1d, 256 KB L2) against balanced semiprimes;
 * `blocks` is 32 KB blocks per side of the interval, `slack` the bits below
 * the full log the threshold is set to (it pays for the unsieved small
 * primes and lets more large-prime partials through). */
static const struct { unsigned bits, fb, blocks, lpmult, slack; } QS_PARAMS[] = {
    {  64,    60,  1,  30, 3 },
    { 100,   200,  1,  40, 4 },
    { 120,   400,  1,  60, 5 },
    { 133,   600,  1, 100, 6 },
    { 166,  2500,  2, 250, 8 },
    { 200,  6000,  2, 250, 8 },
    { 233, 18000,  3, 250, 8 },
    { 266, 34000,  4, 300, 9 },
    { 300, 55000,  6, 350, 10 },
    { 333, 80000,  8, 400, 10 },
    { 366, 110000, 10, 400, 11 },
};
#define QS_NPARAMS (sizeof QS_PARAMS / sizeof QS_PARAMS[0])

static void qs_params(unsigned bits, unsigned* fb, unsigned* blocks, unsigned* lpmult, double* slack) {
    size_t i = 1;
    if (bits <= QS_PARAMS[0].bits) i = 0;
    while (i < QS_NPARAMS && bits > QS_PARAMS[i].bits) i++;
    if (i == 0 || i == QS_NPARAMS) {
        size_t k = i == 0 ? 0 : QS_NPARAMS - 1;
        *fb = QS_PARAMS[k].fb; *blocks = QS_PARAMS[k].blocks;
        *lpmult = QS_PARAMS[k].lpmult; *slack = QS_PARAMS[k].slack;
        return;
    }
    double t = (double)(bits - QS_PARAMS[i - 1].bits) / (QS_PARAMS[i].bits - QS_PARAMS[i - 1].bits);
    *fb = (unsigned)(QS_PARAMS[i - 1].fb + t * ((double)QS_PARAMS[i].fb - QS_PARAMS[i - 1].fb));
    *blocks = t < 0.5 ? QS_PARAMS[i - 1].blocks : QS_PARAMS[i].blocks;
    *lpmult = (unsigned)(QS_PARAMS[i - 1].lpmult + t * ((double)QS_PARAMS[i].lpmult - QS_PARAMS[i - 1].lpmult));
    *slack = QS_PARAMS[i - 1].slack + t * ((double)QS_PARAMS[i].slack - QS_PARAMS[i - 1].slack);
}

/* ------------------------------------------------------------------ */
/* Word arithmetic mod p (p < 2^32)                                    */
/* ------------------------------------------------------------------ */

static inline int ctz32(uint32_t x) {
    int c = 0;
    while (!(x & 1)) { x >>= 1; c++; }
    return c;
}

static inline uint32_t mulmod32(uint32_t a, uint32_t b, uint32_t p) {
    return (uint32_t)(((uint64_t)a * b) % p);
}

static uint32_t powmod32(uint32_t b, uint32_t e, uint32_t p) {
    uint32_t r = 1;
    while (e) {
        if (e & 1) r = mulmod32(r, b, p);
        b = mulmod32(b, b, p);
        e >>= 1;
    }
    return r;
}

static uint32_t invmod32(uint32_t a, uint32_t p) {
    int64_t t0 = 0, t1 = 1, r0 = p, r1 = a % p;
    while (r1) {
        int64_t q = r0 / r1, r2 = r0 - q * r1, t2 = t0 - q * t1;
        r0 = r1; r1 = r2; t0 = t1; t1 = t2;
    }
    return (uint32_t)(t0 < 0 ? t0 + p : t0);
}

/* Square root of a quadratic residue a mod odd prime p (Tonelli–Shanks). */
static uint32_t sqrtmod32(uint32_t a, uint32_t p) {
    a %= p;
    if (a == 0) return 0;
    if (p % 4 == 3) return powmod32(a, (p + 1) / 4, p);
    uint32_t q = p - 1, s = 0;
    while (!(q & 1)) { q >>= 1; s++; }
    uint32_t z = 2;
    while (powmod32(z, (p - 1) / 2, p) != p - 1) z++;
    uint32_t m = s, c = powmod32(z, q, p), t = powmod32(a, q, p), r = powmod32(a, (q + 1) / 2, p);
    while (t != 1) {
        uint32_t i = 0, tt = t;
        while (tt != 1) { tt = mulmod32(tt, tt, p); i++; }
        uint32_t b = c;
        for (uint32_t j = 0; j + i + 1 < m; j++) b = mulmod32(b, b, p);
        m = i;
        c = mulmod32(b, b, p);
        t = mulmod32(t, c, p);
        r = mulmod32(r, b, p);
    }
    return r;
}

/* ------------------------------------------------------------------ */
/* Sieve state                                                         */
/* ------------------------------------------------------------------ */

typedef struct {
    uint32_t* hits;     /* (prime index << 15) | offset in block */
    size_t    n, cap;
} Bucket;

typedef struct {
    mpz_srcptr N;
    mpz_t     kN;
    uint32_t  nfb;          /* base primes are fb[1..nfb] (fb[0]: sign) */
    uint32_t* fb;
    uint32_t* sqrtkn;
    uint8_t*  logp;
    uint32_t* ainv;
    uint32_t* r1;           /* roots as positions in [0, p): x + M == r (mod p) */
    uint32_t* r2;
    uint32_t* pos1;         /* medium primes: next hit within the current block */
    uint32_t* pos2;
    uint32_t* bainv;        /* [s][nfb + 1] */
    uint32_t  med_end;      /* fb indices < med_end have p < QS_BLOCK */
    uint32_t  sieve_start;  /* first sieved index (p >= QS_SMALL_P) */
    uint32_t  M, L, nblocks;
    uint64_t  lp_max;
    uint64_t  dlp_max;      /* 0: single large prime only */
    uint64_t  pmax2;        /* pmax^2 */
    uint64_t  rng;          /* cofactor-splitting curves */
    uint8_t   init;         /* 128 - threshold */
    double    scale;        /* log2 -> byte units */
    /* polynomial */
    int       s;
    uint32_t  qidx[QS_MAXS];
    mpz_t     A, B, C, Bl[QS_MAXS];
    int       sgn[QS_MAXS];
    /* sieve memory */
    uint8_t*  sieve;
    uint16_t* candidx;      /* offset -> candidate index + 1 while scanning */
    uint32_t (*bigbuf)[64]; /* bucket primes per candidate */
    uint8_t*  nbig;
    Bucket*   bucket;
    /* relation buffers */
    uint32_t* col;
    uint32_t* exps;
    mpz_t     g, X, t;
    FiRelLog  log;          /* relations found, replayed by the main thread */
} QS;


/* ------------------------------------------------------------------ */
/* Polynomial selection                                                */
/* ------------------------------------------------------------------ */

/* Choose A = prod q_l with A close to sqrt(2kN)/M: s - 1 random distinct
 * base primes from a window [want/2, 2 want] (want = target^(1/s), widened
 * if too narrow), then the last q_l fitted to the residual.  A key of every
 * A used is kept so no polynomial family repeats. */
static uint32_t fb_index_at_least(const QS* Q, uint32_t lo, uint32_t hi, double v) {
    while (lo < hi) {
        uint32_t mid = lo + (hi - lo) / 2;
        if ((double)Q->fb[mid] < v) lo = mid + 1; else hi = mid;
    }
    return lo;
}

static int choose_A(QS* Q, FiCtx* c, uint64_t* used, size_t* nused) {
    long ex;
    double mant = mpz_get_d_2exp(&ex, Q->kN);
    double target = 0.5 * (log2(2.0 * mant) + (double)ex) - log2((double)Q->M);
    uint32_t lo_i = Q->sieve_start + 1, hi_i = Q->nfb;
    /* About 11-bit q_l, but never larger than the base allows. */
    double qbits = log2((double)Q->fb[hi_i]) - 1.5;
    if (qbits > 11.0) qbits = 11.0;
    int s = (int)ceil(target / qbits);
    if (s < 2) s = 2;
    if (s > QS_MAXS) s = QS_MAXS;
    double each = target / s;
    double ratio = 2.0;
    uint32_t a, b;
    for (;;) {
        a = fb_index_at_least(Q, lo_i, hi_i, pow(2.0, each) / ratio);
        b = fb_index_at_least(Q, lo_i, hi_i, pow(2.0, each) * ratio);
        if (b > hi_i) b = hi_i;
        if (b - a >= (uint32_t)(4 * s) || (a == lo_i && b == hi_i)) break;
        ratio *= 1.5;
    }
    if (b <= a + (uint32_t)s) return 0;
    double lo_bits = log2((double)Q->fb[lo_i]), hi_bits = log2((double)Q->fb[hi_i]);
    for (int tries = 0; tries < 2000; tries++) {
        double acc = 0;
        int ok = 1;
        for (int l = 0; l < s - 1 && ok; l++) {
            uint32_t idx;
            int guard = 0;
            do {
                idx = a + (uint32_t)(fi_rand(c) % (b - a));
                ok = 1;
                for (int j = 0; j < l; j++) if (Q->qidx[j] == idx) ok = 0;
            } while (!ok && ++guard < 64);
            Q->qidx[l] = idx;
            acc += log2((double)Q->fb[idx]);
        }
        if (!ok) continue;
        double rest = target - acc;
        if (rest < lo_bits || rest > hi_bits) continue;
        uint32_t li = fb_index_at_least(Q, lo_i, hi_i, pow(2.0, rest));
        if (li > hi_i) li = hi_i;
        for (int j = 0; j < s - 1; j++) if (Q->qidx[j] == li) ok = 0;
        if (!ok) continue;
        Q->qidx[s - 1] = li;
        Q->s = s;
        mpz_set_ui(Q->A, 1);
        for (int l = 0; l < s; l++) mpz_mul_ui(Q->A, Q->A, Q->fb[Q->qidx[l]]);
        uint64_t key = mth_mpz_get_u64(Q->A) ^ ((uint64_t)mpz_sizeinbase(Q->A, 2) << 56);
        int dup = 0;
        for (size_t u = 0; u < *nused; u++) if (used[u] == key) { dup = 1; break; }
        if (dup) continue;
        used[(*nused)++] = key;
        return 1;
    }
    return 0;
}

/* Set up B_l, the first B, C, A^-1 mod p, 2 B_l A^-1 mod p, and roots. */
static void init_poly(QS* Q) {
    int s = Q->s;
    mpz_set_ui(Q->B, 0);
    for (int l = 0; l < s; l++) {
        uint32_t q = Q->fb[Q->qidx[l]];
        mpz_divexact_ui(Q->t, Q->A, q);                 /* A / q */
        uint32_t aq = (uint32_t)mpz_fdiv_ui(Q->t, q);
        uint32_t gam = mulmod32(Q->sqrtkn[Q->qidx[l]], invmod32(aq, q), q);
        if (gam > q / 2) gam = q - gam;
        mpz_mul_ui(Q->Bl[l], Q->t, gam);
        mpz_add(Q->B, Q->B, Q->Bl[l]);
        Q->sgn[l] = 1;
    }
    /* C = (B^2 - kN) / A, exact because B^2 == kN (mod A). */
    mpz_mul(Q->C, Q->B, Q->B);
    mpz_sub(Q->C, Q->C, Q->kN);
    mpz_divexact(Q->C, Q->C, Q->A);
    uint32_t W = Q->nfb + 1;
    for (uint32_t i = Q->sieve_start; i <= Q->nfb; i++) {
        uint32_t p = Q->fb[i];
        uint32_t am = (uint32_t)mpz_fdiv_ui(Q->A, p);
        if (am == 0) { Q->ainv[i] = 0; continue; }       /* p | A: not sieved */
        uint32_t ai = invmod32(am, p);
        Q->ainv[i] = ai;
        for (int l = 0; l < s; l++) {
            uint32_t bl = (uint32_t)mpz_fdiv_ui(Q->Bl[l], p);
            Q->bainv[(size_t)l * W + i] = mulmod32(mulmod32(2, bl, p), ai, p);
        }
        uint32_t bm = (uint32_t)mpz_fdiv_ui(Q->B, p);
        uint32_t t = Q->sqrtkn[i], Mp = Q->M % p;
        uint32_t x1 = mulmod32(ai, (t + p - bm) % p, p);
        uint32_t x2 = mulmod32(ai, (2 * (uint64_t)p - t - bm) % p, p);
        Q->r1[i] = (uint32_t)(((uint64_t)x1 + Mp) % p);
        Q->r2[i] = (uint32_t)(((uint64_t)x2 + Mp) % p);
    }
}

/* Gray-code step to the next B: flip the sign of term v.  Roots of the
 * medium primes move here; the bucket-sieved primes move inside
 * fill_buckets, in the same pass that sieves them. */
static int next_poly(QS* Q, int v) {
    int e = Q->sgn[v];
    /* B' = B - 2 e B_v; roots move by + e * 2 B_v A^-1. */
    mpz_submul_ui(Q->B, Q->Bl[v], 2);
    if (e < 0) mpz_addmul_ui(Q->B, Q->Bl[v], 4);
    Q->sgn[v] = -e;
    mpz_mul(Q->C, Q->B, Q->B);
    mpz_sub(Q->C, Q->C, Q->kN);
    mpz_divexact(Q->C, Q->C, Q->A);
    const uint32_t* d = Q->bainv + (size_t)v * (Q->nfb + 1);
    for (uint32_t i = Q->sieve_start; i < Q->med_end; i++) {
        if (!Q->ainv[i]) continue;
        uint32_t p = Q->fb[i], dv = d[i];
        if (e > 0) {
            uint32_t a = Q->r1[i] + dv; Q->r1[i] = a >= p ? a - p : a;
            a = Q->r2[i] + dv;          Q->r2[i] = a >= p ? a - p : a;
        } else {
            Q->r1[i] = Q->r1[i] >= dv ? Q->r1[i] - dv : Q->r1[i] + p - dv;
            Q->r2[i] = Q->r2[i] >= dv ? Q->r2[i] - dv : Q->r2[i] + p - dv;
        }
    }
    return e;
}

/* ------------------------------------------------------------------ */
/* Sieving                                                             */
/* ------------------------------------------------------------------ */

/* Move the bucket-sieved primes' roots (when v >= 0: Gray-code term v with
 * sign e) and drop each of their hits in the interval into the bucket of
 * its block. */
static void fill_buckets(QS* Q, int v, int e) {
    const uint32_t* d = v >= 0 ? Q->bainv + (size_t)v * (Q->nfb + 1) : NULL;
    const uint32_t L = Q->L;
    uint32_t* wp[64];                       /* write cursor per block */
    for (uint32_t b = 0; b < Q->nblocks; b++) wp[b] = Q->bucket[b].hits;
    for (uint32_t i = Q->med_end; i <= Q->nfb; i++) {
        if (!Q->ainv[i]) continue;
        uint32_t p = Q->fb[i], r1 = Q->r1[i], r2 = Q->r2[i];
        if (d) {
            uint32_t dv = d[i];
            if (e > 0) {
                r1 += dv; if (r1 >= p) r1 -= p;
                r2 += dv; if (r2 >= p) r2 -= p;
            } else {
                r1 = r1 >= dv ? r1 - dv : r1 + p - dv;
                r2 = r2 >= dv ? r2 - dv : r2 + p - dv;
            }
            Q->r1[i] = r1; Q->r2[i] = r2;
        }
        uint32_t tag = i << QS_BLOCK_BITS;
        for (uint32_t pos = r1; pos < L; pos += p)
            *wp[pos >> QS_BLOCK_BITS]++ = tag | (pos & (QS_BLOCK - 1));
        if (r2 != r1)
            for (uint32_t pos = r2; pos < L; pos += p)
                *wp[pos >> QS_BLOCK_BITS]++ = tag | (pos & (QS_BLOCK - 1));
    }
    for (uint32_t b = 0; b < Q->nblocks; b++) Q->bucket[b].n = (size_t)(wp[b] - Q->bucket[b].hits);
}

/* Each bucket is sized once for the worst case, so filling needs no
 * capacity checks: a prime p >= QS_BLOCK hits a block at most twice (one
 * per root), and in expectation 2 QS_BLOCK / p times. */
static void alloc_buckets(QS* Q) {
    size_t cap = 0;
    for (uint32_t i = Q->med_end; i <= Q->nfb; i++) {
        uint32_t p = Q->fb[i];
        cap += 2 * (QS_BLOCK / p + 1);
    }
    for (uint32_t b = 0; b < Q->nblocks; b++) {
        Q->bucket[b].cap = cap + 16;
        Q->bucket[b].hits = malloc(Q->bucket[b].cap * sizeof(uint32_t));
        Q->bucket[b].n = 0;
    }
}

/* Medium primes keep their next hit relative to the current block in
 * pos1/pos2 (start of polynomial: the roots themselves, since a root < p <
 * QS_BLOCK), so moving to the next block is a subtraction, not a modulo. */
static void start_medium(QS* Q) {
    for (uint32_t i = Q->sieve_start; i < Q->med_end; i++) {
        Q->pos1[i] = Q->r1[i];
        Q->pos2[i] = Q->r2[i];
    }
}

static void sieve_block(QS* Q, uint32_t blk) {
    uint8_t* S = Q->sieve;
    memset(S, Q->init, QS_BLOCK);
    for (uint32_t i = Q->sieve_start; i < Q->med_end; i++) {
        if (!Q->ainv[i]) continue;
        uint32_t p = Q->fb[i];
        uint8_t lg = Q->logp[i];
        uint32_t x1 = Q->pos1[i], x2 = Q->pos2[i];
        if (x1 > x2) { uint32_t t = x1; x1 = x2; x2 = t; }
        /* Both roots in one pass while the larger is in range. */
        uint32_t d = x2 - x1;
        uint32_t lim = QS_BLOCK - d;
        while (x1 + 2 * p < lim) {
            S[x1] += lg; S[x1 + d] += lg; x1 += p;
            S[x1] += lg; S[x1 + d] += lg; x1 += p;
        }
        while (x1 < lim) { S[x1] += lg; S[x1 + d] += lg; x1 += p; }
        x2 = x1 + d;
        if (x1 < QS_BLOCK) { S[x1] += lg; x1 += p; }
        /* d == 0 (p | k): one root, sieved twice above; undo is not
         * needed for a threshold test, the doubled log only biases up. */
        Q->pos1[i] = x1 - QS_BLOCK;
        Q->pos2[i] = x2 - QS_BLOCK;
    }
    const Bucket* bk = &Q->bucket[blk];
    for (size_t k = 0; k < bk->n; k++) {
        uint32_t h = bk->hits[k];
        S[h & (QS_BLOCK - 1)] += Q->logp[h >> QS_BLOCK_BITS];
    }
    (void)blk;
}

/* Sort (col, exp) pairs by column (short lists: insertion sort). */
static void sort_pairs(uint32_t* col, uint32_t* ex, size_t n) {
    for (size_t i = 1; i < n; i++) {
        uint32_t c = col[i], e = ex[i];
        size_t j = i;
        while (j > 0 && col[j - 1] > c) { col[j] = col[j - 1]; ex[j] = ex[j - 1]; j--; }
        col[j] = c; ex[j] = e;
    }
}

/* Trial-divide g(x) for the candidate at absolute position pos and record a
 * relation if it is smooth or a single large-prime partial.  `big` lists
 * the bucket-sieved primes that hit this position. */
static void check_candidate(QS* Q, uint32_t blk, uint32_t off,
                            const uint32_t* big, size_t nbig) {
    uint32_t pos = blk * QS_BLOCK + off;
    long x = (long)pos - (long)Q->M;
    /* X = A x + B;  g = (X^2 - kN) / A = A x^2 + 2 B x + C. */
    mpz_mul_si(Q->X, Q->A, x);
    mpz_add(Q->X, Q->X, Q->B);
    mpz_mul(Q->g, Q->X, Q->X);
    mpz_sub(Q->g, Q->g, Q->kN);
    mpz_divexact(Q->g, Q->g, Q->A);
    size_t cnt = 0;
    if (mpz_sgn(Q->g) < 0) { Q->col[cnt] = 0; Q->exps[cnt++] = 1; mpz_neg(Q->g, Q->g); }
    if (mpz_sgn(Q->g) == 0) return;
    /* Small (unsieved) primes and A's primes: direct division. */
    for (uint32_t i = 1; i < Q->sieve_start; i++) {
        uint32_t p = Q->fb[i], e = 0;
        while (mpz_divisible_ui_p(Q->g, p)) { mpz_divexact_ui(Q->g, Q->g, p); e++; }
        if (e) { Q->col[cnt] = i; Q->exps[cnt++] = e; }
    }
    for (int l = 0; l < Q->s; l++) {
        uint32_t i = Q->qidx[l], p = Q->fb[i], e = 1;    /* A contributes q once */
        while (mpz_divisible_ui_p(Q->g, p)) { mpz_divexact_ui(Q->g, Q->g, p); e++; }
        Q->col[cnt] = i; Q->exps[cnt++] = e;
    }
    /* Medium primes: a root match means p | g(x). */
    for (uint32_t i = Q->sieve_start; i < Q->med_end; i++) {
        if (!Q->ainv[i]) continue;
        uint32_t p = Q->fb[i], m = pos % p;
        if (m != Q->r1[i] && m != Q->r2[i]) continue;
        uint32_t e = 0;
        while (mpz_divisible_ui_p(Q->g, p)) { mpz_divexact_ui(Q->g, Q->g, p); e++; }
        if (e) { Q->col[cnt] = i; Q->exps[cnt++] = e; }
    }
    /* Bucket-sieved primes. */
    for (size_t k = 0; k < nbig; k++) {
        uint32_t i = big[k], p = Q->fb[i], e = 0;
        while (mpz_divisible_ui_p(Q->g, p)) { mpz_divexact_ui(Q->g, Q->g, p); e++; }
        if (e) { Q->col[cnt] = i; Q->exps[cnt++] = e; }
    }
    uint64_t L1 = 1, L2 = 1;
    if (mpz_cmp_ui(Q->g, 1) != 0) {
        if (mpz_sizeinbase(Q->g, 2) > 64) return;
        uint64_t c = mth_mpz_get_u64(Q->g);
        if (c <= Q->lp_max) {
            L2 = c;                              /* < pmax^2: prime */
        } else {
            /* Double large prime: c = L1 L2 with both <= lp_max.  A cofactor
             * below pmax^2 is a prime too large to use; a prime cofactor is
             * useless; otherwise split it in word arithmetic. */
            if (!Q->dlp_max || c > Q->dlp_max || c < Q->pmax2 || fi_u64_is_prime(c)) return;
            uint64_t a = fi_u64_ecm(c, &Q->rng, 8);
            if (!a) a = fi_u64_squfof(c, 0);
            if (!a) return;
            uint64_t b = c / a;
            if (a > Q->lp_max || b > Q->lp_max) return;
            L1 = a; L2 = b;
        }
    }
    sort_pairs(Q->col, Q->exps, cnt);
    mpz_mod(Q->t, Q->X, Q->N);
    fi_rellog_push(&Q->log, 0, Q->t, Q->col, Q->exps, cnt, L1, L2);
}

/* Scan one sieved block: gather the candidate offsets, attach to each the
 * bucket primes that hit it (one pass over the bucket), then check them. */
#define QS_MAXCAND 256
#define QS_MAXBIG  64
static size_t scan_block(QS* Q, uint32_t blk) {
    uint32_t cand[QS_MAXCAND];
    size_t nc = 0;
    const uint64_t* w = (const uint64_t*)Q->sieve;
    for (uint32_t j = 0; j < QS_BLOCK / 8 && nc < QS_MAXCAND; j++) {
        if (!(w[j] & UINT64_C(0x8080808080808080))) continue;
        for (uint32_t b = 0; b < 8 && nc < QS_MAXCAND; b++)
            if (Q->sieve[8 * j + b] & 0x80) cand[nc++] = 8 * j + b;
    }
    if (!nc) return 0;
    uint32_t (*big)[QS_MAXBIG] = Q->bigbuf;
    uint8_t* nbig = Q->nbig;
    for (size_t k = 0; k < nc; k++) { Q->candidx[cand[k]] = (uint16_t)(k + 1); nbig[k] = 0; }
    const Bucket* bk = &Q->bucket[blk];
    for (size_t k = 0; k < bk->n; k++) {
        uint32_t h = bk->hits[k];
        uint16_t ci = Q->candidx[h & (QS_BLOCK - 1)];
        if (ci && nbig[ci - 1] < QS_MAXBIG) big[ci - 1][nbig[ci - 1]++] = h >> QS_BLOCK_BITS;
    }
    for (size_t k = 0; k < nc; k++) {
        Q->candidx[cand[k]] = 0;
        check_candidate(Q, blk, cand[k], big[k], nbig[k]);
    }
    return nc;
}

/* ------------------------------------------------------------------ */
/* Driver                                                              */
/* ------------------------------------------------------------------ */

/* Sieve every polynomial of the family whose A is already chosen. */
static void sieve_family(QS* Q) {
    init_poly(Q);
    uint32_t npoly = 1u << (Q->s - 1);
    for (uint32_t pi = 0; pi < npoly; pi++) {
        int v = -1, e = 0;
        if ((pi & 15u) == 15u && fi_stop_requested()) break;   /* TimeConstrained */
        if (pi) { v = 1 + ctz32(pi); e = next_poly(Q, v); }
        fill_buckets(Q, v, e);
        start_medium(Q);
        for (uint32_t blk = 0; blk < Q->nblocks; blk++) {
            sieve_block(Q, blk);
            scan_block(Q, blk);
        }
    }
}

/* A worker shares the master's read-only factor base and owns the rest. */
static void worker_init(QS* W, const QS* M, uint64_t seed) {
    *W = *M;                                    /* scalars + shared pointers */
    uint32_t nfb = M->nfb;
    W->rng = seed | 1;
    W->ainv = calloc(nfb + 1, sizeof(uint32_t));
    W->r1 = calloc(nfb + 1, sizeof(uint32_t));
    W->r2 = calloc(nfb + 1, sizeof(uint32_t));
    W->pos1 = calloc(nfb + 1, sizeof(uint32_t));
    W->pos2 = calloc(nfb + 1, sizeof(uint32_t));
    W->bainv = malloc((size_t)QS_MAXS * (nfb + 1) * sizeof(uint32_t));
    W->sieve = malloc(QS_BLOCK + 8);
    W->candidx = calloc(QS_BLOCK, sizeof(uint16_t));
    W->bigbuf = malloc(256 * sizeof *W->bigbuf);
    W->nbig = malloc(256);
    W->bucket = calloc(W->nblocks, sizeof(Bucket));
    W->col = malloc((nfb + 4 + QS_MAXS) * sizeof(uint32_t));
    W->exps = malloc((nfb + 4 + QS_MAXS) * sizeof(uint32_t));
    mpz_inits(W->A, W->B, W->C, W->g, W->X, W->t, NULL);
    for (int l = 0; l < QS_MAXS; l++) mpz_init(W->Bl[l]);
    memset(&W->log, 0, sizeof W->log);
    alloc_buckets(W);
}

static void worker_free(QS* W) {
    for (uint32_t b = 0; b < W->nblocks; b++) free(W->bucket[b].hits);
    free(W->bucket);
    free(W->ainv); free(W->r1); free(W->r2); free(W->pos1); free(W->pos2); free(W->bainv);
    free(W->sieve); free(W->candidx); free(W->bigbuf); free(W->nbig);
    free(W->col); free(W->exps);
    mpz_clears(W->A, W->B, W->C, W->g, W->X, W->t, NULL);
    for (int l = 0; l < QS_MAXS; l++) mpz_clear(W->Bl[l]);
    fi_rellog_free(&W->log);
}

static void sieve_family_job(void* p) { sieve_family((QS*)p); }

static int qs_run(mpz_t d, const mpz_t n, FiCtx* c) {
    QS Q;
    memset(&Q, 0, sizeof Q);
    Q.N = n;
    mpz_init(Q.kN);
    unsigned k = fi_ks_multiplier(n);
    mpz_mul_ui(Q.kN, n, k);
    unsigned bits = (unsigned)mpz_sizeinbase(Q.kN, 2);
    unsigned want, nblk, lpmult;
    double slack;
    qs_params(bits, &want, &nblk, &lpmult, &slack);
    FI_TUNE_U("FI_QS_FB", want);
    FI_TUNE_U("FI_QS_BLK", nblk);
    FI_TUNE_U("FI_QS_LP", lpmult);
    FI_TUNE_D("FI_QS_SLACK", slack);

    /* Factor base: 2 (if usable), then primes with (kN/p) != -1. */
    size_t np;
    const uint32_t* P = fi_small_primes(&np);
    Q.fb = malloc((want + 2) * sizeof(uint32_t));
    Q.sqrtkn = malloc((want + 2) * sizeof(uint32_t));
    Q.logp = malloc(want + 2);
    Q.fb[0] = 0;
    uint32_t nfb = 0;
    int found = 0;
    for (size_t i = 0; i < np && nfb < want; i++) {
        uint32_t p = P[i];
        uint32_t r = (uint32_t)mpz_fdiv_ui(Q.kN, p);
        if (mpz_divisible_ui_p(n, p) && mpz_cmp_ui(n, p) != 0) { mpz_set_ui(d, p); found = 1; break; }
        if (p == 2) { Q.fb[++nfb] = 2; Q.sqrtkn[nfb] = r & 1; continue; }
        if (r != 0 && powmod32(r, (p - 1) / 2, p) != 1) continue;
        Q.fb[++nfb] = p;
        Q.sqrtkn[nfb] = sqrtmod32(r, p);
    }
    if (found) goto cleanup_fb;
    Q.nfb = nfb;
    Q.sieve_start = 1;
    unsigned small_p = bits < 120 ? 32u : QS_SMALL_P;
    FI_TUNE_U("FI_QS_SMALL", small_p);
    while (Q.sieve_start <= nfb && Q.fb[Q.sieve_start] < small_p) Q.sieve_start++;
    Q.med_end = Q.sieve_start;
    unsigned med_p = QS_MED_P;
    FI_TUNE_U("FI_QS_MED", med_p);
    while (Q.med_end <= nfb && Q.fb[Q.med_end] < med_p) Q.med_end++;
    if (nfb >= (1u << (32 - QS_BLOCK_BITS)) || Q.sieve_start + 4 > nfb) goto cleanup_fb;

    /* nblk is blocks per side; 0 means a single block in all (tiny inputs). */
    Q.nblocks = nblk ? nblk * 2 : 1;
    if (Q.nblocks > 64) Q.nblocks = 64;
    Q.L = Q.nblocks * QS_BLOCK;
    Q.M = Q.L / 2;
    uint32_t pmax = Q.fb[nfb];
    Q.lp_max = (uint64_t)pmax * lpmult;
    if (Q.lp_max > (uint64_t)pmax * pmax) Q.lp_max = (uint64_t)pmax * pmax;
    Q.pmax2 = (uint64_t)pmax * pmax;
    /* Double large primes pay from ~72 digits: cofactors up to lp_max^1.8
     * are split, and the threshold drops to let them through. */
    unsigned dlp = bits >= QS_DLP_BITS;
    double dlp_exp = 1.8, dlp_slack = 0.7;
    FI_TUNE_U("FI_QS_DLP", dlp);
    FI_TUNE_D("FI_QS_DLPEXP", dlp_exp);
    FI_TUNE_D("FI_QS_DLPSLACK", dlp_slack);
    Q.dlp_max = 0;
    if (dlp) {
        double m = pow((double)Q.lp_max, dlp_exp);
        Q.dlp_max = m >= 9.2e18 ? UINT64_C(9200000000000000000) : (uint64_t)m;
        slack += dlp_slack * log2((double)Q.lp_max);
    }

    /* Threshold: |g| ~ M sqrt(kN / 2) at the ends of the interval; accept
     * when the sieved log is within log2(lp_max) (+ slack for the unsieved
     * small primes) of that. */
    {
        long ex;
        double mant = mpz_get_d_2exp(&ex, Q.kN);
        double lg = log2((double)Q.M) + 0.5 * (log2(mant) + (double)ex) - 0.5;
        double T = lg - log2((double)Q.lp_max) - slack;
        if (T < 10) T = 10;
        Q.scale = T > 100.0 ? 100.0 / T : 1.0;
        Q.init = (uint8_t)(128 - (int)(T * Q.scale + 0.5));
        for (uint32_t i = 1; i <= nfb; i++)
            Q.logp[i] = (uint8_t)(log2((double)Q.fb[i]) * Q.scale + 0.5);
    }

    /* Lazily built tables the workers read must exist before they start. */
    {
        uint32_t lim;
        (void)fi_prime_bitmap(&lim);
        (void)fi_u64_squfof(15, 0);
        (void)fi_u64_is_prime(97);
    }
    int nw = bits >= QS_PAR_BITS ? fi_nthreads() : 1;
    QS* W = malloc((size_t)nw * sizeof(QS));
    for (int w = 0; w < nw; w++) worker_init(&W[w], &Q, fi_rand(c) ^ (uint64_t)(w + 1));

    FiRels* R = fi_rels_new(Q.fb, nfb + 1);
    size_t need = nfb + 1 + 64, extra = 0;
    size_t used_cap = 1 << 16, nused = 0;
    uint64_t* used = malloc(used_cap * sizeof(uint64_t));
    while (!found) {
        int give_up = 0;
        for (int w = 0; w < nw; w++)
            if (nused + 1 >= used_cap || !choose_A(&W[w], c, used, &nused)) { give_up = 1; break; }
        if (give_up) break;
        fi_run_parallel(sieve_family_job, W, sizeof(QS), nw);
        for (int w = 0; w < nw; w++) fi_rellog_replay(&W[w].log, R, n);
        /* Dependencies are guaranteed once relations outnumber the columns
         * actually present; most base primes above the median appear
         * rarely, so this fires well before nfb + 64. */
        if (fi_rels_count(R) >= need || fi_rels_count(R) >= fi_rels_active(R) + 64 + extra) {
            found = fi_rels_solve(R, n, d);
            if (!found) { need += 64; extra += 64; }
        }
    }
    free(used);
    fi_rels_free(R);
    for (int w = 0; w < nw; w++) worker_free(&W[w]);
    free(W);
cleanup_fb:
    free(Q.fb); free(Q.sqrtkn); free(Q.logp);
    mpz_clear(Q.kN);
    return found;
}

int fi_split_qs(mpz_t d, const mpz_t n, FiCtx* c) {
    if (mpz_sizeinbase(n, 2) <= 64) {
        uint64_t w = mth_mpz_get_u64(n), r = 0;
        uint64_t seed = w | 1;
        r = fi_u64_ecm(w, &seed, 256);
        if (!r) r = fi_u64_squfof(w, 0);
        if (r > 1 && r < w) { mth_mpz_set_u64(d, (uint64_t)r); return 1; }
        return 0;
    }
    return qs_run(d, n, c);
}
