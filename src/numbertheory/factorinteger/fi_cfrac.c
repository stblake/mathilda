/*
 * fi_cfrac.c — the continued-fraction method (Morrison & Brillhart, 1975).
 *
 * The partial quotients of the continued fraction of sqrt(kN) give
 * A_{i-1}^2 == (-1)^i Q_i (mod N) with 0 < Q_i < 2 sqrt(kN), far smaller
 * than the residues of Dixon's method, so a useful fraction of them factor
 * over a base of small primes.  Ingredients:
 *   - Knuth–Schroeppel multiplier k, chosen to maximise the expected
 *     contribution of small primes to Q_i (primes with (kN/p) = 1 divide
 *     some Q_i; p | k always does);
 *   - factor base of the primes with (kN/p) != -1, sized by digit count;
 *   - batched trial division (one bignum remainder per word-sized product
 *     of base primes) with early abort when the cofactor after the first
 *     quarter of the base is still too large to plausibly finish;
 *   - the single large-prime variation (partials pair up in fi_gf2.c);
 *   - GF(2) elimination and square roots, shared (fi_gf2.c).
 *
 * fi_fb_factor() is shared with Dixon's method.
 */
#include "fi_internal.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Factor v (> 0) over the factor base fb[1..nfb] (fb[0] is the sign column
 * and not a prime).  Writes the sorted (col, exp) list starting after any
 * entries already in col/exps (count in *cnt), and the large-prime cofactor
 * in *L (1 when v was smooth).  Returns 1 if v is smooth or a partial with
 * L <= lp_max, 0 otherwise.  `abort_bits`: after the first quarter of the
 * base, give up when the cofactor still has more bits than this. */
int fi_fb_factor(mpz_t v, const uint32_t* fb, size_t nfb, uint32_t* col, uint32_t* exps,
                 size_t* cnt, uint64_t* L, uint64_t lp_max, size_t abort_bits) {
    size_t k = *cnt;
    size_t quarter = nfb / 4 + 1;
    mp_bitcnt_t z = mpz_scan1(v, 0);
    size_t c = 1;
    if (fb[1] == 2) {
        if (z) { mpz_tdiv_q_2exp(v, v, z); col[k] = 1; exps[k++] = (uint32_t)z; }
        c = 2;
    }
    while (c <= nfb) {
        /* Batch base primes into a word-sized product. */
        uint64_t M = 1;
        size_t c0 = c;
        while (c <= nfb) {
            uint64_t h, l;
            fi_mul64(M, fb[c], &h, &l);
            if (h) break;
            M = l;
            c++;
        }
        uint64_t r = mth_mpz_fdiv_u64(v, (uint64_t)M);
        for (size_t j = c0; j < c; j++) {
            if (r % fb[j]) continue;
            uint32_t e = 0;
            while (mpz_divisible_ui_p(v, fb[j])) { mpz_divexact_ui(v, v, fb[j]); e++; }
            if (e) { col[k] = (uint32_t)j; exps[k++] = e; }
        }
        if (mpz_cmp_ui(v, 1) == 0) break;
        if (c0 <= quarter && c > quarter && mpz_sizeinbase(v, 2) > abort_bits) { *cnt = k; return 0; }
    }
    *cnt = k;
    if (mpz_cmp_ui(v, 1) == 0) { *L = 1; return 1; }
    if (mpz_sizeinbase(v, 2) <= 64) {
        uint64_t rem = mth_mpz_get_u64(v);
        if (rem <= lp_max) { *L = rem; return 1; }   /* < pmax^2: prime */
    }
    return 0;
}

/* Jacobi symbol (a/m) for odd m > 0, word arithmetic. */
static int jacobi_u32(uint32_t a, uint32_t m) {
    int t = 1;
    a %= m;
    while (a) {
        while (!(a & 1)) {
            a >>= 1;
            uint32_t r = m & 7;
            if (r == 3 || r == 5) t = -t;
        }
        uint32_t tmp = a; a = m; m = tmp;
        if ((a & 3) == 3 && (m & 3) == 3) t = -t;
        a %= m;
    }
    return m == 1 ? t : 0;
}

/* Knuth–Schroeppel multiplier (shared with the quadratic sieve): maximise
 * f(k) = -log(k)/2 + sum_p g(p, kN) log p over the odd primes p < 1000 and
 * p = 2, where g = 2/(p-1) when (kN/p) = 1 and 1/p when p | k.  (kN/p) is
 * computed as (k/p)(N/p), with (N/p) taken once per prime. */
unsigned fi_ks_multiplier(const mpz_t n) {
    static const unsigned K[] = { 1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23,
                                  26, 29, 30, 31, 33, 34, 35, 37, 38, 39, 41, 42, 43, 46, 47,
                                  51, 53, 55, 57, 58, 59, 61, 62, 65, 66, 67, 69, 70, 71, 73 };
    size_t np;
    const uint32_t* P = fi_small_primes(&np);
    enum { MAXP = 200 };                 /* odd primes below ~1230 */
    int legn[MAXP];
    double w[MAXP], wk[MAXP];
    uint32_t pr[MAXP];
    size_t m = 0;
    for (size_t i = 1; i < np && m < MAXP; i++) {
        pr[m] = P[i];
        legn[m] = jacobi_u32((uint32_t)mpz_fdiv_ui(n, P[i]), P[i]);
        w[m] = 2.0 * log((double)P[i]) / (P[i] - 1);
        wk[m] = log((double)P[i]) / P[i];
        m++;
    }
    unsigned n8 = (unsigned)mpz_fdiv_ui(n, 8);
    unsigned best = 1;
    double bs = -1e9;
    for (size_t t = 0; t < sizeof K / sizeof K[0]; t++) {
        unsigned k = K[t];
        double s = -0.5 * log((double)k);
        unsigned m8 = (n8 * k) % 8;
        s += (m8 == 1) ? 2 * log(2.0) : (m8 == 5) ? log(2.0) : 0.5 * log(2.0);
        for (size_t i = 0; i < m; i++) {
            if (k % pr[i] == 0) { s += wk[i]; continue; }
            if (legn[i] * jacobi_u32(k, pr[i]) == 1) s += w[i];
        }
        if (s > bs) { bs = s; best = k; }
    }
    return best;
}

/* Factor-base size by digit count (interpolated). */
size_t fi_cfrac_fb_size(size_t digits) {
    static const struct { size_t d, f; } T[] = {
        { 10, 40 }, { 20, 100 }, { 25, 160 }, { 30, 260 }, { 35, 420 }, { 40, 650 },
        { 45, 1000 }, { 50, 1500 }, { 55, 2200 }, { 60, 3200 }, { 65, 4500 }, { 70, 5800 }
    };
    size_t nt = sizeof T / sizeof T[0];
    if (digits <= T[0].d) return T[0].f;
    for (size_t i = 1; i < nt; i++)
        if (digits <= T[i].d)
            return T[i - 1].f + (T[i].f - T[i - 1].f) * (digits - T[i - 1].d) / (T[i].d - T[i - 1].d);
    return T[nt - 1].f;
}

/* Build fb[1..] = primes p with (kn/p) != -1; returns count. */
size_t fi_build_fb(const mpz_t kn, size_t want, uint32_t* fb) {
    size_t np, c = 1;
    const uint32_t* P = fi_small_primes(&np);
    for (size_t i = 0; i < np && c <= want; i++) {
        uint32_t p = P[i];
        if (p == 2) { fb[c++] = 2; continue; }
        if (mpz_kronecker_ui(kn, p) != -1) fb[c++] = p;
    }
    return c - 1;
}

/* ------------------------------------------------------------------ */
/* The expansion, in batches                                           */
/* ------------------------------------------------------------------ */

/* Candidates per thread per region at full size (a few ms: a timeout is
 * taken between regions), the claim chunk, and the first batch: batches
 * double from there, so a small n is not charged a large batch's overshoot,
 * and the threads pay for themselves from the smallest bignum inputs (equal
 * at 11 digits, 1.7x at 20 on 8 cores). */
#define CFRAC_BATCH 1024
#define CFRAC_CHUNK 64
#define CFRAC_FIRST 64

enum { CF_MORE = 0, CF_PERIOD, CF_STOP };

/* Consecutive terms of the expansion: candidate e is A^2 == (-1)^odd Q (mod N). */
typedef struct {
    size_t   n, ninit, cap;
    mpz_t*   Q;
    mpz_t*   A;
    uint8_t* odd;
    int      end;        /* CF_PERIOD: the last term closes the period;
                            CF_STOP: the expansion ended (Q = 0 or budget) */
} CfBatch;

/* The recurrence state.  Term i: P_i = q_{i-1} Q_{i-1} - P_{i-1},
 * Q_i = (kN - P_i^2) / Q_{i-1}, relation A_{i-1}^2 == (-1)^i Q_i (mod N);
 * then q_i = floor((g + P_i) / Q_i), A_i = q_i A_{i-1} + A_{i-2}. */
typedef struct {
    mpz_srcptr n;
    mpz_t      kn, g, P, Q, Qp, q, t, A, Ap;
    uint64_t   i, budget;
    size_t     want;     /* next batch size, set by the main thread */
    CfBatch*   out;      /* where side() writes the next batch */
} CfWalk;

static void cf_walk_reset(CfWalk* W, unsigned k) {
    mpz_mul_ui(W->kn, W->n, k);
    mpz_sqrt(W->g, W->kn);
    /* i = 0: P0 = 0, Q0 = 1, q0 = g, A_{-1} = 1, A_0 = g. */
    mpz_set_ui(W->P, 0);
    mpz_set_ui(W->Qp, 1);              /* Q_{i-1} */
    mpz_set(W->q, W->g);
    mpz_set_ui(W->Ap, 1);              /* A_{i-2} */
    mpz_mod(W->A, W->g, W->n);         /* A_{i-1} */
    W->i = 1;
}

static void cf_batch_reserve(CfBatch* B, size_t m) {
    if (m > B->cap) {
        B->Q = realloc(B->Q, m * sizeof(mpz_t));
        B->A = realloc(B->A, m * sizeof(mpz_t));
        B->odd = realloc(B->odd, m);
        B->cap = m;
    }
    while (B->ninit < m) { mpz_init(B->Q[B->ninit]); mpz_init(B->A[B->ninit]); B->ninit++; }
}

static void cf_batch_free(CfBatch* B) {
    for (size_t e = 0; e < B->ninit; e++) { mpz_clear(B->Q[e]); mpz_clear(B->A[e]); }
    free(B->Q); free(B->A); free(B->odd);
}

/* Run the expansion for up to W->want terms into W->out.  Stops after a term
 * with Q_i = 1 (the period closes: from there the terms repeat and every
 * dependency would be trivial) and before Q = 0 or the budget. */
static void cf_walk_fill(CfWalk* W) {
    CfBatch* B = W->out;
    cf_batch_reserve(B, W->want);
    B->n = 0;
    B->end = CF_MORE;
    while (B->n < W->want) {
        if (W->i >= W->budget) { B->end = CF_STOP; return; }
        mpz_mul(W->t, W->q, W->Qp);
        mpz_sub(W->P, W->t, W->P);
        mpz_mul(W->t, W->P, W->P);
        mpz_sub(W->t, W->kn, W->t);
        mpz_divexact(W->Q, W->t, W->Qp);
        if (mpz_sgn(W->Q) == 0) { B->end = CF_STOP; return; }
        mpz_set(B->Q[B->n], W->Q);
        mpz_set(B->A[B->n], W->A);
        B->odd[B->n] = (uint8_t)(W->i & 1);
        B->n++;
        mpz_add(W->t, W->g, W->P);
        mpz_fdiv_q(W->q, W->t, W->Q);
        mpz_mul(W->t, W->q, W->A);
        mpz_add(W->t, W->t, W->Ap);
        mpz_set(W->Ap, W->A);
        mpz_mod(W->A, W->t, W->n);
        mpz_swap(W->Qp, W->Q);
        W->i++;
        if (mpz_cmp_ui(W->Qp, 1) == 0) { B->end = CF_PERIOD; return; }
    }
}

/* ------------------------------------------------------------------ */
/* Collection and replay                                               */
/* ------------------------------------------------------------------ */

typedef struct {
    mpz_srcptr n;
    mpz_ptr    d;
    CfBatch*   cur;      /* being trial divided (read-only in the region) */
    FiRels*    R;
    size_t     need;
    unsigned   trivial;
    CfWalk*    W;
} CfState;

/* Side task: worker 0 extends the expansion into the next batch. */
static void cf_side(void* user) { cf_walk_fill(((CfState*)user)->W); }

static int cf_prep(void* user, size_t idx, mpz_t X, mpz_t v, mpz_t tmp,
                   uint32_t* col, uint32_t* exps, size_t* cnt) {
    (void)X; (void)tmp;
    const CfBatch* B = ((CfState*)user)->cur;
    mpz_set(v, B->Q[idx]);
    if (B->odd[idx]) { col[0] = 0; exps[0] = 1; *cnt = 1; }
    return 1;
}

/* Replay one hit, exactly as the serial loop handles it.  1: found;
 * 2: restart on a new multiplier (repeated all-trivial solves). */
static int cf_visit(void* user, const FiLogRec* r, const uint32_t* col, const uint32_t* exps) {
    CfState* C = (CfState*)user;
    uint64_t L = r->L1;
    /* On a small n the large-prime bound can exceed a prime of n: such an L
     * is itself the answer. */
    if (L > 1 && mth_mpz_divisible_u64_p(C->n, (uint64_t)L) && mth_mpz_cmp_u64(C->n, (uint64_t)L) > 0) {
        mth_mpz_set_u64(C->d, (uint64_t)L);
        return 1;
    }
    fi_rels_add(C->R, C->cur->A[r->tag], col, exps, r->cnt, L, C->n);
    if (fi_rels_count(C->R) >= C->need) {
        if (fi_rels_solve(C->R, C->n, C->d)) return 1;
        C->need += 32;                      /* every dependency trivial: gather more */
        C->trivial++;
    }
    /* Repeated all-trivial solves mean the same as a closed period: for some
     * (N, k) every relation of the expansion lies in one square class mod
     * each prime (seen with k = 2 and both primes == 7 mod 8). */
    return C->trivial >= 3 ? 2 : 0;
}

/* Threads: the continued-fraction recurrence is sequential, but it is cheap
 * next to trial division, so the expansion is pipelined.  In each region
 * worker 0 first extends the expansion into the next batch (cf_walk_fill as
 * the side task) and then joins the others trial dividing the current one;
 * the main thread replays the current batch's hits in term order.  A restart
 * discards the batch generated ahead.  Terms, relations, restarts and the
 * divisor are those of the serial loop for any thread count. */
int fi_split_cfrac(mpz_t d, const mpz_t n, FiCtx* ctx) {
    size_t digits = mpz_sizeinbase(n, 10);
    unsigned k = fi_ks_multiplier(n), k0 = k, ki = 0;
    CfWalk W;
    W.n = n;
    mpz_inits(W.kn, W.g, W.P, W.Q, W.Qp, W.q, W.t, W.A, W.Ap, NULL);
    W.budget = ctx->iterations ? ctx->iterations : (UINT64_C(1) << 32);
    cf_walk_reset(&W, k);

    size_t want = fi_cfrac_fb_size(digits);
    uint32_t* fb = malloc((want + 2) * sizeof(uint32_t));
    fb[0] = 0;
    size_t nfb = fi_build_fb(W.kn, want, fb);
    /* A base prime dividing n itself is a factor. */
    int found = 0;
    for (size_t i = 1; i <= nfb && !found; i++)
        if (mpz_divisible_ui_p(n, fb[i]) && mpz_cmp_ui(n, fb[i]) != 0) { mpz_set_ui(d, fb[i]); found = 1; }

    unsigned per = CFRAC_BATCH;
    FI_TUNE_U("FI_CFRAC_BATCH", per);
    int nt = fi_nthreads();
    size_t full = (size_t)per * (size_t)nt, size = CFRAC_FIRST < full ? CFRAC_FIRST : full;

    CfState C;
    C.n = n;
    C.d = d;
    C.R = fi_rels_new(fb, nfb + 1);
    C.need = nfb + 1 + 48;
    C.trivial = 0;
    C.W = &W;
    CfBatch batch[2] = { { 0 }, { 0 } };
    C.cur = &batch[0];
    FiSmooth* S = fi_smooth_new(nt, want);
    uint64_t pmax = fb[nfb];
    uint64_t lp_max = pmax * 128 < pmax * pmax ? pmax * 128 : pmax * pmax;
    size_t qbits = mpz_sizeinbase(W.g, 2) + 1;
    size_t abort_bits = (size_t)(0.70 * qbits) + 8;

    W.want = size;
    W.out = C.cur;
    if (!found) cf_walk_fill(&W);
    while (!found) {
        CfBatch* nxt = C.cur == &batch[0] ? &batch[1] : &batch[0];
        size = 2 * size < full ? 2 * size : full;
        W.want = size;
        W.out = nxt;
        fi_smooth_run(S, C.cur->n, CFRAC_CHUNK, cf_prep, C.cur->end == CF_MORE ? cf_side : NULL,
                      &C, fb, nfb, lp_max, abort_bits);
        int rc = fi_smooth_replay(S, cf_visit, &C);
        if (rc == 1) { found = 1; break; }
        if (rc == 2 || C.cur->end == CF_PERIOD) {
            /* Restart on a new multiplier, with a fresh base and relation
             * store: the old relations all lie in the degenerate class, and
             * elimination would keep returning dependencies built from them. */
            static const unsigned KS[] = { 1, 3, 5, 7, 11, 13, 15, 17, 19, 21, 23, 29, 31, 33, 35, 37 };
            C.trivial = 0;
            do { k = KS[ki++ % (sizeof KS / sizeof KS[0])]; } while (k == k0 && ki < 64);
            if (ki >= 64) break;
            cf_walk_reset(&W, k);
            nfb = fi_build_fb(W.kn, want, fb);
            fi_rels_free(C.R);
            C.R = fi_rels_new(fb, nfb + 1);
            C.need = nfb + 1 + 48;
            pmax = fb[nfb];
            lp_max = pmax * 128 < pmax * pmax ? pmax * 128 : pmax * pmax;
            for (size_t t = 1; t <= nfb && !found; t++)
                if (mpz_divisible_ui_p(n, fb[t]) && mpz_cmp_ui(n, fb[t]) != 0) { mpz_set_ui(d, fb[t]); found = 1; }
            if (found || mpz_perfect_square_p(W.kn)) break;
            W.want = size;
            W.out = C.cur;
            cf_walk_fill(&W);
            continue;
        }
        if (C.cur->end == CF_STOP) break;
        C.cur = nxt;
    }
    fi_smooth_free(S);
    cf_batch_free(&batch[0]);
    cf_batch_free(&batch[1]);
    fi_rels_free(C.R);
    free(fb);
    mpz_clears(W.kn, W.g, W.P, W.Q, W.Qp, W.q, W.t, W.A, W.Ap, NULL);
    return found;
}
