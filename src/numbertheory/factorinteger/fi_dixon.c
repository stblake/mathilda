/*
 * fi_dixon.c — Dixon's random-squares method.
 *
 * Dixon draws random x and keeps those whose square residue z == x^2
 * (mod N) factors over a base of small primes; a GF(2) dependency among the
 * exponent vectors gives X^2 == Y^2 (mod N).  Textbook Dixon takes x
 * uniformly, so z is as large as N and almost never smooth.  Here each x is
 * drawn as floor(sqrt(kN)) + j with random k < 2^16 and random j < 2^20, so
 * z = x^2 - kN is only O(j sqrt(kN)) — the residues are square-root sized,
 * which is what makes the method usable at all, while the choices stay
 * random and every residue is trial divided (no sieve: that refinement is
 * the quadratic sieve).  Since k varies, every prime can divide some z, so
 * the base is all primes up to the bound.
 *
 * Shares batched factor-base trial division with early abort, the single
 * large-prime variation, linear algebra and square roots with CFRAC.
 *
 * Threads (fi_smooth.c): the main thread draws a batch of random words from
 * the context generator, the workers build and trial-divide the candidates,
 * and the hits are replayed in draw order; afterwards the generator is left
 * exactly where the serial loop would have left it.  Relations, solves, the
 * divisor and the random stream are therefore those of the serial method.
 */
#include "fi_internal.h"

#include <stdlib.h>

size_t fi_cfrac_fb_size(size_t digits);

/* Candidates per thread per region (a few ms: a timeout is taken between
 * regions) and the claim chunk.  Dixon's candidates cost about the same at
 * every size, so the threads pay at once (16x at 9 digits on 8 cores). */
#define DIXON_BATCH 2048
#define DIXON_CHUNK 16

typedef struct {
    mpz_srcptr n;
    uint64_t*  draw;        /* the batch's random words */
    /* replay */
    FiRels*    R;
    mpz_ptr    d;
    mpz_t      x;
    size_t     need;
    size_t     stop;        /* index of the candidate that ended the run */
} DixonState;

static int dixon_prep(void* user, size_t idx, mpz_t X, mpz_t v, mpz_t kn,
                      uint32_t* col, uint32_t* exps, size_t* cnt) {
    (void)col; (void)exps; (void)cnt;
    DixonState* D = (DixonState*)user;
    uint64_t r = D->draw[idx];
    unsigned long k = 1 + (unsigned long)(r & 0xffff);
    unsigned long j = 1 + (unsigned long)((r >> 16) & 0xfffff);
    mpz_mul_ui(kn, D->n, k);
    mpz_sqrt(X, kn);
    mpz_add_ui(X, X, j);
    mpz_mul(v, X, X);
    mpz_sub(v, v, kn);                      /* z = x^2 - kN > 0, z == x^2 (mod N) */
    return 1;
}

static int dixon_visit(void* user, const FiLogRec* r, const uint32_t* col, const uint32_t* exps) {
    DixonState* D = (DixonState*)user;
    uint64_t L = r->L1;
    D->stop = (size_t)r->tag;
    if (L > 1 && mth_mpz_divisible_u64_p(D->n, (uint64_t)L) && mth_mpz_cmp_u64(D->n, (uint64_t)L) > 0) {
        mth_mpz_set_u64(D->d, (uint64_t)L);     /* a prime of n as the large prime */
        return 1;
    }
    mpz_mod(D->x, r->X, D->n);
    fi_rels_add(D->R, D->x, col, exps, r->cnt, L, D->n);
    if (fi_rels_count(D->R) >= D->need) {
        if (fi_rels_solve(D->R, D->n, D->d)) return 1;
        D->need += 32;
    }
    return 0;
}

int fi_split_dixon(mpz_t d, const mpz_t n, FiCtx* ctx) {
    size_t digits = mpz_sizeinbase(n, 10);
    /* All primes, so a base about 1.5x CFRAC's covers the same QR primes. */
    size_t want = fi_cfrac_fb_size(digits) * 3 / 2;
    size_t np;
    const uint32_t* P = fi_small_primes(&np);
    if (want > np - 1) want = np - 1;
    uint32_t* fb = malloc((want + 2) * sizeof(uint32_t));
    fb[0] = 0;
    for (size_t i = 0; i < want; i++) fb[i + 1] = P[i];
    size_t nfb = want;

    int found = 0;
    for (size_t i = 1; i <= nfb && !found; i++)
        if (mpz_divisible_ui_p(n, fb[i]) && mpz_cmp_ui(n, fb[i]) != 0) { mpz_set_ui(d, fb[i]); found = 1; }

    unsigned per = DIXON_BATCH;
    FI_TUNE_U("FI_DIXON_BATCH", per);
    int nt = fi_nthreads();
    size_t batch = (size_t)per * (size_t)nt;

    DixonState D;
    D.n = n;
    D.d = d;
    D.R = fi_rels_new(fb, nfb + 1);
    D.need = nfb + 48;
    D.draw = malloc(batch * sizeof(uint64_t));
    mpz_init(D.x);
    FiSmooth* S = fi_smooth_new(nt, nfb);
    uint64_t pmax = fb[nfb];
    uint64_t lp_max = pmax * 128;
    size_t zbits = mpz_sizeinbase(n, 2) / 2 + 40;
    size_t abort_bits = (size_t)(0.70 * zbits) + 8;
    uint64_t budget = ctx->iterations ? ctx->iterations : (UINT64_C(1) << 30);
    for (uint64_t it = 0; it < budget && !found; ) {
        size_t m = budget - it < batch ? (size_t)(budget - it) : batch;
        uint64_t rng0 = ctx->rng;
        for (size_t i = 0; i < m; i++) D.draw[i] = fi_rand(ctx);
        fi_smooth_run(S, m, DIXON_CHUNK, dixon_prep, NULL, &D, fb, nfb, lp_max, abort_bits);
        found = fi_smooth_replay(S, dixon_visit, &D);
        if (found) {
            /* Leave the generator where the serial loop stops: after draw
             * number D.stop of this batch. */
            ctx->rng = rng0;
            for (size_t i = 0; i <= D.stop; i++) (void)fi_rand(ctx);
        }
        it += m;
    }
    fi_smooth_free(S);
    mpz_clear(D.x);
    free(D.draw);
    fi_rels_free(D.R);
    free(fb);
    return found;
}
