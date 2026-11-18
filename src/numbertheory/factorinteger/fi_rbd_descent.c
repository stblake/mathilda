/*
 * fi_rbd_descent.c — the descent kernels of Blake rational base descent.
 *
 * Two descents ([RBD] = rational_base_descent_revised.tex):
 *
 *   iterated (R1-R4)  Q <- floor(b Q / a) from N, testing gcd(Q - j, N) for
 *                     |j| <= J = ceil(2 + a/(a-b)).  The reference algorithm:
 *                     its truncation error E_n < a/(a-b) is what forces the
 *                     base-dependent window.
 *   exact (5.1)       the same recurrence on S N, S = 2^s >= 2a/(a-b); then
 *                     Qhat = floor(Qtilde / S) is within one of the exact
 *                     floor(N b^k / a^k), and the target c M is always one of
 *                     Qhat, ..., Qhat + 3: four candidates for every base.
 *
 * Cost structure of the exact descent, which is the default:
 *
 *   - Deferred testing.  No target c M (c <= cmax) can lie above
 *     T = sqrt(cmax N) + 4, so nothing is tested there.  Because the exact
 *     values are a closed form, the descent does not iterate down to T at
 *     all: one division floor(S N b^k0 / a^k0) jumps to the depth k0 just
 *     above it.  (When b^k0 would be enormous — a base very close to 1 —
 *     the descent instead takes macro-steps of j at once with a^j one word,
 *     which the error analysis allows.)  For a balanced N this skips half
 *     the descent outright.
 *   - Word bases.  With a and b single words every step is mpz_mul_ui +
 *     mpz_tdiv_q_ui on preallocated scratch: two linear passes, no heap.
 *   - Batched gcds.  The four candidates' product Q(Q+1)(Q+2)(Q+3) =
 *     y(y+2), y = Q^2 + 3Q, is folded into an accumulator mod N, and one gcd
 *     is taken per block of B steps.  Only a block whose gcd is non-trivial
 *     is replayed candidate by candidate, so the hits (and the first hit's
 *     (k, j)) are exactly those of the unbatched search.
 *
 * The iterated descent is kept faithful to the paper (no jump; batching only
 * on request) as the reference the tests compare against; record_all mode
 * on it is the truncation-divisibility hit record of [RBD] Section 6.
 */
#include "fi_rbd.h"
#include "expr.h"   /* mth_mpz_add_u64: mpz_add_ui is 32-bit on LLP64 */

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Hit records                                                         */
/* ------------------------------------------------------------------ */

void fi_rbd_hits_init(fi_rbd_hits_t* h) {
    h->v = NULL;
    h->n = h->cap = 0;
}

void fi_rbd_hits_clear(fi_rbd_hits_t* h) {
    for (size_t i = 0; i < h->n; i++) mpz_clear(h->v[i].g);
    free(h->v);
    fi_rbd_hits_init(h);
}

int fi_rbd_hits_add(fi_rbd_hits_t* h, const mpz_t g, uint64_t k, int64_t j) {
    for (size_t i = 0; i < h->n; i++)
        if (mpz_cmp(h->v[i].g, g) == 0) return 0;
    if (h->n == h->cap) {
        size_t nc = h->cap ? 2 * h->cap : 8;
        fi_rbd_hit_t* nv = realloc(h->v, nc * sizeof *nv);
        if (!nv) return 0;
        h->v = nv;
        h->cap = nc;
    }
    mpz_init_set(h->v[h->n].g, g);
    h->v[h->n].k = k;
    h->v[h->n].j = j;
    h->n++;
    return 1;
}

void fi_rbd_free_mpz_array(mpz_t* v, size_t n) {
    for (size_t i = 0; i < n; i++) mpz_clear(v[i]);
    free(v);
}

/* Factor refinement.  The base B is kept pairwise coprime; inserting x
 * either finds it already present, or splits the first element b sharing
 * d = gcd(b, x) > 1 into d, b/d, x/d (re-inserted from a stack), or appends
 * x.  Every value pushed divides N, and the total size strictly decreases
 * along each split, so it terminates; with h hits it is O(h^2) gcds, which
 * is noise next to the descents that produced them. */
mpz_t* fi_rbd_coprime_base(const fi_rbd_hits_t* h, const mpz_t N, size_t* count) {
    *count = 0;
    if (!h || h->n == 0) return NULL;
    size_t cap = 2 * h->n + 4, n = 0;
    mpz_t* B = malloc(cap * sizeof(mpz_t));
    size_t scap = 3 * h->n + 8, sn = 0;
    mpz_t* stk = malloc(scap * sizeof(mpz_t));
    if (!B || !stk) { free(B); free(stk); return NULL; }
    mpz_t d;
    mpz_init(d);

#define RBD_PUSH(x) do {                                                     \
        if (mpz_cmp_ui((x), 1) > 0) {                                        \
            if (sn == scap) { scap *= 2; stk = realloc(stk, scap * sizeof(mpz_t)); } \
            mpz_init_set(stk[sn++], (x));                                    \
        }                                                                    \
    } while (0)

    RBD_PUSH(N);
    for (size_t i = 0; i < h->n; i++) RBD_PUSH(h->v[i].g);
    while (sn) {
        mpz_t x;
        mpz_init_set(x, stk[--sn]);       /* move the top out of the stack */
        mpz_clear(stk[sn]);
        int placed = 0;
        for (size_t i = 0; i < n && !placed; i++) {
            mpz_gcd(d, B[i], x);
            if (mpz_cmp_ui(d, 1) == 0) continue;
            if (mpz_cmp(d, B[i]) == 0 && mpz_cmp(d, x) == 0) { placed = 1; break; }
            /* Remove B[i] (swap with the last element), then re-insert the
             * three refined pieces. */
            mpz_t b;
            mpz_init_set(b, B[i]);
            mpz_swap(B[i], B[n - 1]);
            mpz_clear(B[--n]);
            mpz_divexact(b, b, d);
            mpz_divexact(x, x, d);
            RBD_PUSH(d);
            RBD_PUSH(b);
            RBD_PUSH(x);
            mpz_clear(b);
            placed = 1;
        }
        if (!placed) {
            if (n == cap) { cap *= 2; B = realloc(B, cap * sizeof(mpz_t)); }
            mpz_init_set(B[n++], x);
        }
        mpz_clear(x);
    }
#undef RBD_PUSH
    mpz_clear(d);
    free(stk);
    *count = n;
    return B;
}

/* ------------------------------------------------------------------ */
/* Scratch                                                             */
/* ------------------------------------------------------------------ */

void fi_rbd_scratch_init(fi_rbd_scratch_t* s) {
    mpz_inits(s->Q, s->Qh, s->t, s->u, s->v, s->g, s->Neff, NULL);
    s->blk = NULL;
    s->blk_k = NULL;
    s->blk_cap = 0;
    s->ml = NULL;
    s->ml_cap = 0;
}

void fi_rbd_scratch_clear(fi_rbd_scratch_t* s) {
    mpz_clears(s->Q, s->Qh, s->t, s->u, s->v, s->g, s->Neff, NULL);
    free(s->ml);
    s->ml = NULL;
    s->ml_cap = 0;
    for (size_t i = 0; i < s->blk_cap; i++) mpz_clear(s->blk[i]);
    free(s->blk);
    free(s->blk_k);
    s->blk = NULL;
    s->blk_k = NULL;
    s->blk_cap = 0;
}

static int scratch_reserve(fi_rbd_scratch_t* s, size_t n) {
    if (n <= s->blk_cap) return 1;
    mpz_t* nb = realloc(s->blk, n * sizeof(mpz_t));
    if (!nb) return 0;
    s->blk = nb;
    uint64_t* nk = realloc(s->blk_k, n * sizeof(uint64_t));
    if (!nk) return 0;
    s->blk_k = nk;
    for (size_t i = s->blk_cap; i < n; i++) mpz_init(s->blk[i]);
    s->blk_cap = n;
    return 1;
}

/* ------------------------------------------------------------------ */
/* Window and ceiling                                                  */
/* ------------------------------------------------------------------ */

#define RBD_MAX_WINDOW (UINT64_C(1) << 24)

uint64_t fi_rbd_window(const mpz_t a, const mpz_t b) {
    mpz_t d, q;
    mpz_inits(d, q, NULL);
    mpz_sub(d, a, b);
    mpz_cdiv_q(q, a, d);                   /* ceil(a / (a - b)) */
    uint64_t J = 0;
    if (mpz_sgn(d) > 0 && mpz_cmp_ui(q, (unsigned long)(RBD_MAX_WINDOW - 2)) <= 0)
        J = (uint64_t)mpz_get_ui(q) + 2;
    mpz_clears(d, q, NULL);
    return J;
}

void fi_rbd_ceiling(mpz_t T, const mpz_t N, mpz_srcptr mult, const mpz_t cmax, uint64_t slack) {
    mpz_mul(T, N, cmax);
    if (mult) mpz_mul(T, T, mult);
    mpz_sqrt(T, T);
    mth_mpz_add_u64(T, T, slack);   /* not clamped to a 32-bit `long` (LLP64) */
}

/* ------------------------------------------------------------------ */
/* Descent state                                                       */
/* ------------------------------------------------------------------ */

typedef struct {
    mpz_srcptr N;
    const fi_rbd_opts_t* o;
    fi_rbd_scratch_t* s;
    fi_rbd_hits_t* hits;
    uint64_t J;            /* iterated window radius */
    size_t   B, nblk;      /* batch block length and fill */
    size_t   nbits;        /* bits of N */
    const fi_rbd_modctx_t* mc;   /* accumulator modulus (batch mode) */
    int      found, added;
    mpz_ptr  g;
    uint64_t hk;
    int64_t  hj;
} Desc;

/* Test one candidate cand = Q - j.  Returns 1 when the descent should stop. */
static int test_cand(Desc* D, const mpz_t cand, uint64_t k, int64_t j) {
    if (mpz_sgn(cand) == 0) return 0;      /* gcd(0, N) = N: trivial */
    mpz_gcd(D->s->g, cand, D->N);
    if (mpz_cmp_ui(D->s->g, 1) <= 0 || mpz_cmp(D->s->g, D->N) >= 0) return 0;
    if (!D->found) {
        mpz_set(D->g, D->s->g);
        D->hk = k;
        D->hj = j;
    }
    D->found = 1;
    if (!D->o->record_all) return 1;
    if (D->hits) D->added += fi_rbd_hits_add(D->hits, D->s->g, k, j);
    return 0;
}

/* Test the window around Q (depth k) candidate by candidate, in the
 * reference order: exact j = 0, -1, -2, -3; iterated j = 0, 1, -1, 2, -2, ...
 * (Q - j for j = 1 before Q + j, as in [RBD] Appendix A). */
static int test_window(Desc* D, const mpz_t Q, uint64_t k) {
    mpz_ptr c = D->s->v;
    if (D->o->exact) {
        for (int off = 0; off <= 3; off++) {
            mpz_add_ui(c, Q, (unsigned long)off);
            if (test_cand(D, c, k, -(int64_t)off)) return 1;
        }
        return 0;
    }
    if (test_cand(D, Q, k, 0)) return 1;
    for (uint64_t j = 1; j <= D->J; j++) {
        mpz_sub_ui(c, Q, (unsigned long)j);
        if (test_cand(D, c, k, (int64_t)j)) return 1;
        mpz_add_ui(c, Q, (unsigned long)j);
        if (test_cand(D, c, k, -(int64_t)j)) return 1;
    }
    return 0;
}

/* Fold the factor x into the accumulator (s->ml[0..n)), growing the fold
 * scratch behind it on demand.  Returns 0 on allocation failure, which the
 * caller turns into "replay this block": the hits stay exact either way. */
static int acc_fold(Desc* D, mpz_srcptr x) {
    fi_rbd_scratch_t* s = D->s;
    size_t n = (size_t)D->mc->n, xs = mpz_size(x);
    size_t need = n + fi_rbd_modacc_work(D->mc, xs);
    if (need > s->ml_cap) {
        mp_limb_t* nm = realloc(s->ml, need * sizeof(mp_limb_t));
        if (!nm) return 0;
        s->ml = nm;
        s->ml_cap = need;
    }
    fi_rbd_modacc_fold(D->mc, s->ml, mpz_limbs_read(x), (mp_size_t)xs, s->ml + n);
    return 1;
}

/* Fold the window product of Q into the accumulator.  A failed fold zeroes
 * the accumulator, so the block's gcd is N and it is replayed. */
static void acc_window(Desc* D, const mpz_t Q) {
    fi_rbd_scratch_t* s = D->s;
    size_t n = (size_t)D->mc->n;
    int ok = 1;
    if (D->o->exact) {
        if (mpz_sgn(Q) == 0) {             /* candidates 0..3: 0 is trivial */
            mpz_set_ui(s->t, 6);
            ok = acc_fold(D, s->t);
        } else {
            mpz_add_ui(s->u, Q, 3);
            mpz_mul(s->u, s->u, Q);        /* y = Q (Q + 3) */
            mpz_add_ui(s->t, s->u, 2);     /* y (y + 2) = Q (Q+1) (Q+2) (Q+3) */
            if (2 * mpz_size(s->u) <= n) {
                /* Small y: one product, one fold (same limb count, half the
                 * REDC set-up). */
                mpz_mul(s->t, s->t, s->u);
                ok = acc_fold(D, s->t);
            } else {
                ok = acc_fold(D, s->u) && acc_fold(D, s->t);
            }
        }
    } else {
        /* Iterated window: prod_{|j| <= J} |Q - j|, skipping the trivial 0,
         * folded whenever the partial product outgrows N. */
        mpz_abs(s->t, Q);
        if (mpz_sgn(s->t) == 0) mpz_set_ui(s->t, 1);
        for (uint64_t j = 1; j <= D->J && ok; j++) {
            for (int sg = 0; sg < 2; sg++) {
                if (sg) mpz_add_ui(s->u, Q, (unsigned long)j);
                else    mpz_sub_ui(s->u, Q, (unsigned long)j);
                if (mpz_sgn(s->u) == 0) continue;
                mpz_abs(s->u, s->u);
                mpz_mul(s->t, s->t, s->u);
                if (mpz_size(s->t) > n) {
                    if (!(ok = acc_fold(D, s->t))) break;
                    mpz_set_ui(s->t, 1);
                }
            }
        }
        if (ok) ok = acc_fold(D, s->t);
    }
    if (!ok) memset(s->ml, 0, n * sizeof(mp_limb_t));
}

/* Close the current block: one gcd, and a candidate-by-candidate replay
 * only when it is non-trivial.  Returns 1 when the descent should stop. */
static int flush(Desc* D) {
    if (D->nblk == 0) return 0;
    fi_rbd_scratch_t* s = D->s;
    mpz_t acc;
    mpz_roinit_n(acc, s->ml, D->mc->n);
    mpz_gcd(s->g, acc, D->N);
    int stop = 0;
    if (mpz_cmp_ui(s->g, 1) != 0) {
        for (size_t i = 0; i < D->nblk && !stop; i++) stop = test_window(D, s->blk[i], s->blk_k[i]);
    }
    D->nblk = 0;
    fi_rbd_modacc_one(D->mc, s->ml);
    return stop;
}

/* Test depth k with working value Q.  Returns 1 when the descent should
 * stop. */
static int test_depth(Desc* D, const mpz_t Q, uint64_t k) {
    if (!D->o->batch) return test_window(D, Q, k);
    mpz_set(D->s->blk[D->nblk], Q);
    D->s->blk_k[D->nblk] = k;
    D->nblk++;
    acc_window(D, Q);
    if (D->nblk == D->B) return flush(D);
    return 0;
}

/* log2 of a positive integer, as a double. */
static double lg2(const mpz_t x) {
    signed long e;
    double m = mpz_get_d_2exp(&e, x);
    return (double)e + log2(m);
}

/* One step Q <- floor(b Q / a) (Q >= 0). */
static inline void step(mpz_t Q, const mpz_t a, const mpz_t b, unsigned long aw,
                        unsigned long bw, int word) {
    if (word) {
        mpz_mul_ui(Q, Q, bw);
        mpz_tdiv_q_ui(Q, Q, aw);
    } else {
        mpz_mul(Q, Q, b);
        mpz_tdiv_q(Q, Q, a);
    }
}

/* Position the scaled exact descent just above the deferred ceiling.
 * On entry s->Q = S * Neff; on exit s->Q = Qtilde_k, *k the depth reached
 * (Qtilde_k still above T, so no tested depth is skipped).  Uses one exact
 * division when b^k0 stays affordable, else word macro-steps.  T is the
 * ceiling already scaled by S (a bound on Qtilde, not Qhat). */
static void exact_skip(fi_rbd_scratch_t* s, uint64_t* k, const mpz_t a, const mpz_t b,
                       unsigned long aw, unsigned long bw, int word, mpz_srcptr T) {
    if (!T || mpz_cmp(s->Q, T) <= 0) return;
    double la = lg2(a), lb = mpz_cmp_ui(b, 1) == 0 ? 0.0 : lg2(b);
    double lr = la - lb;                    /* log2(a/b) > 0 */
    if (!(lr > 0)) return;
    mpz_add_ui(s->t, T, 1);
    double lQ = lg2(s->Q), lT = lg2(s->t);
    double kf = floor((lQ - lT) / lr) - 2.0;
    if (kf < 1.0) return;
    size_t qbits = mpz_sizeinbase(s->Q, 2);
    if (kf * la <= 4.0 * (double)qbits + 65536.0 && kf < (double)ULONG_MAX) {
        /* Jump: Qtilde_k0 = floor(S Neff b^k0 / a^k0), exact. */
        unsigned long k0 = (unsigned long)kf;
        mpz_pow_ui(s->t, b, k0);
        mpz_mul(s->Q, s->Q, s->t);
        mpz_pow_ui(s->t, a, k0);
        mpz_tdiv_q(s->Q, s->Q, s->t);
        *k += (uint64_t)k0;
        return;
    }
    if (!word) return;
    /* Macro-steps: Q <- floor(Q b^m / a^m) with a^m one word.  Each
     * truncation contributes < 1 after scaling by (b/a)^(later steps), so the
     * accumulated error stays below a/(a-b) <= S/2 as for single steps. */
    unsigned long ma = aw, mb = bw;
    unsigned m = 1;
    while (ma <= ULONG_MAX / aw) { ma *= aw; mb *= bw; m++; }
    size_t tbits = mpz_sizeinbase(T, 2);
    double drop = m * lr;
    while ((double)mpz_sizeinbase(s->Q, 2) > (double)tbits + drop + 3.0) {
        mpz_mul_ui(s->Q, s->Q, mb);
        mpz_tdiv_q_ui(s->Q, s->Q, ma);
        *k += m;
    }
}

/* S = 2^sh >= 2a/(a-b): sh = bit length of ceil(2a/(a-b)). */
static unsigned scale_shift(const mpz_t a, const mpz_t b) {
    mpz_t d, q;
    mpz_inits(d, q, NULL);
    mpz_sub(d, a, b);
    mpz_mul_2exp(q, a, 1);
    mpz_cdiv_q(q, q, d);
    unsigned sh = (unsigned)mpz_sizeinbase(q, 2);
    mpz_clears(d, q, NULL);
    return sh < 2 ? 2 : sh;
}

int fi_rbd_descend(mpz_t g, uint64_t* kout, int64_t* jout, const mpz_t N, const mpz_t a,
                   const mpz_t b, const fi_rbd_opts_t* o, fi_rbd_hits_t* hits,
                   fi_rbd_scratch_t* s) {
    if (mpz_sgn(b) <= 0 || mpz_cmp(a, b) <= 0 || mpz_cmp_ui(N, 4) < 0) return 0;

    Desc D;
    D.N = N;
    D.o = o;
    D.s = s;
    D.hits = hits;
    D.J = 0;
    D.nbits = mpz_sizeinbase(N, 2);
    D.mc = NULL;
    D.B = D.nbits / 8;
    if (D.B < 16) D.B = 16;
    if (D.B > 256) D.B = 256;
    D.nblk = 0;
    D.found = D.added = 0;
    D.g = g;
    D.hk = 0;
    D.hj = 0;
    if (!o->exact) {
        D.J = fi_rbd_window(a, b);
        if (!D.J) return 0;
        /* A wide iterated window makes the product per step large; keep
         * the block short so a replay stays cheap. */
        if (D.J > 64 && D.B > 32) D.B = 32;
    }
    /* Batching needs the accumulator modulus: the caller's, else a local
     * one.  An even N (Montgomery needs an odd modulus) runs unbatched,
     * which finds exactly the same hits. */
    fi_rbd_opts_t ob;
    fi_rbd_modctx_t mloc;
    int own_mc = 0;
    if (o->batch) {
        if (o->mod) D.mc = o->mod;
        else if (fi_rbd_modctx_init(&mloc, N)) { D.mc = &mloc; own_mc = 1; }
        if (!D.mc || !scratch_reserve(s, D.B)) {
            if (own_mc) fi_rbd_modctx_clear(&mloc);
            ob = *o;
            ob.batch = 0;
            o = D.o = &ob;
            D.mc = NULL;
            own_mc = 0;
        } else {
            size_t need = (size_t)D.mc->n;
            if (need > s->ml_cap) {
                mp_limb_t* nm = realloc(s->ml, need * sizeof(mp_limb_t));
                if (!nm) { if (own_mc) fi_rbd_modctx_clear(&mloc); return 0; }
                s->ml = nm;
                s->ml_cap = need;
            }
            fi_rbd_modacc_one(D.mc, s->ml);
        }
    }

    int word = mpz_fits_ulong_p(a) && mpz_fits_ulong_p(b);
    unsigned long aw = word ? mpz_get_ui(a) : 0, bw = word ? mpz_get_ui(b) : 0;

    if (o->mult) mpz_mul(s->Neff, N, o->mult);
    else         mpz_set(s->Neff, N);

    uint64_t k = 0;
    int stop = 0;
    if (o->exact) {
        unsigned sh = scale_shift(a, b);
        mpz_mul_2exp(s->Q, s->Neff, sh);
        if (o->ceiling) {
            mpz_mul_2exp(s->u, o->ceiling, sh);  /* ceiling on Qtilde */
            exact_skip(s, &k, a, b, aw, bw, word, s->u);
        }
        for (;;) {
            step(s->Q, a, b, aw, bw, word);
            k++;
            mpz_tdiv_q_2exp(s->Qh, s->Q, sh);
            if (o->ceiling && mpz_cmp(s->Qh, o->ceiling) > 0) continue;
            if ((stop = test_depth(&D, s->Qh, k)) != 0) break;
            if (mpz_sgn(s->Qh) == 0) break;
        }
    } else {
        mpz_set(s->Q, s->Neff);
        for (;;) {                          /* R2 descend, R3 search, R4 loop */
            step(s->Q, a, b, aw, bw, word);
            k++;
            if (!(o->ceiling && mpz_cmp(s->Q, o->ceiling) > 0))
                if ((stop = test_depth(&D, s->Q, k)) != 0) break;
            if (mpz_cmp_ui(s->Q, 2) < 0) break;
        }
    }
    if (!stop && o->batch) flush(&D);
    D.nblk = 0;
    if (own_mc) fi_rbd_modctx_clear(&mloc);

    if (!D.found) return 0;
    if (kout) *kout = D.hk;
    if (jout) *jout = D.hj;
    if (o->record_all) return D.added > 0 ? D.added : 1;
    return 1;
}

size_t fi_rbd_band(mpz_t** Qs, uint64_t** ks, const mpz_t N, mpz_srcptr mult,
                   const mpz_t a, const mpz_t b, const mpz_t lo, const mpz_t hi) {
    *Qs = NULL;
    *ks = NULL;
    if (mpz_sgn(b) <= 0 || mpz_cmp(a, b) <= 0) return 0;
    fi_rbd_scratch_t s;
    fi_rbd_scratch_init(&s);
    int word = mpz_fits_ulong_p(a) && mpz_fits_ulong_p(b);
    unsigned long aw = word ? mpz_get_ui(a) : 0, bw = word ? mpz_get_ui(b) : 0;
    unsigned sh = scale_shift(a, b);
    if (mult) mpz_mul(s.Neff, N, mult);
    else      mpz_set(s.Neff, N);
    mpz_mul_2exp(s.Q, s.Neff, sh);
    mpz_mul_2exp(s.u, hi, sh);
    uint64_t k = 0;
    exact_skip(&s, &k, a, b, aw, bw, word, s.u);
    size_t n = 0, cap = 0;
    mpz_t* q = NULL;
    uint64_t* kk = NULL;
    for (;;) {
        step(s.Q, a, b, aw, bw, word);
        k++;
        mpz_tdiv_q_2exp(s.Qh, s.Q, sh);
        if (mpz_cmp(s.Qh, hi) > 0) continue;
        if (mpz_cmp(s.Qh, lo) < 0 || mpz_sgn(s.Qh) == 0) break;
        if (n == cap) {
            cap = cap ? 2 * cap : 64;
            q = realloc(q, cap * sizeof(mpz_t));
            kk = realloc(kk, cap * sizeof(uint64_t));
        }
        mpz_init_set(q[n], s.Qh);
        kk[n] = k;
        n++;
    }
    fi_rbd_scratch_clear(&s);
    *Qs = q;
    *ks = kk;
    return n;
}
