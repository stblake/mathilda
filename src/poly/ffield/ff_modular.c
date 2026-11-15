/*
 * ff_modular.c -- the certified multi-modular driver of the ffield library
 * (see ff_modular.h).
 */
#include "ff_modular.h"

#ifdef USE_FLINT
#include <stdlib.h>
#include <string.h>

#include <flint/fmpz.h>
#include <flint/nmod_poly_factor.h>
#include <flint/fq_nmod.h>
#include <flint/ulong_extras.h>
#include <flint/nf_elem.h>

#define FF_ENTRY(ptr, i, sz) ((gr_ptr)((char*)(ptr) + (i) * (sz)))

/* ------------------------------------------------------------------ */
/*  K[u] elements as t-polynomial coefficient arrays                   */
/* ------------------------------------------------------------------ */

void ff_qtpoly_init(ff_qtpoly_struct* a) { a->len = 0; a->c = NULL; }

void ff_qtpoly_clear(ff_qtpoly_struct* a) {
    for (slong d = 0; d < a->len; d++) fmpq_poly_clear(a->c + d);
    free(a->c);
    a->c = NULL; a->len = 0;
}

void ff_qtpoly_set_length(ff_qtpoly_struct* a, slong len) {
    if (len < a->len) {
        for (slong d = len; d < a->len; d++) fmpq_poly_clear(a->c + d);
    } else if (len > a->len) {
        a->c = realloc(a->c, (size_t)len * sizeof(fmpq_poly_struct));
        for (slong d = a->len; d < len; d++) fmpq_poly_init(a->c + d);
    }
    if (len == 0) { free(a->c); a->c = NULL; }
    a->len = len;
}

void ff_qtpoly_normalise(ff_qtpoly_struct* a) {
    slong len = a->len;
    while (len > 0 && fmpq_poly_is_zero(a->c + len - 1)) len--;
    ff_qtpoly_set_length(a, len);
}

/* ------------------------------------------------------------------ */
/*  Reduction mod p and the factor fields                              */
/* ------------------------------------------------------------------ */

static slong g_ff_mm_max_primes = 2000;
void ff_mm_set_max_primes(slong m) { g_ff_mm_max_primes = m > 0 ? m : 2000; }

/* a (rational t-polynomial, p-integral) mod p */
static void ff_reduce_fmpq_poly(nmod_poly_t out, const fmpq_poly_t a, ulong p) {
    nmod_poly_zero(out);
    slong len = fmpq_poly_length(a);
    if (len == 0) return;
    ulong dinv = n_invmod(fmpz_fdiv_ui(fmpq_poly_denref(a), p), p);
    const fmpz* num = fmpq_poly_numref(a);
    for (slong i = len - 1; i >= 0; i--) {
        ulong c = fmpz_fdiv_ui(num + i, p);
        if (c) nmod_poly_set_coeff_ui(out, i, n_mulmod2(c, dinv, p));
    }
}

static void ff_img_field_init(ff_img_field* G, const nmod_poly_t f, ulong p) {
    nmod_poly_init(G->f, p);
    nmod_poly_set(G->f, f);
    G->deg1 = nmod_poly_degree(f) == 1;
    G->root = 0;
    if (G->deg1) {
        /* f monic: t + c, root -c */
        G->root = n_negmod(nmod_poly_get_coeff_ui(f, 0), p);
        gr_ctx_init_nmod(G->K, p);
    } else {
        /* gr_ctx_clear frees a context made from a reference (as gr_ctx_init_fq_nmod) */
        fq_nmod_ctx_struct* q = flint_malloc(sizeof(fq_nmod_ctx_struct));
        fq_nmod_ctx_init_modulus(q, f, "t");
        _gr_ctx_init_fq_nmod_from_ref(G->K, q);
    }
    G->sz = gr_ctx_sizeof_elem(G->K);
}

static void ff_img_field_clear(ff_img_field* G) {
    gr_ctx_clear(G->K);
    nmod_poly_clear(G->f);
}

void ff_img_set(gr_ptr elem, const nmod_poly_t a, ff_img_field* G, nmod_poly_t scratch) {
    if (G->deg1) {
        *(ulong*)elem = nmod_poly_evaluate_nmod(a, G->root);
    } else {
        nmod_poly_rem(scratch, a, G->f);
        fq_nmod_set_nmod_poly((fq_nmod_struct*)elem, scratch,
                              (const fq_nmod_ctx_struct*)GR_CTX_DATA_AS_PTR(G->K));
    }
}

void ff_img_get(nmod_poly_t out, gr_srcptr elem, ff_img_field* G) {
    if (G->deg1) {
        nmod_poly_zero(out);
        nmod_poly_set_coeff_ui(out, 0, *(const ulong*)elem);
    } else {
        fq_nmod_get_nmod_poly(out, (const fq_nmod_struct*)elem,
                              (const fq_nmod_ctx_struct*)GR_CTX_DATA_AS_PTR(G->K));
    }
}

int ff_img_poly(gr_poly_t P, const nmod_poly_struct* ap, slong len, ff_img_field* G,
                nmod_poly_t scratch) {
    int st = gr_poly_zero(P, G->K);
    if (len == 0) return st == GR_SUCCESS;
    gr_poly_fit_length(P, len, G->K);
    for (slong d = 0; d < len; d++) ff_img_set(FF_ENTRY(P->coeffs, d, G->sz), ap + d, G, scratch);
    _gr_poly_set_length(P, len, G->K);
    _gr_poly_normalise(P, G->K);
    return st == GR_SUCCESS;
}

static void ff_mm_reduced_init(ff_mm_reduced* R, const ff_qtpoly_struct* in, slong nin, ulong p) {
    R->nin = nin;
    R->off = malloc((size_t)(nin > 0 ? nin : 1) * sizeof(slong));
    R->len = malloc((size_t)(nin > 0 ? nin : 1) * sizeof(slong));
    R->total = 0;
    for (slong i = 0; i < nin; i++) { R->off[i] = R->total; R->len[i] = in[i].len; R->total += in[i].len; }
    R->red = malloc((size_t)(R->total > 0 ? R->total : 1) * sizeof(nmod_poly_struct));
    for (slong i = 0; i < nin; i++)
        for (slong d = 0; d < in[i].len; d++) {
            nmod_poly_init(R->red + R->off[i] + d, p);
            ff_reduce_fmpq_poly(R->red + R->off[i] + d, in[i].c + d, p);
        }
}

static void ff_mm_reduced_clear(ff_mm_reduced* R) {
    for (slong i = 0; i < R->total; i++) nmod_poly_clear(R->red + i);
    free(R->red); free(R->off); free(R->len);
}

/* ------------------------------------------------------------------ */
/*  K in characteristic 0                                              */
/* ------------------------------------------------------------------ */

void ff_k0_init(ff_k0* Z, const fmpq_poly_struct* mu) {
    Z->isnf = mu != NULL && fmpq_poly_degree(mu) >= 2;
    if (Z->isnf) {
        nf_init(Z->nf, mu);
        gr_ctx_init_nf(Z->K, mu);
    } else {
        gr_ctx_init_fmpq(Z->K);
    }
    Z->sz = gr_ctx_sizeof_elem(Z->K);
}

void ff_k0_clear(ff_k0* Z) {
    gr_ctx_clear(Z->K);
    if (Z->isnf) nf_clear(Z->nf);
}

int ff_k0_set(gr_ptr c, const fmpq_poly_t tp, ff_k0* Z) {
    if (Z->isnf) { nf_elem_set_fmpq_poly((nf_elem_struct*)c, tp, Z->nf); return 1; }
    fmpq_t q; fmpq_init(q);
    fmpq_poly_get_coeff_fmpq(q, tp, 0);
    int ok = gr_set_fmpq(c, q, Z->K) == GR_SUCCESS;
    fmpq_clear(q);
    return ok;
}

int ff_k0_poly(gr_poly_t P, const ff_qtpoly_struct* a, ff_k0* Z) {
    int st = gr_poly_zero(P, Z->K);
    if (a->len == 0) return st == GR_SUCCESS;
    gr_poly_fit_length(P, a->len, Z->K);
    int ok = 1;
    for (slong d = 0; d < a->len && ok; d++) ok = ff_k0_set(FF_ENTRY(P->coeffs, d, Z->sz), a->c + d, Z);
    _gr_poly_set_length(P, a->len, Z->K);
    _gr_poly_normalise(P, Z->K);
    return ok && st == GR_SUCCESS;
}

/* ------------------------------------------------------------------ */
/*  The driver                                                         */
/* ------------------------------------------------------------------ */

static void ff_lcm_den(fmpz_t L, const ff_qtpoly_struct* a) {
    for (slong d = 0; d < a->len; d++) fmpz_lcm(L, L, fmpq_poly_denref(a->c + d));
}

static slong ff_total(const slong* v, slong k) {
    slong s = 0;
    for (slong i = 0; i < k; i++) s += v[i];
    return s;
}

int ff_mm_solve(const ff_mm_problem* P, ff_mm_stats* stats) {
    if (stats) memset(stats, 0, sizeof *stats);
    slong np = P->nprof;
    if (np < 1) return 0;
    /* mu monic (mu = t over Q) */
    fmpq_poly_t mul; fmpq_poly_init(mul);
    if (P->mu) fmpq_poly_make_monic(mul, P->mu);
    else fmpq_poly_set_coeff_si(mul, 1, 1);
    slong n = fmpq_poly_degree(mul);

    fmpz_t L; fmpz_init_set_ui(L, 1);
    for (slong i = 0; i < P->nin; i++) ff_lcm_den(L, P->in + i);
    fmpz_lcm(L, L, fmpq_poly_denref(mul));

    slong* prof = calloc((size_t)np, sizeof(slong));      /* the accepted profile */
    slong* pprof = calloc((size_t)np, sizeof(slong));     /* this prime's         */
    slong* fprof = calloc((size_t)np, sizeof(slong));     /* this factor's        */
    slong npos = 0, nc = 0, nused = 0, mismatches = 0, next_try = 2;
    fmpz* acc = NULL;
    ulong* primes = NULL;          /* the primes used, and their residues: */
    ulong* rtab = NULL;            /* row k = residues mod primes[k]        */
    slong cap_used = 0;
    fmpq* cand = NULL;
    ulong* res = NULL;
    int have_cand = 0, done = 0;
    fmpz_t M, tmp; fmpz_init(M); fmpz_init(tmp);
    ff_k0 Z; int Z_ready = 0;

    ulong p = UWORD(1) << 62;
    slong max_primes = g_ff_mm_max_primes;
    if (P->height_mult > 0) {
        /* the largest numerator or denominator of the inputs, in bits */
        slong hin = 0;
        for (slong i = 0; i < P->nin; i++)
            for (slong d = 0; d < P->in[i].len; d++) {
                const fmpq_poly_struct* c = P->in[i].c + d;
                slong hb = fmpz_bits(fmpq_poly_denref(c));
                if (hb > hin) hin = hb;
                for (slong j = 0; j < fmpq_poly_length(c); j++) {
                    hb = fmpz_bits(fmpq_poly_numref(c) + j);
                    if (hb > hin) hin = hb;
                }
            }
        slong want = (2 * P->height_mult * (hin + 64)) / 60 + 50;
        if (want > 100000) want = 100000;
        if (want > max_primes) max_primes = want;
    }
    for (slong tried = 0; tried < 2 * max_primes + 50 && !done && nused < max_primes; tried++) {
        p = n_nextprime(p, 1);
        if (stats) stats->primes_tried++;
        if (fmpz_fdiv_ui(L, p) == 0) continue;
        nmod_poly_t mup, dmu, gg; nmod_poly_init(mup, p); nmod_poly_init(dmu, p); nmod_poly_init(gg, p);
        ff_reduce_fmpq_poly(mup, mul, p);
        nmod_poly_derivative(dmu, mup);
        nmod_poly_gcd(gg, mup, dmu);
        /* mu mod p squarefree: K/p is the product of the factor fields */
        int good = nmod_poly_degree(mup) == n && nmod_poly_degree(gg) == 0;
        nmod_poly_factor_t fac; nmod_poly_factor_init(fac);
        if (good) nmod_poly_factor(fac, mup);
        slong s = good ? fac->num : 0;
        ff_mm_reduced R; int R_ready = 0;
        ff_img_field* G = NULL;
        void** img = NULL;
        if (good) {
            ff_mm_reduced_init(&R, P->in, P->nin, p); R_ready = 1;
            G = malloc((size_t)s * sizeof(ff_img_field));
            img = calloc((size_t)s, sizeof(void*));
            for (slong f = 0; f < s; f++) ff_img_field_init(G + f, fac->p + f, p);
            for (slong f = 0; f < s && good; f++) {
                img[f] = P->image(P->ctx, &R, G + f, fprof);
                if (!img[f]) { good = 0; break; }
                if (f == 0) memcpy(pprof, fprof, (size_t)np * sizeof(slong));
                else if (memcmp(pprof, fprof, (size_t)np * sizeof(slong)) != 0) {
                    good = 0;                          /* the factors disagree: unlucky */
                    if (stats) stats->unlucky++;
                }
            }
        }
        /* the profile against the accepted one */
        if (good) {
            slong tp = ff_total(pprof, np), ta = ff_total(prof, np);
            int reset = nused == 0;
            if (!reset && tp < ta) reset = 1;                   /* earlier primes unlucky */
            else if (!reset && tp > ta) { good = 0; if (stats) stats->unlucky++; }
            else if (!reset && memcmp(pprof, prof, (size_t)np * sizeof(slong)) != 0) {
                if (stats) stats->unlucky++;
                if (++mismatches > nused + 2) reset = 1; else good = 0;
            }
            if (good && reset) {
                if (acc) { for (slong i = 0; i < nc; i++) fmpz_clear(acc + i); free(acc); }
                if (cand) { for (slong i = 0; i < nc; i++) fmpq_clear(cand + i); free(cand); }
                free(res);
                memcpy(prof, pprof, (size_t)np * sizeof(slong));
                npos = P->npos(P->ctx, prof);
                nc = npos * n;
                acc = malloc((size_t)(nc > 0 ? nc : 1) * sizeof(fmpz));
                cand = malloc((size_t)(nc > 0 ? nc : 1) * sizeof(fmpq));
                res = malloc((size_t)(nc > 0 ? nc : 1) * sizeof(ulong));
                for (slong i = 0; i < nc; i++) { fmpz_init(acc + i); fmpq_init(cand + i); }
                nused = 0; have_cand = 0; mismatches = 0; next_try = 2;
                free(primes); free(rtab); primes = NULL; rtab = NULL; cap_used = 0;
                fmpz_one(M);
            }
        }
        if (good) {
            /* residues: CRT over the factors to F_p[t]/(mu_p), t-coefficients */
            nmod_poly_t v, w, cof, inv;
            nmod_poly_init(v, p); nmod_poly_init(w, p); nmod_poly_init(cof, p); nmod_poly_init(inv, p);
            nmod_poly_struct* Ef = malloc((size_t)s * sizeof(nmod_poly_struct));
            nmod_poly_struct** ex = malloc((size_t)s * sizeof(nmod_poly_struct*));
            for (slong f = 0; f < s; f++) {
                nmod_poly_init(Ef + f, p);
                if (s > 1) {
                    nmod_poly_div(cof, mup, fac->p + f);
                    nmod_poly_rem(w, cof, fac->p + f);
                    nmod_poly_invmod(inv, w, fac->p + f);
                    nmod_poly_mulmod(Ef + f, cof, inv, mup);
                }
                ex[f] = malloc((size_t)(npos > 0 ? npos : 1) * sizeof(nmod_poly_struct));
                for (slong q = 0; q < npos; q++) nmod_poly_init(ex[f] + q, p);
                P->extract(P->ctx, img[f], G + f, prof, ex[f]);
            }
            for (slong q = 0; q < npos; q++) {
                if (s == 1) nmod_poly_set(v, ex[0] + q);
                else {
                    nmod_poly_zero(v);
                    for (slong f = 0; f < s; f++) {
                        nmod_poly_mulmod(w, ex[f] + q, Ef + f, mup);
                        nmod_poly_add(v, v, w);
                    }
                }
                for (slong j = 0; j < n; j++) res[q * n + j] = nmod_poly_get_coeff_ui(v, j);
            }
            for (slong f = 0; f < s; f++) {
                for (slong q = 0; q < npos; q++) nmod_poly_clear(ex[f] + q);
                free(ex[f]);
                nmod_poly_clear(Ef + f);
            }
            free(ex); free(Ef);
            nmod_poly_clear(v); nmod_poly_clear(w); nmod_poly_clear(cof); nmod_poly_clear(inv);

            /* a candidate that agrees with this image is certified */
            if (have_cand) {
                int agree = 1;
                for (slong i = 0; i < nc && agree; i++) {
                    ulong dm = fmpz_fdiv_ui(fmpq_denref(cand + i), p);
                    if (dm == 0) { agree = 0; break; }
                    ulong nm = fmpz_fdiv_ui(fmpq_numref(cand + i), p);
                    agree = n_mulmod2(nm, n_invmod(dm, p), p) == res[i];
                }
                if (agree) {
                    if (!Z_ready) { ff_k0_init(&Z, P->mu ? mul : NULL); Z_ready = 1; }
                    if (stats) stats->certificates++;
                    if (P->certify(P->ctx, cand, prof, &Z)) done = 1;
                }
                if (!done) have_cand = 0;
            }
            if (!done) {
                /* Keep the residues; combine them only at a reconstruction
                 * attempt, by one product-tree CRT per coefficient (fmpz_comb).
                 * Updating every coefficient on every prime was quadratic in
                 * the number of primes: 1.1 of 3.5 s on a 3600-prime xgcd. */
                if (nused == cap_used) {
                    cap_used = cap_used ? 2 * cap_used : 64;
                    primes = realloc(primes, (size_t)cap_used * sizeof(ulong));
                    rtab = realloc(rtab, (size_t)cap_used * (size_t)(nc > 0 ? nc : 1) * sizeof(ulong));
                }
                primes[nused] = p;
                memcpy(rtab + (size_t)nused * (size_t)nc, res, (size_t)nc * sizeof(ulong));
                if (nused == 0) fmpz_set_ui(M, p); else fmpz_mul_ui(M, M, p);
                nused++;
                if (!have_cand && nused >= next_try) {
                    int okr = 1;
                    if (nused == 1) {
                        for (slong i = 0; i < nc && okr; i++) {
                            fmpz_set_ui(acc + i, rtab[i]);
                            okr = fmpq_reconstruct_fmpz(cand + i, acc + i, M);
                        }
                    } else {
                        fmpz_comb_t comb; fmpz_comb_temp_t ct;
                        fmpz_comb_init(comb, primes, nused);
                        fmpz_comb_temp_init(ct, comb);
                        ulong* col = malloc((size_t)nused * sizeof(ulong));
                        for (slong i = 0; i < nc && okr; i++) {
                            for (slong k = 0; k < nused; k++) col[k] = rtab[(size_t)k * (size_t)nc + (size_t)i];
                            fmpz_multi_CRT_ui(acc + i, col, comb, ct, 0);
                            okr = fmpq_reconstruct_fmpz(cand + i, acc + i, M);
                        }
                        free(col);
                        fmpz_comb_temp_clear(ct);
                        fmpz_comb_clear(comb);
                    }
                    have_cand = okr;
                    next_try = nused + (nused / 2 > 1 ? nused / 2 : 1);
                }
            }
        }
        if (G) {
            for (slong f = 0; f < s; f++) {
                if (img[f]) P->img_clear(P->ctx, img[f], G + f);
                ff_img_field_clear(G + f);
            }
            free(G); free(img);
        }
        if (R_ready) ff_mm_reduced_clear(&R);
        nmod_poly_factor_clear(fac);
        nmod_poly_clear(mup); nmod_poly_clear(dmu); nmod_poly_clear(gg);
    }
    if (stats) stats->primes = nused;
    if (Z_ready) ff_k0_clear(&Z);
    free(primes); free(rtab);
    if (acc) { for (slong i = 0; i < nc; i++) fmpz_clear(acc + i); free(acc); }
    if (cand) { for (slong i = 0; i < nc; i++) fmpq_clear(cand + i); free(cand); }
    free(res);
    free(prof); free(pprof); free(fprof);
    fmpz_clear(M); fmpz_clear(tmp); fmpz_clear(L);
    fmpq_poly_clear(mul);
    return done;
}

#endif /* USE_FLINT */
