/*
 * ff_xgcd.c -- the certified multi-modular extended gcd (see ff_xgcd.h).
 */
#include "ff_xgcd.h"

#ifdef USE_FLINT
#include <stdlib.h>
#include <string.h>

/* gr_poly_xgcd arrived in FLINT 3.2 (see flint_bridge.c). */
#if defined(__FLINT_RELEASE) && __FLINT_RELEASE >= 30200
#define FF_GR_POLY_XGCD gr_poly_xgcd
#else
#define FF_GR_POLY_XGCD gr_poly_xgcd_euclidean
#endif

#define FF_ENTRY(ptr, i, sz) ((gr_ptr)((char*)(ptr) + (i) * (sz)))

typedef struct {
    slong da, db;                  /* deg a, deg b                            */
    const ff_qtpoly_struct* in;    /* a, b                                    */
    ff_qtpoly_struct *g, *s, *t;   /* the answer                              */
} ff_xgcd_ctx;

typedef struct { gr_poly_t G, S, T; } ff_xgcd_img;

/* positions: g (dg + 1 coefficients), s (db - dg), t (da - dg) */
static slong ff_xgcd_npos(void* vctx, const slong* prof) {
    ff_xgcd_ctx* C = vctx;
    slong dg = prof[0];
    return (dg + 1) + (C->db - dg) + (C->da - dg);
}

static void* ff_xgcd_image(void* vctx, const ff_mm_reduced* R, ff_img_field* G, slong* prof) {
    ff_xgcd_ctx* C = vctx;
    nmod_poly_t scratch; nmod_poly_init(scratch, G->f->mod.n);
    gr_poly_t A, B; gr_poly_init(A, G->K); gr_poly_init(B, G->K);
    ff_xgcd_img* I = malloc(sizeof *I);
    gr_poly_init(I->G, G->K); gr_poly_init(I->S, G->K); gr_poly_init(I->T, G->K);
    int ok = ff_img_poly(A, R->red + R->off[0], R->len[0], G, scratch) &&
             ff_img_poly(B, R->red + R->off[1], R->len[1], G, scratch);
    /* the degrees of a and b must survive the reduction */
    ok = ok && A->length == C->da + 1 && B->length == C->db + 1;
    int st = GR_SUCCESS;
    if (ok) {
        st |= FF_GR_POLY_XGCD(I->G, I->S, I->T, A, B, G->K);
        ok = st == GR_SUCCESS && I->G->length > 0;
    }
    if (ok) {
        gr_ptr lc = gr_heap_init(G->K);
        st |= gr_inv(lc, FF_ENTRY(I->G->coeffs, I->G->length - 1, G->sz), G->K);
        st |= gr_poly_mul_scalar(I->G, I->G, lc, G->K);
        st |= gr_poly_mul_scalar(I->S, I->S, lc, G->K);
        st |= gr_poly_mul_scalar(I->T, I->T, lc, G->K);
        gr_heap_clear(lc, G->K);
        slong dg = I->G->length - 1;
        prof[0] = dg;
        /* the normalised cofactors (unique, see ff_xgcd.h) */
        ok = st == GR_SUCCESS && !(dg == C->da && dg == C->db) &&
             I->S->length <= C->db - dg && I->T->length <= C->da - dg;
    }
    gr_poly_clear(A, G->K); gr_poly_clear(B, G->K);
    nmod_poly_clear(scratch);
    if (!ok) {
        gr_poly_clear(I->G, G->K); gr_poly_clear(I->S, G->K); gr_poly_clear(I->T, G->K);
        free(I);
        return NULL;
    }
    return I;
}

static void ff_xgcd_put(nmod_poly_struct* out, const gr_poly_struct* P, slong len, ff_img_field* G) {
    for (slong d = 0; d < len; d++)
        if (d < P->length) ff_img_get(out + d, FF_ENTRY(P->coeffs, d, G->sz), G);
        else nmod_poly_zero(out + d);
}

static void ff_xgcd_extract(void* vctx, void* vimg, ff_img_field* G, const slong* prof, nmod_poly_struct* out) {
    ff_xgcd_ctx* C = vctx;
    ff_xgcd_img* I = vimg;
    slong dg = prof[0];
    ff_xgcd_put(out, I->G, dg + 1, G);
    ff_xgcd_put(out + dg + 1, I->S, C->db - dg, G);
    ff_xgcd_put(out + dg + 1 + (C->db - dg), I->T, C->da - dg, G);
}

static void ff_xgcd_img_clear(void* vctx, void* vimg, ff_img_field* G) {
    (void)vctx;
    ff_xgcd_img* I = vimg;
    gr_poly_clear(I->G, G->K); gr_poly_clear(I->S, G->K); gr_poly_clear(I->T, G->K);
    free(I);
}

static void ff_xgcd_take(ff_qtpoly_struct* x, const fmpq* cand, slong base, slong len, slong n) {
    ff_qtpoly_set_length(x, len);
    for (slong d = 0; d < len; d++) {
        fmpq_poly_zero(x->c + d);
        for (slong j = 0; j < n; j++) fmpq_poly_set_coeff_fmpq(x->c + d, j, cand + (base + d) * n + j);
    }
    ff_qtpoly_normalise(x);
}

static int ff_xgcd_certify(void* vctx, const fmpq* cand, const slong* prof, ff_k0* Z) {
    ff_xgcd_ctx* C = vctx;
    slong dg = prof[0], n = Z->isnf ? fmpq_poly_degree(Z->nf->pol) : 1;
    ff_xgcd_take(C->g, cand, 0, dg + 1, n);
    ff_xgcd_take(C->s, cand, dg + 1, C->db - dg, n);
    ff_xgcd_take(C->t, cand, dg + 1 + (C->db - dg), C->da - dg, n);
    gr_ctx_struct* K = Z->K;
    gr_poly_t A, B, G, S, T, q, r, u;
    gr_poly_init(A, K); gr_poly_init(B, K); gr_poly_init(G, K); gr_poly_init(S, K);
    gr_poly_init(T, K); gr_poly_init(q, K); gr_poly_init(r, K); gr_poly_init(u, K);
    int ok = ff_k0_poly(A, C->in + 0, Z) && ff_k0_poly(B, C->in + 1, Z) &&
             ff_k0_poly(G, C->g, Z) && ff_k0_poly(S, C->s, Z) && ff_k0_poly(T, C->t, Z);
    int st = GR_SUCCESS;
    /* g monic of degree dg */
    ok = ok && G->length == dg + 1 && gr_is_one(FF_ENTRY(G->coeffs, dg, Z->sz), K) == T_TRUE;
    /* g | a, g | b */
    if (ok) { st |= gr_poly_divrem(q, r, A, G, K); ok = st == GR_SUCCESS && r->length == 0; }
    if (ok) { st |= gr_poly_divrem(q, r, B, G, K); ok = st == GR_SUCCESS && r->length == 0; }
    /* s a + t b = g */
    if (ok) {
        st |= gr_poly_mul(u, S, A, K);
        st |= gr_poly_mul(r, T, B, K);
        st |= gr_poly_add(u, u, r, K);
        ok = st == GR_SUCCESS && gr_poly_equal(u, G, K) == T_TRUE;
    }
    /* the degree bounds hold by construction (the positions) */
    gr_poly_clear(A, K); gr_poly_clear(B, K); gr_poly_clear(G, K); gr_poly_clear(S, K);
    gr_poly_clear(T, K); gr_poly_clear(q, K); gr_poly_clear(r, K); gr_poly_clear(u, K);
    return ok;
}

int ff_xgcd_modular(ff_qtpoly_struct* g, ff_qtpoly_struct* s, ff_qtpoly_struct* t,
                    const ff_qtpoly_struct* a, const ff_qtpoly_struct* b,
                    const fmpq_poly_struct* mu, ff_mm_stats* stats) {
    if (stats) memset(stats, 0, sizeof *stats);
    if (a->len == 0 || b->len == 0) return 0;
    if (fmpq_poly_is_zero(a->c + a->len - 1) || fmpq_poly_is_zero(b->c + b->len - 1)) return 0;
    ff_qtpoly_struct in[2] = { *a, *b };       /* shallow, read only */
    ff_xgcd_ctx C = { a->len - 1, b->len - 1, in, g, s, t };
    ff_mm_problem P = { &C, in, 2, mu, 1, ff_xgcd_image, ff_xgcd_npos, ff_xgcd_extract,
                        ff_xgcd_img_clear, ff_xgcd_certify };
    int ok = ff_mm_solve(&P, stats);
    if (!ok) { ff_qtpoly_set_length(g, 0); ff_qtpoly_set_length(s, 0); ff_qtpoly_set_length(t, 0); }
    return ok;
}

#endif /* USE_FLINT */
