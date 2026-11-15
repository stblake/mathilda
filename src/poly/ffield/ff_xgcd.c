/*
 * ff_xgcd.c -- the certified multi-modular extended gcd (see ff_xgcd.h).
 */
#include "ff_xgcd.h"

#ifdef USE_FLINT
#include <stdlib.h>
#include <string.h>

#include <flint/fmpz_mpoly.h>
#include <flint/fmpq_mpoly.h>
#include <flint/nf_elem.h>

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
    const fmpq_poly_struct* mu;    /* monic, or NULL (Q)                      */
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

/* L a as an integer polynomial in (u, t), L the lcm of a's denominators. */
static void ff_qt_to_zmpoly(fmpz_mpoly_t out, fmpz_t L, const ff_qtpoly_struct* a,
                            const fmpz_mpoly_ctx_t ctx) {
    fmpz_one(L);
    for (slong d = 0; d < a->len; d++) fmpz_lcm(L, L, fmpq_poly_denref(a->c + d));
    fmpz_mpoly_zero(out, ctx);
    fmpz_t c; fmpz_init(c);
    ulong e[2];
    for (slong d = 0; d < a->len; d++) {
        const fmpq_poly_struct* cd = a->c + d;
        /* coefficient j of cd is num_j / den; times L: num_j * (L / den) */
        fmpz_t sc; fmpz_init(sc);
        fmpz_divexact(sc, L, fmpq_poly_denref(cd));
        for (slong j = 0; j < fmpq_poly_length(cd); j++) {
            if (fmpz_is_zero(fmpq_poly_numref(cd) + j)) continue;
            fmpz_mul(c, fmpq_poly_numref(cd) + j, sc);
            e[0] = (ulong)d; e[1] = (ulong)j;
            fmpz_mpoly_push_term_fmpz_ui(out, c, e, ctx);
        }
        fmpz_clear(sc);
    }
    fmpz_mpoly_sort_terms(out, ctx);
    fmpz_mpoly_combine_like_terms(out, ctx);
    fmpz_clear(c);
}

/* P == 0 modulo mu(t) (P an integer polynomial in (u, t); mu monic, or NULL
 * or of degree < 2 when K = Q: then P itself must vanish). */
static int ff_zmpoly_zero_mod_mu(const fmpz_mpoly_t P, const fmpq_poly_struct* mu,
                                 const fmpz_mpoly_ctx_t ctx) {
    int ok;
    if (fmpz_mpoly_is_zero(P, ctx)) ok = 1;
    else if (!mu || fmpq_poly_degree(mu) < 2) ok = 0;
    else {
        /* reduce modulo mu(t) over Q */
        fmpq_mpoly_ctx_t qctx; fmpq_mpoly_ctx_init(qctx, 2, ORD_LEX);
        fmpq_mpoly_t E, Mu, Qq, R;
        fmpq_mpoly_init(E, qctx); fmpq_mpoly_init(Mu, qctx); fmpq_mpoly_init(Qq, qctx); fmpq_mpoly_init(R, qctx);
        /* E = P: content 1 over the integer polynomial, then canonicalised */
        fmpz_mpoly_set(E->zpoly, P, qctx->zctx);
        fmpq_one(E->content);
        fmpq_mpoly_reduce(E, qctx);
        fmpq_t c; fmpq_init(c);
        ulong e[2] = { 0, 0 };
        for (slong j = 0; j <= fmpq_poly_degree(mu); j++) {
            fmpq_poly_get_coeff_fmpq(c, mu, j);
            if (fmpq_is_zero(c)) continue;
            e[1] = (ulong)j;
            fmpq_mpoly_push_term_fmpq_ui(Mu, c, e, qctx);
        }
        fmpq_mpoly_sort_terms(Mu, qctx); fmpq_mpoly_combine_like_terms(Mu, qctx);
        fmpq_mpoly_divrem(Qq, R, E, Mu, qctx);
        ok = fmpq_mpoly_is_zero(R, qctx);
        fmpq_clear(c);
        fmpq_mpoly_clear(E, qctx); fmpq_mpoly_clear(Mu, qctx); fmpq_mpoly_clear(Qq, qctx); fmpq_mpoly_clear(R, qctx);
        fmpq_mpoly_ctx_clear(qctx);
    }
    return ok;
}

/* s a + t b == g in K[u] = Q[t]/(mu)[u], decided over Z[u, t]: with each x
 * scaled to the integer polynomial X = L_x x, the identity times
 * L_s L_a L_t L_b L_g is  S A (L_t L_b L_g) + T B (L_s L_a L_g) - G (L_s L_a L_t L_b),
 * reduced modulo mu(t) (a monic divisor in t alone, so the multivariate
 * remainder is the normal form).  FLINT's integer multiplication instead of
 * number-field arithmetic that canonicalised a fraction at every step: 1.5 s
 * -> a few ms on 26000-digit cofactors. */
static int ff_xgcd_identity(const ff_qtpoly_struct* s, const ff_qtpoly_struct* a,
                            const ff_qtpoly_struct* t, const ff_qtpoly_struct* b,
                            const ff_qtpoly_struct* g, const fmpq_poly_struct* mu) {
    fmpz_mpoly_ctx_t ctx; fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);
    fmpz_mpoly_t S, A, T, B, G, P, Q;
    fmpz_mpoly_init(S, ctx); fmpz_mpoly_init(A, ctx); fmpz_mpoly_init(T, ctx);
    fmpz_mpoly_init(B, ctx); fmpz_mpoly_init(G, ctx); fmpz_mpoly_init(P, ctx); fmpz_mpoly_init(Q, ctx);
    fmpz_t Ls, La, Lt, Lb, Lg, f;
    fmpz_init(Ls); fmpz_init(La); fmpz_init(Lt); fmpz_init(Lb); fmpz_init(Lg); fmpz_init(f);
    ff_qt_to_zmpoly(S, Ls, s, ctx); ff_qt_to_zmpoly(A, La, a, ctx);
    ff_qt_to_zmpoly(T, Lt, t, ctx); ff_qt_to_zmpoly(B, Lb, b, ctx);
    ff_qt_to_zmpoly(G, Lg, g, ctx);
    fmpz_mpoly_mul(P, S, A, ctx);
    fmpz_mul(f, Lt, Lb); fmpz_mul(f, f, Lg);
    fmpz_mpoly_scalar_mul_fmpz(P, P, f, ctx);
    fmpz_mpoly_mul(Q, T, B, ctx);
    fmpz_mul(f, Ls, La); fmpz_mul(f, f, Lg);
    fmpz_mpoly_scalar_mul_fmpz(Q, Q, f, ctx);
    fmpz_mpoly_add(P, P, Q, ctx);
    fmpz_mul(f, Ls, La); fmpz_mul(f, f, Lt); fmpz_mul(f, f, Lb);
    fmpz_mpoly_scalar_mul_fmpz(Q, G, f, ctx);
    fmpz_mpoly_sub(P, P, Q, ctx);
    int ok = ff_zmpoly_zero_mod_mu(P, mu, ctx);
    fmpz_clear(Ls); fmpz_clear(La); fmpz_clear(Lt); fmpz_clear(Lb); fmpz_clear(Lg); fmpz_clear(f);
    fmpz_mpoly_clear(S, ctx); fmpz_mpoly_clear(A, ctx); fmpz_mpoly_clear(T, ctx);
    fmpz_mpoly_clear(B, ctx); fmpz_mpoly_clear(G, ctx); fmpz_mpoly_clear(P, ctx); fmpz_mpoly_clear(Q, ctx);
    fmpz_mpoly_ctx_clear(ctx);
    return ok;
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
    /* s a + t b = g, over Z[u, t] modulo mu */
    if (ok) ok = ff_xgcd_identity(C->s, C->in + 0, C->t, C->in + 1, C->g, C->mu);
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
    fmpq_poly_t mum; fmpq_poly_init(mum);
    if (mu) fmpq_poly_make_monic(mum, mu);
    ff_xgcd_ctx C = { a->len - 1, b->len - 1, in, g, s, t, mu ? mum : NULL };
    ff_mm_problem P = { &C, in, 2, mu, 1, (a->len - 1) + (b->len - 1), ff_xgcd_image, ff_xgcd_npos,
                        ff_xgcd_extract, ff_xgcd_img_clear, ff_xgcd_certify };
    int ok = ff_mm_solve(&P, stats);
    if (!ok) { ff_qtpoly_set_length(g, 0); ff_qtpoly_set_length(s, 0); ff_qtpoly_set_length(t, 0); }
    fmpq_poly_clear(mum);
    return ok;
}

/* ------------------------------------------------------------------ */
/*  The gcd alone (ff_gcd_modular)                                     */
/* ------------------------------------------------------------------ */

typedef struct {
    slong da, db;
    const ff_qtpoly_struct* in;
    ff_qtpoly_struct* g;
} ff_gcd_ctx;

static slong ff_gcd_npos(void* vctx, const slong* prof) { (void)vctx; return prof[0] + 1; }

static void* ff_gcd_image(void* vctx, const ff_mm_reduced* R, ff_img_field* G, slong* prof) {
    ff_gcd_ctx* C = vctx;
    nmod_poly_t scratch; nmod_poly_init(scratch, G->f->mod.n);
    gr_poly_t A, B; gr_poly_init(A, G->K); gr_poly_init(B, G->K);
    gr_poly_struct* Gp = malloc(sizeof(gr_poly_struct));
    gr_poly_init(Gp, G->K);
    int ok = ff_img_poly(A, R->red + R->off[0], R->len[0], G, scratch) &&
             ff_img_poly(B, R->red + R->off[1], R->len[1], G, scratch);
    /* leading coefficients survive: the image gcd then has degree >= the true one */
    ok = ok && A->length == C->da + 1 && B->length == C->db + 1;
    ok = ok && gr_poly_gcd(Gp, A, B, G->K) == GR_SUCCESS && Gp->length > 0 &&
         gr_poly_make_monic(Gp, Gp, G->K) == GR_SUCCESS;
    if (ok) prof[0] = Gp->length - 1;
    gr_poly_clear(A, G->K); gr_poly_clear(B, G->K);
    nmod_poly_clear(scratch);
    if (!ok) { gr_poly_clear(Gp, G->K); free(Gp); return NULL; }
    return Gp;
}

static void ff_gcd_extract(void* vctx, void* vimg, ff_img_field* G, const slong* prof, nmod_poly_struct* out) {
    (void)vctx;
    ff_xgcd_put(out, (const gr_poly_struct*)vimg, prof[0] + 1, G);
}

static void ff_gcd_img_clear(void* vctx, void* vimg, ff_img_field* G) {
    (void)vctx;
    gr_poly_clear((gr_poly_struct*)vimg, G->K);
    free(vimg);
}

/* Certificate: g monic of the image degree d, g | a and g | b exactly.  Then g
 * divides gcd(a, b), so deg gcd(a, b) >= d; and the subresultants of a and b
 * that vanish over K vanish under the reduction modulo (p, f) (the leading
 * coefficients survive it), so the image gcd, of degree d, has degree
 * >= deg gcd(a, b).  Hence g = gcd(a, b).  When d = 0 that degree argument
 * alone gives gcd = 1. */
static int ff_gcd_certify(void* vctx, const fmpq* cand, const slong* prof, ff_k0* Z) {
    ff_gcd_ctx* C = vctx;
    slong dg = prof[0], n = Z->isnf ? fmpq_poly_degree(Z->nf->pol) : 1;
    ff_xgcd_take(C->g, cand, 0, dg + 1, n);
    if (C->g->len != dg + 1 || !fmpq_poly_is_one(C->g->c + dg)) return 0;
    if (dg == 0) return 1;
    gr_ctx_struct* K = Z->K;
    gr_poly_t A, B, G, q, r;
    gr_poly_init(A, K); gr_poly_init(B, K); gr_poly_init(G, K); gr_poly_init(q, K); gr_poly_init(r, K);
    int ok = ff_k0_poly(A, C->in + 0, Z) && ff_k0_poly(B, C->in + 1, Z) && ff_k0_poly(G, C->g, Z);
    int st = GR_SUCCESS;
    if (ok) { st |= gr_poly_divrem(q, r, A, G, K); ok = st == GR_SUCCESS && r->length == 0; }
    if (ok) { st |= gr_poly_divrem(q, r, B, G, K); ok = st == GR_SUCCESS && r->length == 0; }
    gr_poly_clear(A, K); gr_poly_clear(B, K); gr_poly_clear(G, K); gr_poly_clear(q, K); gr_poly_clear(r, K);
    return ok;
}

int ff_gcd_modular(ff_qtpoly_struct* g, const ff_qtpoly_struct* a, const ff_qtpoly_struct* b,
                   const fmpq_poly_struct* mu, ff_mm_stats* stats) {
    if (stats) memset(stats, 0, sizeof *stats);
    if (a->len == 0 || b->len == 0) return 0;
    if (fmpq_poly_is_zero(a->c + a->len - 1) || fmpq_poly_is_zero(b->c + b->len - 1)) return 0;
    ff_qtpoly_struct in[2] = { *a, *b };       /* shallow, read only */
    ff_gcd_ctx C = { a->len - 1, b->len - 1, in, g };
    ff_mm_problem P = { &C, in, 2, mu, 1, 0, ff_gcd_image, ff_gcd_npos,
                        ff_gcd_extract, ff_gcd_img_clear, ff_gcd_certify };
    int ok = ff_mm_solve(&P, stats);
    if (!ok) ff_qtpoly_set_length(g, 0);
    return ok;
}

/* q = a / b over K[u] when b divides a exactly (1), else 0.  Division with
 * remainder in exact arithmetic over K. */
int ff_divexact(ff_qtpoly_struct* q, const ff_qtpoly_struct* a, const ff_qtpoly_struct* b,
                const fmpq_poly_struct* mu) {
    if (b->len == 0) return 0;
    fmpq_poly_t mum; fmpq_poly_init(mum);
    if (mu) fmpq_poly_make_monic(mum, mu);
    ff_k0 Z; ff_k0_init(&Z, mu ? mum : NULL);
    gr_ctx_struct* K = Z.K;
    gr_poly_t A, B, Q, R;
    gr_poly_init(A, K); gr_poly_init(B, K); gr_poly_init(Q, K); gr_poly_init(R, K);
    int ok = ff_k0_poly(A, a, &Z) && ff_k0_poly(B, b, &Z) &&
             gr_poly_divrem(Q, R, A, B, K) == GR_SUCCESS && R->length == 0;
    if (ok) {
        ff_qtpoly_set_length(q, Q->length);
        for (slong d = 0; d < Q->length && ok; d++) {
            gr_srcptr c = (const char*)Q->coeffs + d * Z.sz;
            if (Z.isnf) nf_elem_get_fmpq_poly(q->c + d, (const nf_elem_struct*)c, Z.nf);
            else fmpq_poly_set_fmpq(q->c + d, (const fmpq*)c);
        }
        ff_qtpoly_normalise(q);
    }
    gr_poly_clear(A, K); gr_poly_clear(B, K); gr_poly_clear(Q, K); gr_poly_clear(R, K);
    ff_k0_clear(&Z);
    fmpq_poly_clear(mum);
    return ok;
}

#endif /* USE_FLINT */
