/*
 * ff_hnf.c -- Hermite normal forms over K[u]: the direct elimination over any
 * FLINT generic field, and the certified multi-modular form over Q and number
 * fields.  See ff_hnf.h for the contract and the certificate.
 */
#include "ff_hnf.h"

#ifdef USE_FLINT
#include <stdlib.h>
#include <string.h>

#include <flint/fmpq.h>

/* gr_poly_xgcd arrived in FLINT 3.2 (see flint_bridge.c). */
#if defined(__FLINT_RELEASE) && __FLINT_RELEASE >= 30200
#define FF_GR_POLY_XGCD gr_poly_xgcd
#else
#define FF_GR_POLY_XGCD gr_poly_xgcd_euclidean
#endif

#define FF_AT(e, nc, i, j) ((e) + (i) * (nc) + (j))
#define FF_ENTRY(ptr, i, sz) ((gr_ptr)((char*)(ptr) + (i) * (sz)))

/* ------------------------------------------------------------------ */
/*  The direct elimination                                             */
/* ------------------------------------------------------------------ */

int ff_hnf_gr(gr_poly_struct* e, slong nr, slong nc, const gr_poly_struct* D, gr_ctx_t K) {
    slong sz = gr_ctx_sizeof_elem(K);
    int st = GR_SUCCESS;
    gr_poly_t g, s, t, a, b, q, r1, r2, tmp, tmp2;
    gr_poly_init(g, K); gr_poly_init(s, K); gr_poly_init(t, K); gr_poly_init(a, K);
    gr_poly_init(b, K); gr_poly_init(q, K); gr_poly_init(r1, K); gr_poly_init(r2, K);
    gr_poly_init(tmp, K); gr_poly_init(tmp2, K);
    gr_ptr lc = gr_heap_init(K);
    int ok = 1;
    for (slong j = 0; j < nc && ok; j++) {
        /* pivot: least degree in column j among rows j.. */
        slong piv = -1, best = 0;
        for (slong i = j; i < nr; i++) {
            slong len = FF_AT(e, nc, i, j)->length;
            if (len > 0 && (piv < 0 || len < best)) { piv = i; best = len; }
        }
        if (piv < 0) { ok = 0; break; }
        if (piv != j)
            for (slong k = 0; k < nc; k++) gr_poly_swap(FF_AT(e, nc, j, k), FF_AT(e, nc, piv, k), K);
        for (slong i = j + 1; i < nr && ok; i++) {
            if (FF_AT(e, nc, i, j)->length == 0) continue;
            st |= gr_poly_set(a, FF_AT(e, nc, j, j), K);
            st |= gr_poly_set(b, FF_AT(e, nc, i, j), K);
            st |= FF_GR_POLY_XGCD(g, s, t, a, b, K);      /* s a + t b = g */
            st |= gr_poly_divexact(a, a, g, K);           /* a/g */
            st |= gr_poly_divexact(b, b, g, K);           /* b/g */
            for (slong k = j; k < nc; k++) {
                gr_poly_struct* rj = FF_AT(e, nc, j, k);
                gr_poly_struct* ri = FF_AT(e, nc, i, k);
                st |= gr_poly_set(r1, rj, K);
                st |= gr_poly_set(r2, ri, K);
                st |= gr_poly_mul(tmp, s, r1, K);
                st |= gr_poly_mul(tmp2, t, r2, K);
                st |= gr_poly_add(rj, tmp, tmp2, K);
                st |= gr_poly_mul(tmp, a, r2, K);
                st |= gr_poly_mul(tmp2, b, r1, K);
                st |= gr_poly_sub(ri, tmp, tmp2, K);
                if (D && k > j) {
                    st |= gr_poly_rem(rj, rj, D, K);
                    st |= gr_poly_rem(ri, ri, D, K);
                }
            }
            if (st != GR_SUCCESS) ok = 0;
        }
        if (!ok) break;
        /* monic pivot */
        gr_poly_struct* pv = FF_AT(e, nc, j, j);
        if (pv->length == 0) { ok = 0; break; }
        st |= gr_inv(lc, FF_ENTRY(pv->coeffs, pv->length - 1, sz), K);
        for (slong k = j; k < nc; k++)
            st |= gr_poly_mul_scalar(FF_AT(e, nc, j, k), FF_AT(e, nc, j, k), lc, K);
        if (st != GR_SUCCESS) ok = 0;
    }
    /* rows below the form must have vanished (always so modulo D) */
    for (slong i = nc; i < nr && ok; i++)
        for (slong k = 0; k < nc; k++)
            if (FF_AT(e, nc, i, k)->length != 0) { ok = 0; break; }
    /* reduce above the pivots, column by column */
    for (slong j = 0; j < nc && ok; j++) {
        gr_poly_struct* pv = FF_AT(e, nc, j, j);
        for (slong i = 0; i < j; i++) {
            gr_poly_struct* x = FF_AT(e, nc, i, j);
            if (x->length < pv->length) continue;
            st |= gr_poly_divrem(q, tmp, x, pv, K);
            for (slong k = j; k < nc; k++) {
                st |= gr_poly_mul(tmp, q, FF_AT(e, nc, j, k), K);
                st |= gr_poly_sub(FF_AT(e, nc, i, k), FF_AT(e, nc, i, k), tmp, K);
            }
        }
        if (st != GR_SUCCESS) ok = 0;
    }
    gr_heap_clear(lc, K);
    gr_poly_clear(g, K); gr_poly_clear(s, K); gr_poly_clear(t, K); gr_poly_clear(a, K);
    gr_poly_clear(b, K); gr_poly_clear(q, K); gr_poly_clear(r1, K); gr_poly_clear(r2, K);
    gr_poly_clear(tmp, K); gr_poly_clear(tmp2, K);
    return ok;
}

/* ------------------------------------------------------------------ */
/*  The multi-modular Hermite form: a client of ff_mm_solve            */
/* ------------------------------------------------------------------ */

typedef struct {
    slong nr, nc;
    const ff_qtpoly_struct* in;    /* rows (nr*nc entries), then D            */
    ff_qtpoly_struct* H;           /* the answer, nc*nc                        */
    slong* start;                  /* first position of entry (i, k)          */
} ff_hnf_ctx;

typedef struct { gr_poly_struct* e; slong n; } ff_hnf_img;

/* positions of the coefficients of a Hermite form with diagonal degrees deg:
 * entry (i, k), k >= i, has deg[i] + 1 (k == i) or deg[k] (k > i) of them */
static slong ff_entry_len(const slong* deg, slong i, slong k) {
    return k < i ? 0 : (k == i ? deg[i] + 1 : deg[k]);
}

static slong ff_hnf_npos(void* vctx, const slong* deg) {
    ff_hnf_ctx* C = vctx;
    slong pos = 0;
    for (slong i = 0; i < C->nc; i++)
        for (slong k = 0; k < C->nc; k++) {
            C->start[i * C->nc + k] = pos;
            pos += ff_entry_len(deg, i, k);
        }
    return pos;
}

static void* ff_hnf_image(void* vctx, const ff_mm_reduced* R, ff_img_field* G, slong* deg) {
    ff_hnf_ctx* C = vctx;
    slong nr = C->nr, nc = C->nc, rows = nr + nc;
    ff_hnf_img* I = malloc(sizeof *I);
    I->n = rows * nc;
    I->e = malloc((size_t)I->n * sizeof(gr_poly_struct));
    for (slong i = 0; i < I->n; i++) gr_poly_init(I->e + i, G->K);
    nmod_poly_t scratch; nmod_poly_init(scratch, G->f->mod.n);
    gr_poly_t Dg; gr_poly_init(Dg, G->K);
    slong iD = nr * nc;
    int ok = ff_img_poly(Dg, R->red + R->off[iD], R->len[iD], G, scratch);
    ok = ok && Dg->length == C->in[iD].len;       /* D monic: its degree is kept */
    int st = GR_SUCCESS;
    for (slong i = 0; ok && i < nr * nc; i++) {
        ok = ff_img_poly(I->e + i, R->red + R->off[i], R->len[i], G, scratch);
        if (ok && I->e[i].length >= Dg->length) st |= gr_poly_rem(I->e + i, I->e + i, Dg, G->K);
    }
    for (slong k = 0; ok && k < nc; k++) st |= gr_poly_set(FF_AT(I->e, nc, nr + k, k), Dg, G->K);
    ok = ok && st == GR_SUCCESS && ff_hnf_gr(I->e, rows, nc, Dg, G->K);
    for (slong j = 0; ok && j < nc; j++) deg[j] = FF_AT(I->e, nc, j, j)->length - 1;
    gr_poly_clear(Dg, G->K);
    nmod_poly_clear(scratch);
    if (!ok) {
        for (slong i = 0; i < I->n; i++) gr_poly_clear(I->e + i, G->K);
        free(I->e); free(I);
        return NULL;
    }
    return I;
}

static void ff_hnf_extract(void* vctx, void* vimg, ff_img_field* G, const slong* deg, nmod_poly_struct* out) {
    ff_hnf_ctx* C = vctx;
    ff_hnf_img* I = vimg;
    slong nc = C->nc;
    for (slong i = 0; i < nc; i++)
        for (slong k = i; k < nc; k++) {
            const gr_poly_struct* x = FF_AT(I->e, nc, i, k);
            slong len = ff_entry_len(deg, i, k);
            for (slong d = 0; d < len; d++) {
                nmod_poly_struct* o = out + C->start[i * nc + k] + d;
                if (d < x->length) ff_img_get(o, FF_ENTRY(x->coeffs, d, G->sz), G);
                else nmod_poly_zero(o);
            }
        }
}

static void ff_hnf_img_clear(void* vctx, void* vimg, ff_img_field* G) {
    (void)vctx;
    ff_hnf_img* I = vimg;
    for (slong i = 0; i < I->n; i++) gr_poly_clear(I->e + i, G->K);
    free(I->e); free(I);
}

/* (a) the Hermite shape of H (given the degree profile) and (b) every
 * generator of span(rows) + D K[u]^nc reduces to 0 by H.  Exact. */
static int ff_hnf_check(const ff_qtpoly_struct* H, const slong* deg, const ff_qtpoly_struct* rows,
                        slong nr, slong nc, const ff_qtpoly_struct* D, ff_k0* Z) {
    gr_ctx_struct* K = Z->K;
    gr_poly_struct* h = malloc((size_t)(nc * nc) * sizeof(gr_poly_struct));
    gr_poly_struct* r = malloc((size_t)nc * sizeof(gr_poly_struct));
    for (slong i = 0; i < nc * nc; i++) gr_poly_init(h + i, K);
    for (slong k = 0; k < nc; k++) gr_poly_init(r + k, K);
    gr_poly_t q, rem, tmp; gr_poly_init(q, K); gr_poly_init(rem, K); gr_poly_init(tmp, K);
    int ok = 1, st = GR_SUCCESS;
    for (slong i = 0; i < nc * nc && ok; i++) ok = ff_k0_poly(h + i, H + i, Z);
    /* (a) shape */
    for (slong i = 0; i < nc && ok; i++)
        for (slong k = 0; k < nc && ok; k++) {
            const gr_poly_struct* x = FF_AT(h, nc, i, k);
            if (k < i) ok = x->length == 0;
            else if (k == i)
                ok = x->length == deg[i] + 1 &&
                     gr_is_one(FF_ENTRY(x->coeffs, x->length - 1, Z->sz), K) == T_TRUE;
            else ok = x->length <= deg[k];
        }
    /* (b) membership of the generators: the nr rows, then D e_k */
    for (slong g = 0; g < nr + nc && ok; g++) {
        for (slong k = 0; k < nc && ok; k++) {
            if (g < nr) ok = ff_k0_poly(r + k, FF_AT(rows, nc, g, k), Z);
            else if (k == g - nr) ok = ff_k0_poly(r + k, D, Z);
            else st |= gr_poly_zero(r + k, K);
        }
        for (slong j = 0; j < nc && ok; j++) {
            if (r[j].length == 0) continue;
            st |= gr_poly_divrem(q, rem, r + j, FF_AT(h, nc, j, j), K);
            if (st != GR_SUCCESS || rem->length != 0) { ok = 0; break; }
            for (slong k = j; k < nc; k++) {
                st |= gr_poly_mul(tmp, q, FF_AT(h, nc, j, k), K);
                st |= gr_poly_sub(r + k, r + k, tmp, K);
            }
            if (st != GR_SUCCESS || r[j].length != 0) ok = 0;
        }
    }
    gr_poly_clear(q, K); gr_poly_clear(rem, K); gr_poly_clear(tmp, K);
    for (slong i = 0; i < nc * nc; i++) gr_poly_clear(h + i, K);
    for (slong k = 0; k < nc; k++) gr_poly_clear(r + k, K);
    free(h); free(r);
    return ok && st == GR_SUCCESS;
}

static int ff_hnf_certify(void* vctx, const fmpq* cand, const slong* deg, ff_k0* Z) {
    ff_hnf_ctx* C = vctx;
    slong nc = C->nc, n = Z->isnf ? fmpq_poly_degree(Z->nf->pol) : 1;
    for (slong i = 0; i < nc * nc; i++) {
        slong len = ff_entry_len(deg, i / nc, i % nc);
        ff_qtpoly_set_length(C->H + i, len);
        for (slong d = 0; d < len; d++) {
            fmpq_poly_zero(C->H[i].c + d);
            slong base = (C->start[i] + d) * n;
            for (slong j = 0; j < n; j++) fmpq_poly_set_coeff_fmpq(C->H[i].c + d, j, cand + base + j);
        }
        ff_qtpoly_normalise(C->H + i);
    }
    return ff_hnf_check(C->H, deg, C->in, C->nr, nc, C->in + C->nr * nc, Z);
}

int ff_hnf_modular(ff_qtpoly_struct* H, const ff_qtpoly_struct* rows, slong nr, slong nc,
                   const ff_qtpoly_struct* D, const fmpq_poly_struct* mu, ff_mm_stats* stats) {
    if (stats) memset(stats, 0, sizeof *stats);
    if (nc < 1 || nr < 0 || D->len < 2 || !fmpq_poly_is_one(D->c + D->len - 1)) return 0;
    /* the inputs, contiguous (shallow copies: read only) */
    slong nin = nr * nc + 1;
    ff_qtpoly_struct* in = malloc((size_t)nin * sizeof(ff_qtpoly_struct));
    for (slong i = 0; i < nr * nc; i++) in[i] = rows[i];
    in[nr * nc] = *D;
    ff_hnf_ctx C = { nr, nc, in, H, malloc((size_t)(nc * nc) * sizeof(slong)) };
    ff_mm_problem P = { &C, in, nin, mu, nc, ff_hnf_image, ff_hnf_npos, ff_hnf_extract,
                        ff_hnf_img_clear, ff_hnf_certify };
    int ok = ff_mm_solve(&P, stats);
    if (!ok)
        for (slong i = 0; i < nc * nc; i++) ff_qtpoly_set_length(H + i, 0);
    free(C.start);
    free(in);
    return ok;
}

#endif /* USE_FLINT */
