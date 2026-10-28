/*
 * gb_flint.c -- Buchberger's algorithm on FLINT arithmetic.
 *
 * The engine behind GroebnerBasis, Eliminate, Solve's and NSolve's
 * polynomial systems, and the GF(p) Modulus path.  gb_buchberger,
 * gb_groebner_walk (groebner.c, groebnerwalk.c) and gfp_buchberger (gbmod.c)
 * try it first and keep their own engines for what it declines, so every
 * caller of those gets it without a change of its own.
 *
 * Polynomials are fmpz_mpoly (over Q, kept primitive with a positive leading
 * coefficient: no rational coefficient growth) or nmod_mpoly (over GF(p),
 * kept monic) in a FLINT context of the order -- ORD_LEX, ORD_DEGLEX,
 * ORD_DEGREVLEX, variable 0 most significant, which are exactly the
 * GBPoly orders Lexicographic, DegreeLexicographic (the weight matrix the
 * named order builds) and DegreeReverseLexicographic.  Elimination block
 * orders and other weight matrices decline.
 *
 * The algorithm:
 *  - the inputs are reduced against the basis as they are added;
 *  - pairs are chosen by the SUGAR strategy (smallest sugar degree, then the
 *    smallest lcm in the order), which is what keeps Lex on katsura from the
 *    intermediate degree blow-up that the normal strategy suffers;
 *  - the pair set is maintained with the Gebauer-Moeller criteria (the chain
 *    criterion on the new pairs, the product criterion, and criterion B on
 *    the old ones), and a basis element whose leading monomial the new one
 *    divides leaves the reducer set;
 *  - every S-polynomial is brought to full normal form by the reducer set,
 *    the shortest applicable reducer first: over Z fraction-free with the
 *    content divided out as the coefficients grow (gf_reduce_z), over GF(p)
 *    in one heap division (nmod_mpoly_divrem_ideal); then made primitive
 *    (monic);
 *  - the final reducer set is minimal, so one pass of tail reduction makes
 *    it the reduced basis; it is returned primitive over Z (monic over
 *    GF(p)), sorted ascending by leading monomial -- the form gb_buchberger
 *    and gfp_buchberger return.
 *
 * TimeConstrained: the deadline is polled once per pair with the
 * non-jumping tc_deadline_passed; on expiry the FLINT state is freed and only
 * then tc_check_deadline() unwinds, so an aborted run leaks nothing.
 *
 * MATHILDA_GB_FLINT=off makes every entry decline (A/B against the classical
 * engines).
 */
#include "gb_flint.h"

#ifdef USE_FLINT

#include <stdlib.h>
#include <string.h>

#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
#include <flint/nmod_vec.h>
#include <flint/nmod.h>
#include <flint/nmod_mpoly.h>
#include <flint/ulong_extras.h>

#include "core.h"           /* tc_deadline_passed, tc_check_deadline */

static int gf_enabled(void) {
    static int on = -1;
    if (on < 0) {
        const char* e = getenv("MATHILDA_GB_FLINT");
        on = !(e && strcmp(e, "off") == 0);
    }
    return on;
}

/* ------------------------------------------------------------------ */
/*  The monomial order                                                 */
/* ------------------------------------------------------------------ */

/* The FLINT ordering of a GBPoly order, or -1 when there is none. */
static int gf_flint_order(GBOrder order, const GBWeightMatrix* wmat, int n_vars) {
    if (order == GB_ORDER_LEX) return ORD_LEX;
    if (order == GB_ORDER_GREVLEX) return ORD_DEGREVLEX;
    if (order != GB_ORDER_MATRIX || !wmat || wmat->n_vars != n_vars) return -1;
    /* DegreeLexicographic reaches the engine as the weight matrix the named
     * order builds; recognise exactly that matrix. */
    int ds, rev, nl; GBOrder hint;
    if (!gb_classify_named_order("DegreeLexicographic", &ds, &rev, &nl, &hint)) return -1;
    int rows = 0;
    int64_t* w = gb_build_order_matrix(n_vars, ds, rev, nl, &rows);
    int same = w && rows == wmat->n_rows &&
               memcmp(w, wmat->w, (size_t)rows * (size_t)n_vars * sizeof(int64_t)) == 0;
    free(w);
    return same ? ORD_DEGLEX : -1;
}

/* Descending comparison of exponent vectors in the order: < 0 when a is the
 * larger monomial (FLINT's term order). */
static int gf_cmp(const ulong* a, const ulong* b, slong nv, int ord) {
    if (ord != ORD_LEX) {
        ulong da = 0, db = 0;
        for (slong v = 0; v < nv; v++) { da += a[v]; db += b[v]; }
        if (da != db) return da > db ? -1 : 1;
    }
    if (ord == ORD_DEGREVLEX) {
        for (slong v = nv - 1; v >= 0; v--)
            if (a[v] != b[v]) return a[v] < b[v] ? -1 : 1;
        return 0;
    }
    for (slong v = 0; v < nv; v++)
        if (a[v] != b[v]) return a[v] > b[v] ? -1 : 1;
    return 0;
}

static int gf_divides(const ulong* a, const ulong* b, slong nv) {     /* a | b */
    for (slong v = 0; v < nv; v++) if (a[v] > b[v]) return 0;
    return 1;
}

static int gf_coprime(const ulong* a, const ulong* b, slong nv) {
    for (slong v = 0; v < nv; v++) if (a[v] && b[v]) return 0;
    return 1;
}

static void gf_lcm(ulong* l, const ulong* a, const ulong* b, slong nv) {
    for (slong v = 0; v < nv; v++) l[v] = a[v] > b[v] ? a[v] : b[v];
}

static ulong gf_deg(const ulong* a, slong nv) {
    ulong d = 0;
    for (slong v = 0; v < nv; v++) d += a[v];
    return d;
}

/* ------------------------------------------------------------------ */
/*  The engine state                                                   */
/* ------------------------------------------------------------------ */

typedef struct {
    slong i, j;         /* basis indices, i < j */
    ulong* lcm;         /* nv exponents (owned) */
    ulong sugar;
    int dead;
} GfPair;

typedef struct {
    int mod;                    /* 0: Q (fmpz_mpoly), 1: GF(p) (nmod_mpoly) */
    slong nv;
    int ord;
    fmpz_mpoly_ctx_t Z;
    nmod_mpoly_ctx_t M;
    /* basis: every element ever added (pairs refer to them by index) */
    fmpz_mpoly_struct* zb;
    nmod_mpoly_struct* mb;
    ulong* lm;                  /* nb * nv */
    ulong* sugar;
    char* active;               /* still in the reducer set */
    slong nb, cap;
    /* a work budget (0: none): gf_run gives up -- `over` -- past this many
     * pairs or once a new element's coefficients exceed this many bits */
    slong budget_pairs; flint_bitcnt_t budget_bits; int over;
    int expired;                /* the deadline passed inside a reduction */
    GfPair* P;
    slong np, capp;
    /* GF(p) heap-division workspace: a pointer per reducer, the quotients */
    nmod_mpoly_struct** mq; nmod_mpoly_struct** mr;
    slong nq;                   /* quotient slots initialised */
    ulong* tmp;                 /* 2 nv scratch exponents */
} Gf;

static void gf_init(Gf* G, int mod, ulong p, slong nv, int ord) {
    memset(G, 0, sizeof *G);
    G->mod = mod; G->nv = nv; G->ord = ord;
    if (mod) nmod_mpoly_ctx_init(G->M, nv, (ordering_t)ord, p);
    else     fmpz_mpoly_ctx_init(G->Z, nv, (ordering_t)ord);
    G->tmp = flint_malloc(2 * (size_t)nv * sizeof(ulong));
}

static void gf_clear(Gf* G) {
    for (slong k = 0; k < G->nb; k++) {
        if (G->mod) nmod_mpoly_clear(G->mb + k, G->M);
        else        fmpz_mpoly_clear(G->zb + k, G->Z);
    }
    for (slong k = 0; k < G->nq; k++) { nmod_mpoly_clear(G->mq[k], G->M); flint_free(G->mq[k]); }
    for (slong k = 0; k < G->np; k++) flint_free(G->P[k].lcm);
    flint_free(G->zb); flint_free(G->mb); flint_free(G->lm); flint_free(G->sugar);
    flint_free(G->active); flint_free(G->P);
    flint_free(G->mq); flint_free(G->mr);
    flint_free(G->tmp);
    if (G->mod) nmod_mpoly_ctx_clear(G->M);
    else        fmpz_mpoly_ctx_clear(G->Z);
}

static int gf_is_zero(const Gf* G, const void* p) {
    return G->mod ? nmod_mpoly_is_zero((const nmod_mpoly_struct*)p, G->M)
                  : fmpz_mpoly_is_zero((const fmpz_mpoly_struct*)p, G->Z);
}

/* Room for one more basis element (and its quotient slot). */
static void gf_grow(Gf* G) {
    if (G->nb < G->cap) return;
    slong nc = G->cap ? 2 * G->cap : 16;
    if (G->mod) G->mb = flint_realloc(G->mb, (size_t)nc * sizeof(nmod_mpoly_struct));
    else        G->zb = flint_realloc(G->zb, (size_t)nc * sizeof(fmpz_mpoly_struct));
    G->lm = flint_realloc(G->lm, (size_t)(nc * G->nv) * sizeof(ulong));
    G->sugar = flint_realloc(G->sugar, (size_t)nc * sizeof(ulong));
    G->active = flint_realloc(G->active, (size_t)nc);
    if (G->mod) {
        G->mq = flint_realloc(G->mq, (size_t)nc * sizeof(nmod_mpoly_struct*));
        G->mr = flint_realloc(G->mr, (size_t)nc * sizeof(nmod_mpoly_struct*));
    }
    G->cap = nc;
}

/* Make p primitive with a positive leading coefficient (Q), or monic. */
static void gf_normalise(Gf* G, void* p) {
    if (G->mod) {
        nmod_mpoly_struct* q = p;
        if (!nmod_mpoly_is_zero(q, G->M)) nmod_mpoly_make_monic(q, q, G->M);
        return;
    }
    fmpz_mpoly_struct* q = p;
    if (q->length == 0) return;
    fmpz_t c; fmpz_init(c);
    _fmpz_vec_content(c, q->coeffs, q->length);
    if (fmpz_sgn(q->coeffs + 0) < 0) fmpz_neg(c, c);
    if (!fmpz_is_one(c)) fmpz_mpoly_scalar_divexact_fmpz(q, q, c, G->Z);
    fmpz_clear(c);
}

/* The GF(p) reducers -- every basis element (the active ones when
 * only_active) except `skip` -- as the pointer array the heap division
 * takes, shortest first (it reduces by the first divisor that applies; see
 * gf_divisor), with a quotient slot per reducer.  Returns their number. */
static int gf_len_cmp(const void* a, const void* b) {
    slong la = (*(nmod_mpoly_struct* const*)a)->length, lb = (*(nmod_mpoly_struct* const*)b)->length;
    return la < lb ? -1 : la > lb;
}

static slong gf_reducers(Gf* G, slong skip, int only_active) {
    slong n = 0;
    for (slong k = 0; k < G->nb; k++) {
        if ((only_active && !G->active[k]) || k == skip) continue;
        if (n == G->nq) {
            G->mq[n] = flint_malloc(sizeof(nmod_mpoly_struct));
            nmod_mpoly_init(G->mq[n], G->M);
            G->nq++;
        }
        G->mr[n] = G->mb + k;
        n++;
    }
    qsort(G->mr, (size_t)n, sizeof(nmod_mpoly_struct*), gf_len_cmp);
    return n;
}

/* The reducer whose leading monomial divides ex: the SHORTEST one.  Every
 * basis element may reduce -- an element leaves the reducer set of the pair
 * criteria (`active`) when a later one's leading monomial divides its own,
 * but it is still in the ideal, and it is usually the shorter of the two
 * (the inputs, the early elements).  Reducing by the minimal set alone made
 * the remainders snowball (50 000 terms where ~270 suffice on a Lex basis of
 * two quintics).  The final tail reduction (`only_active`) uses the minimal
 * set, which the reduced basis is made of.  -1 when none divides. */
static slong gf_divisor(const Gf* G, const ulong* ex, slong skip, int only_active) {
    slong best = -1, nv = G->nv;
    for (slong k = 0; k < G->nb; k++) {
        if ((only_active && !G->active[k]) || k == skip || !gf_divides(G->lm + k * nv, ex, nv)) continue;
        if (best < 0 || G->zb[k].length < G->zb[best].length) best = k;
    }
    return best;
}

/* The fraction-free normal form over Z.  The irreducible terms collect in
 * `done` (always above everything left in `rest`, so pushing keeps it
 * sorted); one step is rest <- a t - b m g in a single fused pass, a = lc(g)/c
 * and b = lc(t)/c for c their gcd, and multiplies `done` by a.  After every
 * step that scaled (a != 1) the joint content of done and rest is divided out,
 * as a rational engine's fractions are reduced after every operation -- so the
 * coefficients stay the size of the answer.  With `scale`, the factor the
 * result carries is tracked: the true normal form is h / scale.  (FLINT's quasidivrem_ideal
 * scales the whole dividend by every leading coefficient it meets and takes
 * the content only at the end: its coefficients exploded on a Lex basis the
 * rational engine finishes in 0.35 s.) */
static void gf_reduce_z(Gf* G, fmpz_mpoly_struct* h, slong skip, int only_active, fmpq* scale) {
    const fmpz_mpoly_ctx_struct* Z = G->Z;
    slong nv = G->nv;
    ulong* ex = G->tmp; ulong* me = G->tmp + nv;
    fmpz_mpoly_t done, rest, t, mg; fmpz_mpoly_init(done, Z); fmpz_mpoly_init(rest, Z);
    fmpz_mpoly_init(t, Z); fmpz_mpoly_init(mg, Z);
    fmpz_t a, b, c, one; fmpz_init(a); fmpz_init(b); fmpz_init(c); fmpz_init_set_ui(one, 1);
    fmpz_mpoly_swap(rest, h, Z);
    slong steps = 0;
    while (rest->length) {
        /* one normal form over Q can run for seconds; TimeConstrained must
         * not wait for it (the caller sees `expired` and unwinds) */
        if ((++steps & 63) == 0 && tc_deadline_passed()) { G->expired = 1; break; }
        slong k, d = -1;
        for (k = 0; k < rest->length; k++) {
            fmpz_mpoly_get_term_exp_ui(ex, rest, k, Z);
            if ((d = gf_divisor(G, ex, skip, only_active)) >= 0) break;
        }
        for (slong i = 0; i < k; i++) {                    /* irreducible: done */
            fmpz_mpoly_get_term_exp_ui(me, rest, i, Z);
            fmpz_mpoly_push_term_fmpz_ui(done, rest->coeffs + i, me, Z);
        }
        if (d < 0) { fmpz_mpoly_zero(rest, Z); break; }
        fmpz_mpoly_zero(t, Z);
        for (slong i = k; i < rest->length; i++) {
            fmpz_mpoly_get_term_exp_ui(me, rest, i, Z);
            fmpz_mpoly_push_term_fmpz_ui(t, rest->coeffs + i, me, Z);
        }
        const fmpz_mpoly_struct* g = G->zb + d;
        for (slong v = 0; v < nv; v++) me[v] = ex[v] - G->lm[d * nv + v];
        fmpz_mpoly_zero(mg, Z);
        fmpz_mpoly_set_coeff_fmpz_ui(mg, one, me, Z);           /* the monomial m */
        fmpz_mpoly_mul(mg, mg, g, Z);
        fmpz_gcd(c, g->coeffs + 0, t->coeffs + 0);
        fmpz_divexact(a, g->coeffs + 0, c);
        fmpz_divexact(b, t->coeffs + 0, c);
        fmpz_neg(b, b);
        fmpz_mpoly_scalar_fmma(rest, t, a, mg, b, Z);            /* a t - b m g */
        if (fmpz_is_one(a)) continue;                /* nothing was scaled */
        if (scale) fmpq_mul_fmpz(scale, scale, a);
        fmpz_mpoly_scalar_mul_fmpz(done, done, a, Z);
        /* the joint content, with an early exit at 1 */
        fmpz_zero(c);
        for (slong i = 0; i < rest->length && !fmpz_is_one(c); i++) fmpz_gcd(c, c, rest->coeffs + i);
        for (slong i = 0; i < done->length && !fmpz_is_one(c); i++) fmpz_gcd(c, c, done->coeffs + i);
        if (!fmpz_is_one(c) && !fmpz_is_zero(c)) {
            if (scale) fmpq_div_fmpz(scale, scale, c);
            fmpz_mpoly_scalar_divexact_fmpz(rest, rest, c, Z);
            fmpz_mpoly_scalar_divexact_fmpz(done, done, c, Z);
        }
    }
    fmpz_mpoly_swap(h, done, Z);
    fmpz_clear(a); fmpz_clear(b); fmpz_clear(c); fmpz_clear(one);
    fmpz_mpoly_clear(done, Z); fmpz_mpoly_clear(rest, Z); fmpz_mpoly_clear(t, Z); fmpz_mpoly_clear(mg, Z);
}

/* h = its normal form by the reducers (all but `skip`), normalised.  Over
 * GF(p) one heap division (no coefficient growth there); over Z the
 * fraction-free loop above. */
static void gf_reduce(Gf* G, void* h, slong skip, int only_active) {
    if (!G->mod) {
        gf_reduce_z(G, (fmpz_mpoly_struct*)h, skip, only_active, NULL);
        gf_normalise(G, h);
        return;
    }
    slong n = gf_reducers(G, skip, only_active);
    if (n > 0) {
        nmod_mpoly_t R; nmod_mpoly_init(R, G->M);
        nmod_mpoly_divrem_ideal(G->mq, R, (nmod_mpoly_struct*)h, G->mr, n, G->M);
        nmod_mpoly_swap((nmod_mpoly_struct*)h, R, G->M);
        nmod_mpoly_clear(R, G->M);
    }
    gf_normalise(G, h);
}

/* The S-polynomial of basis elements i, j into s (initialised). */
static void gf_spoly(Gf* G, void* s, slong i, slong j, const ulong* lcm) {
    slong nv = G->nv;
    ulong* mi = G->tmp; ulong* mj = G->tmp + nv;
    for (slong v = 0; v < nv; v++) {
        mi[v] = lcm[v] - G->lm[i * nv + v];
        mj[v] = lcm[v] - G->lm[j * nv + v];
    }
    if (G->mod) {
        /* monic elements: s = m_i f_i - m_j f_j */
        nmod_mpoly_t a, b, t; nmod_mpoly_init(a, G->M); nmod_mpoly_init(b, G->M); nmod_mpoly_init(t, G->M);
        nmod_mpoly_set_coeff_ui_ui(a, 1, mi, G->M);
        nmod_mpoly_set_coeff_ui_ui(b, 1, mj, G->M);
        nmod_mpoly_mul(t, G->mb + i, a, G->M);
        nmod_mpoly_mul(a, G->mb + j, b, G->M);
        nmod_mpoly_sub((nmod_mpoly_struct*)s, t, a, G->M);
        nmod_mpoly_clear(a, G->M); nmod_mpoly_clear(b, G->M); nmod_mpoly_clear(t, G->M);
        return;
    }
    /* primitive elements: s = (lc_j/c) m_i f_i - (lc_i/c) m_j f_j, c = gcd */
    const fmpz* li = G->zb[i].coeffs + 0;
    const fmpz* lj = G->zb[j].coeffs + 0;
    fmpz_t c, ci, cj; fmpz_init(c); fmpz_init(ci); fmpz_init(cj);
    fmpz_gcd(c, li, lj);
    fmpz_divexact(ci, lj, c);
    fmpz_divexact(cj, li, c);
    fmpz_mpoly_t a, b, t; fmpz_mpoly_init(a, G->Z); fmpz_mpoly_init(b, G->Z); fmpz_mpoly_init(t, G->Z);
    fmpz_mpoly_set_coeff_fmpz_ui(a, ci, mi, G->Z);
    fmpz_mpoly_set_coeff_fmpz_ui(b, cj, mj, G->Z);
    fmpz_mpoly_mul(t, G->zb + i, a, G->Z);
    fmpz_mpoly_mul(a, G->zb + j, b, G->Z);
    fmpz_mpoly_sub((fmpz_mpoly_struct*)s, t, a, G->Z);
    fmpz_mpoly_clear(a, G->Z); fmpz_mpoly_clear(b, G->Z); fmpz_mpoly_clear(t, G->Z);
    fmpz_clear(c); fmpz_clear(ci); fmpz_clear(cj);
}

/* ------------------------------------------------------------------ */
/*  Pairs: Gebauer-Moeller                                             */
/* ------------------------------------------------------------------ */

static void gf_push_pair(Gf* G, slong i, slong j, ulong* lcm, ulong sugar) {
    if (G->np == G->capp) {
        G->capp = G->capp ? 2 * G->capp : 64;
        G->P = flint_realloc(G->P, (size_t)G->capp * sizeof(GfPair));
    }
    GfPair* q = G->P + G->np++;
    q->i = i; q->j = j; q->lcm = lcm; q->sugar = sugar; q->dead = 0;
}

/* The sugar of the pair (i, j) with lcm l. */
static ulong gf_pair_sugar(const Gf* G, slong i, slong j, const ulong* l) {
    ulong dl = gf_deg(l, G->nv);
    ulong si = G->sugar[i] + dl - gf_deg(G->lm + i * G->nv, G->nv);
    ulong sj = G->sugar[j] + dl - gf_deg(G->lm + j * G->nv, G->nv);
    return si > sj ? si : sj;
}

/* Register basis element h (already stored, its lm set): the Gebauer-Moeller
 * UPDATE of the pair set and of the reducer set. */
static void gf_update(Gf* G, slong h) {
    slong nv = G->nv;
    const ulong* t = G->lm + h * nv;
    /* the candidates (g, h), g active */
    slong nc = 0;
    slong* cg = flint_malloc((size_t)(G->nb + 1) * sizeof(slong));
    ulong* cl = flint_malloc((size_t)((G->nb + 1) * nv) * sizeof(ulong));
    char* keep = flint_malloc((size_t)(G->nb + 1));
    for (slong g = 0; g < h; g++) {
        if (!G->active[g]) continue;
        cg[nc] = g;
        gf_lcm(cl + nc * nv, G->lm + g * nv, t, nv);
        keep[nc] = 1;
        nc++;
    }
    /* chain criterion on the new pairs: drop (g, h) when another candidate's
     * lcm properly divides its lcm; of equal lcms keep one, and none of them
     * when one is coprime (the product criterion would discard it) */
    for (slong a = 0; a < nc; a++) {
        for (slong b = 0; b < nc && keep[a]; b++) {
            if (a == b || !keep[b]) continue;
            const ulong* la = cl + a * nv; const ulong* lb = cl + b * nv;
            if (!gf_divides(lb, la, nv)) continue;
            int equal = memcmp(la, lb, (size_t)nv * sizeof(ulong)) == 0;
            if (!equal || b < a) keep[a] = 0;
        }
    }
    for (slong a = 0; a < nc; a++) {
        if (!keep[a]) continue;
        /* equal-lcm group with a coprime member: all useless */
        int coprime = gf_coprime(G->lm + cg[a] * nv, t, nv);
        for (slong b = 0; b < nc && !coprime; b++) {
            if (b == a) continue;
            if (memcmp(cl + a * nv, cl + b * nv, (size_t)nv * sizeof(ulong)) == 0 &&
                gf_coprime(G->lm + cg[b] * nv, t, nv))
                coprime = 1;
        }
        if (coprime) keep[a] = 0;
    }
    /* criterion B on the old pairs */
    ulong* l1 = G->tmp; ulong* l2 = G->tmp + nv;
    for (slong k = 0; k < G->np; k++) {
        GfPair* q = G->P + k;
        if (q->dead || !gf_divides(t, q->lcm, nv)) continue;
        gf_lcm(l1, G->lm + q->i * nv, t, nv);
        gf_lcm(l2, G->lm + q->j * nv, t, nv);
        if (memcmp(l1, q->lcm, (size_t)nv * sizeof(ulong)) != 0 &&
            memcmp(l2, q->lcm, (size_t)nv * sizeof(ulong)) != 0)
            q->dead = 1;
    }
    /* compact, then add the new pairs */
    slong w = 0;
    for (slong k = 0; k < G->np; k++) {
        if (G->P[k].dead) { flint_free(G->P[k].lcm); continue; }
        G->P[w++] = G->P[k];
    }
    G->np = w;
    for (slong a = 0; a < nc; a++) {
        if (!keep[a]) continue;
        ulong* l = flint_malloc((size_t)nv * sizeof(ulong));
        memcpy(l, cl + a * nv, (size_t)nv * sizeof(ulong));
        gf_push_pair(G, cg[a], h, l, gf_pair_sugar(G, cg[a], h, l));
    }
    flint_free(cg); flint_free(cl); flint_free(keep);
    /* the reducer set: h replaces every element whose lm it divides */
    for (slong g = 0; g < h; g++)
        if (G->active[g] && gf_divides(t, G->lm + g * nv, nv)) G->active[g] = 0;
    G->active[h] = 1;
}

/* The next pair: smallest sugar, then smallest lcm. */
static slong gf_pick(const Gf* G) {
    slong best = -1;
    for (slong k = 0; k < G->np; k++) {
        const GfPair* q = G->P + k;
        if (best < 0) { best = k; continue; }
        const GfPair* b = G->P + best;
        if (q->sugar != b->sugar) { if (q->sugar < b->sugar) best = k; continue; }
        if (gf_cmp(q->lcm, b->lcm, G->nv, G->ord) > 0) best = k;
    }
    return best;
}

/* Store the (nonzero, normalised) polynomial *h as a new basis element with
 * the given sugar, moving it in, and update the pairs. */
static void gf_add(Gf* G, void* h, ulong sugar) {
    gf_grow(G);
    slong k = G->nb;
    if (G->mod) {
        nmod_mpoly_init(G->mb + k, G->M);
        nmod_mpoly_swap(G->mb + k, (nmod_mpoly_struct*)h, G->M);
        nmod_mpoly_get_term_exp_ui(G->lm + k * G->nv, G->mb + k, 0, G->M);
    } else {
        fmpz_mpoly_init(G->zb + k, G->Z);
        fmpz_mpoly_swap(G->zb + k, (fmpz_mpoly_struct*)h, G->Z);
        fmpz_mpoly_get_term_exp_ui(G->lm + k * G->nv, G->zb + k, 0, G->Z);
    }
    G->sugar[k] = sugar;
    G->active[k] = 0;
    G->nb++;
    gf_update(G, k);
}

/* Run to completion.  0 when the deadline passed (the caller unwinds). */
static int gf_run(Gf* G) {
    nmod_mpoly_t ms; fmpz_mpoly_t zs;
    if (G->mod) nmod_mpoly_init(ms, G->M); else fmpz_mpoly_init(zs, G->Z);
    void* s = G->mod ? (void*)ms : (void*)zs;
    int ok = 1;
    slong done_pairs = 0;
    while (G->np > 0) {
        if (tc_deadline_passed()) { ok = 0; break; }
        if (G->budget_pairs && ++done_pairs > G->budget_pairs) { G->over = 1; ok = 0; break; }
        slong k = gf_pick(G);
        GfPair q = G->P[k];
        G->P[k] = G->P[--G->np];
        gf_spoly(G, s, q.i, q.j, q.lcm);
        flint_free(q.lcm);
        gf_reduce(G, s, -1, 0);
        if (G->expired) { ok = 0; break; }
        if (!gf_is_zero(G, s)) {
            /* a reduction step m g can raise the degree past the pair's sugar
             * (not in a degree order, but in Lex); the sugar is at least the
             * remainder's own degree */
            ulong d = G->mod ? (ulong)nmod_mpoly_total_degree_si((nmod_mpoly_struct*)s, G->M)
                             : (ulong)fmpz_mpoly_total_degree_si((fmpz_mpoly_struct*)s, G->Z);
            if (G->budget_bits && !G->mod) {
                const fmpz_mpoly_struct* z = s;
                slong b = _fmpz_vec_max_bits(z->coeffs, z->length);
                if ((flint_bitcnt_t)(b < 0 ? -b : b) > G->budget_bits) { G->over = 1; ok = 0; break; }
            }
            gf_add(G, s, d > q.sugar ? d : q.sugar);
        }
    }
    if (G->mod) nmod_mpoly_clear(ms, G->M); else fmpz_mpoly_clear(zs, G->Z);
    return ok;
}

/* The reduced basis: tail-reduce each reducer by the others (one pass: the
 * set is minimal, so reducibility depends only on the fixed leading
 * monomials), then the indices sorted ascending by leading monomial. */
typedef struct { const ulong* lm; slong nv; int ord; slong k; } GfSort;

static int gf_sort_cmp(const void* pa, const void* pb) {
    const GfSort* a = pa; const GfSort* b = pb;
    return -gf_cmp(a->lm, b->lm, a->nv, a->ord);     /* ascending */
}

static slong gf_final(Gf* G, slong* idx) {
    slong n = 0;
    for (slong k = 0; k < G->nb; k++) {
        if (!G->active[k]) continue;
        gf_reduce(G, G->mod ? (void*)(G->mb + k) : (void*)(G->zb + k), k, 1);
        idx[n++] = k;
    }
    GfSort* S = flint_malloc((size_t)(n ? n : 1) * sizeof(GfSort));
    for (slong i = 0; i < n; i++) {
        S[i].lm = G->lm + idx[i] * G->nv; S[i].nv = G->nv; S[i].ord = G->ord; S[i].k = idx[i];
    }
    qsort(S, (size_t)n, sizeof(GfSort), gf_sort_cmp);
    for (slong i = 0; i < n; i++) idx[i] = S[i].k;
    flint_free(S);
    return n;
}

/* ------------------------------------------------------------------ */
/*  Conversions                                                        */
/* ------------------------------------------------------------------ */

/* A rational GBPoly into Q (denominators cleared, primitive) or GF(p).  0
 * when a denominator is divisible by p. */
static int gf_from_gbpoly(Gf* G, void* out, const GBPoly* g, ulong p) {
    slong nv = G->nv;
    ulong* ex = G->tmp;
    if (G->mod) {
        nmod_mpoly_struct* o = out;
        nmod_mpoly_zero(o, G->M);
        nmod_t mod; nmod_init(&mod, p);
        for (size_t t = 0; t < g->n_terms; t++) {
            ulong d = mpz_fdiv_ui(mpq_denref(g->coefs[t]), p);
            if (d == 0) return 0;
            ulong c = nmod_mul(mpz_fdiv_ui(mpq_numref(g->coefs[t]), p), n_invmod(d, p), mod);
            if (!c) continue;
            for (slong v = 0; v < nv; v++) ex[v] = (ulong)g->exps[t * (size_t)nv + v];
            nmod_mpoly_push_term_ui_ui(o, c, ex, G->M);
        }
        nmod_mpoly_sort_terms(o, G->M);
        nmod_mpoly_combine_like_terms(o, G->M);
        return 1;
    }
    fmpz_mpoly_struct* o = out;
    fmpz_mpoly_zero(o, G->Z);
    mpz_t L, k; mpz_init_set_ui(L, 1); mpz_init(k);
    for (size_t t = 0; t < g->n_terms; t++) mpz_lcm(L, L, mpq_denref(g->coefs[t]));
    fmpz_t c; fmpz_init(c);
    for (size_t t = 0; t < g->n_terms; t++) {
        mpz_divexact(k, L, mpq_denref(g->coefs[t]));
        mpz_mul(k, k, mpq_numref(g->coefs[t]));
        fmpz_set_mpz(c, k);
        for (slong v = 0; v < nv; v++) ex[v] = (ulong)g->exps[t * (size_t)nv + v];
        fmpz_mpoly_push_term_fmpz_ui(o, c, ex, G->Z);
    }
    fmpz_clear(c); mpz_clear(L); mpz_clear(k);
    fmpz_mpoly_sort_terms(o, G->Z);
    fmpz_mpoly_combine_like_terms(o, G->Z);
    return 1;
}

static GBPoly* gf_to_gbpoly(const Gf* G, slong k, const GBPoly* like) {
    slong nv = G->nv;
    GBPoly* r = gb_poly_new(like->n_vars, like->order, like->elim_pivot);
    if (like->wmat) gb_poly_set_wmat(r, like->wmat);
    int* ie = malloc((size_t)(nv ? nv : 1) * sizeof(int));
    ulong* ex = G->tmp;
    mpq_t c; mpq_init(c);
    slong len = G->mod ? G->mb[k].length : G->zb[k].length;
    gb_poly_reserve(r, (size_t)len);
    for (slong t = 0; t < len; t++) {
        if (G->mod) {
            nmod_mpoly_get_term_exp_ui(ex, G->mb + k, t, G->M);
            mpq_set_ui(c, G->mb[k].coeffs[t], 1);
        } else {
            fmpz_mpoly_get_term_exp_ui(ex, G->zb + k, t, G->Z);
            fmpz_get_mpz(mpq_numref(c), G->zb[k].coeffs + t);
            mpz_set_ui(mpq_denref(c), 1);
        }
        for (slong v = 0; v < nv; v++) ie[v] = (int)ex[v];
        gb_poly_push_term(r, ie, c);
    }
    mpq_clear(c);
    free(ie);
    gb_poly_normalize(r);
    return r;
}

/* A GBPoly of the unit ideal's basis {1}, shaped like `like`. */
static GBPoly* gf_one_gbpoly(const GBPoly* like) {
    GBPoly* r = gb_poly_new(like->n_vars, like->order, like->elim_pivot);
    if (like->wmat) gb_poly_set_wmat(r, like->wmat);
    int* z = calloc((size_t)like->n_vars, sizeof(int));
    gb_poly_push_term_si(r, z, 1, 1);
    free(z);
    gb_poly_normalize(r);
    return r;
}

/* Load F into G and run it to the reduced basis: idx[0..*nr) the basis
 * elements in ascending order (idx sized G->nb by the caller after the run:
 * it is allocated here).  0 on a denominator divisible by p or an expired
 * deadline (*timed_out). */
static int gf_engine(Gf* G, GBPoly* const* F, size_t n, ulong p, slong** idx, slong* nr,
                     int* timed_out) {
    nmod_mpoly_t mh; fmpz_mpoly_t zh;
    if (G->mod) nmod_mpoly_init(mh, G->M); else fmpz_mpoly_init(zh, G->Z);
    void* h = G->mod ? (void*)mh : (void*)zh;
    int ok = 1;
    for (size_t i = 0; i < n && ok; i++) {
        if (!F[i]->n_terms) continue;
        if (!(ok = gf_from_gbpoly(G, h, F[i], p))) break;
        ulong sugar = G->mod ? (ulong)nmod_mpoly_total_degree_si(mh, G->M)
                             : (ulong)fmpz_mpoly_total_degree_si(zh, G->Z);
        gf_reduce(G, h, -1, 0);
        if (G->expired) { ok = 0; break; }
        if (gf_is_zero(G, h)) continue;
        gf_add(G, h, sugar);
    }
    if (G->mod) nmod_mpoly_clear(mh, G->M); else fmpz_mpoly_clear(zh, G->Z);
    if (ok) ok = gf_run(G);
    if (ok) {
        *idx = flint_malloc((size_t)(G->nb ? G->nb : 1) * sizeof(slong));
        *nr = gf_final(G, *idx);
        if (G->expired) { ok = 0; flint_free(*idx); *idx = NULL; *nr = 0; }
    }
    *timed_out = !ok && tc_deadline_passed();
    return ok;
}

static int gf_is_unit(const Gf* G, slong k) {
    return G->mod ? nmod_mpoly_is_ui(G->mb + k, G->M) : fmpz_mpoly_is_fmpz(G->zb + k, G->Z);
}

/* ------------------------------------------------------------------ */
/*  Lex over Q, zero-dimensional: DRL, then FGLM                       */
/* ------------------------------------------------------------------ */

/* Lex over Q is where Buchberger's coefficients explode: cyclic-5 does not
 * finish in minutes on fraction-free arithmetic, while its DRL basis takes
 * 4 ms.  For a zero-dimensional ideal the Lex basis is instead read off the
 * quotient Q[X]/I (dimension D, the standard monomials of the DRL basis):
 *
 *  1. the multiplication matrices M_v (x_v times each standard monomial, as a
 *     normal form), exactly over Q;
 *  2. FGLM modulo 62-bit primes: the monomials in increasing Lex order, each
 *     NF vector M_v NF(m') reduced against the accepted ones; a dependency is
 *     a basis element (its leading monomial, coefficients on the accepted
 *     monomials), an independent vector a new standard monomial.  A prime
 *     whose staircase differs from the first is unlucky and dropped;
 *  3. the coefficients combined by CRT and rationally reconstructed, until two
 *     consecutive reconstructions agree;
 *  4. CERTIFIED exactly: every element reduces to 0 by the DRL basis (so it is
 *     in I), and its leading monomials leave exactly D standard monomials.  Then
 *     LT(I) contains <LM(G)> and both have codimension D, so they are equal:
 *     G is the reduced Lex basis of I.  No probabilistic step.
 *
 * Declines (the Lex engine runs directly) for a positive-dimensional ideal,
 * D > GF_FGLM_MAXD, or no agreement within GF_FGLM_PRIMES primes. */
#define GF_FGLM_MAXD   600
#define GF_LEX_BUDGET_PAIRS 64        /* the direct Lex attempt first (gf_basis) */
#define GF_LEX_BUDGET_BITS  256
#define GF_FGLM_PRIMES 400

/* A set of exponent vectors with their indices (open addressing). */
typedef struct { ulong* ex; slong nv, n, cap; slong* tab; slong tcap; } GfMon;

static ulong gf_mon_hash(const ulong* e, slong nv) {
    ulong h = UWORD(1469598103934665603);
    for (slong v = 0; v < nv; v++) h = (h ^ e[v]) * UWORD(1099511628211);
    return h;
}
static void gfm_init(GfMon* S, slong nv) {
    memset(S, 0, sizeof *S); S->nv = nv;
    S->tcap = 64; S->tab = flint_calloc((size_t)S->tcap, sizeof(slong));
}
static void gfm_clear(GfMon* S) { flint_free(S->ex); flint_free(S->tab); }
static slong gfm_find(const GfMon* S, const ulong* e) {
    ulong h = gf_mon_hash(e, S->nv) & (ulong)(S->tcap - 1);
    while (S->tab[h]) {
        slong k = S->tab[h] - 1;
        if (memcmp(S->ex + k * S->nv, e, (size_t)S->nv * sizeof(ulong)) == 0) return k;
        h = (h + 1) & (ulong)(S->tcap - 1);
    }
    return -1;
}
static slong gfm_add(GfMon* S, const ulong* e) {          /* index; new or old */
    slong k = gfm_find(S, e);
    if (k >= 0) return k;
    if (S->n == S->cap) {
        S->cap = S->cap ? 2 * S->cap : 64;
        S->ex = flint_realloc(S->ex, (size_t)(S->cap * S->nv) * sizeof(ulong));
    }
    memcpy(S->ex + S->n * S->nv, e, (size_t)S->nv * sizeof(ulong));
    k = S->n++;
    if (2 * S->n > S->tcap) {                             /* rehash */
        flint_free(S->tab);
        S->tcap *= 2;
        S->tab = flint_calloc((size_t)S->tcap, sizeof(slong));
        for (slong i = 0; i < S->n; i++) {
            ulong h = gf_mon_hash(S->ex + i * S->nv, S->nv) & (ulong)(S->tcap - 1);
            while (S->tab[h]) h = (h + 1) & (ulong)(S->tcap - 1);
            S->tab[h] = i + 1;
        }
    } else {
        ulong h = gf_mon_hash(e, S->nv) & (ulong)(S->tcap - 1);
        while (S->tab[h]) h = (h + 1) & (ulong)(S->tcap - 1);
        S->tab[h] = k + 1;
    }
    return k;
}

static int gf_divisible_by_any(const Gf* G, const slong* idx, slong ng, const ulong* e) {
    for (slong i = 0; i < ng; i++) if (gf_divides(G->lm + idx[i] * G->nv, e, G->nv)) return 1;
    return 0;
}

/* The standard monomials of the reduced basis idx[0..ng) (index 0 is 1), or 0
 * when the ideal is not zero-dimensional or D > maxd. */
static int gf_staircase(const Gf* G, const slong* idx, slong ng, GfMon* S, slong maxd) {
    slong nv = G->nv;
    for (slong v = 0; v < nv; v++) {                      /* a pure power per variable */
        int pure = 0;
        for (slong i = 0; i < ng && !pure; i++) {
            const ulong* l = G->lm + idx[i] * nv;
            pure = l[v] > 0;
            for (slong w = 0; w < nv && pure; w++) if (w != v && l[w]) pure = 0;
        }
        if (!pure) return 0;
    }
    ulong* e = flint_calloc((size_t)nv, sizeof(ulong));
    gfm_add(S, e);
    for (slong i = 0; i < S->n; i++) {
        for (slong v = 0; v < nv; v++) {
            memcpy(e, S->ex + i * nv, (size_t)nv * sizeof(ulong));
            e[v]++;
            if (gf_divisible_by_any(G, idx, ng, e) || gfm_find(S, e) >= 0) continue;
            gfm_add(S, e);
            if (S->n > maxd) { flint_free(e); return 0; }
        }
    }
    flint_free(e);
    return 1;
}

/* The multiplication matrices over Q: column c of matrix v is NF(x_v b_c) on
 * the standard monomials, as numerators num[(v D + c) D + r] over the column
 * denominator den[v D + c]. */
static int gf_mult_matrices(Gf* G, const slong* idx, slong ng, const GfMon* S, fmpz* num, fmpz* den) {
    slong nv = G->nv, D = S->n;
    ulong* e = flint_malloc((size_t)nv * sizeof(ulong));
    ulong* te = flint_malloc((size_t)nv * sizeof(ulong));
    fmpz_mpoly_t m; fmpz_mpoly_init(m, G->Z);
    fmpq_t sc; fmpq_init(sc);
    fmpz_t one; fmpz_init_set_ui(one, 1);
    int ok = 1;
    for (slong v = 0; v < nv && ok; v++) {
        for (slong c = 0; c < D && ok; c++) {
            fmpz* col = num + (v * D + c) * D;
            fmpz* dc = den + v * D + c;
            memcpy(e, S->ex + c * nv, (size_t)nv * sizeof(ulong));
            e[v]++;
            slong k = gfm_find(S, e);
            if (k >= 0) { fmpz_one(col + k); fmpz_one(dc); continue; }
            /* a border monomial: the leading monomial of a basis element is
             * x_v b_c - g / lc(g); anything else is reduced */
            slong b = -1;
            for (slong i = 0; i < ng && b < 0; i++)
                if (memcmp(G->lm + idx[i] * nv, e, (size_t)nv * sizeof(ulong)) == 0) b = idx[i];
            if (b >= 0) {
                const fmpz_mpoly_struct* g = G->zb + b;
                fmpz_set(dc, g->coeffs + 0);
                for (slong t = 1; t < g->length && ok; t++) {
                    fmpz_mpoly_get_term_exp_ui(te, g, t, G->Z);
                    slong r = gfm_find(S, te);
                    if (r < 0) { ok = 0; break; }             /* not reduced: impossible */
                    fmpz_neg(col + r, g->coeffs + t);
                }
                continue;
            }
            fmpz_mpoly_zero(m, G->Z);
            fmpz_mpoly_set_coeff_fmpz_ui(m, one, e, G->Z);
            fmpq_one(sc);
            gf_reduce_z(G, m, -1, 1, sc);                   /* NF = m / sc */
            for (slong t = 0; t < m->length && ok; t++) {
                fmpz_mpoly_get_term_exp_ui(te, m, t, G->Z);
                slong r = gfm_find(S, te);
                if (r < 0) { ok = 0; break; }
                fmpz_mul(col + r, m->coeffs + t, fmpq_denref(sc));
            }
            fmpz_set(dc, fmpq_numref(sc));
        }
    }
    fmpz_clear(one); fmpq_clear(sc);
    fmpz_mpoly_clear(m, G->Z);
    flint_free(e); flint_free(te);
    return ok;
}

/* One FGLM run modulo p.  Out: the accepted (standard) monomials in Lex order
 * (stair, nstair), and the relations: leading monomial rel_ex[k], the
 * coefficients rel_c[k][0..rel_n[k]) on stair[0..rel_n[k]) (LM = sum c s). */
typedef struct {
    GfMon stair;                 /* accepted monomials, Lex ascending */
    ulong* rel_ex; slong nrel, crel;
    slong* rel_n;
    ulong** rel_c;
} GfFglm;

static void gff_clear(GfFglm* R) {
    gfm_clear(&R->stair);
    for (slong k = 0; k < R->nrel; k++) flint_free(R->rel_c[k]);
    flint_free(R->rel_c); flint_free(R->rel_n); flint_free(R->rel_ex);
}

static int gf_fglm_mod(const ulong* M, slong nv, slong D, slong one_idx, nmod_t mod, GfFglm* R) {
    gfm_init(&R->stair, nv);
    R->rel_ex = NULL; R->nrel = R->crel = 0; R->rel_n = NULL; R->rel_c = NULL;
    /* candidates: monomial, its predecessor (accepted index) and variable */
    GfMon cand; gfm_init(&cand, nv);
    slong ccap = 64, *cpre = flint_malloc((size_t)ccap * sizeof(slong)), *cvar = flint_malloc((size_t)ccap * sizeof(slong));
    char* cdone = flint_calloc((size_t)ccap, 1);
    /* accepted: NF vectors; echelon rows (pivot, row, combination) */
    ulong* vec = flint_malloc((size_t)(D * D) * sizeof(ulong));        /* D accepted max */
    ulong* row = flint_malloc((size_t)((D + 1) * D) * sizeof(ulong));   /* + the slot reduced into */
    ulong* comb = flint_calloc((size_t)(D * D), sizeof(ulong));
    slong* piv = flint_malloc((size_t)D * sizeof(slong));
    ulong* v = flint_malloc((size_t)D * sizeof(ulong));
    ulong* lam = flint_malloc((size_t)D * sizeof(ulong));
    ulong* e = flint_calloc((size_t)nv, sizeof(ulong));
    int ok = 1;
    gfm_add(&cand, e); cpre[0] = -1; cvar[0] = -1;
    for (;;) {
        /* the Lex-smallest open candidate */
        slong best = -1;
        for (slong k = 0; k < cand.n; k++) {
            if (cdone[k]) continue;
            if (best < 0 || gf_cmp(cand.ex + k * nv, cand.ex + best * nv, nv, ORD_LEX) > 0) best = k;
        }
        if (best < 0) break;
        cdone[best] = 1;
        const ulong* m = cand.ex + best * nv;
        int divisible = 0;
        for (slong k = 0; k < R->nrel && !divisible; k++) divisible = gf_divides(R->rel_ex + k * nv, m, nv);
        if (divisible) continue;
        /* v = NF(m) = M_var NF(pre), or the unit vector of 1 */
        slong na = R->stair.n;
        if (cpre[best] < 0) {
            memset(v, 0, (size_t)D * sizeof(ulong));
            v[one_idx] = 1;
        } else {
            const ulong* Mv = M + (size_t)cvar[best] * (size_t)(D * D);   /* column-major */
            const ulong* u = vec + cpre[best] * D;
            memset(v, 0, (size_t)D * sizeof(ulong));
            for (slong c = 0; c < D; c++)
                if (u[c]) _nmod_vec_scalar_addmul_nmod(v, Mv + c * D, D, u[c], mod);
        }
        memset(lam, 0, (size_t)D * sizeof(ulong));
        ulong* vr = row + na * D;                        /* reduce into the next row slot */
        memcpy(vr, v, (size_t)D * sizeof(ulong));
        for (slong t = 0; t < na; t++) {
            ulong a = vr[piv[t]];
            if (!a) continue;
            _nmod_vec_scalar_addmul_nmod(vr, row + t * D, D, nmod_neg(a, mod), mod);
            _nmod_vec_scalar_addmul_nmod(lam, comb + t * D, na, a, mod);
        }
        slong pv = -1;
        for (slong r = 0; r < D && pv < 0; r++) if (vr[r]) pv = r;
        if (pv < 0) {                                    /* a relation m = sum lam s */
            if (R->nrel == R->crel) {
                R->crel = R->crel ? 2 * R->crel : 16;
                R->rel_ex = flint_realloc(R->rel_ex, (size_t)(R->crel * nv) * sizeof(ulong));
                R->rel_n = flint_realloc(R->rel_n, (size_t)R->crel * sizeof(slong));
                R->rel_c = flint_realloc(R->rel_c, (size_t)R->crel * sizeof(ulong*));
            }
            memcpy(R->rel_ex + R->nrel * nv, m, (size_t)nv * sizeof(ulong));
            R->rel_n[R->nrel] = na;
            R->rel_c[R->nrel] = flint_malloc((size_t)(na ? na : 1) * sizeof(ulong));
            memcpy(R->rel_c[R->nrel], lam, (size_t)na * sizeof(ulong));
            R->nrel++;
            continue;
        }
        if (na == D) { ok = 0; break; }                   /* rank > D: impossible */
        /* accept m: row normalised at its pivot; comb = (e_m - lam) / pivot */
        ulong inv = n_invmod(vr[pv], mod.n);
        _nmod_vec_scalar_mul_nmod(vr, vr, D, inv, mod);
        ulong* cb = comb + na * D;
        for (slong k = 0; k < na; k++) cb[k] = nmod_mul(nmod_neg(lam[k], mod), inv, mod);
        cb[na] = inv;
        piv[na] = pv;
        memcpy(vec + na * D, v, (size_t)D * sizeof(ulong));
        gfm_add(&R->stair, m);
        slong me = R->stair.n - 1;
        for (slong w = 0; w < nv; w++) {
            memcpy(e, m, (size_t)nv * sizeof(ulong));
            e[w]++;
            if (gfm_find(&cand, e) >= 0) continue;
            slong k = gfm_add(&cand, e);
            if (k >= ccap) {
                slong nc = 2 * ccap;
                cpre = flint_realloc(cpre, (size_t)nc * sizeof(slong));
                cvar = flint_realloc(cvar, (size_t)nc * sizeof(slong));
                cdone = flint_realloc(cdone, (size_t)nc);
                memset(cdone + ccap, 0, (size_t)(nc - ccap));
                ccap = nc;
            }
            cpre[k] = me; cvar[k] = w;
        }
    }
    if (ok && R->stair.n != D) ok = 0;                    /* an unlucky prime */
    flint_free(e); flint_free(lam); flint_free(v); flint_free(piv);
    flint_free(comb); flint_free(row); flint_free(vec);
    flint_free(cdone); flint_free(cvar); flint_free(cpre);
    gfm_clear(&cand);
    if (!ok) gff_clear(R);
    return ok;
}

/* Two FGLM runs with the same staircase and the same leading monomials. */
static int gff_same_shape(const GfFglm* A, const GfFglm* B, slong nv) {
    if (A->stair.n != B->stair.n || A->nrel != B->nrel) return 0;
    if (memcmp(A->stair.ex, B->stair.ex, (size_t)(A->stair.n * nv) * sizeof(ulong))) return 0;
    if (memcmp(A->rel_ex, B->rel_ex, (size_t)(A->nrel * nv) * sizeof(ulong))) return 0;
    for (slong k = 0; k < A->nrel; k++) if (A->rel_n[k] != B->rel_n[k]) return 0;
    return 1;
}

/* The Lex basis of the zero-dimensional ideal whose reduced DRL basis is
 * Dg (idx[0..ng)), as primitive fmpz_mpoly in Dg's (DRL) context, ascending in
 * Lex.  1 on success with *out (nout) allocated; 0 declines. */
static int gf_fglm(Gf* Dg, const slong* idx, slong ng, fmpz_mpoly_struct** out, slong* nout) {
    slong nv = Dg->nv;
    GfMon S; gfm_init(&S, nv);
    if (!gf_staircase(Dg, idx, ng, &S, GF_FGLM_MAXD)) { gfm_clear(&S); return 0; }
    slong D = S.n;
    fmpz* num = _fmpz_vec_init(nv * D * D);
    fmpz* den = _fmpz_vec_init(nv * D);
    int ok = gf_mult_matrices(Dg, idx, ng, &S, num, &den[0]) && !Dg->expired;
    ulong* Mp = ok ? flint_malloc((size_t)(nv * D * D) * sizeof(ulong)) : NULL;
    GfFglm ref; int have = 0;
    memset(&ref, 0, sizeof ref);
    fmpz* crt = NULL; fmpz_t mod_prod; fmpz_init(mod_prod);
    fmpq* rec = NULL; int have_rec = 0, done = 0;
    slong ncoef = 0;
    ulong p = UWORD(1) << 62;
    for (slong it = 0; ok && !done && it < GF_FGLM_PRIMES; it++) {
        if (tc_deadline_passed()) { ok = 0; break; }
        p = n_nextprime(p, 1);
        nmod_t mod; nmod_init(&mod, p);
        int good = 1;
        for (slong k = 0; k < nv * D && good; k++) good = fmpz_fdiv_ui(den + k, p) != 0;
        if (!good) continue;
        for (slong k = 0; k < nv * D; k++) {
            ulong di = n_invmod(fmpz_fdiv_ui(den + k, p), p);
            fmpz* col = num + k * D;
            ulong* mc = Mp + k * D;
            for (slong r = 0; r < D; r++)
                mc[r] = fmpz_is_zero(col + r) ? 0 : nmod_mul(fmpz_fdiv_ui(col + r, p), di, mod);
        }
        GfFglm R;
        if (!gf_fglm_mod(Mp, nv, D, 0, mod, &R)) continue;
        if (!have) {
            ref = R; have = 1;
            for (slong k = 0; k < ref.nrel; k++) ncoef += ref.rel_n[k];
            crt = _fmpz_vec_init(ncoef);
            rec = flint_malloc((size_t)(ncoef ? ncoef : 1) * sizeof(fmpq));
            for (slong k = 0; k < ncoef; k++) fmpq_init(rec + k);
            slong q = 0;
            for (slong k = 0; k < ref.nrel; k++)
                for (slong i = 0; i < ref.rel_n[k]; i++) fmpz_set_ui(crt + q++, ref.rel_c[k][i]);
            fmpz_set_ui(mod_prod, p);
            continue;
        }
        if (!gff_same_shape(&ref, &R, nv)) { gff_clear(&R); continue; }   /* unlucky */
        slong q = 0;
        for (slong k = 0; k < R.nrel; k++)
            for (slong i = 0; i < R.rel_n[k]; i++, q++)
                fmpz_CRT_ui(crt + q, crt + q, mod_prod, R.rel_c[k][i], p, 0);
        fmpz_mul_ui(mod_prod, mod_prod, p);
        gff_clear(&R);
        /* reconstruct; two equal reconstructions in a row are a candidate */
        int all = 1, same = have_rec;
        fmpq_t r; fmpq_init(r);
        for (slong k = 0; k < ncoef && all; k++) {
            all = fmpq_reconstruct_fmpz(r, crt + k, mod_prod);
            if (all && same && !fmpq_equal(r, rec + k)) same = 0;
            if (all) fmpq_set(rec + k, r);
        }
        fmpq_clear(r);
        have_rec = all;
        if (!all || !same) continue;
        /* the candidate, certified: each element in I, staircase size D */
        fmpz_mpoly_struct* B = flint_malloc((size_t)(ref.nrel ? ref.nrel : 1) * sizeof(fmpz_mpoly_struct));
        fmpq_mpoly_t Q; fmpq_mpoly_ctx_t QC; fmpq_mpoly_ctx_init(QC, nv, ORD_LEX);
        fmpq_mpoly_init(Q, QC);
        int cert = ref.stair.n == D;
        slong q2 = 0, nb = 0;
        for (slong k = 0; k < ref.nrel; k++) {
            /* LM - sum c s over Q, primitive over Z, in the DRL context */
            fmpq_mpoly_zero(Q, QC);
            {
                fmpq_t one; fmpq_init(one); fmpq_one(one);
                fmpq_mpoly_set_coeff_fmpq_ui(Q, one, ref.rel_ex + k * nv, QC);
                fmpq_clear(one);
            }
            for (slong i = 0; i < ref.rel_n[k]; i++, q2++) {
                if (fmpq_is_zero(rec + q2)) continue;
                fmpq_t neg; fmpq_init(neg); fmpq_neg(neg, rec + q2);
                fmpq_mpoly_t T; fmpq_mpoly_init(T, QC);
                fmpq_mpoly_set_coeff_fmpq_ui(T, neg, ref.stair.ex + i * nv, QC);
                fmpq_mpoly_add(Q, Q, T, QC);
                fmpq_mpoly_clear(T, QC); fmpq_clear(neg);
            }
            fmpz_mpoly_init(B + nb, Dg->Z);
            /* Q = content * zpoly; take the integer polynomial and move it to
             * the DRL context term by term */
            ulong* te = Dg->tmp;
            for (slong t = 0; t < Q->zpoly->length; t++) {
                fmpz_mpoly_get_term_exp_ui(te, Q->zpoly, t, QC->zctx);
                fmpz_mpoly_push_term_fmpz_ui(B + nb, Q->zpoly->coeffs + t, te, Dg->Z);
            }
            fmpz_mpoly_sort_terms(B + nb, Dg->Z);
            fmpz_mpoly_combine_like_terms(B + nb, Dg->Z);
            nb++;
        }
        for (slong k = 0; k < nb && cert; k++) {
            fmpz_mpoly_t h; fmpz_mpoly_init(h, Dg->Z);
            fmpz_mpoly_set(h, B + k, Dg->Z);
            gf_reduce_z(Dg, h, -1, 1, NULL);
            cert = !Dg->expired && fmpz_mpoly_is_zero(h, Dg->Z);
            fmpz_mpoly_clear(h, Dg->Z);
        }
        fmpq_mpoly_clear(Q, QC); fmpq_mpoly_ctx_clear(QC);
        if (Dg->expired) ok = 0;
        if (cert) {
            /* primitive with a positive leading coefficient in Lex: the Lex
             * leading term is the relation's LM, coefficient positive after
             * clearing (Q's content carried the sign) */
            for (slong k = 0; k < nb; k++) {
                fmpz_t c; fmpz_init(c);
                _fmpz_vec_content(c, B[k].coeffs, B[k].length);
                /* the sign of the LM's coefficient decides */
                slong t0 = -1;
                ulong* te = Dg->tmp;
                for (slong t = 0; t < B[k].length && t0 < 0; t++) {
                    fmpz_mpoly_get_term_exp_ui(te, B + k, t, Dg->Z);
                    if (memcmp(te, ref.rel_ex + k * nv, (size_t)nv * sizeof(ulong)) == 0) t0 = t;
                }
                if (t0 >= 0 && fmpz_sgn(B[k].coeffs + t0) < 0) fmpz_neg(c, c);
                if (!fmpz_is_one(c)) fmpz_mpoly_scalar_divexact_fmpz(B + k, B + k, c, Dg->Z);
                fmpz_clear(c);
            }
            *out = B; *nout = nb; done = 1;
        } else {
            for (slong k = 0; k < nb; k++) fmpz_mpoly_clear(B + k, Dg->Z);
            flint_free(B);
        }
    }
    if (have) gff_clear(&ref);
    if (crt) _fmpz_vec_clear(crt, ncoef);
    if (rec) { for (slong k = 0; k < ncoef; k++) fmpq_clear(rec + k); flint_free(rec); }
    fmpz_clear(mod_prod);
    flint_free(Mp);
    _fmpz_vec_clear(num, nv * D * D);
    _fmpz_vec_clear(den, nv * D);
    gfm_clear(&S);
    return done;
}

/* An fmpz_mpoly (any context of nv variables) as a GBPoly shaped like `like`. */
static GBPoly* gf_zpoly_to_gbpoly(const fmpz_mpoly_t P, const fmpz_mpoly_ctx_t Z, ulong* ex,
                                  const GBPoly* like) {
    slong nv = like->n_vars;
    GBPoly* r = gb_poly_new(like->n_vars, like->order, like->elim_pivot);
    if (like->wmat) gb_poly_set_wmat(r, like->wmat);
    int* ie = malloc((size_t)(nv ? nv : 1) * sizeof(int));
    mpq_t c; mpq_init(c);
    gb_poly_reserve(r, (size_t)P->length);
    for (slong t = 0; t < P->length; t++) {
        fmpz_mpoly_get_term_exp_ui(ex, P, t, Z);
        fmpz_get_mpz(mpq_numref(c), P->coeffs + t);
        mpz_set_ui(mpq_denref(c), 1);
        for (slong v = 0; v < nv; v++) ie[v] = (int)ex[v];
        gb_poly_push_term(r, ie, c);
    }
    mpq_clear(c); free(ie);
    gb_poly_normalize(r);
    return r;
}

/* ------------------------------------------------------------------ */
/*  Entry                                                              */
/* ------------------------------------------------------------------ */

/* The common driver: F rational GBPolys (one order), p = 0 for Q. */
static int gf_basis(GBPoly* const* F, size_t n, ulong p, GBPoly*** out, size_t* out_n) {
    if (!gf_enabled() || n == 0) return 0;
    const GBPoly* f0 = NULL;
    for (size_t i = 0; i < n; i++) if (F[i]->n_terms) { f0 = F[i]; break; }
    if (!f0) { *out = NULL; *out_n = 0; return 1; }          /* the zero ideal */
    if (f0->n_vars < 1) return 0;
    int ord = gf_flint_order(f0->order, f0->wmat, f0->n_vars);
    if (ord < 0) return 0;
    for (size_t i = 0; i < n; i++)
        for (size_t t = 0; t < F[i]->n_terms * (size_t)F[i]->n_vars; t++)
            if (F[i]->exps[t] < 0) return 0;
    int timed_out = 0;
    GBPoly** R = NULL;
    slong nr = 0;
    /* Lex over Q in several variables.  A small system is cheapest computed
     * directly; the run is given a work budget, and only when it is exceeded
     * (the coefficient growth Lex suffers) does the basis come from the DRL
     * basis by FGLM, for a zero-dimensional ideal. */
    if (!p && ord == ORD_LEX && f0->n_vars >= 2 && !getenv("MATHILDA_GB_FGLM_OFF")) {
        Gf L; gf_init(&L, 0, 0, f0->n_vars, ORD_LEX);
        L.budget_pairs = GF_LEX_BUDGET_PAIRS; L.budget_bits = GF_LEX_BUDGET_BITS;
        slong* idx = NULL;
        int ok = gf_engine(&L, F, n, 0, &idx, &nr, &timed_out);
        if (ok) {
            int unit = 0;
            for (slong i = 0; i < nr && !unit; i++) unit = gf_is_unit(&L, idx[i]);
            if (unit) nr = 1;
            R = malloc((size_t)(nr ? nr : 1) * sizeof(GBPoly*));
            for (slong i = 0; i < nr; i++)
                R[i] = unit ? gf_one_gbpoly(f0) : gf_to_gbpoly(&L, idx[i], f0);
        }
        flint_free(idx);
        gf_clear(&L);
        if (timed_out) tc_check_deadline();
        if (ok) { *out = R; *out_n = (size_t)nr; return 1; }
        R = NULL; nr = 0;
    }
    if (!p && ord == ORD_LEX && f0->n_vars >= 2 && !getenv("MATHILDA_GB_FGLM_OFF")) {
        Gf Dg; gf_init(&Dg, 0, 0, f0->n_vars, ORD_DEGREVLEX);
        slong* idx = NULL; slong nd = 0;
        int ok = gf_engine(&Dg, F, n, 0, &idx, &nd, &timed_out);
        int handled = 0;
        if (ok) {
            int unit = 0;
            for (slong i = 0; i < nd && !unit; i++) unit = gf_is_unit(&Dg, idx[i]);
            if (unit) {
                R = malloc(sizeof(GBPoly*)); R[0] = gf_one_gbpoly(f0); nr = 1; handled = 1;
            } else {
                fmpz_mpoly_struct* B = NULL; slong nb = 0;
                if (gf_fglm(&Dg, idx, nd, &B, &nb)) {
                    R = malloc((size_t)(nb ? nb : 1) * sizeof(GBPoly*));
                    for (slong k = 0; k < nb; k++) {
                        R[k] = gf_zpoly_to_gbpoly(B + k, Dg.Z, Dg.tmp, f0);
                        fmpz_mpoly_clear(B + k, Dg.Z);
                    }
                    flint_free(B);
                    nr = nb; handled = 1;
                }
            }
        }
        flint_free(idx);
        gf_clear(&Dg);
        if (timed_out) tc_check_deadline();
        if (handled) { *out = R; *out_n = (size_t)nr; return 1; }
    }
    Gf G; gf_init(&G, p != 0, p, f0->n_vars, ord);
    slong* idx = NULL;
    int ok = gf_engine(&G, F, n, p, &idx, &nr, &timed_out);
    if (ok) {
        int unit = 0;
        for (slong i = 0; i < nr && !unit; i++) unit = gf_is_unit(&G, idx[i]);
        if (unit) nr = 1;
        R = malloc((size_t)(nr ? nr : 1) * sizeof(GBPoly*));
        for (slong i = 0; i < nr; i++)
            R[i] = unit ? gf_one_gbpoly(f0) : gf_to_gbpoly(&G, idx[i], f0);
    }
    flint_free(idx);
    gf_clear(&G);
    if (timed_out) tc_check_deadline();       /* unwinds; nothing is held */
    if (!ok) return 0;
    *out = R;
    *out_n = (size_t)nr;
    return 1;
}

int gb_flint_buchberger(GBPoly* const* F, size_t n, GBPoly*** out, size_t* out_n) {
    return gf_basis(F, n, 0, out, out_n);
}

int gb_flint_buchberger_mod(GBPoly* const* F, size_t n, uint64_t p,
                            GBPoly*** out, size_t* out_n) {
    if (p < 2 || p > (uint64_t)UWORD_MAX || !n_is_prime((ulong)p)) return 0;
    return gf_basis(F, n, (ulong)p, out, out_n);
}

#else  /* !USE_FLINT */

int gb_flint_buchberger(GBPoly* const* F, size_t n, GBPoly*** out, size_t* out_n) {
    (void)F; (void)n; (void)out; (void)out_n; return 0;
}
int gb_flint_buchberger_mod(GBPoly* const* F, size_t n, uint64_t p,
                            GBPoly*** out, size_t* out_n) {
    (void)F; (void)n; (void)p; (void)out; (void)out_n; return 0;
}

#endif /* USE_FLINT */
