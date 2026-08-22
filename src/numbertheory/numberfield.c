/*
 * numberfield.c -- see numberfield.h.
 *
 * Gate 1 of the Thue solver: certify (exactly) whether the equation order
 * Z[theta] of K = Q(theta) is the full ring of integers O_K.
 *
 *   - signature (r1, r2): an exact Sturm real-root count over Q[x]
 *     (FLINT-free, always available);
 *   - disc(f) and its full factorisation (FLINT + the facint helper);
 *   - Dedekind's criterion at each prime p with p^2 | disc(f), decided by
 *     exact finite-field arithmetic (no factoring needed -- the radical is
 *     obtained from gcd(f, f')).
 *
 * Every maximality verdict here is exact; there is no floating point.  On
 * any input the layer cannot fully decide (disc will not factor, a prime
 * too large for a single-word modulus, FLINT absent) the caller is told to
 * DECLINE.
 */
#include "numberfield_internal.h"
#include "facint.h"

#include <stdlib.h>
#include <string.h>

#ifdef USE_FLINT
#include <flint/flint.h>
#include <flint/fmpz_poly_factor.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>
#include <flint/arb.h>
#include <flint/arb_fmpz_poly.h>
/* (fmpz, fmpz_poly, fmpq, fmpq_poly, nf, nf_elem, acb come via the
 *  internal header.) */
#endif

bool nf_available(void) {
#ifdef USE_FLINT
    return true;
#else
    return false;
#endif
}

/* ================================================================== *
 *  Signature via an exact Sturm sequence over Q[x]  (FLINT-free).     *
 * ================================================================== */

typedef struct { int alloc; int deg; mpq_t* c; } QPoly;   /* c[i] = coef of x^i */

static QPoly* qp_new(int alloc) {
    QPoly* p = malloc(sizeof(QPoly));
    p->alloc = alloc < 1 ? 1 : alloc;
    p->c = malloc(sizeof(mpq_t) * (size_t)p->alloc);
    for (int i = 0; i < p->alloc; i++) mpq_init(p->c[i]);
    p->deg = -1;
    return p;
}

static void qp_free(QPoly* p) {
    if (!p) return;
    for (int i = 0; i < p->alloc; i++) mpq_clear(p->c[i]);
    free(p->c);
    free(p);
}

static void qp_normalize(QPoly* p) {
    p->deg = -1;
    for (int i = p->alloc - 1; i >= 0; i--)
        if (mpq_sgn(p->c[i]) != 0) { p->deg = i; break; }
}

/* Fresh r = a mod b (exact Euclidean division in Q[x]); b must be nonzero. */
static QPoly* qp_rem(const QPoly* a, const QPoly* b) {
    QPoly* r = qp_new(a->deg + 1 > 0 ? a->deg + 1 : 1);
    for (int i = 0; i <= a->deg; i++) mpq_set(r->c[i], a->c[i]);
    qp_normalize(r);
    mpq_t factor, t; mpq_init(factor); mpq_init(t);
    while (r->deg >= b->deg && r->deg >= 0) {
        mpq_div(factor, r->c[r->deg], b->c[b->deg]);
        int shift = r->deg - b->deg;
        for (int i = 0; i <= b->deg; i++) {
            mpq_mul(t, factor, b->c[i]);
            mpq_sub(r->c[i + shift], r->c[i + shift], t);
        }
        qp_normalize(r);
    }
    mpq_clear(factor); mpq_clear(t);
    return r;
}

bool nf_signature(const mpz_t* coeffs, int deg, int* r1) {
    if (deg < 1) return false;

    QPoly* g0 = qp_new(deg + 1);
    for (int i = 0; i <= deg; i++) mpq_set_z(g0->c[i], coeffs[i]);
    qp_normalize(g0);

    QPoly* g1 = qp_new(deg);                 /* f' */
    {
        mpq_t iq, tmp; mpq_init(iq); mpq_init(tmp);
        for (int i = 1; i <= deg; i++) {
            mpq_set_ui(iq, (unsigned long)i, 1);
            mpq_set_z(tmp, coeffs[i]);
            mpq_mul(g1->c[i - 1], tmp, iq);
        }
        mpq_clear(iq); mpq_clear(tmp);
    }
    qp_normalize(g1);

    QPoly** chain = malloc(sizeof(QPoly*) * (size_t)(deg + 2));
    int nchain = 0;
    chain[nchain++] = g0;
    chain[nchain++] = g1;

    while (chain[nchain - 1]->deg >= 0) {
        QPoly* r = qp_rem(chain[nchain - 2], chain[nchain - 1]);
        for (int i = 0; i < r->alloc; i++) mpq_neg(r->c[i], r->c[i]);
        qp_normalize(r);
        if (r->deg < 0) { qp_free(r); break; }
        chain[nchain++] = r;
        if (r->deg == 0) break;              /* constant: next remainder is 0 */
    }

    int Vplus = 0, Vminus = 0, prevp = 0, prevm = 0;
    for (int k = 0; k < nchain; k++) {
        int d = chain[k]->deg;
        if (d < 0) continue;                 /* defensive; should not occur */
        int sp = mpq_sgn(chain[k]->c[d]);            /* sign at +inf */
        int sm = (d % 2 == 0) ? sp : -sp;            /* sign at -inf */
        if (sp != 0) { if (prevp != 0 && sp != prevp) Vplus++; prevp = sp; }
        if (sm != 0) { if (prevm != 0 && sm != prevm) Vminus++; prevm = sm; }
    }

    for (int k = 0; k < nchain; k++) qp_free(chain[k]);
    free(chain);

    *r1 = Vminus - Vplus;
    return true;
}

/* ================================================================== *
 *  Dedekind's criterion  (needs FLINT).                               *
 * ================================================================== */

#ifdef USE_FLINT
static void build_fmpz_poly(fmpz_poly_t f, const mpz_t* coeffs, int deg) {
    fmpz_poly_zero(f);
    fmpz_t ci; fmpz_init(ci);
    for (int i = 0; i <= deg; i++) {
        fmpz_set_mpz(ci, coeffs[i]);
        fmpz_poly_set_coeff_fmpz(f, i, ci);
    }
    fmpz_clear(ci);
}

static void lift_nmod_to_fmpz(fmpz_poly_t out, const nmod_poly_t in) {
    fmpz_poly_zero(out);
    slong d = nmod_poly_degree(in);
    for (slong i = 0; i <= d; i++)
        fmpz_poly_set_coeff_ui(out, i, nmod_poly_get_coeff_ui(in, i));
}

/* FLINT 3.0 has no fmpz_poly_is_irreducible; derive it from the factorisation.
 * A degree->=2 integer polynomial is irreducible over Q iff its factorisation
 * is a single primitive factor to the first power spanning the full degree. */
static bool fmpz_poly_irreducible_q(const fmpz_poly_t f) {
    fmpz_poly_factor_t fac;
    fmpz_poly_factor_init(fac);
    fmpz_poly_factor(fac, f);
    bool irr = (fac->num == 1 && fac->exp[0] == 1
                && fmpz_poly_degree(fac->p) == fmpz_poly_degree(f));
    fmpz_poly_factor_clear(fac);
    return irr;
}
#endif

bool nf_dedekind_p_maximal(const mpz_t* coeffs, int deg, const mpz_t p, bool* ok) {
#ifndef USE_FLINT
    (void)coeffs; (void)deg; (void)p;
    *ok = false;
    return false;
#else
    *ok = true;
    if (!mpz_fits_ulong_p(p)) { *ok = false; return false; }  /* modulus too large */
    mp_limb_t pl = (mp_limb_t)mpz_get_ui(p);
    if (pl < 2) { *ok = false; return false; }

    fmpz_poly_t f; fmpz_poly_init(f);
    build_fmpz_poly(f, coeffs, deg);

    nmod_poly_t fbar, hbar, gbar, A, U, Tbar;
    nmod_poly_init(fbar, pl);
    fmpz_poly_get_nmod_poly(fbar, f);

    /* Radical gbar = prod of DISTINCT irreducible factors of fbar.  This must
     * come from the true mod-p factorisation, not gcd(f, f'): in small
     * characteristic (p | e_i) the derivative kills too much -- e.g. over F_2,
     * (x^2)' = 0 -- and gcd(f, f') then over-counts, silently mis-deciding
     * exactly the small primes Dedekind's criterion is run at. */
    nmod_poly_factor_t fac;
    nmod_poly_factor_init(fac);
    nmod_poly_factor(fac, fbar);
    nmod_poly_init(gbar, pl);
    nmod_poly_set_coeff_ui(gbar, 0, 1);                     /* gbar = 1 */
    for (slong i = 0; i < fac->num; i++) {
        nmod_poly_t t; nmod_poly_init(t, pl);
        nmod_poly_mul(t, gbar, &fac->p[i]);
        nmod_poly_swap(gbar, t);
        nmod_poly_clear(t);
    }
    nmod_poly_factor_clear(fac);

    nmod_poly_init(hbar, pl); nmod_poly_div(hbar, fbar, gbar); /* prod g_i^{e_i-1} */

    fmpz_poly_t G, H, GH, T; fmpz_poly_init(G); fmpz_poly_init(H);
    fmpz_poly_init(GH); fmpz_poly_init(T);
    lift_nmod_to_fmpz(G, gbar);
    lift_nmod_to_fmpz(H, hbar);
    fmpz_poly_mul(GH, G, H);
    fmpz_poly_sub(T, f, GH);                    /* f - G*H  (== 0 mod p) */
    fmpz_t pf; fmpz_init(pf); fmpz_set_mpz(pf, p);
    fmpz_poly_scalar_divexact_fmpz(T, T, pf);   /* (f - G*H) / p */

    nmod_poly_init(Tbar, pl); fmpz_poly_get_nmod_poly(Tbar, T);
    nmod_poly_init(A, pl); nmod_poly_gcd(A, gbar, hbar);   /* prod_{e_i>=2} g_i */
    nmod_poly_init(U, pl); nmod_poly_gcd(U, A, Tbar);

    bool maximal = (nmod_poly_degree(U) <= 0);   /* p-maximal iff gcd is a unit */

    nmod_poly_clear(fbar); nmod_poly_clear(hbar);
    nmod_poly_clear(gbar); nmod_poly_clear(A); nmod_poly_clear(U); nmod_poly_clear(Tbar);
    fmpz_poly_clear(G); fmpz_poly_clear(H); fmpz_poly_clear(GH); fmpz_poly_clear(T);
    fmpz_clear(pf); fmpz_poly_clear(f);
    return maximal;
#endif
}

/* ================================================================== *
 *  Whole-field maximality (Gate 1 driver).                            *
 * ================================================================== */

bool nf_is_maximal_order(const mpz_t* coeffs, int deg, bool* certified) {
#ifndef USE_FLINT
    (void)coeffs; (void)deg;
    *certified = false;
    return false;
#else
    *certified = true;

    fmpz_poly_t f; fmpz_poly_init(f);
    build_fmpz_poly(f, coeffs, deg);
    fmpz_t discf; fmpz_init(discf);
    fmpz_poly_discriminant(discf, f);
    mpz_t D; mpz_init(D); fmpz_get_mpz(D, discf);
    fmpz_clear(discf); fmpz_poly_clear(f);

    if (mpz_sgn(D) == 0) { mpz_clear(D); *certified = false; return false; }

    enum { CAP = 128 };
    mpz_t primes[CAP]; int64_t exps[CAP]; int cnt = 0;
    if (!facint_factor_complete(D, primes, exps, CAP, &cnt)) {
        mpz_clear(D); *certified = false; return false;   /* disc won't fully factor */
    }
    mpz_clear(D);

    bool maximal = true;
    for (int i = 0; i < cnt && maximal; i++) {
        if (exps[i] >= 2) {
            bool ok = true;
            bool pm = nf_dedekind_p_maximal(coeffs, deg, primes[i], &ok);
            if (!ok) { *certified = false; maximal = false; }   /* cannot decide */
            else if (!pm) { maximal = false; }                  /* proven non-maximal */
        }
    }
    for (int i = 0; i < cnt; i++) mpz_clear(primes[i]);
    return maximal;
#endif
}

/* ================================================================== *
 *  Construction / accessors.                                          *
 * ================================================================== */

NumberField* nf_field_create(const mpz_t* coeffs, int deg) {
#ifndef USE_FLINT
    (void)coeffs; (void)deg;
    return NULL;
#else
    if (deg < 2) return NULL;
    if (mpz_cmp_si(coeffs[deg], 1) != 0) return NULL;      /* must be monic */

    fmpz_poly_t f; fmpz_poly_init(f);
    build_fmpz_poly(f, coeffs, deg);
    if (!fmpz_poly_irreducible_q(f)) { fmpz_poly_clear(f); return NULL; }

    int r1 = 0;
    if (!nf_signature(coeffs, deg, &r1) || r1 < 0 || r1 > deg || (deg - r1) % 2 != 0) {
        fmpz_poly_clear(f); return NULL;
    }

    /* disc(f) and its full factorisation (drives both the Dedekind maximality
     * test and, when non-monogenic, Round 2). */
    fmpz_t discf; fmpz_init(discf); fmpz_poly_discriminant(discf, f);
    mpz_t Df; mpz_init(Df); fmpz_get_mpz(Df, discf);
    if (mpz_sgn(Df) == 0) { mpz_clear(Df); fmpz_clear(discf); fmpz_poly_clear(f); return NULL; }
    enum { PCAP = 128 };
    mpz_t primes[PCAP]; int64_t exps[PCAP]; int pcnt = 0;
    if (!facint_factor_complete(Df, primes, exps, PCAP, &pcnt)) {   /* disc won't factor -> decline */
        mpz_clear(Df); fmpz_clear(discf); fmpz_poly_clear(f); return NULL;
    }
    mpz_clear(Df);

    /* Maximality: Dedekind at each p with p^2 | disc. */
    bool all_maximal = true, decline = false;
    for (int i = 0; i < pcnt && !decline; i++) if (exps[i] >= 2) {
        bool okd = true;
        bool pm = nf_dedekind_p_maximal(coeffs, deg, primes[i], &okd);
        if (!okd) decline = true;                 /* cannot decide (e.g. p too big) */
        else if (!pm) all_maximal = false;
    }

    /* Integral basis (1/D)*W and d_K: monogenic fast path, else Round 2. */
    mpz_t* okW = malloc(sizeof(mpz_t) * (size_t)deg * (size_t)deg);
    for (int i = 0; i < deg * deg; i++) mpz_init(okW[i]);
    mpz_t okD, dK; mpz_init(okD); mpz_init(dK);
    int64_t index = 1;
    bool ib_ok = false;
    if (!decline) {
        if (all_maximal) {                        /* Z[theta] = O_K (unchanged path) */
            for (int i = 0; i < deg; i++) for (int j = 0; j < deg; j++)
                mpz_set_ui(okW[i*deg + j], (i == j) ? 1u : 0u);
            mpz_set_ui(okD, 1); fmpz_get_mpz(dK, discf); index = 1; ib_ok = true;
        } else {                                  /* Round 2 (Pohst-Zassenhaus) */
            ib_ok = nf_round2_maximal_order(coeffs, deg, (const mpz_t*)primes, exps, pcnt,
                                            okW, okD, dK, &index);
        }
    }
    for (int i = 0; i < pcnt; i++) mpz_clear(primes[i]);
    if (!ib_ok) {
        for (int i = 0; i < deg * deg; i++) mpz_clear(okW[i]);
        free(okW);
        mpz_clear(okD); mpz_clear(dK); fmpz_clear(discf); fmpz_poly_clear(f);
        return NULL;                              /* Gate 1 DECLINE */
    }

    NumberField* K = malloc(sizeof *K);
    K->deg = deg;
    K->r1 = r1;
    K->r2 = (deg - r1) / 2;
    mpz_init_set(K->disc, dK);                    /* d_K (== disc(f) when monogenic) */
    K->coeffs = malloc(sizeof(mpz_t) * (size_t)(deg + 1));
    for (int i = 0; i <= deg; i++) mpz_init_set(K->coeffs[i], coeffs[i]);

    /* O_K integral basis (Z[theta] = identity/1 in the monogenic case). */
    K->ok_num = okW;                              /* transfer ownership */
    mpz_init_set(K->ok_den, okD);
    K->ok_index = index;
    mpz_clear(okD); mpz_clear(dK);

    fmpq_poly_t fq; fmpq_poly_init(fq); fmpq_poly_set_fmpz_poly(fq, f);
    nf_init(K->nf, fq);
    fmpq_poly_clear(fq);
    fmpz_clear(discf); fmpz_poly_clear(f);

    K->prec = 0;
    K->emb = NULL;
    K->roots = NULL;
    if (!nf_ensure_prec(K, 256)) { nf_field_free(K); return NULL; }
    return K;
#endif
}

void nf_field_free(NumberField* K) {
    if (!K) return;
    mpz_clear(K->disc);
    for (int i = 0; i <= K->deg; i++) mpz_clear(K->coeffs[i]);
    free(K->coeffs);
    if (K->ok_num) { for (int i = 0; i < K->deg * K->deg; i++) mpz_clear(K->ok_num[i]); free(K->ok_num); }
    mpz_clear(K->ok_den);
#ifdef USE_FLINT
    if (K->emb)   { _acb_vec_clear(K->emb, K->r1 + K->r2); }
    if (K->roots) { _acb_vec_clear(K->roots, K->deg); }
    nf_clear(K->nf);
#endif
    free(K);
}

int  nf_degree(const NumberField* K)    { return K->deg; }
int  nf_r1(const NumberField* K)        { return K->r1; }
int  nf_r2(const NumberField* K)        { return K->r2; }
int  nf_unit_rank(const NumberField* K) { return K->r1 + K->r2 - 1; }
void nf_disc(const NumberField* K, mpz_t out) { mpz_set(out, K->disc); }

/* ================================================================== *
 *  FLINT-typed internals (embeddings, norm) -- see the internal hdr.  *
 * ================================================================== */
#ifdef USE_FLINT

bool nf_ensure_prec(struct NumberField* K, slong prec) {
    if (K->emb && K->prec >= prec) return true;
    if (prec < 64) prec = 64;

    fmpz_poly_t f; fmpz_poly_init(f);
    { fmpz_t ci; fmpz_init(ci);
      for (int i = 0; i <= K->deg; i++) {
          fmpz_set_mpz(ci, K->coeffs[i]);
          fmpz_poly_set_coeff_fmpz(f, i, ci);
      }
      fmpz_clear(ci); }

    slong wp = prec + 16;
    acb_ptr all = _acb_vec_init(K->deg);
    arb_fmpz_poly_complex_roots(all, f, 0, wp);   /* isolates each root to wp bits */
    fmpz_poly_clear(f);

    acb_ptr emb   = _acb_vec_init(K->r1 + K->r2);
    acb_ptr roots = _acb_vec_init(K->deg);
    int ri = 0, cj = K->r1;
    for (int i = 0; i < K->deg; i++) {
        if (arb_contains_zero(acb_imagref(all + i))) {          /* real embedding */
            if (ri < K->r1) { acb_set(emb + ri, all + i); acb_set(roots + ri, all + i); }
            ri++;
        } else if (arb_is_positive(acb_imagref(all + i))) {     /* complex rep (Im>0) */
            if (cj < K->r1 + K->r2) {
                acb_set(emb + cj, all + i);
                acb_set(roots + cj, all + i);
                acb_conj(roots + K->r2 + cj, all + i);          /* matching conjugate */
            }
            cj++;
        }
    }
    bool okcount = (ri == K->r1 && cj == K->r1 + K->r2);
    _acb_vec_clear(all, K->deg);
    if (!okcount) {                                             /* prec too low to classify */
        _acb_vec_clear(emb, K->r1 + K->r2);
        _acb_vec_clear(roots, K->deg);
        return false;
    }

    if (K->emb)   _acb_vec_clear(K->emb, K->r1 + K->r2);
    if (K->roots) _acb_vec_clear(K->roots, K->deg);
    K->emb = emb;
    K->roots = roots;
    K->prec = prec;
    return true;
}

void nf_embed_int(const struct NumberField* K, const mpz_t* c, int k, acb_t out) {
    slong wp = K->prec + 16;
    acb_zero(out);
    fmpz_t ci; fmpz_init(ci);
    for (int i = K->deg - 1; i >= 0; i--) {              /* Horner in theta^(k) */
        acb_mul(out, out, K->emb + k, wp);
        fmpz_set_mpz(ci, c[i]);
        acb_add_fmpz(out, out, ci, wp);
    }
    fmpz_clear(ci);
}

void nf_log_embedding_int(const struct NumberField* K, const mpz_t* c, arb_ptr out) {
    slong wp = K->prec + 16;
    acb_t v; acb_init(v);
    arb_t ab; arb_init(ab);
    for (int i = 0; i < K->r1 + K->r2; i++) {
        nf_embed_int(K, c, i, v);
        acb_abs(ab, v, wp);
        arb_log(out + i, ab, wp);
        if (i >= K->r1) arb_mul_ui(out + i, out + i, 2, wp);   /* complex weight */
    }
    arb_clear(ab); acb_clear(v);
}

bool nf_norm_int(const struct NumberField* K, const mpz_t* c, mpz_t out) {
    fmpq_poly_t g; fmpq_poly_init(g);
    fmpz_t z; fmpz_init(z);
    for (int i = 0; i < K->deg; i++) {
        fmpz_set_mpz(z, c[i]);
        fmpq_poly_set_coeff_fmpz(g, i, z);
    }
    fmpz_clear(z);
    nf_elem_t e; nf_elem_init(e, K->nf);
    nf_elem_set_fmpq_poly(e, g, K->nf);
    fmpq_t nrm; fmpq_init(nrm);
    nf_elem_norm(nrm, e, K->nf);
    bool isint = fmpz_is_one(fmpq_denref(nrm));
    if (isint) fmpz_get_mpz(out, fmpq_numref(nrm));
    fmpq_clear(nrm); nf_elem_clear(e, K->nf); fmpq_poly_clear(g);
    return isint;
}

const nf_t* nf_get_nf(const struct NumberField* K) { return (const nf_t*)&K->nf; }

void nf_ok_to_theta(const struct NumberField* K, const mpz_t* c, mpz_t* v) {
    int n = K->deg;
    for (int j = 0; j < n; j++) {
        mpz_set_ui(v[j], 0);
        for (int i = 0; i < n; i++)
            if (mpz_sgn(c[i]) != 0) mpz_addmul(v[j], c[i], K->ok_num[i*n + j]);
    }
}

const mpz_t* nf_ok_basis(const struct NumberField* K) { return (const mpz_t*)K->ok_num; }
void nf_ok_denom(const struct NumberField* K, mpz_t out) { mpz_set(out, K->ok_den); }
int64_t nf_ok_index(const struct NumberField* K) { return K->ok_index; }

#endif /* USE_FLINT */
