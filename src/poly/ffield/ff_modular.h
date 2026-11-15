/*
 * ff_modular.h -- certified multi-modular computation over K = Q[t]/(mu) (or
 * Q), part of the ffield library (src/poly/ffield/, Phase 4.5).
 *
 * WHY.  Exact algorithms over Q and number fields (the Hermite form of a
 * polynomial matrix, the extended gcd, norms and resultants) suffer from
 * coefficient swell: the intermediate rationals are far larger than those of
 * the answer.  The answer is computed instead from its images over the finite
 * fields F_p[t]/(f), f an irreducible factor of mu mod p, combined over the
 * factors by the polynomial CRT, over the primes by the integer CRT, and
 * lifted by rational reconstruction.  The number of primes follows the height
 * of the ANSWER, and every image costs word-size arithmetic.
 *
 * TRUST.  A reconstruction is a candidate only.  It must first agree with one
 * more image (cheap), and is then CERTIFIED by an exact check in
 * characteristic 0 supplied by the problem (each kernel's header states its
 * certificate and why it suffices).  A problem whose certificate fails at
 * every attempt within the prime budget is declined: the caller runs its exact
 * direct algorithm.  Nothing unchecked is ever returned.
 *
 * SHAPES AND UNLUCKY PRIMES.  Each image carries a profile (a vector of
 * degrees, e.g. the diagonal degrees of a Hermite form or deg gcd).  For all
 * but finitely many primes the profile is the true one, and an unlucky prime
 * can only raise the profile's total.  The driver keeps the images of the least
 * total seen, restarts when a smaller one appears, and switches between profiles
 * of equal total only on repeated evidence.
 *
 * The driver is generic: a problem supplies its inputs (polynomials in u over
 * K, reduced mod p by the driver), an image routine, the positions of the
 * answer's coefficients, and a certificate.  ff_hnf.c and ff_xgcd.c are its
 * clients; norms and resultants by CRT are the next.
 */
#ifndef FF_MODULAR_H
#define FF_MODULAR_H

#ifdef USE_FLINT
#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/nmod_poly.h>
#include <flint/nf.h>
#include <flint/gr.h>
#include <flint/gr_poly.h>

/* An element of K[u], K = Q[t]/(mu): coefficient d of u is the t-polynomial
 * c[d] (degree < deg mu; a rational constant over Q).  len 0 is zero. */
typedef struct {
    slong len;
    fmpq_poly_struct* c;
} ff_qtpoly_struct;

void ff_qtpoly_init(ff_qtpoly_struct* a);
void ff_qtpoly_clear(ff_qtpoly_struct* a);
/* make a have exactly len coefficients, the new ones zero */
void ff_qtpoly_set_length(ff_qtpoly_struct* a, slong len);
/* strip trailing zero coefficients */
void ff_qtpoly_normalise(ff_qtpoly_struct* a);

/* Diagnostics of the last call (tests, benchmarks). */
typedef struct {
    slong primes;           /* images in the accepted reconstruction         */
    slong primes_tried;     /* every prime looked at                         */
    slong unlucky;          /* images discarded for their profile            */
    slong certificates;     /* exact certificate checks run                  */
} ff_mm_stats;

/* The prime budget of every multi-modular kernel (default 2000, about 37000 digits); a hook for the
 * tests.  m <= 0 restores the default. */
void ff_mm_set_max_primes(slong m);

/* One factor field F = F_p[t]/(f), f monic irreducible: nmod when deg f = 1
 * (elements are residues, t = root), fq_nmod otherwise. */
typedef struct {
    int deg1;
    ulong root;
    nmod_poly_t f;
    gr_ctx_t K;
    slong sz;
} ff_img_field;

/* elem <- a (a polynomial over F_p in t) mod f; out <- elem as a polynomial in t */
void ff_img_set(gr_ptr elem, const nmod_poly_t a, ff_img_field* G, nmod_poly_t scratch);
void ff_img_get(nmod_poly_t out, gr_srcptr elem, ff_img_field* G);
/* the K[u] polynomial with t-coefficients ap[0..len) (already mod p) over G */
int ff_img_poly(gr_poly_t P, const nmod_poly_struct* ap, slong len, ff_img_field* G, nmod_poly_t scratch);

/* The driver's inputs reduced mod p: input i is red[off[i] .. off[i] + len[i]). */
typedef struct {
    slong nin;
    slong* off;
    slong* len;
    slong total;
    nmod_poly_struct* red;
} ff_mm_reduced;

/* K itself, for the certificates: fmpq, or nf (Q[t]/(mu), deg mu >= 2). */
typedef struct {
    gr_ctx_t K;
    int isnf;
    nf_t nf;                /* conversions nf_elem <-> fmpq_poly             */
    slong sz;
} ff_k0;

/* over Q (deg mu <= 1) only the constant term of tp is read: tp must already be
 * reduced mod mu, as every kernel input is */
/* K for mu (monic irreducible, or NULL / degree 1: Q) */
void ff_k0_init(ff_k0* Z, const fmpq_poly_struct* mu);
void ff_k0_clear(ff_k0* Z);
int ff_k0_set(gr_ptr c, const fmpq_poly_t tp, ff_k0* Z);
int ff_k0_poly(gr_poly_t P, const ff_qtpoly_struct* a, ff_k0* Z);

typedef struct {
    void* ctx;
    const ff_qtpoly_struct* in;    /* the inputs                              */
    slong nin;
    const fmpq_poly_struct* mu;    /* monic IRREDUCIBLE over Q (a hypothesis of
                                      every certificate: K must be a field; the
                                      caller checks it), or NULL (Q)           */
    slong nprof;                   /* profile length (>= 1)                   */
    /* Size hint: the answer's height is expected to be at most about
     * height_mult * (input height + 64) bits (e.g. deg a + deg b for an extended
     * gcd, a resultant-type bound).  The prime budget is then raised to cover
     * twice that (rational reconstruction); early termination still stops at
     * the true size.  0: the fixed budget only. */
    slong height_mult;
    /* The image over G: fills prof, returns a handle (NULL when this prime is
     * bad for the problem, e.g. a degree dropped). */
    void* (*image)(void* ctx, const ff_mm_reduced* R, ff_img_field* G, slong* prof);
    /* the number of coefficient positions of an answer with this profile */
    slong (*npos)(void* ctx, const slong* prof);
    /* every coefficient of the image answer as a t-polynomial mod f: out[0..npos) */
    void (*extract)(void* ctx, void* img, ff_img_field* G, const slong* prof, nmod_poly_struct* out);
    void (*img_clear)(void* ctx, void* img, ff_img_field* G);
    /* Build the answer from cand (component j of position q at q*n + j, n =
     * deg mu) and check it exactly; 1 accepts it (the problem keeps it). */
    int (*certify)(void* ctx, const fmpq* cand, const slong* prof, ff_k0* Z);
} ff_mm_problem;

/* Run the driver.  1 when an answer was certified (the problem holds it),
 * 0 to decline.  stats may be NULL. */
int ff_mm_solve(const ff_mm_problem* P, ff_mm_stats* stats);

#endif /* USE_FLINT */
#endif /* FF_MODULAR_H */
