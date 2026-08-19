/*
 * numberfield_internal.h
 * ----------------------
 * FLINT-typed internals of the number-field layer, shared by
 * numberfield.c and the unit engine nfunits.c (and the Thue driver).
 * NOT part of the public API.
 *
 * The struct is defined in BOTH build configurations (so the always-
 * compiled accessors in numberfield.c see a complete type); the FLINT
 * members and helper declarations are guarded by USE_FLINT.  Under
 * USE_FLINT=0 no NumberField is ever constructed (nf_field_create
 * returns NULL), so the accessors are dead but must still compile.
 */
#ifndef NUMBERFIELD_INTERNAL_H
#define NUMBERFIELD_INTERNAL_H

#include "numberfield.h"
#include <gmp.h>

#ifdef USE_FLINT
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/nf.h>
#include <flint/nf_elem.h>
#include <flint/acb.h>
#endif

struct NumberField {
    int    deg;
    int    r1;          /* real embeddings                     */
    int    r2;          /* complex-conjugate embedding pairs   */
    mpz_t  disc;        /* d_K == disc(f) in the monogenic case */
    mpz_t* coeffs;      /* monic defining polynomial, [0..deg] */
#ifdef USE_FLINT
    nf_t     nf;        /* FLINT number field on f             */
    slong    prec;      /* current embedding precision (bits)  */
    acb_ptr  emb;       /* r1+r2 embeddings of theta: indices  */
                        /* 0..r1-1 real, r1..r1+r2-1 one rep   */
                        /* per complex pair (Im > 0).          */
    acb_ptr  roots;     /* ALL n conjugate roots of f, ordered */
                        /* [0..r1-1] real, [r1..r1+r2-1] Im>0, */
                        /* [r1+r2..n-1] the matching conjugates */
                        /* (roots[r1+r2+i] = conj roots[r1+i]).*/
#endif
};

/* Voronoi's algorithm for the fundamental unit of a RANK-1 complex cubic
 * (deg 3, signature (1,1)): the fallback for large-regulator fields whose
 * units exceed nfunits.c's coefficient box.  Writes the unit's power-basis
 * integer coordinates into out[0..2] (pre-init'd mpz) and returns true, or
 * false to DECLINE.  PROPOSAL only -- the caller must certify (p-saturation).
 * Returns false in a non-FLINT build.  (Defined in nfvoronoi.c.) */
bool nf_voronoi_unit_cubic11(struct NumberField* K, mpz_t* out);

#ifdef USE_FLINT

/* Ensure the conjugate embeddings are computed to at least `prec` bits
 * (recompute if the stored precision is lower).  Returns true on success. */
bool nf_ensure_prec(struct NumberField* K, slong prec);

/* Embedding of the algebraic integer  b = sum_{i=0}^{deg-1} c[i] theta^i
 * at conjugate index k in [0, r1+r2).  `out` is a pre-init'd acb; computed
 * at the field's current precision. */
void nf_embed_int(const struct NumberField* K, const mpz_t* c, int k, acb_t out);

/* Log-embedding vector of a unit given by integer coords c[0..deg-1]:
 * out[i] = m_i * log|sigma_i(b)|, m_i = 1 (i < r1) or 2 (complex), for
 * i in [0, r1+r2).  `out` is r1+r2 pre-init'd arb_t.  Precision = field's. */
void nf_log_embedding_int(const struct NumberField* K, const mpz_t* c, arb_ptr out);

/* Exact norm N_{K/Q}(sum c[i] theta^i) as an mpz.  Returns true and sets
 * `out` iff the norm is a rational integer (always so for algebraic
 * integers in Z[theta]). */
bool nf_norm_int(const struct NumberField* K, const mpz_t* c, mpz_t out);

const nf_t* nf_get_nf(const struct NumberField* K);

#endif /* USE_FLINT */
#endif /* NUMBERFIELD_INTERNAL_H */
