/*
 * ff_hnf.h -- Hermite normal forms of polynomial matrices over K[u] (ffield
 * library, src/poly/ffield/; Phase 4.5 of the Bronstein completion campaign).
 *
 * Two kernels:
 *
 *  - ff_hnf_gr: the elimination over any FLINT generic field K (gr_ctx), column
 *    by column with unimodular extended-gcd steps, optionally modulo a monic D.
 *    It is the algorithm of Trager's hnf / hnfModDet (src/internal/trager/
 *    Normalize.m) and of the FLINT` polynomial-matrix builtins.
 *
 *  - ff_hnf_modular: the same Hermite form over K = Q[t]/(mu) (or Q) computed
 *    from images over the finite fields F_p[t]/(f), f an irreducible factor of
 *    mu mod p, combined by CRT and rational reconstruction, and CERTIFIED in
 *    characteristic 0 before it is returned.  Over Q and number fields the
 *    direct elimination modulo D bounds the degrees of the entries but not the
 *    heights of their rational coefficients: the extended-gcd cofactors grow,
 *    and a 9x3 matrix of degree 21 over Q(omega) took 27 s (G-perf-2) where the
 *    form itself has coefficients below 10^12.  The images cost microseconds
 *    each and the number of primes follows the height of the RESULT.
 *
 * The certificate (theory note, section E7).  Let L' = span(rows) + D K[u]^c,
 * D monic.  The candidate H is accepted only when
 *   (a) H has Hermite shape: upper triangular, monic pivots, every entry above
 *       a pivot of lower degree;
 *   (b) every generator of L' (each row, each D e_k) reduces to 0 by H, so
 *       L' is contained in the module L(H) spanned by H, and deg det L(H) <=
 *       deg det L';
 *   (c) for a prime p not dividing a denominator of the rows, of D or of mu,
 *       and an irreducible factor f of mu mod p, the Hermite form of the image
 *       of L' over F_p[t]/(f) has diagonal degrees summing to deg det L(H).
 * The K-dimension of L'/(D) is the rank of the matrix of the u^k r mod D, and a
 * rank only drops under the reduction mod (p, f), so deg det L' <= deg det of
 * the image; with (b) and (c), deg det L' = deg det L(H), hence L(H) = L', and
 * H is its Hermite form by uniqueness.  A reconstruction that fails any check
 * is discarded; the kernel declines (returns 0) rather than return an
 * unchecked form.
 *
 * Without FLINT nothing is declared; callers are compiled under USE_FLINT.
 */
#ifndef FF_HNF_H
#define FF_HNF_H

#ifdef USE_FLINT
#include <flint/flint.h>
#include <flint/fmpq_poly.h>
#include <flint/gr.h>
#include <flint/gr_poly.h>
#include "ff_modular.h"

/* Hermite form of the nr x nc matrix e (row-major gr_poly entries over the
 * field K), in place: the form in the first nc rows, the rows below zero.
 * With D != NULL (monic, degree >= 1) the rows D e_k must already be among
 * the rows, and every entry right of the current column is reduced mod D after
 * each step.  Returns 1, or 0 on rank deficiency or a FLINT failure. */
int ff_hnf_gr(gr_poly_struct* e, slong nr, slong nc, const gr_poly_struct* D, gr_ctx_t K);

/* The Hermite form H (nc*nc entries, row-major, initialised by the caller) of
 * span(rows) + D K[u]^nc over K = Q[t]/(mu), mu monic irreducible over Q, or
 * K = Q when mu is NULL.  rows is nr x nc; D is monic of degree >= 1.  Returns
 * 1 with H certified (see above), 0 to decline (the caller then runs the
 * direct elimination).  stats may be NULL; the prime budget is
 * ff_mm_set_max_primes. */
int ff_hnf_modular(ff_qtpoly_struct* H, const ff_qtpoly_struct* rows, slong nr, slong nc,
                   const ff_qtpoly_struct* D, const fmpq_poly_struct* mu, ff_mm_stats* stats);

#endif /* USE_FLINT */
#endif /* FF_HNF_H */
