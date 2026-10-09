#ifndef REDUCE_SOS_H
#define REDUCE_SOS_H

/*
 * reduce_sos.c -- Positivstellensatz / sum-of-squares emptiness certificate for
 * `Reduce[..., Reals]` (REDUCE_PLAN.md, D5).
 *
 * A LAST-RESORT stage: when the Fourier-Motzkin / CAD / zero-dimensional engines
 * all decline a real conjunction, this tries to PROVE the region empty by
 * exhibiting a certificate that one of its strict inequalities contradicts the
 * rest.  Concretely, for a region {g_i >= 0, h_j == 0, q < 0} it tries to show
 * q >= 0 on {g_i >= 0, h_j == 0} (so the region is empty), via a Putinar
 * representation
 *        q  ==  sigma_0 + sum_i sigma_i * g_i   (mod the ideal <h_j>)
 * with every sigma an exact sum of squares.  The sigmas are found numerically
 * (an SDP by alternating projection onto the PSD cone and the coefficient-
 * matching affine subspace, with numerical facial reduction for a boundary/
 * interior zero), then ROUNDED to rationals and the whole identity + PSD-ness is
 * re-verified in EXACT rational arithmetic.  Only a fully verified certificate
 * returns False; anything else declines (NULL).  Soundness over completeness --
 * it never returns a wrong answer and never emits a solution set.
 *
 * Needs LAPACK for the numeric guide (USE_LAPACK); without it the stage declines.
 */

#include "reduce_form.h"

struct Expr;

/* Returns a freshly-allocated `False` symbol if the region described by `F` (a
 * single conjunction over the variables `vlist[0..nv-1]`, interpreted over the
 * Reals) is PROVEN empty by a verified Positivstellensatz certificate; otherwise
 * NULL (decline).  `F` and `vlist` are borrowed. */
struct Expr* reduce_sos(const RForm* F, struct Expr** vlist, int nv);

#endif /* REDUCE_SOS_H */
