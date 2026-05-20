#ifndef MATLSTSQ_H
#define MATLSTSQ_H

#include "expr.h"

/* Linear least-squares solver.
 *
 *   LeastSquares[m, b]
 *       Returns an x that minimises Norm[m . x - b].  When m has full
 *       column rank the minimiser is unique; when m is rank-deficient
 *       the result is the minimum-norm minimiser (the Moore-Penrose
 *       pseudoinverse solution PseudoInverse[m] . b).
 *
 *   LeastSquares[m, b, Method -> "<name>"]
 *   LeastSquares[m, b, Tolerance -> t]
 *       Optional Method and Tolerance arguments may appear in either
 *       order and may both be given.
 *
 * Supported Method names:
 *   "Automatic"            -- alias for "Direct" (default).
 *   "Direct"               -- Moore-Penrose solve via PseudoInverse.
 *                             Works on every input family (exact /
 *                             rational / symbolic / inexact / complex).
 *   "IterativeRefinement"  -- residual-correction loop on top of Direct
 *                             until ||dx||^2 <= tol^2 or 50 iterations.
 *                             For exact inputs the first correction is
 *                             exactly zero and the answer equals Direct;
 *                             for inexact inputs the iteration drives
 *                             round-off down to the configured Tolerance.
 *   "LSQR"                 -- Paige-Saunders LSQR iterative method.
 *                             Currently dispatched to Direct so the
 *                             user-facing API matches Mathematica; a
 *                             dedicated implementation is a future
 *                             extension.
 *   "Krylov"               -- iterative Krylov method (CGNR on the
 *                             normal equations).  Currently dispatched
 *                             to Direct.
 *
 *  The Tolerance option is forwarded to PseudoInverse and (in a future
 *  refinement) is also used as the convergence threshold for the
 *  iterative methods.
 *
 *  Registration is performed by `matlstsq_init`, called from
 *  `core_init` in `core.c`.
 */
Expr* builtin_leastsquares(Expr* res);
void  matlstsq_init(void);

#endif /* MATLSTSQ_H */
