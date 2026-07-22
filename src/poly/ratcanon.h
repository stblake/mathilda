/* ratcanon.h — unified rational normalization for Together/Cancel (see
 * RATCANON_REWRITE_PLAN.md).  The target design: build ONE differential/
 * algebraic tower for the input, reduce once via FLINT, map back with one
 * output convention — replacing the ~8-path engine zoo.
 *
 * PHASE 1 (current): a throwaway prototype `RatCanonPrototype[expr]` that proves
 * the one-front-end + one-reduction pipeline on the four representative regimes
 * (plain-Q, transcendental-kernel, Q(i), symbolic Sqrt[k]).  It substitutes
 * every kernel / algebraic constant / radical to a fresh free symbol, reduces
 * with the existing plain-Q engine (flint_rational_together), maps back, and
 * lets `eval` apply the algebraic relations (I^2->-1, Sqrt[k]^2->k).  The real
 * builder (rat_canon_build) and reduction (rat_canon_reduce) land in Phases 2-3;
 * this prototype is deleted at the end of Phase 3.
 */
#ifndef MATHILDA_RATCANON_H
#define MATHILDA_RATCANON_H

/* Register the Phase-1 prototype builtin.  Called from core_init(). */
void ratcanon_init(void);

#endif /* MATHILDA_RATCANON_H */
