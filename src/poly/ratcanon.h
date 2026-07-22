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

#include <stddef.h>
#include "expr.h"

/* Register the Phase-1 prototype builtin.  Called from core_init(). */
void ratcanon_init(void);

/* ---- Phase 2: the canonical tower IR + builder ------------------------- *
 *
 * rat_canon_build(e) turns a rational expression into one ordered generator
 * tower: every transcendental kernel (Log/Exp/Tan/inverse-trig) and algebraic
 * generator (I, radicals, roots of unity) is replaced by a fresh indeterminate,
 * yielding num/den as a plain rational function in {base vars, generators}.
 * Transcendental generators are free (algebraically independent — logs expanded,
 * commensurate exponentials collapsed to a common fundamental); algebraic
 * generators carry an explicit relation (relation == 0, e.g. sym^2 + 1 for I,
 * sym^q - radicand for a radical).  Algebraic generators are ordered first
 * (leading) so their relations form a Groebner basis for the Phase-3 reduction.
 * No reduction and no builtin change here — this is the representation only.
 */

typedef enum { RCG_TRANSCENDENTAL, RCG_ALGEBRAIC } RcGenKind;

typedef struct {
    Expr*     kernel;    /* surface form to map back to: Log[u], E^fund, Sqrt[k], I, ... (owned) */
    char*     sym;       /* fresh generator symbol name "$rcgN$"                         (owned) */
    RcGenKind kind;
    Expr*     relation;  /* NULL (transcendental) or the minimal poly, == 0              (owned) */
} RcGen;

typedef struct {
    RcGen*  gens;        /* ordered: algebraic (leading) then transcendental */
    size_t  n, cap;
    Expr*   num;         /* numerator over {base vars, gen syms}   (owned) */
    Expr*   den;         /* denominator                            (owned) */
    Expr*   var;         /* heuristic main variable, or NULL       (owned) */
} RatCanonForm;

/* Build the tower IR for `e`.  Returns NULL only on allocation failure — an
 * expression with no kernels yields a form with n == 0 (plain rational). The
 * caller owns the result; free with rat_canon_free. */
RatCanonForm* rat_canon_build(const Expr* e);
void          rat_canon_free(RatCanonForm* f);

/* Reconstruct the original expression from a form by substituting each
 * generator symbol back to its kernel (algebraic relations applied by eval).
 * Used by the round-trip tests; returns a fresh owned Expr. */
Expr* rat_canon_roundtrip(const RatCanonForm* f);

/* Substitute every generator symbol in `e` back to its kernel (no eval).
 * Exposed for tests (e.g. checking a relation vanishes at its kernel). */
Expr* rat_canon_subst_back(const RatCanonForm* f, const Expr* e);

#endif /* MATHILDA_RATCANON_H */
