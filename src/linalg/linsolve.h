#ifndef MATSOL_H
#define MATSOL_H

#include <stdint.h>
#include "expr.h"

/* RowReduce[m]                                  -- default Method -> Automatic.
 * RowReduce[m, Method -> "<name>"]              -- explicit dispatcher.
 *
 * LinearSolve[m, b]                             -- default Method -> Automatic.
 * LinearSolve[m, b, Method -> "<name>"]         -- explicit dispatcher.
 *
 * Supported method names:
 *   "Automatic"                  -- aliased to "DivisionFreeRowReduction"
 *   "DivisionFreeRowReduction"   -- Bareiss-like fraction-free Gauss-Jordan
 *   "OneStepRowReduction"        -- classical Gauss-Jordan with division
 *   "CofactorExpansion"          -- Cramer's rule (LinearSolve); identity-if-
 *                                   invertible (RowReduce), with fallback to
 *                                   DivisionFreeRowReduction on
 *                                   singular / rectangular input.
 *                                   For Inverse: adjugate / det formula.
 *
 * Method -> Automatic (the symbol) is also accepted, matching Mathematica.
 *
 * The MatsolMethod enum and helpers below are shared with src/matinv.c
 * so that Inverse, RowReduce and LinearSolve all parse the same
 * Method-option grammar and reuse the same Together / exact-division
 * primitives used inside the OneStep workers.
 */
typedef enum {
    MATSOL_AUTOMATIC = 0,
    MATSOL_DIVFREE,
    MATSOL_ONESTEP,
    MATSOL_COFACTOR,
    MATSOL_INVALID
} MatsolMethod;

/* Parse a `Method -> <value>` rule.  Returns MATSOL_INVALID when the
 * argument isn't a Method rule with a recognised RHS. */
MatsolMethod matsol_parse_method_option(Expr* opt);

/* Rate-limit a per-call warning so test loops don't spew.  `key` is
 * hashed to detect repeated invocations of the same call. */
void matsol_warn_once(uint64_t* last_hash, Expr* key, const char* msg);

/* Canonicalise a matrix entry via `Together` so subsequent
 * is_zero_poly checks behave correctly on symbolic rationals. */
Expr* matsol_canon_entry(Expr* e);

/* num / den, canonicalised.  Prefers exact polynomial division via
 * exact_div_wrapper when applicable, else builds Times[num, Power[den,
 * -1]] then Together-canonicalises. */
Expr* matsol_div_entry(Expr* num, Expr* den);

Expr* builtin_rowreduce(Expr* res);
Expr* builtin_linearsolve(Expr* res);
void  matsol_init(void);

#endif /* MATSOL_H */
