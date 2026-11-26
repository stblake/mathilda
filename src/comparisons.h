#ifndef COMPARISONS_H
#define COMPARISONS_H

#include "expr.h"

// Comparison routines
Expr* builtin_sameq(Expr* res);
Expr* builtin_unsameq(Expr* res);
Expr* builtin_equal(Expr* res);
Expr* builtin_unequal(Expr* res);
Expr* builtin_less(Expr* res);
Expr* builtin_greater(Expr* res);
Expr* builtin_lessequal(Expr* res);
Expr* builtin_greaterequal(Expr* res);

#include <stdbool.h>
/* The sign of a real exact numeric constant, decided as the comparison heads
 * decide e vs 0; false when undecidable (see comparisons.c). */
bool comparisons_real_constant_sign(Expr* e, int* sgn);

/* The certified numeric order of a and b (-1, 0, +1) that Less/Greater use:
 * exact for integers and rationals, MPFR at increasing precision for exact
 * quantities, +/-Infinity against every real. *decided is false when the pair
 * has no numeric order (non-numeric, complex, or inseparable). */
int comparisons_numeric_order(Expr* a, Expr* b, bool* decided);

/* SameQ[a, b]: structural identity, approximate numbers within
 * Internal`$SameQTolerance -- what a builtin whose default test is SameQ
 * (SplitBy's key comparison) must use instead of the exact expr_eq. */
bool expr_sameq(const Expr* a, const Expr* b);

/* SameQ classes over g keys that are already pairwise distinct under expr_eq
 * (the representatives an exact hash grouping produced, in first-appearance
 * order): cls[i] receives the class of key i, classes numbered 0, 1, ... in
 * order of their first key, and the class count is returned (g when no two
 * keys are SameQ -- the common case, decided without a sort when fewer than
 * two keys hold an approximate number).  Keys that differ only in approximate
 * numbers within Internal`$SameQTolerance share a class.  The tolerance is
 * not transitive, so it cannot be hashed: the keys holding approximate numbers
 * are sorted by an order that treats SameQ numbers as tied, and neighbours
 * that are SameQ are joined -- O(g log g) on top of the exact grouping.
 * `key(ctx, i)` returns key i, borrowed. */
size_t sameq_merge_classes(const Expr* (*key)(void* ctx, size_t i), void* ctx,
                           size_t g, size_t* cls);

// Initialize comparison builtins
void comparisons_init(void);

#endif // COMPARISONS_H
