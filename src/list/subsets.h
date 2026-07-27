#ifndef SUBSETS_H
#define SUBSETS_H

#include "expr.h"

/* Subsets[list] / Subsets[list, nspec] / Subsets[list, nspec, s]
 *
 * Enumerates sublists of `list` ordered by increasing length, and
 * lexicographically by original element position within each length. See
 * subsets.c for the full semantics and the lazy generation strategy used by
 * the 3-argument (first-`s`) form. */
Expr* builtin_subsets(Expr* res);

#endif /* SUBSETS_H */
