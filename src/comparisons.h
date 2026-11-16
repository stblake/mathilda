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

// Initialize comparison builtins
void comparisons_init(void);

#endif // COMPARISONS_H
