#ifndef STATS_COMMON_H
#define STATS_COMMON_H

/* Cross-cutting helpers shared between the per-builtin translation units of
 * src/stats/ (mirrors the numbertheory_internal.h pattern in
 * src/numbertheory/). Each helper is defined once in stats_common.c. */

#include "expr.h"
#include <stdbool.h>

/* True for a real/rational/complex numeric leaf. On success stores the real
 * approximation in *val (for real inputs) and whether the value is complex in
 * *out_complex; either out-pointer may be NULL. (Named out_complex, not
 * complex, because <complex.h> makes `complex` a reserved macro.) */
bool stats_is_numeric(Expr* e, double* val, bool* out_complex);

/* Apply the builtin named func_name column-wise to a matrix, i.e.
 * Map[func_name, Transpose[matrix]]. Returns NULL if the transpose is not a
 * function. */
Expr* stats_apply_columnwise(const char* func_name, Expr* matrix);

/* True iff e is NumericQ and free of the imaginary unit I (a real number). */
bool stats_is_real_numeric(Expr* e);

#endif
