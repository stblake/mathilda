#ifndef LIST_COMMON_H
#define LIST_COMMON_H

/* Shared infrastructure for the split list.c modules under src/list/.
 *
 * This header pulls in the common project and standard headers that the
 * original monolithic list.c included, so each per-function module only
 * needs to `#include "list_common.h"` plus its own public header. It also
 * declares the handful of small predicate/constructor helpers that are used
 * by more than one module (they live in list_common.c). */

#include "expr.h"
#include "common.h"
#include "symtab.h"
#include "eval.h"
#include "iter.h"
#include "core.h"
#include "arithmetic.h"
#include "print.h"
#include "attr.h"
#include "sym_intern.h"
#include "sym_names.h"

#include <inttypes.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

/* True for Overflow[]. */
bool is_overflow(Expr* e);
/* True for a List[...] expression. */
bool is_listq(Expr* e);
/* True for the symbol Infinity. */
bool is_infinity(Expr* e);
/* True for Times[-1, Infinity] (i.e. -Infinity). */
bool is_minus_infinity(Expr* e);
/* Build a fresh -Infinity expression (Times[-1, Infinity]). */
Expr* make_minus_infinity(void);
/* True for a real-valued number: integer, real, bigint, rational, MPFR, or a
 * complex with zero imaginary part. */
bool is_real_numeric(Expr* e);

#endif /* LIST_COMMON_H */
