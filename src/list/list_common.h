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

/* ------------------------------------------------------------------------- *
 * Exact numeric ordering for list builtins.
 *
 * These exist because the obvious public helpers are wrong for ordering by
 * VALUE, each in a way that is a silent wrong answer rather than a failure:
 *
 *   expr_compare        a canonical total order, not a numeric one. It settles
 *                       a value tie between different ExprTypes on the type
 *                       enum (sort.c:376), so 1 and 1.0 do not compare equal;
 *                       and for atoms not both integer-like it falls back to
 *                       get_numeric_value() doubles (sort.c:372-377), so two
 *                       exact rationals differing below double resolution do.
 *   is_real_numeric     routes rationals through is_rational
 *                       (arithmetic.c:105-117), which requires int64 numerator
 *                       AND denominator, so an ordinary exact 1/10^25 --
 *                       Rational[1, BigInt] -- is rejected.
 *   expr_numeric_sign   returns a bare 0 for both "genuinely zero" and "I do
 *                       not recognise this", and recognises neither MPFR nor a
 *                       bigint-component Rational.
 *
 * Promoted from src/list/nearest.c, where they were written for Nearest and
 * are now needed by a second caller.
 * ------------------------------------------------------------------------- */

/* True for a real number the helpers below can order exactly. Bigint-aware,
 * unlike is_real_numeric; NaN is rejected. Note this is not a superset of
 * is_real_numeric -- it also rejects a real-valued Complex, which that
 * accepts -- so the two are not interchangeable. */
bool list_real_number_q(Expr* e);

/* Sign of a real number: -1, 0, +1. Sets *ok false for anything no sign can be
 * read from, which callers turn into an unevaluated result. The flag is the
 * whole point: an unreadable value must not be reported as zero. */
int list_numeric_sign(Expr* e, bool* ok);

/* Order two reals by NUMERIC VALUE: -1, 0, +1. Sets *ok false when the
 * comparison cannot be decided.
 *
 * Exact throughout. Same-exactness pairs subtract (exact for exact operands);
 * a MIXED exact/inexact pair is lifted into mpq and compared there, because
 * subtracting would widen the pair to a double, lose the exact operand, and
 * make the comparator INTRANSITIVE -- with plain int64 values, 2^60 == 2.0^60
 * and 2^60 + 1 == 2.0^60 while 2^60 < 2^60 + 1 all held at once. */
int list_numeric_cmp(Expr* a, Expr* b, bool* ok);

#endif /* LIST_COMMON_H */
