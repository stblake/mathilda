#ifndef MATHILDA_ASSOC_QUERY_H
#define MATHILDA_ASSOC_QUERY_H

/* ---------------------------------------------------------------------------
 * Query and Dataset -- the Wolfram Language's query language over nested
 * lists and associations (src/assoc_query.c).
 *
 *   Query[op1, op2, ...][data]   applies op_i at level i of data
 *   Dataset[data][op1, ...]      runs the same query and re-wraps a
 *                                List/Association result as a Dataset
 *
 * Both are *callable objects*: the head of the call is the compound
 * expression Query[...] / Dataset[...].  The evaluator reaches them through
 * query_callable_probe / query_callable_apply, in the same composite-head
 * chain as Function[...][args] and <|...|>[key].
 * -------------------------------------------------------------------------- */

#include "expr.h"
#include <stdbool.h>
#include <stddef.h>

/* True when `head` is Query[...], Dataset[...] or GeneralUtilities`Slice[...]
 * -- a compound head this module knows how to apply.  Cheap (pointer compare). */
bool query_callable_probe(const Expr* head);

/* Apply a Query / Dataset / Slice head to `args` (borrowed).  Returns an owned
 * result, or NULL to leave the call unevaluated (wrong arity). */
Expr* query_callable_apply(const Expr* head, Expr** args, size_t nargs);

/* The REPL text form of a Dataset: a small aligned table for tabular data
 * (a list or association of records, a list or association of scalars),
 * falling back to Dataset[<expr>] otherwise.  Heap string (caller frees), or
 * NULL when `ds` is not a well-formed Dataset[data]. */
char* dataset_format(const Expr* ds);

/* Registers Query and Dataset and installs the Dataset-transparent wrappers
 * on Normal/Length/Keys/Values/Dimensions/Select/SortBy/Map/First/Last/Take/
 * Part/Reverse.  Must run after every module that registers those builtins. */
void assoc_query_init(void);

#endif /* MATHILDA_ASSOC_QUERY_H */
