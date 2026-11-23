#ifndef MATHILDA_DS_OPS_H
#define MATHILDA_DS_OPS_H

/* ---------------------------------------------------------------------------
 * ds_ops -- one table of operations per data-structure kind.
 *
 * Mathilda's non-function payload kinds (packed arrays, sparse arrays,
 * LinkedLists, and the kinds still to come) used to be wired into every generic
 * builtin by hard-coded `type == EXPR_X` branches: Length, Dimensions, Head and
 * the pattern matcher's head test each repeated the same three-way split.  A new
 * structure then meant editing every one of those sites (and, separately, about
 * 25 Compile[] sites) before it behaved like a sequence at all.
 *
 * This table is the single place a kind states what it can do.  A generic
 * builtin asks `ds_ops(e->type)` and calls the operation if the kind provides it
 * (a NULL slot means "this kind does not have that notion"), so registering a
 * kind makes it work in every builtin that has been migrated, with no edit to
 * the builtin.  This slice migrates head / length / dimensions; element access,
 * iteration and rebuild follow in later slices, then the Compile[] registry.
 *
 * Behaviour is deliberately identical to the branches it replaces.
 * -------------------------------------------------------------------------- */

#include "expr.h"
#include <stdbool.h>
#include <stdint.h>

typedef struct DsOps {
    ExprType type;           /* the payload kind this row describes */
    const char* name;        /* for diagnostics and tests */

    /* The Head[] of a value of this kind (owned symbol), never NULL. */
    Expr* (*head)(const Expr* e);
    /* Does `e` match the blank head symbol name `hn` (interned)?  e.g. _SparseArray. */
    bool (*head_is)(const Expr* e, const char* hn);
    /* Leading-axis length, as Length[] reports it. */
    int64_t (*length)(const Expr* e);
    /* Write up to `max` leading dimensions into `out`; return how many (the rank
     * when it fits).  NULL = the kind has no rectangular shape. */
    int (*dims)(const Expr* e, int64_t* out, int max);
    /* Element / selection access for Part[e, i, ...].  Returns the result (NULL
     * with *st == DS_PART_DONE = leave Part unevaluated).  *st is set to
     * DS_PART_DEGRADE when the spec must be met on the boxed List path
     * (caller delists and re-evaluates) and DS_PART_DECLINE when this kind does
     * not take the spec (caller falls through, e.g. [[0]] head extraction). */
    Expr* (*part)(Expr* e, Expr** idx, size_t n, int* st);
    /* Materialise as an ordinary expression (what Normal[] returns): a new owned
     * value, NULL on failure.  NULL slot = the kind has no such conversion. */
    Expr* (*to_list)(const Expr* e);
} DsOps;

enum { DS_PART_DONE = 0, DS_PART_DEGRADE = 1, DS_PART_DECLINE = 2 };

/* Dispatch Part through the row; *st as above.  Returns NULL with
 * *st == DS_PART_DECLINE when the kind has no row or no part slot. */
Expr* ds_part(Expr* e, Expr** idx, size_t n, int* st);

/* Normal[] through the row.  *handled is false (and NULL returned) when the
 * kind has no to_list slot, so the caller can fall through. */
Expr* ds_to_list(const Expr* e, bool* handled);

/* The ops row for a payload kind, or NULL for kinds with none (atoms, function
 * nodes).  O(1): a small array indexed by ExprType. */
const DsOps* ds_ops(ExprType t);

/* Register or replace the row for `ops->type` (the hook for new structure
 * kinds).  The row must outlive the process; returns false if the type is out of
 * range of the table. */
bool ds_register(const DsOps* ops);

#endif /* MATHILDA_DS_OPS_H */
