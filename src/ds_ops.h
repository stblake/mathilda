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
} DsOps;

/* The ops row for a payload kind, or NULL for kinds with none (atoms, function
 * nodes).  O(1): a small array indexed by ExprType. */
const DsOps* ds_ops(ExprType t);

/* Register or replace the row for `ops->type` (the hook for new structure
 * kinds).  The row must outlive the process; returns false if the type is out of
 * range of the table. */
bool ds_register(const DsOps* ops);

#endif /* MATHILDA_DS_OPS_H */
