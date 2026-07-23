#ifndef ITER_H
#define ITER_H

#include <stdbool.h>
#include "expr.h"
#include "symtab.h"   /* Rule, SymbolDef for the localization helpers */

Expr* builtin_do(Expr* res);
Expr* builtin_for(Expr* res);
Expr* builtin_while(Expr* res);
Expr* builtin_break(Expr* res);
Expr* builtin_continue(Expr* res);

void iter_init(void);

/*
 * ----------------------------------------------------------------------------
 *  Shared Wolfram iterator-spec helpers
 * ----------------------------------------------------------------------------
 *
 * A single structural parser + Block-style localization shared by every
 * iterator-driven builtin (Table, Do, Sum, ...).  This replaces the
 * near-identical spec-parsing blocks that previously lived inside each
 * builtin.  Each consumer decides for itself whether a parsed bound is
 * usable (numeric range vs. symbolic vs. Infinity), so the parser stays
 * policy-free: it only classifies and evaluates.
 */

typedef enum {
    ITER_KIND_COUNT,   /* bare n or {n}: repeat-count, var == NULL          */
    ITER_KIND_RANGE,   /* {i,imax} / {i,imin,imax} / {i,imin,imax,di}       */
    ITER_KIND_LIST     /* {i, {a,b,...}}: iterate over explicit list values */
} IterKind;

typedef struct {
    Expr*    var;    /* iterator symbol (owned copy); NULL for COUNT        */
    IterKind kind;
    Expr*    imin;   /* owned; RANGE only                                   */
    Expr*    imax;   /* owned; RANGE upper bound, or COUNT count            */
    Expr*    di;     /* owned; RANGE step (defaults to Integer 1)           */
    Expr*    list;   /* owned List[...] of values; LIST only                */
} IterSpec;

/*
 * Parse a held iterator spec into `out`.  Bound expressions are evaluate()d
 * so that symbolic bounds survive as symbols and constant arithmetic folds.
 * Returns true with `out` fully populated (caller must iter_spec_free), or
 * false for a malformed spec (variable slot not a symbol, or arity > 4), in
 * which case `out` is left freed/zeroed.
 */
bool iter_spec_parse(Expr* spec, IterSpec* out);

/* Free every Expr* owned by `s` and zero it.  NULL-safe. */
void iter_spec_free(IterSpec* s);

/*
 * Resolve a COUNT/RANGE spec to machine doubles for in-loop iteration.
 * `allow_inf` permits an Infinity upper bound (Do/Sum: true; Table: false).
 * Returns false when a bound is non-numeric (and not an allowed Infinity),
 * or the step is zero.  LIST specs are not numeric and return false.
 */
bool iter_spec_resolve_numeric(const IterSpec* s, bool allow_inf,
                               double* min_val, double* max_val,
                               double* di_val, bool* is_real, bool* is_inf);

/*
 * Block-style localization of the iterator variable.  shadow() saves and
 * clears the variable's current OwnValues (returning them); bind each
 * iteration's value with symtab_add_own_value(); restore() frees the
 * iteration-created OwnValue chain and reinstates the saved values.
 * `var` may be NULL (COUNT specs): both are no-ops then.
 */
Rule* iter_spec_shadow(Expr* var);
void  iter_spec_restore(Expr* var, Rule* saved_own);

#endif
