#ifndef SUM_INTERNAL_H
#define SUM_INTERNAL_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "expr.h"

/*
 * Shared helpers for the Sum sub-algorithms (Sum`Polynomial, Sum`Geometric,
 * Sum`Gosper).  Each takes/returns owned Expr* per the usual contract and
 * routes through the evaluator so existing builtins (ReplaceAll, Factor,
 * Together, ...) do the heavy lifting.
 */

/* Build head[args...] from `n` owned args, evaluate, free the call. */
Expr* sum_eval(const char* head, Expr** args, size_t n);

/* Evaluate e /. var -> val  (val and e are copied, not consumed). */
Expr* sum_subst(Expr* e, Expr* var, Expr* val);

/* Factor[e]; falls back to a copy of e if Factor leaves itself unevaluated.
 * e is not consumed. */
Expr* sum_factor(Expr* e);

/* True if e is free of var (FreeQ[e, var]). */
bool sum_free_of(Expr* e, Expr* var);

/* True if e contains an opaque indexed accessor (v[i] for an undefined v) whose
 * argument depends on var.  Used to bail out of Gosper and to skip the
 * closed-form cascade on a short finite range, where such a body has no closed
 * form and only churns Simplify.  Defined in sum_gosper.c. */
bool sum_body_has_opaque_index(Expr* e, Expr* var);

/* Convenience: integer node. */
Expr* sum_int(int64_t v);

/* a - b  (both copied), evaluated. */
Expr* sum_sub(Expr* a, Expr* b);

/* Parse the (f, i[, imin, imax]) argument shape shared by every stage.
 * Returns true and sets f and var (and, when definite, imin and imax) for
 * head[f,i] (indefinite, argc==2) or head[f,i,imin,imax] (definite, argc==4).
 * The out pointers alias res's args (not copied). */
bool sum_stage_args(Expr* res, Expr** f, Expr** var,
                    Expr** imin, Expr** imax, bool* definite);

#endif
