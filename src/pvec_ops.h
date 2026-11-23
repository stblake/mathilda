/* pvec_ops.h -- PVector at the expression level: constructor, the evaluator
 * gate, and the builtins (pvec_ops.c). The persistent-vector core is pvec.h. */
#ifndef MATHILDA_PVEC_OPS_H
#define MATHILDA_PVEC_OPS_H

#include "expr.h"
#include <stdbool.h>
#include <stdint.h>

void pvec_init(void);

/* A new PVector expression holding the elements of a List or packed List `e`
 * (borrowed); NULL when `e` is not one. New reference. */
Expr* pvec_make(const Expr* e);

/* The evaluator gate (eval.c).  `res` is an evaluated call whose head is
 * `head` (interned name) with attribute bits `attrs`.  If an argument is a
 * PVector and the head is one this layer handles, returns the evaluated answer
 * (a new reference; `res` is untouched), else NULL.  Native O(1) / O(log n)
 * operations run on the vector itself; structure-preserving heads (Map, Select,
 * Take, Join, Reverse, ..., listable arithmetic) run on the equivalent List and
 * the result is wrapped back into a PVector, so every answer equals PVector of
 * the List answer by construction; other sequence heads (Total, Max, Count, ...)
 * answer exactly as they do for the List. */
Expr* pvec_gate(Expr* res, const char* head, uint32_t attrs);

#endif
