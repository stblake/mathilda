/*
 * besseljzero.c -- BesselJZero[n, k], the k-th positive zero of BesselJ[n, x].
 *
 * Currently a symbolic placeholder: it stays unevaluated for all arguments so
 * that the Hadamard-product recogniser (Product`BesselZero) can match the
 * canonical infinite product
 *
 *   Product[1 - x^2/BesselJZero[n,k]^2, {k,1,Inf}] = Gamma[n+1] (2/x)^n BesselJ[n,x].
 *
 * A future numeric path can evaluate BesselJZero[n, k] for exact numeric n and
 * positive-integer k via McMahon asymptotics + Newton refinement on BesselJ.
 *
 * Memory contract: takes ownership of res but must not free it; returns NULL
 * (leave unevaluated) or an owned closed form.
 */

#include "besseljzero.h"
#include "symtab.h"
#include "expr.h"
#include "attr.h"
#include <stdlib.h>

Expr* builtin_besseljzero(Expr* res) {
    /* Two-argument form only; stay symbolic for all arguments (for now). */
    (void)res;
    return NULL;
}

void besseljzero_init(void) {
    symtab_add_builtin("BesselJZero", builtin_besseljzero);
    symtab_get_def("BesselJZero")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED | ATTR_READPROTECTED);
    symtab_set_docstring("BesselJZero",
        "BesselJZero[n, k] gives the k-th positive zero of BesselJ[n, x]. "
        "Stays symbolic for symbolic arguments.");
}
