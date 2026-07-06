#ifndef STRINGS_INTERNAL_H
#define STRINGS_INTERNAL_H

/*
 * strings_internal.h - helpers shared between the per-builtin translation
 * units under src/strings/. Kept as static inline functions so no separate
 * object file is needed; each user simply includes this header.
 */

#include <stdbool.h>
#include <stdint.h>
#include "expr.h"
#include "sym_names.h"

/*
 * strings_is_upto:
 * Returns true if e is UpTo[n] where n is an integer, and writes the integer
 * value to *out. Shared by StringTake and StringDrop, both of which accept an
 * UpTo[n] count specification.
 */
static inline bool strings_is_upto(Expr* e, int64_t* out) {
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_UpTo &&
        e->data.function.arg_count == 1 &&
        e->data.function.args[0]->type == EXPR_INTEGER) {
        *out = e->data.function.args[0]->data.integer;
        return true;
    }
    return false;
}

#endif // STRINGS_INTERNAL_H
