/*
 * stringlength.c - StringLength builtin for Mathilda
 *
 * StringLength["string"] - gives the number of characters in a string
 */

#include "picostrings.h"
#include <string.h>

/*
 * builtin_stringlength:
 * StringLength["string"] returns the number of characters in the string
 * as an integer. Returns NULL (unevaluated) for non-string arguments.
 */
Expr* builtin_stringlength(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1)
        return NULL;

    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_STRING)
        return NULL;

    int64_t len = (int64_t)strlen(arg->data.string);
    Expr* result = expr_new_integer(len);

    return result;
}
