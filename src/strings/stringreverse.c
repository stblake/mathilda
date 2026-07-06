/*
 * stringreverse.c - StringReverse builtin for Mathilda
 *
 * StringReverse["string"] - reverses the order of the characters in a string
 */

#include "picostrings.h"
#include <string.h>
#include <stdio.h>

/*
 * sr_emit_argx:
 * Emit `StringReverse::argx: StringReverse called with N arguments; 1
 * argument is expected.` to stderr and return NULL so the call is left
 * unevaluated. Matches Mathematica's diagnostic for a wrong argument count.
 */
static Expr* sr_emit_argx(size_t argc) {
    fprintf(stderr,
            "StringReverse::argx: StringReverse called with %zu argument%s; "
            "1 argument is expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

/*
 * builtin_stringreverse:
 * StringReverse["string"] reverses the order of the characters in "string",
 * returning a new string. StringReverse is Listable, so the evaluator threads
 * it element-wise over a list argument before this builtin is reached; each
 * threaded call therefore sees a single string.
 *
 * A call with a number of arguments other than one emits StringReverse::argx
 * and is left unevaluated. A single non-string argument also leaves the call
 * unevaluated (NULL) so symbolic arguments flow through unchanged.
 */
Expr* builtin_stringreverse(Expr* res) {
    if (res->type != EXPR_FUNCTION)
        return NULL;

    size_t argc = res->data.function.arg_count;
    if (argc != 1)
        return sr_emit_argx(argc);

    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_STRING)
        return NULL;

    const char* str = arg->data.string;
    size_t len = strlen(str);

    char* buf = malloc(len + 1);
    if (!buf)
        return NULL;

    for (size_t i = 0; i < len; i++)
        buf[i] = str[len - 1 - i];
    buf[len] = '\0';

    Expr* result = expr_new_string(buf);
    free(buf);
    return result;
}
