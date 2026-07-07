/*
 * characters.c - Characters builtin for Mathilda
 *
 * Characters["string"] - gives a list of the individual characters in a string
 */

#include "picostrings.h"
#include "sym_names.h"
#include <string.h>

/*
 * builtin_characters:
 * Characters["string"] returns a List of single-character strings.
 * Each character in the input string becomes a length-1 string element.
 * Returns NULL (unevaluated) for non-string arguments.
 */
Expr* builtin_characters(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1)
        return NULL;

    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_STRING)
        return NULL;

    const char* str = arg->data.string;
    size_t len = strlen(str);

    Expr** chars = malloc(sizeof(Expr*) * len);
    if (!chars) {

        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    }

    char buf[2];
    buf[1] = '\0';
    for (size_t i = 0; i < len; i++) {
        buf[0] = str[i];
        chars[i] = expr_new_string(buf);
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), chars, len);
    free(chars);

    return result;
}
