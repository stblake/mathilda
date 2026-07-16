/*
 * stringjoin.c - StringJoin builtin for Mathilda
 *
 * StringJoin["s1", "s2", ...] - concatenates strings together
 */

#include "picostrings.h"
#include "sym_names.h"
#include <string.h>

/*
 * collect_strings:
 * Recursively collects all strings from an expression, flattening any List
 * wrappers encountered. Appends string pointers (borrowed, not copied) to
 * the output array. Returns true if all leaves are strings, false otherwise.
 */
static bool collect_strings(Expr* e, const char*** strs, size_t* count, size_t* cap) {
    if (e->type == EXPR_STRING) {
        if (*count >= *cap) {
            *cap *= 2;
            *strs = realloc(*strs, sizeof(const char*) * (*cap));
        }
        (*strs)[(*count)++] = e->data.string;
        return true;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (!collect_strings(e->data.function.args[i], strs, count, cap))
                return false;
        }
        return true;
    }
    return false;
}

/*
 * builtin_stringjoin:
 * StringJoin[s1, s2, ...] concatenates all string arguments into a single
 * string. Any List wrappers in the arguments are recursively flattened so
 * that all enclosed strings participate in the join.
 *
 * StringJoin[] with no arguments returns the empty string "".
 * Returns NULL (unevaluated) if any non-string, non-list leaf is found.
 */
Expr* builtin_stringjoin(Expr* res) {
    if (res->type != EXPR_FUNCTION)
        return NULL;

    size_t argc = res->data.function.arg_count;

    /* StringJoin[] -> "" */
    if (argc == 0)
        return expr_new_string("");

    /* Collect all strings, flattening lists */
    size_t cap = 16, count = 0;
    const char** strs = malloc(sizeof(const char*) * cap);

    for (size_t i = 0; i < argc; i++) {
        if (!collect_strings(res->data.function.args[i], &strs, &count, &cap)) {
            free(strs);
            return NULL;
        }
    }

    /* Compute total length */
    size_t total_len = 0;
    for (size_t i = 0; i < count; i++) {
        total_len += strlen(strs[i]);
    }

    /* Build concatenated string */
    char* buf = malloc(total_len + 1);
    buf[0] = '\0';
    size_t offset = 0;
    for (size_t i = 0; i < count; i++) {
        size_t slen = strlen(strs[i]);
        memcpy(buf + offset, strs[i], slen);
        offset += slen;
    }
    buf[total_len] = '\0';

    free(strs);

    Expr* result = expr_new_string(buf);
    free(buf);
    return result;
}
