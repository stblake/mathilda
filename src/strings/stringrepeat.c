/*
 * stringrepeat.c - StringRepeat builtin for Mathilda
 *
 * StringRepeat["str", n]        - "str" concatenated n times.
 * StringRepeat["str", n, max]   - up to n copies of "str", truncated so the
 *                                 total length is at most max (a partial final
 *                                 copy is allowed).
 *
 * Strings are treated as raw byte arrays (consistent with StringTake/StringDrop
 * and StringPartition across this subsystem); no UTF-8 codepoint decoding is
 * performed, so lengths count bytes.
 */

#include "picostrings.h"
#include "common.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/*
 * builtin_stringrepeat:
 * See file header for the accepted forms. Returns NULL (unevaluated) for a
 * non-string first argument, a non-integer or negative count/max, or when the
 * requested length would overflow. A wrong argument count emits
 * StringRepeat::argt (via builtin_arg_error) and leaves the call unevaluated.
 */
Expr* builtin_stringrepeat(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;

    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc < 2 || argc > 3)
        return builtin_arg_error("StringRepeat", argc, 2, 3);

    Expr* arg0 = args[0];
    if (arg0->type != EXPR_STRING) return NULL;

    /* Count n: non-negative integer. */
    if (args[1]->type != EXPR_INTEGER) return NULL;
    int64_t n = args[1]->data.integer;
    if (n < 0) return NULL;

    /* Optional maximum total length: non-negative integer. */
    bool have_max = false;
    int64_t max = 0;
    if (argc == 3) {
        if (args[2]->type != EXPR_INTEGER) return NULL;
        max = args[2]->data.integer;
        if (max < 0) return NULL;
        have_max = true;
    }

    const char* str = arg0->data.string;
    size_t len = strlen(str);

    /* Determine the output length. An empty base string or zero copies yields
     * the empty string; guarding len == 0 also avoids the i % len below. */
    size_t out_len;
    if (n == 0 || len == 0) {
        out_len = 0;
    } else {
        /* full = n * len, guarded against size_t overflow. */
        size_t un = (size_t)n;
        if (un > SIZE_MAX / len) {
            /* Overflow: only tolerable if a max cap keeps the result finite. */
            if (!have_max) return NULL;
            out_len = (size_t)max;  /* n*len exceeds SIZE_MAX, so cap wins. */
        } else {
            size_t full = un * len;
            out_len = (have_max && (size_t)max < full) ? (size_t)max : full;
        }
    }

    char* buf = malloc(out_len + 1);
    if (!buf) return NULL;

    /* Fill cyclically so a truncated final copy is handled naturally. */
    for (size_t i = 0; i < out_len; i++)
        buf[i] = str[i % len];
    buf[out_len] = '\0';

    Expr* result = expr_new_string(buf);
    free(buf);
    return result;
}
