/*
 * stringpad.c - StringPadLeft and StringPadRight builtins for Mathilda
 *
 * StringPadLeft["str", n]        - "str" of length n, padded on the left with
 *                                  spaces or truncated (keeping the last n chars).
 * StringPadLeft["str", n, "p"]   - as above, padded by repeating copies of "p".
 * StringPadLeft["str"]           - returns "str" unchanged (n = length of str).
 * StringPadLeft[{s1, s2, ...}]   - pads each string with spaces to the length of
 *                                  the longest, so all become the same length.
 * StringPadLeft[{s1, ...}, n]    - pads or truncates each string to length n.
 * StringPadLeft[{s1, ...}, n, "p"] - as above, using padding string "p".
 *
 * StringPadRight is identical except padding is added on the right and
 * truncation keeps the first n characters.
 *
 * Padding is laid down as cyclic copies of the pad string read left-to-right
 * (pad[i] = p[i % plen]) on both sides, truncated when the target width is
 * reached; this matches the Wolfram Language documentation and all documented
 * examples (e.g. StringPadLeft["abcde", 10, "."] -> ".....abcde").
 *
 * Strings are treated as raw byte arrays (consistent with StringRepeat /
 * StringTake / StringPartition across this subsystem); no UTF-8 codepoint
 * decoding is performed, so lengths count bytes.
 */

#include "picostrings.h"
#include "common.h"       /* builtin_arg_error */
#include "sym_names.h"    /* SYM_List */
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

/*
 * pad_one:
 * Returns a fresh EXPR_STRING that is s padded (or truncated) to length n.
 *   left == true : pad on the left / truncate keeping the last n bytes.
 *   left == false: pad on the right / truncate keeping the first n bytes.
 * Padding bytes are cyclic copies of p (p[i % plen]). Returns NULL on
 * allocation failure or when padding is required (n > strlen(s)) but the pad
 * string is empty (plen == 0).
 */
static Expr* pad_one(const char* s, int64_t n, const char* p, size_t plen,
                     bool left) {
    size_t slen = strlen(s);

    if ((size_t)n <= slen) {
        /* Truncation (or exact fit): copy the appropriate n-byte window. */
        size_t start = left ? (slen - (size_t)n) : 0;
        char* buf = malloc((size_t)n + 1);
        if (!buf) return NULL;
        memcpy(buf, s + start, (size_t)n);
        buf[n] = '\0';
        Expr* result = expr_new_string(buf);
        free(buf);
        return result;
    }

    /* Padding: pad region of width w precedes (left) or follows (right) s. */
    size_t w = (size_t)n - slen;
    if (plen == 0) return NULL;  /* cannot build padding from an empty string */

    char* buf = malloc((size_t)n + 1);
    if (!buf) return NULL;

    if (left) {
        for (size_t i = 0; i < w; i++) buf[i] = p[i % plen];
        memcpy(buf + w, s, slen);
    } else {
        memcpy(buf, s, slen);
        for (size_t i = 0; i < w; i++) buf[slen + i] = p[i % plen];
    }
    buf[n] = '\0';

    Expr* result = expr_new_string(buf);
    free(buf);
    return result;
}

/*
 * pad_parse_n:
 * Reads the target length from args[1] into *out. Returns true if args[1] is a
 * non-negative EXPR_INTEGER, false otherwise (caller leaves the call
 * unevaluated).
 */
static bool pad_parse_n(Expr* arg, int64_t* out) {
    if (arg->type != EXPR_INTEGER) return false;
    if (arg->data.integer < 0) return false;
    *out = arg->data.integer;
    return true;
}

/*
 * pad_dispatch:
 * Shared implementation for StringPadLeft (left == true) and StringPadRight
 * (left == false). See file header for the accepted forms. Returns NULL to
 * leave the call unevaluated on any argument that does not match a supported
 * form; a wrong argument count emits head::argb via builtin_arg_error.
 */
static Expr* pad_dispatch(Expr* res, const char* head, bool left) {
    if (res->type != EXPR_FUNCTION) return NULL;

    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc < 1 || argc > 3)
        return builtin_arg_error(head, argc, 1, 3);

    /* Optional padding string (defaults to a single space). A List padding
     * (Wolfram's list form) is not supported -> leave unevaluated. */
    const char* p = " ";
    size_t plen = 1;
    if (argc == 3) {
        if (args[2]->type != EXPR_STRING) return NULL;
        p = args[2]->data.string;
        plen = strlen(p);
    }

    Expr* arg0 = args[0];

    /* --- List of strings ------------------------------------------------- */
    if (arg0->type == EXPR_FUNCTION &&
        arg0->data.function.head->type == EXPR_SYMBOL &&
        arg0->data.function.head->data.symbol.name == SYM_List) {

        size_t m = arg0->data.function.arg_count;
        Expr** elems = arg0->data.function.args;

        /* Every element must be a string. */
        for (size_t i = 0; i < m; i++)
            if (elems[i]->type != EXPR_STRING) return NULL;

        /* Target length: explicit n when given, else the longest element. */
        int64_t n;
        if (argc >= 2) {
            if (!pad_parse_n(args[1], &n)) return NULL;
        } else {
            size_t maxlen = 0;
            for (size_t i = 0; i < m; i++) {
                size_t li = strlen(elems[i]->data.string);
                if (li > maxlen) maxlen = li;
            }
            n = (int64_t)maxlen;
        }

        Expr** out = malloc(sizeof(Expr*) * (m ? m : 1));
        if (!out) return NULL;
        for (size_t i = 0; i < m; i++) {
            out[i] = pad_one(elems[i]->data.string, n, p, plen, left);
            if (!out[i]) {
                for (size_t j = 0; j < i; j++) expr_free(out[j]);
                free(out);
                return NULL;
            }
        }
        Expr* result = expr_new_function(expr_new_symbol(SYM_List), out, m);
        free(out);
        return result;
    }

    /* --- Single string --------------------------------------------------- */
    if (arg0->type != EXPR_STRING) return NULL;

    const char* s = arg0->data.string;
    int64_t n;
    if (argc >= 2) {
        if (!pad_parse_n(args[1], &n)) return NULL;
    } else {
        n = (int64_t)strlen(s);  /* 1-arg single string: no change */
    }

    return pad_one(s, n, p, plen, left);
}

/*
 * builtin_stringpadleft:
 * StringPadLeft pads on the left / truncates keeping the last n characters.
 */
Expr* builtin_stringpadleft(Expr* res) {
    return pad_dispatch(res, "StringPadLeft", true);
}

/*
 * builtin_stringpadright:
 * StringPadRight pads on the right / truncates keeping the first n characters.
 */
Expr* builtin_stringpadright(Expr* res) {
    return pad_dispatch(res, "StringPadRight", false);
}
