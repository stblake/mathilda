/*
 * stringpad.c - StringPadLeft and StringPadRight builtins for Mathilda
 *
 * StringPadLeft["str", n]        - "str" of length n, padded on the left with
 *                                  spaces or truncated (keeping the last n chars).
 * StringPadLeft["str", n, "p"]   - as above, padded by repeating copies of "p".
 * StringPadLeft["str"]           - unevaluated (StringPadLeft::strlist).
 * StringPadLeft[{s1, s2, ...}]   - pads each string with spaces to the length of
 *                                  the longest, so all become the same length.
 * StringPadLeft[{s1, ...}, n]    - pads or truncates each string to length n.
 * StringPadLeft[{s1, ...}, n, "p"] - as above, using padding string "p".
 *
 * StringPadRight is identical except padding is added on the right and
 * truncation keeps the first n characters.
 *
 * Padding is the pad string repeated cyclically and ALIGNED to the result, as
 * in Mathematica: on the right, result character i (0-based) is p[i mod plen];
 * on the left, it is p[(i - n) mod plen], i.e. the cycle is anchored at the
 * result's end. So StringPadRight["ab", 6, "xyz"] is "abzxyz" and
 * StringPadLeft["ab", 7, "xyz"] is "zxyzxab".
 *
 * Lengths count characters (Unicode code points); the pad string is cycled a
 * character at a time.
 */

#include "picostrings.h"
#include "strings_internal.h"
#include "utf8.h"
#include "common.h"       /* builtin_arg_error */
#include "sym_names.h"    /* SYM_List */
#include "linkedlist.h"   /* a LinkedList of strings */
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

/*
 * pad_one:
 * Returns a fresh EXPR_STRING that is s padded (or truncated) to n characters.
 *   left == true : pad on the left / truncate keeping the last n characters.
 *   left == false: pad on the right / truncate keeping the first n characters.
 * `pv` is the pad string's character view (pv->nchars >= 1). Returns NULL on
 * allocation failure.
 */
static Expr* pad_one(const char* s, int64_t n, const U8View* pv, bool left) {
    size_t sbytes = strlen(s);
    size_t slen = u8_length(s, sbytes);

    if ((size_t)n <= slen) {
        /* Truncation (or exact fit): copy the appropriate n-character window. */
        size_t from = left ? u8_offset(s, sbytes, slen - (size_t)n) : 0;
        size_t to   = left ? sbytes : u8_offset(s, sbytes, (size_t)n);
        return expr_new_string_n(s + from, to - from);
    }

    /* Padding: w characters before (left) or after (right) s. */
    size_t w = (size_t)n - slen;
    size_t plen = pv->nchars;
    size_t pbytes = 0;                         /* bytes in the pad region */
    size_t start = left ? (size_t)(((int64_t)plen - (n % (int64_t)plen)) % (int64_t)plen)
                        : slen % plen;         /* cycle index of the first pad char */
    for (size_t k = 0; k < w; k++) {
        size_t c = (start + k) % plen;
        pbytes += u8_view_off(pv, c + 1) - u8_view_off(pv, c);
    }

    char* buf = malloc(sbytes + pbytes + 1);
    if (!buf) return NULL;
    char* w0 = left ? buf : buf + sbytes;      /* where the padding goes */
    size_t o = 0;
    for (size_t k = 0; k < w; k++) {
        size_t c = (start + k) % plen;
        size_t a = u8_view_off(pv, c), b = u8_view_off(pv, c + 1);
        memcpy(w0 + o, pv->s + a, b - a);
        o += b - a;
    }
    memcpy(left ? buf + pbytes : buf, s, sbytes);
    buf[sbytes + pbytes] = '\0';
    return expr_new_string_adopt(buf);
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

    /* Optional padding string (defaults to a single space). */
    const char* p = " ";
    if (argc == 3) {
        if (args[2]->type != EXPR_STRING || args[2]->data.string[0] == '\0')
            return strings_arg_message(head, "stringnz", "String of non-zero length",
                                       3, res);
        p = args[2]->data.string;
    }
    if (argc >= 2) {
        int64_t probe;
        if (!pad_parse_n(args[1], &probe))
            return strings_arg_message(head, "intnm",
                                       "Non-negative machine-sized integer", 2, res);
    }
    U8View pv;
    if (u8_view_init(&pv, p, strlen(p)) != 0) return NULL;

    Expr* arg0 = args[0];
    Expr* result = NULL;

    LLElems sv;
    if (ll_elems_open(arg0, &sv)) {
        /* --- List (or LinkedList) of strings: pad each, to n or to the longest */
        size_t m = sv.n;
        Expr** elems = sv.v;
        for (size_t i = 0; i < m; i++)
            if (elems[i]->type != EXPR_STRING) {
                ll_elems_close(&sv);
                u8_view_free(&pv);
                return strings_arg_message(head, "strlist", "List of strings", 1, res);
            }
        if (m == 0) {   /* Mathematica: StringPadLeft[{}, n] is "", StringPadLeft[{}] is {} */
            ll_elems_close(&sv);
            u8_view_free(&pv);
            return argc >= 2 ? expr_new_string("")
                             : ll_seq_like(arg0, NULL, 0);
        }

        int64_t n = 0;
        if (argc >= 2) {
            pad_parse_n(args[1], &n);
        } else {
            size_t maxlen = 0;
            for (size_t i = 0; i < m; i++) {
                const char* si = elems[i]->data.string;
                size_t li = u8_length(si, strlen(si));
                if (li > maxlen) maxlen = li;
            }
            n = (int64_t)maxlen;
        }

        Expr** out = malloc(sizeof(Expr*) * m);
        size_t built = 0;
        for (; out && built < m; built++) {
            out[built] = pad_one(elems[built]->data.string, n, &pv, left);
            if (!out[built]) break;
        }
        if (out && built == m) result = ll_seq_like(arg0, out, m);
        else for (size_t k = 0; out && k < built; k++) expr_free(out[k]);
        free(out);
        ll_elems_close(&sv);
    } else if (arg0->type == EXPR_STRING && argc >= 2) {
        /* --- Single string ---------------------------------------------- */
        int64_t n = 0;
        pad_parse_n(args[1], &n);
        result = pad_one(arg0->data.string, n, &pv, left);
    } else {
        /* A single string needs a target length; anything else is not a
         * list of strings. */
        result = strings_arg_message(head, "strlist", "List of strings", 1, res);
    }
    u8_view_free(&pv);
    return result;
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
