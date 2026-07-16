/*
 * stringtake.c - StringTake builtin for Mathilda
 *
 * StringTake["string", spec] - extracts substrings by count, range, or step
 */

#include "picostrings.h"
#include "strings_internal.h"
#include "eval.h"
#include "sym_names.h"
#include <string.h>

/*
 * stringtake_substring:
 * Helper to extract a substring from str[start-1..end-1] (1-based indices,
 * already resolved to positive). Returns a new EXPR_STRING.
 */
static Expr* stringtake_substring(const char* str, int64_t start, int64_t end) {
    size_t slen = (size_t)(end - start + 1);
    char* buf = malloc(slen + 1);
    if (!buf) return NULL;
    memcpy(buf, str + start - 1, slen);
    buf[slen] = '\0';
    Expr* result = expr_new_string(buf);
    free(buf);
    return result;
}

/*
 * builtin_stringtake:
 * StringTake["string", n]       - first n characters
 * StringTake["string", -n]      - last n characters
 * StringTake["string", UpTo[n]] - first n chars, or as many as available
 * StringTake["string", {n}]     - nth character as a string
 * StringTake["string", {m,n}]   - characters m through n
 * StringTake["string", {m,n,s}] - characters m through n in steps of s
 * StringTake[{s1,s2,...}, spec]  - maps over list of strings
 *
 * Returns NULL (unevaluated) for invalid arguments.
 */
Expr* builtin_stringtake(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;

    Expr* arg0 = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    /* StringTake[{s1, s2, ...}, spec] - map over list of strings */
    if (arg0->type == EXPR_FUNCTION &&
        arg0->data.function.head->type == EXPR_SYMBOL &&
        arg0->data.function.head->data.symbol.name == SYM_List) {
        size_t n = arg0->data.function.arg_count;
        Expr** results = malloc(sizeof(Expr*) * n);
        if (!results) return NULL;
        for (size_t i = 0; i < n; i++) {
            Expr** inner_args = malloc(sizeof(Expr*) * 2);
            inner_args[0] = expr_copy(arg0->data.function.args[i]);
            inner_args[1] = expr_copy(spec);
            Expr* call = expr_new_function(expr_new_symbol(SYM_StringTake),
                                           inner_args, 2);
            free(inner_args);
            Expr* r = evaluate(call);
            expr_free(call);
            results[i] = r;
        }
        Expr* result = expr_new_function(expr_new_symbol(SYM_List), results, n);
        free(results);
        return result;
    }

    if (arg0->type != EXPR_STRING) return NULL;

    const char* str = arg0->data.string;
    int64_t len = (int64_t)strlen(str);

    /* StringTake["string", n] - first n chars (n>0) or last |n| chars (n<0) */
    if (spec->type == EXPR_INTEGER) {
        int64_t n = spec->data.integer;
        if (n == 0) return expr_new_string("");
        if (n > 0) {
            if (n > len) return NULL;
            return stringtake_substring(str, 1, n);
        } else {
            int64_t take = -n;
            if (take > len) return NULL;
            return stringtake_substring(str, len - take + 1, len);
        }
    }

    /* StringTake["string", UpTo[n]] */
    int64_t upto_n;
    if (strings_is_upto(spec, &upto_n)) {
        if (upto_n <= 0) return NULL;
        int64_t take = upto_n < len ? upto_n : len;
        if (take == 0) return expr_new_string("");
        return stringtake_substring(str, 1, take);
    }

    /* StringTake["string", {spec}] - list spec */
    if (spec->type == EXPR_FUNCTION &&
        spec->data.function.head->type == EXPR_SYMBOL &&
        spec->data.function.head->data.symbol.name == SYM_List) {
        size_t spec_argc = spec->data.function.arg_count;

        /* StringTake["string", {n}] - single character */
        if (spec_argc == 1) {
            Expr* a1 = spec->data.function.args[0];
            if (a1->type != EXPR_INTEGER) return NULL;
            int64_t k = a1->data.integer;
            if (k < 0) k = len + k + 1;
            if (k < 1 || k > len) return NULL;
            char buf[2] = { str[k - 1], '\0' };
            return expr_new_string(buf);
        }

        /* StringTake["string", {m, n}] - range */
        if (spec_argc == 2) {
            Expr* a1 = spec->data.function.args[0];
            Expr* a2 = spec->data.function.args[1];
            if (a1->type != EXPR_INTEGER || a2->type != EXPR_INTEGER)
                return NULL;
            int64_t m = a1->data.integer;
            int64_t n = a2->data.integer;
            if (m < 0) m = len + m + 1;
            if (n < 0) n = len + n + 1;
            if (m < 1 || m > len || n < 1 || n > len || m > n) return NULL;
            return stringtake_substring(str, m, n);
        }

        /* StringTake["string", {m, n, s}] - range with step */
        if (spec_argc == 3) {
            Expr* a1 = spec->data.function.args[0];
            Expr* a2 = spec->data.function.args[1];
            Expr* a3 = spec->data.function.args[2];
            if (a1->type != EXPR_INTEGER || a2->type != EXPR_INTEGER ||
                a3->type != EXPR_INTEGER)
                return NULL;
            int64_t m = a1->data.integer;
            int64_t n = a2->data.integer;
            int64_t s = a3->data.integer;
            if (s == 0) return NULL;
            if (m < 0) m = len + m + 1;
            if (n < 0) n = len + n + 1;

            /* Calculate count */
            int64_t count = 0;
            if (s > 0) {
                if (m >= 1 && n <= len && m <= n)
                    count = (n - m) / s + 1;
            } else {
                if (m <= len && n >= 1 && m >= n)
                    count = (m - n) / (-s) + 1;
            }
            if (count <= 0) return expr_new_string("");

            char* buf = malloc((size_t)count + 1);
            if (!buf) return NULL;
            int64_t current = m;
            for (int64_t i = 0; i < count; i++) {
                buf[i] = str[current - 1];
                current += s;
            }
            buf[count] = '\0';
            Expr* result = expr_new_string(buf);
            free(buf);
            return result;
        }
    }

    return NULL;
}
