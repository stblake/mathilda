/*
 * stringpart.c - StringPart builtin for Mathilda
 *
 * StringPart["string", spec] - extracts characters by index, list, or Span
 */

#include "picostrings.h"
#include "eval.h"
#include "sym_names.h"
#include <string.h>

/*
 * stringpart_single:
 * Helper to extract a single character from a string by 1-based index.
 * Negative indices count from the end. Returns a new single-character
 * string Expr*, or NULL if the index is out of bounds or not an integer.
 */
static Expr* stringpart_single(const char* str, int64_t len, Expr* idx) {
    if (idx->type != EXPR_INTEGER) return NULL;
    int64_t k = idx->data.integer;
    if (k < 0) k = len + k + 1;
    if (k < 1 || k > len) return NULL;
    char buf[2] = { str[k - 1], '\0' };
    return expr_new_string(buf);
}

/*
 * builtin_stringpart:
 * StringPart["string", n]        - gives the nth character as a string
 * StringPart["string", {n1,n2,...}] - gives a list of characters
 * StringPart["string", m;;n]     - gives characters m through n
 * StringPart["string", m;;n;;s]  - gives characters m through n in steps of s
 * StringPart[{s1,s2,...}, spec]   - gives the list of results for each si
 *
 * Negative indices count from the end. Returns NULL (unevaluated) for
 * invalid arguments.
 */
Expr* builtin_stringpart(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;

    Expr* arg0 = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    /* StringPart[{s1, s2, ...}, spec] - map over list of strings */
    if (arg0->type == EXPR_FUNCTION &&
        arg0->data.function.head->type == EXPR_SYMBOL &&
        arg0->data.function.head->data.symbol == SYM_List) {
        size_t n = arg0->data.function.arg_count;
        Expr** results = malloc(sizeof(Expr*) * n);
        if (!results) return NULL;
        for (size_t i = 0; i < n; i++) {
            /* Build StringPart[si, spec] and evaluate */
            Expr** inner_args = malloc(sizeof(Expr*) * 2);
            inner_args[0] = expr_copy(arg0->data.function.args[i]);
            inner_args[1] = expr_copy(spec);
            Expr* call = expr_new_function(expr_new_symbol(SYM_StringPart),
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

    /* StringPart["string", n] - single integer index */
    if (spec->type == EXPR_INTEGER) {
        return stringpart_single(str, len, spec);
    }

    /* StringPart["string", {n1, n2, ...}] - list of indices */
    if (spec->type == EXPR_FUNCTION &&
        spec->data.function.head->type == EXPR_SYMBOL &&
        spec->data.function.head->data.symbol == SYM_List) {
        size_t n = spec->data.function.arg_count;
        Expr** results = malloc(sizeof(Expr*) * n);
        if (!results) return NULL;
        for (size_t i = 0; i < n; i++) {
            Expr* ch = stringpart_single(str, len, spec->data.function.args[i]);
            if (!ch) {
                for (size_t j = 0; j < i; j++) expr_free(results[j]);
                free(results);
                return NULL;
            }
            results[i] = ch;
        }
        Expr* result = expr_new_function(expr_new_symbol(SYM_List), results, n);
        free(results);
        return result;
    }

    /* StringPart["string", m;;n] or StringPart["string", m;;n;;s] - Span */
    if (spec->type == EXPR_FUNCTION &&
        spec->data.function.head->type == EXPR_SYMBOL &&
        spec->data.function.head->data.symbol == SYM_Span) {
        int64_t start = 1, end = len, step = 1;
        size_t span_argc = spec->data.function.arg_count;

        if (span_argc >= 1) {
            Expr* a1 = spec->data.function.args[0];
            if (a1->type == EXPR_INTEGER) {
                start = a1->data.integer;
                if (start < 0) start = len + start + 1;
            } else if (a1->type == EXPR_SYMBOL &&
                       a1->data.symbol == SYM_All) {
                start = 1;
            } else return NULL;
        }
        if (span_argc >= 2) {
            Expr* a2 = spec->data.function.args[1];
            if (a2->type == EXPR_INTEGER) {
                end = a2->data.integer;
                if (end < 0) end = len + end + 1;
            } else if (a2->type == EXPR_SYMBOL &&
                       a2->data.symbol == SYM_All) {
                end = len;
            } else return NULL;
        }
        if (span_argc >= 3) {
            Expr* a3 = spec->data.function.args[2];
            if (a3->type == EXPR_INTEGER) step = a3->data.integer;
            else return NULL;
            if (step == 0) return NULL;
        }

        /* Calculate number of elements */
        int64_t count = 0;
        if (step > 0) {
            if (start >= 1 && end <= len && start <= end)
                count = (end - start) / step + 1;
        } else {
            if (start <= len && end >= 1 && start >= end)
                count = (start - end) / (-step) + 1;
        }
        if (count < 0) count = 0;

        Expr** results = NULL;
        if (count > 0) {
            results = malloc(sizeof(Expr*) * (size_t)count);
            if (!results) return NULL;
        }

        int64_t current = start;
        for (int64_t i = 0; i < count; i++) {
            char buf[2] = { str[current - 1], '\0' };
            results[i] = expr_new_string(buf);
            current += step;
        }

        Expr* result = expr_new_function(expr_new_symbol(SYM_List),
                                         results, (size_t)count);
        free(results);
        return result;
    }

    return NULL;
}
