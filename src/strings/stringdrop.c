/*
 * stringdrop.c - StringDrop builtin for Mathilda
 *
 * StringDrop["string", spec] - removes characters by count, range, or step
 */

#include "picostrings.h"
#include "strings_internal.h"
#include "eval.h"
#include "sym_names.h"
#include <string.h>
#include <stdio.h>

/*
 * sd_emit_argrx:
 * Emit `StringDrop::argrx: StringDrop called with N arguments; 2 arguments
 * are expected.` to stderr and return NULL so the call is left unevaluated.
 * Matches Mathematica's diagnostic for a wrong argument count.
 */
static Expr* sd_emit_argrx(size_t argc) {
    fprintf(stderr,
            "StringDrop::argrx: StringDrop called with %zu argument%s; "
            "2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

/*
 * stringdrop_build_kept:
 * Rebuild a string from str keeping only the positions i (0-based) for which
 * keep[i] is true, in order. Returns a new EXPR_STRING (possibly ""). The
 * complementary "keep" mask is how StringDrop is expressed: the sequence spec
 * marks characters to drop, and everything else survives.
 */
static Expr* stringdrop_build_kept(const char* str, int64_t len, const bool* keep) {
    char* buf = malloc((size_t)len + 1);
    if (!buf) return NULL;
    int64_t j = 0;
    for (int64_t i = 0; i < len; i++) {
        if (keep[i]) buf[j++] = str[i];
    }
    buf[j] = '\0';
    Expr* result = expr_new_string(buf);
    free(buf);
    return result;
}

/*
 * builtin_stringdrop:
 * StringDrop["string", n]       - drop the first n characters
 * StringDrop["string", -n]      - drop the last n characters
 * StringDrop["string", UpTo[n]] - drop n characters, or as many as available
 * StringDrop["string", {n}]     - drop the nth character
 * StringDrop["string", {m, n}]  - drop characters m through n
 * StringDrop["string", {m, n, s}] - drop characters m through n in steps of s
 * StringDrop[{s1, s2, ...}, spec] - maps over the list of strings
 *
 * StringDrop is the complement of StringTake: the same standard Wolfram
 * Language sequence specification selects the characters to remove, and the
 * remaining characters are concatenated in order. Negative indices count from
 * the end. A call whose argument count is not two emits StringDrop::argrx and
 * is left unevaluated; an out-of-range position or a symbolic argument leaves
 * the call unevaluated (NULL) so it flows through the evaluator unchanged.
 */
Expr* builtin_stringdrop(Expr* res) {
    if (res->type != EXPR_FUNCTION)
        return NULL;

    size_t argc = res->data.function.arg_count;
    if (argc != 2)
        return sd_emit_argrx(argc);

    Expr* arg0 = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    /* StringDrop[{s1, s2, ...}, spec] - map over list of strings */
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
            Expr* call = expr_new_function(expr_new_symbol(SYM_StringDrop),
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

    /* keep[i] tracks whether position i (0-based) survives the drop. */
    bool* keep = malloc(sizeof(bool) * (len > 0 ? (size_t)len : 1));
    if (!keep) return NULL;
    for (int64_t i = 0; i < len; i++) keep[i] = true;

    /* StringDrop["string", n] - drop first n (n>0) or last |n| (n<0) chars */
    if (spec->type == EXPR_INTEGER) {
        int64_t n = spec->data.integer;
        if (n == 0) {
            /* drop nothing */
        } else if (n > 0) {
            if (n > len) { free(keep); return NULL; }
            for (int64_t i = 0; i < n; i++) keep[i] = false;
        } else {
            int64_t take = -n;
            if (take > len) { free(keep); return NULL; }
            for (int64_t i = len - take; i < len; i++) keep[i] = false;
        }
        Expr* result = stringdrop_build_kept(str, len, keep);
        free(keep);
        return result;
    }

    /* StringDrop["string", UpTo[n]] - drop n, or as many as are available */
    int64_t upto_n;
    if (strings_is_upto(spec, &upto_n)) {
        if (upto_n <= 0) { free(keep); return NULL; }
        int64_t take = upto_n < len ? upto_n : len;
        for (int64_t i = 0; i < take; i++) keep[i] = false;
        Expr* result = stringdrop_build_kept(str, len, keep);
        free(keep);
        return result;
    }

    /* StringDrop["string", {spec}] - list spec */
    if (spec->type == EXPR_FUNCTION &&
        spec->data.function.head->type == EXPR_SYMBOL &&
        spec->data.function.head->data.symbol.name == SYM_List) {
        size_t spec_argc = spec->data.function.arg_count;

        /* StringDrop["string", {n}] - drop the single character n */
        if (spec_argc == 1) {
            Expr* a1 = spec->data.function.args[0];
            if (a1->type != EXPR_INTEGER) { free(keep); return NULL; }
            int64_t k = a1->data.integer;
            if (k < 0) k = len + k + 1;
            if (k < 1 || k > len) { free(keep); return NULL; }
            keep[k - 1] = false;
            Expr* result = stringdrop_build_kept(str, len, keep);
            free(keep);
            return result;
        }

        /* StringDrop["string", {m, n}] - drop the range m through n */
        if (spec_argc == 2) {
            Expr* a1 = spec->data.function.args[0];
            Expr* a2 = spec->data.function.args[1];
            if (a1->type != EXPR_INTEGER || a2->type != EXPR_INTEGER) {
                free(keep); return NULL;
            }
            int64_t m = a1->data.integer;
            int64_t n = a2->data.integer;
            if (m < 0) m = len + m + 1;
            if (n < 0) n = len + n + 1;
            /* A decreasing range (m > n) is empty and drops nothing. */
            if (m <= n) {
                if (m < 1 || m > len || n < 1 || n > len) {
                    free(keep); return NULL;
                }
                for (int64_t i = m; i <= n; i++) keep[i - 1] = false;
            }
            Expr* result = stringdrop_build_kept(str, len, keep);
            free(keep);
            return result;
        }

        /* StringDrop["string", {m, n, s}] - drop the range m..n in steps of s */
        if (spec_argc == 3) {
            Expr* a1 = spec->data.function.args[0];
            Expr* a2 = spec->data.function.args[1];
            Expr* a3 = spec->data.function.args[2];
            if (a1->type != EXPR_INTEGER || a2->type != EXPR_INTEGER ||
                a3->type != EXPR_INTEGER) {
                free(keep); return NULL;
            }
            int64_t m = a1->data.integer;
            int64_t n = a2->data.integer;
            int64_t s = a3->data.integer;
            if (s == 0) { free(keep); return NULL; }
            if (m < 0) m = len + m + 1;
            if (n < 0) n = len + n + 1;

            for (int64_t p = m; s > 0 ? p <= n : p >= n; p += s) {
                if (p < 1 || p > len) { free(keep); return NULL; }
                keep[p - 1] = false;
            }
            Expr* result = stringdrop_build_kept(str, len, keep);
            free(keep);
            return result;
        }
    }

    free(keep);
    return NULL;
}
