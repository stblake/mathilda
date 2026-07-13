/*
 * stringinsert.c - StringInsert builtin for Mathilda
 *
 * StringInsert["string", "snew", spec] - inserts a string at given positions
 */

#include "picostrings.h"
#include "eval.h"
#include "sym_names.h"
#include <string.h>
#include <stdio.h>

/*
 * si_emit_argrx:
 * Emit `StringInsert::argrx: StringInsert called with N arguments; 3 arguments
 * are expected.` to stderr and return NULL so the call is left unevaluated.
 * Matches Mathematica's diagnostic for a wrong argument count.
 */
static Expr* si_emit_argrx(size_t argc) {
    fprintf(stderr,
            "StringInsert::argrx: StringInsert called with %zu argument%s; "
            "3 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

/*
 * si_pos_to_offset:
 * Convert one StringInsert position specification p (which must be a nonzero
 * integer) into a 0-based insertion offset in [0, len] measured against the
 * ORIGINAL string. StringInsert["string", "snew", n] makes the first character
 * of "snew" the nth character of the result, i.e. it inserts before original
 * character n, so a positive n maps to offset n - 1. A negative -k makes the
 * last character of "snew" the kth character from the end, i.e. it inserts
 * before original character len - k + 2, so a negative n maps to offset
 * len + n + 1. Returns true and writes *out on success, false if p is not an
 * integer, is zero, or resolves out of range.
 */
static bool si_pos_to_offset(Expr* p, int64_t len, int64_t* out) {
    if (p->type != EXPR_INTEGER) return false;
    int64_t n = p->data.integer;
    int64_t off;
    if (n > 0) off = n - 1;
    else if (n < 0) off = len + n + 1;
    else return false;              /* position 0 is invalid */
    if (off < 0 || off > len) return false;
    *out = off;
    return true;
}

/*
 * builtin_stringinsert:
 * StringInsert["string", "snew", n]           - insert "snew" before position n
 * StringInsert["string", "snew", -n]          - insert n positions from the end
 * StringInsert["string", "snew", {n1, n2, ...}] - insert a copy at each position
 * StringInsert[{s1, s2, ...}, "snew", spec]   - maps over the list of strings
 *
 * All positions refer to the ORIGINAL string, before any insertion is done, so
 * the offsets are resolved up front and every copy of "snew" is placed relative
 * to the untouched input. Negative positions count from the end (see
 * si_pos_to_offset). A call whose argument count is not three emits
 * StringInsert::argrx and is left unevaluated; a non-string first/second
 * argument, a non-integer position, or an out-of-range position leaves the call
 * unevaluated (NULL) so symbolic arguments flow through the evaluator unchanged.
 */
Expr* builtin_stringinsert(Expr* res) {
    if (res->type != EXPR_FUNCTION)
        return NULL;

    size_t argc = res->data.function.arg_count;
    if (argc != 3)
        return si_emit_argrx(argc);

    Expr* arg0 = res->data.function.args[0];
    Expr* snew = res->data.function.args[1];
    Expr* spec = res->data.function.args[2];

    /* StringInsert[{s1, s2, ...}, "snew", spec] - map over list of strings */
    if (arg0->type == EXPR_FUNCTION &&
        arg0->data.function.head->type == EXPR_SYMBOL &&
        arg0->data.function.head->data.symbol.name == SYM_List) {
        size_t n = arg0->data.function.arg_count;
        Expr** results = malloc(sizeof(Expr*) * n);
        if (!results) return NULL;
        for (size_t i = 0; i < n; i++) {
            Expr** inner_args = malloc(sizeof(Expr*) * 3);
            inner_args[0] = expr_copy(arg0->data.function.args[i]);
            inner_args[1] = expr_copy(snew);
            inner_args[2] = expr_copy(spec);
            Expr* call = expr_new_function(expr_new_symbol(SYM_StringInsert),
                                           inner_args, 3);
            free(inner_args);
            Expr* r = evaluate(call);
            expr_free(call);
            results[i] = r;
        }
        Expr* result = expr_new_function(expr_new_symbol(SYM_List), results, n);
        free(results);
        return result;
    }

    if (arg0->type != EXPR_STRING || snew->type != EXPR_STRING)
        return NULL;

    const char* str = arg0->data.string;
    const char* ins = snew->data.string;
    int64_t len = (int64_t)strlen(str);
    size_t inslen = strlen(ins);

    /* Gather the position specs into a flat array of Expr* to resolve. */
    Expr** pos_specs;
    size_t npos;
    if (spec->type == EXPR_INTEGER) {
        npos = 1;
        pos_specs = &spec;
    } else if (spec->type == EXPR_FUNCTION &&
               spec->data.function.head->type == EXPR_SYMBOL &&
               spec->data.function.head->data.symbol.name == SYM_List) {
        npos = spec->data.function.arg_count;
        pos_specs = spec->data.function.args;
    } else {
        return NULL;
    }

    /* counts[i] = number of copies of "snew" to insert before original
     * character i (0-based); index len is "insert at the very end". */
    int64_t* counts = calloc((size_t)len + 1, sizeof(int64_t));
    if (!counts) return NULL;

    for (size_t i = 0; i < npos; i++) {
        int64_t off;
        if (!si_pos_to_offset(pos_specs[i], len, &off)) {
            free(counts);
            return NULL;
        }
        counts[off]++;
    }

    /* Build the result: total length is the original plus one "snew" per
     * insertion. npos is bounded and inslen small, so this cannot overflow in
     * any realistic input. */
    size_t total = (size_t)len + npos * inslen;
    char* buf = malloc(total + 1);
    if (!buf) { free(counts); return NULL; }

    size_t pos = 0;
    for (int64_t i = 0; i <= len; i++) {
        for (int64_t c = 0; c < counts[i]; c++) {
            memcpy(buf + pos, ins, inslen);
            pos += inslen;
        }
        if (i < len) buf[pos++] = str[i];
    }
    buf[pos] = '\0';

    free(counts);
    Expr* result = expr_new_string(buf);
    free(buf);
    return result;
}
