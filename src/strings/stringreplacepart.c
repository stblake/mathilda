/*
 * stringreplacepart.c - StringReplacePart builtin for Mathilda
 *
 * StringReplacePart["string", "snew", part] - replaces character ranges by new
 * strings, with support for multiple ranges and the operator form.
 */

#include "picostrings.h"
#include "eval.h"
#include "sym_names.h"
#include <string.h>
#include <stdio.h>

/*
 * srp_is_list:
 * Returns true if e is a List[...] expression.
 */
static bool srp_is_list(Expr* e) {
    return e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

/*
 * srp_emit_argt:
 * Emit `StringReplacePart::argt: StringReplacePart called with N arguments;
 * 2 or 3 arguments are expected.` to stderr and return NULL so the call is
 * left unevaluated. Matches Mathematica's diagnostic for a wrong argument
 * count.
 */
static Expr* srp_emit_argt(size_t argc) {
    fprintf(stderr,
            "StringReplacePart::argt: StringReplacePart called with %zu "
            "argument%s; 2 or 3 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

/*
 * srp_emit_ovlp:
 * Emit the StringReplacePart::ovlp diagnostic when a position specification
 * overlaps a previously accepted one. The offending copy of the new string is
 * not inserted. The message reports the position as originally given (before
 * negative-index resolution), matching Mathematica.
 */
static void srp_emit_ovlp(int64_t m, int64_t n, const char* repl) {
    fprintf(stderr,
            "StringReplacePart::ovlp: Position {%lld,%lld} overlaps previous "
            "positions; new string %s will not be inserted.\n",
            (long long)m, (long long)n, repl);
}

/*
 * SrpRange: one resolved replacement instruction. m_orig/n_orig retain the
 * position as written (for the ovlp message); start/end are the 1-based
 * inclusive character range after negative-index resolution; repl is a
 * borrowed pointer to the replacement string for this range.
 */
typedef struct {
    int64_t m_orig, n_orig;
    int64_t start, end;
    const char* repl;
} SrpRange;

/*
 * srp_parse_range:
 * Parse a single position specification {m, n} (a List of two integers) and
 * resolve negative indices against a string of length len. Writes m_orig,
 * n_orig, start, and end into *out. Returns false if r is not a two-integer
 * list.
 */
static bool srp_parse_range(Expr* r, int64_t len, SrpRange* out) {
    if (!srp_is_list(r) || r->data.function.arg_count != 2)
        return false;
    Expr* a = r->data.function.args[0];
    Expr* b = r->data.function.args[1];
    if (a->type != EXPR_INTEGER || b->type != EXPR_INTEGER)
        return false;
    int64_t m = a->data.integer, n = b->data.integer;
    out->m_orig = m;
    out->n_orig = n;
    out->start = m < 0 ? len + m + 1 : m;
    out->end   = n < 0 ? len + n + 1 : n;
    return true;
}

/*
 * builtin_stringreplacepart:
 * StringReplacePart["string", "snew", {m, n}]
 *     replaces characters m through n in "string" by "snew".
 * StringReplacePart["string", "snew", {{m1, n1}, {m2, n2}, ...}]
 *     inserts a copy of "snew" at each of the given ranges.
 * StringReplacePart["string", {"snew1", "snew2", ...}, {{m1, n1}, ...}]
 *     replaces the i-th range by the i-th new string; the two lists must be
 *     the same length.
 * StringReplacePart[{s1, s2, ...}, snew, part]
 *     maps over a list of strings.
 * StringReplacePart[new, part][old]
 *     is the operator form, equivalent to StringReplacePart[old, new, part].
 *
 * Position specifications use the form returned by StringPosition: each is a
 * pair {m, n} of first/last character positions, with negative positions
 * counting from the end. All positions refer to the ORIGINAL string, before
 * any replacement is done. Overlapping positions are not allowed: a range that
 * overlaps a previously accepted one triggers StringReplacePart::ovlp and its
 * new string is not inserted. An empty replacement string deletes the selected
 * characters.
 *
 * A call whose argument count is not 2 or 3 emits StringReplacePart::argt and
 * is left unevaluated. A non-string argument, a malformed or out-of-range
 * position, or a new-string/position length mismatch leaves the call
 * unevaluated (NULL) so symbolic arguments flow through the evaluator
 * unchanged.
 */
Expr* builtin_stringreplacepart(Expr* res) {
    if (res->type != EXPR_FUNCTION)
        return NULL;

    size_t argc = res->data.function.arg_count;

    /* Operator form: StringReplacePart[new, part][old] is realised as the pure
     * function Function[StringReplacePart[#1, new, part]], mirroring how the
     * other operator-form builtins (Cases, ...) are implemented. */
    if (argc == 2) {
        Expr* slot_args[1] = { expr_new_integer(1) };
        Expr* slot = expr_new_function(expr_new_symbol(SYM_Slot), slot_args, 1);
        Expr* inner_args[3] = {
            slot,
            expr_copy(res->data.function.args[0]),
            expr_copy(res->data.function.args[1])
        };
        Expr* inner = expr_new_function(expr_new_symbol(SYM_StringReplacePart),
                                        inner_args, 3);
        Expr* func_args[1] = { inner };
        return expr_new_function(expr_new_symbol(SYM_Function), func_args, 1);
    }

    if (argc != 3)
        return srp_emit_argt(argc);

    Expr* arg0 = res->data.function.args[0];
    Expr* snew = res->data.function.args[1];
    Expr* part = res->data.function.args[2];

    /* StringReplacePart[{s1, s2, ...}, snew, part] - map over list of strings */
    if (srp_is_list(arg0)) {
        size_t n = arg0->data.function.arg_count;
        Expr** results = malloc(sizeof(Expr*) * (n ? n : 1));
        if (!results) return NULL;
        for (size_t i = 0; i < n; i++) {
            Expr** inner_args = malloc(sizeof(Expr*) * 3);
            inner_args[0] = expr_copy(arg0->data.function.args[i]);
            inner_args[1] = expr_copy(snew);
            inner_args[2] = expr_copy(part);
            Expr* call = expr_new_function(
                expr_new_symbol(SYM_StringReplacePart), inner_args, 3);
            free(inner_args);
            results[i] = evaluate(call);
            expr_free(call);
        }
        Expr* result = expr_new_function(expr_new_symbol(SYM_List), results, n);
        free(results);
        return result;
    }

    if (arg0->type != EXPR_STRING)
        return NULL;

    /* The position spec must be a list. */
    if (!srp_is_list(part))
        return NULL;

    const char* str = arg0->data.string;
    int64_t len = (int64_t)strlen(str);
    size_t part_n = part->data.function.arg_count;

    /* An empty position list makes no replacements. */
    if (part_n == 0) {
        if (snew->type == EXPR_STRING) return expr_new_string(str);
        if (srp_is_list(snew) && snew->data.function.arg_count == 0)
            return expr_new_string(str);
        return NULL; /* {snew...} length must match: 0 positions, N strings */
    }

    /* Distinguish the flat single-range form {m, n} from the nested list of
     * ranges {{m1, n1}, ...} by whether the first element is itself a list. */
    bool nested = srp_is_list(part->data.function.args[0]);
    size_t nranges = nested ? part_n : 1;
    Expr** range_exprs = nested ? part->data.function.args : &part;

    /* Resolve the replacement string for each range: a single string is
     * broadcast to every range; a list of strings must match the number of
     * ranges one-for-one. */
    const char** repls = malloc(sizeof(const char*) * nranges);
    if (!repls) return NULL;
    if (snew->type == EXPR_STRING) {
        for (size_t i = 0; i < nranges; i++) repls[i] = snew->data.string;
    } else if (srp_is_list(snew)) {
        if (snew->data.function.arg_count != nranges) { free(repls); return NULL; }
        for (size_t i = 0; i < nranges; i++) {
            Expr* s = snew->data.function.args[i];
            if (s->type != EXPR_STRING) { free(repls); return NULL; }
            repls[i] = s->data.string;
        }
    } else {
        free(repls);
        return NULL;
    }

    /* Parse and validate every range against the original string. */
    SrpRange* ranges = malloc(sizeof(SrpRange) * nranges);
    if (!ranges) { free(repls); return NULL; }
    for (size_t i = 0; i < nranges; i++) {
        if (!srp_parse_range(range_exprs[i], len, &ranges[i])) {
            free(repls); free(ranges); return NULL;
        }
        ranges[i].repl = repls[i];
        if (ranges[i].start < 1 || ranges[i].end > len ||
            ranges[i].start > ranges[i].end) {
            free(repls); free(ranges); return NULL;
        }
    }
    free(repls);

    /* Reject overlapping ranges in the order given: the first claimant of a
     * character wins; a later range that touches an already-claimed position
     * triggers StringReplacePart::ovlp and is dropped. */
    bool* covered = calloc((size_t)(len > 0 ? len : 1), sizeof(bool));
    if (!covered) { free(ranges); return NULL; }
    SrpRange* acc = malloc(sizeof(SrpRange) * nranges);
    if (!acc) { free(covered); free(ranges); return NULL; }
    size_t nacc = 0;
    for (size_t i = 0; i < nranges; i++) {
        bool overlap = false;
        for (int64_t p = ranges[i].start; p <= ranges[i].end; p++) {
            if (covered[p - 1]) { overlap = true; break; }
        }
        if (overlap) {
            srp_emit_ovlp(ranges[i].m_orig, ranges[i].n_orig, ranges[i].repl);
            continue;
        }
        for (int64_t p = ranges[i].start; p <= ranges[i].end; p++)
            covered[p - 1] = true;
        acc[nacc++] = ranges[i];
    }
    free(covered);
    free(ranges);

    /* Sort the accepted, non-overlapping ranges by start position so the
     * output can be assembled in a single left-to-right pass. Insertion sort
     * is ample: the number of ranges is small in practice. */
    for (size_t i = 1; i < nacc; i++) {
        SrpRange key = acc[i];
        size_t j = i;
        while (j > 0 && acc[j - 1].start > key.start) {
            acc[j] = acc[j - 1];
            j--;
        }
        acc[j] = key;
    }

    /* First pass: compute the exact output length. */
    size_t total = 0;
    int64_t cursor = 1;
    for (size_t k = 0; k < nacc; k++) {
        if (acc[k].start > cursor) total += (size_t)(acc[k].start - cursor);
        total += strlen(acc[k].repl);
        cursor = acc[k].end + 1;
    }
    if (cursor <= len) total += (size_t)(len - cursor + 1);

    /* Second pass: emit literal spans and replacement strings. */
    char* buf = malloc(total + 1);
    if (!buf) { free(acc); return NULL; }
    size_t off = 0;
    cursor = 1;
    for (size_t k = 0; k < nacc; k++) {
        if (acc[k].start > cursor) {
            size_t span = (size_t)(acc[k].start - cursor);
            memcpy(buf + off, str + cursor - 1, span);
            off += span;
        }
        size_t rl = strlen(acc[k].repl);
        memcpy(buf + off, acc[k].repl, rl);
        off += rl;
        cursor = acc[k].end + 1;
    }
    if (cursor <= len) {
        size_t span = (size_t)(len - cursor + 1);
        memcpy(buf + off, str + cursor - 1, span);
        off += span;
    }
    buf[off] = '\0';

    free(acc);
    Expr* result = expr_new_string(buf);
    free(buf);
    return result;
}
