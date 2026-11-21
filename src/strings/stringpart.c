/*
 * stringpart.c - StringPart builtin for Mathilda
 *
 * StringPart["string", spec] - extracts characters by index, list, or Span.
 * Positions count characters (code points), not UTF-8 bytes.
 */

#include "picostrings.h"
#include "strings_internal.h"
#include "eval.h"
#include "sym_names.h"
#include <string.h>

/* StringPart::partw: Part spec of string does not exist. (OutputForm both) */
static void sp_emit_partw(Expr* spec, Expr* str) {
    char* a = expr_to_string_outputform(spec);
    mth_message("StringPart", "partw", "Part %s of %s does not exist.",
                a ? a : "", str->data.string);
    free(a);
}

/* StringPart::take: Cannot take positions m through n in string. The Span
 * bounds are printed as written (they may be All, a symbol, a bignum). */
static void sp_emit_take(Expr* span, Expr* str) {
    char* a = expr_to_string_outputform(span->data.function.args[0]);
    char* b = expr_to_string_outputform(span->data.function.args[1]);
    mth_message("StringPart", "take", "Cannot take positions %s through %s in %s.",
                a ? a : "", b ? b : "", str->data.string);
    free(a);
    free(b);
}

/* StringPart::pkspec1: The expression spec cannot be used as a part spec. */
static Expr* sp_emit_pkspec1(Expr* spec) {
    char* a = expr_to_string_outputform(spec);
    mth_message("StringPart", "pkspec1",
                "The expression %s cannot be used as a part specification.",
                a ? a : "");
    free(a);
    return NULL;
}

/* Resolve a 1-based / negative index against len; false if out of range. */
static bool sp_index(int64_t k, int64_t len, int64_t* pos) {
    if (k < 0) k = len + 1 + k;          /* k >= INT64_MIN + 1 + len: no overflow */
    if (k < 1 || k > len) return false;
    *pos = k;
    return true;
}

/*
 * SpKind: the spec, classified once (it does not depend on the string).
 * A Span's bounds are parsed here too; its range is resolved per string.
 */
typedef enum { SPK_INVALID, SPK_INDEX, SPK_LIST, SPK_ALL, SPK_SPAN, SPK_SPAN_BAD } SpKind;

static SpKind sp_classify(Expr* spec, int64_t* m, int64_t* n, int64_t* s) {
    if (spec->type == EXPR_INTEGER || spec->type == EXPR_BIGINT) return SPK_INDEX;
    if (strings_is_seq(spec)) return SPK_LIST;   /* a List or a LinkedList */
    if (strings_is_sym(spec, SYM_All)) return SPK_ALL;
    if (strings_has_head(spec, SYM_Span) &&
        (spec->data.function.arg_count == 2 || spec->data.function.arg_count == 3)) {
        Expr** a = spec->data.function.args;
        /* All as the start is 1; as the end it is -1 (the last character). */
        if (strings_is_sym(a[0], SYM_All)) *m = 1;
        else if (!strings_machine_int(a[0], true, m)) return SPK_SPAN_BAD;
        if (strings_is_sym(a[1], SYM_All)) *n = -1;
        else if (!strings_machine_int(a[1], true, n)) return SPK_SPAN_BAD;
        *s = 1;
        if (spec->data.function.arg_count == 3 &&
            (!strings_machine_int(a[2], true, s) || *s == 0))
            return SPK_SPAN_BAD;
        return SPK_SPAN;
    }
    return SPK_INVALID;
}

/*
 * part_scalar:
 * StringPart[str, spec] for a classified spec. Returns the character (n) or
 * list of characters, or NULL after emitting partw / take. A single index or
 * a Span only scans the bytes up to the characters it selects (a U8View of
 * just the Span's byte range); an index list views the whole string.
 */
static Expr* part_scalar(Expr* str, Expr* spec, SpKind kind,
                         int64_t m, int64_t n, int64_t s) {
    /* A positive single index needs only the characters up to it. */
    StrChars c;
    int64_t prefix = -1;
    if (kind == SPK_INDEX && spec->type == EXPR_INTEGER && spec->data.integer > 0)
        prefix = spec->data.integer;
    strings_chars_init(&c, str->data.string, prefix);
    int64_t len = c.len;
    size_t b0, b1;

    if (kind == SPK_INDEX) {
        int64_t pos;
        if (spec->type != EXPR_INTEGER || !sp_index(spec->data.integer, len, &pos)) {
            sp_emit_partw(spec, str);
            return NULL;
        }
        strings_chars_span(&c, pos, pos, &b0, &b1);
        return expr_new_string_n(c.s + b0, b1 - b0);
    }

    if (kind == SPK_LIST) {
        /* Every index must be a machine integer in range, else the whole
         * list is reported. Validate first, then view once. */
        LLElems ie;                       /* a List or a LinkedList of indices */
        ll_elems_open(spec, &ie);
        size_t k = ie.n;
        Expr** idx = ie.v;
        int64_t pos;
        for (size_t i = 0; i < k; i++) {
            if (idx[i]->type != EXPR_INTEGER || !sp_index(idx[i]->data.integer, len, &pos)) {
                ll_elems_close(&ie);
                sp_emit_partw(spec, str);
                return NULL;
            }
        }
        U8View v;
        if (u8_view_init(&v, c.s, c.nbytes) != 0) { ll_elems_close(&ie); return NULL; }
        Expr** items = malloc(sizeof(Expr*) * (k ? k : 1));
        if (!items) { ll_elems_close(&ie); u8_view_free(&v); return NULL; }
        for (size_t i = 0; i < k; i++) {
            sp_index(idx[i]->data.integer, len, &pos);
            items[i] = strings_char_at(&v, pos);
        }
        ll_elems_close(&ie);
        Expr* out = ll_seq_like(spec, items, k);
        free(items);
        u8_view_free(&v);
        return out;
    }

    /* All, or a Span under the Take rules. */
    int64_t first = 1, count = len;
    if (kind == SPK_ALL) s = 1;
    if (kind == SPK_SPAN_BAD ||
        (kind == SPK_SPAN && !strings_seq_range(m, n, s, len, &first, &count))) {
        sp_emit_take(spec, str);
        return NULL;
    }
    Expr** items = malloc(sizeof(Expr*) * (count ? (size_t)count : 1));
    if (!items) return NULL;
    if (count > 0) {
        int64_t last = first + (count - 1) * s;
        int64_t lo = s > 0 ? first : last, hi = s > 0 ? last : first;
        strings_chars_span(&c, lo, hi, &b0, &b1);
        U8View v;                                /* characters lo..hi only */
        if (u8_view_init(&v, c.s + b0, b1 - b0) != 0) { free(items); return NULL; }
        for (int64_t i = 0; i < count; i++)
            items[i] = strings_char_at(&v, first + i * s - lo + 1);
        u8_view_free(&v);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)count);
    free(items);
    return out;
}

/*
 * builtin_stringpart:
 * StringPart["string", n]        - gives the nth character as a string
 * StringPart["string", {n1,n2,...}] - gives a list of characters
 * StringPart["string", m;;n]     - gives characters m through n
 * StringPart["string", m;;n;;s]  - gives characters m through n in steps of s
 * StringPart["string", All]      - gives every character
 * StringPart[{s1,s2,...}, spec]   - gives the list of results for each si
 *
 * Negative indices count from the end; Span ranges follow the Take rules
 * (strings_seq_range). Messages (call left unevaluated): argr/argrx, strse,
 * partw for a nonexistent index, take for an invalid Span, pkspec1 for a spec
 * that is none of the above. Over a list of strings a failing element stays
 * as the unevaluated StringPart[si, spec].
 */
Expr* builtin_stringpart(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return strings_argr_error("StringPart", argc, 2);

    Expr* arg0 = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    bool threaded = strings_is_seq(arg0);
    if (!(arg0->type == EXPR_STRING || (threaded && strings_seq_all_strings(arg0))))
        return strings_msg_strse("StringPart", res, 1);

    int64_t m = 1, n = -1, s = 1;
    SpKind kind = sp_classify(spec, &m, &n, &s);
    if (kind == SPK_INVALID) return sp_emit_pkspec1(spec);

    if (!threaded) return part_scalar(arg0, spec, kind, m, n, s);

    LLElems sv;                           /* a List or a LinkedList of strings */
    ll_elems_open(arg0, &sv);
    size_t cnt = sv.n;
    Expr** results = malloc(sizeof(Expr*) * (cnt ? cnt : 1));
    if (!results) { ll_elems_close(&sv); return NULL; }
    for (size_t i = 0; i < cnt; i++) {
        Expr* si = sv.v[i];
        results[i] = part_scalar(si, spec, kind, m, n, s);
        if (!results[i]) results[i] = strings_uneval(SYM_StringPart, si, spec, NULL);
    }
    ll_elems_close(&sv);
    Expr* out = ll_seq_like(arg0, results, cnt);
    free(results);
    return out;
}
