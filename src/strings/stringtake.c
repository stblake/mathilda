/*
 * stringtake.c - StringTake builtin for Mathilda
 *
 * StringTake["string", spec] - extracts substrings by count, range, or step.
 * Positions count characters (code points), not UTF-8 bytes.
 */

#include "picostrings.h"
#include "strings_internal.h"
#include "eval.h"
#include "sym_names.h"
#include <string.h>
#include <inttypes.h>

/*
 * st_emit_take:
 * StringTake::take: Cannot take positions m through n in "string".
 * Positions are reported as written ({m, n} of the spec); the string in
 * InputForm, quoted, exactly as Mathematica does.
 */
static void st_emit_take(const StrSeq* q, Expr* str) {
    char* s = expr_to_string(str);
    mth_message("StringTake", "take",
                "Cannot take positions %" PRId64 " through %" PRId64 " in %s.",
                q->m, q->n, s ? s : "");
    free(s);
}

/* StringTake::mseqs - the spec is not a sequence specification at all. */
static Expr* st_emit_mseqs(Expr* call) {
    char* c = expr_to_string_outputform(call);
    mth_message("StringTake", "mseqs",
                "Sequence specification (+n, -n, {+n}, {-n}, {m, n}, or {m, n, s}) "
                "or a list of sequence specifications expected at position 2 in %s.",
                c ? c : "");
    free(c);
    return NULL;
}

/*
 * take_scalar:
 * StringTake[str, q] for one parsed sequence spec. Returns the new String, or
 * NULL after emitting StringTake::take when the positions are out of range
 * (or on allocation failure). Only the byte span of the selected characters
 * is ever scanned: a forward contiguous range is one memcpy of it, any other
 * step gathers characters through a U8View of just that span.
 */
Expr* strings_take_seq(Expr* str, const StrSeq* q) {
    if (q->kind == SSEQ_ALL) return expr_copy(str);          /* no scan needed */
    if (q->kind == SSEQ_NONE) return expr_new_string_n("", 0);
    StrChars c;
    strings_chars_init(&c, str->data.string, strings_seq_prefix(q));

    int64_t first, step, count;
    if (!strings_seq_resolve(q, c.len, &first, &step, &count)) {
        st_emit_take(q, str);
        return NULL;
    }
    if (count == 0) return expr_new_string_n("", 0);

    int64_t last = first + (count - 1) * step;   /* a valid position */
    int64_t lo = step > 0 ? first : last, hi = step > 0 ? last : first;
    size_t b_lo, b_hi;
    strings_chars_span(&c, lo, hi, &b_lo, &b_hi);
    if (step == 1) return expr_new_string_n(c.s + b_lo, b_hi - b_lo);

    U8View v;                                    /* characters lo..hi only */
    if (u8_view_init(&v, c.s + b_lo, b_hi - b_lo) != 0) return NULL;
    /* Positions are distinct, so the result is at most the span. */
    char* buf = malloc(v.nbytes + 1);
    if (!buf) { u8_view_free(&v); return NULL; }
    size_t o = 0;
    for (int64_t i = 0; i < count; i++) {
        int64_t k = first + i * step - lo;       /* 0-based within the span */
        size_t b0 = STR_OFF(&v, k), b1 = STR_OFF(&v, k + 1);
        memcpy(buf + o, v.s + b0, b1 - b0);
        o += b1 - b0;
    }
    buf[o] = '\0';
    u8_view_free(&v);
    return expr_new_string_adopt(buf);
}

/*
 * StSpec: the whole second argument, parsed once (not per string). Either a
 * single sequence spec, or a list of them (`list` set, one StrSeq per element
 * of `elems`).
 */
typedef struct {
    bool    list;
    StrSeq  one;
    StrSeq* seqs;
    Expr**  elems;
    size_t  n;
} StSpec;

/*
 * st_parse_spec:
 * Parse StringTake's second argument. A List of 1-3 integers is a single
 * spec; any other List is a list of specs, each of which must itself be a
 * sequence spec. Mathematica's two quirks are kept: {m, n, 0} silently becomes
 * the list of specs {m, n, 0}, and a list of four or more integers does too,
 * after a StringTake::ambgsntx warning. {} takes nothing. Returns false (the
 * caller emits mseqs) for anything unparseable.
 */
static bool st_parse_spec(Expr* spec, StSpec* out) {
    out->list = false; out->seqs = NULL; out->elems = NULL; out->n = 0;
    if (strings_seq_parse(spec, &out->one)) return true;
    if (!strings_is_list(spec)) return false;

    size_t k = spec->data.function.arg_count;
    Expr** v = spec->data.function.args;
    if (k == 0) { out->one.kind = SSEQ_NONE; return true; }

    bool all_ints = true;
    int64_t tmp;
    for (size_t i = 0; i < k && all_ints; i++)
        all_ints = strings_machine_int(v[i], true, &tmp);
    if (all_ints && k >= 4)
        mth_message("StringTake", "ambgsntx",
                    "Warning: interpreting list of integers as a list of sequence specifications.");

    out->seqs = malloc(sizeof(StrSeq) * k);
    if (!out->seqs) return false;
    for (size_t i = 0; i < k; i++) {
        if (strings_is_list(v[i]) && v[i]->data.function.arg_count == 0) {
            out->seqs[i].kind = SSEQ_NONE;          /* {} takes nothing */
            continue;
        }
        if (!strings_seq_parse(v[i], &out->seqs[i])) {
            free(out->seqs);
            out->seqs = NULL;
            return false;
        }
    }
    out->list = true;
    out->elems = v;
    out->n = k;
    return true;
}

/*
 * take_one_string:
 * Apply a parsed spec to one string. For a single spec a failure returns NULL
 * when `top` (the call stays unevaluated), else the unevaluated element
 * StringTake[str, spec]; a list of specs always yields a List whose failed
 * entries are unevaluated StringTake[str, spec_i] (Mathematica's behaviour).
 */
static Expr* take_one_string(Expr* str, Expr* spec, const StSpec* p, bool top) {
    if (!p->list) {
        Expr* r = strings_take_seq(str, &p->one);
        if (!r && !top) r = strings_uneval(SYM_StringTake, str, spec, NULL);
        return r;
    }
    Expr** items = malloc(sizeof(Expr*) * p->n);
    if (!items) return NULL;
    for (size_t i = 0; i < p->n; i++) {
        items[i] = strings_take_seq(str, &p->seqs[i]);
        if (!items[i])
            items[i] = strings_uneval(SYM_StringTake, str, p->elems[i], NULL);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, p->n);
    free(items);
    return out;
}

/*
 * builtin_stringtake:
 * StringTake["string", n]       - first n characters
 * StringTake["string", -n]      - last n characters
 * StringTake["string", {n}]     - nth character as a string
 * StringTake["string", {m,n}]   - characters m through n
 * StringTake["string", {m,n,s}] - characters m through n in steps of s
 * StringTake["string", UpTo[n]] - first n chars, or as many as available
 * StringTake["string", All]     - the whole string; None gives ""
 * StringTake["string", {spec1, spec2, ...}] - a list of substrings
 * StringTake[{s1,s2,...}, spec] - maps over a list of strings
 *
 * Messages (call left unevaluated): argr/argrx for the wrong argument count,
 * strse for a first argument that is not a string or list of strings, mseqs
 * for an invalid spec, take for out-of-range positions. Over a list of
 * strings, an element whose positions are out of range is left as the
 * unevaluated StringTake[si, spec] inside the result list.
 */
Expr* builtin_stringtake(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return strings_argr_error("StringTake", argc, 2);

    Expr* arg0 = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    bool threaded = strings_is_seq(arg0);
    if (!(arg0->type == EXPR_STRING || (threaded && strings_seq_all_strings(arg0))))
        return strings_msg_strse("StringTake", res, 1);

    StSpec p;
    if (!st_parse_spec(spec, &p)) return st_emit_mseqs(res);

    Expr* out;
    if (!threaded) {
        out = take_one_string(arg0, spec, &p, true);
    } else {
        LLElems sv;                       /* a List or a LinkedList of strings */
        ll_elems_open(arg0, &sv);
        size_t n = sv.n;
        Expr** results = malloc(sizeof(Expr*) * (n ? n : 1));
        if (!results) { ll_elems_close(&sv); free(p.seqs); return NULL; }
        for (size_t i = 0; i < n; i++)
            results[i] = take_one_string(sv.v[i], spec, &p, false);
        ll_elems_close(&sv);
        out = ll_seq_like(arg0, results, n);
        free(results);
    }
    free(p.seqs);
    return out;
}
