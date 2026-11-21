/*
 * stringdrop.c - StringDrop builtin for Mathilda
 *
 * StringDrop["string", spec] - removes characters by count, range, or step.
 * Positions count characters (code points), not UTF-8 bytes.
 */

#include "picostrings.h"
#include "strings_internal.h"
#include "eval.h"
#include "sym_names.h"
#include <string.h>
#include <inttypes.h>

/*
 * sd_emit_drop:
 * StringDrop::drop: Cannot drop positions m through n in "string".
 */
static void sd_emit_drop(const StrSeq* q, Expr* str) {
    char* s = expr_to_string(str);
    mth_message("StringDrop", "drop",
                "Cannot drop positions %" PRId64 " through %" PRId64 " in %s.",
                q->m, q->n, s ? s : "");
    free(s);
}

/* StringDrop::seqs - the spec is not a sequence specification. Unlike
 * StringTake, StringDrop takes no list of specifications. */
static Expr* sd_emit_seqs(Expr* call) {
    char* c = expr_to_string_outputform(call);
    mth_message("StringDrop", "seqs",
                "Sequence specification (+n, -n, {+n}, {-n}, {m, n}, or {m, n, s}) "
                "expected at position 2 in %s.", c ? c : "");
    free(c);
    return NULL;
}

/*
 * drop_scalar:
 * StringDrop[str, q]. The dropped positions are first, first+step, ...
 * (count of them), spanning characters lo..hi; everything else is kept in
 * order. The prefix before lo and the suffix after hi are copied whole; a
 * contiguous run (|step| == 1, or one position) drops the span outright, a
 * strided one walks a U8View of just the span, testing membership
 * arithmetically (no per-character mask). Returns NULL after StringDrop::drop
 * when the positions are out of range.
 */
Expr* strings_drop_seq(Expr* str, const StrSeq* q) {
    if (q->kind == SSEQ_ALL) return expr_new_string_n("", 0); /* no scan needed */
    if (q->kind == SSEQ_NONE) return expr_copy(str);
    StrChars c;
    strings_chars_init(&c, str->data.string, strings_seq_prefix(q));

    int64_t first, step, count;
    if (!strings_seq_resolve(q, c.len, &first, &step, &count)) {
        sd_emit_drop(q, str);
        return NULL;
    }
    if (count == 0) return expr_copy(str);

    int64_t last = first + (count - 1) * step;
    int64_t lo = first < last ? first : last;
    int64_t hi = first < last ? last : first;
    size_t b_lo, b_hi;
    strings_chars_span(&c, lo, hi, &b_lo, &b_hi);

    char* buf = malloc(c.nbytes + 1);
    if (!buf) return NULL;
    memcpy(buf, c.s, b_lo);
    size_t o = b_lo;
    if (count >= 2 && step != 1 && step != -1) {
        U8View v;                                /* characters lo..hi only */
        if (u8_view_init(&v, c.s + b_lo, b_hi - b_lo) != 0) { free(buf); return NULL; }
        /* count >= 2 bounds |step| by len, so negating it is safe. */
        int64_t a = step > 0 ? step : -step;
        for (int64_t k = 0; k < (int64_t)v.nchars; k++) {
            if (k % a == 0) continue;            /* lo, lo+a, ..., hi dropped */
            size_t b0 = STR_OFF(&v, k), b1 = STR_OFF(&v, k + 1);
            memcpy(buf + o, v.s + b0, b1 - b0);
            o += b1 - b0;
        }
        u8_view_free(&v);
    }
    memcpy(buf + o, c.s + b_hi, c.nbytes - b_hi);
    o += c.nbytes - b_hi;
    buf[o] = '\0';
    return expr_new_string_adopt(buf);
}

/*
 * builtin_stringdrop:
 * StringDrop["string", n]       - drop the first n characters
 * StringDrop["string", -n]      - drop the last n characters
 * StringDrop["string", {n}]     - drop the nth character
 * StringDrop["string", {m, n}]  - drop characters m through n
 * StringDrop["string", {m, n, s}] - drop characters m through n in steps of s
 * StringDrop["string", UpTo[n]] - drop n characters, or as many as available
 * StringDrop["string", All]     - drop everything; None drops nothing
 * StringDrop[{s1, s2, ...}, spec] - maps over the list of strings
 *
 * StringDrop is the complement of StringTake: the same standard sequence
 * specification selects the characters to remove, with the same range rules
 * (strings_seq_range). Messages (call left unevaluated): argr/argrx, strse for
 * a first argument that is not a string or list of strings, seqs for an
 * invalid spec, drop for out-of-range positions. Over a list of strings, an
 * element whose positions are out of range stays as the unevaluated
 * StringDrop[si, spec] inside the result.
 */
Expr* builtin_stringdrop(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return strings_argr_error("StringDrop", argc, 2);

    Expr* arg0 = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    bool threaded = strings_is_seq(arg0);
    if (!(arg0->type == EXPR_STRING || (threaded && strings_seq_all_strings(arg0))))
        return strings_msg_strse("StringDrop", res, 1);

    StrSeq q;
    if (strings_is_seq(spec) && ll_seq_length(spec) == 0)
        q.kind = SSEQ_NONE;                       /* {} drops nothing */
    else if (!strings_seq_parse(spec, &q))
        return sd_emit_seqs(res);

    if (!threaded) return strings_drop_seq(arg0, &q);

    LLElems sv;                           /* a List or a LinkedList of strings */
    ll_elems_open(arg0, &sv);
    size_t n = sv.n;
    Expr** results = malloc(sizeof(Expr*) * (n ? n : 1));
    if (!results) { ll_elems_close(&sv); return NULL; }
    for (size_t i = 0; i < n; i++) {
        Expr* si = sv.v[i];
        results[i] = strings_drop_seq(si, &q);
        if (!results[i]) results[i] = strings_uneval(SYM_StringDrop, si, spec, NULL);
    }
    ll_elems_close(&sv);
    Expr* out = ll_seq_like(arg0, results, n);
    free(results);
    return out;
}
