#ifndef STRINGS_INTERNAL_H
#define STRINGS_INTERNAL_H

/*
 * strings_internal.h - helpers shared between the per-builtin translation
 * units under src/strings/. Kept as static inline functions so no separate
 * object file is needed; each user simply includes this header.
 *
 * The positional builtins (StringTake, StringDrop, StringPart, StringInsert,
 * StringReplacePart, StringPartition) index CHARACTERS (code points), as
 * Mathematica does, not bytes: every one of them resolves character positions
 * through a U8View (utf8.h), which costs nothing for pure-ASCII strings.
 */

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "expr.h"
#include "sym_names.h"
#include "print.h"
#include "message.h"
#include "utf8.h"
#include "linkedlist.h"
#include "eval.h"   /* eval_mark_fixed_point: strings_uneval */

/* True if e is the symbol whose interned name is `name`. */
static inline bool strings_is_sym(const Expr* e, const char* name) {
    return e->type == EXPR_SYMBOL && e->data.symbol.name == name;
}

/* True if e is head[...] for the interned symbol name `head`. */
static inline bool strings_has_head(const Expr* e, const char* head) {
    return e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == head;
}

static inline bool strings_is_list(const Expr* e) {
    return strings_has_head(e, SYM_List);
}

/* True if e is a List whose every element is a String. */
static inline bool strings_all_strings(const Expr* e) {
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (e->data.function.args[i]->type != EXPR_STRING) return false;
    return true;
}

/* A List or a LinkedList: a sequence of subjects the string heads thread over.
 * A LinkedList is read through LLReader (linkedlist.h) -- never converted --
 * and the answer is built in its image with ll_seq_like. */
static inline bool strings_is_seq(const Expr* e) {
    return strings_is_list(e) || e->type == EXPR_LINKEDLIST;
}

/* True if the sequence e (strings_is_seq) holds only Strings. A typed
 * LinkedList holds numbers, so only a boxed (or empty) one can qualify. */
static inline bool strings_seq_all_strings(const Expr* e) {
    const LinkedList* ll = expr_linkedlist(e);
    if (!ll) return strings_all_strings(e);
    int64_t n = ll_length(ll);
    if (n == 0) return true;
    if (ll_kind(ll) != LL_BOXED) return false;
    for (int64_t i = 0; i < n; i++)
        if (ll_peek_boxed(ll, i)->type != EXPR_STRING) return false;
    return true;
}

/*
 * strings_machine_int:
 * True (writing *out) if e is a machine-sized integer. When allow_real is set,
 * an integer-valued machine Real (1., -3.) is accepted too: the sequence
 * specifications of StringTake/StringDrop and the Span bounds of StringPart
 * take them, while position specifications (StringPart n, StringInsert,
 * StringReplacePart) do not. A bignum is never machine-sized.
 */
static inline bool strings_machine_int(const Expr* e, bool allow_real,
                                       int64_t* out) {
    if (e->type == EXPR_INTEGER) { *out = e->data.integer; return true; }
    if (allow_real && e->type == EXPR_REAL) {
        double d = e->data.real;
        /* [-2^63, 2^63) are exactly the doubles that convert without UB. */
        if (d == floor(d) && d >= -9223372036854775808.0 &&
            d < 9223372036854775808.0) {
            *out = (int64_t)d;
            return true;
        }
    }
    return false;
}

/*
 * strings_is_upto:
 * Returns true if e is UpTo[n] where n is an integer, and writes the integer
 * value to *out.
 */
static inline bool strings_is_upto(Expr* e, int64_t* out) {
    if (strings_has_head(e, SYM_UpTo) &&
        e->data.function.arg_count == 1 &&
        e->data.function.args[0]->type == EXPR_INTEGER) {
        *out = e->data.function.args[0]->data.integer;
        return true;
    }
    return false;
}

/* True if e is Infinity (the bare symbol or DirectedInfinity[1]). */
static inline bool strings_is_infinity(const Expr* e) {
    if (strings_is_sym(e, SYM_Infinity)) return true;
    return strings_has_head(e, SYM_DirectedInfinity) &&
           e->data.function.arg_count == 1 &&
           e->data.function.args[0]->type == EXPR_INTEGER &&
           e->data.function.args[0]->data.integer == 1;
}

/* ------------------------------------------------------------ messages */

/*
 * Wrong-argument-count diagnostics, with Mathematica's tag choice: a fixed
 * arity builtin called with one argument says `argr`, any other count `argrx`;
 * a two-or-three-argument builtin says `argtu` for one argument, else `argt`.
 * All return NULL so the call stays unevaluated.
 */
static inline Expr* strings_argr_error(const char* head, size_t argc,
                                       size_t expected) {
    mth_message(head, argc == 1 ? "argr" : "argrx",
                "%s called with %zu argument%s; %zu arguments are expected.",
                head, argc, argc == 1 ? "" : "s", expected);
    return NULL;
}

static inline Expr* strings_argt_error(const char* head, size_t argc) {
    mth_message(head, argc == 1 ? "argtu" : "argt",
                "%s called with %zu argument%s; 2 or 3 arguments are expected.",
                head, argc, argc == 1 ? "" : "s");
    return NULL;
}

/* Head::strse: String or list of strings expected at position pos in call. */
static inline Expr* strings_msg_strse(const char* head, Expr* call, int pos) {
    char* c = expr_to_string_outputform(call);
    mth_message(head, "strse",
                "String or list of strings expected at position %d in %s.",
                pos, c ? c : "");
    free(c);
    return NULL;
}

/* Head::string: String expected at position pos in call. */
static inline Expr* strings_msg_string(const char* head, Expr* call, int pos) {
    char* c = expr_to_string_outputform(call);
    mth_message(head, "string", "String expected at position %d in %s.",
                pos, c ? c : "");
    free(c);
    return NULL;
}

/* ------------------------------------------- standard sequence specs */

/*
 * StrSeq: one parsed standard sequence specification, as taken by StringTake
 * and StringDrop (and, through Span, StringPart). RANGE carries {m, n, s} as
 * written (unresolved, s != 0); UPTO carries the cap in m (INT64_MAX for
 * UpTo[Infinity]).
 */
typedef enum { SSEQ_RANGE, SSEQ_ALL, SSEQ_NONE, SSEQ_UPTO } StrSeqKind;
typedef struct {
    StrSeqKind kind;
    int64_t    m, n, s;
} StrSeq;

/*
 * strings_seq_parse:
 * Parse one sequence specification: n, -n, {n}, {m, n}, {m, n, s} (s != 0),
 * All, None, UpTo[k] (k a non-negative machine integer or Infinity). Bounds
 * must be machine-sized (integer-valued Reals accepted, as Mathematica does).
 * Returns false for anything else, including {m, n, 0}; the caller decides
 * whether that is an error or (StringTake) a list of specifications.
 */
static inline bool strings_seq_parse(Expr* spec, StrSeq* out) {
    int64_t a, b, c;
    if (strings_machine_int(spec, true, &a)) {
        out->kind = SSEQ_RANGE;
        if (a >= 0) { out->m = 1; out->n = a; }
        else        { out->m = a; out->n = -1; }
        out->s = 1;
        return true;
    }
    if (strings_is_sym(spec, SYM_All))  { out->kind = SSEQ_ALL;  return true; }
    if (strings_is_sym(spec, SYM_None)) { out->kind = SSEQ_NONE; return true; }
    if (strings_has_head(spec, SYM_UpTo) && spec->data.function.arg_count == 1) {
        Expr* k = spec->data.function.args[0];
        if (k->type == EXPR_INTEGER && k->data.integer >= 0) {
            out->kind = SSEQ_UPTO; out->m = k->data.integer; return true;
        }
        if (strings_is_infinity(k)) {
            out->kind = SSEQ_UPTO; out->m = INT64_MAX; return true;
        }
        return false;
    }
    /* {n}, {m, n} or {m, n, s}: a List, or a LinkedList read from its slots. */
    const LinkedList* sl = expr_linkedlist(spec);
    if (!sl && !strings_is_list(spec)) return false;
    size_t k = sl ? (size_t)ll_length(sl) : spec->data.function.arg_count;
    if (k < 1 || k > 3) return false;
    Expr* tmp[3] = { NULL, NULL, NULL };
    Expr** v = spec->data.function.args;
    if (sl) {
        for (size_t i = 0; i < k; i++) tmp[i] = ll_get(sl, (int64_t)i);
        v = tmp;
    }
    bool ok = strings_machine_int(v[0], true, &a);
    b = a; c = 1;
    if (ok && k >= 2) ok = strings_machine_int(v[1], true, &b);
    if (ok && k == 3) ok = strings_machine_int(v[2], true, &c) && c != 0;
    for (size_t i = 0; i < 3; i++) if (tmp[i]) expr_free(tmp[i]);
    if (!ok) return false;
    out->kind = SSEQ_RANGE;
    out->m = a; out->n = b; out->s = c;
    return true;
}

/*
 * strings_seq_range:
 * Resolve the range {m, n, s} against a string of `len` characters, using
 * Mathematica's Take rules (verified exhaustively against Mathematica 14 for
 * lengths 0, 3 and 5): negative positions count from the end (k -> len+1+k),
 * a literal m == 0 means "just past the end" (len+1), and the result may be
 * empty only when the bounds are adjacent (n == m-1 for s > 0, n == m+1 for
 * s < 0). Valid positions are 1 <= M <= len+1, 0 <= N <= len for s > 0, and
 * 0 <= M <= len, 1 <= N <= len+1 for s < 0. On success writes the 1-based first
 * position *first and the element count *count (positions first, first+s, ...).
 * Every bound is validated before any arithmetic that could overflow.
 */
static inline bool strings_seq_range(int64_t m, int64_t n, int64_t s,
                                     int64_t len, int64_t* first,
                                     int64_t* count) {
    int64_t M = m < 0 ? len + 1 + m : (m == 0 ? len + 1 : m);
    int64_t N = n < 0 ? len + 1 + n : n;
    if (s > 0) {
        if (M < 1 || M > len + 1 || N < 0 || N > len || N < M - 1) return false;
        *count = N < M ? 0 : (N - M) / s + 1;
    } else {
        if (M < 0 || M > len || N < 1 || N > len + 1 || M < N - 1) return false;
        *count = N > M ? 0 : (N - M) / s + 1;
    }
    *first = M;
    return true;
}

/*
 * strings_seq_resolve:
 * Resolve any StrSeq against `len` characters into (first, step, count).
 * Returns false only for an out-of-range RANGE.
 */
static inline bool strings_seq_resolve(const StrSeq* q, int64_t len,
                                       int64_t* first, int64_t* step,
                                       int64_t* count) {
    *first = 1; *step = 1;
    switch (q->kind) {
    case SSEQ_ALL:  *count = len; return true;
    case SSEQ_NONE: *count = 0;   return true;
    case SSEQ_UPTO: *count = q->m < len ? q->m : len; return true;
    case SSEQ_RANGE:
        *step = q->s;
        return strings_seq_range(q->m, q->n, q->s, len, first, count);
    }
    return false;
}

/*
 * strings_u8_advance:
 * Advance from byte i of s[0, n) over up to k characters; returns the new byte
 * index and writes the number of characters actually passed to *done (< k
 * only when the string ends first). ASCII bytes are stepped inline; only a
 * non-ASCII lead byte pays the u8_char_len call, so the split is exactly
 * utf8.h's.
 */
static inline size_t strings_u8_advance(const char* s, size_t n, size_t i,
                                        int64_t k, int64_t* done) {
    int64_t d = 0;
    while (d < k && i < n) {
        if (!((unsigned char)s[i] & 0x80)) i++;
        else i += u8_char_len(s, n, i);
        d++;
    }
    *done = d;
    return i;
}

/*
 * StrChars: the character geometry of one string WITHOUT an offset table --
 * its byte and character lengths and whether it is pure ASCII. For a builtin
 * that touches one contiguous run of characters (StringTake n / {m, n},
 * StringDrop, StringPart n) this avoids a U8View's O(n) allocation: the byte
 * span of the run costs one O(prefix + run) scan, and nothing for ASCII. A
 * builtin that needs arbitrary positions inside the run builds a U8View over
 * just that byte span (the greedy UTF-8 split of a slice that starts and ends
 * on character boundaries is the same as the whole string's).
 *
 * `prefix`, when >= 0, says the caller only needs the first `prefix`
 * characters to resolve its positions (a forward spec anchored at the start);
 * the length scan then stops there and len is min(prefix, true length) --
 * exact when the string is shorter, and a sufficient stand-in otherwise (see
 * strings_seq_prefix). So StringTake[s, 5] on a long non-ASCII string never
 * scans past character 5.
 */
typedef struct {
    const char* s;
    size_t      nbytes;
    int64_t     len;       /* characters (see `prefix` above) */
    bool        ascii;
} StrChars;

static inline void strings_chars_init(StrChars* c, const char* s, int64_t prefix) {
    c->s = s;
    c->nbytes = strlen(s);
    c->ascii = u8_is_ascii(s, c->nbytes) != 0;
    if (c->ascii) c->len = (int64_t)c->nbytes;
    else if (prefix >= 0) (void)strings_u8_advance(s, c->nbytes, 0, prefix, &c->len);
    else c->len = (int64_t)u8_length(s, c->nbytes);
}

/*
 * strings_seq_prefix:
 * The number of leading characters that decide a sequence spec, or -1 if the
 * whole length is needed. A RANGE with 1 <= m, 0 <= n, s > 0 resolves the same
 * against any len >= n (M <= N+1 <= len+1 and N <= len hold already), and
 * UpTo[k] needs only min(k, len); everything else (negative or zero bounds,
 * negative steps, All) depends on the true length.
 */
static inline int64_t strings_seq_prefix(const StrSeq* q) {
    if (q->kind == SSEQ_UPTO) return q->m;
    if (q->kind == SSEQ_RANGE && q->m >= 1 && q->n >= 0 && q->s > 0) return q->n;
    return -1;
}

/* Byte span [*b0, *b1) of characters lo..hi (1-based inclusive, lo <= hi+1). */
static inline void strings_chars_span(const StrChars* c, int64_t lo, int64_t hi,
                                      size_t* b0, size_t* b1) {
    if (c->ascii) { *b0 = (size_t)(lo - 1); *b1 = (size_t)hi; return; }
    int64_t d;
    *b0 = strings_u8_advance(c->s, c->nbytes, 0, lo - 1, &d);
    *b1 = strings_u8_advance(c->s, c->nbytes, *b0, hi - lo + 1, &d);
}

/* Byte offset of 0-based character k in the viewed string. */
#define STR_OFF(v, k) u8_view_off((v), (size_t)(k))

/* The 1-based character `pos` of the viewed string as a new String. */
static inline Expr* strings_char_at(const U8View* v, int64_t pos) {
    size_t b0 = STR_OFF(v, pos - 1);
    return expr_new_string_n(v->s + b0, STR_OFF(v, pos) - b0);
}

/* head[a, b] / head[a, b, c] from copies, for an unevaluated element: the
 * call the head just declined (after its message) on evaluated arguments.
 * It is stamped as a fixed point, so the evaluator's walk of the threaded
 * result does not re-ask the head and print the message a second time
 * (StringPart[{"abc", "de"}, 3] printed StringPart::partw twice). */
static inline Expr* strings_uneval(const char* head, Expr* a, Expr* b, Expr* c) {
    Expr* args[3] = { expr_copy(a), expr_copy(b), c ? expr_copy(c) : NULL };
    Expr* call = expr_new_function(expr_new_symbol(head), args, c ? 3 : 2);
    eval_mark_fixed_point(call);
    return call;
}

/* A String holding the character s[0, n) (one UTF-8 sequence). Single ASCII
 * characters are shared, refcounted Strings, so this allocates nothing for
 * them. Defined in characters.c. */
Expr* strings_char_expr(const char* s, size_t n);

/* Mathematica's argument diagnostic "head::tag: <what> expected at position
 * <pos> in <call>." (e.g. tag "string", what "String"), routed through
 * mth_message. Returns NULL. Defined in strings_common.c. */
Expr* strings_arg_message(const char* head, const char* tag, const char* what,
                          size_t pos, Expr* call);

/* StringTake / StringDrop of one string by one parsed sequence spec (see
 * stringtake.c / stringdrop.c); NULL after StringTake::take / StringDrop::drop
 * when the positions are out of range. Shared with StringTakeDrop. */
Expr* strings_take_seq(Expr* str, const StrSeq* q);
Expr* strings_drop_seq(Expr* str, const StrSeq* q);

#endif // STRINGS_INTERNAL_H
