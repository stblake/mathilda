/*
 * stringreplacepart.c - StringReplacePart builtin for Mathilda
 *
 * StringReplacePart["string", "snew", part] - replaces character ranges by new
 * strings, with support for multiple ranges and the operator form. Positions
 * count characters (code points), not UTF-8 bytes.
 */

#include "picostrings.h"
#include "strings_internal.h"
#include "eval.h"
#include "sym_names.h"
#include <string.h>
#include <inttypes.h>

/*
 * SrpRange: one replacement instruction. m_orig/n_orig keep the position as
 * written (for messages); start/end are the resolved 1-based inclusive
 * character range, with end == start - 1 for an empty range (a pure
 * insertion before character `start`); repl is the replacement String
 * (borrowed); order is the index in the position list.
 */
typedef struct {
    int64_t m_orig, n_orig;
    int64_t start, end;
    Expr*   repl;
    size_t  order;
} SrpRange;

/* StringReplacePart::psl for the scalar spec `what` inside `where`. */
static Expr* srp_emit_psl(Expr* what, Expr* where) {
    char* a = expr_to_string_outputform(what);
    char* b = expr_to_string_outputform(where);
    mth_message("StringReplacePart", "psl",
                "Position specification %s in %s is not a machine-sized integer "
                "or a list of machine-sized integers.", a ? a : "", b ? b : "");
    free(a);
    free(b);
    return NULL;
}

/* StringReplacePart::spos: element `index` of `list` is not a pair. */
static Expr* srp_emit_spos(size_t index, Expr* list) {
    char* b = expr_to_string_outputform(list);
    mth_message("StringReplacePart", "spos",
                "Position specification {m, n} is expected at position %zu in %s.",
                index, b ? b : "");
    free(b);
    return NULL;
}

/* StringReplacePart::plen: the new strings and the positions differ in length. */
static Expr* srp_emit_plen(Expr* snew, Expr* list) {
    char* a = expr_to_string_outputform(snew);
    char* b = expr_to_string_outputform(list);
    mth_message("StringReplacePart", "plen", "%s and %s should have the same length.",
                a ? a : "", b ? b : "");
    free(a);
    free(b);
    return NULL;
}

/*
 * srp_cmp_sorted:
 * Mathematica orders the accepted ranges as the triples {start, end, "new"}
 * in canonical order, so an insertion sorts before a replacement starting at
 * the same character, and equal positions fall back to the new strings.
 */
static int srp_cmp_sorted(const void* x, const void* y) {
    const SrpRange* a = (const SrpRange*)x;
    const SrpRange* b = (const SrpRange*)y;
    if (a->start != b->start) return a->start < b->start ? -1 : 1;
    if (a->end != b->end) return a->end < b->end ? -1 : 1;
    int c = expr_compare(a->repl, b->repl);
    if (c) return c;
    return (a->order > b->order) - (a->order < b->order);
}

/*
 * srp_overlaps:
 * Phase 1 overlap test of the nonempty range [s, e] against the accepted
 * nonempty ranges, kept disjoint and sorted by start in acc[0..nacc).
 * Binary search for the first range ending at or after s.
 */
static bool srp_overlaps(const SrpRange* acc, size_t nacc, int64_t s, int64_t e,
                         size_t* slot) {
    size_t lo = 0, hi = nacc;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (acc[mid].end < s) lo = mid + 1; else hi = mid;
    }
    *slot = lo;
    return lo < nacc && acc[lo].start <= e;
}

/*
 * srp_apply:
 * Validate the parsed ranges against `str` and build the result. Order of
 * work, each step matching Mathematica 14:
 *   1. every range must be valid ({m, n} under the Take rules, step 1),
 *      else StringReplacePart::repart and the call stays unevaluated;
 *   2. nonempty ranges are checked for overlap in the order GIVEN: the first
 *      claimant of a character wins, a later one gets ::ovlp and is dropped;
 *   3. the survivors are sorted canonically (srp_cmp_sorted) and walked: an
 *      insertion that falls inside the preceding replacement (start <= its
 *      end) is dropped silently and that replacement keeps its last character
 *      (its end moves back by one). This reproduces Mathematica's output
 *      exactly, quirk included;
 *   4. one pass emits the untouched spans and the new strings.
 */
static Expr* srp_apply(Expr* str, SrpRange* r, size_t nr) {
    const char* s = str->data.string;
    U8View v;
    if (u8_view_init(&v, s, strlen(s)) != 0) return NULL;
    int64_t len = (int64_t)v.nchars;

    for (size_t i = 0; i < nr; i++) {
        int64_t first, count;
        if (!strings_seq_range(r[i].m_orig, r[i].n_orig, 1, len, &first, &count)) {
            char* q = expr_to_string(str);
            mth_message("StringReplacePart", "repart",
                        "Cannot replace positions %" PRId64 " through %" PRId64 " in %s.",
                        r[i].m_orig, r[i].n_orig, q ? q : "");
            free(q);
            u8_view_free(&v);
            return NULL;
        }
        r[i].start = first;
        r[i].end = first + count - 1;
    }

    /* Phase 1: overlap among nonempty ranges, in the order given. kept[]
     * collects every survivor (insertions included); acc[] the nonempty
     * survivors sorted by start, for the binary search. */
    SrpRange* kept = malloc(sizeof(SrpRange) * (nr ? nr : 1));
    SrpRange* acc = malloc(sizeof(SrpRange) * (nr ? nr : 1));
    if (!kept || !acc) { free(kept); free(acc); u8_view_free(&v); return NULL; }
    size_t nk = 0, nacc = 0;
    for (size_t i = 0; i < nr; i++) {
        if (r[i].end >= r[i].start) {
            size_t slot;
            if (srp_overlaps(acc, nacc, r[i].start, r[i].end, &slot)) {
                char* q = expr_to_string_outputform(r[i].repl);
                mth_message("StringReplacePart", "ovlp",
                            "Position {%" PRId64 ", %" PRId64 "} overlaps previous "
                            "positions; new string %s will not be inserted.",
                            r[i].m_orig, r[i].n_orig, q ? q : "");
                free(q);
                continue;
            }
            memmove(acc + slot + 1, acc + slot, sizeof(SrpRange) * (nacc - slot));
            acc[slot] = r[i];
            nacc++;
        }
        kept[nk++] = r[i];
    }
    free(acc);

    /* Phase 2: canonical order; drop insertions inside a replacement. */
    qsort(kept, nk, sizeof(SrpRange), srp_cmp_sorted);
    size_t nout = 0;
    for (size_t i = 0; i < nk; i++) {
        bool insertion = kept[i].end < kept[i].start;
        if (insertion && nout > 0 && kept[i].start <= kept[nout - 1].end) {
            kept[nout - 1].end--;
            continue;
        }
        kept[nout++] = kept[i];
    }

    /* Output length, then one emitting pass. */
    size_t total = 0;
    int64_t cursor = 1;
    for (size_t k = 0; k < nout; k++) {
        total += STR_OFF(&v, kept[k].start - 1) - STR_OFF(&v, cursor - 1);
        total += strlen(kept[k].repl->data.string);
        cursor = kept[k].end + 1;
    }
    total += v.nbytes - STR_OFF(&v, cursor - 1);

    char* buf = malloc(total + 1);
    if (!buf) { free(kept); u8_view_free(&v); return NULL; }
    size_t o = 0;
    cursor = 1;
    for (size_t k = 0; k < nout; k++) {
        size_t b0 = STR_OFF(&v, cursor - 1), b1 = STR_OFF(&v, kept[k].start - 1);
        memcpy(buf + o, s + b0, b1 - b0);
        o += b1 - b0;
        size_t rl = strlen(kept[k].repl->data.string);
        memcpy(buf + o, kept[k].repl->data.string, rl);
        o += rl;
        cursor = kept[k].end + 1;
    }
    size_t bt = STR_OFF(&v, cursor - 1);
    memcpy(buf + o, s + bt, v.nbytes - bt);
    o += v.nbytes - bt;
    buf[o] = '\0';

    free(kept);
    u8_view_free(&v);
    return expr_new_string_adopt(buf);
}

/*
 * builtin_stringreplacepart:
 * StringReplacePart["string", "snew", {m, n}]
 *     replaces characters m through n in "string" by "snew"; {m, m-1}
 *     inserts "snew" before character m.
 * StringReplacePart["string", "snew", n]
 *     replaces the first n characters (-n: the last n), as StringTake's n.
 * StringReplacePart["string", "snew", {{m1, n1}, {m2, n2}, ...}]
 *     inserts a copy of "snew" at each of the given ranges.
 * StringReplacePart["string", {"snew1", "snew2", ...}, {{m1, n1}, ...}]
 *     replaces the i-th range by the i-th new string.
 * StringReplacePart[new, part][old]
 *     is the operator form, equivalent to StringReplacePart[old, new, part].
 *
 * All positions refer to the ORIGINAL string. Unlike StringTake, the first
 * argument does not thread over a list (Mathematica emits ::string).
 * Messages (call left unevaluated unless noted): argtu/argt for the wrong
 * argument count, string (position 1), strse (position 2), plen, spos, psl,
 * repart; ovlp drops only the overlapping replacement.
 */
Expr* builtin_stringreplacepart(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    /* Operator form: StringReplacePart[new, part][old] is realised as the pure
     * function Function[StringReplacePart[#1, new, part]], mirroring how the
     * other operator-form builtins (Cases, ...) are implemented. */
    /* StringReplacePart[new, part] is the operator form: it stays inert
     * (Head StringReplacePart, not a pure Function) and is applied by the
     * generic operator-form table, src/opform.c:
     * StringReplacePart[new, part][old] == StringReplacePart[old, new, part]. */
    if (argc == 2) return NULL;
    if (argc != 3) return strings_argt_error("StringReplacePart", argc);

    Expr* arg0 = res->data.function.args[0];
    Expr* snew = res->data.function.args[1];
    Expr* part = res->data.function.args[2];

    if (arg0->type != EXPR_STRING)
        return strings_msg_string("StringReplacePart", res, 1);
    bool snew_list = strings_is_list(snew);
    if (!(snew->type == EXPR_STRING || (snew_list && strings_all_strings(snew))))
        return strings_msg_strse("StringReplacePart", res, 2);

    /* Normalise the positions to a list of ranges. A bare integer n is the
     * range of StringTake's n; a flat {m, n} is one range; a list whose first
     * element is a list is a list of ranges. */
    int64_t one_m = 0, one_n = 0;
    Expr** ranges;
    size_t nr;
    Expr* wrapped = NULL;           /* {part}, for messages about a flat range */
    if (part->type == EXPR_INTEGER) {
        int64_t k = part->data.integer;
        if (k >= 0) { one_m = 1; one_n = k; } else { one_m = k; one_n = -1; }
        ranges = NULL;
        nr = 1;
    } else if (strings_is_list(part)) {
        if (part->data.function.arg_count == 0) return expr_copy(arg0);
        if (strings_is_list(part->data.function.args[0])) {
            ranges = part->data.function.args;
            nr = part->data.function.arg_count;
        } else {
            ranges = &part;
            nr = 1;
            Expr* w[1] = { expr_copy(part) };
            wrapped = expr_new_function(expr_new_symbol(SYM_List), w, 1);
        }
    } else {
        return srp_emit_psl(part, res);
    }
    Expr* plist = wrapped ? wrapped : part;

    Expr* out = NULL;
    SrpRange* r = NULL;
    if (snew_list && snew->data.function.arg_count != nr) {
        srp_emit_plen(snew, plist);
        goto done;
    }
    r = malloc(sizeof(SrpRange) * nr);
    if (!r) goto done;
    for (size_t i = 0; i < nr; i++) {
        r[i].repl = snew_list ? snew->data.function.args[i] : snew;
        r[i].order = i;
        if (!ranges) { r[i].m_orig = one_m; r[i].n_orig = one_n; continue; }
        Expr* e = ranges[i];
        if (!strings_is_list(e) || e->data.function.arg_count != 2) {
            srp_emit_spos(i + 1, plist);
            goto done;
        }
        Expr* a = e->data.function.args[0];
        Expr* b = e->data.function.args[1];
        if (a->type != EXPR_INTEGER || b->type != EXPR_INTEGER) {
            srp_emit_psl(e, plist);
            goto done;
        }
        r[i].m_orig = a->data.integer;
        r[i].n_orig = b->data.integer;
    }
    out = srp_apply(arg0, r, nr);

done:
    free(r);
    if (wrapped) expr_free(wrapped);
    return out;
}
