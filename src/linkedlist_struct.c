/*
 * linkedlist_struct.c -- the copying heads over LinkedList storage.
 *
 * Take, Drop, Most, Rest, Reverse, Append, Prepend, Join, Insert, Delete,
 * Select and Listable threading. Each returns a NEW LinkedList and leaves its
 * argument alone (only the *To forms, Push/Pop/PopFirst and Part assignment
 * mutate). Map is in funcprog.c and ReplacePart in replace.c, next to the
 * level and position machinery they share with the List versions.
 *
 * Two techniques, and no List is ever built of the whole argument:
 *
 *   Structural heads copy RANGES of the storage. A run of typed slots moves a
 *   chunk span at a time (ll_append_range: memcpy on the back run), so
 *   Take/Drop/Most/Rest/Join of a typed list cost one pass of memcpy and never
 *   create an Expr per element. Drop and Delete copy the ranges between the
 *   removed positions; Insert copies the ranges between insertion points.
 *
 *   Element-computing heads (Select, Listable threading, and Map in
 *   funcprog.c) STREAM the list through the existing List machinery in blocks
 *   of LL_STREAM_BLOCK elements: ll_block hands a typed block over as a packed
 *   array gathered straight from the slots, so Sin[ll], ll + 1, Select[ll,
 *   EvenQ] run the very packed kernels and compiled predicates a packed List
 *   would, and agree with them to the bit; each result block is appended
 *   natively (memcpy when it comes back packed). See LL_STREAM_BLOCK for the
 *   block size.
 *
 * Messages: an out-of-range position is reported with the head's own tag
 * (Take::take, Drop::drop, Insert::ins, Delete::partw, Most::nomost,
 * Rest::norest, Join::heads, Thread::tdlen) and the call stays unevaluated.
 */
#include "linkedlist.h"

#include "common.h"        /* head_is */
#include "eval.h"          /* evaluate */
#include "message.h"
#include "ndarray.h"       /* is_ndarray, is_packed_list */
#include "pack.h"
#include "part.h"          /* builtin_insert, builtin_delete */
#include "print.h"         /* expr_to_string */
#include "sym_names.h"
#include "list/join.h"
#include "list/reverse.h"
#include "list/take_drop.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static Expr* ll_wrap(LinkedList* ll) { return expr_new_linkedlist(ll); }

static Expr* sym_call(const char* head, Expr** args, size_t n) {
    return expr_new_function(expr_new_symbol(head), args, n);
}

static bool is_plain_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && head_is(e, SYM_List);
}

/* "LinkedList[<n>]"-style rendering for a message argument. Caller frees. */
static char* shown(const Expr* e) { return expr_to_string((Expr*)e); }

/* ------------------------------------------------------------------------
 *  Take / Drop sequence specifications
 * ---------------------------------------------------------------------- */

/* A resolved sequence specification: `count` positions starting at the
 * 0-based index `first` and stepping by `step`. */
typedef struct { int64_t first, count, step; } LLSeq;

static void seq_message(const char* head, bool take, const Expr* res, int pos,
                        const Expr* ll, int64_t m, int64_t n) {
    (void)res; (void)pos;
    char* ls = shown(ll);
    mth_message(head, take ? "take" : "drop", "Cannot %s positions %lld through %lld in %s.",
                take ? "take" : "drop", (long long)m, (long long)n, ls ? ls : "");
    free(ls);
}

static void seqs_message(const char* head, const Expr* res, int pos) {
    char* rs = shown(res);
    mth_message(head, "seqs",
                "Sequence specification (+n, -n, {+n}, {-n}, {m, n}, or {m, n, s}) "
                "expected at position %d in %s.", pos, rs ? rs : "");
    free(rs);
}

/* Resolve `spec` against a list of length `len` the way Take and Drop read it:
 * n, -n, {m}, {m, n}, {m, n, s}, All, None, UpTo[n]. False after the message
 * (out of range, or not a specification at all). */
static bool ll_seq_resolve(const char* head, bool take, const Expr* res, int pos,
                           const Expr* llx, const Expr* spec, int64_t len, LLSeq* out) {
    int64_t m, n, s = 1;            /* 1-based, inclusive, as written */
    int64_t shown_m, shown_n;       /* the values the message names */
    if (spec->type == EXPR_SYMBOL && spec->data.symbol.name == SYM_All) {
        out->first = 0; out->count = len; out->step = 1;
        return true;
    }
    if (spec->type == EXPR_SYMBOL && spec->data.symbol.name == SYM_None) {
        out->first = 0; out->count = 0; out->step = 1;
        return true;
    }
    if (spec->type == EXPR_INTEGER) {
        int64_t k = spec->data.integer;
        if (k >= 0) { m = 1; n = k; shown_m = 1; shown_n = k; }
        else        { m = len + k + 1; n = len; shown_m = k; shown_n = -1; }
    } else if (spec->type == EXPR_FUNCTION && head_is(spec, SYM_UpTo) &&
               spec->data.function.arg_count == 1 &&
               spec->data.function.args[0]->type == EXPR_INTEGER) {
        int64_t k = spec->data.function.args[0]->data.integer;
        if (k >= 0) { m = 1; n = k > len ? len : k; }
        else        { m = len + k + 1; if (m < 1) m = 1; n = len; }
        shown_m = m; shown_n = n;
    } else if (is_plain_list(spec) && spec->data.function.arg_count >= 1 &&
               spec->data.function.arg_count <= 3) {
        size_t c = spec->data.function.arg_count;
        for (size_t i = 0; i < c; i++)
            if (spec->data.function.args[i]->type != EXPR_INTEGER) {
                seqs_message(head, res, pos);
                return false;
            }
        shown_m = spec->data.function.args[0]->data.integer;
        shown_n = c >= 2 ? spec->data.function.args[1]->data.integer : shown_m;
        if (c == 3) s = spec->data.function.args[2]->data.integer;
        if (s == 0) { seqs_message(head, res, pos); return false; }
        m = shown_m < 0 ? len + shown_m + 1 : shown_m;
        n = shown_n < 0 ? len + shown_n + 1 : shown_n;
    } else {
        seqs_message(head, res, pos);
        return false;
    }
    int64_t count = 0;
    if (s > 0 && m <= n) count = (n - m) / s + 1;
    if (s < 0 && m >= n) count = (m - n) / (-s) + 1;
    if (count > 0) {
        int64_t last = m + (count - 1) * s;
        if (m < 1 || m > len || last < 1 || last > len) {
            seq_message(head, take, res, pos, llx, shown_m, shown_n);
            return false;
        }
    } else if (s > 0 && (m < 1 || m > len + 1 || n > len)) {
        /* An empty forward range must still sit inside the list ({len+1, len}
         * is the empty tail; {1, 7} of a 3-list is out of range). */
        seq_message(head, take, res, pos, llx, shown_m, shown_n);
        return false;
    }
    out->first = m - 1;
    out->count = count;
    out->step = s;
    return true;
}

/* Apply the remaining Take/Drop specifications to one element. NULL when the
 * element cannot take them (the whole call then stays unevaluated). */
static Expr* take_drop_inner(bool take, Expr* elem, Expr** specs, size_t nspecs) {
    Expr** a = malloc(sizeof(Expr*) * (nspecs + 1));
    a[0] = elem;
    for (size_t i = 0; i < nspecs; i++) a[i + 1] = expr_copy(specs[i]);
    Expr* call = sym_call(take ? SYM_Take : SYM_Drop, a, nspecs + 1);
    free(a);
    Expr* r = take ? builtin_take(call) : builtin_drop(call);
    expr_free(call);
    return r;
}

/* Push element i of src to dst, with the deeper specifications applied. */
static bool push_inner(LinkedList* dst, const LinkedList* src, int64_t i, bool take,
                       Expr** specs, size_t nspecs) {
    if (nspecs == 0) { ll_push_back_from(dst, src, i); return true; }
    Expr* r = take_drop_inner(take, ll_get(src, i), specs, nspecs);
    if (!r) return false;
    ll_push_back(dst, r);
    return true;
}

Expr* ll_take_drop(Expr* res, bool take) {
    const char* head = take ? "Take" : "Drop";
    Expr* llx = res->data.function.args[0];
    const LinkedList* src = llx->data.llist;
    int64_t len = ll_length(src);
    Expr** rest = res->data.function.args + 2;
    size_t nrest = res->data.function.arg_count - 2;
    LLSeq q;
    if (!ll_seq_resolve(head, take, res, 2, llx, res->data.function.args[1], len, &q))
        return NULL;
    LinkedList* out = ll_new();
    bool ok = true;
    if (take) {
        if (nrest == 0 && q.step == 1) {
            ll_append_range(out, src, q.first, q.first + q.count);
        } else if (nrest == 0 && q.step == -1) {
            ll_append_range_rev(out, src, q.first - q.count + 1, q.first + 1);
        } else {
            for (int64_t t = 0; t < q.count && ok; t++)
                ok = push_inner(out, src, q.first + t * q.step, true, rest, nrest);
        }
    } else {
        /* Normalise to an ascending run of removed positions. */
        int64_t lo = q.first, step = q.step;
        if (q.count > 0 && step < 0) { lo = q.first + (q.count - 1) * step; step = -step; }
        if (nrest == 0 && (q.count <= 1 || step == 1)) {
            ll_append_range(out, src, 0, lo);
            ll_append_range(out, src, q.count ? lo + q.count * step - step + 1 : lo, len);
        } else {
            int64_t next = q.count ? lo : -1, left = q.count;
            for (int64_t i = 0; i < len && ok; i++) {
                if (i == next) { if (--left > 0) next += step; else next = -1; continue; }
                ok = push_inner(out, src, i, false, rest, nrest);
            }
        }
    }
    if (!ok) { ll_release(out); return NULL; }
    return ll_wrap(out);
}

/* ------------------------------------------------------------------------
 *  Most / Rest / Reverse / Append / Prepend
 * ---------------------------------------------------------------------- */

Expr* ll_most_rest(Expr* res, bool most) {
    Expr* llx = res->data.function.args[0];
    const LinkedList* src = llx->data.llist;
    int64_t len = ll_length(src);
    if (len == 0) {
        char* ls = shown(llx);
        mth_message(most ? "Most" : "Rest", most ? "nomost" : "norest",
                    "Cannot take %s of expression %s with length zero.",
                    most ? "Most" : "Rest", ls ? ls : "");
        free(ls);
        return NULL;
    }
    LinkedList* out = ll_new();
    ll_append_range(out, src, most ? 0 : 1, most ? len - 1 : len);
    return ll_wrap(out);
}

/* Does a Reverse level specification name level `lv`? (An integer or a List
 * of integers, as reverse.c reads it.) */
static bool rev_names_level(const Expr* spec, int64_t lv) {
    if (!spec) return lv == 1;
    if (spec->type == EXPR_INTEGER) return spec->data.integer == lv;
    if (is_plain_list(spec))
        for (size_t i = 0; i < spec->data.function.arg_count; i++) {
            const Expr* a = spec->data.function.args[i];
            if (a->type == EXPR_INTEGER && a->data.integer == lv) return true;
        }
    return false;
}

/* The level specification one level down (every level minus one), or NULL
 * when it names no level below 1 -- the elements are then untouched. */
static Expr* rev_spec_below(const Expr* spec) {
    if (!spec) return NULL;
    if (spec->type == EXPR_INTEGER)
        return spec->data.integer >= 2 ? expr_new_integer(spec->data.integer - 1) : NULL;
    if (!is_plain_list(spec)) return NULL;
    size_t n = spec->data.function.arg_count, k = 0;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        const Expr* x = spec->data.function.args[i];
        if (x->type == EXPR_INTEGER && x->data.integer >= 2)
            a[k++] = expr_new_integer(x->data.integer - 1);
    }
    Expr* r = k ? sym_call(SYM_List, a, k) : NULL;
    free(a);
    return r;
}

Expr* ll_reverse(Expr* res) {
    const LinkedList* src = res->data.function.args[0]->data.llist;
    const Expr* spec = res->data.function.arg_count == 2 ? res->data.function.args[1] : NULL;
    int64_t len = ll_length(src);
    bool top = rev_names_level(spec, 1);
    Expr* below = rev_spec_below(spec);
    LinkedList* out = ll_new();
    if (!below) {
        if (top) ll_append_range_rev(out, src, 0, len);
        else     ll_append_range(out, src, 0, len);
        return ll_wrap(out);
    }
    for (int64_t t = 0; t < len; t++) {
        int64_t i = top ? len - 1 - t : t;
        Expr* a[2] = { ll_get(src, i), expr_copy(below) };
        Expr* call = sym_call(SYM_Reverse, a, 2);
        Expr* r = builtin_reverse(call);
        if (!r) r = expr_copy(call->data.function.args[0]);    /* an atom: unchanged */
        expr_free(call);
        ll_push_back(out, r);
    }
    expr_free(below);
    return ll_wrap(out);
}

Expr* ll_append_prepend(Expr* res, bool prepend) {
    LinkedList* out = ll_copy(res->data.function.args[0]->data.llist);
    Expr* e = expr_copy(res->data.function.args[1]);
    if (prepend) ll_push_front(out, e);
    else         ll_push_back(out, e);
    return ll_wrap(out);
}

/* ------------------------------------------------------------------------
 *  Join
 * ---------------------------------------------------------------------- */

/* The elements of a Join operand as a LinkedList: a LinkedList is borrowed
 * (retained), a List or packed List is copied in (typed buffers by memcpy).
 * NULL for anything else. */
static LinkedList* join_operand(const Expr* a) {
    if (a->type == EXPR_LINKEDLIST) return ll_retain(a->data.llist);
    if (is_plain_list(a) || is_packed_list(a)) {
        LinkedList* t = ll_new();
        ll_append_elements(t, a);
        return t;
    }
    return NULL;
}

static Expr* join_heads_message(const Expr* res, size_t bad) {
    const Expr* a = res->data.function.args[bad];
    Expr* h = a->type == EXPR_FUNCTION ? a->data.function.head : NULL;
    char* hs = h ? shown(h) : shown(a);
    mth_message("Join", "heads",
                "Heads LinkedList and %s at positions 1 and %zu are expected to be the same.",
                hs ? hs : "", bad + 1);
    free(hs);
    return NULL;
}

Expr* ll_join(Expr* res) {
    size_t argc = res->data.function.arg_count, nl = argc;
    int64_t level = 1;
    const Expr* last = res->data.function.args[argc - 1];
    if (argc >= 2 && last->type == EXPR_INTEGER) {
        level = last->data.integer;
        if (level < 1) return NULL;
        nl = argc - 1;
    }
    LinkedList** ops = malloc(sizeof(LinkedList*) * nl);
    for (size_t i = 0; i < nl; i++) {
        ops[i] = join_operand(res->data.function.args[i]);
        if (!ops[i]) {
            for (size_t j = 0; j < i; j++) ll_release(ops[j]);
            free(ops);
            return join_heads_message(res, i);
        }
    }
    LinkedList* out = ll_new();
    bool ok = true;
    if (level == 1) {
        for (size_t i = 0; i < nl; i++) ll_append_range(out, ops[i], 0, ll_length(ops[i]));
    } else {
        /* Join[a, b, n]: position k of the result joins the k-th elements of
         * every operand that has one at level n - 1 (ragged operands simply
         * contribute nothing past their end), as join.c does for Lists. */
        int64_t maxlen = 0;
        for (size_t i = 0; i < nl; i++)
            if (ll_length(ops[i]) > maxlen) maxlen = ll_length(ops[i]);
        Expr** subs = malloc(sizeof(Expr*) * (nl + 1));
        for (int64_t k = 0; k < maxlen && ok; k++) {
            size_t ns = 0;
            for (size_t i = 0; i < nl; i++)
                if (k < ll_length(ops[i])) subs[ns++] = ll_get(ops[i], k);
            if (ns == 1) { ll_push_back(out, subs[0]); continue; }
            subs[ns] = expr_new_integer(level - 1);
            Expr* call = sym_call(SYM_Join, subs, ns + 1);
            Expr* r = builtin_join(call);
            expr_free(call);
            if (!r) { ok = false; break; }
            ll_push_back(out, r);
        }
        free(subs);
    }
    for (size_t i = 0; i < nl; i++) ll_release(ops[i]);
    free(ops);
    if (!ok) { ll_release(out); return NULL; }
    return ll_wrap(out);
}

/* ------------------------------------------------------------------------
 *  Insert / Delete
 * ---------------------------------------------------------------------- */

/* A position specification, split into its first index and the rest of the
 * path. `path` borrows the spec's own nodes. */
typedef struct {
    int64_t idx;         /* 0-based slot (Insert: in [0, len]; Delete: [0, len)) */
    Expr**  rest;        /* the deeper indices, borrowed */
    size_t  nrest;
} LLPos;

static int lpos_cmp(const void* a, const void* b) {
    const LLPos* x = a;
    const LLPos* y = b;
    return x->idx < y->idx ? -1 : x->idx > y->idx ? 1 : 0;
}

static bool is_path(const Expr* p) {
    if (!is_plain_list(p)) return false;
    for (size_t i = 0; i < p->data.function.arg_count; i++)
        if (p->data.function.args[i]->type != EXPR_INTEGER) return false;
    return true;
}

/* Split a position spec into its paths: n, {i, j, ...} (one path) or
 * {{...}, {...}} (several). *paths borrows the spec's nodes; a bare integer is
 * returned through `one`. False for anything that is not a position spec. */
static bool split_paths(Expr* spec, Expr** one, Expr*** paths, size_t* np) {
    if (spec->type == EXPR_INTEGER) { *one = spec; *paths = one; *np = 1; return true; }
    if (!is_plain_list(spec)) return false;
    size_t n = spec->data.function.arg_count;
    if (n == 0) { *paths = NULL; *np = 0; return true; }        /* no positions */
    if (is_path(spec)) { *one = spec; *paths = one; *np = 1; return true; }
    for (size_t i = 0; i < n; i++)
        if (!is_path(spec->data.function.args[i])) return false;
    *paths = spec->data.function.args;
    *np = n;
    return true;
}

static void pos_message(bool insert, const Expr* llx, const Expr* path) {
    char* ps = shown(path);
    char* ls = shown(llx);
    if (insert)
        mth_message("Insert", "ins", "Cannot insert at position %s in %s.", ps ? ps : "", ls ? ls : "");
    else
        mth_message("Delete", "partw", "Part %s of %s does not exist.", ps ? ps : "", ls ? ls : "");
    free(ps);
    free(ls);
}

/* Resolve every path's first index. Insert reads a single-index path as a
 * slot in [1, len + 1] (negative from the end: -1 is after the last
 * element); a deeper path, and every Delete path, names an element. */
static LLPos* resolve_paths(bool insert, const Expr* llx, Expr** paths, size_t np,
                            int64_t len) {
    LLPos* v = malloc(sizeof(LLPos) * (np ? np : 1));
    for (size_t i = 0; i < np; i++) {
        Expr* p = paths[i];
        int64_t k;
        Expr** rest = NULL;
        size_t nrest = 0;
        if (p->type == EXPR_INTEGER) {
            k = p->data.integer;
        } else {
            k = p->data.function.args[0]->data.integer;
            rest = p->data.function.args + 1;
            nrest = p->data.function.arg_count - 1;
        }
        bool slot = insert && nrest == 0;
        int64_t top = slot ? len + 1 : len;
        int64_t j = k > 0 ? k : top + k + 1;
        if (k == 0 || j < 1 || j > top) {
            Expr* shownp = p->type == EXPR_INTEGER
                ? sym_call(SYM_List, (Expr*[]){ expr_copy(p) }, 1) : expr_copy(p);
            pos_message(insert, llx, shownp);
            expr_free(shownp);
            free(v);
            return NULL;
        }
        v[i].idx = j - 1;
        v[i].rest = rest;
        v[i].nrest = nrest;
    }
    qsort(v, np, sizeof(LLPos), lpos_cmp);
    return v;
}

/* Apply the deeper paths v[lo, hi) (all naming element `i`) to that element
 * through the head itself, so a nested LinkedList or List takes them its own
 * way. NULL when the element rejects one (after its message). */
static Expr* apply_deep(bool insert, Expr* elem, Expr* ins, const LLPos* v, size_t lo, size_t hi) {
    size_t m = hi - lo;
    Expr** ps = malloc(sizeof(Expr*) * m);
    for (size_t t = 0; t < m; t++) {
        Expr** a = malloc(sizeof(Expr*) * v[lo + t].nrest);
        for (size_t r = 0; r < v[lo + t].nrest; r++) a[r] = expr_copy(v[lo + t].rest[r]);
        ps[t] = sym_call(SYM_List, a, v[lo + t].nrest);
        free(a);
    }
    Expr* spec = m == 1 ? ps[0] : sym_call(SYM_List, ps, m);
    free(ps);
    Expr* call;
    if (insert) {
        Expr* a[3] = { elem, expr_copy(ins), spec };
        call = sym_call("Insert", a, 3);
    } else {
        Expr* a[2] = { elem, spec };
        call = sym_call("Delete", a, 2);
    }
    Expr* r = insert ? builtin_insert(call) : builtin_delete(call);
    expr_free(call);
    return r;
}

/* Insert[ll, e, n], Insert[ll, e, {i, j, ...}], Insert[ll, e, {{...}, ...}]
 * and Delete[ll, n | {i, ...} | {{...}, ...}], built range by range. */
static Expr* ll_insert_delete(Expr* res, bool insert) {
    Expr* llx = res->data.function.args[0];
    const LinkedList* src = llx->data.llist;
    Expr* ins = insert ? res->data.function.args[1] : NULL;
    Expr* spec = res->data.function.args[insert ? 2 : 1];
    int64_t len = ll_length(src);
    Expr* one = NULL;
    Expr** paths;
    size_t np;
    if (!split_paths(spec, &one, &paths, &np)) return NULL;
    if (np == 0) return ll_wrap(ll_copy(src));                 /* {}: no position */
    /* Delete[ll, {{}}] removes the list itself; Delete[ll, 0] its head. */
    for (size_t i = 0; i < np && !insert; i++) {
        Expr* p = paths[i];
        if (is_plain_list(p) && p->data.function.arg_count == 0)
            return sym_call(SYM_Sequence, NULL, 0);
        if ((p->type == EXPR_INTEGER && p->data.integer == 0) ||
            (is_plain_list(p) && p->data.function.arg_count == 1 &&
             p->data.function.args[0]->data.integer == 0)) {
            Expr** a = malloc(sizeof(Expr*) * (len ? (size_t)len : 1));
            for (int64_t t = 0; t < len; t++) a[t] = ll_get(src, t);
            Expr* r = sym_call(SYM_Sequence, a, (size_t)len);
            free(a);
            return r;
        }
    }
    if (insert && is_plain_list(spec) && spec->data.function.arg_count == 0)
        return ll_wrap(ll_copy(src));
    LLPos* v = resolve_paths(insert, llx, paths, np, len);
    if (!v) return NULL;

    LinkedList* out = ll_new();
    int64_t done = 0;            /* elements [0, done) are in `out` */
    bool ok = true;
    for (size_t a = 0; a < np && ok; ) {
        int64_t i = v[a].idx;
        size_t b = a;
        while (b < np && v[b].idx == i) b++;
        /* Within the group: top-level entries (insertion slots, deletions)
         * and deeper paths into element i. */
        size_t tops = 0, deep_lo = b;
        for (size_t t = a; t < b; t++) {
            if (v[t].nrest == 0) tops++;
            else if (deep_lo == b) deep_lo = t;
        }
        ll_append_range(out, src, done, i);
        done = i;
        if (insert) {
            for (size_t t = 0; t < tops; t++) ll_push_back(out, expr_copy(ins));
        }
        bool deleted = !insert && tops > 0;
        bool has_deep = false;
        for (size_t t = a; t < b; t++) if (v[t].nrest) { has_deep = true; break; }
        if (deleted) {
            done = i + 1;
        } else if (has_deep) {
            /* Gather this group's deep paths contiguously. */
            LLPos* d = malloc(sizeof(LLPos) * (b - a));
            size_t nd = 0;
            for (size_t t = a; t < b; t++) if (v[t].nrest) d[nd++] = v[t];
            Expr* r = apply_deep(insert, ll_get(src, i), ins, d, 0, nd);
            free(d);
            if (!r) { ok = false; break; }
            ll_push_back(out, r);
            done = i + 1;
        }
        (void)deep_lo;
        a = b;
    }
    free(v);
    if (!ok) { ll_release(out); return NULL; }
    ll_append_range(out, src, done, len);
    return ll_wrap(out);
}

Expr* ll_insert(Expr* res) { return ll_insert_delete(res, true); }
Expr* ll_delete(Expr* res) { return ll_insert_delete(res, false); }

/* ------------------------------------------------------------------------
 *  Streaming: Select and Listable threading
 * ---------------------------------------------------------------------- */

/* Append one evaluated result block. False when the result is not a list --
 * an abort or a throw passing through -- which the caller then returns. */
static bool stream_take(LinkedList* out, const Expr* r) {
    if (!ll_is_list_source(r)) return false;
    ll_append_elements(out, r);
    return true;
}

Expr* ll_select(Expr* res) {
    const LinkedList* src = res->data.function.args[0]->data.llist;
    Expr* crit = res->data.function.args[1];
    int64_t want = -1;
    if (res->data.function.arg_count == 3) {
        const Expr* n = res->data.function.args[2];
        if (n->type != EXPR_INTEGER || n->data.integer < 0) return NULL;
        want = n->data.integer;
    }
    int64_t len = ll_length(src);
    LinkedList* out = ll_new();
    for (int64_t lo = 0; lo < len && want != 0; lo += LL_STREAM_BLOCK) {
        int64_t n = len - lo < LL_STREAM_BLOCK ? len - lo : LL_STREAM_BLOCK;
        Expr* a[3] = { ll_block(src, lo, n), expr_copy(crit),
                       want > 0 ? expr_new_integer(want) : NULL };
        Expr* call = sym_call(SYM_Select, a, want > 0 ? 3 : 2);
        Expr* r = evaluate(call);
        expr_free(call);
        int64_t before = ll_length(out);
        if (!stream_take(out, r)) { ll_release(out); return r; }
        expr_free(r);
        if (want > 0) want -= ll_length(out) - before;
    }
    return ll_wrap(out);
}

/* The leading length of a threaded operand of a Listable call: a LinkedList,
 * a List, or an array (packed or visible). -1 for a scalar operand. */
static int64_t thread_len(const Expr* a) {
    if (a->type == EXPR_LINKEDLIST) return ll_length(a->data.llist);
    if (is_plain_list(a)) return (int64_t)a->data.function.arg_count;
    if (is_ndarray(a) && a->data.ndarray.rank >= 1) return a->data.ndarray.dims[0];
    return -1;
}

/* Elements [lo, lo + n) of a threaded operand, as the List machinery takes
 * them: a LinkedList block (packed when typed), a List of the shared
 * elements, or the array's leading-axis slice (one memcpy of whole rows). */
static Expr* thread_slice(const Expr* a, int64_t lo, int64_t n) {
    if (a->type == EXPR_LINKEDLIST) return ll_block(a->data.llist, lo, n);
    if (is_plain_list(a)) {
        Expr** v = malloc(sizeof(Expr*) * (size_t)n);
        for (int64_t i = 0; i < n; i++) v[i] = expr_copy(a->data.function.args[lo + i]);
        Expr* r = sym_call(SYM_List, v, (size_t)n);
        free(v);
        return r;
    }
    const NDArrayData* d = &a->data.ndarray;
    size_t row = ndt_elem_size(d->dtype);
    for (int k = 1; k < d->rank; k++) row *= (size_t)d->dims[k];
    int64_t dims[NDARRAY_MAX_RANK];
    for (int k = 0; k < d->rank; k++) dims[k] = d->dims[k];
    dims[0] = n;
    void* buf = malloc(row * (size_t)n);
    if (!buf) { fprintf(stderr, "LinkedList: out of memory\n"); abort(); }
    memcpy(buf, (const unsigned char*)d->data + row * (size_t)lo, row * (size_t)n);
    return expr_new_ndarray_like(a, d->rank, dims, buf, d->dtype);
}

/* A Listable call with at least one LinkedList argument threads over it, as
 * over a List, and answers with a LinkedList. Every threaded operand
 * (LinkedList, List, array) must have the same length, else Thread::tdlen and
 * the call settles as it stands (expr_copy). NULL: no LinkedList argument. */
Expr* ll_thread_listable(Expr* res) {
    size_t argc = res->data.function.arg_count;
    int64_t len = -1;
    bool any_ll = false;
    for (size_t i = 0; i < argc; i++) {
        const Expr* a = res->data.function.args[i];
        if (a->type == EXPR_LINKEDLIST) any_ll = true;
        int64_t li = thread_len(a);
        if (li < 0) continue;
        if (len >= 0 && li != len) {
            char* s = expr_to_string_outputform(res);
            mth_message("Thread", "tdlen", "Objects of unequal length in %s cannot be combined.",
                        s ? s : "");
            free(s);
            return expr_copy(res);
        }
        len = li;
    }
    if (!any_ll) return NULL;
    { Expr* d = ll_listable_direct(res); if (d) return d; }
    LinkedList* out = ll_new();
    Expr** a = malloc(sizeof(Expr*) * (argc ? argc : 1));
    for (int64_t lo = 0; lo < len; lo += LL_STREAM_BLOCK) {
        int64_t n = len - lo < LL_STREAM_BLOCK ? len - lo : LL_STREAM_BLOCK;
        for (size_t i = 0; i < argc; i++) {
            const Expr* x = res->data.function.args[i];
            a[i] = thread_len(x) >= 0 ? thread_slice(x, lo, n) : expr_copy((Expr*)x);
        }
        Expr* call = expr_new_function(expr_copy(res->data.function.head), a, argc);
        Expr* r = evaluate(call);
        expr_free(call);
        if (!stream_take(out, r)) { free(a); ll_release(out); return r; }
        expr_free(r);
    }
    free(a);
    return ll_wrap(out);
}
