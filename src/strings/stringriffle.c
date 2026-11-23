/*
 * stringriffle.c - StringRiffle builtin for Mathilda
 *
 * StringRiffle assembles a string from a (possibly nested) list of elements by
 * inserting separators between them - the inverse of StringSplit. Non-string
 * leaves are converted with expr_to_string (so the integer 27 -> "27"); string
 * leaves are used verbatim.
 *
 * Forms:
 *   StringRiffle[list]                       - default separator scheme: a
 *       single space at the innermost level, one extra newline per level going
 *       up (2-D: rows by "\n", cells by " "; 3-D: blocks by "\n\n", ...).
 *   StringRiffle[list, sep]                  - string sep between the top-level
 *       elements; deeper levels fall back to the default scheme.
 *   StringRiffle[list, {"l","sep","r"}]      - a 3-string list is a delimiter
 *       triple: wrap the join with "l"..."r", joining with "sep".
 *   StringRiffle[list, sep1, sep2, ...]      - sep_i (a string or a delimiter
 *       triple) between elements at level i (1 = outermost); deeper levels use
 *       the default scheme.
 *
 * Strings are treated as raw byte arrays (consistent with the rest of the
 * string subsystem); no UTF-8 decoding is performed.
 */

#include "picostrings.h"
#include "strings_internal.h"
#include "regex_common.h"   /* RegexBuf */
#include "common.h"
#include "sym_names.h"
#include "linkedlist.h"   /* a LinkedList riffles as its elements */
#include "print.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

/* True iff e is a List[...] expression or a LinkedList (riffled as the List
 * of its elements, read from the slots). */
static bool is_list(const Expr* e) {
    return e->type == EXPR_LINKEDLIST
        || (e->type == EXPR_FUNCTION
            && e->data.function.head->type == EXPR_SYMBOL
            && e->data.function.head->data.symbol.name == SYM_List);
}

/* Element count and element i of a list (is_list): a List's argument, or a
 * LinkedList's slot viewed through ll_elem_view (release with seq_done). */
static size_t seq_len(const Expr* e) {
    return e->type == EXPR_LINKEDLIST ? (size_t)ll_length(e->data.llist)
                                      : e->data.function.arg_count;
}
static Expr* seq_at(const Expr* e, size_t i, bool* owned) {
    if (e->type == EXPR_LINKEDLIST) return ll_elem_view(e->data.llist, (int64_t)i, owned);
    *owned = false;
    return e->data.function.args[i];
}
static void seq_done(Expr* x, bool owned) { ll_elem_done(x, owned); }

/*
 * SepSpec:
 * A resolved separator for one nesting level. left/right are the delimiters
 * that wrap the whole level ("" for a plain separator); sep goes between
 * successive elements. All three point into borrowed argument strings that
 * remain valid for the duration of the call.
 */
typedef struct {
    const char* left;
    const char* sep;
    const char* right;
} SepSpec;

/*
 * parse_sep:
 * Validates a user-supplied separator argument and fills *out. A plain string
 * becomes {"", string, ""}; a 3-element List of strings becomes its
 * {left, sep, right}. Returns false (leaving the call unevaluated) for anything
 * else.
 */
static bool parse_sep(Expr* e, SepSpec* out) {
    if (e->type == EXPR_STRING) {
        out->left = "";
        out->sep = e->data.string;
        out->right = "";
        return true;
    }
    if (is_list(e) && seq_len(e) == 3) {
        /* Strings are boxed: a LinkedList's slots are borrowed in place. */
        const Expr* a[3];
        for (size_t i = 0; i < 3; i++) {
            bool own;
            Expr* x = seq_at(e, i, &own);
            bool str = x->type == EXPR_STRING;
            seq_done(x, own);
            if (!str) return false;
            a[i] = e->type == EXPR_LINKEDLIST ? ll_peek_boxed(e->data.llist, (int64_t)i)
                                              : e->data.function.args[i];
        }
        out->left = a[0]->data.string;
        out->sep = a[1]->data.string;
        out->right = a[2]->data.string;
        return true;
    }
    return false;
}

/*
 * RiffleFrame / RiffleStack:
 * One open list on the explicit traversal stack shared by depth_from_bottom
 * and riffle_build, which are iterative so a list nested to any depth (e.g.
 * Nest[List, "a", 100000]) is riffled without exhausting the C stack (audit
 * StringRiffle#1). `e` is the list, `i` the next element to visit; `own` says
 * whether `e` is a box ll_elem_view lent (released when the frame pops). The
 * remaining fields are riffle_build's resolved separator for this level and
 * depth_from_bottom's deepest child so far.
 */
typedef struct {
    Expr* e;
    bool own;
    size_t i, n, level, df;
    const char* sep;     /* NULL => sep_len newlines */
    size_t sep_len;
    const char* right;
    size_t max_child;
} RiffleFrame;

typedef struct { RiffleFrame* v; size_t n, cap; } RiffleStack;

/* Push a frame for list e; returns the new frame, or NULL on allocation
 * failure (the caller still owns e then). */
static RiffleFrame* rs_push(RiffleStack* st, Expr* e, bool own) {
    if (st->n == st->cap) {
        size_t nc = st->cap ? 2 * st->cap : 32;
        RiffleFrame* g = realloc(st->v, nc * sizeof(RiffleFrame));
        if (!g) return NULL;
        st->v = g; st->cap = nc;
    }
    RiffleFrame* f = &st->v[st->n++];
    memset(f, 0, sizeof *f);
    f->e = e;
    f->own = own;
    f->n = seq_len(e);
    return f;
}

/* Release every frame still on the stack (an early exit) and the stack. */
static void rs_free(RiffleStack* st) {
    for (size_t k = st->n; k-- > 0;) seq_done(st->v[k].e, st->v[k].own);
    free(st->v);
}

/*
 * depth_from_bottom:
 * Distance from e to its deepest leaf: 0 for a leaf, 1 + max child depth for a
 * list (an empty list counts as 1). Used to pick the default separator (space
 * at depth 1, one extra newline per level above). Returns SIZE_MAX on
 * allocation failure.
 */
static size_t depth_from_bottom(Expr* e) {
    if (!is_list(e)) return 0;
    RiffleStack st = { NULL, 0, 0 };
    if (!rs_push(&st, e, false)) return SIZE_MAX;
    size_t d = 0;
    while (st.n) {
        RiffleFrame* f = &st.v[st.n - 1];
        if (f->i == f->n) {
            d = f->max_child + 1;
            seq_done(f->e, f->own);
            if (--st.n && d > st.v[st.n - 1].max_child) st.v[st.n - 1].max_child = d;
            continue;
        }
        bool own;
        Expr* x = seq_at(f->e, f->i++, &own);
        if (!is_list(x)) { seq_done(x, own); continue; }
        if (!rs_push(&st, x, own)) { seq_done(x, own); rs_free(&st); return SIZE_MAX; }
    }
    free(st.v);
    return d;
}

/* Append a leaf as it appears in the result: a string verbatim, anything else
 * as ToString would show it (OutputForm). Returns 0, or -1 on failure. */
static int riffle_leaf(Expr* e, RegexBuf* out) {
    if (e->type == EXPR_STRING)
        return regexbuf_add(out, e->data.string, strlen(e->data.string));
    char* t = expr_to_string_outputform(e);
    if (!t) return -1;
    int rc = regexbuf_add(out, t, strlen(t));
    free(t);
    return rc;
}

/* Open list e at nesting `level` (0-based from the top) with distance-from-
 * bottom df: resolve its separator, emit its left delimiter and push its frame.
 * An explicit separator wins; else the default scheme (space at the bottom,
 * df-1 newlines above). Returns 0, or -1 on failure (e then still the
 * caller's). */
static int riffle_open(RiffleStack* st, Expr* e, bool own, const SepSpec* seps,
                       size_t nseps, size_t level, size_t df, RegexBuf* out) {
    const char* left = "";
    const char* right = "";
    const char* sep = " ";
    size_t sep_len = 1;
    if (level < nseps) {
        left = seps[level].left;
        sep = seps[level].sep;
        right = seps[level].right;
        sep_len = strlen(sep);
    } else if (df > 1) {
        sep = NULL;                    /* df-1 newlines between elements */
        sep_len = df - 1;
    }
    if (regexbuf_add(out, left, strlen(left))) return -1;
    RiffleFrame* f = rs_push(st, e, own);
    if (!f) return -1;
    f->level = level;
    f->df = df;
    f->sep = sep;
    f->sep_len = sep_len;
    f->right = right;
    return 0;
}

/*
 * riffle_build:
 * Appends the riffled form of e to `out` -- one buffer for the whole result, so
 * every leaf is copied exactly once whatever the nesting depth. `depth` is e's
 * distance from bottom (depth_from_bottom), which picks the default separator
 * where no explicit separator covers a level. Returns 0, or -1 on allocation
 * failure.
 */
static int riffle_build(Expr* e, const SepSpec* seps, size_t nseps, size_t depth,
                        RegexBuf* out) {
    if (!is_list(e)) return riffle_leaf(e, out);
    RiffleStack st = { NULL, 0, 0 };
    if (riffle_open(&st, e, false, seps, nseps, 0, depth, out)) { free(st.v); return -1; }
    while (st.n) {
        RiffleFrame* f = &st.v[st.n - 1];
        if (f->i == f->n) {
            int rc = regexbuf_add(out, f->right, strlen(f->right));
            seq_done(f->e, f->own);
            st.n--;
            if (rc) { rs_free(&st); return -1; }
            continue;
        }
        if (f->i > 0) {
            int rc = 0;
            if (f->sep) {
                rc = regexbuf_add(out, f->sep, f->sep_len);
            } else {
                for (size_t k = 0; k < f->sep_len && !rc; k++)
                    rc = regexbuf_add(out, "\n", 1);
            }
            if (rc) { rs_free(&st); return -1; }
        }
        bool own;
        Expr* x = seq_at(f->e, f->i++, &own);
        if (!is_list(x)) {
            int rc = riffle_leaf(x, out);
            seq_done(x, own);
            if (rc) { rs_free(&st); return -1; }
            continue;
        }
        size_t level = f->level + 1, df = f->df ? f->df - 1 : 0;
        if (riffle_open(&st, x, own, seps, nseps, level, df, out)) {
            seq_done(x, own);
            rs_free(&st);
            return -1;
        }
    }
    free(st.v);
    return 0;
}

/*
 * builtin_stringriffle:
 * See file header for the accepted forms. Emits StringRiffle::list /
 * StringRiffle::string and returns NULL (unevaluated) for a
 * first argument that is neither a list nor a string, a malformed separator, or
 * allocation failure. StringRiffle[] with no arguments emits StringRiffle::argm
 * (via builtin_arg_error) and leaves the call unevaluated.
 */
Expr* builtin_stringriffle(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;

    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc == 0)
        return builtin_arg_error("StringRiffle", argc, 1, SIZE_MAX);

    /* The first argument is the data; the rest are per-level separators. */
    size_t nseps = argc - 1;
    SepSpec* seps = malloc(nseps ? nseps * sizeof(SepSpec) : 1);
    if (!seps) return NULL;
    for (size_t i = 0; i < nseps; i++) {
        if (!parse_sep(args[i + 1], &seps[i])) {
            free(seps);
            return strings_arg_message("StringRiffle", "string", "String", i + 2, res);
        }
    }

    /* Only a list or a bare string can be riffled. */
    if (!is_list(args[0]) && args[0]->type != EXPR_STRING) {
        free(seps);
        return strings_arg_message("StringRiffle", "list", "List", 1, res);
    }

    size_t depth = depth_from_bottom(args[0]);
    if (depth == SIZE_MAX) { free(seps); return NULL; }
    RegexBuf out = { NULL, 0, 0 };
    int rc = riffle_build(args[0], seps, nseps, depth, &out);
    free(seps);
    if (rc) { free(out.p); return NULL; }
    if (!out.p) return expr_new_string("");
    return expr_new_string_adopt(out.p);
}
