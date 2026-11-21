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
 * depth_from_bottom:
 * Distance from e to its deepest leaf: 0 for a leaf, 1 + max child depth for a
 * list (an empty list counts as 1). Used to pick the default separator (space
 * at depth 1, one extra newline per level above).
 */
static size_t depth_from_bottom(const Expr* e) {
    if (!is_list(e)) return 0;
    size_t max_child = 0;
    for (size_t i = 0, n = seq_len(e); i < n; i++) {
        bool own;
        Expr* x = seq_at(e, i, &own);
        size_t d = depth_from_bottom(x);
        seq_done(x, own);
        if (d > max_child) max_child = d;
    }
    return max_child + 1;
}

/*
 * riffle_build:
 * Appends the riffled form of e to `out` -- one buffer for the whole result, so
 * every leaf is copied exactly once whatever the nesting depth. level is the
 * 0-based nesting level from the top; df is the distance from bottom (used for
 * the default separator when no explicit separator covers this level). Leaves
 * that are not strings appear as ToString would show them (OutputForm).
 * Returns 0, or -1 on allocation failure.
 */
static int riffle_build(Expr* e, const SepSpec* seps, size_t nseps,
                        size_t level, size_t df, RegexBuf* out) {
    if (!is_list(e)) {
        if (e->type == EXPR_STRING)
            return regexbuf_add(out, e->data.string, strlen(e->data.string));
        char* t = expr_to_string_outputform(e);
        if (!t) return -1;
        int rc = regexbuf_add(out, t, strlen(t));
        free(t);
        return rc;
    }

    /* Resolve the separator for this level. An explicit separator wins; else
     * fall back to the default scheme (space at the bottom, newlines above). */
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
    size_t n = seq_len(e);
    for (size_t i = 0; i < n; i++) {
        if (i > 0) {
            if (sep) {
                if (regexbuf_add(out, sep, sep_len)) return -1;
            } else {
                for (size_t k = 0; k < sep_len; k++)
                    if (regexbuf_add(out, "\n", 1)) return -1;
            }
        }
        bool own;
        Expr* x = seq_at(e, i, &own);
        int rc = riffle_build(x, seps, nseps, level + 1, df ? df - 1 : 0, out);
        seq_done(x, own);
        if (rc) return -1;
    }
    return regexbuf_add(out, right, strlen(right));
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
    RegexBuf out = { NULL, 0, 0 };
    int rc = riffle_build(args[0], seps, nseps, 0, depth, &out);
    free(seps);
    if (rc) { free(out.p); return NULL; }
    if (!out.p) return expr_new_string("");
    return expr_new_string_adopt(out.p);
}
