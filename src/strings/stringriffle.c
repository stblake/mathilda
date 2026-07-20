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
#include "common.h"
#include "sym_names.h"
#include "print.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

/* True iff e is a List[...] expression. */
static bool is_list(const Expr* e) {
    return e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List;
}

/*
 * leaf_to_str:
 * Renders a non-list element as a freshly malloc'd, NUL-terminated C string.
 * String elements are copied verbatim (no surrounding quotes); every other
 * expression is rendered via expr_to_string (StandardForm). Returns NULL on
 * allocation failure.
 */
static char* leaf_to_str(Expr* e) {
    if (e->type == EXPR_STRING) {
        size_t len = strlen(e->data.string);
        char* out = malloc(len + 1);
        if (!out) return NULL;
        memcpy(out, e->data.string, len + 1);
        return out;
    }
    /* expr_to_string returns a heap string the caller must free. */
    return expr_to_string(e);
}

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
    if (is_list(e) && e->data.function.arg_count == 3) {
        Expr** a = e->data.function.args;
        if (a[0]->type != EXPR_STRING || a[1]->type != EXPR_STRING
            || a[2]->type != EXPR_STRING)
            return false;
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
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        size_t d = depth_from_bottom(e->data.function.args[i]);
        if (d > max_child) max_child = d;
    }
    return max_child + 1;
}

/*
 * riffle_build:
 * Recursively assembles the riffled string for e. level is the 0-based nesting
 * level from the top; df is the distance from bottom (used for the default
 * separator when no explicit separator covers this level). Returns a freshly
 * malloc'd string, or NULL on failure (propagated up).
 */
static char* riffle_build(Expr* e, const SepSpec* seps, size_t nseps,
                          size_t level, size_t df) {
    if (!is_list(e))
        return leaf_to_str(e);

    size_t n = e->data.function.arg_count;

    /* Build each child string first. */
    char** parts = malloc(n ? n * sizeof(char*) : 1);
    if (!parts) return NULL;
    for (size_t i = 0; i < n; i++) {
        parts[i] = riffle_build(e->data.function.args[i], seps, nseps,
                                level + 1, df ? df - 1 : 0);
        if (!parts[i]) {
            for (size_t j = 0; j < i; j++) free(parts[j]);
            free(parts);
            return NULL;
        }
    }

    /* Resolve the separator for this level. An explicit separator wins; else
     * fall back to the default scheme (space at the bottom, newlines above). */
    SepSpec spec;
    char* gen = NULL;                 /* owned buffer for a default newline run */
    if (level < nseps) {
        spec = seps[level];
    } else {
        spec.left = "";
        spec.right = "";
        if (df <= 1) {
            spec.sep = " ";
        } else {
            size_t nl = df - 1;       /* df-1 newlines between elements */
            gen = malloc(nl + 1);
            if (!gen) {
                for (size_t i = 0; i < n; i++) free(parts[i]);
                free(parts);
                return NULL;
            }
            memset(gen, '\n', nl);
            gen[nl] = '\0';
            spec.sep = gen;
        }
    }

    /* Two-pass assemble: left + p0 + sep + p1 + ... + right. */
    size_t sep_len = strlen(spec.sep);
    size_t total = strlen(spec.left) + strlen(spec.right);
    for (size_t i = 0; i < n; i++)
        total += strlen(parts[i]);
    if (n > 1)
        total += (n - 1) * sep_len;

    char* buf = malloc(total + 1);
    if (!buf) {
        for (size_t i = 0; i < n; i++) free(parts[i]);
        free(parts);
        free(gen);
        return NULL;
    }

    size_t off = 0;
    size_t ll = strlen(spec.left);
    memcpy(buf + off, spec.left, ll); off += ll;
    for (size_t i = 0; i < n; i++) {
        if (i > 0) { memcpy(buf + off, spec.sep, sep_len); off += sep_len; }
        size_t pl = strlen(parts[i]);
        memcpy(buf + off, parts[i], pl); off += pl;
    }
    size_t rl = strlen(spec.right);
    memcpy(buf + off, spec.right, rl); off += rl;
    buf[off] = '\0';

    for (size_t i = 0; i < n; i++) free(parts[i]);
    free(parts);
    free(gen);
    return buf;
}

/*
 * builtin_stringriffle:
 * See file header for the accepted forms. Returns NULL (unevaluated) for a
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
            return NULL;   /* malformed separator: leave unevaluated */
        }
    }

    /* Only a list or a bare string can be riffled; anything else stays
     * symbolic (e.g. StringRiffle[x]). */
    if (!is_list(args[0]) && args[0]->type != EXPR_STRING) {
        free(seps);
        return NULL;
    }

    size_t depth = depth_from_bottom(args[0]);
    char* buf = riffle_build(args[0], seps, nseps, 0, depth);
    free(seps);
    if (!buf) return NULL;

    Expr* out = expr_new_string(buf);
    free(buf);
    return out;
}
