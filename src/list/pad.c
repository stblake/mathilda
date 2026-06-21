#include "list_common.h"
#include "pad.h"

/* =====================================================================
 * PadRight / PadLeft  --  pad a (possibly ragged, possibly nested) array
 * on the right (resp. left) to a requested length / set of dimensions.
 *
 * The two are exact mirrors and share one recursive engine, pr_build,
 * selected by a pad_left flag and parameterised by:
 *   - dimv[k] : target length at each of k levels (negative => pad on the
 *               opposite side, i.e. on the left, at that level);
 *   - marv[k] : margin of padding to leave on the left at each level
 *               (negative margin truncates leading elements);
 *   - P       : the padding "block" (a scalar, or a nested List that is
 *               indexed cyclically per level).
 *
 * pr_pad_at(P, coords) yields a padding leaf/subtree by cyclically
 * indexing into P with the per-level coordinates accumulated so far.
 * For the default scalar padding 0 this just returns 0 regardless of
 * coordinates.
 * ===================================================================== */

/* Cyclically index into the padding block P using the coordinate path
 * coords[0..n-1].  A non-List P (or running out of coordinates) is a leaf
 * and is returned as a fresh copy. */
static Expr* pr_pad_at(Expr* P, const int64_t* coords, size_t n) {
    if (n == 0 || P->type != EXPR_FUNCTION ||
        P->data.function.head->type != EXPR_SYMBOL ||
        P->data.function.head->data.symbol != SYM_List ||
        P->data.function.arg_count == 0) {
        return expr_copy(P);
    }
    int64_t len = (int64_t)P->data.function.arg_count;
    int64_t idx = ((coords[0] % len) + len) % len;   /* floored modulo */
    return pr_pad_at(P->data.function.args[idx], coords + 1, n - 1);
}

/* Walk a ragged array, recording the maximum width seen at each level.
 * Only EXPR_FUNCTION nodes contribute depth; atoms terminate a branch.
 * dims is grown as needed; *depth is the number of populated levels. */
static void pr_scan_dims(Expr* node, size_t level,
                         int64_t** dims, size_t* depth, size_t* cap) {
    if (node->type != EXPR_FUNCTION) return;
    if (level >= *cap) {
        size_t newcap = *cap ? *cap * 2 : 4;
        *dims = realloc(*dims, sizeof(int64_t) * newcap);
        for (size_t i = *cap; i < newcap; i++) (*dims)[i] = 0;
        *cap = newcap;
    }
    int64_t w = (int64_t)node->data.function.arg_count;
    if (w > (*dims)[level]) (*dims)[level] = w;
    if (level + 1 > *depth) *depth = level + 1;
    for (size_t i = 0; i < node->data.function.arg_count; i++) {
        pr_scan_dims(node->data.function.args[i], level + 1, dims, depth, cap);
    }
}

/* Margin to leave at a given level.  margins may be a single integer
 * (applies to level 0 only) or a List of per-level integers. */
static int64_t pr_margin_at(Expr* margins, size_t level) {
    if (!margins) return 0;
    if (margins->type == EXPR_INTEGER) return level == 0 ? margins->data.integer : 0;
    if (margins->type == EXPR_FUNCTION &&
        margins->data.function.head->type == EXPR_SYMBOL &&
        margins->data.function.head->data.symbol == SYM_List) {
        if (level < margins->data.function.arg_count) {
            Expr* m = margins->data.function.args[level];
            if (m->type == EXPR_INTEGER) return m->data.integer;
        }
    }
    return 0;
}

/* Recursive padding engine, shared by PadRight (pad_left == 0) and PadLeft
 * (pad_left == 1).  node is the original sub-expression occupying this slot
 * (NULL for a pure-padding slot).  coords[0..level-1] is already filled; this
 * call fills coords[level] as it lays out the level. */
static Expr* pr_build(Expr* node, const int64_t* dimv, Expr* margins,
                      size_t k, size_t level, Expr* P, int64_t* coords,
                      int pad_left, int empty_anchor) {
    if (level == k) {
        /* Leaf of the padded structure: keep the original element (whatever
         * its remaining depth) or synthesise a padding leaf. */
        return node ? expr_copy(node) : pr_pad_at(P, coords, k);
    }

    int64_t n = dimv[level];
    size_t  N = (size_t)(n < 0 ? -n : n);
    int64_t m = pr_margin_at(margins, level);

    /* Determine this level's source row: its head, length and a way to
     * fetch child j. */
    Expr*  head_src = NULL;     /* function whose head we copy */
    int    atom_row = 0;        /* node is a non-list atom -> 1-element row */
    size_t L = 0;
    if (node && node->type == EXPR_FUNCTION) {
        head_src = node;
        L = node->data.function.arg_count;
    } else if (node) {
        atom_row = 1;           /* e.g. {a, {b,c}} -> a becomes a 1-row */
        L = 1;
    }

    /* Which side gets the padding?  PadRight pads on the left only for a
     * negative length; PadLeft pads on the left for a non-negative one. */
    int pad_on_left = pad_left ? (n >= 0) : (n < 0);

    /* Placement of the original row inside the N output slots.
     *   pad on left  -> row right-aligned, margin m on the right;
     *   pad on right -> row left-aligned, margin m on the left. */
    int64_t list_start = pad_on_left ? ((int64_t)N - (int64_t)L - m) : m;

    /* Padding-phase anchor.  For right padding the first original element
     * sits on P[0] (coord = i - list_start).  For left padding the last
     * original element sits on P[s-1]; with pr_pad_at's floored modulo that
     * is coord = i - (list_start + L).  For left padding L cancels against
     * list_start, so empty pure-padding sub-rows still share the data rows'
     * tiling phase.  Only a wholly empty *top-level* list (empty_anchor)
     * falls back to laying the sequence from position 0. */
    int64_t coord_base = (L == 0 && empty_anchor)
        ? 0
        : (pad_on_left ? list_start + (int64_t)L : list_start);

    Expr* head = head_src ? expr_copy(head_src->data.function.head)
                          : expr_new_symbol(SYM_List);
    if (N == 0) {
        return expr_new_function(head, NULL, 0);
    }

    Expr** out = malloc(sizeof(Expr*) * N);
    for (size_t i = 0; i < N; i++) {
        int64_t j = (int64_t)i - list_start;     /* index into the original row */
        Expr* child = NULL;
        if (j >= 0 && j < (int64_t)L) {
            child = atom_row ? node : head_src->data.function.args[j];
        }
        coords[level] = (int64_t)i - coord_base;
        out[i] = pr_build(child, dimv, margins, k, level + 1, P, coords,
                          pad_left, /*empty_anchor=*/0);
    }
    Expr* result = expr_new_function(head, out, N);
    free(out);
    return result;
}

/* `<name>::argb: <name> called with N arguments; between 1 and 4 arguments
 * are expected.` */
static Expr* pad_emit_argb(const char* name, size_t argc) {
    fprintf(stderr,
            "%s::argb: %s called with %zu argument%s; "
            "between 1 and 4 arguments are expected.\n",
            name, name, argc, argc == 1 ? "" : "s");
    return NULL;
}

/* Shared driver for PadRight (pad_left == 0) and PadLeft (pad_left == 1). */
static Expr* pad_dispatch(Expr* res, int pad_left, const char* name) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 4) return pad_emit_argb(name, argc);

    Expr* list = res->data.function.args[0];
    if (list->type != EXPR_FUNCTION) return NULL;   /* not an array: leave alone */

    /* ---- dimensions (arg 1) ---- */
    int64_t* dimv = NULL;
    size_t   k = 0;
    int      own_dimv = 0;       /* dimv heap-owned (full-dims path) */
    Expr* nspec = (argc >= 2) ? res->data.function.args[1] : NULL;

    int64_t scalar_dim[1];
    if (!nspec ||
        (nspec->type == EXPR_SYMBOL && nspec->data.symbol == SYM_Automatic)) {
        /* Pad*[list] or Pad*[list, Automatic, ...]: pad to full. */
        size_t cap = 0, depth = 0;
        pr_scan_dims(list, 0, &dimv, &depth, &cap);
        if (depth == 0) {                 /* atom array already handled above */
            free(dimv);
            return expr_copy(list);
        }
        k = depth;
        own_dimv = 1;
    } else if (nspec->type == EXPR_INTEGER) {
        scalar_dim[0] = nspec->data.integer;
        dimv = scalar_dim;
        k = 1;
    } else if (nspec->type == EXPR_FUNCTION &&
               nspec->data.function.head->type == EXPR_SYMBOL &&
               nspec->data.function.head->data.symbol == SYM_List) {
        size_t kk = nspec->data.function.arg_count;
        if (kk == 0) return expr_copy(list);
        dimv = malloc(sizeof(int64_t) * kk);
        for (size_t i = 0; i < kk; i++) {
            Expr* d = nspec->data.function.args[i];
            if (d->type != EXPR_INTEGER) { free(dimv); return NULL; }
            dimv[i] = d->data.integer;
        }
        k = kk;
        own_dimv = 1;
    } else {
        return NULL;   /* unrecognised length spec: leave unevaluated */
    }

    /* ---- padding (arg 2) and margins (arg 3) ---- */
    Expr* zero = NULL;
    Expr* P = (argc >= 3) ? res->data.function.args[2] : (zero = expr_new_integer(0));
    Expr* margins = (argc >= 4) ? res->data.function.args[3] : NULL;

    int64_t* coords = malloc(sizeof(int64_t) * k);
    int empty_anchor = (list->data.function.arg_count == 0);
    Expr* out = pr_build(list, dimv, margins, k, 0, P, coords, pad_left, empty_anchor);

    free(coords);
    if (own_dimv) free(dimv);
    if (zero) expr_free(zero);
    return out;
}

Expr* builtin_padright(Expr* res) {
    return pad_dispatch(res, /*pad_left=*/0, "PadRight");
}

Expr* builtin_padleft(Expr* res) {
    return pad_dispatch(res, /*pad_left=*/1, "PadLeft");
}
