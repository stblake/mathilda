#include "list_common.h"
#include "array_pad.h"
#include "pad_schemes.h"
#include "options.h"
#include "ndarray.h"     /* is_ndarray, ndarray_to_nested_list */
#include "ndstruct.h"    /* ndstruct_arraypad */
#include "arithmetic.h"  /* is_infinity_sym */

/* ArrayPad[array, m]            pad m elements of padding on every side/level.
 * ArrayPad[array, {m, n}]       m before, n after, on every dimension.
 * ArrayPad[array, {{m1,n1},...}] per-level amounts ({mi} means {mi,mi}).
 * ArrayPad[array, amounts, padding]  use the given padding scheme (default 0).
 *
 * A negative amount removes elements from that side.  Padding may be a constant,
 * a cyclic list of constants, or one of the named schemes handled by
 * pad_schemes.c ("Fixed"/"Periodic"/"Reflected"/"Reversed"/"ReversedNegation"/
 * "ReflectedDifferences"/"ReversedDifferences"/"Extrapolated", the last taking
 * InterpolationOrder). */

#define AP_MAX_RANK 64

/* Nesting depth of a nested-List array (0 for an atom). {{},{}} has depth 2. */
static size_t ap_depth(const Expr* e) {
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL
        || e->data.function.head->data.symbol.name != SYM_List)
        return 0;
    size_t best = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        size_t d = ap_depth(e->data.function.args[i]);
        if (d > best) best = d;
    }
    return 1 + best;
}

/* Rectangular dimension along each level, following the first-child spine. */
static void ap_dims(const Expr* e, int64_t* dims, size_t num_levels, size_t level) {
    if (level >= num_levels) return;
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL
        || e->data.function.head->data.symbol.name != SYM_List) {
        dims[level] = 0; return;
    }
    dims[level] = (int64_t)e->data.function.arg_count;
    if (e->data.function.arg_count > 0)
        ap_dims(e->data.function.args[0], dims, num_levels, level + 1);
    else
        for (size_t k = level + 1; k < num_levels; k++) dims[k] = 0;
}

/* True for a List[...] expression. */
static bool ap_is_list(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List;
}

/* ---- amounts parsing --------------------------------------------------- *
 * Fill lo[]/hi[] for each of *num_levels levels.  Returns false on a bad spec. */
static bool ap_parse_amounts(const Expr* spec, size_t rank,
                             int64_t* lo, int64_t* hi, size_t* num_levels) {
    if (spec->type == EXPR_INTEGER) {
        for (size_t k = 0; k < rank; k++) { lo[k] = hi[k] = spec->data.integer; }
        *num_levels = rank;
        return true;
    }
    if (!ap_is_list(spec)) return false;
    size_t n = spec->data.function.arg_count;
    if (n == 0) return false;
    bool list_of_lists = ap_is_list(spec->data.function.args[0]);
    if (list_of_lists) {
        if (n > AP_MAX_RANK) return false;
        for (size_t k = 0; k < n; k++) {
            const Expr* in = spec->data.function.args[k];
            if (!ap_is_list(in)) return false;
            size_t ic = in->data.function.arg_count;
            const Expr* a = (ic >= 1) ? in->data.function.args[0] : NULL;
            const Expr* b = (ic == 2) ? in->data.function.args[1] : a;
            if (ic < 1 || ic > 2 || a->type != EXPR_INTEGER
                || b->type != EXPR_INTEGER) return false;
            lo[k] = a->data.integer; hi[k] = b->data.integer;
        }
        *num_levels = n;
        return true;
    }
    /* flat integer list: {m} -> {m,m}, {m,n} -> {m,n}, applied to every dim */
    if (n > 2) return false;
    const Expr* a = spec->data.function.args[0];
    const Expr* b = (n == 2) ? spec->data.function.args[1] : a;
    if (a->type != EXPR_INTEGER || b->type != EXPR_INTEGER) return false;
    for (size_t k = 0; k < rank; k++) { lo[k] = a->data.integer; hi[k] = b->data.integer; }
    *num_levels = rank;
    return true;
}

/* ---- constant / cyclic builder ----------------------------------------- *
 * Cyclically index the padding block by the accumulated coordinates (matches
 * PadLeft/PadRight's pr_pad_at). A scalar block ignores coords. */
static Expr* ap_pad_at(const Expr* block, const int64_t* coords, size_t n) {
    if (!block) return expr_new_integer(0);
    if (n == 0 || !ap_is_list(block) || block->data.function.arg_count == 0)
        return expr_copy((Expr*)block);
    int64_t len = (int64_t)block->data.function.arg_count;
    int64_t idx = ((coords[0] % len) + len) % len;
    return ap_pad_at(block->data.function.args[idx], coords + 1, n - 1);
}

/* Build the padded array for a constant/cyclic block. `node` is the original
 * sub-expression at this position (NULL for a pure-padding subtree). Output
 * width at each level is orig_dim[level]+lo+hi (rectangular target). */
static Expr* ap_build_const(const Expr* node, const int64_t* lo, const int64_t* hi,
                            const int64_t* orig_dim, size_t num_levels, size_t level,
                            const Expr* block, int64_t* coords) {
    if (level == num_levels)
        return node ? expr_copy((Expr*)node) : ap_pad_at(block, coords, num_levels);

    int64_t N = orig_dim[level] + lo[level] + hi[level];
    if (N < 0) N = 0;
    int64_t L = (node && ap_is_list(node)) ? (int64_t)node->data.function.arg_count : 0;
    Expr* head = (node && ap_is_list(node))
                 ? expr_copy(node->data.function.head) : expr_new_symbol(SYM_List);
    if (N == 0) return expr_new_function(head, NULL, 0);

    Expr** out = malloc(sizeof(Expr*) * (size_t)N);
    for (int64_t i = 0; i < N; i++) {
        int64_t src = i - lo[level];               /* original element index */
        const Expr* child = (node && ap_is_list(node) && src >= 0 && src < L)
                            ? node->data.function.args[src] : NULL;
        coords[level] = src;
        out[i] = ap_build_const(child, lo, hi, orig_dim, num_levels, level + 1,
                                block, coords);
    }
    Expr* r = expr_new_function(head, out, (size_t)N);
    free(out);
    return r;
}

/* ---- value-dependent builder (outer-to-inner) -------------------------- *
 * At each level, extend the fiber of sub-nodes by the scheme, then recurse into
 * every resulting sub-node to pad the deeper levels. Returns NULL on failure
 * (e.g. a too-short axis for a difference scheme). */
static Expr* ap_pad_valdep(const Expr* node, const int64_t* lo, const int64_t* hi,
                           size_t num_levels, size_t level, PadScheme sc,
                           int64_t order) {
    if (level == num_levels || !ap_is_list(node)) return expr_copy((Expr*)node);

    int64_t L = (int64_t)node->data.function.arg_count;
    int64_t l = lo[level], h = hi[level];
    int64_t start = l < 0 ? -l : 0;
    int64_t end = L - (h < 0 ? -h : 0);
    if (end < start) end = start;
    int64_t src_len = end - start;
    int64_t front = l > 0 ? l : 0;
    int64_t back = h > 0 ? h : 0;
    /* A value-dependent scheme cannot derive padding from an empty axis, so
     * such a dimension is left unpadded (matching Wolfram's ArrayPad). */
    if (src_len == 0) { front = 0; back = 0; }

    const Expr** src = (const Expr**)malloc(sizeof(Expr*) * (size_t)(src_len == 0 ? 1 : src_len));
    for (int64_t j = 0; j < src_len; j++) src[j] = node->data.function.args[start + j];

    Expr** cur;            /* owned when padded, else borrowed slice */
    size_t cur_len;
    bool cur_owned;
    if (front > 0 || back > 0) {
        if (!pad_scheme_extend(src, (size_t)src_len, front, back, NULL, sc, order,
                               &cur, &cur_len)) {
            free(src);
            return NULL;
        }
        cur_owned = true;
    } else {
        cur = (Expr**)src;           /* borrowed slice */
        cur_len = (size_t)src_len;
        cur_owned = false;
    }

    Expr** out = malloc(sizeof(Expr*) * (cur_len == 0 ? 1 : cur_len));
    bool failed = false;
    size_t built = 0;
    for (size_t i = 0; i < cur_len; i++) {
        out[i] = ap_pad_valdep(cur[i], lo, hi, num_levels, level + 1, sc, order);
        if (!out[i]) { failed = true; break; }
        built++;
    }

    Expr* r = NULL;
    if (!failed)
        r = expr_new_function(expr_copy(node->data.function.head), out, cur_len);
    else
        for (size_t i = 0; i < built; i++) expr_free(out[i]);
    free(out);
    if (cur_owned) { for (size_t i = 0; i < cur_len; i++) expr_free(cur[i]); free(cur); }
    free(src);
    return r;
}

/* Emit `ArrayPad::mindimsize` when a value-difference scheme is asked for
 * positive padding on an axis of length < 2. */
static Expr* ap_mindimsize(const Expr* scheme, const Expr* amounts, const Expr* arr) {
    /* Wolfram prints the scheme name without quotes in this message. */
    char* ss = (scheme && scheme->type == EXPR_STRING)
               ? NULL : expr_to_string((Expr*)scheme);
    const char* sname = ss ? ss : scheme->data.string;
    char* as = expr_to_string((Expr*)amounts);
    char* xs = expr_to_string((Expr*)arr);
    fprintf(stderr,
            "ArrayPad::mindimsize: With padding %s, the padding amount %s should "
            "specify positive padding only in dimensions of length at least 2 in "
            "array %s.\n", sname, as, xs);
    free(ss); free(as); free(xs);
    return NULL;
}

Expr* builtin_array_pad(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;

    /* Strip a trailing InterpolationOrder option. */
    const Expr* ord_opt = NULL;
    bool ord_given = false;
    const OptEntry ents[1] = { { "InterpolationOrder", &ord_opt, &ord_given } };
    size_t argc = res->data.function.arg_count;
    if (!options_extract(res, "ArrayPad", ents, 1, &argc)) return NULL;
    if (argc < 2 || argc > 3) return builtin_arg_error("ArrayPad", argc, 2, 3);

    int64_t order = 1;                     /* default linear extrapolation */
    if (ord_given && ord_opt) {
        if (is_infinity_sym((Expr*)ord_opt)) order = -1;
        else if (ord_opt->type == EXPR_INTEGER && ord_opt->data.integer >= 0)
            order = ord_opt->data.integer;
        else return NULL;                  /* bad InterpolationOrder value */
    }

    Expr* array = res->data.function.args[0];
    const Expr* amounts = res->data.function.args[1];
    const Expr* padding = (argc == 3) ? res->data.function.args[2] : NULL;

    /* NDArray input: rank-1 constant fast path, else materialize once. */
    Expr* owned_array = NULL;
    if (is_ndarray(array)) {
        Expr* fast = ndstruct_arraypad(res);
        if (fast) return fast;
        owned_array = ndarray_to_nested_list(array);
        if (!owned_array) return NULL;
        array = owned_array;
    }

    if (!ap_is_list(array)) { if (owned_array) expr_free(owned_array); return NULL; }

    size_t rank = ap_depth(array);
    if (rank == 0 || rank > AP_MAX_RANK) { if (owned_array) expr_free(owned_array); return NULL; }

    int64_t lo[AP_MAX_RANK], hi[AP_MAX_RANK];
    for (size_t k = 0; k < AP_MAX_RANK; k++) { lo[k] = hi[k] = 0; }
    size_t num_levels = 0;
    if (!ap_parse_amounts(amounts, rank, lo, hi, &num_levels)) {
        if (owned_array) expr_free(owned_array);
        return NULL;
    }

    PadScheme sc = pad_scheme_classify(padding);
    if (sc == PS_INVALID) { if (owned_array) expr_free(owned_array); return NULL; }

    Expr* result;
    if (!pad_scheme_value_dependent(sc)) {
        int64_t orig_dim[AP_MAX_RANK];
        for (size_t k = 0; k < AP_MAX_RANK; k++) orig_dim[k] = 0;
        ap_dims(array, orig_dim, num_levels, 0);
        int64_t coords[AP_MAX_RANK];
        Expr* block = padding ? (Expr*)padding : NULL;
        result = ap_build_const(array, lo, hi, orig_dim, num_levels, 0, block, coords);
    } else {
        /* mindimsize guard for the difference/reflected schemes. */
        if (pad_scheme_needs_two(sc)) {
            int64_t dm[AP_MAX_RANK];
            for (size_t k = 0; k < AP_MAX_RANK; k++) dm[k] = 0;
            ap_dims(array, dm, num_levels, 0);
            for (size_t k = 0; k < num_levels; k++)
                if ((lo[k] > 0 || hi[k] > 0) && dm[k] < 2) {
                    if (owned_array) expr_free(owned_array);
                    return ap_mindimsize(padding, amounts, array);
                }
        }
        result = ap_pad_valdep(array, lo, hi, num_levels, 0, sc, order);
    }

    if (owned_array) expr_free(owned_array);
    return result;   /* may be NULL -> leave unevaluated */
}
