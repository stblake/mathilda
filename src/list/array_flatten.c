#include "list_common.h"
#include "array_flatten.h"

/* ArrayFlatten[a]      -> a single flattened matrix from a matrix of matrices
 *                         (block matrix). Equivalent to ArrayFlatten[a, 2].
 * ArrayFlatten[a, r]   -> flatten out r pairs of levels of a rank-2r array,
 *                         yielding a rank-r array.
 *
 * `a` is treated as a rank-r *grid* of blocks. Each block is either:
 *   - an "array block": a nested list of array depth >= r, whose first r
 *     dimensions give the block's size along each grid axis (any deeper
 *     structure is carried along verbatim as element content); or
 *   - a "scalar": an atom, a non-List head, or a list shallower than r. A
 *     scalar is replicated to fill a rank-r block of the size demanded by its
 *     grid position (this is how a 0 becomes a zero block).
 *
 * Block sizes must fit: all array blocks sharing grid index i along axis k
 * must agree on their k-th dimension. The output size along axis k is the sum
 * over grid positions i of that shared length (positions with no array block
 * contribute length 1). If two array blocks in the same row/column disagree,
 * the blocks do not fit and the call is left unevaluated.
 *
 * ArrayFlatten[a]    == Flatten[a, {{1,3},{2,4}}]
 * ArrayFlatten[a, r] == Flatten[a, {{1,r+1},{2,r+2},...,{r,2r}}]
 * (both implemented directly here; the list-of-lists Flatten levelspec does
 * not exist in this codebase).
 */

/* Maximum r (number of level pairs). Keeps the fixed-size dims/index buffers
 * below well within the 64-slot recursion contract of af_tensor_dims. */
#define AF_MAX_R 32

/* Local copy of linalg/util.c's get_tensor_dims, kept here so this module has
 * no cross-module link dependency. Returns the rank of a rectangular
 * nested-List tensor (0 for a non-List atom/head), filling dims[0..rank-1], or
 * -1 if the structure is jagged (non-rectangular). dims must hold >= 64 slots. */
static int af_tensor_dims(Expr* e, int64_t* dims) {
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL
        || e->data.function.head->data.symbol.name != SYM_List) {
        return 0; /* rank 0 */
    }
    int64_t len = (int64_t)e->data.function.arg_count;
    dims[0] = len;
    if (len == 0) return 1;

    int sub_rank = af_tensor_dims(e->data.function.args[0], dims + 1);
    for (int64_t i = 1; i < len; i++) {
        int64_t cur_dims[64];
        int cur_rank = af_tensor_dims(e->data.function.args[i], cur_dims);
        if (cur_rank != sub_rank) return -1; /* jagged */
        for (int j = 0; j < sub_rank; j++) {
            if (cur_dims[j] != dims[j + 1]) return -1; /* jagged */
        }
    }
    return sub_rank + 1;
}

/* Recursively walk the block grid down r List levels, collecting borrowed
 * pointers to each block in row-major order. Records the grid dimension at each
 * level in D and enforces that the grid is rectangular. Returns 1 on success,
 * 0 if the structure is not a rectangular rank-r List grid. */
static int af_collect_grid(Expr* node, int level, int r, int64_t* D,
                           Expr*** blocks, size_t* nblocks, size_t* cap) {
    if (level == r) {
        if (*nblocks == *cap) {
            *cap *= 2;
            *blocks = realloc(*blocks, sizeof(Expr*) * (*cap));
        }
        (*blocks)[(*nblocks)++] = node;
        return 1;
    }
    if (node->type != EXPR_FUNCTION || node->data.function.head->type != EXPR_SYMBOL
        || node->data.function.head->data.symbol.name != SYM_List) {
        return 0;
    }
    int64_t len = (int64_t)node->data.function.arg_count;
    if (D[level] == -1) D[level] = len;
    else if (D[level] != len) return 0; /* ragged grid */
    for (int64_t i = 0; i < len; i++) {
        if (!af_collect_grid(node->data.function.args[i], level + 1, r,
                             D, blocks, nblocks, cap)) {
            return 0;
        }
    }
    return 1;
}

/* Assembled per-flatten context passed to af_build. */
typedef struct {
    int r;
    int64_t* O;          /* output dims, length r */
    int64_t** map_i;     /* per axis: output index -> grid index, length O[k] */
    int64_t** map_j;     /* per axis: output index -> within-block index */
    Expr** blocks;       /* row-major grid of borrowed block pointers */
    int* is_arr;         /* per block: 1 if an array block, 0 if scalar */
    int64_t* gstride;    /* grid strides, length r */
    int64_t* o_idx;      /* scratch: current output multi-index, length r */
} AFCtx;

/* Build the rank-r output nested List by recursion over the output axes. */
static Expr* af_build(AFCtx* c, int level) {
    if (level == c->r) {
        int64_t bflat = 0;
        int64_t j[AF_MAX_R];
        for (int k = 0; k < c->r; k++) {
            int64_t o = c->o_idx[k];
            bflat += c->map_i[k][o] * c->gstride[k];
            j[k] = c->map_j[k][o];
        }
        Expr* block = c->blocks[bflat];
        if (c->is_arr[bflat]) {
            Expr* cur = block;
            for (int k = 0; k < c->r; k++) {
                cur = cur->data.function.args[j[k]];
            }
            return expr_copy(cur);
        }
        return expr_copy(block); /* scalar replication */
    }

    int64_t n = c->O[level];
    if (n == 0) return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);

    Expr** args = malloc(sizeof(Expr*) * (size_t)n);
    for (int64_t o = 0; o < n; o++) {
        c->o_idx[level] = o;
        args[o] = af_build(c, level + 1);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), args, (size_t)n);
    free(args);
    return out;
}

Expr* builtin_array_flatten(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;

    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) {
        return builtin_arg_error("ArrayFlatten", argc, 1, 2);
    }

    Expr* a = res->data.function.args[0];
    int64_t r = 2;
    if (argc == 2) {
        Expr* rr = res->data.function.args[1];
        if (rr->type != EXPR_INTEGER || rr->data.integer < 1
            || rr->data.integer > AF_MAX_R) {
            return NULL; /* leave unevaluated on a bad level count */
        }
        r = rr->data.integer;
    }

    /* Collect the rank-r block grid. */
    int64_t D[AF_MAX_R];
    for (int k = 0; k < (int)r; k++) D[k] = -1;

    size_t nblocks = 0, cap = 16;
    Expr** blocks = malloc(sizeof(Expr*) * cap);
    if (!af_collect_grid(a, 0, (int)r, D, &blocks, &nblocks, &cap)) {
        free(blocks);
        return NULL; /* not a rectangular rank-r List grid */
    }

    /* Grid strides for row-major flat block indexing. */
    int64_t gstride[AF_MAX_R];
    gstride[r - 1] = 1;
    for (int k = (int)r - 2; k >= 0; k--) gstride[k] = gstride[k + 1] * D[k + 1];

    /* Per-axis, per-position block length; -1 = not yet determined. */
    int64_t* len[AF_MAX_R];
    for (int k = 0; k < (int)r; k++) {
        len[k] = (D[k] > 0) ? malloc(sizeof(int64_t) * (size_t)D[k]) : NULL;
        for (int64_t i = 0; i < D[k]; i++) len[k][i] = -1;
    }

    int* is_arr = malloc(sizeof(int) * (nblocks > 0 ? nblocks : 1));

    /* Classify each block and derive the fitted length along every axis. */
    int mismatch = 0;
    for (size_t f = 0; f < nblocks && !mismatch; f++) {
        int64_t bdims[64];
        int rank = af_tensor_dims(blocks[f], bdims);
        is_arr[f] = (rank >= (int)r);
        if (!is_arr[f]) continue;
        for (int k = 0; k < (int)r; k++) {
            int64_t idx = ((int64_t)f / gstride[k]) % D[k];
            if (len[k][idx] == -1) len[k][idx] = bdims[k];
            else if (len[k][idx] != bdims[k]) { mismatch = 1; break; }
        }
    }

    if (mismatch) {
        for (int k = 0; k < (int)r; k++) free(len[k]);
        free(is_arr);
        free(blocks);
        return NULL; /* blocks do not fit together */
    }

    /* Positions with no array block contribute length 1; sum to output dims. */
    int64_t O[AF_MAX_R];
    for (int k = 0; k < (int)r; k++) {
        O[k] = 0;
        for (int64_t i = 0; i < D[k]; i++) {
            if (len[k][i] == -1) len[k][i] = 1;
            O[k] += len[k][i];
        }
    }

    /* Per-axis maps from an output index to its (grid index, within-block
     * index). */
    int64_t* map_i[AF_MAX_R];
    int64_t* map_j[AF_MAX_R];
    for (int k = 0; k < (int)r; k++) {
        map_i[k] = (O[k] > 0) ? malloc(sizeof(int64_t) * (size_t)O[k]) : NULL;
        map_j[k] = (O[k] > 0) ? malloc(sizeof(int64_t) * (size_t)O[k]) : NULL;
        int64_t o = 0;
        for (int64_t i = 0; i < D[k]; i++) {
            for (int64_t jj = 0; jj < len[k][i]; jj++) {
                map_i[k][o] = i;
                map_j[k][o] = jj;
                o++;
            }
        }
    }

    int64_t o_idx[AF_MAX_R];
    AFCtx ctx = { (int)r, O, map_i, map_j, blocks, is_arr, gstride, o_idx };
    Expr* result = af_build(&ctx, 0);

    for (int k = 0; k < (int)r; k++) {
        free(len[k]);
        free(map_i[k]);
        free(map_j[k]);
    }
    free(is_arr);
    free(blocks);
    return result;
}
