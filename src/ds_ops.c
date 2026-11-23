/* ds_ops.c -- the registry and the rows for the built-in payload kinds.
 * See ds_ops.h.  Each function body is the code that used to be inlined at the
 * hard-coded site, moved here unchanged. */
#include "ds_ops.h"
#include "sym_names.h"
#include "ndarray.h"
#include "linkedlist.h"
#include "sparse/sparse.h"
#include "part.h"

#define DS_MAX_TYPES 64
static const DsOps* g_rows[DS_MAX_TYPES];
static bool g_init = false;

/* ---- packed arrays (EXPR_NDARRAY) --------------------------------------- */
static Expr* nd_head(const Expr* e) {
    /* A packed list IS a List and must say so (see expr_head). */
    return expr_new_symbol(is_packed_list(e) ? SYM_List : SYM_NDArray);
}
static int64_t nd_length(const Expr* e) { return e->data.ndarray.dims[0]; }   /* rank >= 1 */
static int nd_dims(const Expr* e, int64_t* out, int max) {
    int d = e->data.ndarray.rank;
    if (d > max) d = max;
    for (int i = 0; i < d; i++) out[i] = e->data.ndarray.dims[i];
    return d;
}

/* ---- sparse arrays (EXPR_SPARSE) ---------------------------------------- */
static Expr* sp_head(const Expr* e) { (void)e; return expr_new_symbol(SYM_SparseArray); }
static bool sp_head_is(const Expr* e, const char* hn) { (void)e; return hn == SYM_SparseArray; }
static int64_t sp_length(const Expr* e) { return e->data.sparse->dims[0]; }
static int sp_dims(const Expr* e, int64_t* out, int max) {
    int d = e->data.sparse->rank;
    if (d > max) d = max;
    for (int i = 0; i < d; i++) out[i] = e->data.sparse->dims[i];
    return d;
}

/* ---- LinkedList (EXPR_LINKEDLIST) --------------------------------------- */
static Expr* ll_head_sym(const Expr* e) { (void)e; return expr_new_symbol(SYM_LinkedList); }
static bool ll_head_is(const Expr* e, const char* hn) { (void)e; return hn == SYM_LinkedList; }
static int64_t ll_len(const Expr* e) { return ll_length(e->data.llist); }       /* O(1) */

/* ---- element access ------------------------------------------------------ */
static Expr* nd_part_row(Expr* e, Expr** idx, size_t n, int* st) {
    *st = DS_PART_DONE;
    if (n > 0 && idx[0]->type == EXPR_INTEGER && idx[0]->data.integer == 0) {
        *st = DS_PART_DECLINE;       /* [[0]]: head extraction via expr_part */
        return NULL;
    }
    bool degrade = false;
    Expr* r = ndarray_part(e, idx, n, &degrade);
    if (!r && degrade) *st = DS_PART_DEGRADE;
    return r;
}
static Expr* sp_part_row(Expr* e, Expr** idx, size_t n, int* st) {
    bool err = false;
    *st = DS_PART_DONE;
    return sp_part(e, idx, n, &err);
}
static Expr* ll_part_row(Expr* e, Expr** idx, size_t n, int* st) {
    *st = DS_PART_DONE;
    return expr_part(e, idx, n);
}

/* A packed list/array matches _List or _NDArray by its presentation; the
 * matcher keeps its own NDArray head logic, so no head_is is given here. */
static const DsOps ROW_NDARRAY    = { EXPR_NDARRAY,    "NDArray",    nd_head,      NULL,       nd_length, nd_dims, nd_part_row };
static const DsOps ROW_SPARSE     = { EXPR_SPARSE,     "SparseArray", sp_head,     sp_head_is, sp_length, sp_dims, sp_part_row };
static const DsOps ROW_LINKEDLIST = { EXPR_LINKEDLIST, "LinkedList", ll_head_sym,  ll_head_is, ll_len,    NULL, ll_part_row };

static void ds_init(void) {
    g_init = true;
    g_rows[ROW_NDARRAY.type]    = &ROW_NDARRAY;
    g_rows[ROW_SPARSE.type]     = &ROW_SPARSE;
    g_rows[ROW_LINKEDLIST.type] = &ROW_LINKEDLIST;
}

const DsOps* ds_ops(ExprType t) {
    if (!g_init) ds_init();
    return ((unsigned)t < DS_MAX_TYPES) ? g_rows[t] : NULL;
}

bool ds_register(const DsOps* ops) {
    if (!g_init) ds_init();
    if (!ops || (unsigned)ops->type >= DS_MAX_TYPES) return false;
    g_rows[ops->type] = ops;
    return true;
}

Expr* ds_part(Expr* e, Expr** idx, size_t n, int* st) {
    const DsOps* o = ds_ops(e->type);
    if (!o || !o->part) { *st = DS_PART_DECLINE; return NULL; }
    return o->part(e, idx, n, st);
}
