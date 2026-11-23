/* ds_ops.c -- the registry and the rows for the built-in payload kinds.
 * See ds_ops.h.  Each function body is the code that used to be inlined at the
 * hard-coded site, moved here unchanged. */
#include "ds_ops.h"
#include "sym_names.h"
#include "ndarray.h"
#include "linkedlist.h"
#include "sparse/sparse.h"
#include "part.h"
#include "patterns.h"
#include "pvec.h"
#include "pvec_ops.h"

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
    if (n > 0 && idx[0]->type == EXPR_INTEGER && idx[0]->data.integer == 0) {
        *st = DS_PART_DECLINE;          /* [[0]]: head extraction */
        return NULL;
    }
    return ll_part_expr(e, idx, n);
}

/* ---- materialisation (Normal) ------------------------------------------- */
static Expr* nd_to_list(const Expr* e) { return ndarray_to_nested_list(e); }
static Expr* sp_to_list(const Expr* e) { return sp_normal(e->data.sparse); }
static Expr* ll_to_list_row(const Expr* e) { return ll_to_list(e->data.llist); }

/* ---- PVector (EXPR_PVEC) -------------------------------------------------- */
static Expr* pv_head_row(const Expr* e) { (void)e; return expr_new_symbol(SYM_PVector); }
static bool pv_head_is(const Expr* e, const char* hn) { (void)e; return hn == SYM_PVector; }
static int64_t pv_len_row(const Expr* e) { return pv_length(e->data.pvec); }
static int pv_dims_row(const Expr* e, int64_t* out, int max) {
    if (max < 1) return 0;
    out[0] = pv_length(e->data.pvec);
    return 1;
}
static Expr* pv_part_row(Expr* e, Expr** idx, size_t n, int* st) {
    *st = DS_PART_DONE;
    if (n > 0 && idx[0]->type == EXPR_INTEGER && idx[0]->data.integer == 0) {
        *st = DS_PART_DECLINE;          /* [[0]]: head extraction */
        return NULL;
    }
    return pvec_part(e, idx, n);
}
static Expr* pv_to_list_row(const Expr* e) { return pv_to_list(e->data.pvec); }

/* ---- element access and rebuild -------------------------------------------- */
static Expr* ll_get_row(const Expr* e, int64_t i) { return ll_get(e->data.llist, i); }
static Expr* pv_get_row(const Expr* e, int64_t i) { return pv_get(e->data.pvec, i); }
static Expr* ll_from_list_row(const Expr* list) {
    if (!ll_is_list_source(list)) return NULL;
    LinkedList* ll = ll_new();
    if (!ll_append_elements(ll, list)) { ll_release(ll); return NULL; }
    return expr_new_linkedlist(ll);
}
static Expr* pv_from_list_row(const Expr* list) { return pvec_make(list); }

/* ---- iteration -------------------------------------------------------------- */
static bool ll_each_row(const Expr* e, bool (*fn)(Expr*, void*), void* ud) {
    const LinkedList* ll = e->data.llist;
    int64_t n = ll_length(ll);
    bool boxed = ll_kind(ll) == LL_BOXED;
    for (int64_t i = 0; i < n; i++) {
        Expr* owned = boxed ? NULL : ll_get(ll, i);
        Expr* x = owned ? owned : (Expr*)ll_peek_boxed(ll, i);
        bool go = fn(x, ud);
        if (owned) expr_free(owned);
        if (!go) return true;
    }
    return true;
}
static bool ll_atoms_only(const Expr* e) { return ll_kind(e->data.llist) != LL_BOXED; }
static bool pv_each_row(const Expr* e, bool (*fn)(Expr*, void*), void* ud) {
    const PVec* v = e->data.pvec;
    int64_t n = pv_length(v);
    for (int64_t i = 0; i < n; i++) {
        Expr* x = pv_get(v, i);
        bool go = fn(x, ud);
        expr_free(x);
        if (!go) return true;
    }
    return true;
}
static bool pv_atoms_only(const Expr* e) { PVKind k = pv_kind(e->data.pvec); return k == PV_I64 || k == PV_REAL; }

/* ---- native heads ---------------------------------------------------------- */
static Expr* ll_native(const char* h, Expr* res, bool* handled) {
    *handled = true;
    size_t n = res->data.function.arg_count;
    if (h == SYM_Append)      return ll_append_prepend(res, false);
    if (h == SYM_Prepend)     return ll_append_prepend(res, true);
    if (h == SYM_Reverse)     return ll_reverse(res);
    if (h == SYM_Join)        return ll_join(res);
    if (h == SYM_Take)        return ll_take_drop(res, true);
    if (h == SYM_Drop)        return ll_take_drop(res, false);
    if (h == SYM_Most)        return ll_most_rest(res, true);
    if (h == SYM_Rest)        return ll_most_rest(res, false);
    if (h == SYM_Insert)      return ll_insert(res);
    if (h == SYM_Delete)      return ll_delete(res);
    if (h == SYM_ReplacePart) return ll_replace_part(res);
    if (h == SYM_Map)         return ll_map(res);
    if (h == SYM_Select)      return ll_select(res);
    if (h == SYM_First || h == SYM_Last)
        return ll_first_last_expr(res->data.function.args[0], n == 2 ? res->data.function.args[1] : NULL,
                                  h == SYM_Last);
    *handled = false;
    return NULL;
}

/* PVector: the O(1) / O(log n) heads that need no materialisation. Anything
 * else (or an empty vector, an option argument) is not handled here and the
 * evaluator gate answers it through the equivalent List. */
static Expr* pv_native(const char* h, Expr* res, bool* handled) {
    *handled = false;
    size_t argc = res->data.function.arg_count;
    Expr* a0 = argc ? res->data.function.args[0] : NULL;
    if (!a0 || a0->type != EXPR_PVEC) return NULL;
    PVec* pv = a0->data.pvec;
    int64_t len = pv_length(pv);
    if (h == SYM_Append && argc == 2) {
        *handled = true;
        return expr_new_pvec(pv_push(pv_retain(pv), expr_copy(res->data.function.args[1])));
    }
    if (h == SYM_Prepend && argc == 2) {
        *handled = true;
        PVec* v = pv_push(pv_new(), expr_copy(res->data.function.args[1]));
        for (int64_t i = 0; i < len; i++) v = pv_push(v, pv_get(pv, i));
        return expr_new_pvec(v);
    }
    if (argc == 1 && len > 0) {
        if (h == SYM_First) { *handled = true; return pv_get(pv, 0); }
        if (h == SYM_Last)  { *handled = true; return pv_get(pv, len - 1); }
        if (h == SYM_Most)  { *handled = true; return expr_new_pvec(pv_pop(pv_retain(pv), NULL)); }
        if (h == SYM_Rest) {
            *handled = true;
            PVec* v = pv_new();
            for (int64_t i = 1; i < len; i++) v = pv_push(v, pv_get(pv, i));
            return expr_new_pvec(v);
        }
    }
    return NULL;
}

/* A packed list/array matches _List or _NDArray by its presentation; the
 * matcher keeps its own NDArray head logic, so no head_is is given here. */
static const DsOps ROW_NDARRAY    = { EXPR_NDARRAY,    "NDArray",    nd_head,      NULL,       nd_length, nd_dims, nd_part_row, nd_to_list };
static const DsOps ROW_SPARSE     = { EXPR_SPARSE,     "SparseArray", sp_head,     sp_head_is, sp_length, sp_dims, sp_part_row, sp_to_list };
static const DsOps ROW_LINKEDLIST = { EXPR_LINKEDLIST, "LinkedList", ll_head_sym,  ll_head_is, ll_len,    NULL, ll_part_row, ll_to_list_row, ll_get_row, ll_from_list_row, ll_each_row, ll_atoms_only, ll_native, ll_descend_expr };

static const DsOps ROW_PVEC       = { EXPR_PVEC,       "PVector",    pv_head_row,  pv_head_is, pv_len_row, pv_dims_row, pv_part_row, pv_to_list_row, pv_get_row, pv_from_list_row, pv_each_row, pv_atoms_only, pv_native, NULL };

static void ds_init(void) {
    g_init = true;
    g_rows[ROW_NDARRAY.type]    = &ROW_NDARRAY;
    g_rows[ROW_SPARSE.type]     = &ROW_SPARSE;
    g_rows[ROW_LINKEDLIST.type] = &ROW_LINKEDLIST;
    g_rows[ROW_PVEC.type]       = &ROW_PVEC;
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

Expr* ds_to_list(const Expr* e, bool* handled) {
    const DsOps* o = ds_ops(e->type);
    if (!o || !o->to_list) { *handled = false; return NULL; }
    *handled = true;
    return o->to_list(e);
}

Expr* ds_get(const Expr* e, int64_t i) {
    const DsOps* o = ds_ops(e->type);
    return (o && o->get) ? o->get(e, i) : NULL;
}

Expr* ds_rebuild_like(const Expr* like, const Expr* list) {
    const DsOps* o = ds_ops(like->type);
    return (o && o->from_list) ? o->from_list(list) : NULL;
}

bool ds_has_each(ExprType t) { const DsOps* o = ds_ops(t); return o && o->each; }
bool ds_each(const Expr* e, bool (*fn)(Expr*, void*), void* ud) {
    const DsOps* o = ds_ops(e->type);
    return (o && o->each) ? o->each(e, fn, ud) : false;
}
bool ds_atoms_only(const Expr* e) {
    const DsOps* o = ds_ops(e->type);
    return o && o->atoms_only && o->atoms_only(e);
}

Expr* ds_native(const Expr* kind_val, const char* head, Expr* res, bool* handled) {
    const DsOps* o = ds_ops(kind_val->type);
    if (!o || !o->native) { *handled = false; return NULL; }
    return o->native(head, res, handled);
}

int64_t (*ds_descender(ExprType t))(Expr*, int64_t, int64_t, int64_t, bool, Expr*, bool) {
    const DsOps* o = ds_ops(t);
    return o ? o->descend : NULL;
}
