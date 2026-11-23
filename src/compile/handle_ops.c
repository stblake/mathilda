/* CompileHandleOps registry: the rows and lookup.  See handle_ops.h. */
#include "handle_ops.h"
#include "compile_internal.h"      /* OP_LL_RETAIN / OP_PV_RETAIN */
#include "../assoc.h"              /* is_association / assoc_prebuild_index */
#include "../linkedlist.h"         /* expr_linkedlist / ll_length / ll_kind */
#include "../pvec.h"               /* pv_length / pv_kind */

/* ---- boundary hooks (moved verbatim from compiled_function.c / compile_vm.c) */

static bool array_load_arg(const CompileValue* v, CompileType want, Expr** arr) {
    Expr* x = CT_IS_ARRAY(v->type) ? v->v.a : NULL;
    if (!x || x->type != EXPR_NDARRAY || x->data.ndarray.rank != CT_RANK(want))
        return false;
    *arr = x;
    return true;
}

/* A LinkedList argument is the caller's own list, BORROWED: the program
 * mutates it in place, the interpreter's reference semantics.  Its element
 * kind is checked when the VM binds it (llist_load_arg). */
static bool llist_borrow_arg(const Expr* e, CompileType t) {
    (void)t;
    return e->type == EXPR_LINKEDLIST;
}
/* A PVector argument is an immutable value, BORROWED; a boxed vector, a List or
 * anything else declines.  The element kind is checked at load_arg. */
static bool pvec_borrow_arg(const Expr* e, CompileType t) {
    (void)t;
    return e->type == EXPR_PVEC;
}
/* A read-only Association bag is BORROWED like an NDArray argument.  The
 * single-key index is built HERE, at the single-threaded marshalling boundary,
 * so a parallel VM run never triggers the lazy build (which mutates the shared
 * node) from a worker thread. */
static bool assoc_borrow_arg(const Expr* e, CompileType t) {
    (void)t;
    if (!is_association(e)) return false;
    assoc_prebuild_index(e);
    return true;
}

/* A read-only Association bag rides in `.arr` of a SCALAR-bank argument
 * register (arg registers are < arr_base), so arr_sweep -- which only walks
 * [arr_base, tile_base) -- never frees this borrowed handle. */
static bool assoc_load_arg(const CompileValue* v, CompileType want, Expr** arr) {
    (void)want;
    Expr* x = CT_IS_ASSOC(v->type) ? v->v.a : NULL;
    if (!x || !is_association(x)) return false;
    *arr = x;
    return true;
}
/* BORROWED, in a scalar-bank argument register (never swept), and the caller's
 * own list.  Its storage must already be of the promised element kind -- or
 * empty, which any typed push retypes. */
static bool llist_load_arg(const CompileValue* v, CompileType want, Expr** arr) {
    Expr* x = CT_IS_LLIST(v->type) ? v->v.a : NULL;
    LinkedList* ll = expr_linkedlist(x);
    if (!ll) return false;
    CompileType e = CT_LLIST_ELEM(want);
    LLKind k = e == CT_INT ? LL_I64 : e == CT_REAL ? LL_REAL : LL_CPLX;
    if (ll_length(ll) > 0 && ll_kind(ll) != k) return false;
    *arr = x;
    return true;
}
/* BORROWED (scalar-bank register, never swept) and a value, so the program never
 * changes it: an update starts from a retained reference. */
static bool pvec_load_arg(const CompileValue* v, CompileType want, Expr** arr) {
    Expr* x = CT_IS_PVEC(v->type) ? v->v.a : NULL;
    if (!x || x->type != EXPR_PVEC) return false;
    CompileType e = CT_PVEC_ELEM(want);
    PVKind k = e == CT_INT ? PV_I64 : PV_REAL;
    if (pv_length(x->data.pvec) > 0 && pv_kind(x->data.pvec) != k) return false;
    *arr = x;
    return true;
}

/* ---- the rows ------------------------------------------------------------ */

static const CompileHandleOps g_rows[] = {
    { "Array", CT_ARR, 4 * CT_MAX_RANK, CHF_NUMERIC_OPERAND | CHF_SCOPE_LOCAL,
      0, 0, NULL, NULL, array_load_arg },
    { "Association", CT_ASSOC, 4, 0,
      0, 0, NULL, assoc_borrow_arg, assoc_load_arg },
    { "LinkedList", CT_LLIST, 4, CHF_SCOPE_LOCAL,
      OP_LL_RETAIN, 3, "a LinkedList argument needs a machine element type",
      llist_borrow_arg, llist_load_arg },
    { "PVector", CT_PVEC, 4, CHF_VALUE_SEMANTICS | CHF_SCOPE_LOCAL,
      OP_PV_RETAIN, 2, "a PVector argument needs a machine element type",
      pvec_borrow_arg, pvec_load_arg },
};

const CompileHandleOps* ct_handle_ops(CompileType t) {
    int v = (int)t;
    for (size_t i = 0; i < sizeof g_rows / sizeof g_rows[0]; i++)
        if (v >= (int)g_rows[i].lo && v < (int)g_rows[i].lo + g_rows[i].span) return &g_rows[i];
    return NULL;
}

CompileType ct_handle_elem(CompileType t) {
    const CompileHandleOps* h = ct_handle_ops(t);
    return (CompileType)(((int)t - (int)(h ? h->lo : CT_ARR)) & 3);
}

const char* ct_handle_name(CompileType t) {
    const CompileHandleOps* h = ct_handle_ops(t);
    return h ? h->name : NULL;
}

bool ct_is_opaque_handle(CompileType t) {
    const CompileHandleOps* h = ct_handle_ops(t);
    return h && !(h->flags & CHF_NUMERIC_OPERAND);
}

bool ct_has_value_semantics(CompileType t) {
    const CompileHandleOps* h = ct_handle_ops(t);
    return h && (h->flags & CHF_VALUE_SEMANTICS);
}

bool ct_is_scope_local(CompileType t) {
    const CompileHandleOps* h = ct_handle_ops(t);
    return h && (h->flags & CHF_SCOPE_LOCAL);
}

uint16_t ct_retain_op(CompileType t) {
    const CompileHandleOps* h = ct_handle_ops(t);
    return h ? h->retain_op : 0;
}
