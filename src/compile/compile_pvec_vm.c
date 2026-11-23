/* Mathilda — Compile[] engine: the PVector opcodes at run time.
 *
 * Each op works on the vector through the same pv_* helpers the interpreter
 * uses (src/pvec.h, pvec_ops.c), so a compiled push, read or update is the
 * interpreted one by construction: the same trie, the same copy-on-write rule,
 * no Expr built per element on the typed paths. Total of reals and Normal go
 * through the interpreter's own List answer (see below), bit for bit.
 *
 * A vector register holds an EXPR_PVEC node in `.arr`. The emitter guarantees
 * the storage kind matches the op's element type (PVF_ELEM): a vector the body
 * created only ever receives stores of that type, and an argument was checked
 * at the boundary (load_arg). A vector that still turns out to be of another
 * kind fails the op.
 *
 * VALUE semantics: nothing here is undone on failure. A failing op returns
 * false, the VM abandons the call, the register bank is swept (every owned
 * node released), and the interpreter re-runs the body from the start. An
 * argument is only ever read; every update either works on a node the program
 * owns (in place when its PVec is held once, a path copy when it is shared) or
 * builds a new node from a retained reference.
 */
#include "compile_internal.h"
#include "../pvec.h"
#include "../pvec_ops.h"
#include "../eval.h"
#include "../ndarray.h"
#include "../sym_names.h"

#include <stdint.h>
#include <string.h>

static PVec* reg_pv(const Slot* s) {
    return (s->arr && s->arr->type == EXPR_PVEC) ? s->arr->data.pvec : NULL;
}

static PVKind elem_kind(CompileType e) { return e == CT_INT ? PV_I64 : PV_REAL; }

/* The vector holds elements of kind `e` (an empty one holds any). */
static bool kind_ok(const PVec* pv, CompileType e) {
    return pv_length(pv) == 0 || pv_kind(pv) == elem_kind(e);
}

/* A 1-based Part index (negative from the end) to a 0-based slot, or -1. */
static int64_t resolve_index(const PVec* pv, long long k) {
    int64_t n = pv_length(pv);
    int64_t j = k > 0 ? (int64_t)k - 1 : n + (int64_t)k;
    return (k == 0 || j < 0 || j >= n) ? -1 : j;
}

/* Replace the node in a destination register (a reused register may still hold
 * last call's value). */
static void set_node(Slot* d, Expr* node) {
    expr_free(d->arr);
    d->arr = node;
}

/* The register's node, made PRIVATE to the register.  expr_copy shares a node
 * (a reference count), so `old = v` leaves one EXPR_PVEC node in two registers;
 * an in-place update would then show through the other name. A shared node is
 * replaced by a fresh one over the same PVec (retained, so the PVec's own
 * copy-on-write keeps the other holder's value); an unshared node is returned
 * as is. */
static Expr* own_node(Slot* s) {
    Expr* n = s->arr;
    if (n->refcount > 1) {
        Expr* f = expr_new_pvec(pv_retain(n->data.pvec));
        expr_free(n);                        /* drops this register's share */
        s->arr = f;
        n = f;
    }
    return n;
}

/* The exact sum of an Integer vector, or false past int64 (the interpreter's
 * answer is then a bignum, which no machine register holds). */
static bool sum_i64(const PVec* pv, int64_t* out) {
    int64_t n = pv_length(pv), acc = 0;
    for (int64_t i = 0; i < n; ) {
        const void* p;
        int64_t k = pv_span_at(pv, i, &p);
        const int64_t* w = (const int64_t*)p;
        for (int64_t t = 0; t < k; t++) {
            int64_t x = w[t];
            if ((x > 0 && acc > INT64_MAX - x) || (x < 0 && acc < INT64_MIN - x)) return false;
            acc += x;
        }
        i += k;
    }
    *out = acc;
    return true;
}

bool vm_pvec_op(const Instr* c, Slot* d, Slot* a, Slot* b) {
    const unsigned f = c->flags;
    const CompileType e = PVF_ELEM(f);
    switch (c->op) {
        case OP_PV_NEW:
            set_node(d, expr_new_pvec(pv_new()));
            return true;

        case OP_PV_FROMARR: {
            /* The interpreter's own constructor, over the packed array. */
            Expr* src = a->arr;
            if (!src || src->type != EXPR_NDARRAY || src->data.ndarray.rank != 1
                || src->data.ndarray.dtype != (e == CT_INT ? NDT_INT64 : NDT_FLOAT64)) return false;
            Expr* r = pvec_make(src);
            if (!r) return false;
            if (!kind_ok(r->data.pvec, e)) { expr_free(r); return false; }
            if (f & PVF_FREE_A) { expr_free(a->arr); a->arr = NULL; }
            set_node(d, r);
            return true;
        }

        case OP_PV_RETAIN: {
            if (!reg_pv(a)) return false;
            set_node(d, expr_copy(a->arr));
            return true;
        }

        case OP_PV_STRUCT: {
            PVec* src = reg_pv(a);
            if (!src || !kind_ok(src, e)) return false;
            const bool most = c->imm.i == PVS_MOST;
            if (most && pv_length(src) == 0) return false;            /* Most::nomost */
            if (f & PVF_INPLACE) {
                if (d != a) return false;
                /* The node is this register's own: consume its PVec and put the
                 * update back. Held once: updated in place; shared with another
                 * name or an earlier snapshot: pv_* copies the touched nodes and
                 * leaves the other holder's value alone. */
                Expr* node = own_node(a);
                PVec* p = node->data.pvec;
                if (most) p = pv_pop(p, NULL);
                else if (e == CT_INT) p = pv_push_i64(p, (int64_t)b->i);
                else p = pv_push_real(p, b->r);
                node->data.pvec = p;
                expr_invalidate_hash(node);
                return true;
            }
            PVec* np = pv_retain(src);
            if (most) np = pv_pop(np, NULL);
            else if (e == CT_INT) np = pv_push_i64(np, (int64_t)b->i);
            else np = pv_push_real(np, b->r);
            set_node(d, expr_new_pvec(np));
            return true;
        }

        case OP_PV_SET: {
            PVec* pv = reg_pv(d);
            if (!pv || !kind_ok(pv, e)) return false;
            int64_t j = resolve_index(pv, a->i);
            if (j < 0) return false;                             /* Set::partw */
            Expr* node = own_node(d);
            Expr* x = e == CT_INT ? expr_new_integer((int64_t)b->i) : expr_new_real(b->r);
            node->data.pvec = pv_set(node->data.pvec, j, x);
            expr_invalidate_hash(node);
            return true;
        }

        case OP_PV_GET: {
            const PVec* pv = reg_pv(a);
            if (!pv) return false;
            int64_t j = resolve_index(pv, b->i);
            if (j < 0) return false;                             /* Part::partw, First of {} */
            if (e == CT_INT) { int64_t v; if (!pv_get_i64(pv, j, &v)) return false; d->i = (long long)v; return true; }
            return pv_get_real(pv, j, &d->r);
        }

        case OP_PV_LEN: {
            const PVec* pv = reg_pv(a);
            if (!pv) return false;
            d->i = (long long)pv_length(pv);
            return true;
        }

        case OP_PV_TOTAL: {
            const PVec* pv = reg_pv(a);
            if (!pv) return false;
            if (pv_length(pv) == 0) {
                if (e != CT_INT) return false;       /* Total[{}] is the Integer 0 */
                d->i = 0;
                return true;
            }
            if (!kind_ok(pv, e)) return false;
            if (e == CT_INT) {
                int64_t s;
                if (!sum_i64(pv, &s)) return false;
                d->i = (long long)s;
                return true;
            }
            /* Reals: the interpreter's own Total over the List the vector
             * stands for -- the order of summation is the evaluator's (packed
             * pairwise above its size threshold), so reproduce it exactly
             * rather than re-deriving it. */
            Expr* arg[1] = { pv_to_list(pv) };
            Expr* r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Total), arg, 1));
            bool ok = r && r->type == EXPR_REAL;
            if (ok) d->r = r->data.real;
            if (r) expr_free(r);
            return ok;
        }

        case OP_PV_NORMAL: {
            const PVec* pv = reg_pv(a);
            if (!pv || pv_length(pv) == 0 || !kind_ok(pv, e)) return false;   /* empty: {} */
            Expr* l = pv_to_list(pv);
            Expr* r = ndarray_from_nested_list(l, e == CT_INT ? NDT_INT64 : NDT_FLOAT64);
            expr_free(l);
            if (!r) return false;
            set_node(d, r);
            return true;
        }

        default:
            return false;
    }
}
