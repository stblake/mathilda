/* Mathilda — Compile[] engine: lowering the PVector heads.
 *
 * What compiles (elements: machine Integer or Real, one kind per vector):
 *
 *   creating   PVector[], PVector[array]
 *   reading    Length, v[[i]], First, Last, Total, Normal
 *   updating   v = Append[v, x], AppendTo[v, x], v[[i]] = x, v = Most[v]
 *   copying    Most[v], Append[v, x] (a new value; the operand is untouched)
 *
 * where a vector comes from a Module local (`Module[{v = PVector[]}, ...]`) or
 * an argument declared `{v, _PVector}` / `{v, _PVector, _Integer|_Real}`. A
 * vector of any other kind (boxed, mixed, strings), a deeper Part, a sequence
 * spec, a With constant, a global vector: the body bails and the interpreter
 * runs it. A PVector result is handed back as the interpreter's own EXPR_PVEC.
 *
 * VALUE semantics.  A PVector is an immutable value (src/pvec.h): an update
 * yields a new vector and every other holder of the old one keeps it. So this
 * lowering needs none of what compile_llist.c needs for a mutable shared list:
 * no undo frame (a failed run simply drops its registers; an argument is only
 * ever BORROWED and never changed), no aliasing hazard, no global folding.
 *
 * The register holds an EXPR_PVEC node in `.arr`, owned in the array bank (a
 * vector the body made, or a retained reference to an argument), borrowed in an
 * argument register. A plain copy `w = v` retains, so both names share one
 * PVec at refcount 2 and the next update through either path-copies instead of
 * disturbing the other. An update the program makes to a vector it holds ONCE
 * (`v = Append[v, x]`, AppendTo, `v[[i]] = x`, `v = Most[v]`, an operand that
 * is a temporary) is done IN PLACE on that register's PVec -- the exact
 * `v = pv_push(v, x)` loop that makes the interpreter's append O(1).
 *
 * Element type.  A compiled vector has ONE element type for the whole program
 * and every value stored into it must have that type: the interpreter keeps
 * {1, 2.5} as a boxed vector, so storing an Integer into a Real vector cannot
 * be "widened" without changing the answer. A Module local initialised to
 * PVector[] has no type until something is stored, so it is GUESSED -- Real
 * first -- and when a store of the other machine type turns up the store
 * reports it (the ll_want_* fields shared with LinkedList) and the enclosing
 * Module lowers its body again with that type (llist_retry). A second,
 * different type means a mixed vector: the body bails.
 *
 * Evaluation order.  The interpreter evaluates the vector operand BEFORE the
 * other operands. A register read happens at the op, i.e. after them, so an
 * other operand that assigns (v[[(v = Most[v]; 1)]]) could change what is read:
 * such a body bails (has_assign), as the one place the two orders can differ.
 */
#include "compile_emit.h"
#include "../pvec.h"
#include "../symtab.h"
#include "../sym_names.h"
#include "../sym_intern.h"

#include <stdlib.h>
#include <string.h>

static bool pv_elem_ok(CompileType e) { return e == CT_INT || e == CT_REAL; }

static uint16_t pvf(CompileType elem, unsigned extra) {
    return (uint16_t)(((unsigned)elem & 3u) | extra);
}

/* Does `e` contain an assignment-like head anywhere? */
static bool has_assign(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        static const char* names[] = {
            "Set", "SetDelayed", "AddTo", "SubtractFrom", "TimesBy", "DivideBy",
            "Increment", "Decrement", "PreIncrement", "PreDecrement", "AppendTo",
            "PrependTo", "Unset", "Push", "Pop", "PopFirst"
        };
        const char* h = e->data.function.head->data.symbol.name;
        for (size_t i = 0; i < sizeof names / sizeof *names; i++)
            if (h == intern_symbol(names[i])) return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_assign(e->data.function.args[i])) return true;
    return has_assign(e->data.function.head);
}

/* ------------------------------------------------------------------------
 *  Operands
 * ---------------------------------------------------------------------- */

/* Lower `x` as a vector operand. 1: `v` holds it; 0: `x` is not a vector
 * (another lowering may take the head); -1: it is, but cannot be used here.
 * A MUTATING use needs a Module local (an owned register): an argument is the
 * caller's value and is never updated in place. `name` receives the symbol
 * (when `x` is one) for the element-type retry. */
static int pv_operand(Ctx* c, const Expr* x, bool mutate, Val* v, const char** name) {
    if (name) *name = NULL;
    if (x->type == EXPR_SYMBOL) {
        const char* nm = x->data.symbol.name;
        if (name) *name = nm;
        CompileType st; bool sb = false;
        int sr = scope_find(c, nm, &st, &sb);
        if (sr >= 0) {
            if (!CT_IS_PVEC(st)) return 0;
            v->reg = sr; v->tmp = false; v->type = st; v->built = false;
            return 1;
        }
        int k = arg_find(c, nm);
        if (k >= 0) {
            if (!CT_IS_PVEC(c->arg_types[k])) return 0;
            if (mutate) return -1;
            c->argdep[k] = 1;
            v->reg = k; v->tmp = false; v->type = c->arg_types[k]; v->built = false;
            return 1;
        }
        return 0;
    }
    CompileType t;
    if (!infer_type(c, x, &t) || !CT_IS_PVEC(t)) return 0;
    if (mutate) return -1;
    if (!emit(c, x, v)) return -1;
    return CT_IS_PVEC(v->type) ? 1 : -1;
}

static bool pv_operand_type(Ctx* c, const Expr* x, CompileType* t) {
    return infer_type(c, x, t) && CT_IS_PVEC(*t);
}

/* A store of a value of type `vt` into a vector of element type `elem`: equal,
 * or the body bails -- after naming the wanted type, so that a guessed vector
 * can be retyped by the enclosing Module. */
static bool store_type_ok(Ctx* c, const char* name, CompileType elem, CompileType vt) {
    if (vt == elem) return true;
    if (name && pv_elem_ok(vt) && !c->ll_want_name) {
        c->ll_want_name = name;
        c->ll_want_elem = vt;
    }
    c->ok = false;
    return false;
}

/* ------------------------------------------------------------------------
 *  Lowering
 * ---------------------------------------------------------------------- */

/* Checked by pointer BEFORE any operand is typed (both hooks run at every
 * node; typing an operand recurses through the hook again). */
static bool pv_head(const char* h) {
    return h == SYM_PVector || h == SYM_AppendTo || h == SYM_Set || h == SYM_Length
        || h == SYM_Total || h == SYM_Part || h == SYM_First || h == SYM_Last
        || h == SYM_Normal || h == SYM_Most || h == SYM_Append;
}

/* Is every element of the List literal `x` of one machine type? (A literal
 * {1, 2.5} lowers to a Real array where the interpreter keeps a boxed vector.) */
static bool literal_list_uniform(Ctx* c, const Expr* x) {
    if (x->type != EXPR_FUNCTION || x->data.function.head->type != EXPR_SYMBOL
        || x->data.function.head->data.symbol.name != SYM_List) return true;
    CompileType first = CT_ERR;
    for (size_t i = 0; i < x->data.function.arg_count; i++) {
        CompileType t;
        if (!infer_type(c, x->data.function.args[i], &t) || !pv_elem_ok(t)) return false;
        if (i == 0) first = t; else if (t != first) return false;
    }
    return true;
}

static int emit_new(Ctx* c, Expr** A, size_t na, Val* out) {
    Slot z; memset(&z, 0, sizeof z);
    if (na == 0) {
        CompileType e = llist_guess(c, c->ll_init_name);
        if (!pv_elem_ok(e)) return -1;
        int r = alloc_arr(c);
        ins_f(c, OP_PV_NEW, pvf(e, 0), (uint32_t)r, 0, 0, z);
        out->reg = r; out->tmp = true; out->type = CT_PVEC_TYPE(e); out->built = false;
        return c->ok ? 1 : -1;
    }
    if (na != 1 || !literal_list_uniform(c, A[0])) return -1;
    Val av;
    if (!emit(c, A[0], &av)) return -1;
    if (!CT_IS_ARRAY(av.type) || CT_RANK(av.type) != 1 || !pv_elem_ok(CT_ELEM(av.type))) return -1;
    CompileType e = CT_ELEM(av.type);
    unsigned f = 0;
    int r;
    if (av.tmp) { r = av.reg; f |= PVF_FREE_A; }
    else r = alloc_arr(c);
    ins_f(c, OP_PV_FROMARR, pvf(e, f), (uint32_t)r, (uint32_t)av.reg, 0, z);
    out->reg = r; out->tmp = true; out->type = CT_PVEC_TYPE(e); out->built = false;
    return c->ok ? 1 : -1;
}

/* v = Append[v, x] / AppendTo[v, x]: the local's register updated IN PLACE. */
static int emit_push(Ctx* c, const Expr* target, const Expr* valx, Val* out) {
    Val lv, vv;
    const char* nm;
    int k = pv_operand(c, target, true, &lv, &nm);
    if (k <= 0) return -1;
    CompileType e = CT_PVEC_ELEM(lv.type);
    if (has_assign(valx)) return -1;
    if (!emit(c, valx, &vv)) return -1;
    if (!store_type_ok(c, nm, e, vv.type)) return -1;
    Slot z; memset(&z, 0, sizeof z); z.i = PVS_APPEND;
    ins_f(c, OP_PV_STRUCT, pvf(e, PVF_INPLACE), (uint32_t)lv.reg, (uint32_t)lv.reg, (uint32_t)vv.reg, z);
    free_if_tmp(c, vv);
    *out = lv;            /* AppendTo's value is the vector itself (non-tmp) */
    return c->ok ? 1 : -1;
}

static int emit_set_part(Ctx* c, const Expr* part, const Expr* rhs, Val* out) {
    if (part->data.function.arg_count != 2) return 0;
    CompileType lt;
    if (!pv_operand_type(c, part->data.function.args[0], &lt)) return 0;
    Val lv, iv, vv;
    const char* nm;
    int k = pv_operand(c, part->data.function.args[0], true, &lv, &nm);
    if (k <= 0) return -1;
    CompileType e = CT_PVEC_ELEM(lv.type);
    if (has_assign(part->data.function.args[1]) || has_assign(rhs)) return -1;
    if (!emit(c, part->data.function.args[1], &iv)) return -1;
    if (iv.type != CT_INT) return -1;
    if (!emit(c, rhs, &vv)) return -1;
    if (!store_type_ok(c, nm, e, vv.type)) return -1;
    Slot z; memset(&z, 0, sizeof z);
    ins_f(c, OP_PV_SET, pvf(e, 0), (uint32_t)lv.reg, (uint32_t)iv.reg, (uint32_t)vv.reg, z);
    /* Set's value is the right-hand side. */
    if (iv.tmp) {
        ins(c, OP_MOVE, (uint32_t)iv.reg, (uint32_t)vv.reg, 0, z);
        pop_tmp(c, vv);
        *out = iv; out->type = e;
    } else {
        *out = vv;
    }
    return c->ok ? 1 : -1;
}

/* A scalar result of a vector read: the op, then the vector operand released. */
static int pv_scalar_read(Ctx* c, uint16_t op, Val lv, int ib, CompileType rt, Val* out) {
    Slot z; memset(&z, 0, sizeof z);
    int r = alloc_temp(c);
    ins_f(c, op, pvf(CT_PVEC_ELEM(lv.type), 0), (uint32_t)r, (uint32_t)lv.reg, (uint32_t)ib, z);
    free_if_tmp(c, lv);
    out->reg = r; out->tmp = true; out->type = rt; out->built = false;
    return c->ok ? 1 : -1;
}

static int emit_index_read(Ctx* c, const Expr* vx, const Expr* ix, long long fixed, Val* out) {
    Val lv, iv;
    int k = pv_operand(c, vx, false, &lv, NULL);
    if (k <= 0) return k;
    Slot z; memset(&z, 0, sizeof z);
    if (ix) {
        if (has_assign(ix)) return -1;
        if (!emit(c, ix, &iv)) return -1;
        if (iv.type != CT_INT) return -1;
    } else {
        Slot s; memset(&s, 0, sizeof s); s.i = fixed;
        iv = emit_const(c, s, CT_INT);
    }
    pop_tmp(c, iv);
    return pv_scalar_read(c, OP_PV_GET, lv, iv.reg, CT_PVEC_ELEM(lv.type), out);
}

/* Most[v] / Append[v, x]: a NEW value. A temporary operand is consumed and
 * updated in place; a named one is left untouched (the op retains it). */
static int emit_struct(Ctx* c, int sub, Expr** A, size_t na, Val* out) {
    Val lv, bv;
    const char* nm = NULL;
    int k = pv_operand(c, A[0], false, &lv, &nm);
    if (k <= 0) return k;
    CompileType e = CT_PVEC_ELEM(lv.type);
    bool has_b = false;
    if (na == 2) {
        if (!lv.tmp && has_assign(A[1])) return -1;
        if (!emit(c, A[1], &bv)) return -1;
        has_b = true;
        if (!store_type_ok(c, nm, e, bv.type)) return -1;
    }
    unsigned f = 0;
    int dst;
    if (lv.tmp) { f |= PVF_INPLACE; dst = lv.reg; }
    else dst = alloc_arr(c);
    Slot im; memset(&im, 0, sizeof im); im.i = sub;
    ins_f(c, OP_PV_STRUCT, pvf(e, f), (uint32_t)dst, (uint32_t)lv.reg,
          has_b ? (uint32_t)bv.reg : 0u, im);
    if (has_b) free_if_tmp(c, bv);
    out->reg = dst; out->tmp = true; out->type = lv.type; out->built = false;
    return c->ok ? 1 : -1;
}

static int emit_pvec(Ctx* c, const char* h, Expr** A, size_t na, Val* out) {
    if (h == SYM_PVector) return emit_new(c, A, na, out);
    if (na == 0) return 0;

    if (h == SYM_AppendTo && na == 2) {
        CompileType lt;
        if (!pv_operand_type(c, A[0], &lt)) return 0;
        return emit_push(c, A[0], A[1], out);
    }
    if (h == SYM_Set && na == 2 && A[0]->type == EXPR_FUNCTION
        && A[0]->data.function.head->type == EXPR_SYMBOL
        && A[0]->data.function.head->data.symbol.name == SYM_Part)
        return emit_set_part(c, A[0], A[1], out);
    if (h == SYM_Set) return 0;      /* `v = rhs`: the assignment lowering (pvec_emit_assign) */

    /* Everything below reads (or copies) a vector operand. */
    CompileType lt;
    if (!pv_operand_type(c, A[0], &lt)) return 0;
    if (h == SYM_Length && na == 1) {
        Val lv;
        int k = pv_operand(c, A[0], false, &lv, NULL);
        if (k <= 0) return k ? k : -1;
        return pv_scalar_read(c, OP_PV_LEN, lv, 0, CT_INT, out);
    }
    if (h == SYM_Total && na == 1) {
        Val lv;
        int k = pv_operand(c, A[0], false, &lv, NULL);
        if (k <= 0) return k ? k : -1;
        return pv_scalar_read(c, OP_PV_TOTAL, lv, 0, CT_PVEC_ELEM(lv.type), out);
    }
    if (h == SYM_Part && na == 2) return emit_index_read(c, A[0], A[1], 0, out);
    if (h == SYM_First && na == 1) return emit_index_read(c, A[0], NULL, 1, out);
    if (h == SYM_Last && na == 1) return emit_index_read(c, A[0], NULL, -1, out);
    if (h == SYM_Normal && na == 1) {
        Val lv;
        int k = pv_operand(c, A[0], false, &lv, NULL);
        if (k <= 0) return k ? k : -1;
        CompileType e = CT_PVEC_ELEM(lv.type);
        Slot z; memset(&z, 0, sizeof z);
        int r = lv.tmp ? lv.reg : alloc_arr(c);
        ins_f(c, OP_PV_NORMAL, pvf(e, 0), (uint32_t)r, (uint32_t)lv.reg, 0, z);
        /* Normal[v] is a plain List in the interpreter (pv_to_list), packed or
         * not whatever the vector came from: pv_plain makes the boundary present
         * an array result that way (cf_apply), instead of pack_offer's rule. */
        out->reg = r; out->tmp = true; out->type = CT_ARRAY(e, 1); out->built = false;
        c->pv_plain = true;
        return c->ok ? 1 : -1;
    }
    if (h == SYM_Most && na == 1)   return emit_struct(c, PVS_MOST, A, na, out);
    if (h == SYM_Append && na == 2) return emit_struct(c, PVS_APPEND, A, na, out);
    /* Any other head over a vector (Map, Select, Sin[v], Rest, a sequence spec
     * ...) has no compiled form: the body bails. */
    return -1;
}

bool try_emit_pvec(Ctx* c, const char* h, Expr** A, size_t na, Val* out) {
    if (!pv_head(h)) return false;
    int r = emit_pvec(c, h, A, na, out);
    if (r == 0) return false;
    if (r < 0) c->ok = false;
    return true;
}

/* `v = rhs` where v is a PVector local.  Append[v, x] / Most[v] on the SAME
 * local update its register in place (the interpreter's `v = pv_push(v, x)`);
 * any other right-hand side is adopted (a temporary) or retained (a name). */
int pvec_emit_assign(Ctx* c, const Expr* target, const Expr* rhs, Val* out) {
    if (target->type != EXPR_SYMBOL) return 0;
    const char* nm = target->data.symbol.name;
    CompileType vt; bool sb = false;
    int vreg = scope_find(c, nm, &vt, &sb);
    if (vreg < 0 || !CT_IS_PVEC(vt)) return 0;
    Slot z; memset(&z, 0, sizeof z);

    if (rhs->type == EXPR_FUNCTION && rhs->data.function.head->type == EXPR_SYMBOL) {
        const char* rh = rhs->data.function.head->data.symbol.name;
        Expr** RA = rhs->data.function.args;
        size_t rn = rhs->data.function.arg_count;
        bool same = rn >= 1 && RA[0]->type == EXPR_SYMBOL && RA[0]->data.symbol.name == nm;
        if (same && rh == SYM_Append && rn == 2) return emit_push(c, target, RA[1], out);
        if (same && rh == SYM_Most && rn == 1) {
            CompileType e = CT_PVEC_ELEM(vt);
            Slot im; memset(&im, 0, sizeof im); im.i = PVS_MOST;
            ins_f(c, OP_PV_STRUCT, pvf(e, PVF_INPLACE), (uint32_t)vreg, (uint32_t)vreg, 0, im);
            out->reg = vreg; out->tmp = false; out->type = vt; out->built = false;
            return c->ok ? 1 : -1;
        }
    }

    Val rv;
    if (!emit(c, rhs, &rv)) return -1;
    if (!CT_IS_PVEC(rv.type)) { c->ok = false; return -1; }
    if (rv.type != vt) {
        /* the other machine element type: say so, for a retry of the Module */
        if (!c->ll_want_name && pv_elem_ok(CT_PVEC_ELEM(rv.type))) {
            c->ll_want_name = nm;
            c->ll_want_elem = CT_PVEC_ELEM(rv.type);
        }
        c->ok = false; return -1;
    }
    if (rv.tmp && reg_is_owned_arr(rv.reg)) {
        ins(c, OP_A_XFER, (uint32_t)vreg, (uint32_t)rv.reg, 0, z);   /* move the handle */
        pop_tmp(c, rv);
    } else {
        ins(c, OP_PV_RETAIN, (uint32_t)vreg, (uint32_t)rv.reg, 0, z);
    }
    out->reg = vreg; out->tmp = false; out->type = vt; out->built = false;
    return c->ok ? 1 : -1;
}

/* ------------------------------------------------------------------------
 *  Inference (the twin of emit_pvec: same heads, same result types)
 * ---------------------------------------------------------------------- */

bool try_infer_pvec(Ctx* c, const char* h, Expr** A, size_t na, CompileType* out) {
    if (!pv_head(h)) return false;
    CompileType lt, at;
    if (h == SYM_PVector) {
        if (na == 0) {
            CompileType e = llist_guess(c, c->ll_init_name);
            if (!pv_elem_ok(e)) return false;
            *out = CT_PVEC_TYPE(e);
            return true;
        }
        if (na == 1 && literal_list_uniform(c, A[0]) && infer_type(c, A[0], &at)
            && CT_IS_ARRAY(at) && CT_RANK(at) == 1 && pv_elem_ok(CT_ELEM(at))) {
            *out = CT_PVEC_TYPE(CT_ELEM(at));
            return true;
        }
        return false;
    }
    if (na == 0) return false;
    if (h == SYM_Set && na == 2 && A[0]->type == EXPR_FUNCTION
        && A[0]->data.function.head->type == EXPR_SYMBOL
        && A[0]->data.function.head->data.symbol.name == SYM_Part
        && A[0]->data.function.arg_count == 2) {
        if (!pv_operand_type(c, A[0]->data.function.args[0], &lt)) return false;
        *out = CT_PVEC_ELEM(lt);
        return true;
    }
    if (h == SYM_Set) return false;        /* a plain assignment: typed by infer_type itself */
    if (!pv_operand_type(c, A[0], &lt)) return false;
    CompileType e = CT_PVEC_ELEM(lt);
    if (h == SYM_AppendTo && na == 2) { *out = lt; return true; }
    if (h == SYM_Length && na == 1) { *out = CT_INT; return true; }
    if (h == SYM_Total && na == 1) { *out = e; return true; }
    if ((h == SYM_Part && na == 2) || ((h == SYM_First || h == SYM_Last) && na == 1)) { *out = e; return true; }
    if (h == SYM_Normal && na == 1) { *out = CT_ARRAY(e, 1); return true; }
    if ((h == SYM_Most && na == 1) || (h == SYM_Append && na == 2)) { *out = lt; return true; }
    return false;
}
