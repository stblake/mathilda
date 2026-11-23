/* Mathilda — Compile[] engine: lowering the LinkedList heads.
 *
 * What compiles (elements: machine Integer, Real or Complex, one kind per list):
 *
 *   creating   LinkedList[], LinkedList[array]
 *   mutating   AppendTo / PrependTo / Push [l, v], Pop / PopFirst [l],
 *              l[[i]] = v
 *   reading    Length, l[[i]], First, Last, Total, Normal
 *   copying    Most, Rest, Reverse, Take[l, k], Drop[l, k], Append / Prepend
 *              [l, v], Join[l1, l2]
 *
 * where a list comes from a Module local (`Module[{b = LinkedList[]}, ...]`), an
 * argument declared `{b, _LinkedList}` or `{b, _LinkedList, _Integer}`, or -- in
 * auto-compiled code only -- a global whose value is a LinkedList. Anything
 * else (a boxed list, a deeper Part, a sequence spec) bails and the interpreter
 * runs the body.
 *
 * Element type.  A compiled list has ONE element type for the whole program,
 * and it must be the type of every value pushed into it: the interpreter keeps
 * {1, 2.5} as a mixed (boxed) list, so pushing an Integer into a Real list
 * cannot be "coerced" without changing the answer. A Module local initialised
 * to LinkedList[] has no type until something is pushed, so it is GUESSED --
 * Real first -- and when a push of another machine type turns up, the push
 * reports it (ll_want_*) and the enclosing Module lowers its body again with
 * that type (llist_retry). A second, different push type means a mixed list:
 * the body bails, which is the right answer, slowly.
 *
 * Who may mutate a list from outside.  A compiled call that fails half way is
 * abandoned and the interpreter runs the body from the start, so a list the run
 * did not create must look untouched afterwards. compiled_eval runs every
 * program with list ops inside an undo frame for exactly that (linkedlist.c,
 * "Undo log"). That covers one program run once. It cannot cover an
 * auto-compiled PER-ELEMENT program (Table, Map, Sum) whose later element fails
 * after earlier ones succeeded: the interpreter would then redo the earlier
 * elements too. So a folded global list may be MUTATED only under
 * COMPILE_LL_MUTATE_GLOBALS, which only the whole-loop runner (Do / While /
 * For, autocompile_run_closed_stmt) sets; elsewhere a global list is read-only.
 * A CompiledFunction with a list parameter is never inlined into, or called
 * from, another program (compile.c), so it always runs as its own call.
 */
#include "compile_emit.h"
#include "../linkedlist.h"
#include "../symtab.h"
#include "../sym_names.h"

#include <stdlib.h>
#include <string.h>

static bool ll_elem_ok(CompileType e) { return e == CT_INT || e == CT_REAL || e == CT_COMPLEX; }

/* ------------------------------------------------------------------------
 *  Element-type guesses (Module locals initialised to LinkedList[])
 * ---------------------------------------------------------------------- */

CompileType llist_guess(const Ctx* c, const char* name) {
    if (name)
        for (int i = 0; i < c->nll_guess; i++)
            if (c->ll_guess[i].name == name) return c->ll_guess[i].elem;
    return CT_REAL;
}

static bool guess_set(Ctx* c, const char* name, CompileType e) {
    for (int i = 0; i < c->nll_guess; i++)
        if (c->ll_guess[i].name == name) { c->ll_guess[i].elem = e; return true; }
    if (c->nll_guess == (int)(sizeof c->ll_guess / sizeof c->ll_guess[0])) return false;
    c->ll_guess[c->nll_guess].name = name;
    c->ll_guess[c->nll_guess].elem = e;
    c->nll_guess++;
    return true;
}

/* Is `spec` the local declaration `name = LinkedList[]`? */
static bool declares_empty_list(const Expr* spec, const char* name) {
    if (spec->type != EXPR_FUNCTION || spec->data.function.head->type != EXPR_SYMBOL
        || spec->data.function.head->data.symbol.name != SYM_Set
        || spec->data.function.arg_count != 2) return false;
    const Expr* v = spec->data.function.args[0];
    const Expr* init = spec->data.function.args[1];
    return v->type == EXPR_SYMBOL && v->data.symbol.name == name
        && init->type == EXPR_FUNCTION && init->data.function.head->type == EXPR_SYMBOL
        && (init->data.function.head->data.symbol.name == SYM_LinkedList
            || init->data.function.head->data.symbol.name == SYM_PVector)   /* PVector[] shares the retry */
        && init->data.function.arg_count == 0;
}

/* After a failed lowering of a Module with local list `locals`: if the failure
 * was a push of another element type into one of ITS LinkedList[] locals whose
 * guess that type is not yet, adopt it and say "lower the body again". */
bool llist_retry(Ctx* c, const Expr* locals) {
    const char* nm = c->ll_want_name;
    CompileType want = c->ll_want_elem;
    if (!nm || !ll_elem_ok(want) || locals->type != EXPR_FUNCTION) return false;
    bool mine = false;
    for (size_t i = 0; i < locals->data.function.arg_count && !mine; i++)
        mine = declares_empty_list(locals->data.function.args[i], nm);
    if (!mine || llist_guess(c, nm) == want) return false;
    if (!guess_set(c, nm, want)) return false;
    c->ll_want_name = NULL;
    return true;
}

/* ------------------------------------------------------------------------
 *  Folded global lists
 * ---------------------------------------------------------------------- */

static const Expr* global_list(const Ctx* c, const char* nm) {
    if (!(c->flags & COMPILE_FOLD_GLOBALS)) return NULL;
    for (Rule* r = symtab_get_own_values(nm); r; r = r->next) {
        if (!r->pattern || r->pattern->type != EXPR_SYMBOL) continue;
        if (strcmp(r->pattern->data.symbol.name, nm) != 0) continue;
        return (r->replacement && r->replacement->type == EXPR_LINKEDLIST) ? r->replacement : NULL;
    }
    return NULL;
}

bool llist_global_type(const Ctx* c, const char* nm, CompileType* out) {
    const Expr* x = global_list(c, nm);
    if (!x) return false;
    const LinkedList* ll = x->data.llist;
    CompileType e;
    switch (ll_length(ll) ? ll_kind(ll) : LL_EMPTY) {
        case LL_I64:  e = CT_INT; break;
        case LL_REAL: e = CT_REAL; break;
        case LL_CPLX: e = CT_COMPLEX; break;
        case LL_EMPTY:
            e = (c->flags & COMPILE_LL_EMPTY_INT)  ? CT_INT
              : (c->flags & COMPILE_LL_EMPTY_CPLX) ? CT_COMPLEX : CT_REAL;
            break;
        default: return false;                    /* boxed: not a machine list */
    }
    *out = CT_LLIST_TYPE(e);
    return true;
}

/* A register holding a reference to the global list `nm` (LL_REF). */
static bool emit_global_ref(Ctx* c, const char* nm, Val* out) {
    CompileType t;
    const Expr* x = global_list(c, nm);
    if (!x || !llist_global_type(c, nm, &t)) return false;
    if (c->nllrefs == c->llrefs_cap) {
        int nc = c->llrefs_cap ? 2 * c->llrefs_cap : 4;
        Expr** nv = realloc(c->llrefs, (size_t)nc * sizeof *nv);
        if (!nv) { c->ok = false; return false; }
        c->llrefs = nv;
        c->llrefs_cap = nc;
    }
    Expr* ref = expr_copy((Expr*)x);
    c->llrefs[c->nllrefs++] = ref;
    Slot k; memset(&k, 0, sizeof k); k.p = ref;
    int r = alloc_arr(c);
    ins(c, OP_LL_REF, (uint32_t)r, 0, 0, k);
    c->uses_llist = true;
    out->reg = r; out->tmp = true; out->type = t; out->built = false;
    return c->ok;
}

/* A bare global list as a VALUE (bound to a local, returned): it may then be
 * mutated through another name, so only when mutating globals is allowed. */
bool llist_global_value(Ctx* c, const char* nm, Val* out) {
    if (!(c->flags & COMPILE_LL_MUTATE_GLOBALS)) return false;
    return emit_global_ref(c, nm, out);
}

/* ------------------------------------------------------------------------
 *  Operands
 * ---------------------------------------------------------------------- */

/* Lower `x` as a list operand. 1: `v` holds it; 0: `x` is not a list (another
 * lowering may take the head); -1: it is a list but cannot be used here. A
 * MUTATING use of a folded global needs COMPILE_LL_MUTATE_GLOBALS. `name`
 * receives the symbol, when `x` is one, for the element-type retry. */
static int ll_operand(Ctx* c, const Expr* x, bool mutate, Val* v, const char** name) {
    if (name) *name = NULL;
    if (x->type == EXPR_SYMBOL) {
        const char* nm = x->data.symbol.name;
        if (name) *name = nm;
        CompileType st; bool sb = false;
        int sr = scope_find(c, nm, &st, &sb);
        if (sr >= 0) {
            if (!CT_IS_LLIST(st)) return 0;
            v->reg = sr; v->tmp = false; v->type = st; v->built = false;
            return 1;
        }
        int k = arg_find(c, nm);
        if (k >= 0) {
            if (!CT_IS_LLIST(c->arg_types[k])) return 0;
            c->argdep[k] = 1;
            v->reg = k; v->tmp = false; v->type = c->arg_types[k]; v->built = false;
            return 1;
        }
        CompileType gt;
        if (!llist_global_type(c, nm, &gt)) return 0;
        if (mutate && !(c->flags & COMPILE_LL_MUTATE_GLOBALS)) return -1;
        return emit_global_ref(c, nm, v) ? 1 : -1;
    }
    CompileType t;
    if (!infer_type(c, x, &t) || !CT_IS_LLIST(t)) return 0;
    if (!emit(c, x, v)) return -1;
    return CT_IS_LLIST(v->type) ? 1 : -1;
}

/* The type a list operand would have, without lowering it. */
static bool ll_operand_type(Ctx* c, const Expr* x, CompileType* t) {
    return infer_type(c, x, t) && CT_IS_LLIST(*t);
}

/* A push of a value of type `vt` into a list of element type `elem`: equal, or
 * the body bails -- after naming the wanted type, so that a guessed list can be
 * retyped by the enclosing Module. */
static bool push_type_ok(Ctx* c, const char* list_name, CompileType elem, CompileType vt) {
    if (vt == elem) return true;
    if (list_name && ll_elem_ok(vt) && !c->ll_want_name) {
        c->ll_want_name = list_name;
        c->ll_want_elem = vt;
    }
    c->ok = false;
    return false;
}

static uint16_t llf(CompileType elem, unsigned extra) {
    return (uint16_t)(((unsigned)elem & 3u) | extra);
}

/* A scalar result of a list read: the op, then the list operand released. */
static int ll_scalar_read(Ctx* c, uint16_t op, Val lv, int ib, CompileType rt, Val* out) {
    Slot z; memset(&z, 0, sizeof z);
    int r = alloc_temp(c);
    ins_f(c, op, llf(CT_LLIST_ELEM(lv.type), 0), (uint32_t)r, (uint32_t)lv.reg, (uint32_t)ib, z);
    free_if_tmp(c, lv);
    c->uses_llist = true;
    out->reg = r; out->tmp = true; out->type = rt; out->built = false;
    return c->ok ? 1 : -1;
}

/* The destination of an op that yields a NEW handle from list `lv` and an
 * optional second handle `bv` (consumed by the op when temporary). Reusing an
 * operand's own register keeps the array bank LIFO. */
static int handle_dst(Ctx* c, Val lv, const Val* bv, unsigned* flags) {
    if (bv && bv->tmp && CT_IS_HANDLE(bv->type)) *flags |= LLF_FREE_B;
    if (lv.tmp) return lv.reg;
    if (bv && bv->tmp && CT_IS_HANDLE(bv->type)) return bv->reg;
    return alloc_arr(c);
}

static void handle_pop_b(Ctx* c, Val lv, const Val* bv) {
    /* bv's register is the result's when lv was borrowed; otherwise it is
     * above lv's and is released here (the op already freed its handle). */
    if (bv && bv->tmp && CT_IS_HANDLE(bv->type) && lv.tmp) pop_tmp(c, *bv);
}

/* ------------------------------------------------------------------------
 *  Lowering
 * ---------------------------------------------------------------------- */

/* The heads this file can lower. Checked by pointer BEFORE any operand is
 * typed: both hooks run at every node, and typing an operand recurses through
 * the hook again, so an ungated check would cost time exponential in depth. */
static bool ll_head(const char* h) {
    return h == SYM_LinkedList || h == SYM_AppendTo || h == SYM_PrependTo || h == SYM_Push
        || h == SYM_Pop || h == SYM_PopFirst || h == SYM_Set || h == SYM_Length
        || h == SYM_Total || h == SYM_Part || h == SYM_First || h == SYM_Last
        || h == SYM_Normal || h == SYM_Most || h == SYM_Rest || h == SYM_Reverse
        || h == SYM_Take || h == SYM_Drop || h == SYM_Append || h == SYM_Prepend
        || h == SYM_Join;
}

static int emit_new(Ctx* c, Expr** A, size_t na, Val* out) {
    Slot z; memset(&z, 0, sizeof z);
    if (na == 0) {
        CompileType e = llist_guess(c, c->ll_init_name);
        int r = alloc_arr(c);
        ins_f(c, OP_LL_NEW, llf(e, 0), (uint32_t)r, 0, 0, z);
        c->uses_llist = true;
        out->reg = r; out->tmp = true; out->type = CT_LLIST_TYPE(e); out->built = false;
        return c->ok ? 1 : -1;
    }
    if (na != 1) return -1;
    Val av;
    if (!emit(c, A[0], &av)) return -1;
    if (!CT_IS_ARRAY(av.type) || CT_RANK(av.type) != 1 || !ll_elem_ok(CT_ELEM(av.type))) return -1;
    CompileType e = CT_ELEM(av.type);
    unsigned f = 0;
    int r;
    if (av.tmp) { r = av.reg; f |= LLF_FREE_A; }
    else r = alloc_arr(c);
    ins_f(c, OP_LL_FROMARR, llf(e, f), (uint32_t)r, (uint32_t)av.reg, 0, z);
    c->uses_llist = true;
    out->reg = r; out->tmp = true; out->type = CT_LLIST_TYPE(e); out->built = false;
    return c->ok ? 1 : -1;
}

static int emit_push(Ctx* c, Expr** A, bool front, Val* out) {
    Val lv, vv;
    const char* nm;
    int k = ll_operand(c, A[0], true, &lv, &nm);
    if (k <= 0) return k;
    CompileType e = CT_LLIST_ELEM(lv.type);
    if (!emit(c, A[1], &vv)) return -1;
    if (!push_type_ok(c, nm, e, vv.type)) return -1;
    Slot z; memset(&z, 0, sizeof z);
    ins_f(c, OP_LL_PUSH, llf(e, front ? LLF_FRONT : 0), (uint32_t)lv.reg, (uint32_t)vv.reg, 0, z);
    free_if_tmp(c, vv);
    c->uses_llist = true;
    *out = lv;            /* AppendTo's value is the list itself */
    return c->ok ? 1 : -1;
}

static int emit_set_part(Ctx* c, const Expr* part, const Expr* rhs, Val* out) {
    if (part->data.function.arg_count != 2) return 0;
    CompileType lt;
    if (!ll_operand_type(c, part->data.function.args[0], &lt)) return 0;
    Val lv, iv, vv;
    const char* nm;
    int k = ll_operand(c, part->data.function.args[0], true, &lv, &nm);
    if (k <= 0) return k ? k : -1;
    CompileType e = CT_LLIST_ELEM(lv.type);
    if (!emit(c, part->data.function.args[1], &iv)) return -1;
    if (iv.type != CT_INT) return -1;
    if (!emit(c, rhs, &vv)) return -1;
    if (!push_type_ok(c, nm, e, vv.type)) return -1;
    Slot z; memset(&z, 0, sizeof z);
    ins_f(c, OP_LL_SET, llf(e, 0), (uint32_t)lv.reg, (uint32_t)iv.reg, (uint32_t)vv.reg, z);
    c->uses_llist = true;
    /* Set's value is the right-hand side. */
    if (iv.tmp) {
        ins(c, OP_MOVE, (uint32_t)iv.reg, (uint32_t)vv.reg, 0, z);
        pop_tmp(c, vv);
        *out = iv; out->type = e;
    } else {
        *out = vv;
    }
    free_if_tmp(c, lv);
    return c->ok ? 1 : -1;
}

static int emit_index_read(Ctx* c, const Expr* lx, const Expr* ix, long long fixed, Val* out) {
    Val lv, iv;
    int k = ll_operand(c, lx, false, &lv, NULL);
    if (k <= 0) return k;
    Slot z; memset(&z, 0, sizeof z);
    if (ix) {
        if (!emit(c, ix, &iv)) return -1;
        if (iv.type != CT_INT) return -1;
    } else {
        Slot s; memset(&s, 0, sizeof s); s.i = fixed;
        iv = emit_const(c, s, CT_INT);
    }
    pop_tmp(c, iv);
    return ll_scalar_read(c, OP_LL_GET, lv, iv.reg, CT_LLIST_ELEM(lv.type), out);
}

static int emit_struct(Ctx* c, int sub, Expr** A, size_t na, Val* out) {
    Val lv, bv;
    bool has_b = false;
    int k = ll_operand(c, A[0], false, &lv, NULL);
    if (k <= 0) return k;
    CompileType e = CT_LLIST_ELEM(lv.type);
    if (na == 2) {
        if (!emit(c, A[1], &bv)) return -1;
        has_b = true;
        if (sub == LLS_TAKE || sub == LLS_DROP) { if (bv.type != CT_INT) return -1; }
        else if (sub == LLS_JOIN) { if (bv.type != lv.type) return -1; }
        else if (bv.type != e) return -1;   /* Append/Prepend of another type: mixed */
    }
    unsigned f = 0;
    int dst = handle_dst(c, lv, has_b ? &bv : NULL, &f);
    Slot im; memset(&im, 0, sizeof im); im.i = sub;
    ins_f(c, OP_LL_STRUCT, llf(e, f), (uint32_t)dst, (uint32_t)lv.reg,
          has_b ? (uint32_t)bv.reg : 0u, im);
    c->uses_llist = true;
    if (has_b && !CT_IS_HANDLE(bv.type)) pop_tmp(c, bv);
    else if (has_b) handle_pop_b(c, lv, &bv);
    out->reg = dst; out->tmp = true; out->type = lv.type; out->built = false;
    return c->ok ? 1 : -1;
}

static int emit_llist(Ctx* c, const char* h, Expr** A, size_t na, Val* out) {
    if (h == SYM_LinkedList) return emit_new(c, A, na, out);
    if (na == 0) return 0;

    if ((h == SYM_AppendTo || h == SYM_Push) && na == 2) return emit_push(c, A, false, out);
    if (h == SYM_PrependTo && na == 2) return emit_push(c, A, true, out);
    if (h == SYM_Set && na == 2 && A[0]->type == EXPR_FUNCTION
        && A[0]->data.function.head->type == EXPR_SYMBOL
        && A[0]->data.function.head->data.symbol.name == SYM_Part)
        return emit_set_part(c, A[0], A[1], out);

    /* Everything below reads (or copies) a list operand. */
    CompileType lt;
    if (!ll_operand_type(c, A[0], &lt)) {
        /* Pop/PopFirst exist only for lists; let the operand decide the rest. */
        return 0;
    }
    if ((h == SYM_Pop || h == SYM_PopFirst) && na == 1) {
        Val lv;
        int k = ll_operand(c, A[0], true, &lv, NULL);
        if (k <= 0) return k ? k : -1;
        unsigned f = h == SYM_PopFirst ? LLF_FRONT : 0;
        Slot z; memset(&z, 0, sizeof z);
        int r = alloc_temp(c);
        ins_f(c, OP_LL_POP, llf(CT_LLIST_ELEM(lv.type), f), (uint32_t)r, (uint32_t)lv.reg, 0, z);
        free_if_tmp(c, lv);
        c->uses_llist = true;
        out->reg = r; out->tmp = true; out->type = CT_LLIST_ELEM(lv.type); out->built = false;
        return c->ok ? 1 : -1;
    }
    if (h == SYM_Length && na == 1) {
        Val lv;
        int k = ll_operand(c, A[0], false, &lv, NULL);
        if (k <= 0) return k ? k : -1;
        return ll_scalar_read(c, OP_LL_LEN, lv, 0, CT_INT, out);
    }
    if (h == SYM_Total && na == 1) {
        Val lv;
        int k = ll_operand(c, A[0], false, &lv, NULL);
        if (k <= 0) return k ? k : -1;
        return ll_scalar_read(c, OP_LL_TOTAL, lv, 0, CT_LLIST_ELEM(lv.type), out);
    }
    if (h == SYM_Part && na == 2) return emit_index_read(c, A[0], A[1], 0, out);
    if (h == SYM_First && na == 1) return emit_index_read(c, A[0], NULL, 1, out);
    if (h == SYM_Last && na == 1) return emit_index_read(c, A[0], NULL, -1, out);
    if (h == SYM_Normal && na == 1) {
        Val lv;
        int k = ll_operand(c, A[0], false, &lv, NULL);
        if (k <= 0) return k ? k : -1;
        CompileType e = CT_LLIST_ELEM(lv.type);
        Slot z; memset(&z, 0, sizeof z);
        int r = lv.tmp ? lv.reg : alloc_arr(c);
        ins_f(c, OP_LL_NORMAL, llf(e, 0), (uint32_t)r, (uint32_t)lv.reg, 0, z);
        c->uses_llist = true;
        /* Built: the interpreter's Normal[ll] is a List (packed above the size
         * threshold), whatever the list came from. */
        out->reg = r; out->tmp = true; out->type = CT_ARRAY(e, 1); out->built = true;
        return c->ok ? 1 : -1;
    }
    if (h == SYM_Most && na == 1)    return emit_struct(c, LLS_MOST, A, na, out);
    if (h == SYM_Rest && na == 1)    return emit_struct(c, LLS_REST, A, na, out);
    if (h == SYM_Reverse && na == 1) return emit_struct(c, LLS_REVERSE, A, na, out);
    if (h == SYM_Take && na == 2)    return emit_struct(c, LLS_TAKE, A, na, out);
    if (h == SYM_Drop && na == 2)    return emit_struct(c, LLS_DROP, A, na, out);
    if (h == SYM_Append && na == 2)  return emit_struct(c, LLS_APPEND, A, na, out);
    if (h == SYM_Prepend && na == 2) return emit_struct(c, LLS_PREPEND, A, na, out);
    if (h == SYM_Join && na == 2)    return emit_struct(c, LLS_JOIN, A, na, out);
    /* Any other head over a list (Map, Select, Sin[l], a sequence spec ...)
     * has no compiled form: the body bails. */
    return -1;
}

bool try_emit_llist(Ctx* c, const char* h, Expr** A, size_t na, Val* out) {
    if (!ll_head(h)) return false;
    int r = emit_llist(c, h, A, na, out);
    if (r == 0) return false;
    if (r < 0) c->ok = false;
    return true;
}

/* ------------------------------------------------------------------------
 *  Inference (the twin of emit_llist: same heads, same result types)
 * ---------------------------------------------------------------------- */

bool try_infer_llist(Ctx* c, const char* h, Expr** A, size_t na, CompileType* out) {
    if (!ll_head(h)) return false;
    CompileType lt, at;
    if (h == SYM_LinkedList) {
        if (na == 0) { *out = CT_LLIST_TYPE(llist_guess(c, c->ll_init_name)); return true; }
        if (na == 1 && infer_type(c, A[0], &at) && CT_IS_ARRAY(at) && CT_RANK(at) == 1
            && ll_elem_ok(CT_ELEM(at))) { *out = CT_LLIST_TYPE(CT_ELEM(at)); return true; }
        return false;
    }
    if (na == 0) return false;
    if (h == SYM_Set && na == 2 && A[0]->type == EXPR_FUNCTION
        && A[0]->data.function.head->type == EXPR_SYMBOL
        && A[0]->data.function.head->data.symbol.name == SYM_Part
        && A[0]->data.function.arg_count == 2) {
        if (!ll_operand_type(c, A[0]->data.function.args[0], &lt)) return false;
        *out = CT_LLIST_ELEM(lt);
        return true;
    }
    if (!ll_operand_type(c, A[0], &lt)) return false;
    CompileType e = CT_LLIST_ELEM(lt);
    if ((h == SYM_AppendTo || h == SYM_PrependTo || h == SYM_Push) && na == 2) { *out = lt; return true; }
    if ((h == SYM_Pop || h == SYM_PopFirst) && na == 1) { *out = e; return true; }
    if (h == SYM_Length && na == 1) { *out = CT_INT; return true; }
    if (h == SYM_Total && na == 1) { *out = e; return true; }
    if ((h == SYM_Part && na == 2) || ((h == SYM_First || h == SYM_Last) && na == 1)) { *out = e; return true; }
    if (h == SYM_Normal && na == 1) { *out = CT_ARRAY(e, 1); return true; }
    if (((h == SYM_Most || h == SYM_Rest || h == SYM_Reverse) && na == 1)
        || ((h == SYM_Take || h == SYM_Drop || h == SYM_Append || h == SYM_Prepend
             || h == SYM_Join) && na == 2)) { *out = lt; return true; }
    return false;
}
