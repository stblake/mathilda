/* Mathilda — Compile[]: the Association read-op lowering (Compile[] B1–B5).
 *
 * A declared _Association argument is a read-only parameter bag: Lookup /
 * KeyExistsQ / Length / Values (B1), a runtime-varying key (B2), the pure
 * set-algebra KeyDrop / KeyTake / Counts (B3), the higher-order Map / Select
 * (B4) and the functional Append update (B5) all lower to dedicated opcodes
 * whose immediate borrows a program-owned AssocSpec / AssocCalleeSpec.  This
 * file is the EMIT + inference side (try_emit_assoc / try_infer_assoc and the
 * spec builders); the matching runtime (vm_assoc_*) lives in compile_vm.c and
 * meets this file at the three shared coercions declared in compile_internal.h
 * (assoc_elem_ndt / assoc_value_to_slot / assoc_slot_to_value). */
#include "compile_emit.h"        /* Ctx / Val + the promoted emit-core primitives */
#include "../arithmetic.h"       /* make_complex / expr_new_* */
#include "../symtab.h"           /* symtab_get_own_values, Rule */
#include "../assoc.h"            /* is_association / assoc_lookup_value / assoc_values_list */
#include "../ndarray.h"          /* NDType / NDT_* for assoc_elem_ndt */
#include "../sym_names.h"        /* SYM_Rule / SYM_RuleDelayed */
#include <string.h>
#include <stdlib.h>

/* ================================================================== *
 *  Association read-only parameter bag (Compile[] B1)                 *
 * ================================================================== */

/* NDType for a packed vector of Association values of element type `e`. */
NDType assoc_elem_ndt(CompileType e) {
    return e == CT_COMPLEX ? NDT_COMPLEX64 : e == CT_INT ? NDT_INT64 : NDT_FLOAT64;
}

/* True when `e` contains no argument / local (scope) symbol — i.e. it is a
 * compile-time constant.  Association keys and Lookup defaults must be constant
 * in B1 (a runtime-varying key is B2); this is the gate. */
static bool expr_is_compile_const(const Ctx* c, const Expr* e) {
    if (!e) return true;
    switch (e->type) {
        case EXPR_SYMBOL: {
            const char* nm = e->data.symbol.name;
            CompileType st;
            if (scope_find(c, nm, &st, NULL) >= 0) return false;   /* a local  */
            if (arg_find(c, nm) >= 0)             return false;   /* an arg   */
            return true;                                          /* global   */
        }
        case EXPR_FUNCTION: {
            if (!expr_is_compile_const(c, e->data.function.head)) return false;
            for (size_t i = 0; i < e->data.function.arg_count; i++)
                if (!expr_is_compile_const(c, e->data.function.args[i])) return false;
            return true;
        }
        default: return true;   /* Integer / Real / String / BigInt / ... */
    }
}

/* Common machine element type across a constant association's VALUES, via the
 * same numeric classifier the whole compiler uses (literal).  Fails on an empty,
 * non-canonical, or non-numeric-valued association. */
static bool assoc_const_values_elem(const Expr* assoc, CompileType* out) {
    size_t n = assoc->data.function.arg_count;
    if (n == 0) return false;
    CompileType elem = CT_INT;
    for (size_t i = 0; i < n; i++) {
        const Expr* entry = assoc->data.function.args[i];
        if (entry->type != EXPR_FUNCTION || entry->data.function.arg_count != 2) return false;
        Slot imm; CompileType vt;
        if (!literal(entry->data.function.args[1], &imm, &vt)) return false;
        elem = num_common(elem, vt);
        if ((int)elem < 0 || CT_IS_ARRAY(elem) || CT_IS_ASSOC(elem)) return false;
    }
    *out = elem;
    return true;
}

/* Under COMPILE_FOLD_GLOBALS, a non-argument symbol whose current OwnValue is an
 * association (`g = <|...|>`) — the mirror of global_const for association bags. */
static const Expr* global_assoc(const Ctx* c, const char* nm) {
    if (!(c->flags & COMPILE_FOLD_GLOBALS)) return NULL;
    for (Rule* r = symtab_get_own_values(nm); r; r = r->next) {
        if (!r->pattern || r->pattern->type != EXPR_SYMBOL) continue;
        if (strcmp(r->pattern->data.symbol.name, nm) != 0) continue;
        return (r->replacement && is_association(r->replacement)) ? r->replacement : NULL;
    }
    return NULL;
}

typedef enum { ASSOC_NONE, ASSOC_ARG, ASSOC_CONST, ASSOC_EXPR } AssocOpKind;
typedef struct {
    AssocOpKind kind;
    int         reg;      /* ASSOC_ARG: the operand register             */
    CompileType valtype;  /* ASSOC_ARG: declared value element type      */
    const Expr* assoc;    /* ASSOC_CONST: the constant association node   */
    const Expr* expr;     /* ASSOC_EXPR: a produced-association subexpr   */
} AssocOperand;

/* A head that PRODUCES an association value from an association-typed body, so a
 * `KeyDrop[p, k]` / `Counts[v]` may itself be the operand of another association
 * op (composition). */
static bool assoc_producer_head(const char* h, size_t na) {
    if ((strcmp(h, "KeyDrop") == 0 || strcmp(h, "KeyTake") == 0) && na == 2) return true;
    if (strcmp(h, "Counts") == 0 && na == 1) return true;
    return false;
}

/* THE shared resolver — used by both infer_type and emit_node so they can never
 * disagree about what an association operand is.  A declared `_Association`
 * argument -> ASSOC_ARG (register + value type); a literal `<|...|>` or a folded
 * global -> ASSOC_CONST; a produced-association subexpression (KeyDrop/Counts/…)
 * -> ASSOC_EXPR (emitted/inferred by the caller); anything else -> ASSOC_NONE
 * (the head bails, or, for a head shared with the array lowerings, falls through
 * to them). */
static void resolve_assoc_operand(Ctx* c, const Expr* e, AssocOperand* o) {
    o->kind = ASSOC_NONE; o->reg = -1; o->valtype = CT_REAL; o->assoc = NULL; o->expr = NULL;
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        const char* nm = e->data.symbol.name;
        CompileType st;
        if (scope_find(c, nm, &st, NULL) >= 0) return;   /* a local is not a bag */
        int k = arg_find(c, nm);
        if (k >= 0) {
            if (CT_IS_ASSOC(c->arg_types[k])) {
                o->kind = ASSOC_ARG; o->reg = k;
                o->valtype = CT_ASSOC_VALTYPE(c->arg_types[k]);
                c->argdep[k] = 1;
            }
            return;   /* an arg of another type is not an association bag */
        }
        const Expr* g = global_assoc(c, nm);
        if (g) { o->kind = ASSOC_CONST; o->assoc = g; }
        return;
    }
    if (is_association(e)) { o->kind = ASSOC_CONST; o->assoc = e; return; }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && assoc_producer_head(e->data.function.head->data.symbol.name, e->data.function.arg_count)) {
        o->kind = ASSOC_EXPR; o->expr = e;
    }
}

/* Materialise an association SOURCE operand for a consumer: a borrowed argument
 * register, a constant node (carried in the spec), or a freshly EMITTED produced
 * association (an owned array-bank temp).  Returns false (bailing via c->ok) if
 * `e` is not an association.  `*owned` is true only for a produced temp — the
 * consumer must free it (an array-producer via its free-source flag with LIFO
 * slot reuse; a scalar reader via free_if_tmp after the read). */
typedef struct { int reg; bool owned; const Expr* cst; CompileType valtype; } AssocSrc;
static bool materialize_assoc_src(Ctx* c, const AssocOperand* ao, AssocSrc* s) {
    s->reg = -1; s->owned = false; s->cst = NULL; s->valtype = CT_REAL;
    if (ao->kind == ASSOC_ARG)  { s->reg = ao->reg; s->valtype = ao->valtype; return true; }
    if (ao->kind == ASSOC_CONST) {
        s->cst = ao->assoc;
        if (!assoc_const_values_elem(ao->assoc, &s->valtype)) s->valtype = CT_REAL;
        return true;
    }
    if (ao->kind == ASSOC_EXPR) {
        Val av;
        if (!emit(c, ao->expr, &av) || !CT_IS_ASSOC(av.type)) { c->ok = false; return false; }
        s->reg = av.reg; s->owned = av.tmp; s->valtype = CT_ASSOC_VALTYPE(av.type);
        return true;
    }
    return false;
}

/* Program-owned AssocSpec, freed in compiled_free (declared in compile_internal.h). */
void compile_assocspec_free(AssocSpec* p) {
    if (!p) return;
    expr_free(p->key); expr_free(p->deflt); expr_free(p->assoc);
    free(p);
}

static bool expr_eq_or_null(const Expr* a, const Expr* b) {
    if (!a && !b) return true;
    return a && b && expr_eq(a, b);
}
static bool assocspec_eq(const AssocSpec* a, const AssocSpec* b) {
    return expr_eq_or_null(a->key, b->key)
        && expr_eq_or_null(a->deflt, b->deflt)
        && expr_eq_or_null(a->assoc, b->assoc);
}

static bool ctx_own_assocspec(Ctx* c, AssocSpec* p) {
    if (c->nassocs == c->assocs_cap) {
        int nc = c->assocs_cap ? c->assocs_cap * 2 : 4;
        AssocSpec** np = realloc(c->assocs, (size_t)nc * sizeof *np);
        if (!np) { compile_assocspec_free(p); c->ok = false; return false; }
        c->assocs = np; c->assocs_cap = nc;
    }
    c->assocs[c->nassocs++] = p;
    return true;
}

/* Build (or reuse) a program-owned AssocSpec.  The three Exprs are ref-counted
 * keep-alives (expr_copy) so the spec outlives the body; identical specs are
 * interned so two identical lookups share one pointer and thus CSE (K_KERN1's
 * imm_eq compares imm.p).  Returns NULL and sets c->ok=false on failure. */
static AssocSpec* emit_assocspec(Ctx* c, const Expr* key, const Expr* deflt, const Expr* assoc) {
    AssocSpec* sp = calloc(1, sizeof *sp);
    if (!sp) { c->ok = false; return NULL; }
    if (key)   sp->key   = expr_copy((Expr*)key);
    if (deflt) sp->deflt = expr_copy((Expr*)deflt);
    if (assoc) sp->assoc = expr_copy((Expr*)assoc);
    if ((key && !sp->key) || (deflt && !sp->deflt) || (assoc && !sp->assoc)) {
        compile_assocspec_free(sp); c->ok = false; return NULL;
    }
    for (int i = 0; i < c->nassocs; i++)
        if (assocspec_eq(c->assocs[i], sp)) { compile_assocspec_free(sp); return c->assocs[i]; }
    if (!ctx_own_assocspec(c, sp)) return NULL;   /* frees sp, sets c->ok */
    return sp;
}

/* --- B4 higher-order transform callee pool ------------------------------- */
void compile_assoccallee_free(AssocCalleeSpec* p) {
    if (!p) return;
    compiled_free(p->callee);   /* the callee is a self-contained program */
    expr_free(p->assoc);
    free(p);
}
static bool ctx_own_assoccallee(Ctx* c, AssocCalleeSpec* p) {
    if (c->ncallees == c->callees_cap) {
        int nc = c->callees_cap ? c->callees_cap * 2 : 4;
        AssocCalleeSpec** np = realloc(c->callees, (size_t)nc * sizeof *np);
        if (!np) { compile_assoccallee_free(p); c->ok = false; return false; }
        c->callees = np; c->callees_cap = nc;
    }
    c->callees[c->ncallees++] = p;
    return true;
}

/* B4 higher-order callee compiler: turns an embedded one-arg function into a
 * standalone callee program.  Defined below fn_resolve (which it needs);
 * forward-declared here for try_emit_assoc. */
struct CompiledProgram* compile_value_callee(Ctx* c, const Expr* fn, CompileType in);

/* Coerce an association value Expr to a machine Slot of type `t` (the same
 * coercions cf_box uses for a scalar argument).  False if it does not fit. */
bool assoc_value_to_slot(const Expr* v, CompileType t, Slot* out) {
    switch (t) {
        case CT_INT:     return cf_to_ll(v, &out->i);
        case CT_REAL:    return cf_to_double(v, &out->r);
        case CT_COMPLEX: { double re, im; if (!cf_to_complex(v, &re, &im)) return false;
                           out->z = re + im * I; return true; }
        default: return false;
    }
}
/* Box a machine Slot result of type `t` back to an Expr (mirrors cf_unbox). */
Expr* assoc_slot_to_value(Slot s, CompileType t) {
    switch (t) {
        case CT_INT:  return expr_new_integer((int64_t)s.i);
        case CT_REAL: return expr_new_real(s.r);
        case CT_COMPLEX: {
            double re = creal(s.z), im = cimag(s.z);
            return im == 0.0 ? expr_new_real(re) : make_complex(expr_new_real(re), expr_new_real(im));
        }
        default: return NULL;
    }
}

/* Emit the Association read ops (B1).  Returns true if `h` is an association op
 * that was HANDLED — lowered (out set) or cleanly bailed (c->ok=false).  Returns
 * false if `h` is not one, or is a head shared with the array lowerings
 * (Length/Values) whose operand is not an association — so the caller keeps
 * looking. */
bool try_emit_assoc(Ctx* c, const char* h, Expr** A, size_t na, Val* out) {
    bool is_lookup  = strcmp(h, "Lookup") == 0 && (na == 2 || na == 3);
    bool is_exists  = (strcmp(h, "KeyExistsQ") == 0 || strcmp(h, "KeyMemberQ") == 0) && na == 2;
    bool is_free    = strcmp(h, "KeyFreeQ") == 0 && na == 2;
    bool is_len     = strcmp(h, "Length") == 0 && na == 1;
    bool is_values  = strcmp(h, "Values") == 0 && na == 1;
    bool is_keydrop = strcmp(h, "KeyDrop") == 0 && na == 2;
    bool is_keytake = strcmp(h, "KeyTake") == 0 && na == 2;
    bool is_counts  = strcmp(h, "Counts") == 0 && na == 1;
    bool is_map     = strcmp(h, "Map") == 0 && na == 2;      /* Map[f, assoc]     */
    bool is_select  = strcmp(h, "Select") == 0 && na == 2;   /* Select[assoc, p]  */
    bool is_append  = strcmp(h, "Append") == 0 && na == 2;   /* Append[assoc,k->v]*/
    if (!is_lookup && !is_exists && !is_free && !is_len && !is_values
        && !is_keydrop && !is_keytake && !is_counts && !is_map && !is_select && !is_append) return false;

    /* B4: Map[f, assoc] / Select[assoc, pred] — a higher-order transform via a
     * compiled callee (the embedded function, one machine value in).  Map's
     * association operand is A[1], Select's is A[0].  Handled before the
     * association-operand resolver; if the operand is NOT an association we
     * return false so the ordinary array Map/Select lowering takes it. */
    if (is_map || is_select) {
        const Expr* fn_expr    = is_map ? A[0] : A[1];
        const Expr* assoc_expr = is_map ? A[1] : A[0];
        AssocOperand mo; resolve_assoc_operand(c, assoc_expr, &mo);
        if (mo.kind == ASSOC_NONE) return false;     /* not an assoc -> array path */
        AssocSrc s;
        if (!materialize_assoc_src(c, &mo, &s)) { c->ok = false; return true; }
        struct CompiledProgram* callee = compile_value_callee(c, fn_expr, s.valtype);
        if (!callee) { c->ok = false; return true; }
        CompileType rtype = compiled_result_type(callee);
        bool ok_result = is_map ? (rtype == CT_INT || rtype == CT_REAL || rtype == CT_COMPLEX)
                                : (rtype == CT_BOOL);
        if (!ok_result || compiled_num_args(callee) != 1) {
            compiled_free(callee); c->ok = false; return true;
        }
        AssocCalleeSpec* sp = calloc(1, sizeof *sp);
        if (!sp) { compiled_free(callee); c->ok = false; return true; }
        sp->callee = callee;
        if (s.cst) {
            sp->assoc = expr_copy((Expr*)s.cst);
            if (!sp->assoc) { compile_assoccallee_free(sp); c->ok = false; return true; }
        }
        if (!ctx_own_assoccallee(c, sp)) return true;    /* frees sp+callee, sets c->ok */
        CompileType out_valtype = is_map ? rtype : s.valtype;   /* Select keeps values */
        uint16_t f = (uint16_t)((unsigned)s.valtype | ((unsigned)rtype << 4));
        uint32_t areg;
        if (s.cst) areg = 0;
        else if (s.owned) {
            Val sv = { s.reg, true, CT_ASSOC_TYPE(s.valtype), false };
            pop_tmp(c, sv); f |= 0x100u; areg = (uint32_t)s.reg;
        } else areg = (uint32_t)s.reg;
        Slot ip; memset(&ip, 0, sizeof ip); ip.p = sp;
        int dst = alloc_arr(c);
        ins_f(c, is_map ? OP_ASSOC_MAP : OP_ASSOC_SELECT, f, (uint32_t)dst, areg, 0, ip);
        out->reg = dst; out->tmp = true; out->type = CT_ASSOC_TYPE(out_valtype); out->built = true;
        return true;
    }

    /* Counts is the one association op whose OPERAND is a machine array, not an
     * association -> <|element -> count|> (integer values).  Handled before the
     * association-operand resolver below. */
    if (is_counts) {
        Val av;
        if (!emit(c, A[0], &av)) { c->ok = false; return true; }
        if (!CT_IS_ARRAY(av.type) || CT_RANK(av.type) != 1) { c->ok = false; return true; }
        uint16_t f = 0;
        uint32_t areg = (uint32_t)av.reg;
        if (av.tmp) { pop_tmp(c, av); f |= 1u; }   /* free the produced array temp */
        Slot z; memset(&z, 0, sizeof z);
        int dst = alloc_arr(c);
        ins_f(c, OP_ASSOC_COUNTS, f, (uint32_t)dst, areg, 0, z);
        out->reg = dst; out->tmp = true; out->type = CT_ASSOC_TYPE(CT_INT); out->built = true;
        return true;
    }

    AssocOperand ao; resolve_assoc_operand(c, A[0], &ao);
    if (ao.kind == ASSOC_NONE) {
        /* Lookup / KeyExistsQ / KeyMemberQ / KeyFreeQ / KeyDrop / KeyTake are
         * association-only: bail.  Length / Values / Append are shared with the
         * array path — let it try. */
        if (is_len || is_values || is_append) return false;
        c->ok = false; return true;
    }

    /* Unified source resolution.  A produced-association operand (ASSOC_EXPR:
     * `KeyDrop[p,k]`, `Counts[v]`, …) is EMITTED here into an owned array-bank
     * temp; an argument bag stays a borrowed register; a constant rides in the
     * spec.  From here every consumer works off (src_reg | src_cst) + src_owned. */
    AssocSrc src;
    if (!materialize_assoc_src(c, &ao, &src)) { c->ok = false; return true; }
    int         src_reg   = src.reg;      /* ARG / produced register, else -1     */
    const Expr* src_cst   = src.cst;      /* constant source node, else NULL      */
    bool        src_owned = src.owned;    /* produced temp — this op must free it  */
    CompileType src_valtype = src.valtype;

    /* B3: KeyDrop / KeyTake -> an OWNED association (array bank).  An owned
     * produced source is consumed in place: pop it so the result reuses its slot
     * and set the free-source flag (bit1) so the VM frees R[a] after reading it. */
    if (is_keydrop || is_keytake) {
        const Expr* keyspec = A[1];
        if (!expr_is_compile_const(c, keyspec)) { c->ok = false; return true; }
        AssocSpec* sp = emit_assocspec(c, keyspec, NULL, src_cst);
        if (!sp) return true;
        Slot ip; memset(&ip, 0, sizeof ip); ip.p = sp;
        uint16_t f = (uint16_t)(is_keytake ? 1u : 0u);   /* bit0 = take */
        uint32_t areg;
        if (src_cst) areg = 0;                            /* source in the spec   */
        else if (src_owned) {
            Val sv = { src_reg, true, CT_ASSOC_TYPE(src_valtype), false };
            pop_tmp(c, sv); f |= 2u; areg = (uint32_t)src_reg;   /* free-source, reuse slot */
        } else areg = (uint32_t)src_reg;                  /* borrowed argument bag */
        int dst = alloc_arr(c);
        ins_f(c, OP_ASSOC_KEYSEL, f, (uint32_t)dst, areg, 0, ip);
        out->reg = dst; out->tmp = true; out->type = CT_ASSOC_TYPE(src_valtype); out->built = true;
        return true;
    }

    /* B5: Append[assoc, key -> value] -> a new OWNED association with the key set
     * (functional insert/replace).  The key is a compile-time constant; the value
     * is a runtime machine scalar coerced to the bag's value type (so the result
     * keeps that type).  The mutating AssociateTo[sym, …] is deliberately NOT
     * compiled — rebinding a symbol's OwnValue is a side effect the register VM
     * does not model — and stays in the interpreter. */
    if (is_append) {
        const Expr* rule = A[1];
        if (!(rule->type == EXPR_FUNCTION && rule->data.function.head->type == EXPR_SYMBOL
              && rule->data.function.head->data.symbol.name == SYM_Rule
              && rule->data.function.arg_count == 2)) { c->ok = false; return true; }
        const Expr* key = rule->data.function.args[0];
        if (!expr_is_compile_const(c, key)) { c->ok = false; return true; }
        Val vv;
        if (!emit(c, rule->data.function.args[1], &vv)) { c->ok = false; return true; }
        coerce(c, &vv, src_valtype);   /* value keeps the bag's value type */
        if (!c->ok) return true;
        if (src_valtype != CT_INT && src_valtype != CT_REAL && src_valtype != CT_COMPLEX) {
            pop_tmp(c, vv); c->ok = false; return true;
        }
        AssocSpec* sp = emit_assocspec(c, key, NULL, src_cst);
        if (!sp) { pop_tmp(c, vv); return true; }
        Slot ip; memset(&ip, 0, sizeof ip); ip.p = sp;
        uint16_t f = (uint16_t)((unsigned)src_valtype);
        pop_tmp(c, vv);
        uint32_t areg;
        if (src_cst) areg = 0;
        else if (src_owned) {
            Val sv = { src_reg, true, CT_ASSOC_TYPE(src_valtype), false };
            pop_tmp(c, sv); f |= 0x100u; areg = (uint32_t)src_reg;
        } else areg = (uint32_t)src_reg;
        int dst = alloc_arr(c);
        ins_f(c, OP_ASSOC_SET, f, (uint32_t)dst, areg, (uint32_t)vv.reg, ip);
        out->reg = dst; out->tmp = true; out->type = CT_ASSOC_TYPE(src_valtype); out->built = true;
        return true;
    }

    if (is_lookup) {
        const Expr* key   = A[1];
        const Expr* deflt = (na == 3) ? A[2] : NULL;
        if (deflt && !expr_is_compile_const(c, deflt)) { c->ok = false; return true; }

        if (!expr_is_compile_const(c, key)) {
            /* B2: a runtime-varying INTEGER or REAL key.  Emit the key as a
             * machine scalar and probe Part A's index at runtime (O(1)); the op
             * is pure (K_KERN2) so it CSEs by (bag, key, spec) and is hoisted
             * only when the key register is itself loop-invariant. */
            Val kv;
            if (!emit(c, key, &kv)) { c->ok = false; return true; }
            CompileType kt = kv.type;
            if (kt != CT_INT && kt != CT_REAL) { pop_tmp(c, kv); c->ok = false; return true; }
            CompileType rt = src_valtype;
            if (src_cst && !assoc_const_values_elem(src_cst, &rt)) { pop_tmp(c, kv); c->ok = false; return true; }
            if (deflt) {
                Slot d; CompileType dt;
                if (!literal(deflt, &d, &dt)) { pop_tmp(c, kv); c->ok = false; return true; }
                rt = num_common(rt, dt);
                if ((int)rt < 0 || CT_IS_ARRAY(rt) || CT_IS_ASSOC(rt)) { pop_tmp(c, kv); c->ok = false; return true; }
            }
            if (rt != CT_INT && rt != CT_REAL && rt != CT_COMPLEX) { pop_tmp(c, kv); c->ok = false; return true; }
            AssocSpec* sp = emit_assocspec(c, NULL, deflt, src_cst);
            if (!sp) { pop_tmp(c, kv); return true; }
            Slot ip; memset(&ip, 0, sizeof ip); ip.p = sp;
            uint16_t f = (uint16_t)((unsigned)rt | ((unsigned)kt << 4));
            /* For a register bag `a` is that register; for a constant bag the
             * handle is in the spec and `a` is a harmless read of the key. */
            uint32_t bagreg = src_cst ? (uint32_t)kv.reg : (uint32_t)src_reg;
            pop_tmp(c, kv);
            int dst = alloc_temp(c);
            ins_f(c, OP_ASSOC_LOOKUP_DYN, f, (uint32_t)dst, bagreg, (uint32_t)kv.reg, ip);
            if (src_owned) free_if_tmp(c, (Val){ src_reg, true, CT_ASSOC_TYPE(src_valtype), false });
            out->reg = dst; out->tmp = true; out->type = rt; out->built = false;
            return true;
        }
        if (src_cst) {
            /* Fold at compile time via the native helper (never evaluate()). */
            Expr* v = assoc_lookup_value(src_cst, key);
            if (!v) v = (Expr*)deflt;
            Slot imm; CompileType vt;
            if (!v || !literal(v, &imm, &vt)) { c->ok = false; return true; }
            *out = emit_const(c, imm, vt);
            return true;
        }
        CompileType rt = src_valtype;
        if (deflt) {
            Slot d; CompileType dt;
            if (!literal(deflt, &d, &dt)) { c->ok = false; return true; }
            rt = num_common(rt, dt);
            if ((int)rt < 0 || CT_IS_ARRAY(rt) || CT_IS_ASSOC(rt)) { c->ok = false; return true; }
        }
        if (rt != CT_INT && rt != CT_REAL && rt != CT_COMPLEX) { c->ok = false; return true; }
        AssocSpec* sp = emit_assocspec(c, key, deflt, NULL);
        if (!sp) return true;   /* c->ok already false */
        Slot ip; memset(&ip, 0, sizeof ip); ip.p = sp;
        int dst = alloc_temp(c);
        ins_f(c, OP_ASSOC_LOOKUP, (uint16_t)rt, (uint32_t)dst, (uint32_t)src_reg, 0, ip);
        if (src_owned) free_if_tmp(c, (Val){ src_reg, true, CT_ASSOC_TYPE(src_valtype), false });
        out->reg = dst; out->tmp = true; out->type = rt; out->built = false;
        return true;
    }

    if (is_exists || is_free) {
        const Expr* key = A[1];
        if (!expr_is_compile_const(c, key)) { c->ok = false; return true; }
        if (src_cst) {
            int present = assoc_lookup_value(src_cst, key) != NULL;
            Slot imm; imm.i = is_free ? !present : present;
            *out = emit_const(c, imm, CT_BOOL);
            return true;
        }
        AssocSpec* sp = emit_assocspec(c, key, NULL, NULL);
        if (!sp) return true;
        Slot ip; memset(&ip, 0, sizeof ip); ip.p = sp;
        int dst = alloc_temp(c);
        ins_f(c, OP_ASSOC_HASKEY, (uint16_t)(is_free ? 1 : 0), (uint32_t)dst, (uint32_t)src_reg, 0, ip);
        if (src_owned) free_if_tmp(c, (Val){ src_reg, true, CT_ASSOC_TYPE(src_valtype), false });
        out->reg = dst; out->tmp = true; out->type = CT_BOOL; out->built = false;
        return true;
    }

    if (is_len) {
        if (src_cst) {
            Slot imm; imm.i = (long long)src_cst->data.function.arg_count;
            *out = emit_const(c, imm, CT_INT);
            return true;
        }
        Slot z; memset(&z, 0, sizeof z);
        int dst = alloc_temp(c);
        ins(c, OP_ASSOC_LEN, (uint32_t)dst, (uint32_t)src_reg, 0, z);
        if (src_owned) free_if_tmp(c, (Val){ src_reg, true, CT_ASSOC_TYPE(src_valtype), false });
        out->reg = dst; out->tmp = true; out->type = CT_INT; out->built = false;
        return true;
    }

    /* Values: an owned packed vector.  A constant source rides in the spec; a
     * produced source is consumed in place (pop + reuse slot + free-source flag,
     * bit 0x100), an argument bag stays borrowed. */
    {
        CompileType elem = src_valtype;
        if (src_cst && !assoc_const_values_elem(src_cst, &elem)) { c->ok = false; return true; }
        if (elem != CT_INT && elem != CT_REAL && elem != CT_COMPLEX) { c->ok = false; return true; }
        AssocSpec* sp = emit_assocspec(c, NULL, NULL, src_cst);
        if (!sp) return true;
        Slot ip; memset(&ip, 0, sizeof ip); ip.p = sp;
        uint16_t f = (uint16_t)elem;
        uint32_t areg;
        if (src_cst) areg = 0;
        else if (src_owned) {
            Val sv = { src_reg, true, CT_ASSOC_TYPE(src_valtype), false };
            pop_tmp(c, sv); f |= 0x100u; areg = (uint32_t)src_reg;   /* free-source, reuse slot */
        } else areg = (uint32_t)src_reg;
        int dst = alloc_arr(c);
        ins_f(c, OP_ASSOC_VALUES, f, (uint32_t)dst, areg, 0, ip);
        out->reg = dst; out->tmp = true; out->type = CT_ARRAY(elem, 1); out->built = true;
        return true;
    }
}

/* infer_type twin of try_emit_assoc: returns true and sets *out when `h` is an
 * association op with an association operand whose result type is known; false
 * otherwise (so a shared head falls through to the array branch, and an
 * association-only head with no valid operand ultimately bails). */
bool try_infer_assoc(Ctx* c, const char* h, Expr** A, size_t na, CompileType* out) {
    bool is_lookup = strcmp(h, "Lookup") == 0 && (na == 2 || na == 3);
    bool is_exists = (strcmp(h, "KeyExistsQ") == 0 || strcmp(h, "KeyMemberQ") == 0
                      || strcmp(h, "KeyFreeQ") == 0) && na == 2;
    bool is_len    = strcmp(h, "Length") == 0 && na == 1;
    bool is_values = strcmp(h, "Values") == 0 && na == 1;
    bool is_keysel = (strcmp(h, "KeyDrop") == 0 || strcmp(h, "KeyTake") == 0) && na == 2;
    bool is_counts = strcmp(h, "Counts") == 0 && na == 1;
    if (!is_lookup && !is_exists && !is_len && !is_values && !is_keysel && !is_counts) return false;

    if (is_counts) {
        CompileType t;
        if (!infer_type(c, A[0], &t) || !CT_IS_ARRAY(t) || CT_RANK(t) != 1) return false;
        *out = CT_ASSOC_TYPE(CT_INT); return true;
    }

    AssocOperand ao; resolve_assoc_operand(c, A[0], &ao);
    if (ao.kind == ASSOC_NONE) return false;

    /* Unified source resolution (type-only twin of materialize_assoc_src). */
    CompileType src_valtype = CT_REAL;
    const Expr* src_cst = NULL;
    if (ao.kind == ASSOC_ARG) src_valtype = ao.valtype;
    else if (ao.kind == ASSOC_CONST) src_cst = ao.assoc;   /* valtype computed per-op */
    else { /* ASSOC_EXPR */
        CompileType t;
        if (!infer_type(c, ao.expr, &t) || !CT_IS_ASSOC(t)) return false;
        src_valtype = CT_ASSOC_VALTYPE(t);
    }

    if (is_exists) { *out = CT_BOOL; return true; }
    if (is_len)    { *out = CT_INT;  return true; }
    if (is_keysel) {
        if (!expr_is_compile_const(c, A[1])) return false;
        CompileType vt = src_valtype;
        if (src_cst && !assoc_const_values_elem(src_cst, &vt)) vt = CT_REAL;
        *out = CT_ASSOC_TYPE(vt); return true;
    }
    if (is_lookup) {
        const Expr* key   = A[1];
        const Expr* deflt = (na == 3) ? A[2] : NULL;
        if (deflt && !expr_is_compile_const(c, deflt)) return false;
        if (!expr_is_compile_const(c, key)) {
            /* B2 runtime key: an int/real machine scalar; result is the value type
             * (widened by the default), regardless of the specific key. */
            CompileType kt;
            if (!infer_type(c, key, &kt) || (kt != CT_INT && kt != CT_REAL)) return false;
            CompileType rt = src_valtype;
            if (src_cst && !assoc_const_values_elem(src_cst, &rt)) return false;
            if (deflt) {
                Slot d; CompileType dt;
                if (!literal(deflt, &d, &dt)) return false;
                rt = num_common(rt, dt);
                if ((int)rt < 0 || CT_IS_ARRAY(rt) || CT_IS_ASSOC(rt)) return false;
            }
            if (rt != CT_INT && rt != CT_REAL && rt != CT_COMPLEX) return false;
            *out = rt; return true;
        }
        if (src_cst) {
            Expr* v = assoc_lookup_value(src_cst, key);
            if (!v) v = (Expr*)deflt;
            Slot imm; CompileType vt;
            if (!v || !literal(v, &imm, &vt)) return false;
            *out = vt; return true;
        }
        CompileType rt = src_valtype;
        if (deflt) {
            Slot d; CompileType dt;
            if (!literal(deflt, &d, &dt)) return false;
            rt = num_common(rt, dt);
            if ((int)rt < 0 || CT_IS_ARRAY(rt) || CT_IS_ASSOC(rt)) return false;
        }
        if (rt != CT_INT && rt != CT_REAL && rt != CT_COMPLEX) return false;
        *out = rt; return true;
    }
    /* Values */
    {
        CompileType elem = src_valtype;
        if (src_cst && !assoc_const_values_elem(src_cst, &elem)) return false;
        if (elem != CT_INT && elem != CT_REAL && elem != CT_COMPLEX) return false;
        *out = CT_ARRAY(elem, 1); return true;
    }
}

