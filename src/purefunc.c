#include "purefunc.h"
#include "modular.h"   /* scoping_avoid_values: capture avoidance */
#include "symsum.h"
#include "symtab.h"
#include "eval.h"
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#include "attr.h"
#include "sym_names.h"
#include "assoc.h"
#include "message.h"
#include "print.h"
#include <stdio.h>

static uint32_t parse_pure_attr(Expr* attr_expr) {
    if (attr_expr->type == EXPR_SYMBOL) {
        const char* name = attr_expr->data.symbol.name;
        if (name == SYM_Listable) return ATTR_LISTABLE;
        if (name == SYM_Flat) return ATTR_FLAT;
        if (name == SYM_Orderless) return ATTR_ORDERLESS;
        if (name == SYM_NumericFunction) return ATTR_NUMERICFUNCTION;
        if (name == SYM_OneIdentity) return ATTR_ONEIDENTITY;
        if (name == SYM_HoldFirst) return ATTR_HOLDFIRST;
        if (name == SYM_HoldRest) return ATTR_HOLDREST;
        if (name == SYM_HoldAll) return ATTR_HOLDALL;
        if (name == SYM_HoldAllComplete) return ATTR_HOLDALLCOMPLETE;
        if (name == SYM_SequenceHold) return ATTR_SEQUENCEHOLD;
        if (name == SYM_NHoldRest) return ATTR_NHOLDREST;
        if (name == SYM_Protected) return ATTR_PROTECTED;
    }
    return ATTR_NONE;
}

uint32_t pure_function_attributes(Expr* func) {
    if (func->type == EXPR_FUNCTION && func->data.function.head->type == EXPR_SYMBOL &&
        func->data.function.head->data.symbol.name == SYM_Function) {
        /* Mathematica default: Function[params, body] has no Hold attributes,
         * so its arguments are evaluated normally before being substituted
         * into body. The 3-arg form Function[params, body, attrs] can opt in
         * to HoldAll / HoldFirst / HoldRest / HoldAllComplete and a handful
         * of structural attributes. */
        uint32_t attrs = ATTR_NONE;
        if (func->data.function.arg_count >= 3) {
            Expr* attr_spec = func->data.function.args[2];
            if (attr_spec->type == EXPR_SYMBOL) {
                attrs |= parse_pure_attr(attr_spec);
            } else if (attr_spec->type == EXPR_FUNCTION &&
                       attr_spec->data.function.head->type == EXPR_SYMBOL &&
                       attr_spec->data.function.head->data.symbol.name == SYM_List) {
                for (size_t i = 0; i < attr_spec->data.function.arg_count; i++) {
                    attrs |= parse_pure_attr(attr_spec->data.function.args[i]);
                }
            }
        }
        return attrs;
    }
    return ATTR_NONE;
}

void purefunc_init(void) {
    symtab_add_builtin("Slot", builtin_slot);
    symtab_add_builtin("SlotSequence", builtin_slotsequence);
    /* NHoldAll (as the reference): N[(2 + #) &, 40] must not turn Slot[1]
     * into Slot[1.`40], a slot index no call can fill. */
    symtab_get_def("Slot")->attributes |= ATTR_PROTECTED | ATTR_NHOLDALL;
    symtab_get_def("SlotSequence")->attributes |= ATTR_PROTECTED | ATTR_NHOLDALL;
}

Expr* builtin_slot(Expr* res) {
    (void)res;
    // Slot stays unevaluated. 
    // It's only replaced during apply_pure_function.
    return NULL;
}

Expr* builtin_slotsequence(Expr* res) {
    (void)res;
    // SlotSequence stays unevaluated.
    return NULL;
}

/* A slot-form Function (Function[body], or Function[Null, ...]) establishes its
 * own # scope; a named Function (Function[x, ...] / Function[{...}, ...]) does
 * not rebind slots. */
static int function_is_slot_form(Expr* fn) {
    size_t fargc = fn->data.function.arg_count;
    if (fargc <= 1) return 1;                              /* Function[body] */
    Expr* p = fn->data.function.args[0];
    if (p->type == EXPR_SYMBOL && p->data.symbol.name == SYM_Null) return 1;
    return 0;                                              /* named parameter(s) */
}

/* Report a named slot (#name) that cannot be filled, Mathematica-style:
 *   Function::slota  -- the first argument is an association without the key
 *   Function::slot1  -- the first argument is not an association at all
 * The slot is then left in place (WL returns e.g. `1 + #b`). */
static void named_slot_message(Expr* fn, const char* name, Expr** args, size_t arg_count) {
    if (mth_msg_suppressed()) { mth_msg_note_named("Function", arg_count >= 1 && is_association(args[0]) ? "slota" : "slot1"); return; }
    char* fs = fn ? expr_to_string(fn) : NULL;
    if (arg_count >= 1 && is_association(args[0])) {
        char* as = expr_to_string(args[0]);
        mth_message("Function", "slota", "Named slot %s in %s cannot be filled from %s.", name, fs ? fs : "Function[...]", as ? as : "<|...|>");
        free(as);
    } else {
        Expr** cp = malloc(sizeof(Expr*) * (arg_count ? arg_count : 1));
        for (size_t i = 0; i < arg_count; i++) cp[i] = expr_copy(args[i]);
        Expr* call = expr_new_function(expr_new_symbol(SYM_List), cp, arg_count);
        free(cp);
        char* cs = expr_to_string(call);          /* "{a1, a2}" -> "[a1, a2]" */
        expr_free(call);
        size_t cl = cs ? strlen(cs) : 0;
        if (cl >= 2) { cs[0] = '['; cs[cl - 1] = ']'; }
        mth_message("Function", "slot1", "(%s)%s is expected to have an Association as the first argument.", fs ? fs : "Function[...]", cl >= 2 ? cs : "[]");
        free(cs);
    }
    free(fs);
}

/* Copy-on-write substitution over the children of e: sub(child) is called on
 * the head and every argument, and a new node is built only when one of them
 * came back as a different pointer (expr_copy is a refcount bump, so an
 * unchanged child is returned as itself).  Otherwise e itself is shared --
 * with no allocation on the way -- keeping its fixed-point stamps: a fresh
 * node would carry none, so a large slot-free subtree of the body (a
 * 10,000-term sum a pattern bound into (# === x) &) would be evaluated again
 * on every application, as in FeynCalc's Isolate. */
typedef Expr* (*SubstFn)(Expr* e, void* ctx);

static Expr* subst_children(Expr* e, SubstFn sub, void* ctx) {
    size_t n = e->data.function.arg_count;
    Expr* new_head = sub(e->data.function.head, ctx);
    Expr** new_args = NULL;
    if (new_head != e->data.function.head) {
        new_args = malloc(sizeof(Expr*) * (n ? n : 1));
        for (size_t j = 0; j < n; j++) new_args[j] = NULL;
    }
    for (size_t i = 0; i < n; i++) {
        Expr* c = sub(e->data.function.args[i], ctx);
        if (!new_args) {
            if (c == e->data.function.args[i]) { expr_free(c); continue; }
            new_args = malloc(sizeof(Expr*) * (n ? n : 1));
            for (size_t j = 0; j < n; j++) new_args[j] = NULL;
        }
        new_args[i] = c;
    }
    if (!new_args) {                     /* nothing below changed */
        expr_free(new_head);
        return expr_copy(e);
    }
    for (size_t j = 0; j < n; j++)
        if (!new_args[j]) new_args[j] = expr_copy(e->data.function.args[j]);
    Expr* res = expr_new_function(new_head, new_args, n);
    free(new_args);
    return res;
}

/* Binder-summary context of one application (symsum.h): the mask of the
 * names the walk substitutes -- the parameters, or Slot and SlotSequence for
 * a slot-form function -- and whether summaries may be consulted. */
typedef struct { uint64_t mask; bool fast; int* renames; } PfSum;

/* Capture avoidance (audit Function#2, Function#3): before substituting into
 * a nested binding scope (Function with parameters, Module, With), rename
 * every local of it that occurs in a value being substituted, so
 * Function[x, Function[y, x + y]][y] is Function[y$n, y + y$n], not
 * Function[y, y + y].  NULL when nothing clashes. */
/* Cheap pre-test for a lexical binding scope (Function / Module / With):
 * the walk asks it of every node, so it stays a few pointer compares. */
static inline bool pf_scope_head(const Expr* e) {
    const Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL) return false;
    const char* n = h->data.symbol.name;
    return n == SYM_Function || n == SYM_Module || n == SYM_With;
}

static Expr* pf_avoid_capture(Expr* e, Expr** vals, size_t count, const PfSum* s) {
    Expr* safe = scoping_avoid_values(e, vals, count, NULL, 0);
    if (safe && s->renames) (*s->renames)++;
    return safe;
}

typedef struct { Expr** args; size_t arg_count; Expr* fn; const PfSum* s; bool frozen; } SlotCtx;
static Expr* slots_rec(Expr* e, Expr** args, size_t arg_count, Expr* fn, const PfSum* s, bool frozen);
static Expr* slots_sub(Expr* e, void* ctx) {
    SlotCtx* c = ctx;
    return slots_rec(e, c->args, c->arg_count, c->fn, c->s, c->frozen);
}

/* `fn` is the Function being applied (for diagnostics only; may be NULL). */
static Expr* slots_rec(Expr* e, Expr** args, size_t arg_count, Expr* fn, const PfSum* s, bool frozen) {
    if (!e) return NULL;

    if (e->type == EXPR_FUNCTION) {
        /* No Slot / SlotSequence below: nothing to substitute. */
        uint64_t sb;
        if (s->fast && symsum_get(e, &sb) && !(sb & s->mask)) { SCOPE_STAT(pf_skip); return expr_copy(e); }
        bool cfrozen = frozen || e->refcount >= 2;
        SCOPE_STAT(pf_nodes);
        if (pf_scope_head(e) && expr_is_binding_scope(e)) {
            Expr* safe = pf_avoid_capture(e, args, arg_count, s);
            if (safe) {
                Expr* r = slots_rec(safe, args, arg_count, fn, s, false);
                expr_free(safe);
                return r;
            }
        }
        // A nested *slot-form* Function establishes its own # scope, so the
        // enclosing function's slots must not be substituted into it. A nested
        // *named* Function does not rebind slots, so an enclosing slot-form's #
        // must still flow into its body (A10) -- fall through to the standard
        // recursion below in that case.
        if (e->data.function.head->type == EXPR_SYMBOL &&
            e->data.function.head->data.symbol.name == SYM_Function &&
            function_is_slot_form(e)) {
            return expr_copy(e);
        }
        
        // Check if this function is Slot[n] or SlotSequence[n]
        if (e->data.function.head->type == EXPR_SYMBOL) {
            const char* head = e->data.function.head->data.symbol.name;
            
            if (head == SYM_Slot && e->data.function.arg_count == 1) {
                Expr* arg = e->data.function.args[0];
                if (arg->type == EXPR_INTEGER) {
                    int64_t idx = arg->data.integer;
                    if (idx >= 1 && (size_t)idx <= arg_count) {
                        return expr_copy(args[idx - 1]);
                    }
                    /* #0 is the pure function itself: the anonymous recursion
                     * If[# == 0, 1, # #0[# - 1]] & (audit Function#1, Slot#1). */
                    if (idx == 0 && fn) return expr_copy(fn);
                }
                /* Named slot #name == Slot["name"]: the value stored under the
                 * key "name" in the first argument, which must be an
                 * association (O(1) via the key index). */
                if (arg->type == EXPR_STRING) {
                    Expr* v = (arg_count >= 1 && is_association(args[0]))
                                  ? assoc_lookup_value(args[0], arg) : NULL;
                    if (v) return expr_copy(v);
                    named_slot_message(fn, arg->data.string, args, arg_count);
                    return expr_copy(e);
                }
            }
            
            if (head == SYM_SlotSequence && e->data.function.arg_count == 1) {
                Expr* arg = e->data.function.args[0];
                if (arg->type == EXPR_INTEGER) {
                    int64_t idx = arg->data.integer;
                    /* ##n is the arguments from the n-th on, which is the
                     * EMPTY sequence when exactly n - 1 were supplied:
                     * {##}&[] is {} (audit SlotSequence#1). */
                    if (idx >= 1 && (size_t)idx <= arg_count + 1) {
                        size_t seq_count = arg_count - (size_t)idx + 1;
                        Expr** seq_args = malloc(sizeof(Expr*) * (seq_count ? seq_count : 1));
                        for (size_t i = 0; i < seq_count; i++) {
                            seq_args[i] = expr_copy(args[(size_t)idx - 1 + i]);
                        }
                        Expr* seq = expr_new_function(expr_new_symbol(SYM_Sequence), seq_args, seq_count);
                        free(seq_args);
                        return seq;
                    }
                }
            }
        }
        
        // Standard recursion for other functions
        SlotCtx c = { args, arg_count, fn, s, cfrozen };
        Expr* r = subst_children(e, slots_sub, &c);
        if (r == e && s->fast && cfrozen) { symsum_fill(e); SCOPE_STAT(pf_fill); }
        return r;
    }
    
    return expr_copy(e);
}

static uint64_t g_slot_mask = 0;   /* bits of Slot | SlotSequence */

/* `frozen`: `e` is part of a shared (hence immutable) Function node. */
static Expr* substitute_slots(Expr* e, Expr** args, size_t arg_count, Expr* fn, bool frozen) {
    if (!e) return NULL;
    SCOPE_STAT(pf_top);
    int renames = 0;
    PfSum s = { 0, false, &renames };
#ifndef MTH_NO_SYMSUM
    if (g_symsum_on) {
        if (!g_slot_mask) {
            SymbolDef* a = intern_symbol_node(SYM_Slot);
            SymbolDef* b = intern_symbol_node(SYM_SlotSequence);
            symsum_register_binder(a);
            symsum_register_binder(b);
            g_slot_mask = symsum_def_bit(a) | symsum_def_bit(b);
        }
        s.mask = g_slot_mask;
        s.fast = g_symsum_on != 0;
    }
#endif
    Expr* r = slots_rec(e, args, arg_count, fn, &s, frozen || e->refcount >= 2);
#ifndef MTH_NO_SYMSUM
    /* Not re-run when a named slot could message (#name over an association). */
    if (s.fast && g_symsum_diff && !(arg_count >= 1 && is_association(args[0]))) {
        int renames2 = 0;
        PfSum slow = { 0, false, &renames2 };
        Expr* r2 = slots_rec(e, args, arg_count, fn, &slow, false);
        if (renames2 != renames || (renames == 0 && !expr_eq(r, r2))) {
            fprintf(stderr, "MATHILDA_SYMSUM_DIFF: substitute_slots differs\n");
            abort();
        }
        expr_free(r2);
    }
#endif
    return r;
}

/*
 * function_binds_name:
 * True if the nested Function expression `fn` (head == Function) rebinds
 * `name` as one of its own NAMED parameters -- i.e. `name` is shadowed inside
 * `fn`'s body. Slot-form functions (`Function[body]`, `Function[Null, body,
 * ...]`) bind only slots, never names, so they shadow nothing here.
 */
static int function_binds_name(Expr* fn, const char* name) {
    size_t fargc = fn->data.function.arg_count;
    if (fargc <= 1) return 0;                       /* Function[body] -- slot form */
    Expr* p = fn->data.function.args[0];
    if (p->type == EXPR_SYMBOL) {
        if (p->data.symbol.name == SYM_Null) return 0;   /* Function[Null, ...] -- slot form */
        return strcmp(p->data.symbol.name, name) == 0;   /* Function[x, ...] */
    }
    if (p->type == EXPR_FUNCTION &&
        p->data.function.head->type == EXPR_SYMBOL &&
        p->data.function.head->data.symbol.name == SYM_List) {
        for (size_t i = 0; i < p->data.function.arg_count; i++) {
            Expr* v = p->data.function.args[i];
            if (v->type == EXPR_SYMBOL && strcmp(v->data.symbol.name, name) == 0)
                return 1;                                /* Function[{..., x, ...}, ...] */
        }
    }
    return 0;
}

/*
 * substitute_names:
 * Walks `e` and replaces free occurrences of any parameter symbol in
 * `names` with the corresponding expression in `vals`. The substitution is
 * lexical (no evaluation).
 *
 * A nested Function IS descended into -- Mathematica's Function is a lexical
 * closure, so an outer named parameter is visible inside a nested pure
 * function unless that nested Function rebinds the SAME name. Only the
 * shadowed names are dropped from the substitution set for that subtree; a
 * slot-form inner function (`(... # ...)&`) rebinds no names, so every outer
 * name still flows into its body (this is what makes
 * `Function[x, (x + #)&][a]` yield `(a + #)&`, matching the Wolfram Language).
 * Slots themselves are handled separately by substitute_slots and are claimed
 * by the innermost enclosing Function of any form, so they are never touched
 * here.
 *
 * The returned expression is freshly allocated; the caller owns it and
 * `vals` is not consumed (copies are made as needed).
 */
typedef struct { char** names; Expr** vals; size_t count; const PfSum* s; bool frozen; } NameCtx;
static Expr* names_rec(Expr* e, char** names, Expr** vals, size_t count, const PfSum* s, bool frozen);
static Expr* names_sub(Expr* e, void* ctx) {
    NameCtx* c = ctx;
    return names_rec(e, c->names, c->vals, c->count, c->s, c->frozen);
}

/* The name a Module / With binding-list entry binds (a Symbol, or
 * Set / SetDelayed[sym, init]), or NULL. */
static const char* pf_entry_local(const Expr* b) {
    if (b->type == EXPR_SYMBOL) return b->data.symbol.name;
    if (b->type == EXPR_FUNCTION && b->data.function.arg_count == 2
        && b->data.function.head->type == EXPR_SYMBOL
        && (b->data.function.head->data.symbol.name == SYM_Set
            || b->data.function.head->data.symbol.name == SYM_SetDelayed)
        && b->data.function.args[0]->type == EXPR_SYMBOL)
        return b->data.function.args[0]->data.symbol.name;
    return NULL;
}

/* A Module / With that binds one of the parameters shadows it: the
 * parameter is substituted into the init values (evaluated outside the
 * scope) but not into the body, so Function[x, Module[{x = x + 1}, x]][5]
 * is 6.  It used to be substituted everywhere, binding list included
 * (Module[{5 = 5 + 1}, 5]).  NULL when no parameter is shadowed. */
static Expr* names_into_local_scope(Expr* e, char** names, Expr** vals, size_t count,
                                    const PfSum* s, bool frozen) {
    Expr* vars = e->data.function.args[0];
    if (!(vars->type == EXPR_FUNCTION && vars->data.function.head->type == EXPR_SYMBOL
          && vars->data.function.head->data.symbol.name == SYM_List))
        return NULL;
    size_t nb = vars->data.function.arg_count;
    size_t kept = 0;
    for (size_t i = 0; i < count; i++) {
        bool sh = false;
        for (size_t j = 0; j < nb && !sh && names[i]; j++)
            sh = pf_entry_local(vars->data.function.args[j]) == names[i];
        if (!sh) kept++;
    }
    if (kept == count) return NULL;
    char** knames = malloc(sizeof(char*) * (kept ? kept : 1));
    Expr** kvals = malloc(sizeof(Expr*) * (kept ? kept : 1));
    size_t k = 0;
    for (size_t i = 0; i < count; i++) {
        bool sh = false;
        for (size_t j = 0; j < nb && !sh && names[i]; j++)
            sh = pf_entry_local(vars->data.function.args[j]) == names[i];
        if (sh) continue;
        knames[k] = names[i]; kvals[k] = vals[i]; k++;
    }
    /* binding list: inits see every parameter */
    Expr** nbv = malloc(sizeof(Expr*) * (nb ? nb : 1));
    for (size_t j = 0; j < nb; j++) {
        Expr* b = vars->data.function.args[j];
        if (b->type == EXPR_FUNCTION && pf_entry_local(b)) {
            Expr* ba[2] = { expr_copy(b->data.function.args[0]),
                            names_rec(b->data.function.args[1], names, vals, count, s, frozen) };
            nbv[j] = expr_new_function(expr_copy(b->data.function.head), ba, 2);
        } else {
            nbv[j] = expr_copy(b);
        }
    }
    size_t n = e->data.function.arg_count;
    Expr** na = malloc(sizeof(Expr*) * n);
    na[0] = expr_new_function(expr_copy(vars->data.function.head), nbv, nb);
    free(nbv);
    for (size_t i = 1; i < n; i++)
        na[i] = kept ? names_rec(e->data.function.args[i], knames, kvals, kept, s, frozen)
                     : expr_copy(e->data.function.args[i]);
    Expr* r = expr_new_function(expr_copy(e->data.function.head), na, n);
    free(na);
    free(knames);
    free(kvals);
    return r;
}

static Expr* names_rec(Expr* e, char** names, Expr** vals, size_t count, const PfSum* s, bool frozen) {
    if (!e) return NULL;

    if (e->type == EXPR_SYMBOL) {
        if (s->fast && !(symsum_sym_bit(e) & s->mask)) return expr_copy(e);
        for (size_t i = 0; i < count; i++) {
            const char* nm = e->data.symbol.name;
            /* names are interned, so the pointer test decides almost always;
             * the first-character guard keeps the fallback strcmp rare */
            if (names[i] && (names[i] == nm || (names[i][0] == nm[0] && strcmp(nm, names[i]) == 0))) {
                return expr_copy(vals[i]);
            }
        }
        return expr_copy(e);
    }

    if (e->type == EXPR_FUNCTION) {
        /* No parameter name below: nothing to substitute.  (The mask is the
         * outermost set's; a nested Function's shadowing only shrinks it.) */
        uint64_t sb;
        if (s->fast && symsum_get(e, &sb) && !(sb & s->mask)) { SCOPE_STAT(pf_skip); return expr_copy(e); }
        bool cfrozen = frozen || e->refcount >= 2;
        SCOPE_STAT(pf_nodes);
        if (pf_scope_head(e) && expr_is_binding_scope(e)) {
            Expr* safe = pf_avoid_capture(e, vals, count, s);
            if (safe) {
                Expr* r = names_rec(safe, names, vals, count, s, false);
                expr_free(safe);
                return r;
            }
            if (e->data.function.head->data.symbol.name != SYM_Function) {
                Expr* r = names_into_local_scope(e, names, vals, count, s, cfrozen);
                if (r) return r;
            }
        }
        /* For a nested Function, drop from the substitution set any name it
         * rebinds as a named parameter (a slot-form inner function rebinds
         * none). Non-Function heads keep the full set. */
        char** use_names = names;
        Expr** use_vals = vals;
        size_t use_count = count;
        char** alloc_names = NULL;
        Expr** alloc_vals = NULL;
        if (e->data.function.head->type == EXPR_SYMBOL &&
            e->data.function.head->data.symbol.name == SYM_Function) {
            size_t kept = 0;
            for (size_t i = 0; i < count; i++)
                if (!(names[i] && function_binds_name(e, names[i]))) kept++;
            if (kept == 0) return expr_copy(e);          /* every name shadowed */
            if (kept != count) {
                alloc_names = malloc(sizeof(char*) * kept);
                alloc_vals = malloc(sizeof(Expr*) * kept);
                size_t j = 0;
                for (size_t i = 0; i < count; i++) {
                    if (names[i] && function_binds_name(e, names[i])) continue;
                    alloc_names[j] = names[i];
                    alloc_vals[j] = vals[i];
                    j++;
                }
                use_names = alloc_names;
                use_vals = alloc_vals;
                use_count = kept;
            }
        }

        NameCtx c = { use_names, use_vals, use_count, s, cfrozen };
        Expr* res = subst_children(e, names_sub, &c);
        free(alloc_names);
        free(alloc_vals);
        if (res == e && s->fast && cfrozen) { symsum_fill(e); SCOPE_STAT(pf_fill); }
        return res;
    }

    return expr_copy(e);
}

/* Replace the free occurrences of the parameters `names` by `vals` in a
 * Function body (lexical; nested Functions shadow).  A subtree whose binder
 * summary holds none of the parameters' bits is returned shared unvisited. */
static Expr* substitute_names(Expr* e, char** names, Expr** vals, size_t count, bool frozen) {
    if (!e) return NULL;
    SCOPE_STAT(pf_top);
    int renames = 0;
    PfSum s = { 0, false, &renames };
#ifndef MTH_NO_SYMSUM
    if (g_symsum_on) {
        for (size_t i = 0; i < count; i++) {
            if (!names[i]) continue;
            SymbolDef* d = intern_symbol_node(names[i]);
            symsum_register_binder(d);
            s.mask |= symsum_def_bit(d);
        }
        s.fast = g_symsum_on != 0;
    }
#endif
    Expr* r = names_rec(e, names, vals, count, &s, frozen || e->refcount >= 2);
#ifndef MTH_NO_SYMSUM
    if (s.fast && g_symsum_diff) {
        int renames2 = 0;
        PfSum slow = { 0, false, &renames2 };
        Expr* r2 = names_rec(e, names, vals, count, &slow, false);
        if (renames2 != renames || (renames == 0 && !expr_eq(r, r2))) {
            fprintf(stderr, "MATHILDA_SYMSUM_DIFF: substitute_names differs\n");
            abort();
        }
        expr_free(r2);
    }
#endif
    return r;
}

/*
 * apply_pure_function:
 * Applies a Function[...] expression to the given arguments.
 *
 * Four input shapes are supported:
 *   Function[body]                    -- body uses # / ## slots
 *   Function[Null, body, attrs]       -- slot form with attributes
 *   Function[x, body]                 -- single named parameter
 *   Function[{x1, ...}, body]         -- list of named parameters
 *   Function[{x1, ...}, body, attrs]  -- named parameters with attributes
 *
 * Parameter binding is performed by lexical substitution (not by mutating
 * the global symbol table) so that references wrapped in Unevaluated still
 * see the substituted expression. Substitution descends into nested Function
 * expressions as a lexical closure would, stopping at a name only where a
 * nested Function rebinds it (see substitute_names).
 */
/*
 * trap_return:
 * Strip a Return[v, Function] marker. A plain Return[v] / Return[] is NOT
 * consumed: as in the reference it passes out of the pure function to the
 * enclosing definition (or Do / Scan), so Map[(Return[#]; 0) &, {1, 2}] is
 * {Return[1], Return[2]} and (Return[1]; 2) &[] inside a definition leaves
 * the definition. Return[v, h] for another h is handed upward unchanged.
 *
 * Takes ownership of `result` and returns either the original pointer
 * (when the marker is absent or propagating) or a freshly-owned
 * unwrapped value (when the marker is consumed). Frees `result` on the
 * consume path.
 */
static Expr* apply_pure_function_sorted(Expr* head, Expr** args, size_t arg_count);

static Expr* trap_return(Expr* result) {
    Expr* rv = NULL;
    EvalReturnAction ra = eval_classify_return_targeted(result, SYM_Function, &rv);
    if (ra == EVAL_RETURN_CONSUME) {
        expr_free(result);
        return rv;
    }
    return result;
}

/* Function::fpct -- a named-parameter Function called with fewer arguments
 * than it has parameters; the call stays unevaluated, as in the reference
 * (audit Function#5: Function[{x, y}, x + y][1] gave 1 + y). */
static Expr* pf_fpct(Expr* head, Expr* params, Expr** args, size_t arg_count) {
    char* ps = expr_to_string(params);
    char* hs = expr_to_string(head);
    size_t cap = 64 + (hs ? strlen(hs) : 0);
    char** as = malloc(sizeof(char*) * (arg_count ? arg_count : 1));
    for (size_t i = 0; i < arg_count; i++) { as[i] = expr_to_string(args[i]); cap += strlen(as[i]) + 2; }
    char* cs = malloc(cap);
    size_t o = (size_t)snprintf(cs, cap, "%s[", hs ? hs : "");
    for (size_t i = 0; i < arg_count; i++) {
        o += (size_t)snprintf(cs + o, cap - o, "%s%s", i ? ", " : "", as[i]);
        free(as[i]);
    }
    snprintf(cs + o, cap - o, "]");
    mth_message("Function", "fpct", "Too many parameters in %s to be filled from %s.",
                ps ? ps : "", cs);
    free(as); free(cs); free(ps); free(hs);
    return NULL;
}

Expr* apply_pure_function(Expr* head, Expr** args, size_t arg_count) {
    if (head->type != EXPR_FUNCTION) return NULL;

    size_t head_argc = head->data.function.arg_count;

    /* Function[params, body, Orderless]: the arguments are taken in canonical
     * order, as for an Orderless symbol head (audit Function#4).  Sorted in a
     * private copy of the pointer array; the call itself is not mutated. */
    if (head_argc >= 3 && arg_count > 1 &&
        (pure_function_attributes(head) & ATTR_ORDERLESS)) {
        Expr** sorted = malloc(sizeof(Expr*) * arg_count);
        memcpy(sorted, args, sizeof(Expr*) * arg_count);
        expr_orderless_sort(sorted, arg_count);
        Expr* r = apply_pure_function_sorted(head, sorted, arg_count);
        free(sorted);
        return r;
    }
    return apply_pure_function_sorted(head, args, arg_count);
}

static Expr* apply_pure_function_sorted(Expr* head, Expr** args, size_t arg_count) {
    size_t head_argc = head->data.function.arg_count;

    /* Function[body] -- slot form with 1 argument */
    if (head_argc == 1) {
        Expr* body = head->data.function.args[0];
        Expr* substituted = substitute_slots(body, args, arg_count, head, head->refcount >= 2);
        Expr* result = evaluate(substituted);
        expr_free(substituted);
        return trap_return(result);
    }

    if (head_argc < 2) return NULL;

    Expr* params = head->data.function.args[0];
    Expr* body = head->data.function.args[1];

    /* Function[Null, body, ...] -- slot form with attributes.
     * (Null indicates "use slots" rather than named parameters.) */
    if (params->type == EXPR_SYMBOL && params->data.symbol.name == SYM_Null) {
        Expr* substituted = substitute_slots(body, args, arg_count, head, head->refcount >= 2);
        Expr* result = evaluate(substituted);
        expr_free(substituted);
        return trap_return(result);
    }

    /* Function[x, body] -- single named parameter */
    if (params->type == EXPR_SYMBOL) {
        if (arg_count < 1) return pf_fpct(head, params, args, arg_count);
        char* name = params->data.symbol.name;
        Expr* val = args[0];
        Expr* substituted = substitute_names(body, &name, &val, 1, head->refcount >= 2);
        Expr* result = evaluate(substituted);
        expr_free(substituted);
        return trap_return(result);
    }

    /* Function[{x1, x2, ...}, body] -- list of named parameters */
    if (params->type == EXPR_FUNCTION &&
        params->data.function.head->type == EXPR_SYMBOL &&
        params->data.function.head->data.symbol.name == SYM_List) {
        size_t var_count = params->data.function.arg_count;
        if (arg_count < var_count) return pf_fpct(head, params, args, arg_count);
        if (var_count == 0) {
            Expr* result = evaluate(body);
            return trap_return(result);
        }

        char** names = malloc(sizeof(char*) * var_count);
        Expr** vals = malloc(sizeof(Expr*) * var_count);
        for (size_t i = 0; i < var_count; i++) {
            Expr* var = params->data.function.args[i];
            if (var->type == EXPR_SYMBOL) {
                names[i] = var->data.symbol.name;
                /* Missing arguments leave the parameter symbolic, which is
                 * the same behaviour as unbound slots. */
                vals[i] = (i < arg_count) ? args[i] : var;
            } else {
                names[i] = NULL;
                vals[i] = var;
            }
        }

        Expr* substituted = substitute_names(body, names, vals, var_count, head->refcount >= 2);
        Expr* result = evaluate(substituted);
        expr_free(substituted);
        free(names);
        free(vals);
        return trap_return(result);
    }

    return NULL;
}
