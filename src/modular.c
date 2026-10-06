#include "modular.h"
#include "symtab.h"
#include "eval.h"
#include "attr.h"
#include "sym_names.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

static int64_t module_number = 1;

void modular_init(void) {
    symtab_add_builtin("Module", builtin_module);
    symtab_add_builtin("Block", builtin_block);
    symtab_add_builtin("With", builtin_with);
    symtab_add_builtin("Unique", builtin_unique);

    // Initial value for $ModuleNumber
    Expr* mn = expr_new_integer(module_number);
    Expr* sym_mn = expr_new_symbol(SYM_DollarModuleNumber);
    symtab_add_own_value("$ModuleNumber", sym_mn, mn);
    expr_free(mn);
    expr_free(sym_mn);
}

typedef struct ScopingEnv {
    const char* old_name;
    Expr* replacement;
    struct ScopingEnv* next;
} ScopingEnv;

/* The ITERATOR family: head[body, {i, ...}, ...] -- body in arg 0, one or more
 * iterator specs after it, each binding its first element. Table, Do, Sum and
 * Product are structurally identical here and bind their iterator the same way,
 * so capture-avoidance must treat them alike. Only Table used to be listed,
 * which is exactly MATHILDA_DIVERGENCES.md A18: `g[v_] := Module[{s = 0},
 * Do[s += v, {k, 2}]; s]; g[k]` gave 3 where Mathematica gives 2 k, because the
 * caller's `k` was injected into a body whose `Do` then bound it. Sum and
 * Product had it too (`Sum[v, {k, 2}]` → 3, `Product` → 2). */
static bool is_iterator_scope_head(const char* h) {
    return h == SYM_Table || h == SYM_Do || h == SYM_Sum || h == SYM_Product;
}

static bool is_scoping_construct(Expr* e) {
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    return h == SYM_Module || h == SYM_Block || h == SYM_With ||
           h == SYM_Function || is_iterator_scope_head(h);
}

// Does this construct carry its bound names in argument 0?
//
//   Module/Block/With[{vars}, body]  yes -- arg 0 is the binding list.
//   Function[params, body, attrs]    yes -- arg 0 is the parameter spec,
//                                    either a bare symbol or a List of them.
//   Function[body]                   NO -- the single argument is the *body*
//                                    of a Slot form (`(# + 1) &`). Nothing is
//                                    bound: `#` is not a symbol, so there is
//                                    no name to shadow. Treating arg 0 as a
//                                    binding list here copies the body
//                                    verbatim and leaves every enclosing local
//                                    free inside it -- so `Module[{c = 0},
//                                    Scan[(c = c + 1) &, ...]]` incremented the
//                                    *global* c (self-referentially, hence
//                                    $RecursionLimit) and left the local at 0.
//   Table[body, {i, ...}]            NO -- it binds in the iterators, args 1..
static bool scoping_binds_in_arg0(Expr* e) {
    if (!is_scoping_construct(e)) return false;
    const char* h = e->data.function.head->data.symbol.name;
    if (is_iterator_scope_head(h)) return false;   /* arg 0 is the body */
    if (h == SYM_Function) return e->data.function.arg_count >= 2;
    return e->data.function.arg_count >= 1;
}

static Expr* substitute_scoping(Expr* e, ScopingEnv* env) {
    if (!e) return NULL;
    if (e->type == EXPR_SYMBOL) {
        ScopingEnv* curr = env;
        while (curr) {
            if (strcmp(e->data.symbol.name, curr->old_name) == 0) {
                return expr_copy(curr->replacement);
            }
            curr = curr->next;
        }
        return expr_copy(e);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy(e);

    // The iterator family (Table/Do/Sum/Product) binds its iterator variables in
    // the *iterator* specs (args 1..), each of the form {var, ...}; arg 0 is the
    // body. Every other scoping construct (Module/Block/With/Function) binds in
    // arg 0. They therefore need completely different substitution handling, so
    // detect the iterator shape here.
    bool binds_in_iterators = is_scoping_construct(e)
        && is_iterator_scope_head(e->data.function.head->data.symbol.name);

    // Handle shadowing in scoping constructs
    ScopingEnv* filtered_env = env;

    if (binds_in_iterators || scoping_binds_in_arg0(e)) {
        // Collect the names this construct binds, so they are removed from the
        // env we push into the body (lexical shadowing).
        const char* shadow_buf[64];
        size_t nshadow = 0;
        if (binds_in_iterators) {
            for (size_t k = 1; k < e->data.function.arg_count && nshadow < 64; k++) {
                Expr* it = e->data.function.args[k];
                if (it->type == EXPR_FUNCTION
                    && it->data.function.head->type == EXPR_SYMBOL
                    && it->data.function.head->data.symbol.name == SYM_List
                    && it->data.function.arg_count >= 2   /* {n} is a count, binds nothing */
                    && it->data.function.args[0]->type == EXPR_SYMBOL) {
                    shadow_buf[nshadow++] = it->data.function.args[0]->data.symbol.name;
                }
            }
        } else {
            Expr* vars = e->data.function.args[0];
            if (vars->type == EXPR_SYMBOL) {
                // Function[x, body] -- a single bare parameter still binds `x`,
                // so it must shadow an enclosing local of the same name.
                shadow_buf[nshadow++] = vars->data.symbol.name;
            } else if (vars->type == EXPR_FUNCTION
                       && vars->data.function.head->type == EXPR_SYMBOL
                       && vars->data.function.head->data.symbol.name == SYM_List) {
                for (size_t i = 0; i < vars->data.function.arg_count && nshadow < 64; i++) {
                    Expr* v = vars->data.function.args[i];
                    const char* nm = NULL;
                    if (v->type == EXPR_SYMBOL) nm = v->data.symbol.name;
                    else if (v->type == EXPR_FUNCTION && v->data.function.head->data.symbol.name == SYM_Set && v->data.function.arg_count == 2) {
                        if (v->data.function.args[0]->type == EXPR_SYMBOL) nm = v->data.function.args[0]->data.symbol.name;
                    }
                    if (nm) shadow_buf[nshadow++] = nm;
                }
            }
        }

        if (nshadow > 0) {
            // Rebuild the env list skipping every shadowed name.
            ScopingEnv* new_env = NULL;
            for (ScopingEnv* curr = env; curr; curr = curr->next) {
                bool shadowed = false;
                for (size_t s = 0; s < nshadow; s++)
                    if (strcmp(curr->old_name, shadow_buf[s]) == 0) { shadowed = true; break; }
                if (!shadowed) {
                    ScopingEnv* node = malloc(sizeof(ScopingEnv));
                    node->old_name = curr->old_name;
                    node->replacement = curr->replacement;
                    node->next = new_env;
                    new_env = node;
                }
            }
            filtered_env = new_env;
        }
    }

    Expr** new_args = malloc(sizeof(Expr*) * e->data.function.arg_count);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        // Table: arg 0 is the body (substitute normally with the shadowed
        // env); args 1.. are iterator specs {var, lim...} where `var` is a
        // binding occurrence (copied) and the limits are substituted.
        if (binds_in_iterators) {
            if (i == 0) {
                new_args[i] = substitute_scoping(e->data.function.args[i], filtered_env);
            } else {
                Expr* it = e->data.function.args[i];
                // A length-1 spec {n} is an iteration COUNT, not a binding: `n`
                // must be substituted (With[{n=3}, Table[x, {n}]] -> {x,x,x}).
                // Only {i, lim, ...} (length >= 2) binds `i` as an iterator var.
                if (it->type == EXPR_FUNCTION
                    && it->data.function.head->type == EXPR_SYMBOL
                    && it->data.function.head->data.symbol.name == SYM_List
                    && it->data.function.arg_count >= 2
                    && it->data.function.args[0]->type == EXPR_SYMBOL) {
                    size_t na = it->data.function.arg_count;
                    Expr** nb = malloc(sizeof(Expr*) * na);
                    nb[0] = expr_copy(it->data.function.args[0]);
                    for (size_t j = 1; j < na; j++)
                        nb[j] = substitute_scoping(it->data.function.args[j], filtered_env);
                    Expr* lhead = expr_copy(it->data.function.head);
                    new_args[i] = expr_new_function(lhead, nb, na);
                    free(nb);
                } else {
                    new_args[i] = substitute_scoping(it, filtered_env);
                }
            }
            continue;
        }
        // First argument of scoping constructs is the variable list:
        // substitute into each binding's RHS (which sees the outer
        // scope, i.e. the original env -- the rebound names propagate
        // back into bindings under simultaneous-binding semantics) but
        // NOT into the LHS name (which is a binding occurrence).
        // Without this, `With[{q = 12}, With[{k = q}, k]]` would leave
        // the inner `q` in `k = q` as a free symbol, and the CRC table
        // pattern `With[{q = 4 a c - b^2, k = (4 c)/q}, ...]` (which
        // relies on the outer-substituted k getting the value of q)
        // would fall apart on every q-dependent recursion.
        if (i == 0 && scoping_binds_in_arg0(e)) {
            Expr* vars_list = e->data.function.args[i];
            if (vars_list->type == EXPR_FUNCTION
                && vars_list->data.function.head
                && vars_list->data.function.head->type == EXPR_SYMBOL
                && vars_list->data.function.head->data.symbol.name == SYM_List) {
                size_t nb = vars_list->data.function.arg_count;
                Expr** new_b = malloc(sizeof(Expr*) * (nb > 0 ? nb : 1));
                for (size_t bi = 0; bi < nb; bi++) {
                    Expr* b = vars_list->data.function.args[bi];
                    if (b->type == EXPR_FUNCTION && b->data.function.head
                        && b->data.function.head->type == EXPR_SYMBOL
                        && (b->data.function.head->data.symbol.name == SYM_Set
                            || b->data.function.head->data.symbol.name == SYM_SetDelayed)
                        && b->data.function.arg_count == 2) {
                        Expr* lhs = expr_copy(b->data.function.args[0]);
                        Expr* rhs = substitute_scoping(b->data.function.args[1], env);
                        Expr* head_copy = expr_copy(b->data.function.head);
                        Expr* bargs[2] = { lhs, rhs };
                        new_b[bi] = expr_new_function(head_copy, bargs, 2);
                    } else {
                        new_b[bi] = expr_copy(b);
                    }
                }
                Expr* lhead_copy = expr_copy(vars_list->data.function.head);
                new_args[i] = expr_new_function(lhead_copy, new_b, nb);
                free(new_b);
            } else {
                new_args[i] = expr_copy(vars_list);
            }
        } else {
            new_args[i] = substitute_scoping(e->data.function.args[i], filtered_env);
        }
    }
    Expr* new_head = substitute_scoping(e->data.function.head, filtered_env);
    Expr* res = expr_new_function(new_head, new_args, e->data.function.arg_count);
    free(new_args);

    if (filtered_env != env) {
        ScopingEnv* tmp = filtered_env;
        while (tmp) {
            ScopingEnv* next = tmp->next;
            free(tmp);
            tmp = next;
        }
    }

    return res;
}

/* Re-register the user-visible $ModuleNumber OwnValue from the static counter,
 * mirroring builtin_module so Unique and Module share one monotone source. */
static void unique_sync_module_number(void) {
    Expr* mn_sym = expr_new_symbol(SYM_DollarModuleNumber);
    Expr* val = expr_new_integer(module_number);
    symtab_add_own_value("$ModuleNumber", mn_sym, val);
    expr_free(mn_sym);
    expr_free(val);
}

/* A Unique name-prefix source: a Symbol contributes its name, a String its
 * characters; Unique[] (no source) uses "$". Anything else -> NULL. */
static const char* unique_prefix_of(Expr* e) {
    if (!e) return "$";
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name;
    if (e->type == EXPR_STRING) return e->data.string;
    return NULL;
}

/* Build a batch of fresh symbols, one per prefix, all sharing a single counter
 * value chosen so that none of prefix<n> already names a symbol -- guaranteeing
 * the whole batch is fresh and distinct (distinct prefixes, or distinct via the
 * uniqueness scan). Each is registered Temporary. Returns count on success, or
 * -1 if any prefix is unusable (caller then returns NULL to leave unevaluated).
 * out[] must have room for `count` Expr*. */
static int unique_make_batch(const char** prefixes, size_t count, Expr** out) {
    char buf[512];
    for (size_t i = 0; i < count; i++) {
        if (!prefixes[i]) return -1;
    }
    /* Advance the shared counter until this value is clear for every prefix. */
    for (;;) {
        int64_t n = module_number;
        bool clash = false;
        for (size_t i = 0; i < count && !clash; i++) {
            snprintf(buf, sizeof(buf), "%s%lld", prefixes[i], (long long)n);
            if (symtab_lookup(buf) != NULL) clash = true;
        }
        if (!clash) break;
        module_number++;
    }
    int64_t n = module_number++;
    for (size_t i = 0; i < count; i++) {
        snprintf(buf, sizeof(buf), "%s%lld", prefixes[i], (long long)n);
        out[i] = expr_new_symbol(buf);
        symtab_get_def(buf)->attributes |= ATTR_TEMPORARY;
    }
    unique_sync_module_number();
    return (int)count;
}

/* Unique[] / Unique["x"] / Unique[x] / Unique[{a, b, ...}] -- generate one or
 * more fresh, never-before-used symbols with the Temporary attribute, drawing
 * their numeric suffix from the shared $ModuleNumber counter. A list argument
 * yields a list of symbols sharing one suffix (Wolfram behaviour). */
Expr* builtin_unique(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    /* Unique[] -> a single "$"<n> symbol. */
    if (argc == 0) {
        const char* pfx = "$";
        Expr* out[1];
        if (unique_make_batch(&pfx, 1, out) < 0) return NULL;
        return out[0];
    }
    if (argc != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    /* Unique[{a, b, ...}] -> list of fresh symbols, one per element. */
    if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL
        && arg->data.function.head->data.symbol.name == SYM_List) {
        size_t n = arg->data.function.arg_count;
        const char** prefixes = malloc(n * sizeof(char*));
        for (size_t i = 0; i < n; i++) {
            prefixes[i] = unique_prefix_of(arg->data.function.args[i]);
        }
        Expr** syms = malloc(n * sizeof(Expr*));
        int made = unique_make_batch(prefixes, n, syms);
        free(prefixes);
        if (made < 0) { free(syms); return NULL; }
        Expr* list = expr_new_function(expr_new_symbol(SYM_List), NULL, n);
        for (size_t i = 0; i < n; i++) list->data.function.args[i] = syms[i];
        free(syms);
        return list;
    }

    /* Unique["x"] / Unique[x] -> one fresh symbol. */
    const char* pfx = unique_prefix_of(arg);
    if (!pfx) return NULL;
    Expr* out[1];
    if (unique_make_batch(&pfx, 1, out) < 0) return NULL;
    return out[0];
}

/* ------------------------------------------------------------------------
 * Capture-avoiding pattern substitution support (used by replace_bindings).
 *
 * A rule's RHS may hold a scoping construct (Module/Block/With/Function/Table).
 * When a pattern binding carries a free symbol whose name equals one of that
 * construct's bound locals, the naive substitution injects the symbol into the
 * body, where the construct's later localization captures it:
 *   g[v_] := Module[{e = 1}, v + e];  g[e + 1]   gave 3, must give 2 + e.
 * The fix renames every colliding local to a fresh symbol *before* the binding
 * is substituted, so the injected symbol and the (renamed) local stay distinct.
 * See MATHILDA_DIVERGENCES.md A11.
 * ------------------------------------------------------------------------ */

/* True when `e` is a scoping construct that actually binds names.
 *
 * Block is deliberately EXCLUDED: its scope is DYNAMIC, not lexical, so a
 * symbol arriving from a caller's value is precisely what it means to rebind,
 * and renaming the local to avoid "capturing" it destroys the construct's whole
 * purpose. Mathematica agrees -- with the A11 example read over Block,
 *
 *     g[v_] := Block[{e = 1}, v + e];   g[e + 1]
 *
 * Mathematica gives 3 (the injected e sees the Block's value), where the same
 * shape over Module gives 2 + e. Renaming here broke the standard idiom of
 * Block-ing a function symbol to install a temporary rewrite hook and then
 * evaluating a held body that calls it: the hook was installed on a renamed
 * symbol and the body kept reaching the original. */
bool expr_is_binding_scope(Expr* e) {
    if (!is_scoping_construct(e)) return false;
    const char* h = e->data.function.head->data.symbol.name;
    if (h == SYM_Block) return false;                          /* dynamic scope */
    if (is_iterator_scope_head(h)) return e->data.function.arg_count >= 2; /* body + iterator(s) */
    return scoping_binds_in_arg0(e);
}

/* Collect the names `e` binds into out[0..*n) (bounded by cap). Mirrors the
 * shadow-name collection in substitute_scoping. */
static void scoping_collect_locals(Expr* e, const char** out, size_t* n, size_t cap) {
    const char* h = e->data.function.head->data.symbol.name;
    if (is_iterator_scope_head(h)) {
        for (size_t k = 1; k < e->data.function.arg_count && *n < cap; k++) {
            Expr* it = e->data.function.args[k];
            if (it->type == EXPR_FUNCTION
                && it->data.function.head->type == EXPR_SYMBOL
                && it->data.function.head->data.symbol.name == SYM_List
                && it->data.function.arg_count >= 2
                && it->data.function.args[0]->type == EXPR_SYMBOL) {
                out[(*n)++] = it->data.function.args[0]->data.symbol.name;
            }
        }
        return;
    }
    Expr* vars = e->data.function.args[0];
    if (vars->type == EXPR_SYMBOL) {                    /* Function[x, body] */
        if (*n < cap) out[(*n)++] = vars->data.symbol.name;
    } else if (vars->type == EXPR_FUNCTION
               && vars->data.function.head->type == EXPR_SYMBOL
               && vars->data.function.head->data.symbol.name == SYM_List) {
        for (size_t i = 0; i < vars->data.function.arg_count && *n < cap; i++) {
            Expr* v = vars->data.function.args[i];
            if (v->type == EXPR_SYMBOL) out[(*n)++] = v->data.symbol.name;
            else if (v->type == EXPR_FUNCTION
                     && v->data.function.head->data.symbol.name == SYM_Set
                     && v->data.function.arg_count == 2
                     && v->data.function.args[0]->type == EXPR_SYMBOL)
                out[(*n)++] = v->data.function.args[0]->data.symbol.name;
        }
    }
}

/* Rename a bare binding-occurrence symbol per the rename env, else copy it. */
static Expr* rename_bound_symbol(Expr* sym, ScopingEnv* ren) {
    if (sym->type == EXPR_SYMBOL)
        for (ScopingEnv* c = ren; c; c = c->next)
            if (strcmp(sym->data.symbol.name, c->old_name) == 0)
                return expr_copy(c->replacement);
    return expr_copy(sym);
}

/* One binding-list entry: bare `L`, or `Set[L, init]`/`SetDelayed[L, init]`.
 * Rename L per `ren`; substitute `ren` into the init value (which sees the
 * renamed locals). */
static Expr* rename_binding_entry(Expr* b, ScopingEnv* ren) {
    if (b->type == EXPR_FUNCTION && b->data.function.head->type == EXPR_SYMBOL
        && (b->data.function.head->data.symbol.name == SYM_Set
            || b->data.function.head->data.symbol.name == SYM_SetDelayed)
        && b->data.function.arg_count == 2
        && b->data.function.args[0]->type == EXPR_SYMBOL) {
        Expr* lhs = rename_bound_symbol(b->data.function.args[0], ren);
        Expr* rhs = substitute_scoping(b->data.function.args[1], ren);
        Expr* bargs[2] = { lhs, rhs };
        return expr_new_function(expr_copy(b->data.function.head), bargs, 2);
    }
    return rename_bound_symbol(b, ren);
}

/* Rebuild `e` with its own bound locals alpha-renamed per `ren`. Bodies, init
 * values and iterator limits are rewritten with substitute_scoping (which
 * respects any inner scope that re-binds the same name); binding occurrences
 * are renamed directly. */
static Expr* scoping_apply_rename(Expr* e, ScopingEnv* ren) {
    const char* h = e->data.function.head->data.symbol.name;
    size_t argc = e->data.function.arg_count;
    Expr** na = malloc(sizeof(Expr*) * (argc > 0 ? argc : 1));

    if (is_iterator_scope_head(h)) {
        na[0] = substitute_scoping(e->data.function.args[0], ren);   /* body */
        for (size_t i = 1; i < argc; i++) {
            Expr* it = e->data.function.args[i];
            if (it->type == EXPR_FUNCTION
                && it->data.function.head->type == EXPR_SYMBOL
                && it->data.function.head->data.symbol.name == SYM_List
                && it->data.function.arg_count >= 2
                && it->data.function.args[0]->type == EXPR_SYMBOL) {
                size_t nn = it->data.function.arg_count;
                Expr** nb = malloc(sizeof(Expr*) * nn);
                nb[0] = rename_bound_symbol(it->data.function.args[0], ren);
                for (size_t j = 1; j < nn; j++)
                    nb[j] = substitute_scoping(it->data.function.args[j], ren);
                na[i] = expr_new_function(expr_copy(it->data.function.head), nb, nn);
                free(nb);
            } else {
                na[i] = substitute_scoping(it, ren);
            }
        }
    } else if (h == SYM_Function) {
        Expr* params = e->data.function.args[0];
        if (params->type == EXPR_FUNCTION
            && params->data.function.head->type == EXPR_SYMBOL
            && params->data.function.head->data.symbol.name == SYM_List) {
            size_t nn = params->data.function.arg_count;
            Expr** nb = malloc(sizeof(Expr*) * (nn > 0 ? nn : 1));
            for (size_t j = 0; j < nn; j++)
                nb[j] = rename_bound_symbol(params->data.function.args[j], ren);
            na[0] = expr_new_function(expr_copy(params->data.function.head), nb, nn);
            free(nb);
        } else {
            na[0] = rename_bound_symbol(params, ren);
        }
        for (size_t i = 1; i < argc; i++)
            na[i] = substitute_scoping(e->data.function.args[i], ren);
    } else {   /* Module / Block / With : args[0] = List[bindings], args[1] = body */
        Expr* vars = e->data.function.args[0];
        if (vars->type == EXPR_FUNCTION
            && vars->data.function.head->type == EXPR_SYMBOL
            && vars->data.function.head->data.symbol.name == SYM_List) {
            size_t nn = vars->data.function.arg_count;
            Expr** nb = malloc(sizeof(Expr*) * (nn > 0 ? nn : 1));
            for (size_t j = 0; j < nn; j++)
                nb[j] = rename_binding_entry(vars->data.function.args[j], ren);
            na[0] = expr_new_function(expr_copy(vars->data.function.head), nb, nn);
            free(nb);
        } else {
            na[0] = substitute_scoping(vars, ren);
        }
        for (size_t i = 1; i < argc; i++)
            na[i] = substitute_scoping(e->data.function.args[i], ren);
    }

    Expr* r = expr_new_function(expr_copy(e->data.function.head), na, argc);
    free(na);
    return r;
}

/* If substituting values carrying any name in danger[0..ndanger) into `e`'s body
 * could capture one of e's bound locals, return a fresh copy of `e` with the
 * colliding locals alpha-renamed; otherwise return NULL (caller keeps `e`). */
Expr* scoping_capture_avoid(Expr* e, const char** danger, size_t ndanger) {
    if (ndanger == 0 || !expr_is_binding_scope(e)) return NULL;

    const char* locals[64]; size_t nloc = 0;
    scoping_collect_locals(e, locals, &nloc, 64);

    const char* clash[64]; size_t nclash = 0;
    for (size_t i = 0; i < nloc && nclash < 64; i++)
        for (size_t d = 0; d < ndanger; d++)
            if (strcmp(locals[i], danger[d]) == 0) { clash[nclash++] = locals[i]; break; }
    if (nclash == 0) return NULL;

    /* A fresh symbol per colliding local: prefix "L$" so unique_make_batch
     * mints "L$<n>", matching Module's own naming. */
    char* store[64]; const char* prefixes[64];
    for (size_t i = 0; i < nclash; i++) {
        size_t len = strlen(clash[i]);
        char* p = malloc(len + 2);
        memcpy(p, clash[i], len); p[len] = '$'; p[len + 1] = '\0';
        store[i] = p; prefixes[i] = p;
    }
    Expr* fresh[64];
    int made = unique_make_batch(prefixes, nclash, fresh);
    for (size_t i = 0; i < nclash; i++) free(store[i]);
    if (made < 0) return NULL;

    ScopingEnv* ren = NULL;
    for (size_t i = 0; i < nclash; i++) {
        ScopingEnv* node = malloc(sizeof(ScopingEnv));
        node->old_name = clash[i];
        node->replacement = fresh[i];   /* freed with the env below */
        node->next = ren;
        ren = node;
    }

    Expr* result = scoping_apply_rename(e, ren);

    while (ren) {
        ScopingEnv* nx = ren->next;
        expr_free(ren->replacement);
        free(ren);
        ren = nx;
    }
    return result;
}

Expr* builtin_module(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* vars = res->data.function.args[0];
    Expr* body = res->data.function.args[1];

    if (vars->type != EXPR_FUNCTION || vars->data.function.head->data.symbol.name != SYM_List) return NULL;

    size_t var_count = vars->data.function.arg_count;
    ScopingEnv* env = NULL;

    // Increment $ModuleNumber
    Expr* mn_sym = expr_new_symbol(SYM_DollarModuleNumber);
    Expr* mn_val_expr = evaluate(mn_sym);
    if (mn_val_expr->type == EXPR_INTEGER) {
        module_number = mn_val_expr->data.integer;
    }
    expr_free(mn_val_expr);
    int64_t current_mn = module_number++;
    Expr* next_mn = expr_new_integer(module_number);
    symtab_add_own_value("$ModuleNumber", mn_sym, next_mn);
    expr_free(mn_sym);
    expr_free(next_mn);

    typedef struct {
        char* new_name;
        Expr* init_val;
    } VarInfo;
    VarInfo* info = calloc(var_count, sizeof(VarInfo));

    for (size_t i = 0; i < var_count; i++) {
        Expr* v = vars->data.function.args[i];
        const char* orig_name = NULL;
        Expr* init_val = NULL;

        if (v->type == EXPR_SYMBOL) {
            orig_name = v->data.symbol.name;
        } else if (v->type == EXPR_FUNCTION && v->data.function.head->data.symbol.name == SYM_Set && v->data.function.arg_count == 2) {
            if (v->data.function.args[0]->type == EXPR_SYMBOL) {
                orig_name = v->data.function.args[0]->data.symbol.name;
                init_val = evaluate(v->data.function.args[1]);
            }
        }

        if (orig_name) {
            char buf[256];
            snprintf(buf, sizeof(buf), "%s$%lld", orig_name, (long long)current_mn);
            info[i].new_name = mathilda_strdup(buf);
            info[i].init_val = init_val;

            ScopingEnv* new_node = malloc(sizeof(ScopingEnv));
            new_node->old_name = orig_name;
            new_node->replacement = expr_new_symbol(buf);
            new_node->next = env;
            env = new_node;

            // Register temporary symbol
            symtab_get_def(buf)->attributes |= ATTR_TEMPORARY;
            if (init_val) {
                symtab_add_own_value(buf, new_node->replacement, init_val);
            }
        }
    }

    Expr* substituted_body = substitute_scoping(body, env);
    Expr* final_res = evaluate(substituted_body);
    expr_free(substituted_body);

    /* Trap Return targeting this Module. Return[v] (1-arg) is consumed
     * unconditionally; Return[v, h] is consumed only when h == Module,
     * else the marker is handed upward unchanged so an enclosing
     * boundary with the matching head can claim it. */
    {
        Expr* rv = NULL;
        EvalReturnAction ra = eval_classify_return(final_res, SYM_Module, &rv);
        if (ra == EVAL_RETURN_CONSUME) {
            expr_free(final_res);
            final_res = rv;
        }
        /* PROPAGATE / NONE: final_res is returned unchanged. */
    }

    // Cleanup env and info
    while (env) {
        ScopingEnv* next = env->next;
        expr_free(env->replacement);
        free(env);
        env = next;
    }
    for (size_t i = 0; i < var_count; i++) {
        if (info[i].new_name) free(info[i].new_name);
        if (info[i].init_val) expr_free(info[i].init_val);
    }
    free(info);

    return final_res;
}

/* ---------------------------------------------------------------- Block ----
 * Block[{x, y = v}, body] gives its locals DYNAMIC scope: the body sees them
 * with no value, and whatever they had is restored on the way out. "Whatever
 * they had" is BOTH rule lists, not just OwnValues:
 *
 *     gg[a_] := "orig";
 *     Block[{gg}, gg[a_] := "patched"; gg[1]]      (* "patched"            *)
 *     gg[1]                                        (* "orig" -- restored   *)
 *
 * Restoring only own_values (as this did before) left the DownValues written
 * inside the Block installed for the rest of the session. That is not a missing
 * nicety: Block over a function symbol is the standard way to install a
 * temporary rewrite hook (ParallelMixedSpecial's ExtendedBounds does exactly
 * this over four of ParallelMixed's bound-decision symbols), so the first such
 * call would silently and permanently repoint them.
 *
 * DefaultValues (`default_options`, backing Options[f]) are deliberately NOT
 * cleared: Mathematica's Block does clear them, but Options here are registered
 * once at module-init time and nothing in the tree rebinds them under Block, so
 * clearing would be pure blast radius.
 *
 * Non-local exit. Throw/Return travel as ordinary sentinel return values, so
 * they reach the restore below normally. A TimeConstrained timeout does not --
 * it siglongjmps clean past every C frame between the deadline and
 * tc_run_guarded. Frames are therefore also threaded on a global stack so that
 * unwind can drain whatever it jumped over; see mth_block_depth_save /
 * mth_block_depth_unwind, used by tc_run_guarded exactly as it already saves
 * and restores the async-defer count and the message-suppression depth. */

typedef struct BlockSavedVar {
    char*    name;
    Rule*    old_own;
    Rule*    old_down;
#if UP_VALUES
    Rule*    old_up;
#endif
    uint32_t old_attrs;
} BlockSavedVar;

typedef struct BlockFrame {
    BlockSavedVar*     saved;
    size_t             count;
    struct BlockFrame* next;
} BlockFrame;

static BlockFrame* block_stack = NULL;
static int         block_depth = 0;

/* Free one rule list wholesale (the values written inside the Block). */
static void blk_free_rules(Rule* r) {
    while (r) {
        Rule* next = r->next;
        expr_free(r->pattern);
        expr_free(r->replacement);
        free(r);
        r = next;
    }
}

/* Put every saved symbol back and release the frame. Idempotent per frame:
 * only ever called once, either by builtin_block or by the unwind drain. */
static void blk_restore_frame(BlockFrame* f) {
    for (size_t i = 0; i < f->count; i++) {
        if (!f->saved[i].name) continue;
        SymbolDef* def = symtab_get_def(f->saved[i].name);
        blk_free_rules(def->own_values);
        blk_free_rules(def->down_values);
        def->own_values  = f->saved[i].old_own;
        def->down_values = f->saved[i].old_down;
#if UP_VALUES
        /* Free the upvalues written inside the Block, decrementing the global
         * count, then restore the saved ones (which stayed counted while shadowed). */
        { size_t k = 0; for (Rule* r = def->up_values; r; r = r->next) k++;
          blk_free_rules(def->up_values);
          symtab_up_value_count = (symtab_up_value_count >= k) ? symtab_up_value_count - k : 0; }
        def->up_values = f->saved[i].old_up;
#endif
        def->attributes  = f->saved[i].old_attrs;
        free(f->saved[i].name);
    }
    free(f->saved);
    free(f);
}

int mth_block_depth_save(void) { return block_depth; }

/* Drain the Block stack back down to `depth`, restoring each frame. Called
 * after a timeout siglongjmp, where the frames between here and the deadline
 * never got to run their own restore. */
void mth_block_depth_unwind(int depth) {
    while (block_depth > depth && block_stack) {
        BlockFrame* f = block_stack;
        block_stack = f->next;
        block_depth--;
        blk_restore_frame(f);
    }
    /* Defensive: a mismatched count must not leave the counter drifting. */
    if (!block_stack) block_depth = 0;
}

Expr* builtin_block(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* vars = res->data.function.args[0];
    Expr* body = res->data.function.args[1];

    if (vars->type != EXPR_FUNCTION || vars->data.function.head->data.symbol.name != SYM_List) return NULL;

    size_t var_count = vars->data.function.arg_count;

    BlockFrame* frame = malloc(sizeof(BlockFrame));
    frame->saved = calloc(var_count ? var_count : 1, sizeof(BlockSavedVar));
    frame->count = var_count;
    frame->next  = block_stack;
    block_stack  = frame;
    block_depth++;

    for (size_t i = 0; i < var_count; i++) {
        Expr* v = vars->data.function.args[i];
        const char* name = NULL;
        Expr* init_val = NULL;

        if (v->type == EXPR_SYMBOL) {
            name = v->data.symbol.name;
        } else if (v->type == EXPR_FUNCTION && v->data.function.head->data.symbol.name == SYM_Set && v->data.function.arg_count == 2) {
            if (v->data.function.args[0]->type == EXPR_SYMBOL) {
                name = v->data.function.args[0]->data.symbol.name;
                init_val = evaluate(v->data.function.args[1]);
            }
        }

        if (name) {
            SymbolDef* def = symtab_get_def(name);
            frame->saved[i].name      = mathilda_strdup(name);
            frame->saved[i].old_own   = def->own_values;
            frame->saved[i].old_down  = def->down_values;
            frame->saved[i].old_attrs = def->attributes;
            def->own_values  = NULL;   /* the body sees the symbol unset ... */
            def->down_values = NULL;   /* ... in BOTH rule lists */
#if UP_VALUES
            /* Shadow upvalues too. The detached list stays counted (it is only
             * homeless, not freed) -- the count gate is conservative, so an
             * over-count merely skips the early-out, never a real rule. */
            frame->saved[i].old_up = def->up_values;
            def->up_values = NULL;
#endif

            if (init_val) {
                symtab_add_own_value(name, (v->type == EXPR_SYMBOL ? v : v->data.function.args[0]), init_val);
                expr_free(init_val);
            }
        }
    }

    Expr* final_res = evaluate(body);

    /* Trap Return targeting this Block. Symmetric to the Module path. */
    {
        Expr* rv = NULL;
        EvalReturnAction ra = eval_classify_return(final_res, SYM_Block, &rv);
        if (ra == EVAL_RETURN_CONSUME) {
            expr_free(final_res);
            final_res = rv;
        }
    }

    /* Pop this frame (an inner Block always pops before its outer one, so the
     * top of the stack is ours) and restore. */
    if (block_stack == frame) {
        block_stack = frame->next;
        block_depth--;
        blk_restore_frame(frame);
    }

    return final_res;
}

Expr* builtin_with(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* vars = res->data.function.args[0];
    Expr* body = res->data.function.args[1];

    if (vars->type != EXPR_FUNCTION || vars->data.function.head->data.symbol.name != SYM_List) return NULL;

    size_t var_count = vars->data.function.arg_count;
    ScopingEnv* env = NULL;

    for (size_t i = 0; i < var_count; i++) {
        Expr* v = vars->data.function.args[i];
        const char* name = NULL;
        Expr* val = NULL;

        if (v->type == EXPR_FUNCTION && v->data.function.arg_count == 2) {
            const char* h = v->data.function.head->data.symbol.name;
            if (h == SYM_Set || h == SYM_SetDelayed) {
                if (v->data.function.args[0]->type == EXPR_SYMBOL) {
                    name = v->data.function.args[0]->data.symbol.name;
                    if (h == SYM_Set) {
                        val = evaluate(v->data.function.args[1]);
                    } else {
                        val = expr_copy(v->data.function.args[1]);
                    }
                }
            }
        }

        if (name && val) {
            ScopingEnv* new_node = malloc(sizeof(ScopingEnv));
            new_node->old_name = name;
            new_node->replacement = val;
            new_node->next = env;
            env = new_node;
        }
    }

    Expr* substituted_body = substitute_scoping(body, env);
    Expr* final_res = evaluate(substituted_body);
    expr_free(substituted_body);

    /* Trap Return targeting this With. Symmetric to Module / Block. */
    {
        Expr* rv = NULL;
        EvalReturnAction ra = eval_classify_return(final_res, SYM_With, &rv);
        if (ra == EVAL_RETURN_CONSUME) {
            expr_free(final_res);
            final_res = rv;
        }
    }

    while (env) {
        ScopingEnv* next = env->next;
        expr_free(env->replacement);
        free(env);
        env = next;
    }

    return final_res;
}
