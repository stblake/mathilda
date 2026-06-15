#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
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

static bool is_scoping_construct(Expr* e) {
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol;
    return h == SYM_Module || h == SYM_Block || h == SYM_With || 
           h == SYM_Function || h == SYM_Table;
}

static Expr* substitute_scoping(Expr* e, ScopingEnv* env) {
    if (!e) return NULL;
    if (e->type == EXPR_SYMBOL) {
        ScopingEnv* curr = env;
        while (curr) {
            if (strcmp(e->data.symbol, curr->old_name) == 0) {
                return expr_copy(curr->replacement);
            }
            curr = curr->next;
        }
        return expr_copy(e);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy(e);

    // Handle shadowing in scoping constructs
    ScopingEnv* filtered_env = env;
    
    if (is_scoping_construct(e) && e->data.function.arg_count >= 1) {
        Expr* vars = e->data.function.args[0];
        if (vars->type == EXPR_FUNCTION && vars->data.function.head->data.symbol == SYM_List) {
            // Create a filtered env that doesn't contain variables redefined here
            for (size_t i = 0; i < vars->data.function.arg_count; i++) {
                Expr* v = vars->data.function.args[i];
                const char* shadowed_name = NULL;
                if (v->type == EXPR_SYMBOL) shadowed_name = v->data.symbol;
                else if (v->type == EXPR_FUNCTION && v->data.function.head->data.symbol == SYM_Set && v->data.function.arg_count == 2) {
                    if (v->data.function.args[0]->type == EXPR_SYMBOL) shadowed_name = v->data.function.args[0]->data.symbol;
                }

                if (shadowed_name) {
                    // This variable is shadowed, so we remove it from the env we pass down
                    // For simplicity, we just rebuild the env list skipping shadowed names
                    ScopingEnv* new_env = NULL;
                    ScopingEnv* curr = filtered_env;
                    while (curr) {
                        if (strcmp(curr->old_name, shadowed_name) != 0) {
                            ScopingEnv* node = malloc(sizeof(ScopingEnv));
                            node->old_name = curr->old_name;
                            node->replacement = curr->replacement;
                            node->next = new_env;
                            new_env = node;
                        }
                        curr = curr->next;
                    }
                    if (filtered_env != env) {
                        // Free the intermediate filtered env
                        ScopingEnv* tmp = filtered_env;
                        while (tmp) {
                            ScopingEnv* next = tmp->next;
                            free(tmp);
                            tmp = next;
                        }
                    }
                    filtered_env = new_env;
                }
            }
        }
    }
    
    Expr** new_args = malloc(sizeof(Expr*) * e->data.function.arg_count);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
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
        if (i == 0 && is_scoping_construct(e)) {
            Expr* vars_list = e->data.function.args[i];
            if (vars_list->type == EXPR_FUNCTION
                && vars_list->data.function.head
                && vars_list->data.function.head->type == EXPR_SYMBOL
                && vars_list->data.function.head->data.symbol == SYM_List) {
                size_t nb = vars_list->data.function.arg_count;
                Expr** new_b = malloc(sizeof(Expr*) * (nb > 0 ? nb : 1));
                for (size_t bi = 0; bi < nb; bi++) {
                    Expr* b = vars_list->data.function.args[bi];
                    if (b->type == EXPR_FUNCTION && b->data.function.head
                        && b->data.function.head->type == EXPR_SYMBOL
                        && (b->data.function.head->data.symbol == SYM_Set
                            || b->data.function.head->data.symbol == SYM_SetDelayed)
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

Expr* builtin_module(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* vars = res->data.function.args[0];
    Expr* body = res->data.function.args[1];

    if (vars->type != EXPR_FUNCTION || vars->data.function.head->data.symbol != SYM_List) return NULL;

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
            orig_name = v->data.symbol;
        } else if (v->type == EXPR_FUNCTION && v->data.function.head->data.symbol == SYM_Set && v->data.function.arg_count == 2) {
            if (v->data.function.args[0]->type == EXPR_SYMBOL) {
                orig_name = v->data.function.args[0]->data.symbol;
                init_val = evaluate(v->data.function.args[1]);
            }
        }

        if (orig_name) {
            char buf[256];
            snprintf(buf, sizeof(buf), "%s$%lld", orig_name, (long long)current_mn);
            info[i].new_name = strdup(buf);
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

Expr* builtin_block(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* vars = res->data.function.args[0];
    Expr* body = res->data.function.args[1];

    if (vars->type != EXPR_FUNCTION || vars->data.function.head->data.symbol != SYM_List) return NULL;

    size_t var_count = vars->data.function.arg_count;
    
    typedef struct {
        char* name;
        Rule* old_own;
        uint32_t old_attrs;
    } SavedVar;
    SavedVar* saved = calloc(var_count, sizeof(SavedVar));

    for (size_t i = 0; i < var_count; i++) {
        Expr* v = vars->data.function.args[i];
        const char* name = NULL;
        Expr* init_val = NULL;

        if (v->type == EXPR_SYMBOL) {
            name = v->data.symbol;
        } else if (v->type == EXPR_FUNCTION && v->data.function.head->data.symbol == SYM_Set && v->data.function.arg_count == 2) {
            if (v->data.function.args[0]->type == EXPR_SYMBOL) {
                name = v->data.function.args[0]->data.symbol;
                init_val = evaluate(v->data.function.args[1]);
            }
        }

        if (name) {
            SymbolDef* def = symtab_get_def(name);
            saved[i].name = strdup(name);
            saved[i].old_own = def->own_values;
            saved[i].old_attrs = def->attributes;
            def->own_values = NULL; // Clear values
            
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

    // Restore
    for (size_t i = 0; i < var_count; i++) {
        if (saved[i].name) {
            SymbolDef* def = symtab_get_def(saved[i].name);
            // Free temporary values assigned in block
            Rule* curr = def->own_values;
            while (curr) {
                Rule* next = curr->next;
                expr_free(curr->pattern);
                expr_free(curr->replacement);
                free(curr);
                curr = next;
            }
            def->own_values = saved[i].old_own;
            def->attributes = saved[i].old_attrs;
            free(saved[i].name);
        }
    }
    free(saved);

    return final_res;
}

Expr* builtin_with(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* vars = res->data.function.args[0];
    Expr* body = res->data.function.args[1];

    if (vars->type != EXPR_FUNCTION || vars->data.function.head->data.symbol != SYM_List) return NULL;

    size_t var_count = vars->data.function.arg_count;
    ScopingEnv* env = NULL;

    for (size_t i = 0; i < var_count; i++) {
        Expr* v = vars->data.function.args[i];
        const char* name = NULL;
        Expr* val = NULL;

        if (v->type == EXPR_FUNCTION && v->data.function.arg_count == 2) {
            const char* h = v->data.function.head->data.symbol;
            if (h == SYM_Set || h == SYM_SetDelayed) {
                if (v->data.function.args[0]->type == EXPR_SYMBOL) {
                    name = v->data.function.args[0]->data.symbol;
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
