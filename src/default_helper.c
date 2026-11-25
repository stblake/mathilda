#include "expr.h"
#include "eval.h"
#include "sym_names.h"
#include "symtab.h"
#include <string.h>

Expr* get_default_value(Expr* pat_head, int pos, int total_pats) {
    if (!pat_head || pat_head->type != EXPR_SYMBOL) return NULL;
    const char* hn = pat_head->data.symbol.name;

    /* Fast path: when no user-defined Default rules exist, the three
     * evaluate(Default[...]) calls below would all return the constructed
     * expression unchanged, so we'd just fall through to the built-in
     * fallback anyway. Skip the evaluations entirely. This matters
     * because the matcher invokes get_default_value once per Optional
     * pattern (`x_.`) attempt, and each evaluate call is non-trivial. */
    /* The interned-pointer cache, not a by-name lookup: this runs per attempt.
     * Reading Default's rules (or their absence) from C is a dependency of a
     * failed match the evaluator cannot see, so note it. */
    SymbolDef* def_sym = intern_symbol_node(SYM_Default);
    eval_dep_note(def_sym);
    if (!def_sym || def_sym->down_values == NULL) {
        if (hn == SYM_Plus) return expr_new_integer(0);
        if (hn == SYM_Times || hn == SYM_Power) return expr_new_integer(1);
        return NULL;
    }

    // Construct Default[f, pos, total_pats]
    Expr* args3[3] = { expr_copy(pat_head), expr_new_integer(pos), expr_new_integer(total_pats) };
    Expr* def3 = expr_new_function(expr_new_symbol(SYM_Default), args3, 3);
    Expr* res3 = evaluate(def3);
    if (!expr_eq(res3, def3)) {
        expr_free(def3);
        return res3;
    }
    expr_free(res3);
    expr_free(def3);

    // Construct Default[f, pos]
    Expr* args2[2] = { expr_copy(pat_head), expr_new_integer(pos) };
    Expr* def2 = expr_new_function(expr_new_symbol(SYM_Default), args2, 2);
    Expr* res2 = evaluate(def2);
    if (!expr_eq(res2, def2)) {
        expr_free(def2);
        return res2;
    }
    expr_free(res2);
    expr_free(def2);

    // Construct Default[f]
    Expr* args1[1] = { expr_copy(pat_head) };
    Expr* def1 = expr_new_function(expr_new_symbol(SYM_Default), args1, 1);
    Expr* res1 = evaluate(def1);
    if (!expr_eq(res1, def1)) {
        expr_free(def1);
        return res1;
    }
    expr_free(res1);
    expr_free(def1);

    // Built-in fallbacks
    if (hn == SYM_Plus) return expr_new_integer(0);
    if (hn == SYM_Times || hn == SYM_Power) return expr_new_integer(1);

    return NULL;
}

/* Default[f] / Default[f, i]: the built-in defaults the matcher uses for an
 * optional argument x_. -- Default[Plus] (any position) is 0, Default[Times]
 * is 1, Default[Power, 2] is 1 (audit Default#1).  A user's Default[f, ...]
 * values are DownValues and apply for every other f; NULL leaves the call. */
Expr* builtin_default(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    Expr* f = res->data.function.args[0];
    if (f->type != EXPR_SYMBOL) return NULL;
    if (argc == 2 && res->data.function.args[1]->type != EXPR_INTEGER) return NULL;
    const char* h = f->data.symbol.name;
    if (h == SYM_Plus) return expr_new_integer(0);
    if (h == SYM_Times) return expr_new_integer(1);
    if (h == SYM_Power && argc == 2 && res->data.function.args[1]->data.integer == 2)
        return expr_new_integer(1);
    return NULL;
}
