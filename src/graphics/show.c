/* show.c — Show[graphics, opts...] plus the no-Raylib graphics_show()
 * stub. When USE_GRAPHICS is compiled in, render.c provides the real
 * graphics_show(); this file's stub is excluded then (see the #ifndef
 * below) so there's exactly one definition either way. */

#include "show.h"
#include "sym_names.h"
#include <stdlib.h>

static bool is_rule_arg(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return false;
    return (h->data.symbol == SYM_Rule || h->data.symbol == SYM_RuleDelayed)
        && e->data.function.arg_count == 2;
}

static bool is_graphics(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Graphics
        && e->data.function.arg_count >= 1;
}

Expr* builtin_show(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;

    Expr* g = res->data.function.args[0];
    if (!is_graphics(g)) return NULL;

    for (size_t i = 1; i < argc; i++) {
        if (!is_rule_arg(res->data.function.args[i])) return NULL;
    }

    if (argc == 1) {
        Expr* copy = expr_copy(g);
        graphics_show(copy);
        return copy;
    }

    /* Merge: start from g's own options, then let each Show[] option
     * override (by option name) or append, matching Mathematica's
     * rightmost/later-wins option semantics. */
    size_t gargc = g->data.function.arg_count;
    size_t max_opts = (gargc - 1) + (argc - 1);
    Expr** opts = malloc(sizeof(Expr*) * (max_opts > 0 ? max_opts : 1));
    size_t n = 0;
    for (size_t i = 1; i < gargc; i++) opts[n++] = expr_copy(g->data.function.args[i]);

    for (size_t i = 1; i < argc; i++) {
        Expr* newopt = res->data.function.args[i];
        Expr* newlhs = newopt->data.function.args[0];
        bool replaced = false;
        for (size_t j = 0; j < n; j++) {
            Expr* lhs = opts[j]->data.function.args[0];
            if (lhs->type == EXPR_SYMBOL && newlhs->type == EXPR_SYMBOL
                && lhs->data.symbol == newlhs->data.symbol) {
                expr_free(opts[j]);
                opts[j] = expr_copy(newopt);
                replaced = true;
                break;
            }
        }
        if (!replaced) opts[n++] = expr_copy(newopt);
    }

    Expr** gargs = malloc(sizeof(Expr*) * (1 + n));
    gargs[0] = expr_copy(g->data.function.args[0]);
    for (size_t i = 0; i < n; i++) gargs[1 + i] = opts[i];
    free(opts);

    Expr* merged = expr_new_function(expr_new_symbol(SYM_Graphics), gargs, 1 + n);
    free(gargs);

    graphics_show(merged);
    return merged;
}

#ifndef USE_GRAPHICS
#include <stdio.h>
void graphics_show(const Expr* graphics_expr) {
    (void)graphics_expr;
    printf("Show: graphics support not compiled in (install raylib and rebuild).\n");
}
#endif
