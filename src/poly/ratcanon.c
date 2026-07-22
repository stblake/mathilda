/* ratcanon.c — unified rational normalization (see RATCANON_REWRITE_PLAN.md).
 *
 * PHASE 1 prototype only.  Proves the one-front-end + one-reduction pipeline on
 * the four representative regimes.  Throwaway — replaced by rat_canon_build /
 * rat_canon_reduce in Phases 2-3.
 */
#include "ratcanon.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "sym_names.h"
#include "arithmetic.h"
#include "flint_bridge.h"
#include "core.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/* ---- generator substitution map ---------------------------------------- */

typedef struct { Expr* kernel; char* name; } RcpGen;   /* kernel <-> $rcpN$ */
typedef struct {
    RcpGen* g; size_t n, cap;
    int ctr;
    const char* iname;   /* shared name for the imaginary unit, or NULL */
} RcpMap;

static void rcp_map_free(RcpMap* m) {
    for (size_t i = 0; i < m->n; i++) { expr_free(m->g[i].kernel); free(m->g[i].name); }
    free(m->g);
}

/* Return the fresh symbol name assigned to kernel `k` (dedup by expr_eq),
 * creating it on first sight.  `k` is copied. */
static const char* rcp_bind(RcpMap* m, const Expr* k) {
    for (size_t i = 0; i < m->n; i++)
        if (expr_eq(m->g[i].kernel, (Expr*)k)) return m->g[i].name;
    if (m->n == m->cap) {
        m->cap = m->cap ? m->cap * 2 : 8;
        m->g = realloc(m->g, m->cap * sizeof(RcpGen));
    }
    char buf[32];
    snprintf(buf, sizeof buf, "$rcp%d$", m->ctr++);
    char* nm = malloc(strlen(buf) + 1);
    strcpy(nm, buf);
    m->g[m->n].kernel = expr_copy((Expr*)k);
    m->g[m->n].name = nm;
    return m->g[m->n++].name;
}

/* Shared fresh symbol for the imaginary unit (kernel = I, so map-back + eval
 * applies I^2 -> -1). */
static const char* rcp_i_gen(RcpMap* m) {
    if (m->iname) return m->iname;
    Expr* i = expr_new_symbol("I");
    m->iname = rcp_bind(m, i);
    expr_free(i);
    return m->iname;
}

static int rcp_is_number(const Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return 1;
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Rational) return 1;
    return 0;
}

static int rcp_is_indep_kernel_head(const char* h) {
    static const char* const ks[] = {
        "Log", "ArcSin","ArcCos","ArcTan","ArcCot","ArcSec","ArcCsc",
        "ArcSinh","ArcCosh","ArcTanh","ArcCoth","ArcSech","ArcCsch", NULL };
    for (int i = 0; ks[i]; i++) if (strcmp(h, ks[i]) == 0) return 1;
    return 0;
}

/* Substitute every kernel / algebraic constant / radical to a fresh free symbol.
 * Kernels are treated as opaque atoms (arguments not descended into — Phase 2
 * adds recursive argument normalization). */
static Expr* rcp_forward(const Expr* e, RcpMap* m) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    const Expr* h = e->data.function.head;
    size_t ac = e->data.function.arg_count;
    Expr** av = e->data.function.args;
    if (h->type == EXPR_SYMBOL) {
        const char* hn = h->data.symbol.name;
        /* Complex[a,b] -> a + b*Igen */
        if (hn == SYM_Complex && ac == 2) {
            const char* ig = rcp_i_gen(m);
            return expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){
                expr_copy(av[0]),
                expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
                    expr_copy(av[1]), expr_new_symbol(ig) }, 2) }, 2);
        }
        /* Sqrt[symbolic] / symbolic^(p/q) -> fresh algebraic symbol */
        if (hn == SYM_Sqrt && ac == 1 && !rcp_is_number(av[0]))
            return expr_new_symbol(rcp_bind(m, e));
        if (hn == SYM_Power && ac == 2) {
            int64_t p, q;
            if (is_rational(av[1], &p, &q) && q != 1 && !rcp_is_number(av[0]))
                return expr_new_symbol(rcp_bind(m, e));
        }
        /* Log / inverse-trig -> fresh transcendental symbol */
        if (rcp_is_indep_kernel_head(hn) && ac == 1)
            return expr_new_symbol(rcp_bind(m, e));
    }
    /* recurse */
    Expr* nh = rcp_forward(h, m);
    Expr** na = malloc(sizeof(Expr*) * (ac ? ac : 1));
    for (size_t i = 0; i < ac; i++) na[i] = rcp_forward(av[i], m);
    Expr* r = expr_new_function(nh, na, ac);
    free(na);
    return r;
}

/* Substitute the fresh $rcpN$ symbols back to their kernels. */
static Expr* rcp_backward(const Expr* e, RcpMap* m) {
    if (!e) return NULL;
    if (e->type == EXPR_SYMBOL) {
        for (size_t i = 0; i < m->n; i++)
            if (strcmp(e->data.symbol.name, m->g[i].name) == 0)
                return expr_copy(m->g[i].kernel);
        return expr_copy((Expr*)e);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t ac = e->data.function.arg_count;
    Expr* nh = rcp_backward(e->data.function.head, m);
    Expr** na = malloc(sizeof(Expr*) * (ac ? ac : 1));
    for (size_t i = 0; i < ac; i++) na[i] = rcp_backward(e->data.function.args[i], m);
    Expr* r = expr_new_function(nh, na, ac);
    free(na);
    return r;
}

/* RatCanonPrototype[expr]: substitute all kernels -> fresh free symbols, reduce
 * over Q with the existing plain-Q engine, map back, eval to apply the algebraic
 * relations.  Returns NULL (unevaluated) if the reduction declines. */
static Expr* builtin_ratcanon_prototype(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    RcpMap m; memset(&m, 0, sizeof m);
    Expr* sub = rcp_forward(arg, &m);

    Expr* reduced = flint_rational_together(sub);   /* one reduction over Q */
    expr_free(sub);
    if (!reduced) { rcp_map_free(&m); return NULL; }

    Expr* back = rcp_backward(reduced, &m);
    expr_free(reduced);
    rcp_map_free(&m);
    return eval_and_free(back);   /* applies I^2->-1, Sqrt[k]^2->k, ... */
}

void ratcanon_init(void) {
    symtab_add_builtin("RatCanonPrototype", builtin_ratcanon_prototype);
    symtab_get_def("RatCanonPrototype")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("RatCanonPrototype",
        "RatCanonPrototype[expr] (Phase-1 prototype) reduces a rational function "
        "over the differential/algebraic tower of expr via one FLINT reduction.");
}
