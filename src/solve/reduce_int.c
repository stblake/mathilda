/*
 * reduce_int.c
 *
 * Integer / Rationals domain for `Reduce` (REDUCE_PLAN.md, Phase 5).  See
 * reduce_int.h.  Reuses `Solve[..., dom]` for the heavy lifting and reformats
 * its solution list; falls back to bounded integer enumeration (reduce_univar_
 * integers) for the univariate inequality case Solve declines.
 */
#include "reduce_int.h"
#include "reduce_univar.h"

#include "eval.h"
#include "sym_names.h"

#include <stdlib.h>
#include <string.h>

static bool is_head(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == name;
}

/* Collect the distinct generated parameters C[k] appearing in `e`.  The Solve
 * Diophantine engine always emits the head `C`; a non-default
 * GeneratedParameters head is applied afterwards by rename_param_head. */
static void collect_params(const Expr* e, Expr*** list, int* n) {
    if (!e || e->type != EXPR_FUNCTION) return;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && strcmp(h->data.symbol.name, "C") == 0) {
        for (int i = 0; i < *n; i++) if (expr_eq((*list)[i], e)) return;
        *list = realloc(*list, (size_t)(*n + 1) * sizeof(Expr*));
        (*list)[(*n)++] = expr_copy((Expr*)e);
        return;
    }
    collect_params(h, list, n);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        collect_params(e->data.function.args[i], list, n);
}

/* Rewrite every generated parameter C[k] (head C with a single argument) to
 * head[k], returning a fresh tree; used to honor GeneratedParameters -> head.
 * A no-op tree copy when head is "C". */
static Expr* rename_param_head(const Expr* e, const char* head) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && strcmp(h->data.symbol.name, "C") == 0
        && e->data.function.arg_count == 1) {
        return expr_new_function(expr_new_symbol(head),
            (Expr*[]){ expr_copy(e->data.function.args[0]) }, 1);
    }
    size_t n = e->data.function.arg_count;
    Expr* nh = rename_param_head(h, head);
    Expr** args = malloc(n * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) args[i] = rename_param_head(e->data.function.args[i], head);
    Expr* r = expr_new_function(nh, args, n);
    free(args);
    return r;
}

/* One solution tuple (a List of Rule[var,val]) -> a conjunction Expr, or NULL
 * on an unexpected shape.  Sets *all_values when the tuple is empty ({{}}). */
static Expr* tuple_to_conj(const Expr* tuple, const Expr* dom, bool* all_values) {
    if (!is_head(tuple, SYM_List)) return NULL;
    size_t nr = tuple->data.function.arg_count;
    if (nr == 0) { *all_values = true; return NULL; }

    Expr** params = NULL; int npar = 0;
    Expr** eqs = malloc(nr * sizeof(Expr*));
    for (size_t i = 0; i < nr; i++) {
        Expr* rule = tuple->data.function.args[i];
        if (!(is_head(rule, SYM_Rule) || is_head(rule, SYM_RuleDelayed))
            || rule->data.function.arg_count != 2) {
            for (size_t j = 0; j < i; j++) expr_free(eqs[j]);
            for (int j = 0; j < npar; j++) expr_free(params[j]);
            free(eqs); free(params);
            return NULL;
        }
        Expr* var = rule->data.function.args[0];
        Expr* val = rule->data.function.args[1];
        collect_params(val, &params, &npar);
        eqs[i] = expr_new_function(expr_new_symbol(SYM_Equal),
            (Expr*[]){ expr_copy(var), expr_copy(val) }, 2);
    }

    /* Element[C[k], dom] conjuncts first, then the var == val atoms. */
    int total = npar + (int)nr;
    Expr** parts = malloc((size_t)total * sizeof(Expr*));
    int p = 0;
    for (int i = 0; i < npar; i++)
        parts[p++] = expr_new_function(expr_new_symbol(SYM_Element),
            (Expr*[]){ params[i], expr_copy((Expr*)dom) }, 2);
    for (size_t i = 0; i < nr; i++) parts[p++] = eqs[i];
    free(eqs); free(params);

    Expr* conj = (total == 1) ? parts[0]
        : expr_new_function(expr_new_symbol(SYM_And), parts, (size_t)total);
    free(parts);
    return conj;
}

/* Reformat a Solve solution list into a logical Expr, or NULL to decline. */
static Expr* reformat_solutions(const Expr* sols, const Expr* dom) {
    if (!is_head(sols, SYM_List)) return NULL;   /* Solve declined */
    size_t nt = sols->data.function.arg_count;
    if (nt == 0) return expr_new_symbol(SYM_False);

    Expr** disj = malloc(nt * sizeof(Expr*));
    int nd = 0;
    for (size_t i = 0; i < nt; i++) {
        bool all_values = false;
        Expr* conj = tuple_to_conj(sols->data.function.args[i], dom,
                                   &all_values);
        if (all_values) {                        /* {{}} -> every value */
            for (int j = 0; j < nd; j++) expr_free(disj[j]);
            free(disj);
            return expr_new_symbol(SYM_True);
        }
        if (!conj) {                             /* unexpected shape -> decline */
            for (int j = 0; j < nd; j++) expr_free(disj[j]);
            free(disj);
            return NULL;
        }
        disj[nd++] = conj;
    }

    Expr* out = (nd == 1) ? disj[0]
        : expr_new_function(expr_new_symbol(SYM_Or), disj, (size_t)nd);
    free(disj);
    return eval_and_free(out);
}

Expr* reduce_integers(const Expr* expr, const Expr* vars_expr, const RForm* F,
                      Expr** vars, int nv, const Expr* dom,
                      const ReduceOpts* opts) {
    const char* head = (opts && opts->param_head) ? opts->param_head : "C";
    /* Reuse the Solve integer/rational engine (which always emits the parameter
     * head `C`), then rename C[k] -> head[k] to honor GeneratedParameters. */
    Expr* call = expr_new_function(expr_new_symbol(SYM_Solve),
        (Expr*[]){ expr_copy((Expr*)expr), expr_copy((Expr*)vars_expr),
                   expr_copy((Expr*)dom) }, 3);
    Expr* sols = eval_and_free(call);
    Expr* out = reformat_solutions(sols, dom);
    expr_free(sols);
    if (out && strcmp(head, "C") != 0) {
        Expr* renamed = rename_param_head(out, head);
        expr_free(out);
        out = renamed;
    }
    if (out) return out;

    /* Solve declined (e.g. a pure bounded inequality).  Univariate integers get
     * the bounded-enumeration fallback. */
    if (nv == 1 && dom->type == EXPR_SYMBOL && dom->data.symbol.name == SYM_Integers)
        return reduce_univar_integers(F, vars[0], vars, nv, opts);
    return NULL;
}

Expr* reduce_modular(const Expr* expr, const Expr* vars_expr,
                     const ReduceOpts* opts) {
    if (!opts || !opts->modulus) return NULL;
    /* Solve[expr, vars, Modulus -> p].  Solve's modular pre-pass enumerates
     * residues (and declines for a symbolic / out-of-range modulus or a
     * non-modular statement, leaving the result unevaluated). */
    Expr* modrule = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){ expr_new_symbol(SYM_Modulus), expr_copy((Expr*)opts->modulus) }, 2);
    Expr* call = expr_new_function(expr_new_symbol(SYM_Solve),
        (Expr*[]){ expr_copy((Expr*)expr), expr_copy((Expr*)vars_expr), modrule }, 3);
    Expr* sols = eval_and_free(call);
    /* Modular solutions carry no generated parameters, so `dom` is unused by
     * tuple_to_conj here; an empty list -> False, a non-List (Solve declined)
     * -> NULL (unevaluated). */
    Expr* out = reformat_solutions(sols, NULL);
    expr_free(sols);
    return out;
}
