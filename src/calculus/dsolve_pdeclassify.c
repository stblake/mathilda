/*
 * dsolve_pdeclassify.c — PDEClassify[eqn, u, {v1, v2}]: classify a second-order
 * linear PDE by the discriminant of its principal part.
 *
 * For  A u_{v1 v1} + B u_{v1 v2} + C u_{v2 v2} + (lower order) == 0  the type is
 * set by  Δ = B² − 4 A C  (only the highest-order terms matter):
 *
 *     Δ > 0  →  "Hyperbolic"   (e.g. the wave equation u_tt == c² u_xx)
 *     Δ = 0  →  "Parabolic"    (e.g. the heat equation u_t == u_xx)
 *     Δ < 0  →  "Elliptic"     (e.g. Laplace u_xx + u_yy == 0)
 *
 * A, B, C are read from the equation's second-order terms exactly as
 * DSolve`PDELinearSecondOrder reads them.  The sign of Δ is decided for a
 * definite constant discriminant; a variable-coefficient / parameter-dependent Δ
 * whose sign is not a decidable constant (a mixed-type equation such as the
 * Tricomi equation y u_xx + u_yy == 0) leaves the call unevaluated — an honest
 * decline rather than a region-blind label.  A first cut: linear in the
 * second-order terms, over two independent variables.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>
#include <stdint.h>

/* Derivative[o1,o2][u][v1,v2] */
static Expr* pdec_deriv(const char* u, int o1, int o2, const char* v1, const char* v2) {
    Expr* d = expr_new_function(expr_new_symbol(SYM_Derivative),
                  (Expr*[]){ expr_new_integer(o1), expr_new_integer(o2) }, 2);
    Expr* du = expr_new_function(d, (Expr*[]){ expr_new_symbol(u) }, 1);
    return expr_new_function(du, (Expr*[]){ expr_new_symbol(v1), expr_new_symbol(v2) }, 2);
}

/* Decide the sign of the (already Simplify-reduced) discriminant: +1 / -1 / 0,
 * or 2 when undecidable.  `disc` is borrowed. */
static int pdec_sign(const Expr* disc) {
    if (ds_is_zero(disc)) return 0;
    if (disc->type == EXPR_INTEGER) return disc->data.integer > 0 ? 1 : -1;
    if (disc->type == EXPR_REAL)    return disc->data.real > 0 ? 1 : -1;
    Expr* s = eval_and_free(ds_call1("Sign", expr_copy((Expr*)disc)));
    int r = 2;
    if (s->type == EXPR_INTEGER) {
        if (s->data.integer > 0) r = 1;
        else if (s->data.integer < 0) r = -1;
        else r = 0;
    }
    expr_free(s);
    return r;
}

/* Is f the unknown u or one of its derivatives Derivative[..][u]?  On a
 * derivative, *order receives the total order. */
static bool pdec_is_u_head(const Expr* h, const char* u, int64_t* order) {
    *order = 0;
    if (h->type == EXPR_SYMBOL) return h->data.symbol.name == u;
    if (h->type == EXPR_FUNCTION && h->data.function.arg_count == 1 &&
        h->data.function.args[0]->type == EXPR_SYMBOL &&
        h->data.function.args[0]->data.symbol.name == u &&
        h->data.function.head->type == EXPR_FUNCTION &&
        h->data.function.head->data.function.head->type == EXPR_SYMBOL &&
        h->data.function.head->data.function.head->data.symbol.name == SYM_Derivative) {
        const Expr* d = h->data.function.head;
        for (size_t i = 0; i < d->data.function.arg_count; i++) {
            if (d->data.function.args[i]->type != EXPR_INTEGER) { *order = INT64_MAX; return true; }
            *order += d->data.function.args[i]->data.integer;
        }
        return true;
    }
    return false;
}

/* Walk the residual: the highest derivative order of u, and whether u is
 * applied with its arguments in the order (v2, v1) -- the order the caller
 * listed the variables need not be u's own. */
static void pdec_scan(const Expr* e, const char* u, const char* v1, const char* v2,
                      int64_t* max_order, bool* swapped) {
    if (!e || e->type != EXPR_FUNCTION) return;
    int64_t ord;
    if (pdec_is_u_head(e->data.function.head, u, &ord)) {
        if (ord > *max_order) *max_order = ord;
        if (e->data.function.arg_count == 2 &&
            e->data.function.args[0]->type == EXPR_SYMBOL &&
            e->data.function.args[1]->type == EXPR_SYMBOL &&
            e->data.function.args[0]->data.symbol.name == v2 &&
            e->data.function.args[1]->data.symbol.name == v1)
            *swapped = true;
    }
    pdec_scan(e->data.function.head, u, v1, v2, max_order, swapped);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        pdec_scan(e->data.function.args[i], u, v1, v2, max_order, swapped);
}

static Expr* builtin_pdeclassify(Expr* res) {
    DSolveProblem P;
    if (!dsolve_parse(res, &P)) { dsolve_problem_free(&P); return NULL; }
    if (!P.is_pde || P.nfun != 1 || P.nind != 2 || P.neq != 1) {
        dsolve_problem_free(&P); return NULL;
    }
    const char* u  = P.fun_names[0];
    const char* v1 = P.ind_names[0];
    const char* v2 = P.ind_names[1];

    /* The discriminant classifies a SECOND-order equation: a third- or
     * fourth-order one is not typed by its second-order part (audit
     * PDEClassify#1).  Variables listed in another order than u's own
     * arguments are read in u's order -- the type is symmetric in the
     * swap (audit PDEClassify#2). */
    {
        int64_t max_order = 0;
        bool swapped = false;
        pdec_scan(P.eq_residuals[0], u, v1, v2, &max_order, &swapped);
        if (max_order > 2) { dsolve_problem_free(&P); return NULL; }
        if (swapped) { const char* t = v1; v1 = v2; v2 = t; }
    }

    const char* sxx = intern_symbol("DSolve`pdecUxx");
    const char* sxy = intern_symbol("DSolve`pdecUxy");
    const char* syy = intern_symbol("DSolve`pdecUyy");

    Expr* R = expr_copy(P.eq_residuals[0]);
    R = ds_subst(R, pdec_deriv(u, 2, 0, v1, v2), expr_new_symbol(sxx));
    R = ds_subst(R, pdec_deriv(u, 1, 1, v1, v2), expr_new_symbol(sxy));
    R = ds_subst(R, pdec_deriv(u, 0, 2, v1, v2), expr_new_symbol(syy));

    Expr* A = ds_d(expr_copy(R), expr_new_symbol(sxx));
    Expr* B = ds_d(expr_copy(R), expr_new_symbol(sxy));
    Expr* C = ds_d(expr_copy(R), expr_new_symbol(syy));
    expr_free(R);

    /* linear in the second-order terms, and genuinely second order */
    bool ok = ds_free_of(A, sxx) && ds_free_of(A, sxy) && ds_free_of(A, syy)
           && ds_free_of(B, sxx) && ds_free_of(B, sxy) && ds_free_of(B, syy)
           && ds_free_of(C, sxx) && ds_free_of(C, sxy) && ds_free_of(C, syy)
           && !(ds_is_zero(A) && ds_is_zero(B) && ds_is_zero(C));
    if (!ok) { expr_free(A); expr_free(B); expr_free(C); dsolve_problem_free(&P); return NULL; }

    /* Δ = B² − 4 A C */
    Expr* B2  = eval_and_free(ds_call2(SYM_Power, B, expr_new_integer(2)));
    Expr* AC4 = eval_and_free(ds_call2(SYM_Times, expr_new_integer(4),
                    eval_and_free(ds_call2(SYM_Times, A, C))));
    Expr* disc = ds_simplify(eval_and_free(ds_call2(SYM_Subtract, B2, AC4)));

    int sg = pdec_sign(disc);
    expr_free(disc);
    dsolve_problem_free(&P);

    if (sg == 1)  return expr_new_string("Hyperbolic");
    if (sg == -1) return expr_new_string("Elliptic");
    if (sg == 0)  return expr_new_string("Parabolic");
    return NULL;  /* undecidable sign (mixed type / parameter-dependent): decline */
}

void dsolve_pdeclassify_init(void) {
    symtab_add_builtin("PDEClassify", builtin_pdeclassify);
    symtab_get_def("PDEClassify")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("PDEClassify",
        "PDEClassify[eqn, u, {v1, v2}] classifies a second-order linear PDE by the "
        "discriminant Δ = B² − 4 A C of its principal part A u_{v1 v1} + B u_{v1 v2} "
        "+ C u_{v2 v2}: \"Hyperbolic\" (Δ > 0, e.g. the wave equation), \"Parabolic\" "
        "(Δ == 0, e.g. the heat equation), or \"Elliptic\" (Δ < 0, e.g. Laplace's "
        "equation). Only the highest-order terms determine the type. A discriminant "
        "whose sign is not a decidable constant (a mixed-type / parameter-dependent "
        "equation such as Tricomi's y u_xx + u_yy == 0) leaves the call unevaluated.");
}
