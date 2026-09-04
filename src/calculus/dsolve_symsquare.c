/*
 * dsolve_symsquare.c — DSolve`SymmetricSquare.
 *
 * Solves a third-order linear homogeneous ODE whose solution space is the
 * symmetric square of a second-order equation: if u1, u2 span the solutions of
 * u'' = q u (normal form, no u' term), then {u1^2, u1 u2, u2^2} span a
 * third-order equation.  Concretely, after depressing y''' + p2 y'' + p1 y' +
 * p0 y == 0 to w''' + P1 w' + P0 w == 0 via y = w Exp[-Integrate[p2/3, x]], the
 * symmetric-square signature is
 *     P0 == P1'/2,   with underlying   u'' + (P1/4) u == 0.
 * We solve that second-order equation (Airy / Bessel / ... via the existing
 * recognizers), build w = C1 u1^2 + C2 u1 u2 + C3 u2^2, and undo the depression.
 *
 * Depression coefficients (closed form, avoiding an Exp[-.]/Exp[.] cancellation):
 *     P1 = p1 - p2^2/3 - p2'
 *     P0 = p0 - p1 p2/3 + 2 p2^3/27 - p2''/3
 *
 * Examples:
 *   y''' - 4(x+2)y' - 2y == 0            -> Airy:   u''=(x+2)u,
 *       y = C1 AiryAi[x+2]^2 + C2 AiryAi[x+2]AiryBi[x+2] + C3 AiryBi[x+2]^2.
 *   x^3 y''' + 3x^2 y'' + (4x^3-11x)y' + 4x^2 y == 0  -> Bessel (order sqrt 3).
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include <stdlib.h>

/* Extract the RHS of {{u[x] -> expr}} (applied form). */
static Expr* ss_extract_applied(Expr* r, const char* ufun) {
    if (!head_is(r, SYM_List) || r->data.function.arg_count < 1) return NULL;
    Expr* inner = r->data.function.args[0];
    if (!head_is(inner, SYM_List)) return NULL;
    for (size_t k = 0; k < inner->data.function.arg_count; k++) {
        Expr* rule = inner->data.function.args[k];
        if (head_is(rule, SYM_Rule) && rule->data.function.arg_count == 2) {
            Expr* lhs = rule->data.function.args[0];
            if (lhs->type == EXPR_FUNCTION && lhs->data.function.head->type == EXPR_SYMBOL
                && lhs->data.function.head->data.symbol.name == ufun)
                return expr_copy(rule->data.function.args[1]);
        }
    }
    return NULL;
}

/* a/b (simplified); a, b consumed. */
static Expr* ss_div(Expr* a, Expr* b) {
    return ds_simplify(ds_call2(SYM_Times, a,
        expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ b, expr_new_integer(-1) }, 2)));
}
/* n*x (n integer literal); x consumed. */
static Expr* ss_int_times(long k, Expr* x) { return ds_call2(SYM_Times, expr_new_integer(k), x); }
/* Coefficient[e, v]; e borrowed, v consumed. */
static Expr* ss_coeff(const Expr* e, Expr* v) {
    return eval_and_free(ds_call2("Coefficient", expr_copy((Expr*)e), v));
}

Expr** dsolve_symsquare_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 3) return NULL;
    const char* xvar = P->ind_names[0];

    Expr** a = NULL; Expr* g = NULL; int n = 0;
    if (!dsolve_linear_coeffs(P, &a, &g, &n)) return NULL;
    bool ok = (n == 3) && !ds_is_zero(a[3]) && ds_is_zero(g);   /* linear, order 3, homogeneous */
    if (!ok) { for (int k = 0; k <= n; k++) expr_free(a[k]); free(a); if (g) expr_free(g); return NULL; }
    expr_free(g);

    /* normalise: p_k = a_k / a_3 */
    Expr* p2 = ss_div(expr_copy(a[2]), expr_copy(a[3]));
    Expr* p1 = ss_div(expr_copy(a[1]), expr_copy(a[3]));
    Expr* p0 = ss_div(expr_copy(a[0]), expr_copy(a[3]));
    for (int k = 0; k <= 3; k++) expr_free(a[k]);
    free(a);

    Expr* p2p = ds_d(expr_copy(p2), expr_new_symbol(xvar));       /* p2'  */
    Expr* p2pp = ds_d(expr_copy(p2p), expr_new_symbol(xvar));     /* p2'' */

    /* P1 = p1 - p2^2/3 - p2' */
    Expr* P1 = ds_simplify(ds_call2(SYM_Subtract,
                   ds_call2(SYM_Subtract, expr_copy(p1),
                       ss_div(ds_call2(SYM_Power, expr_copy(p2), expr_new_integer(2)), expr_new_integer(3))),
                   expr_copy(p2p)));
    /* P0 = p0 - p1 p2/3 + 2 p2^3/27 - p2''/3 */
    Expr* P0 = ds_simplify(ds_call2(SYM_Plus,
                   ds_call2(SYM_Subtract, expr_copy(p0),
                       ss_div(ds_call2(SYM_Times, expr_copy(p1), expr_copy(p2)), expr_new_integer(3))),
                   ds_call2(SYM_Subtract,
                       ss_div(ss_int_times(2, ds_call2(SYM_Power, expr_copy(p2), expr_new_integer(3))),
                              expr_new_integer(27)),
                       ss_div(expr_copy(p2pp), expr_new_integer(3)))));

    /* symmetric-square condition: P0 == P1'/2 */
    Expr* P1p = ds_d(expr_copy(P1), expr_new_symbol(xvar));
    Expr* cond = ds_simplify(ds_call2(SYM_Subtract, expr_copy(P0),
                     ss_div(expr_copy(P1p), expr_new_integer(2))));
    bool is_symsq = ds_is_zero(cond);
    expr_free(cond); expr_free(P1p); expr_free(P0);
    expr_free(p1); expr_free(p0); expr_free(p2p); expr_free(p2pp);

    Expr* body = NULL;
    if (is_symsq) {
        /* mu = Exp[-Integrate[p2/3, x]] (the depression factor; y = mu w) */
        Expr* ip2 = ds_integrate(ss_div(expr_copy(p2), expr_new_integer(3)), expr_new_symbol(xvar));
        if (!ds_has_head(ip2, SYM_Integrate)) {
            Expr* mu = eval_and_free(ds_call1("Exp",
                           ds_call2(SYM_Times, expr_new_integer(-1), ip2)));
            /* underlying second-order normal form: u'' + (P1/4) u == 0 */
            const char* ufun = intern_symbol("DSolve`ssu");
            Expr* qc = ss_div(expr_copy(P1), expr_new_integer(4));
            Expr* ueq = expr_new_function(expr_new_symbol(SYM_Equal),
                            (Expr*[]){ ds_call2(SYM_Plus, ds_make_funcapp(ufun, 2, xvar),
                                          ds_call2(SYM_Times, qc, ds_make_funcapp(ufun, 0, xvar))),
                                       expr_new_integer(0) }, 2);
            Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                             (Expr*[]){ ueq, ds_make_funcapp(ufun, 0, xvar), expr_new_symbol(xvar) }, 3);
            Expr* ubody = ss_extract_applied(eval_and_free(call), ufun);
            if (ubody) {
                Expr* u1 = ss_coeff(ubody, ds_const(1));
                Expr* u2 = ss_coeff(ubody, ds_const(2));
                if (!ds_is_zero(u1) && !ds_is_zero(u2)) {
                    /* w = C1 u1^2 + C2 u1 u2 + C3 u2^2 */
                    Expr* w = eval_and_free(ds_call2(SYM_Plus,
                                  ds_call2(SYM_Times, ds_const(1),
                                      ds_call2(SYM_Power, expr_copy(u1), expr_new_integer(2))),
                                  ds_call2(SYM_Plus,
                                      ds_call2(SYM_Times, ds_const(2),
                                          ds_call2(SYM_Times, expr_copy(u1), expr_copy(u2))),
                                      ds_call2(SYM_Times, ds_const(3),
                                          ds_call2(SYM_Power, expr_copy(u2), expr_new_integer(2))))));
                    body = eval_and_free(ds_call2(SYM_Times, expr_copy(mu), w));
                }
                expr_free(u1); expr_free(u2); expr_free(ubody);
            }
            expr_free(mu);
        } else {
            expr_free(ip2);
        }
    }
    expr_free(p2); expr_free(P1);

    if (!body) return NULL;
    Expr** out = malloc(sizeof(Expr*));
    out[0] = body;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_symsquare(Expr* res) {
    return dsolve_method_builtin(res, dsolve_symsquare_try);
}

void dsolve_symsquare_init(void) {
    symtab_add_builtin("DSolve`SymmetricSquare", builtin_dsolve_symsquare);
    symtab_get_def("DSolve`SymmetricSquare")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`SymmetricSquare",
        "DSolve`SymmetricSquare[eqn, y, x] solves a third-order linear ODE whose "
        "solution space is the symmetric square of a second-order equation "
        "u'' = q u: it recovers q, solves for u1, u2 (Airy/Bessel/...), and "
        "returns C1 u1^2 + C2 u1 u2 + C3 u2^2 (times the depression factor).");
}
