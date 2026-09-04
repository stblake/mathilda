/*
 * dsolve_eigenvalue.c — DSolve`EigenvalueProblem (Sturm-Liouville, first cut).
 *
 * A PINNED-ONLY method (never in the automatic cascade, so it cannot misfire on an
 * ordinary IVP/BVP DSolve call).  It solves the constant-coefficient eigenvalue
 * problem
 *     y'' + lambda y == 0,   with two HOMOGENEOUS boundary conditions at two
 *                            distinct points a, b (Dirichlet y==0, Neumann y'==0,
 *                            or mixed),
 * returning the discrete eigenvalue family and its eigenfunctions.  With the ansatz
 * y = A Cos[w(x-a)] + B Sin[w(x-a)]  (lambda = w^2), the boundary condition at `a`
 * kills one component (Dirichlet kills the cosine, Neumann kills the sine), leaving a
 * single term; the condition at `b` becomes Sin[w L]==0 or Cos[w L]==0 (L = b-a):
 *
 *     order(a) == order(b)  ->  Sin[w L] == 0  ->  w_n = n Pi / L
 *     order(a) != order(b)  ->  Cos[w L] == 0  ->  w_n = (2n-1) Pi / (2 L)
 *
 * and the eigenfunction is Sin[w_n(x-a)] when the condition at `a` is Dirichlet
 * (order 0), Cos[w_n(x-a)] when it is Neumann (order 1).  The whole result is then
 * VERIFIED by back-substitution under the assumption that the family index n =
 * C[1] is an integer, so a mis-derivation can only decline (return NULL / stay
 * symbolic), never ship a wrong answer.
 *
 * Output:
 *   {{ lambda -> ConditionalExpression[w_n^2, Element[C[1],Integers] && C[1] >= 1],
 *      y      -> Function[{x}, C[2] Sin/Cos[w_n (x-a)]] }}
 * with C[1] the integer family index (matching the C[1]-in-Integers convention Solve
 * already uses for periodic families) and C[2] the free amplitude, both renamed by
 * the GeneratedParameters option.
 *
 * Out of scope (future): general Sturm-Liouville (p(x) y')' + (q + lambda w) y == 0
 * with non-constant p, q, weight w; Robin/periodic BCs; the lambda == 0 (constant)
 * Neumann eigenpair; auto-cascade dispatch; DEigensystem/DEigenvalues surfaces.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

static Expr* sym(const char* s) { return expr_new_symbol(s); }
static Expr* mul(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Times, a, b)); }
static Expr* add(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Plus, a, b)); }
static Expr* powi(Expr* base, long e) {
    return expr_new_function(sym(SYM_Power), (Expr*[]){ base, expr_new_integer(e) }, 2);
}

/* Report whether the residual `e` is zero under the assumption that the family
 * index is an integer.  `e` is consumed.  The output uses the generated constant
 * C[1] as the index (collision-safe, matching Solve's C[k]-in-Integers family
 * convention), but Simplify's assumption machinery only keys on plain SYMBOLS, not
 * on a compound C[1] — so here we first substitute C[1] -> the bare (context-
 * qualified, user-collision-free) symbol DSolve`eigN and assume THAT integer.
 * A boundary residual such as Sin[C[1] Pi] or Cos[(2 C[1]-1) Pi/2] reduces to 0
 * only under this integer index; the ODE residual is 0 without it. */
static bool zero_under_integer_index(Expr* e) {
    const char* nidx = intern_symbol("DSolve`eigN");
    e = ds_subst(e, ds_const(1), expr_new_symbol(nidx));
    /* TrigReduce first: it normalizes a factored argument like Cos[Pi(n-1/2)]
     * (which Simplify alone leaves alone) into a form the integer assumption then
     * collapses to 0.  Harmless on the already-reduced Sin[n Pi] forms. */
    e = eval_and_free(ds_call1("TrigReduce", e));
    Expr* assum = expr_new_function(sym("Element"),
                    (Expr*[]){ expr_new_symbol(nidx), sym("Integers") }, 2);
    Expr* s = eval_and_free(ds_call2("Simplify", e, assum));
    bool z = ds_is_zero(s);
    expr_free(s);
    return z;
}

static Expr* dsolve_eigenvalue_build(DSolveProblem* P) {
    if (P->is_pde || P->nfun != 1 || P->neq != 1 || P->ncond != 2) return NULL;
    const char* xvar  = P->ind_names[0];
    const char* yname = P->fun_names[0];

    /* linear structure: c2 y'' + c1 y' + c0 y == g */
    Expr** coeffs = NULL; Expr* forcing = NULL; int order = 0;
    if (!dsolve_linear_coeffs(P, &coeffs, &forcing, &order)) return NULL;

    bool okform = (order == 2) && ds_is_zero(forcing)
               && ds_is_zero(coeffs[1]) && ds_free_of(coeffs[2], xvar);
    Expr* lam = NULL;
    if (okform) {
        /* eigenparameter lambda = c0 / c2 must be a single free symbol (!= x) */
        lam = ds_simplify(eval_and_free(ds_call2(SYM_Times, expr_copy(coeffs[0]),
                  powi(expr_copy(coeffs[2]), -1))));
        if (!(lam->type == EXPR_SYMBOL && lam->data.symbol.name != xvar
              && ds_free_of(coeffs[2], lam->data.symbol.name)))
            okform = false;
    }
    for (int k = 0; k <= order; k++) if (coeffs[k]) expr_free(coeffs[k]);
    free(coeffs);
    if (forcing) expr_free(forcing);
    if (!okform) { if (lam) expr_free(lam); return NULL; }

    /* two homogeneous BCs on y (fi 0), each order 0 or 1, value 0, distinct points */
    for (size_t c = 0; c < 2; c++)
        if (P->conds[c].fi != 0 || (P->conds[c].order != 0 && P->conds[c].order != 1)
            || !ds_is_zero(P->conds[c].value)) { expr_free(lam); return NULL; }
    Expr* a = P->conds[0].point;                    /* borrowed */
    int oa = P->conds[0].order, ob = P->conds[1].order;
    Expr* diff = eval_and_free(ds_call2(SYM_Subtract,
                     expr_copy(P->conds[1].point), expr_copy(a)));   /* b - a */
    bool distinct = ds_is_nonzero(diff);
    Expr* L = ds_simplify(diff);                     /* L = b - a (consumes diff) */
    if (!distinct) { expr_free(lam); expr_free(L); return NULL; }

    /* w_n = n Pi / L (same-type BCs) or (2n-1) Pi / (2 L) (mixed); n = C[1] */
    Expr* wn;
    if (oa == ob) {
        wn = mul(mul(ds_const(1), sym("Pi")), powi(expr_copy(L), -1));
    } else {
        Expr* num = add(mul(expr_new_integer(2), ds_const(1)), expr_new_integer(-1));
        wn = mul(mul(num, sym("Pi")), powi(mul(expr_new_integer(2), expr_copy(L)), -1));
    }
    wn = ds_simplify(wn);

    /* eigenfunction: Sin[w(x-a)] if Dirichlet at a (oa==0), else Cos[w(x-a)] */
    const char* kind = (oa == 0) ? "Sin" : "Cos";
    Expr* xa   = eval_and_free(ds_call2(SYM_Subtract, sym(xvar), expr_copy(a)));
    Expr* yk   = eval_and_free(ds_call1(kind, mul(expr_copy(wn), xa)));
    Expr* yn   = ds_simplify(mul(ds_const(2), yk));      /* C[2] amplitude */
    Expr* lamn = ds_simplify(powi(expr_copy(wn), 2));    /* lambda_n = w_n^2 */
    expr_free(wn);

    /* ---- verify: ODE residual (structural) and each BC residual (under n integer) ---- */
    Expr* yxx = ds_d(ds_d(expr_copy(yn), sym(xvar)), sym(xvar));
    bool ok = zero_under_integer_index(add(yxx, mul(expr_copy(lamn), expr_copy(yn))));
    for (size_t c = 0; c < 2 && ok; c++) {
        Expr* bexpr = expr_copy(yn);
        for (int d = 0; d < P->conds[c].order; d++) bexpr = ds_d(bexpr, sym(xvar));
        bexpr = ds_subst(bexpr, sym(xvar), expr_copy(P->conds[c].point));
        ok = zero_under_integer_index(bexpr);           /* value is 0 */
    }
    if (!ok) { expr_free(lam); expr_free(L); expr_free(yn); expr_free(lamn); return NULL; }

    /* ---- assemble {{ lambda -> ConditionalExpression[...], y -> Function[...] }} ---- */
    Expr* cond = expr_new_function(sym("And"), (Expr*[]){
        expr_new_function(sym("Element"),
            (Expr*[]){ ds_const(1), sym("Integers") }, 2),
        expr_new_function(sym("GreaterEqual"),
            (Expr*[]){ ds_const(1), expr_new_integer(1) }, 2) }, 2);
    Expr* lamce = expr_new_function(sym("ConditionalExpression"),
                      (Expr*[]){ lamn, cond }, 2);        /* consumes lamn, cond */
    Expr* lam_rhs = ds_rename_param(lamce, P->param_head);
    expr_free(lamce);
    Expr* lam_rule = expr_new_function(sym(SYM_Rule),
                         (Expr*[]){ expr_copy(lam), lam_rhs }, 2);
    expr_free(lam);

    Expr* yn_r = ds_rename_param(yn, P->param_head);
    expr_free(yn);
    Expr* ylhs; Expr* yrhs;
    if (P->applied) {
        ylhs = expr_new_function(sym(yname), (Expr*[]){ sym(xvar) }, 1);
        yrhs = yn_r;
    } else {
        ylhs = sym(yname);
        Expr* plist = expr_new_function(sym(SYM_List), (Expr*[]){ sym(xvar) }, 1);
        yrhs = expr_new_function(sym(SYM_Function), (Expr*[]){ plist, yn_r }, 2);
    }
    Expr* y_rule = expr_new_function(sym(SYM_Rule), (Expr*[]){ ylhs, yrhs }, 2);

    Expr* branch = expr_new_function(sym(SYM_List), (Expr*[]){ lam_rule, y_rule }, 2);
    Expr* result = expr_new_function(sym(SYM_List), (Expr*[]){ branch }, 1);
    expr_free(L);
    return result;
}

static Expr* builtin_dsolve_eigenvalue(Expr* res) {
    DSolveProblem P;
    if (!dsolve_parse(res, &P)) { dsolve_problem_free(&P); return NULL; }
    Expr* r = P.is_pde ? NULL : dsolve_eigenvalue_build(&P);
    dsolve_problem_free(&P);
    return r;
}

void dsolve_eigenvalue_init(void) {
    symtab_add_builtin("DSolve`EigenvalueProblem", builtin_dsolve_eigenvalue);
    symtab_get_def("DSolve`EigenvalueProblem")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`EigenvalueProblem",
        "DSolve`EigenvalueProblem[{y''[x] + lambda y[x] == 0, bc1, bc2}, y, x] solves the "
        "constant-coefficient Sturm-Liouville eigenvalue problem on [a,b] with two "
        "HOMOGENEOUS boundary conditions (Dirichlet y[.]==0, Neumann y'[.]==0, or mixed) at "
        "two distinct points, returning the eigenvalue family and eigenfunctions: "
        "{{lambda -> ConditionalExpression[w_n^2, Element[C[1],Integers] && C[1]>=1], "
        "y -> Function[{x}, C[2] Sin/Cos[w_n (x-a)]]}}, where w_n = n Pi/(b-a) for same-type "
        "conditions and (2n-1)Pi/(2(b-a)) for mixed, n = C[1] the integer index and C[2] the "
        "amplitude. The eigenparameter lambda is detected as the coefficient ratio c0/c2. "
        "Every eigenpair is verified by back-substitution under n integer. Pinned-only (not "
        "in the automatic cascade). Non-constant weight, Robin/periodic BCs, and the "
        "lambda==0 Neumann mode are future work.");
}
