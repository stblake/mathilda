/*
 * dsolve_reduce_order.c — DSolve`ReductionOfOrder.
 *
 * Reduces a second-order ODE that is missing the dependent variable y — i.e.
 * y''[x] == F(x, y'[x]) with F free of y — to a first-order equation in p = y':
 *     p'[x] == F(x, p),
 * solved by recursing into the scalar engine.  The second solution then comes
 * from y = Integrate[p, x] + C[2].
 *
 * (The autonomous case y'' == F(y, y') missing x, and the energy integral
 * y'' == f(y), are left for a later pass — they typically yield elliptic /
 * implicit forms.)
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include "integrate.h"      /* g_integrate_quiet */
#include <stdlib.h>

/* Extract the RHS of {{p[x] -> expr}} (applied form). */
static Expr* extract_applied(Expr* r, const char* pfun) {
    if (!head_is(r, SYM_List) || r->data.function.arg_count < 1) return NULL;
    Expr* inner = r->data.function.args[0];
    if (!head_is(inner, SYM_List)) return NULL;
    for (size_t k = 0; k < inner->data.function.arg_count; k++) {
        Expr* rule = inner->data.function.args[k];
        if (head_is(rule, SYM_Rule) && rule->data.function.arg_count == 2) {
            Expr* lhs = rule->data.function.args[0];
            if (lhs->type == EXPR_FUNCTION && lhs->data.function.head->type == EXPR_SYMBOL
                && lhs->data.function.head->data.symbol.name == pfun)
                return expr_copy(rule->data.function.args[1]);
        }
    }
    return NULL;
}

/* The general solution y = C[2] + Integrate[p, x] when the quadrature does not
 * close: each additive term of p is integrated on its own (bounded, muted), its
 * x-free factors (the C[1] the first-order solve carries, parameters) pulled out,
 * and a term that does not close is held as Inactive[Integrate] -- the reference's
 * shape for y'' + f[x] y' == 0 (C[2] + Inactive[Integrate][C[1] E^-Inactive[
 * Integrate][f, x], x]).  General solutions only (an indefinite inert integral
 * cannot be fitted at a point).  The body is exact by the fundamental theorem
 * (checked: D[y, x] - p must cancel structurally once each inert integral is
 * differentiated).  MATHILDA_NO_RO_INERT disables it.  pbody consumed. */
static Expr* reduce_order_inert(DSolveProblem* P, Expr* pbody, const char* xvar) {
    if (P->ncond > 0 || getenv("MATHILDA_NO_RO_INERT")) { expr_free(pbody); return NULL; }
    Expr* pe = head_is(pbody, SYM_Plus) ? expr_copy(pbody)
             : eval_and_free(ds_call1("Expand", expr_copy(pbody)));
    size_t nt = head_is(pe, SYM_Plus) ? pe->data.function.arg_count : 1;
    Expr** parts = malloc((nt + 1) * sizeof(Expr*));
    for (size_t i = 0; i < nt; i++) {
        Expr* t = (nt == 1 && !head_is(pe, SYM_Plus)) ? expr_copy(pe)
                                                       : expr_copy(pe->data.function.args[i]);
        /* split t = k * g with k free of x */
        Expr* k = expr_new_integer(1);
        Expr* g = t;
        if (head_is(t, SYM_Times)) {
            size_t n = t->data.function.arg_count, nk = 0, ng = 0;
            Expr** ka = malloc(n * sizeof(Expr*));
            Expr** ga = malloc(n * sizeof(Expr*));
            for (size_t j = 0; j < n; j++) {
                Expr* f = t->data.function.args[j];
                if (!ds_contains(f, intern_symbol(xvar))) ka[nk++] = expr_copy(f);
                else ga[ng++] = expr_copy(f);
            }
            expr_free(k);
            k = (nk == 0) ? expr_new_integer(1)
                          : eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), ka, nk));
            g = (ng == 0) ? expr_new_integer(1)
                          : eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), ga, ng));
            free(ka); free(ga);
            expr_free(t);
        }
        Expr* gi = ds_integrate_or_inert(g, xvar, 3);          /* consumes g */
        parts[i] = ds_call2(SYM_Times, k, gi);
    }
    parts[nt] = ds_const(2);
    expr_free(pe);
    Expr* body = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), parts, nt + 1));
    free(parts);
    Expr* chk = eval_and_free(ds_call1("Together", ds_call2(SYM_Subtract,
                    ds_d(expr_copy(body), expr_new_symbol(xvar)), pbody)));
    bool ok = ds_is_zero(chk);
    expr_free(chk);
    if (!ok) { expr_free(body); return NULL; }
    return body;
}

static Expr** reduce_order_wrap(Expr* body, size_t* nbranch) {
    if (!body) return NULL;
    Expr** out = malloc(sizeof(Expr*));
    out[0] = body;
    *nbranch = 1;
    return out;
}

Expr** dsolve_reduce_order_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 2) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar = P->ind_names[0];

    Expr* F = dsolve_solve_top_derivative(P, 2);   /* y'' == F(x, y, y') */
    if (!F) return NULL;
    /* missing y: F must be free of the order-0 funcapp y[x].  (ds_contains on the
     * bare symbol y would false-positive on the y inside y'[x] = Derivative[1][y][x],
     * so detect the funcapp by substituting it for a marker.) */
    const char* ychk = intern_symbol("DSolve`ropY");
    Expr* Fchk = ds_subst(expr_copy(F), ds_make_funcapp(yname, 0, xvar), expr_new_symbol(ychk));
    bool has_y = ds_contains(Fchk, ychk);
    expr_free(Fchk);
    if (has_y) { expr_free(F); return NULL; }

    /* replace y'[x] with p[x], forming the reduced equation p'[x] == F(x, p) */
    const char* pfun = intern_symbol("DSolve`rop");
    Expr* Fp = ds_subst(F, ds_make_funcapp(yname, 1, xvar), ds_make_funcapp(pfun, 0, xvar));
    Expr* reduced = expr_new_function(expr_new_symbol(SYM_Equal),
                        (Expr*[]){ ds_make_funcapp(pfun, 1, xvar), Fp }, 2);

    /* solve the first-order equation for p (applied form) */
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ reduced, ds_make_funcapp(pfun, 0, xvar), expr_new_symbol(xvar) }, 3);
    Expr* r = eval_and_free(call);
    Expr* pbody = extract_applied(r, pfun);
    expr_free(r);
    if (!pbody) return NULL;

    /* y = Integrate[p, x] + C[2] (p already carries C[1]).  The quadrature is an
     * internal probe: a non-elementary p must not surface Integrate::nonelem to the
     * user (§2.2.59-5880 leaked it), so the integrator is muted around it. */
    if (ds_has_inactive_integrate(pbody)) {
        /* The first-order solve already held part of p inert (y'' + f[x] y' == 0
         * gives p = C[1] E^-Inactive[Integrate][f, x]); integrating it again cannot
         * close, so do not pay the integrator for the refusal. */
        return reduce_order_wrap(reduce_order_inert(P, pbody, xvar), nbranch);
    }
    g_integrate_quiet++;
    Expr* yint = ds_integrate(expr_copy(pbody), expr_new_symbol(xvar));
    g_integrate_quiet--;
    if (ds_has_head(yint, SYM_Integrate)) {
        expr_free(yint);
        return reduce_order_wrap(reduce_order_inert(P, pbody, xvar), nbranch);
    }
    /* Guard against a wrong antiderivative (the Integrate engine returns 0 for
     * some rational integrands with symbolic parameters, which would silently
     * yield the degenerate y = const): require D[yint, x] == p.  The correct
     * antiderivative is routinely left in a form zero_test / Simplify cannot
     * reduce — e.g. the first-order solve returns p = -Tan[1/5(-C[1]-5x)] and
     * Integrate gives a multi-Log form whose derivative differs from p only by a
     * trig identity Simplify does not crack.  Fall back to PossibleZeroQ's
     * numeric sampling, exactly as the derivative-divides integrator gates its
     * own antiderivatives.  A genuinely wrong antiderivative (Integrate -> 0)
     * leaves D[yint,x]-p == -p, which PossibleZeroQ rejects, and the final
     * back-substitution verify in dsolve_run is the authoritative backstop. */
    Expr* chk = eval_and_free(ds_call2(SYM_Subtract,
                    ds_d(expr_copy(yint), expr_new_symbol(xvar)), expr_copy(pbody)));
    bool intok = ds_is_zero(chk);
    if (!intok) {
        Expr* pz = eval_and_free(ds_call1("PossibleZeroQ", expr_copy(chk)));
        intok = (pz->type == EXPR_SYMBOL && pz->data.symbol.name == SYM_True);
        expr_free(pz);
    }
    expr_free(chk); expr_free(pbody);
    if (!intok) { expr_free(yint); return NULL; }
    Expr* body = eval_and_free(ds_call2(SYM_Plus, yint, ds_const(2)));

    Expr** out = malloc(sizeof(Expr*));
    out[0] = body;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_reduce_order(Expr* res) {
    return dsolve_method_builtin(res, dsolve_reduce_order_try);
}

void dsolve_reduce_order_init(void) {
    symtab_add_builtin("DSolve`ReductionOfOrder", builtin_dsolve_reduce_order);
    symtab_get_def("DSolve`ReductionOfOrder")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`ReductionOfOrder",
        "DSolve`ReductionOfOrder[eqn, y, x] reduces a second-order ODE missing y "
        "(y'' == F(x, y')) to first order in p = y', solves it, and integrates.");
}
