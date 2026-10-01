/*
 * dsolve_exactode.c — DSolve`ExactODE.
 *
 * Solves a linear ODE of order n >= 2 whose left side is a total derivative,
 *
 *     L[y] = Sum_{k=0}^{n} a_k(x) y^(k) == g(x),     L[y] == d/dx( M[y] ),
 *
 * where M is an order-(n-1) linear operator.  Such an equation is EXACT and
 * integrates once to the first integral
 *
 *     M[y] = Sum_{j=0}^{n-1} b_j(x) y^(j) == Integrate[g, x] + C[k],
 *
 * a lower-order ODE handed back to the scalar cascade (recurse into DSolve).
 * The first-integral coefficients follow from matching
 * Sum a_k y^(k) = d/dx( Sum b_j y^(j) ):
 *
 *     b_{n-1} = a_n,   b_{k-1} = a_k - b_k'   (k = n-1 .. 1),
 *
 * and the exactness condition is the leftover a_0 == b_0'  (equivalently
 * Sum (-1)^k a_k^(k) == 0).  The integration constant is numbered one past the
 * largest C[k] the sub-solve came back with (ex_max_const), not C[n]: for a
 * single application those coincide — the reduced order-(n-1) equation
 * contributes C[1..n-1] — but they do NOT when the recursion nests (below).
 * Iterated exactness needs no special handling: the recursion re-enters ExactODE
 * on the reduced equation, so a doubly-exact 3rd-order equation reduces twice.
 *
 * Nesting and the first-integral constant.  That constant is carried as a PLAIN
 * symbol rather than a C[k] (see the DiffUnderInt note at the call site), and it
 * must be UNIQUE PER NESTING LEVEL.  With one fixed name the outer and inner
 * constants of a doubly-exact equation were the SAME symbol, so they merged in
 * the inner sub-solve — before any rename could tell them apart — and
 * x y''' + 2 y'' == A x came back as C[1] + C[2] x + A x^3/18 + C[2] Log[x]:
 * two constants for a third-order ODE, i.e. an incomplete general solution.
 * Nothing downstream can catch that, because the residual of an incomplete
 * family is still exactly zero (the corpus harness scored it PASS).  Hence the
 * per-level name table below; the nesting is bounded (an exact reduction drops
 * the order by one each time, so the depth cannot exceed the order) and a
 * deeper nesting than the table declines rather than aliasing.
 *
 * Scope: order >= 2, genuinely exact, LINEAR by coefficient matching (above) or
 * NONLINEAR by the jet peel (further down).  First-order exact M + N y' == 0 is
 * DSolve`Exact.  Integrating-factor (adjoint) exactness is future work.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../sym_intern.h"
#include "../eval.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include <stdlib.h>

/* Extract the RHS of {{y[x] -> expr}} (applied form) from a DSolve result. */
static Expr* extract_applied(Expr* r, const char* yfun) {
    if (!head_is(r, SYM_List) || r->data.function.arg_count < 1) return NULL;
    Expr* inner = r->data.function.args[0];
    if (!head_is(inner, SYM_List)) return NULL;
    for (size_t k = 0; k < inner->data.function.arg_count; k++) {
        Expr* rule = inner->data.function.args[k];
        if (head_is(rule, SYM_Rule) && rule->data.function.arg_count == 2) {
            Expr* lhs = rule->data.function.args[0];
            if (lhs->type == EXPR_FUNCTION && lhs->data.function.head->type == EXPR_SYMBOL
                && lhs->data.function.head->data.symbol.name == yfun)
                return expr_copy(rule->data.function.args[1]);
        }
    }
    return NULL;
}

/* Largest index k of any C[k] in e (0 if none) — the sub-solve's own constants. */
static void ex_scan_const(const Expr* e, const char* cname, int* mx) {
    if (!e || e->type != EXPR_FUNCTION) return;
    Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && h->data.symbol.name == cname
        && e->data.function.arg_count == 1
        && e->data.function.args[0]->type == EXPR_INTEGER) {
        int idx = (int)e->data.function.args[0]->data.integer;
        if (idx > *mx) *mx = idx;
    }
    ex_scan_const(h, cname, mx);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        ex_scan_const(e->data.function.args[i], cname, mx);
}

static int ex_max_const(const Expr* e) {
    int mx = 0;
    ex_scan_const(e, intern_symbol("C"), &mx);
    return mx;
}

/* ---- nonlinear exactness ---------------------------------------------------
 *
 * The linear path above matches coefficients; a NONLINEAR left side needs the
 * general construction, over the jet variables j_k = y^(k) treated as independent
 * symbols.  If  L = dF/dx  for some F(x, j_0, ..., j_{n-1}), then expanding the
 * total derivative
 *
 *     dF/dx = F_x + SUM_{k=0}^{n-1} F_{j_k} j_{k+1}
 *
 * shows dL/dj_n = F_{j_{n-1}}, so F's dependence on the top jet is recovered by
 * one integration, and subtracting that piece's total derivative leaves a shorter
 * expression to which the same step applies.  Peeling k = n down to 1 either
 * exhausts L — the equation is exact — or leaves a remainder that still carries a
 * derivative, which proves it is not.  Everything is then checked by one exact
 * total-derivative comparison, so a wrong F can never be emitted.
 *
 * This reaches the `_nth_order _exact _nonlinear` family.  2 y y''' + 2(y+3y')y''
 * + 2 y'^2 == Sin[x] is the §2.2.35-3497 member: the peel gives
 * F = 2 y y'' + 2 y y' + 2 y'^2 (which is (y^2)'' + (y^2)'), so the first integral
 * is 2 y y'' + 2 y y' + 2 y'^2 == -Cos[x] + C[3], handed back to the cascade like
 * any other reduction.  Strictly additive: it runs only where
 * dsolve_linear_coeffs declines, so the linear path is untouched. */
#define DS_EXJET_MAX 9
static const char* const ds_exjet[DS_EXJET_MAX] = {
    "DSolve`exj0", "DSolve`exj1", "DSolve`exj2", "DSolve`exj3", "DSolve`exj4",
    "DSolve`exj5", "DSolve`exj6", "DSolve`exj7", "DSolve`exj8"
};

/* dF/dx with F = F(x, j_0..j_{n-1}): the chain rule over the jet variables. */
static Expr* exj_total_d(const Expr* F, int n, const char* xvar) {
    Expr* acc = ds_d(expr_copy((Expr*)F), expr_new_symbol(xvar));
    for (int k = 0; k <= n - 1; k++) {
        Expr* p = ds_d(expr_copy((Expr*)F), expr_new_symbol(ds_exjet[k]));
        if (ds_is_zero(p)) { expr_free(p); continue; }
        acc = eval_and_free(ds_call2(SYM_Plus, acc,
                  ds_call2(SYM_Times, p, expr_new_symbol(ds_exjet[k + 1]))));
    }
    return eval_and_free(ds_call1("Expand", acc));
}

/* The first integral F of a nonlinear exact operator part L (order n), in the jet
 * variables, or NULL when L is not a total derivative.  `L` borrowed. */
static Expr* exj_first_integral(const Expr* L, int n, const char* yname,
                                const char* xvar) {
    if (n < 2 || n >= DS_EXJET_MAX) return NULL;
    Expr* R = expr_copy((Expr*)L);
    for (int k = n; k >= 0; k--)
        R = ds_subst(R, ds_make_funcapp(yname, k, xvar),
                     expr_new_symbol(ds_exjet[k]));
    R = eval_and_free(ds_call1("Expand", R));

    Expr* F = expr_new_integer(0);
    bool ok = true;
    for (int k = n; k >= 1 && ok; k--) {
        Expr* c = ds_d(expr_copy(R), expr_new_symbol(ds_exjet[k]));
        if (ds_is_zero(c)) { expr_free(c); continue; }
        /* F_{j_{k-1}} cannot itself depend on j_k: L must be linear in its top
         * jet, which is exactly the condition for the peel to be well posed. */
        if (ds_contains(c, intern_symbol(ds_exjet[k]))) { expr_free(c); ok = false; break; }
        Expr* Fk = ds_integrate(expr_copy(c), expr_new_symbol(ds_exjet[k - 1]));
        expr_free(c);
        if (ds_has_head(Fk, SYM_Integrate)) { expr_free(Fk); ok = false; break; }
        Expr* dFk = exj_total_d(Fk, n, xvar);
        R = eval_and_free(ds_call1("Expand", ds_call2(SYM_Subtract, R, dFk)));
        F = eval_and_free(ds_call2(SYM_Plus, F, Fk));
    }
    /* What is left may only be a function of x: a term carrying any jet cannot be
     * produced by dF/dx once every jet coefficient has been peeled. */
    for (int k = 0; k <= n && ok; k++)
        if (ds_contains(R, intern_symbol(ds_exjet[k]))) ok = false;
    if (ok && !ds_is_zero(R)) {
        Expr* G = ds_integrate(expr_copy(R), expr_new_symbol(xvar));
        if (ds_has_head(G, SYM_Integrate)) { expr_free(G); ok = false; }
        else F = eval_and_free(ds_call2(SYM_Plus, F, G));
    }
    expr_free(R);
    if (!ok) { expr_free(F); return NULL; }

    /* Certificate: dF/dx must be L exactly (in the jet variables). */
    Expr* chk = exj_total_d(F, n, xvar);
    Expr* Lj = expr_copy((Expr*)L);
    for (int k = n; k >= 0; k--)
        Lj = ds_subst(Lj, ds_make_funcapp(yname, k, xvar), expr_new_symbol(ds_exjet[k]));
    Expr* diff = eval_and_free(ds_call1("Expand", ds_call2(SYM_Subtract, chk, Lj)));
    bool exact = ds_is_zero(diff);
    if (!exact) {           /* one Simplify pass for a rational-coefficient L */
        Expr* d2 = ds_simplify(expr_copy(diff));
        exact = ds_is_zero(d2);
        expr_free(d2);
    }
    expr_free(diff);
    if (!exact) { expr_free(F); return NULL; }

    /* Back to y, y', ... — the form the cascade re-enters on. */
    for (int k = 0; k <= n - 1; k++)
        F = ds_subst(F, expr_new_symbol(ds_exjet[k]), ds_make_funcapp(yname, k, xvar));
    return F;
}

/* One first-integral-constant name per nesting level (see the header note).  The
 * depth of an exact reduction is bounded by the ODE order, so eight is ample;
 * exceeding it declines instead of reusing a name. */
#define DS_EXFIC_NEST_MAX 8
static const char* const ds_exfic_names[DS_EXFIC_NEST_MAX] = {
    "DSolve`exFIC1", "DSolve`exFIC2", "DSolve`exFIC3", "DSolve`exFIC4",
    "DSolve`exFIC5", "DSolve`exFIC6", "DSolve`exFIC7", "DSolve`exFIC8"
};
static int ds_exfic_nest = 0;

/* mode: 1 = linear coefficient-matching only, 2 = nonlinear jet peel only, 3 = both
 * (the pinned builtin).  The two paths occupy DIFFERENT cascade slots -- see the
 * wrappers at the bottom of this file for why. */
static Expr** exactode_core(DSolveProblem* P, size_t* nbranch, int mode) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar = P->ind_names[0];

    Expr** c = NULL; Expr* g = NULL; int n = 0;
    Expr** b = NULL;            /* linear path: first-integral coefficients */
    Expr* Mnl = NULL;           /* nonlinear path: the first integral itself */
    bool exact = false;

    if ((mode & 1) && dsolve_linear_coeffs(P, &c, &g, &n)) {
        if (n < 2) {                                         /* first-order is DSolve`Exact */
            for (int k = 0; k <= n; k++) expr_free(c[k]);
            free(c); expr_free(g);
            return NULL;
        }
        /* First-integral coefficients b[0..n-1]:  b[n-1]=a_n, b[k-1]=a_k - b[k]'. */
        b = malloc((size_t)n * sizeof(Expr*));
        b[n - 1] = expr_copy(c[n]);
        for (int k = n - 1; k >= 1; k--) {
            Expr* db = ds_d(expr_copy(b[k]), expr_new_symbol(xvar));
            b[k - 1] = eval_and_free(ds_call2(SYM_Subtract, expr_copy(c[k]), db));
        }
        /* Exactness: a_0 - b_0' == 0 (Together first to help zero_test on rationals). */
        Expr* db0 = ds_d(expr_copy(b[0]), expr_new_symbol(xvar));
        Expr* chk = eval_and_free(ds_call2(SYM_Subtract, expr_copy(c[0]), db0));
        Expr* chk_t = eval_and_free(ds_call1("Together", expr_copy(chk)));
        exact = ds_is_zero(chk_t);
        expr_free(chk); expr_free(chk_t);
    } else if (mode & 2) {
        /* NONLINEAR: peel the total derivative over the jet variables (see above).
         * The forcing is what survives setting every jet to 0; the operator part is
         * the rest.
         *
         * GENERAL SOLUTION ONLY.  ExactODE sits early in the cascade, well ahead of
         * the nonlinear-2nd-order specialists, so a new claimant here preempts them
         * — and for an IVP that is a regression, not a gain: the reduction's
         * constants are introduced across two levels of recursion and the fit has to
         * invert the composition, where `AutonomousReduction` fits its stage-1
         * constant directly from the point conditions and answers
         * `y''+2yy'==0, y(0)=0, y'(0)=1` with `Tanh[x]` (M61).  Measured: without
         * this gate the peel claimed §2.2.34-3345/3347 and lost both.  The linear
         * path above is unaffected — it has always run for IVPs and its constants
         * are fitted by the substrate. */
        if (P->ncond > 0) return NULL;
        n = P->max_order[0];
        if (n < 2 || n >= DS_EXJET_MAX) return NULL;
        Expr* R = expr_copy((Expr*)P->eq_residuals[0]);
        Expr* g0 = expr_copy(R);
        for (int k = n; k >= 0; k--)
            g0 = ds_subst(g0, ds_make_funcapp(yname, k, xvar), expr_new_integer(0));
        /* R = L - g with L the operator part and g the forcing, so g0 = R|jets=0
         * is -g and L = R - g0.  Getting this sign wrong still PASSES the
         * dF/dx == L certificate below, because that certificate is checked against
         * the same wrong L: `R + g0` produced a first integral off by 2 Cos[x] on
         * §2.2.35-3497 and was only noticed as garbage two recursion levels further
         * down.  A certificate has to be anchored to the ORIGINAL residual. */
        g = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(g0)));
        Expr* L = eval_and_free(ds_call2(SYM_Subtract, R, g0));
        Mnl = exj_first_integral(L, n, yname, xvar);
        expr_free(L);
        if (!Mnl) { expr_free(g); return NULL; }
        exact = true;
    } else {
        /* mode asked for the other path only, and that path declined. */
        if (c) { for (int k = 0; k <= n; k++) expr_free(c[k]); free(c); }
        expr_free(g);
        return NULL;
    }

    Expr* body = NULL;
    if (exact && ds_exfic_nest < DS_EXFIC_NEST_MAX) {
        /* Reduced forcing: Integrate[g, x] + <first-integral constant>.  Decline on
         * a non-elementary antiderivative (an unevaluated Integrate would leak into
         * the sub-solve).  The first-integral constant is a PLAIN symbol, not C[n]:
         * a generated C[n] on the RHS drives the reduced first-order engine's
         * integrating-factor quadrature into the DiffUnderInt escalation (a C[_]
         * head reads as a parametric function), which hangs on the non-elementary
         * Integrate[C[2] E^(-Cos[x]), x] of §2.2.14-1384 -- a plain symbol returns
         * the inert quadrature instantly.  Renamed to a C[k] in the body, one past
         * the largest the sub-solve used (the per-level name keeps a nested
         * reduction's two constants distinct -- see the header note). */
        Expr* Gint = ds_integrate(expr_copy(g), expr_new_symbol(xvar));
        if (!ds_has_head(Gint, SYM_Integrate)) {
            const char* fic = ds_exfic_names[ds_exfic_nest];  /* per-level plain symbol */
            Expr* rhs = eval_and_free(ds_call2(SYM_Plus, Gint, expr_new_symbol(fic)));

            /* M[y] = Sum_{j=0}^{n-1} b[j] * y^(j), or the peeled nonlinear one. */
            Expr* M;
            if (b) {
                Expr** terms = malloc((size_t)n * sizeof(Expr*));
                for (int j = 0; j < n; j++)
                    terms[j] = ds_call2(SYM_Times, expr_copy(b[j]), ds_make_funcapp(yname, j, xvar));
                M = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)n));
                free(terms);
            } else {
                M = expr_copy(Mnl);
            }

            Expr* reduced = expr_new_function(expr_new_symbol(SYM_Equal),
                                (Expr*[]){ M, rhs }, 2);

            /* Recurse into the scalar cascade on the order-(n-1) equation.  The
             * NONLINEAR path bounds that sub-solve, the standard kit every recursing
             * method here carries (ChangeOfVariable, SecondOrderSymmetry, SolvableForY,
             * NthAlgebraic): the reduction is a full DSolve whose cost is data-dependent,
             * and this path runs last, so time spent here is time the caller has already
             * spent on every cheaper method.  Where the peel is the route that closes an
             * equation it is fast (2.1 s end to end on §2.2.35-3497), so a few seconds is
             * the right bound.  The LINEAR path is left unbounded: it is the historical
             * behaviour, it runs early because its first integral is the cheapest
             * reduction available, and nothing else claims those equations. */
            Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                             (Expr*[]){ reduced, ds_make_funcapp(yname, 0, xvar),
                                        expr_new_symbol(xvar) }, 3);
            if (!b)   /* nonlinear path */
                call = expr_new_function(expr_new_symbol(SYM_TimeConstrained),
                           (Expr*[]){ call, expr_new_integer(3),
                                      expr_new_symbol(intern_symbol("$Aborted")) }, 3);
            ds_exfic_nest++;
            Expr* r = eval_and_free(call);
            ds_exfic_nest--;
            body = extract_applied(r, yname);
            /* A reduced sub-solve that leaves an unevaluated Integrate is not a
             * usable closed form: the first-order linear engine returns the inert
             * integrating-factor quadrature quickly (Integrate[exFIC E^(-Cos[x]),x]
             * for y''+Sin[x]y'+Cos[x]y==0, §2.2.14-1384; Integrate[x^(-3/2)E^(-x^2/2),x]
             * for 2x y''+(1-2x^2)y'-4x y==0).  Decline BEFORE the exFIC->C[k] rename:
             * that rename re-evaluates the body, and renaming exFIC->C[2] INSIDE the
             * inert Integrate would re-trigger the very C[2]-DiffUnderInt hang the
             * plain symbol avoided.  Declining sends the cascade to the Frobenius
             * series (which fits the ICs at the ordinary point x=0). */
            if (body && ds_has_head(body, SYM_Integrate)) { expr_free(body); body = NULL; }
            if (body)
                body = ds_subst(body, expr_new_symbol(fic),
                                ds_const(ex_max_const(body) + 1));
            expr_free(r);
        } else {
            expr_free(Gint);
        }
    }

    if (c) { for (int k = 0; k <= n; k++) expr_free(c[k]); free(c); }
    expr_free(g);
    if (b) { for (int j = 0; j < n; j++) expr_free(b[j]); free(b); }
    if (Mnl) expr_free(Mnl);

    if (!body) return NULL;
    Expr** out = malloc(sizeof(Expr*));
    out[0] = body;
    *nbranch = 1;
    return out;
}

/* The LINEAR path, in its historical early cascade slot: coefficient matching is
 * microseconds and its first integral is the cleanest available reduction. */
Expr** dsolve_exactode_try(DSolveProblem* P, size_t* nbranch) {
    return exactode_core(P, nbranch, 1);
}

/* The NONLINEAR path, in a LATE slot after the nonlinear-second-order specialists.
 * It is a general backstop, and a general backstop placed early is a liability: the
 * peel claims any equation whose left side happens to be a total derivative, then
 * recurses on a first integral the specialists would have answered directly.
 * Measured on §2.1.2-1143 (y y''' == y' y'', whose first integral is y y'' - y'^2):
 * from the early slot it reached the SAME answer in 9.96 s where
 * `AutonomousReduction` takes 0.08 s, i.e. a 125x latency regression for no new
 * answer.  From here it can only turn a decline into an answer. */
Expr** dsolve_exactode_nl_try(DSolveProblem* P, size_t* nbranch) {
    return exactode_core(P, nbranch, 2);
}

static Expr* builtin_dsolve_exact_ode(Expr* res) {
    /* Pinned: both paths, so DSolve`ExactODE[...] is complete on its own. */
    Expr* r = dsolve_method_builtin(res, dsolve_exactode_try);
    if (!r) r = dsolve_method_builtin(res, dsolve_exactode_nl_try);
    return r;
}

void dsolve_exactode_init(void) {
    symtab_add_builtin("DSolve`ExactODE", builtin_dsolve_exact_ode);
    symtab_get_def("DSolve`ExactODE")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`ExactODE",
        "DSolve`ExactODE[eqn, y, x] solves a linear ODE of order >= 2 whose left "
        "side is a total derivative d/dx(M[y]) (an exact equation): it integrates "
        "once to the first integral M[y] == Integrate[g, x] + C[k] and recurses "
        "into the scalar solver.");
}
