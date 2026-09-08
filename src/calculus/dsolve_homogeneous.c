/*
 * dsolve_homogeneous.c — DSolve`Homogeneous.
 *
 * Solves the homogeneous first-order ODE  y'[x] == F(y/x)  (F of degree 0) by
 * the substitution y = v x.  Then y' = v + x v' = F, so x v' = F(v) - v, which
 * is separable:  Integrate[1/(F(v) - v), v] == Log[x] + C[1].  Solving for v and
 * multiplying by x recovers y.
 *
 * Degree-0 homogeneity is detected by substituting Y -> v x into F(x, Y) and
 * requiring the result to be free of x.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../parse.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Reject an explicit body that leaked an internal radical-rationalisation
 * placeholder (e.g. $radu0$ from radrat): it is not a real closed form, so the
 * verified implicit first integral is the correct fallback. */
static bool homog_has_internal_placeholder(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL)
        return e->data.symbol.name && strstr(e->data.symbol.name, "$rad") != NULL;
    if (e->type == EXPR_FUNCTION) {
        if (homog_has_internal_placeholder(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (homog_has_internal_placeholder(e->data.function.args[i])) return true;
    }
    return false;
}

/* Is the explicit body DEMONSTRABLY a wrong solution of the ODE?  The symbolic
 * dsolve_verify_body keeps a residual zero_test cannot decide (Solve policy), so
 * a spurious inverse-function branch — e.g. the -2 I C[1] Pi term an inverse-Log
 * Solve produces by reusing the integration constant C[1] as the branch index
 * (problem 112) — survives it.  Sample the ODE residual at a few (x, C[1], C[2])
 * points; reject ONLY when it is nonzero at the majority of points where it
 * numericizes (never on an undecidable/complex one, matching the keep-the-
 * undecidable policy).  A rejected body falls through to the implicit entry. */
static bool homog_num_wrong(const DSolveProblem* P, const Expr* body) {
    const char* yname = P->fun_names[0];
    const char* xvar  = P->ind_names[0];
    int maxord = P->max_order[0];
    Expr* R = expr_copy(P->eq_residuals[0]);
    for (int k = maxord; k >= 1; k--) {
        Expr* dk = expr_copy((Expr*)body);
        for (int i = 0; i < k; i++) dk = ds_d(dk, expr_new_symbol(xvar));
        R = ds_subst(R, ds_make_funcapp(yname, k, xvar), dk);
    }
    R = ds_subst(R, ds_make_funcapp(yname, 0, xvar), expr_copy((Expr*)body));
    static const double xs[]  = { 1.3, 1.7, 2.3, 0.7, 3.1 };
    static const double c1s[] = { 1.181, 1.4, 2.25, 0.7, 1.9 };
    static const double c2s[] = { 0.714, 0.75, 1.6, 1.3, 0.4 };
    int small = 0, big = 0;
    for (int i = 0; i < 5; i++) {
        Expr* e = expr_copy(R);
        e = ds_subst(e, ds_const(1), expr_new_real(c1s[i]));
        e = ds_subst(e, ds_const(2), expr_new_real(c2s[i]));
        e = ds_subst(e, expr_new_symbol(xvar), expr_new_real(xs[i]));
        e = eval_and_free(ds_call1("Abs", e));
        double m = (e && e->type == EXPR_REAL)    ? e->data.real
                 : (e && e->type == EXPR_INTEGER) ? (double)e->data.integer : NAN;
        expr_free(e);
        if (isnan(m) || !isfinite(m)) continue;      /* complex/singular -> skip */
        if (m < 1e-6) small++; else if (m > 1e-3) big++;
    }
    expr_free(R);
    return big > 0 && big >= small;
}

/*
 * Fallback inversion for the separable relation Integrate[1/(F(v)-v), v] ==
 * Log[x] + C[1] when Solve cannot invert the log form directly.  For a rational
 * F the antiderivative is a sum of logarithms (rational coefficients), so
 * exponentiating turns the relation into an algebraic one,
 *     Product g_i(v)^{c_i} == C[1] x,
 * whose fractional exponents are cleared by raising both sides to a small power
 * d.  Solve then returns explicit (possibly Root) branches; the per-branch
 * back-substitution verify in dsolve_run keeps the genuine ones and drops any
 * spurious root from the exponentiation.  The transcendental family (an ArcTan
 * term -> an E^ArcTan factor that no d clears, e.g. y'==(x+y)/(x-y)) leaves the
 * exponentiated form non-algebraic, so every d fails and this declines.
 */
static Expr** homog_exp_log_invert(const Expr* intV, const char* vn,
                                   const char* xvar, size_t* nv_out) {
    Expr* r1 = parse_expression("E^(a_ + b_) :> E^a E^b");
    Expr* r2 = parse_expression("E^(c_. Log[g_]) :> g^c");
    if (!r1 || !r2) { if (r1) expr_free(r1); if (r2) expr_free(r2); return NULL; }
    Expr* rules = expr_new_function(expr_new_symbol(SYM_List),
                      (Expr*[]){ r1, r2 }, 2);
    Expr* expForm = expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ expr_new_symbol(intern_symbol("E")),
                                   expr_copy((Expr*)intV) }, 2);
    Expr* form = eval_and_free(ds_call2("ReplaceRepeated", expForm, rules));

    static const int dtry[] = { 1, 2, 3, 4, 6 };
    Expr** vs = NULL; size_t nv = 0;
    for (size_t i = 0; i < sizeof(dtry) / sizeof(dtry[0]) && !vs; i++) {
        int d = dtry[i];
        /* lhs = PowerExpand[form^d] (integer exponents); rhs = (C[1] x)^d */
        Expr* lhs = eval_and_free(ds_call1("PowerExpand",
                        expr_new_function(expr_new_symbol(SYM_Power),
                            (Expr*[]){ expr_copy(form), expr_new_integer(d) }, 2)));
        Expr* cx = ds_call2(SYM_Times, ds_const(1), expr_new_symbol(xvar));
        Expr* rhs = expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ cx, expr_new_integer(d) }, 2);
        Expr* eqn = expr_new_function(expr_new_symbol(SYM_Equal),
                        (Expr*[]){ lhs, rhs }, 2);
        Expr* sol = ds_solve(eqn, expr_new_symbol(vn));
        size_t n = 0;
        Expr** got = dsolve_extract_solutions(sol, vn, &n);
        if (sol) expr_free(sol);
        if (got && n > 0) { vs = got; nv = n; }
        else if (got) free(got);
    }
    expr_free(form);
    *nv_out = nv;
    return vs;
}

Expr** dsolve_homogeneous_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar = P->ind_names[0];

    Expr* F = dsolve_solve_top_derivative(P, 1);
    if (!F) return NULL;
    const char* Yn = intern_symbol("DSolve`Y");
    const char* vn = intern_symbol("DSolve`v");
    Expr* FY = ds_subst(F, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Yn));

    /* Detection: Fvx = FY /. Y -> v x must be free of x for degree-0 homogeneity. */
    Expr* vx = ds_call2(SYM_Times, expr_new_symbol(vn), expr_new_symbol(xvar));
    Expr* Fvx = ds_subst(expr_copy(FY), expr_new_symbol(Yn), vx);
    int homog = ds_free_of(Fvx, xvar);
    expr_free(Fvx);
    if (!homog) { expr_free(FY); return NULL; }

    /* Computation: reduced RHS via Fv1 = FY /. {Y -> v, x -> 1}.  Degree-0
     * homogeneity makes F(x, v x) == F(1, v) exactly, but substituting x -> 1
     * textually eliminates x (radicals/exponentials collapse), whereas Y -> v x
     * leaves an uncancelled x for a radical RHS (Sqrt[x^2 + v^2 x^2] does not
     * reduce without x > 0) that leaks a spurious power of x into the v-integral
     * and yields a WRONG explicit solution.  denom = F(v) - v. */
    Expr* Fv1 = ds_subst(expr_copy(FY), expr_new_symbol(Yn), expr_new_symbol(vn));
    Fv1 = ds_subst(Fv1, expr_new_symbol(xvar), expr_new_integer(1));
    expr_free(FY);
    Expr* denom = eval_and_free(ds_call2(SYM_Subtract, Fv1, expr_new_symbol(vn)));
    if (ds_is_zero(denom)) { expr_free(denom); return NULL; }         /* y = C x, handled elsewhere */

    /* Integrate[1/(F(v)-v), v] == Log[x] + C[1] */
    Expr* intV = ds_integrate(
        expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ denom, expr_new_integer(-1) }, 2),
        expr_new_symbol(vn));
    if (ds_has_head(intV, SYM_Integrate)) { expr_free(intV); return NULL; }

    /* Try the direct log-form inversion first; keep a copy of the antiderivative
     * for the algebraic (exponentiated) fallback below. */
    Expr* intVcopy = expr_copy(intV);
    Expr* rhs = eval_and_free(ds_call2(SYM_Plus,
                    ds_call1("Log", expr_new_symbol(xvar)), ds_const(1)));
    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){ intV, rhs }, 2);
    Expr* solres = ds_solve(eq, expr_new_symbol(vn));
    size_t nv = 0;
    Expr** vs = dsolve_extract_solutions(solres, vn, &nv);
    if (solres) expr_free(solres);
    bool from_direct = (vs != NULL);
    /* The exponentiated-log fallback only applies to a log-sum antiderivative
     * (rational F).  For a non-log intV (denom = E^v gives intV = -E^-v) it
     * fabricates a spurious branch, so restrict it to the case it was built for
     * and otherwise decline to the verified implicit first integral. */
    if (!vs && ds_contains(intVcopy, intern_symbol("Log")))
        vs = homog_exp_log_invert(intVcopy, vn, xvar, &nv);
    expr_free(intVcopy);
    if (!vs) return NULL;

    /* Drop the whole branch set if any body leaked an internal placeholder
     * ($rad...) -> decline so the implicit first-integral fallback runs. */
    for (size_t i = 0; i < nv; i++)
        if (homog_has_internal_placeholder(vs[i])) {
            for (size_t j = 0; j < nv; j++) expr_free(vs[j]);
            free(vs);
            return NULL;
        }

    /* y = x v */
    for (size_t i = 0; i < nv; i++)
        vs[i] = eval_and_free(ds_call2(SYM_Times, expr_new_symbol(xvar), vs[i]));

    /* Drop any DIRECT-solve body that numerically fails the ODE (a spurious
     * inverse-function branch the symbolic verify cannot reject, e.g. the
     * -2 I C[1] Pi inverse-Log branch of E^(y/x) homogeneous, problem 112); if
     * none survive, decline so the verified implicit first integral is returned.
     * The exponentiate-and-clear (homog_exp_log_invert) Root branches are NOT
     * checked here: dsolve_run already back-substitution-verifies them, and a
     * Root can look numerically nonzero at a mis-indexed sample point — checking
     * them wrongly rejected the rational log-family (x+2y)/(2x+y). */
    if (from_direct) {
        size_t keep = 0;
        for (size_t i = 0; i < nv; i++) {
            if (homog_num_wrong(P, vs[i])) expr_free(vs[i]);
            else vs[keep++] = vs[i];
        }
        if (keep == 0) { free(vs); return NULL; }
        nv = keep;
    }
    *nbranch = nv;
    return vs;
}

/*
 * Implicit-solution path.  When the substitution v = y/x gives an antiderivative
 * that cannot be inverted for y (a transcendental relation — the ArcTan
 * log-spiral family, e.g. y'==(x+y)/(x-y)), return the general integral
 *     G = (Integrate[1/(F(v)-v), v] /. v -> y[x]/x) - Log[x],
 * meaning the solution G == C[1].  dsolve_run_implicit verifies it by implicit
 * differentiation and assembles {{ G == C[1] }}.
 */
Expr** dsolve_homogeneous_implicit_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar = P->ind_names[0];

    Expr* F = dsolve_solve_top_derivative(P, 1);
    if (!F) return NULL;
    const char* Yn = intern_symbol("DSolve`Y");
    const char* vn = intern_symbol("DSolve`v");
    Expr* FY = ds_subst(F, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Yn));
    Expr* vx = ds_call2(SYM_Times, expr_new_symbol(vn), expr_new_symbol(xvar));
    Expr* Fvx = ds_subst(expr_copy(FY), expr_new_symbol(Yn), vx);
    int homog = ds_free_of(Fvx, xvar);
    expr_free(Fvx);
    if (!homog) { expr_free(FY); return NULL; }
    /* Reduced RHS via Fv1 = FY /. {Y -> v, x -> 1} (see dsolve_homogeneous_try). */
    Expr* Fv1 = ds_subst(expr_copy(FY), expr_new_symbol(Yn), expr_new_symbol(vn));
    Fv1 = ds_subst(Fv1, expr_new_symbol(xvar), expr_new_integer(1));
    expr_free(FY);
    Expr* denom = eval_and_free(ds_call2(SYM_Subtract, Fv1, expr_new_symbol(vn)));
    if (ds_is_zero(denom)) { expr_free(denom); return NULL; }
    Expr* intV = ds_integrate(
        expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ denom, expr_new_integer(-1) }, 2),
        expr_new_symbol(vn));
    if (ds_has_head(intV, SYM_Integrate)) { expr_free(intV); return NULL; }

    /* G = (intV /. v -> y[x]/x) - Log[x] */
    Expr* yx = ds_call2(SYM_Times, ds_make_funcapp(yname, 0, xvar),
                   expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_new_symbol(xvar), expr_new_integer(-1) }, 2));
    Expr* G = ds_subst(intV, expr_new_symbol(vn), yx);
    G = eval_and_free(ds_call2(SYM_Subtract, G, ds_call1("Log", expr_new_symbol(xvar))));

    Expr** out = malloc(sizeof(Expr*));
    out[0] = G;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_homogeneous(Expr* res) {
    /* explicit inversion first; fall back to the implicit first-integral form */
    Expr* r = dsolve_method_builtin(res, dsolve_homogeneous_try);
    if (!r) r = dsolve_method_builtin_implicit(res, dsolve_homogeneous_implicit_try);
    return r;
}

void dsolve_homogeneous_init(void) {
    symtab_add_builtin("DSolve`Homogeneous", builtin_dsolve_homogeneous);
    symtab_get_def("DSolve`Homogeneous")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`Homogeneous",
        "DSolve`Homogeneous[eqn, y, x] solves y'[x] == F(y/x) via the substitution "
        "y = v x, which separates variables in v and x.");
}
