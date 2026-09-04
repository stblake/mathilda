/*
 * dsolve_pdequasi.c — DSolve`PDEQuasilinear (Lagrange's method of characteristics).
 *
 * Solves the first-order quasilinear PDE
 *
 *     P(v1,v2,u) u_{v1} + Q(v1,v2,u) u_{v2} == R(v1,v2,u)
 *
 * (linear in the first derivatives; the coefficients may depend on u) by the
 * characteristic system  dv1/P = dv2/Q = du/R.  Two first integrals phi1, phi2
 * give the general solution phi2 == F(phi1)  (F arbitrary).  This is the
 * canonical generalization of the constant-coefficient DSolve`PDELinearFirstOrder
 * (dsolve_pde1.c), which handles only P, Q constant.
 *
 * First cut — two tractable, bounded classes; anything else declines (the
 * verifier is the correctness gate, so a decline is the only other outcome):
 *
 *   (A) Semilinear: P, Q free of u, and R affine in u (R = c(v1,v2) u + d).
 *       The base characteristic dv2/dv1 = Q/P decouples in (v1,v2); recurse the
 *       scalar ODE cascade to get its first integral xi(v1,v2).  u along a
 *       characteristic then solves the linear ODE u' + (P u-term) = ... via the
 *       integrating factor, exactly as pde1 does but with a computed xi and, when
 *       needed, y expressed along the characteristic (Solve xi == const for v2).
 *       Output is the EXPLICIT u == body(v1,v2) with the arbitrary function
 *       C[1][xi] (routed to the ordinary explicit PDE verify/assemble).
 *
 *   (B) Conservation law (R == 0, P or Q genuinely depending on u): u is constant
 *       along characteristics, so phi2 = u; phi1 solves Q dv1 - P dv2 = 0 with u
 *       held as a parameter (recurse the scalar cascade).  Output is the IMPLICIT
 *       relation phi1(v1,v2,u) == C[1][u] (inviscid Burgers u u_x + u_y == 0 ->
 *       u == C[1][... ], traffic flow), verified by the implicit-function rule.
 *
 * Declines: P or Q depends on u AND R != 0 (needs the full characteristic system
 * — future); a characteristic integral that stays inert (non-elementary).
 * Reuses: recursion into DSolve (characteristic ODEs), ds_solve/ds_d/ds_integrate,
 * and the dsolve_run_pde_implicit substrate (PDEExplicit / PDEImplicit wrappers).
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include <stdlib.h>

static Expr* ev1(const char* h, Expr* a) { return eval_and_free(ds_call1(h, a)); }
static Expr* mul(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Times, a, b)); }
static Expr* add(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Plus, a, b)); }
static Expr* inv(Expr* a) { return expr_new_function(expr_new_symbol(SYM_Power),
                                (Expr*[]){ a, expr_new_integer(-1) }, 2); }

/* Derivative[o1,o2][u][v1,v2]. */
static Expr* pde_deriv(const char* u, int o1, int o2, const char* v1, const char* v2) {
    Expr* d = expr_new_function(expr_new_symbol(SYM_Derivative),
                  (Expr*[]){ expr_new_integer(o1), expr_new_integer(o2) }, 2);
    Expr* du = expr_new_function(d, (Expr*[]){ expr_new_symbol(u) }, 1);
    return expr_new_function(du, (Expr*[]){ expr_new_symbol(v1), expr_new_symbol(v2) }, 2);
}

/* A first integral xi(v1,v2[,u]) of the characteristic base ODE
 * a(v1,v2) y'[v1] == b(v1,v2)  (u carried as a symbolic constant when a/b depend
 * on it).  Recurse the scalar DSolve cascade, convert the single branch to a
 * relation in v1, v2 and the constant kk, and Solve for kk.  NULL on decline.
 * a, b borrowed. */
static Expr* pdequasi_invariant(const Expr* a, const Expr* b,
                                const char* v1, const char* v2) {
    const char* chy = intern_symbol("DSolve`chy");
    const char* kk  = intern_symbol("DSolve`chk");

    Expr* yapp = expr_new_function(expr_new_symbol(chy),
                     (Expr*[]){ expr_new_symbol(v1) }, 1);
    Expr* a_s = ds_subst(expr_copy((Expr*)a), expr_new_symbol(v2), expr_copy(yapp));
    Expr* b_s = ds_subst(expr_copy((Expr*)b), expr_new_symbol(v2), yapp); /* consumes yapp */

    Expr* yder = expr_new_function(expr_new_symbol(SYM_Derivative),
                     (Expr*[]){ expr_new_integer(1) }, 1);         /* Derivative[1]     */
    yder = expr_new_function(yder, (Expr*[]){ expr_new_symbol(chy) }, 1);  /* [chy]     */
    yder = expr_new_function(yder, (Expr*[]){ expr_new_symbol(v1) }, 1);   /* [v1]      */
    Expr* ode = eval_and_free(ds_call2(SYM_Subtract,
                    ds_call2(SYM_Times, a_s, yder), b_s));
    Expr* eqn = expr_new_function(expr_new_symbol(SYM_Equal),
                    (Expr*[]){ ode, expr_new_integer(0) }, 2);
    Expr* chyapp = expr_new_function(expr_new_symbol(chy),
                       (Expr*[]){ expr_new_symbol(v1) }, 1);
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ eqn, chyapp, expr_new_symbol(v1) }, 3);
    Expr* r = eval_and_free(call);

    Expr* rel = NULL;
    if (r && head_is(r, SYM_List) && r->data.function.arg_count >= 1) {
        Expr* branch = r->data.function.args[0];
        if (head_is(branch, SYM_List) && branch->data.function.arg_count == 1) {
            Expr* el = branch->data.function.args[0];
            if (head_is(el, SYM_Rule) && el->data.function.arg_count == 2) {
                rel = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){
                    expr_new_symbol(v2), expr_copy(el->data.function.args[1]) }, 2);
            } else if (head_is(el, SYM_Equal)) {
                rel = expr_copy(el);
            }
        }
    }
    if (r) expr_free(r);
    if (!rel) return NULL;

    Expr* chyv1 = expr_new_function(expr_new_symbol(chy),
                      (Expr*[]){ expr_new_symbol(v1) }, 1);
    rel = ds_subst(rel, chyv1, expr_new_symbol(v2));      /* chy[v1] -> v2 (implicit G) */
    rel = ds_subst(rel, ds_const(1), expr_new_symbol(kk));/* C[1] -> kk                 */

    Expr* sol = ds_solve(rel, expr_new_symbol(kk));
    size_t n = 0;
    Expr** vals = dsolve_extract_solutions(sol, kk, &n);
    if (sol) expr_free(sol);
    Expr* xi = NULL;
    if (vals && n >= 1 && ds_free_of(vals[0], kk)) xi = expr_copy(vals[0]);
    if (vals) { for (size_t i = 0; i < n; i++) expr_free(vals[i]); free(vals); }
    return xi;
}

/* Build the single-element bodies[] array carrying `wrap`. */
static Expr** one(Expr* wrap) {
    Expr** out = malloc(sizeof(Expr*));
    out[0] = wrap;
    return out;
}

Expr** dsolve_pdequasi_solve(DSolveProblem* P) {
    if (P->nfun != 1 || P->nind != 2 || P->neq != 1) return NULL;
    const char* uname = P->fun_names[0];
    const char* v1 = P->ind_names[0];
    const char* v2 = P->ind_names[1];
    const char* sUx = intern_symbol("DSolve`pdeUx");
    const char* sUy = intern_symbol("DSolve`pdeUy");
    const char* sU  = intern_symbol("DSolve`pdeU");

    /* algebraic residual: u_{v1} -> sUx, u_{v2} -> sUy, u -> sU */
    Expr* R = expr_copy(P->eq_residuals[0]);
    R = ds_subst(R, pde_deriv(uname, 1, 0, v1, v2), expr_new_symbol(sUx));
    R = ds_subst(R, pde_deriv(uname, 0, 1, v1, v2), expr_new_symbol(sUy));
    R = ds_subst(R, expr_new_function(expr_new_symbol(uname),
                     (Expr*[]){ expr_new_symbol(v1), expr_new_symbol(v2) }, 2),
                 expr_new_symbol(sU));

    Expr* a = ds_d(expr_copy(R), expr_new_symbol(sUx));
    Expr* b = ds_d(expr_copy(R), expr_new_symbol(sUy));
    /* quasilinear: linear in the first derivatives */
    bool linear = ds_free_of(a, sUx) && ds_free_of(a, sUy)
               && ds_free_of(b, sUx) && ds_free_of(b, sUy);
    if (!linear) { expr_free(R); expr_free(a); expr_free(b); return NULL; }

    /* Rhs = -(R with the derivative symbols zeroed) */
    Expr* rest = expr_copy(R);
    rest = ds_subst(rest, expr_new_symbol(sUx), expr_new_integer(0));
    rest = ds_subst(rest, expr_new_symbol(sUy), expr_new_integer(0));
    expr_free(R);
    Expr* Rhs = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), rest));

    /* need a != 0; else swap the two variable roles */
    if (ds_is_zero(a)) {
        const char* t = v1; v1 = v2; v2 = t;
        Expr* tt = a; a = b; b = tt;
    }
    if (ds_is_zero(a)) { expr_free(a); expr_free(b); expr_free(Rhs); return NULL; }

    bool a_has_u = !ds_free_of(a, sU);
    bool b_has_u = !ds_free_of(b, sU);

    /* -------- class (B): conservation law R == 0, coefficients touch u -------- */
    if (a_has_u || b_has_u) {
        if (!ds_is_zero(Rhs)) {          /* genuinely quasilinear, non-conservation */
            expr_free(a); expr_free(b); expr_free(Rhs); return NULL;
        }
        expr_free(Rhs);
        Expr* xi = pdequasi_invariant(a, b, v1, v2);  /* u carried as a constant */
        expr_free(a); expr_free(b);
        if (!xi) return NULL;
        Expr* wrap = expr_new_function(expr_new_symbol(intern_symbol("DSolve`PDEImplicit")),
                         (Expr*[]){ xi, expr_new_symbol(sU) }, 2);
        return one(wrap);
    }

    /* -------- class (A): semilinear, P, Q free of u, R affine in u -------- */
    Expr* c = ds_d(expr_copy(Rhs), expr_new_symbol(sU));   /* coefficient of u   */
    if (!ds_free_of(c, sU)) {                               /* R nonlinear in u   */
        expr_free(a); expr_free(b); expr_free(Rhs); expr_free(c); return NULL;
    }
    Expr* d = ds_subst(Rhs, expr_new_symbol(sU), expr_new_integer(0));  /* consumes Rhs */

    Expr* xi = pdequasi_invariant(a, b, v1, v2);
    if (!xi) { expr_free(a); expr_free(b); expr_free(c); expr_free(d); return NULL; }

    /* u-ODE along a characteristic: u' - (c/a) u == d/a  =>  Pcoef=-c/a, Qcoef=d/a */
    Expr* Pcoef = mul(mul(expr_new_integer(-1), c), inv(expr_copy(a)));   /* consumes c */
    Expr* Qcoef = mul(d, inv(expr_copy(a)));                              /* consumes d */

    /* express v2 along the characteristic when the coefficients need it */
    const char* xis = intern_symbol("DSolve`pdexi");
    if (!ds_free_of(Pcoef, v2) || !ds_free_of(Qcoef, v2)) {
        Expr* xieq = expr_new_function(expr_new_symbol(SYM_Equal),
                         (Expr*[]){ expr_copy(xi), expr_new_symbol(xis) }, 2);
        Expr* ysol = ds_solve(xieq, expr_new_symbol(v2));
        size_t ny = 0;
        Expr** ys = dsolve_extract_solutions(ysol, v2, &ny);
        if (ysol) expr_free(ysol);
        if (!ys || ny < 1) {
            if (ys) { for (size_t i = 0; i < ny; i++) expr_free(ys[i]); free(ys); }
            expr_free(a); expr_free(b); expr_free(xi); expr_free(Pcoef); expr_free(Qcoef);
            return NULL;
        }
        Pcoef = ds_subst(Pcoef, expr_new_symbol(v2), expr_copy(ys[0]));
        Qcoef = ds_subst(Qcoef, expr_new_symbol(v2), expr_copy(ys[0]));
        for (size_t i = 0; i < ny; i++) expr_free(ys[i]);
        free(ys);
    }
    expr_free(a); expr_free(b);

    /* integrating factor mu = Exp[Integrate[Pcoef, v1]] */
    Expr* Pint = ds_integrate(Pcoef, expr_new_symbol(v1));
    if (ds_has_head(Pint, SYM_Integrate)) {
        expr_free(Pint); expr_free(Qcoef); expr_free(xi); return NULL;
    }
    Expr* mu = ev1("Exp", Pint);

    /* inner = C[1][xi_sym] + Integrate[mu Qcoef, v1] */
    Expr* arb = expr_new_function(ds_const(1), (Expr*[]){ expr_new_symbol(xis) }, 1);
    Expr* inner = arb;
    if (!ds_is_zero(Qcoef)) {
        Expr* Qint = ds_integrate(mul(expr_copy(mu), Qcoef), expr_new_symbol(v1));
        if (ds_has_head(Qint, SYM_Integrate)) {
            expr_free(Qint); expr_free(inner); expr_free(mu); expr_free(xi); return NULL;
        }
        inner = add(inner, Qint);
    } else expr_free(Qcoef);

    /* body = inner / mu, then xi_sym -> xi(v1,v2) */
    Expr* body_xi = mul(inner, inv(expr_copy(mu)));
    expr_free(mu);
    Expr* body = ds_subst(body_xi, expr_new_symbol(xis), xi);   /* consumes xi */

    Expr* wrap = expr_new_function(expr_new_symbol(intern_symbol("DSolve`PDEExplicit")),
                     (Expr*[]){ body }, 1);
    return one(wrap);
}

static Expr* builtin_dsolve_pdequasi(Expr* res) {
    return dsolve_method_builtin_pde_implicit(res, dsolve_pdequasi_solve);
}

void dsolve_pdequasi_init(void) {
    symtab_add_builtin("DSolve`PDEQuasilinear", builtin_dsolve_pdequasi);
    symtab_get_def("DSolve`PDEQuasilinear")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`PDEQuasilinear",
        "DSolve`PDEQuasilinear[eqn, u, {v1, v2}] solves a first-order quasilinear "
        "PDE P u_{v1} + Q u_{v2} == R (linear in the first derivatives; P, Q, R may "
        "depend on u) by Lagrange's method of characteristics.  First cut: the "
        "semilinear class (P, Q free of u; explicit u == body with C[1][xi]) and the "
        "conservation law R == 0 (implicit phi1 == C[1][u], e.g. inviscid Burgers).  "
        "Declines coefficients depending on u with R != 0, and non-elementary "
        "characteristic integrals.");
}
