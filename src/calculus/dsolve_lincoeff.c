/*
 * dsolve_lincoeff.c — DSolve`LinearCoefficients.
 *
 * Solves y' == (a1 x + b1 y + c1) / (a2 x + b2 y + c2) — a ratio of two affine
 * forms — by removing the constant terms:
 *
 *   det = a1 b2 - a2 b1 != 0 : the two lines meet at (x0, y0); the shift
 *       u = x - x0, w = y - y0 gives the HOMOGENEOUS equation
 *       w'(u) == (a1 u + b1 w)/(a2 u + b2 w), solved by recursing DSolve.
 *   det == 0 (parallel lines): the substitution v = a1 x + b1 y makes the
 *       equation autonomous (SEPARABLE), v' == a1 + b1 (v + c1)/((a2/a1) v + c2).
 *
 * Two cascade entries mirror DSolve`Homogeneous: the explicit try recurses and
 * shifts an explicit body back, the implicit try recurses and shifts the implicit
 * first integral G == C[1] back.  Mirrors SymPy's `linear_coefficients`.  (The Lie
 * `linear` heuristic also reaches this class; this deterministic method claims it
 * earlier and with a closed form where one exists.)
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include <stdlib.h>

/* Coefficient[e, var] (var a symbol); e borrowed, result owned. */
static Expr* coeff1(const Expr* e, const char* var) {
    return eval_and_free(ds_call2("Coefficient", expr_copy((Expr*)e), expr_new_symbol(var)));
}
/* e /. {x->0, Y->0} — the constant term; e borrowed. */
static Expr* const_term(const Expr* e, const char* xvar, const char* Yn) {
    Expr* t = ds_subst(expr_copy((Expr*)e), expr_new_symbol(xvar), expr_new_integer(0));
    return ds_subst(t, expr_new_symbol(Yn), expr_new_integer(0));
}

/* Detect the affine-ratio form; on true fill the six owned coefficients. */
static bool lc_setup(DSolveProblem* P, const char* Yn,
                     Expr** a1,Expr** b1,Expr** c1, Expr** a2,Expr** b2,Expr** c2) {
    const char* yname = P->fun_names[0];
    const char* xvar  = P->ind_names[0];
    Expr* F = dsolve_solve_top_derivative(P, 1);        /* y' == F */
    if (!F) return false;
    Expr* Ft  = eval_and_free(ds_call1("Together", F));
    Expr* num = eval_and_free(ds_call1("Numerator",   expr_copy(Ft)));
    Expr* den = eval_and_free(ds_call1("Denominator",  Ft));
    /* funcapp y[x] -> plain symbol Y */
    num = ds_subst(num, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Yn));
    den = ds_subst(den, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Yn));

    *a1 = coeff1(num, xvar); *b1 = coeff1(num, Yn); *c1 = const_term(num, xvar, Yn);
    *a2 = coeff1(den, xvar); *b2 = coeff1(den, Yn); *c2 = const_term(den, xvar, Yn);

    /* require both genuinely affine: a x + b Y + c reconstructs num / den */
    Expr* rn = eval_and_free(ds_call2(SYM_Plus,
                   ds_call2(SYM_Plus, ds_call2(SYM_Times, expr_copy(*a1), expr_new_symbol(xvar)),
                                      ds_call2(SYM_Times, expr_copy(*b1), expr_new_symbol(Yn))),
                   expr_copy(*c1)));
    Expr* rd = eval_and_free(ds_call2(SYM_Plus,
                   ds_call2(SYM_Plus, ds_call2(SYM_Times, expr_copy(*a2), expr_new_symbol(xvar)),
                                      ds_call2(SYM_Times, expr_copy(*b2), expr_new_symbol(Yn))),
                   expr_copy(*c2)));
    Expr* dn = eval_and_free(ds_call2(SYM_Subtract, rn, num));   /* consumes num */
    Expr* dd = eval_and_free(ds_call2(SYM_Subtract, rd, den));   /* consumes den */
    bool ok = ds_is_zero(dn) && ds_is_zero(dd);
    expr_free(dn); expr_free(dd);
    /* the equation must genuinely couple x and y in the affine forms, and at least
     * one constant term must be present (else it is already homogeneous/separable
     * and owned by an earlier method) */
    if (ok) {
        bool has_const = !ds_is_zero(*c1) || !ds_is_zero(*c2);
        bool has_y = !ds_is_zero(*b1) || !ds_is_zero(*b2);
        if (!has_const || !has_y) ok = false;
    }
    if (!ok) { expr_free(*a1);expr_free(*b1);expr_free(*c1);
               expr_free(*a2);expr_free(*b2);expr_free(*c2); }
    return ok;
}

/* Build the reduced sub-equation and, for det!=0, the shift (x0,y0).
 * Returns the sub-equation DSolve[...] result R (owned) via recursion into the
 * scalar engine, or NULL.  On det!=0 x0o and y0o receive the shift (owned) and the
 * sub-solution is in the fresh function w(u); on det==0 x0o and y0o are NULL and the
 * sub-solution is in v(x) with b1o the b1 needed to recover y=(v-a1 x)/b1. */
static Expr* lc_recurse(DSolveProblem* P, const char* Yn, const char* un, const char* wn,
                        bool* parallel, Expr** x0o, Expr** y0o, Expr** a1o, Expr** b1o) {
    Expr *a1,*b1,*c1,*a2,*b2,*c2;
    if (!lc_setup(P, Yn, &a1,&b1,&c1,&a2,&b2,&c2)) return NULL;
    const char* xvar = P->ind_names[0];

    Expr* det = eval_and_free(ds_call2(SYM_Subtract,
                    ds_call2(SYM_Times, expr_copy(a1), expr_copy(b2)),
                    ds_call2(SYM_Times, expr_copy(a2), expr_copy(b1))));
    bool par = ds_is_zero(det);
    expr_free(det);
    *parallel = par;
    Expr* R = NULL;

    if (!par) {
        /* intersection (x0, y0) of a1 x+b1 Y+c1==0, a2 x+b2 Y+c2==0 */
        Expr* eq1 = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){
            eval_and_free(ds_call2(SYM_Plus,
                ds_call2(SYM_Plus, ds_call2(SYM_Times, expr_copy(a1), expr_new_symbol(xvar)),
                                   ds_call2(SYM_Times, expr_copy(b1), expr_new_symbol(Yn))),
                expr_copy(c1))), expr_new_integer(0) }, 2);
        Expr* eq2 = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){
            eval_and_free(ds_call2(SYM_Plus,
                ds_call2(SYM_Plus, ds_call2(SYM_Times, expr_copy(a2), expr_new_symbol(xvar)),
                                   ds_call2(SYM_Times, expr_copy(b2), expr_new_symbol(Yn))),
                expr_copy(c2))), expr_new_integer(0) }, 2);
        Expr* sys = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ eq1, eq2 }, 2);
        Expr* vars = expr_new_function(expr_new_symbol(SYM_List),
                         (Expr*[]){ expr_new_symbol(xvar), expr_new_symbol(Yn) }, 2);
        Expr* sol = eval_and_free(ds_call2("Solve", sys, vars));
        Expr* x0 = NULL; Expr* y0 = NULL;
        if (head_is(sol, SYM_List) && sol->data.function.arg_count >= 1) {
            Expr* br = sol->data.function.args[0];
            if (head_is(br, SYM_List))
                for (size_t k = 0; k < br->data.function.arg_count; k++) {
                    Expr* rl = br->data.function.args[k];
                    if (head_is(rl, SYM_Rule) && rl->data.function.arg_count == 2
                        && rl->data.function.args[0]->type == EXPR_SYMBOL) {
                        const char* nm = rl->data.function.args[0]->data.symbol.name;
                        if (nm == xvar) x0 = expr_copy(rl->data.function.args[1]);
                        else if (nm == Yn) y0 = expr_copy(rl->data.function.args[1]);
                    }
                }
        }
        expr_free(sol);
        if (x0 && y0) {
            /* homogeneous w'(u) == (a1 u + b1 w)/(a2 u + b2 w) */
            Expr* wu = ds_make_funcapp(wn, 0, un);
            Expr* nn = ds_call2(SYM_Plus, ds_call2(SYM_Times, expr_copy(a1), expr_new_symbol(un)),
                                          ds_call2(SYM_Times, expr_copy(b1), expr_copy(wu)));
            Expr* dd = ds_call2(SYM_Plus, ds_call2(SYM_Times, expr_copy(a2), expr_new_symbol(un)),
                                          ds_call2(SYM_Times, expr_copy(b2), wu));
            Expr* rhs = eval_and_free(ds_call2(SYM_Times, nn,
                            expr_new_function(expr_new_symbol(SYM_Power),
                                (Expr*[]){ dd, expr_new_integer(-1) }, 2)));
            Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                           (Expr*[]){ ds_make_funcapp(wn, 1, un), rhs }, 2);
            Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                             (Expr*[]){ eq, ds_make_funcapp(wn, 0, un), expr_new_symbol(un) }, 3);
            R = eval_and_free(call);
            *x0o = x0; *y0o = y0;
        } else { if (x0) expr_free(x0); if (y0) expr_free(y0); }
        *a1o = NULL; *b1o = NULL;
    } else {
        /* parallel: v = a1 x + b1 y, v' == a1 + b1 (v+c1)/((a2/a1) v + c2)
         * (requires a1 != 0; else decline) */
        if (!ds_is_zero(a1)) {
            Expr* vf = ds_make_funcapp("DSolve`lcv", 0, xvar);  /* v[x] */
            Expr* ratio = ds_call2(SYM_Times, expr_copy(a2),
                              expr_new_function(expr_new_symbol(SYM_Power),
                                  (Expr*[]){ expr_copy(a1), expr_new_integer(-1) }, 2));  /* a2/a1 */
            Expr* denv = ds_call2(SYM_Plus, ds_call2(SYM_Times, ratio, expr_copy(vf)), expr_copy(c2));
            Expr* numv = ds_call2(SYM_Plus, expr_copy(vf), expr_copy(c1));
            Expr* frac = eval_and_free(ds_call2(SYM_Times, numv,
                             expr_new_function(expr_new_symbol(SYM_Power),
                                 (Expr*[]){ denv, expr_new_integer(-1) }, 2)));
            Expr* rhs = eval_and_free(ds_call2(SYM_Plus, expr_copy(a1),
                            ds_call2(SYM_Times, expr_copy(b1), frac)));
            Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                           (Expr*[]){ ds_make_funcapp("DSolve`lcv", 1, xvar), rhs }, 2);
            Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                             (Expr*[]){ eq, ds_make_funcapp("DSolve`lcv", 0, xvar),
                                        expr_new_symbol(xvar) }, 3);
            R = eval_and_free(call);
            expr_free(vf);
            *a1o = expr_copy(a1); *b1o = expr_copy(b1);
        }
        *x0o = NULL; *y0o = NULL;
    }

    expr_free(a1);expr_free(b1);expr_free(c1);expr_free(a2);expr_free(b2);expr_free(c2);
    return R;
}

/* Explicit branches. */
Expr** dsolve_lincoeff_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1 || P->max_order[0] != 1) return NULL;
    const char* xvar = P->ind_names[0];
    const char* Yn = intern_symbol("DSolve`lcY");
    const char* un = intern_symbol("DSolve`lcU");
    const char* wn = intern_symbol("DSolve`lcW");
    bool par; Expr* x0=NULL; Expr* y0=NULL; Expr* a1=NULL; Expr* b1=NULL;
    Expr* R = lc_recurse(P, Yn, un, wn, &par, &x0, &y0, &a1, &b1);
    if (!R) return NULL;

    Expr** acc = NULL; size_t nacc = 0;
    if (!par) {
        /* explicit w(u) bodies -> y = (w /. u->x-x0) + y0 */
        size_t nb = 0;
        Expr** ws = dsolve_extract_applied_bodies(R, wn, &nb);
        if (ws) {
            Expr* xm = eval_and_free(ds_call2(SYM_Subtract, expr_new_symbol(xvar), expr_copy(x0)));
            acc = malloc(nb * sizeof(Expr*));
            for (size_t i = 0; i < nb; i++) {
                Expr* body = ds_subst(ws[i], expr_new_symbol(un), expr_copy(xm));
                body = eval_and_free(ds_call2(SYM_Plus, body, expr_copy(y0)));
                acc[nacc++] = body;
            }
            expr_free(xm); free(ws);
        }
    } else if (a1 && b1) {
        /* explicit v(x) -> y = (v - a1 x)/b1 */
        size_t nb = 0;
        Expr** vs = dsolve_extract_applied_bodies(R, "DSolve`lcv", &nb);
        if (vs) {
            acc = malloc(nb * sizeof(Expr*));
            for (size_t i = 0; i < nb; i++) {
                Expr* yb = eval_and_free(ds_call2(SYM_Times,
                    ds_call2(SYM_Subtract, vs[i],
                             ds_call2(SYM_Times, expr_copy(a1), expr_new_symbol(xvar))),
                    expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ expr_copy(b1), expr_new_integer(-1) }, 2)));
                acc[nacc++] = yb;
            }
            free(vs);
        }
    }
    expr_free(R);
    if (x0) expr_free(x0);
    if (y0) expr_free(y0);
    if (a1) expr_free(a1);
    if (b1) expr_free(b1);
    if (nacc == 0) { if (acc) free(acc); return NULL; }
    *nbranch = nacc;
    return acc;
}

/* Extract the first-integral LHS G from a sub-result {{ G == C[1] }}; owned/NULL. */
static Expr* first_integral_lhs(Expr* R) {
    if (!head_is(R, SYM_List) || R->data.function.arg_count < 1) return NULL;
    Expr* br = R->data.function.args[0];
    if (!head_is(br, SYM_List) || br->data.function.arg_count < 1) return NULL;
    Expr* eq = br->data.function.args[0];
    if (head_is(eq, SYM_Equal) && eq->data.function.arg_count == 2)
        return expr_copy(eq->data.function.args[0]);
    return NULL;
}

/* Implicit first integral.  det != 0 (log-spiral subset): the homogeneous
 * sub-solve returns {{ G(u, w[u]) == C[1] }}; shift w[u]->y[x]-y0, u->x-x0.
 * det == 0 (parallel): the separable sub-solve returns {{ G(v[x]) == C[1] }};
 * substitute v[x] -> a1 x + b1 y[x]. */
Expr** dsolve_lincoeff_implicit_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1 || P->max_order[0] != 1) return NULL;
    const char* xvar  = P->ind_names[0];
    const char* yname = P->fun_names[0];
    const char* Yn = intern_symbol("DSolve`lcY");
    const char* un = intern_symbol("DSolve`lcU");
    const char* wn = intern_symbol("DSolve`lcW");
    bool par; Expr* x0=NULL; Expr* y0=NULL; Expr* a1=NULL; Expr* b1=NULL;
    Expr* R = lc_recurse(P, Yn, un, wn, &par, &x0, &y0, &a1, &b1);
    if (!R) {
        if (x0) expr_free(x0);
        if (y0) expr_free(y0);
        if (a1) expr_free(a1);
        if (b1) expr_free(b1);
        return NULL;
    }
    Expr* G = first_integral_lhs(R);
    expr_free(R);

    if (G && !par && x0 && y0) {
        /* w[u] -> y[x] - y0 ;  then u -> x - x0 */
        Expr* wsub = eval_and_free(ds_call2(SYM_Subtract,
                         ds_make_funcapp(yname, 0, xvar), expr_copy(y0)));
        G = ds_subst(G, ds_make_funcapp(wn, 0, un), wsub);
        Expr* xm = eval_and_free(ds_call2(SYM_Subtract,
                       expr_new_symbol(xvar), expr_copy(x0)));
        G = ds_subst(G, expr_new_symbol(un), xm);
    } else if (G && par && a1 && b1) {
        /* v[x] -> a1 x + b1 y[x] */
        Expr* vsub = eval_and_free(ds_call2(SYM_Plus,
                         ds_call2(SYM_Times, expr_copy(a1), expr_new_symbol(xvar)),
                         ds_call2(SYM_Times, expr_copy(b1), ds_make_funcapp(yname, 0, xvar))));
        G = ds_subst(G, ds_make_funcapp("DSolve`lcv", 0, xvar), vsub);
    } else {
        if (G) expr_free(G);
        G = NULL;
    }
    if (x0) expr_free(x0);
    if (y0) expr_free(y0);
    if (a1) expr_free(a1);
    if (b1) expr_free(b1);
    if (!G) return NULL;

    Expr** out = malloc(sizeof(Expr*));
    out[0] = G;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_lincoeff(Expr* res) {
    Expr* r = dsolve_method_builtin(res, dsolve_lincoeff_try);
    if (!r) r = dsolve_method_builtin_implicit(res, dsolve_lincoeff_implicit_try);
    return r;
}

void dsolve_lincoeff_init(void) {
    symtab_add_builtin("DSolve`LinearCoefficients", builtin_dsolve_lincoeff);
    symtab_get_def("DSolve`LinearCoefficients")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`LinearCoefficients",
        "DSolve`LinearCoefficients[eqn, y, x] solves y' == (a1 x+b1 y+c1)/(a2 x+b2 y+c2): "
        "for det != 0 a shift to the lines' intersection gives a homogeneous equation, "
        "for parallel lines the substitution v = a1 x+b1 y gives a separable one.");
}
