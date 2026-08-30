/*
 * dsolve_specialform.c — DSolve`SpecialFunctionForm.
 *
 * Recognises a homogeneous second-order linear ODE whose solutions are named
 * special functions, by matching the normalised form  y'' + P(x) y' + Q(x) y = 0
 * against a table:
 *
 *   Airy            P = 0,      Q = -(A x + B)  ->  AiryAi[u], AiryBi[u],
 *                                                   u = A^(1/3)(x + B/A)
 *   Bessel          P = 1/x,    Q = 1 - v^2/x^2 ->  BesselJ[v, x], BesselY[v, x]
 *   modified Bessel P = 1/x,    Q = -1 - v^2/x^2 -> BesselI[v, x], BesselK[v, x]
 *
 * These heads exist in Mathilda, so the substrate still back-substitution
 * verifies the result.  Equations whose solutions are functions Mathilda does
 * not have (Mathieu, Kelvin, Weierstrass, ...) are left for a later pass.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

/* C[1] b0 + C[2] b1 ; b0,b1 consumed */
static Expr* combo(Expr* b0, Expr* b1) {
    return eval_and_free(ds_call2(SYM_Plus,
        ds_call2(SYM_Times, ds_const(1), b0),
        ds_call2(SYM_Times, ds_const(2), b1)));
}
/* a^(p/q) ; a borrowed */
static Expr* powrat(const Expr* a, int p, int q) {
    Expr* rat = eval_and_free(ds_call2(SYM_Times, expr_new_integer(p),
                    expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ expr_new_integer(q), expr_new_integer(-1) }, 2)));  /* p/q */
    return eval_and_free(ds_call2(SYM_Power, expr_copy((Expr*)a), rat));
}

Expr** dsolve_specialform_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 2) return NULL;
    const char* xvar = P->ind_names[0];

    Expr** c; Expr* g; int n;
    if (!dsolve_linear_coeffs(P, &c, &g, &n)) return NULL;
    bool homog = ds_is_zero(g);
    expr_free(g);
    if (n != 2 || !homog) { for (int k = 0; k <= n; k++) expr_free(c[k]); free(c); return NULL; }

    /* normalised P = c1/c2, Q = c0/c2 */
    Expr* Pc = ds_simplify(ds_call2(SYM_Times, expr_copy(c[1]),
                   expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_copy(c[2]), expr_new_integer(-1) }, 2)));
    Expr* Qc = ds_simplify(ds_call2(SYM_Times, expr_copy(c[0]),
                   expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_copy(c[2]), expr_new_integer(-1) }, 2)));
    for (int k = 0; k <= 2; k++) expr_free(c[k]);
    free(c);

    Expr* general = NULL;

    /* ---- Airy: P == 0, Q = -(A x + B), A = -dQ/dx constant, B = -Q(0) ---- */
    if (!general && ds_is_zero(Pc)) {
        Expr* dQ = ds_d(expr_copy(Qc), expr_new_symbol(xvar));    /* Q' = -A */
        if (ds_free_of(dQ, xvar) && !ds_is_zero(dQ)) {
            Expr* Q0 = ds_subst(expr_copy(Qc), expr_new_symbol(xvar), expr_new_integer(0));
            /* require Q exactly linear: Q == dQ*x + Q0 */
            Expr* lin = eval_and_free(ds_call2(SYM_Subtract, expr_copy(Qc),
                            ds_call2(SYM_Plus, ds_call2(SYM_Times, expr_copy(dQ), expr_new_symbol(xvar)), expr_copy(Q0))));
            if (ds_is_zero(lin)) {
                Expr* A = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(dQ)));   /* A = -Q' */
                Expr* B = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(Q0)));   /* B = -Q0 */
                /* u = A^(1/3) (x + B/A) */
                Expr* cbrtA = powrat(A, 1, 3);
                Expr* BoverA = eval_and_free(ds_call2(SYM_Times, expr_copy(B),
                                   expr_new_function(expr_new_symbol(SYM_Power),
                                       (Expr*[]){ expr_copy(A), expr_new_integer(-1) }, 2)));
                Expr* u = eval_and_free(ds_call2(SYM_Times, cbrtA,
                              ds_call2(SYM_Plus, expr_new_symbol(xvar), BoverA)));
                general = combo(ds_call1("AiryAi", expr_copy(u)), ds_call1("AiryBi", expr_copy(u)));
                expr_free(u); expr_free(A); expr_free(B);
            }
            expr_free(lin); expr_free(Q0);
        }
        expr_free(dQ);
    }

    /* ---- Bessel / modified Bessel: P == 1/x, Q = s - v^2/x^2 ---- */
    if (!general) {
        Expr* oneOverX = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ expr_new_symbol(xvar), expr_new_integer(-1) }, 2));
        Expr* Pdiff = eval_and_free(ds_call2(SYM_Subtract, expr_copy(Pc), oneOverX));
        if (ds_is_zero(Pdiff)) {
            for (int s = 1; s >= -1 && !general; s -= 2) {
                /* nu^2 = x^2 (s - Q) must be free of x */
                Expr* nu2 = ds_simplify(ds_call2(SYM_Times,
                                expr_new_function(expr_new_symbol(SYM_Power),
                                    (Expr*[]){ expr_new_symbol(xvar), expr_new_integer(2) }, 2),
                                ds_call2(SYM_Subtract, expr_new_integer(s), expr_copy(Qc))));
                if (ds_free_of(nu2, xvar)) {
                    Expr* nu = powrat(nu2, 1, 2);   /* Sqrt[nu^2] */
                    const char* fJ = (s == 1) ? "BesselJ" : "BesselI";
                    const char* fY = (s == 1) ? "BesselY" : "BesselK";
                    Expr* b0 = expr_new_function(expr_new_symbol(fJ),
                                   (Expr*[]){ expr_copy(nu), expr_new_symbol(xvar) }, 2);
                    Expr* b1 = expr_new_function(expr_new_symbol(fY),
                                   (Expr*[]){ expr_copy(nu), expr_new_symbol(xvar) }, 2);
                    general = combo(b0, b1);
                    expr_free(nu);
                }
                expr_free(nu2);
            }
        }
        expr_free(Pdiff);
    }

    expr_free(Pc); expr_free(Qc);
    if (!general) return NULL;
    Expr** out = malloc(sizeof(Expr*));
    out[0] = general;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_specialform(Expr* res) {
    return dsolve_method_builtin(res, dsolve_specialform_try);
}

void dsolve_specialform_init(void) {
    symtab_add_builtin("DSolve`SpecialFunctionForm", builtin_dsolve_specialform);
    symtab_get_def("DSolve`SpecialFunctionForm")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("DSolve`SpecialFunctionForm",
        "DSolve`SpecialFunctionForm[eqn, y, x] recognises second-order linear ODEs "
        "whose solutions are Airy (y'' == (A x + B) y) or Bessel / modified Bessel "
        "functions (x^2 y'' + x y' +- (x^2 -+ v^2) y == 0).");
}
