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
 *   Kummer (1F1)    P = b/x-1,  Q = -a/x        ->  Hypergeometric1F1[a,b,x],
 *                                                   x^(1-b) 1F1[a-b+1, 2-b, x]
 *   Gauss (2F1)     W = x(1-x), Q = -a b/W,
 *                   P = (c-(a+b+1)x)/W          ->  Hypergeometric2F1[a,b,c,x],
 *                                                   x^(1-c) 2F1[a-c+1,b-c+1,2-c,x]
 *
 * These heads exist in Mathilda, so the substrate still back-substitution
 * verifies the result.  (For the hypergeometric families the residual is a
 * contiguous-relation identity that Simplify will NOT discharge symbolically;
 * the substrate keeps a branch unless zero_test PROVES it nonzero, so the
 * genuinely-zero residual survives -- the same "structurally exact" acceptance
 * Airy/Bessel rely on.  Do not "strengthen" verify to require a positive proof,
 * or these branches will be dropped.)  The Kummer/Gauss second solution carries
 * the factor x^(1-b) / x^(1-c); it is emitted only when that exponent parameter
 * (b resp. c) is a NUMBER and not an integer.  Two reasons: an integer makes the
 * two solutions dependent (or the pFq lower parameter singular), and a SYMBOLIC
 * exponent makes the verify residual a symbolic-power + pFq sum on which
 * zero_test currently hangs.  Both cases decline to the Frobenius series
 * fallback rather than emit.  The other parameters (a for Kummer; a, b for
 * Gauss) may stay symbolic.  Equations whose solutions are functions Mathilda
 * does not have (Mathieu, Kelvin, Weierstrass, LegendreQ, ...) are a later pass.
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

/* True iff e is an explicit number (NumberQ).  This guards the hypergeometric
 * second solution x^(1-b) / x^(1-c): with a SYMBOLIC exponent the back-
 * substitution residual is a sum of symbolic powers times pFq's, on which
 * zero_test / PossibleZeroQ presently hangs (a pre-existing limitation), so
 * verification would never terminate.  A numeric exponent keeps the power
 * concrete and verification bounded; a symbolic-exponent equation instead
 * declines to the Frobenius series fallback (correct, just not closed form). */
static bool specialform_is_number(Expr* e) {
    Expr* nq = eval_and_free(ds_call1("NumberQ", expr_copy(e)));
    bool r = (nq->type == EXPR_SYMBOL && nq->data.symbol.name == SYM_True);
    expr_free(nq);
    return r;
}

/* Roots {*ra, *rb} of  t^2 - S t + Pr  (owned outputs; S, Pr borrowed).
 * Prefers the clean, exact roots read off FactorList's linear factors (a
 * linear-factor Solve is radical-free, so symbolic a+b / a b factor back to
 * a and b), and falls back to the radical roots of the whole quadratic when it
 * does not split.  Returns true iff two roots were recovered. */
static bool specialform_quad_roots(Expr* S, Expr* Pr, Expr** ra, Expr** rb) {
    const char* t = intern_symbol("DSolve`hgt");
    Expr* quad = eval_and_free(ds_call2(SYM_Plus,
                     ds_call2(SYM_Power, expr_new_symbol(t), expr_new_integer(2)),
                     ds_call2(SYM_Plus,
                         ds_call2(SYM_Times, expr_new_integer(-1),
                             ds_call2(SYM_Times, expr_copy(S), expr_new_symbol(t))),
                         expr_copy(Pr))));
    Expr* got[2]; int nr = 0;

    /* clean path: exact roots from the linear factors of the quadratic */
    Expr* fl = eval_and_free(ds_call1("FactorList", expr_copy(quad)));
    if (fl && ds_has_head(fl, SYM_List)) {
        for (size_t i = 0; i < fl->data.function.arg_count && nr < 2; i++) {
            Expr* pair = fl->data.function.args[i];
            if (!ds_has_head(pair, SYM_List) || pair->data.function.arg_count != 2) continue;
            Expr* fac = pair->data.function.args[0];
            if (ds_free_of(fac, t)) continue;                       /* constant factor */
            Expr* dfac = ds_d(expr_copy(fac), expr_new_symbol(t));
            bool linear = ds_free_of(dfac, t) && !ds_is_zero(dfac);
            expr_free(dfac);
            if (!linear) { nr = 0; break; }                         /* irreducible -> fallback */
            Expr* me = pair->data.function.args[1];
            long m = (me->type == EXPR_INTEGER) ? (long)me->data.integer : 1;
            Expr* sol = ds_solve(ds_call2(SYM_Equal, expr_copy(fac), expr_new_integer(0)),
                                 expr_new_symbol(t));
            size_t k = 0;
            Expr** rs = dsolve_extract_solutions(sol, t, &k);
            if (sol) expr_free(sol);
            if (rs && k >= 1) for (long j = 0; j < m && nr < 2; j++) got[nr++] = expr_copy(rs[0]);
            if (rs) { for (size_t j = 0; j < k; j++) expr_free(rs[j]); free(rs); }
        }
    }
    if (fl) expr_free(fl);

    /* radical fallback: roots of the whole quadratic */
    if (nr < 2) {
        for (int j = 0; j < nr; j++) expr_free(got[j]);
        nr = 0;
        Expr* sol = ds_solve(ds_call2(SYM_Equal, expr_copy(quad), expr_new_integer(0)),
                             expr_new_symbol(t));
        size_t k = 0;
        Expr** rs = dsolve_extract_solutions(sol, t, &k);
        if (sol) expr_free(sol);
        if (rs) {
            for (size_t j = 0; j < k && nr < 2; j++) got[nr++] = expr_copy(rs[j]);
            for (size_t j = 0; j < k; j++) expr_free(rs[j]);
            free(rs);
        }
    }
    expr_free(quad);

    if (nr == 2) { *ra = got[0]; *rb = got[1]; return true; }
    for (int j = 0; j < nr; j++) expr_free(got[j]);
    return false;
}

Expr** dsolve_specialform_try(DSolveProblem* P, size_t* nbranch) {
    /* normalised second-order form y'' + Pc y' + Qc y == 0 (homogeneous only) */
    Expr* Pc; Expr* Qc;
    if (!dsolve_second_order_PQ(P, &Pc, &Qc)) return NULL;
    const char* xvar = P->ind_names[0];

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

    /* ---- Kummer (confluent hypergeometric 1F1): P == b/x - 1, Q == -a/x ----
     * y = C[1] 1F1[a,b,x] + C[2] x^(1-b) 1F1[a-b+1, 2-b, x].  Read a,b directly:
     * x(P+1) free of x is b, -x Q free of x is a.  The second basis needs b not
     * an integer (else 2-b is a non-positive integer -> singular pFq lower
     * parameter); an integer b declines to the Frobenius fallback. */
    if (!general) {
        Expr* kb = ds_simplify(ds_call2(SYM_Times, expr_new_symbol(xvar),
                       ds_call2(SYM_Plus, expr_copy(Pc), expr_new_integer(1))));   /* b */
        if (ds_free_of(kb, xvar)) {
            Expr* ka = ds_simplify(ds_call2(SYM_Times, expr_new_integer(-1),
                           ds_call2(SYM_Times, expr_new_symbol(xvar), expr_copy(Qc))));  /* a */
            if (ds_free_of(ka, xvar)) {
                Expr* iq = eval_and_free(ds_call1("IntegerQ", expr_copy(kb)));
                bool bint = (iq->type == EXPR_SYMBOL && iq->data.symbol.name == SYM_True);
                expr_free(iq);
                if (specialform_is_number(kb) && !bint) {
                    Expr* b0 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Hypergeometric1F1),
                                   (Expr*[]){ expr_copy(ka), expr_copy(kb), expr_new_symbol(xvar) }, 3));
                    Expr* xpow = eval_and_free(ds_call2(SYM_Power, expr_new_symbol(xvar),
                                     ds_call2(SYM_Subtract, expr_new_integer(1), expr_copy(kb))));  /* x^(1-b) */
                    Expr* a2 = eval_and_free(ds_call2(SYM_Plus,
                                   ds_call2(SYM_Subtract, expr_copy(ka), expr_copy(kb)),
                                   expr_new_integer(1)));                          /* a-b+1 */
                    Expr* b2 = eval_and_free(ds_call2(SYM_Subtract, expr_new_integer(2), expr_copy(kb)));  /* 2-b */
                    Expr* h2 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Hypergeometric1F1),
                                   (Expr*[]){ a2, b2, expr_new_symbol(xvar) }, 3));
                    Expr* b1 = eval_and_free(ds_call2(SYM_Times, xpow, h2));
                    general = combo(b0, b1);
                }
            }
            expr_free(ka);
        }
        expr_free(kb);
    }

    /* ---- Gauss (hypergeometric 2F1): W = x(1-x), Q == -a b/W,
     *      P == (c - (a+b+1)x)/W ----
     * y = C[1] 2F1[a,b,c,x] + C[2] x^(1-c) 2F1[a-c+1, b-c+1, 2-c, x].  W Q free
     * of x gives -a b (the tight gate that rejects Bessel's -v^2/x^2 term); W P
     * linear in x gives c = value at 0 and a+b = -(slope)-1; a,b are the roots
     * of t^2 - (a+b) t + a b.  Integer c declines to Frobenius. */
    if (!general) {
        Expr* W = eval_and_free(ds_call2(SYM_Times, expr_new_symbol(xvar),
                      ds_call2(SYM_Subtract, expr_new_integer(1), expr_new_symbol(xvar))));  /* x(1-x) */
        Expr* negprod = ds_simplify(ds_call2(SYM_Times, expr_copy(W), expr_copy(Qc)));       /* -a b */
        if (ds_free_of(negprod, xvar)) {
            Expr* L = ds_simplify(ds_call2(SYM_Times, expr_copy(W), expr_copy(Pc)));          /* c-(a+b+1)x */
            Expr* dL = ds_d(expr_copy(L), expr_new_symbol(xvar));                            /* -(a+b+1) */
            if (ds_free_of(dL, xvar)) {
                Expr* gc = ds_subst(expr_copy(L), expr_new_symbol(xvar), expr_new_integer(0));  /* c = L(0) */
                Expr* iq = eval_and_free(ds_call1("IntegerQ", expr_copy(gc)));
                bool cint = (iq->type == EXPR_SYMBOL && iq->data.symbol.name == SYM_True);
                expr_free(iq);
                if (specialform_is_number(gc) && !cint) {
                    Expr* S = eval_and_free(ds_call2(SYM_Subtract,
                                  ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(dL)),
                                  expr_new_integer(1)));                            /* a+b = -dL - 1 */
                    Expr* Pr = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1),
                                   expr_copy(negprod)));                            /* a b */
                    Expr* ga = NULL; Expr* gb = NULL;
                    if (specialform_quad_roots(S, Pr, &ga, &gb)) {
                        Expr* b0 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Hypergeometric2F1),
                                       (Expr*[]){ expr_copy(ga), expr_copy(gb), expr_copy(gc),
                                                  expr_new_symbol(xvar) }, 4));
                        Expr* xpow = eval_and_free(ds_call2(SYM_Power, expr_new_symbol(xvar),
                                         ds_call2(SYM_Subtract, expr_new_integer(1), expr_copy(gc))));  /* x^(1-c) */
                        Expr* a2 = eval_and_free(ds_call2(SYM_Plus,
                                       ds_call2(SYM_Subtract, expr_copy(ga), expr_copy(gc)),
                                       expr_new_integer(1)));                       /* a-c+1 */
                        Expr* b2 = eval_and_free(ds_call2(SYM_Plus,
                                       ds_call2(SYM_Subtract, expr_copy(gb), expr_copy(gc)),
                                       expr_new_integer(1)));                       /* b-c+1 */
                        Expr* c2 = eval_and_free(ds_call2(SYM_Subtract, expr_new_integer(2), expr_copy(gc)));  /* 2-c */
                        Expr* h2 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Hypergeometric2F1),
                                       (Expr*[]){ a2, b2, c2, expr_new_symbol(xvar) }, 4));
                        Expr* b1 = eval_and_free(ds_call2(SYM_Times, xpow, h2));
                        general = combo(b0, b1);
                        expr_free(ga); expr_free(gb);
                    }
                    expr_free(S); expr_free(Pr);
                }
                expr_free(gc);
            }
            expr_free(dL); expr_free(L);
        }
        expr_free(negprod); expr_free(W);
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
    symtab_get_def("DSolve`SpecialFunctionForm")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`SpecialFunctionForm",
        "DSolve`SpecialFunctionForm[eqn, y, x] recognises second-order linear ODEs "
        "whose solutions are named special functions: Airy (y'' == (A x + B) y), "
        "Bessel / modified Bessel (x^2 y'' + x y' +- (x^2 -+ v^2) y == 0), Kummer "
        "confluent hypergeometric (x y'' + (b - x) y' - a y == 0 -> "
        "Hypergeometric1F1), and Gauss hypergeometric "
        "(x(1-x) y'' + (c - (a+b+1) x) y' - a b y == 0 -> Hypergeometric2F1). The "
        "hypergeometric second solution is emitted only when b (resp. c) is not an "
        "integer; otherwise it declines to the series fallback.");
}
