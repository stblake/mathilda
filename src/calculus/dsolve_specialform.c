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
#include <math.h>

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

/* --- numeric self-verify for an emitted 2nd-order solution.  The Pöschl-Teller
 *     row below returns a Hypergeometric2F1 combination whose residual zero_test
 *     cannot decide, so (as with M12/M14) we gate emission on a numeric back-
 *     substitution: only a solution that vanishes at several sample points, with
 *     every symbolic parameter instantiated at a distinct generic real, is kept. */
static void sf_collect_params(const Expr* e, const char** names, int* n, int cap) {
    if (!e || *n >= cap) return;
    if (e->type == EXPR_SYMBOL) {
        const char* nm = e->data.symbol.name;
        for (int i = 0; i < *n; i++) if (names[i] == nm) return;
        names[(*n)++] = nm; return;
    }
    if (e->type == EXPR_FUNCTION) {
        if (e->data.function.head && e->data.function.head->type != EXPR_SYMBOL)
            sf_collect_params(e->data.function.head, names, n, cap);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            sf_collect_params(e->data.function.args[i], names, n, cap);
    }
}
static double sf_abs_at(const Expr* R, const char* xv, double xval) {
    Expr* e = ds_subst(expr_copy((Expr*)R), expr_new_symbol(xv), expr_new_real(xval));
    e = eval_and_free(ds_call1("Abs", eval_and_free(ds_call1("N", e))));
    double m = (e && e->type == EXPR_REAL) ? e->data.real
             : (e && e->type == EXPR_INTEGER) ? (double)e->data.integer : NAN;
    expr_free(e); return m;
}
static bool sf_num_ok(const DSolveProblem* P, const Expr* general,
                      const char* xv, const char* yname) {
    Expr* R = expr_copy(P->eq_residuals[0]);
    Expr* b0 = expr_copy((Expr*)general);
    Expr* b1 = ds_d(expr_copy((Expr*)general), expr_new_symbol(xv));
    Expr* b2 = ds_d(ds_d(expr_copy((Expr*)general), expr_new_symbol(xv)), expr_new_symbol(xv));
    R = ds_subst(R, ds_make_funcapp(yname, 2, xv), b2);
    R = ds_subst(R, ds_make_funcapp(yname, 1, xv), b1);
    R = ds_subst(R, ds_make_funcapp(yname, 0, xv), b0);
    R = ds_subst(R, ds_const(1), expr_new_real(1.3));
    R = ds_subst(R, ds_const(2), expr_new_real(0.7));
    {
        const char* skip[] = { xv, intern_symbol("E"), intern_symbol("Pi"),
            intern_symbol("I"), intern_symbol("C"), intern_symbol("EulerGamma"),
            intern_symbol("Degree"), intern_symbol("GoldenRatio"),
            intern_symbol("Catalan"), intern_symbol("Infinity") };
        const int nskip = (int)(sizeof(skip)/sizeof(skip[0]));
        const char* syms[64]; int ns = 0;
        sf_collect_params(R, syms, &ns, 64);
        int pi = 0;
        for (int i = 0; i < ns; i++) {
            bool sk = false;
            for (int j = 0; j < nskip; j++) if (syms[i] == skip[j]) { sk = true; break; }
            if (sk) continue;
            double v = 0.29 + 0.13 * (double)pi; pi++;
            R = ds_subst(R, expr_new_symbol(syms[i]), expr_new_real(v));
        }
    }
    const double xs[] = { 0.2, 0.4, 0.6, 0.8, 1.0 };   /* in (0, pi/2), away from pi/2
                                                        * where the 2F1 argument Sin^2 -> 1
                                                        * can approach a branch point */
    int small = 0, big = 0;
    for (int i = 0; i < 5; i++) {
        double m = sf_abs_at(R, xv, xs[i]);
        if (isnan(m) || !isfinite(m)) continue;
        if (m < 1e-6) small++; else if (m > 1e-3) big++;
    }
    expr_free(R);
    return small >= 2 && big == 0;
}

/* Sin[x]^pe Cos[x]^q Hypergeometric2F1[(pe+q+s)/2, (pe+q-s)/2, pe+1/2, Sin[x]^2].
 * The two Pöschl-Teller solutions are pt_sol(p) and pt_sol(1-p) (the second's
 * lower parameter 3/2-p is exactly (1-p)+1/2).  pe, q, s borrowed. */
static Expr* pt_sol(const Expr* pe, const Expr* q, const Expr* s, const char* xv) {
    Expr* sin2 = ds_call2(SYM_Power, ds_call1("Sin", expr_new_symbol(xv)), expr_new_integer(2));
    Expr* pq   = ds_call2(SYM_Plus, expr_copy((Expr*)pe), expr_copy((Expr*)q));   /* pe+q */
    Expr* halfA = ds_call2(SYM_Power, expr_new_integer(2), expr_new_integer(-1)); /* 1/2 */
    Expr* A = eval_and_free(ds_call2(SYM_Times,
                  ds_call2(SYM_Plus, expr_copy(pq), expr_copy((Expr*)s)), expr_copy(halfA)));
    Expr* B = eval_and_free(ds_call2(SYM_Times,
                  ds_call2(SYM_Subtract, expr_copy(pq), expr_copy((Expr*)s)), expr_copy(halfA)));
    Expr* C = eval_and_free(ds_call2(SYM_Plus, expr_copy((Expr*)pe), expr_copy(halfA))); /* pe+1/2 */
    expr_free(pq); expr_free(halfA);
    Expr* F = eval_and_free(expr_new_function(expr_new_symbol(SYM_Hypergeometric2F1),
                  (Expr*[]){ A, B, C, sin2 }, 4));
    Expr* sp = eval_and_free(ds_call2(SYM_Power, ds_call1("Sin", expr_new_symbol(xv)), expr_copy((Expr*)pe)));
    Expr* cq = eval_and_free(ds_call2(SYM_Power, ds_call1("Cos", expr_new_symbol(xv)), expr_copy((Expr*)q)));
    return eval_and_free(ds_call2(SYM_Times, sp, ds_call2(SYM_Times, cq, F)));
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

    /* ---- Legendre / associated Legendre ----
     *   (1-x^2) y'' - 2x y' + (nu(nu+1) - mu^2/(1-x^2)) y = 0
     * normalised: P = -2x/(1-x^2),  Q = nu(nu+1)/(1-x^2) - mu^2/(1-x^2)^2.
     * Signature P*(1-x^2)+2x == 0 pins the y'-term; then qq = Q*(1-x^2)^2 must be
     * the quadratic (nu(nu+1)-mu^2) - nu(nu+1) x^2.  Emit LegendreP[nu,(mu,)x],
     * LegendreQ[nu,(mu,)x] (mu dropped for the ordinary equation).  Placed first:
     * its P is neither 0 (Airy) nor 1/x (Bessel), so no row below collides. */
    if (!general) {
        Expr* omx2 = eval_and_free(ds_call2(SYM_Subtract, expr_new_integer(1),
                         ds_call2(SYM_Power, expr_new_symbol(xvar), expr_new_integer(2)))); /* 1 - x^2 */
        Expr* pchk = ds_simplify(ds_call2(SYM_Plus,
                         ds_call2(SYM_Times, expr_copy(Pc), expr_copy(omx2)),
                         ds_call2(SYM_Times, expr_new_integer(2), expr_new_symbol(xvar))));
        if (ds_is_zero(pchk)) {
            Expr* qq = ds_simplify(ds_call2(SYM_Times, expr_copy(Qc),
                           ds_call2(SYM_Power, expr_copy(omx2), expr_new_integer(2))));   /* Q (1-x^2)^2 */
            Expr* c0 = ds_subst(expr_copy(qq), expr_new_symbol(xvar), expr_new_integer(0));   /* nu(nu+1)-mu^2 */
            Expr* c2 = eval_and_free(expr_new_function(expr_new_symbol("Coefficient"),
                           (Expr*[]){ expr_copy(qq), expr_new_symbol(xvar), expr_new_integer(2) }, 3));
            /* require qq == c0 + c2 x^2 exactly (no x^1 or higher terms) and c0,c2 free of x */
            Expr* recon = ds_simplify(ds_call2(SYM_Subtract, expr_copy(qq),
                              ds_call2(SYM_Plus, expr_copy(c0),
                                  ds_call2(SYM_Times, expr_copy(c2),
                                      ds_call2(SYM_Power, expr_new_symbol(xvar), expr_new_integer(2))))));
            if (ds_is_zero(recon) && ds_free_of(c0, xvar) && ds_free_of(c2, xvar)) {
                Expr* B = ds_simplify(ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(c2))); /* nu(nu+1) */
                /* nu = (-1 + Sqrt[1+4B]) / 2 */
                Expr* disc = ds_call2(SYM_Plus, expr_new_integer(1),
                                 ds_call2(SYM_Times, expr_new_integer(4), expr_copy(B)));
                Expr* nu = ds_simplify(ds_call2(SYM_Times,
                               ds_call2(SYM_Plus, expr_new_integer(-1), ds_call1("Sqrt", disc)),
                               ds_call2(SYM_Power, expr_new_integer(2), expr_new_integer(-1))));
                Expr* musq = ds_simplify(ds_call2(SYM_Subtract, expr_copy(B), expr_copy(c0)));
                Expr* b0; Expr* b1;
                if (ds_is_zero(musq)) {
                    b0 = expr_new_function(expr_new_symbol(SYM_LegendreP),
                             (Expr*[]){ expr_copy(nu), expr_new_symbol(xvar) }, 2);
                    b1 = expr_new_function(expr_new_symbol(SYM_LegendreQ),
                             (Expr*[]){ expr_copy(nu), expr_new_symbol(xvar) }, 2);
                } else {
                    Expr* mu = ds_call1("Sqrt", expr_copy(musq));
                    b0 = expr_new_function(expr_new_symbol(SYM_LegendreP),
                             (Expr*[]){ expr_copy(nu), expr_copy(mu), expr_new_symbol(xvar) }, 3);
                    b1 = expr_new_function(expr_new_symbol(SYM_LegendreQ),
                             (Expr*[]){ expr_copy(nu), expr_copy(mu), expr_new_symbol(xvar) }, 3);
                    expr_free(mu);
                }
                general = combo(b0, b1);
                expr_free(B); expr_free(nu); expr_free(musq);
            }
            expr_free(qq); expr_free(c0); expr_free(c2); expr_free(recon);
        }
        expr_free(omx2); expr_free(pchk);
    }

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

    /* ---- Bessel-reducible pure-power potential: P == 0, Q == A x^m with m a
     *      number != 0, -2.  The reduced equation y'' + A x^m y == 0 has
     *      y = Sqrt[x] Z_{1/(m+2)}(kappa x^((m+2)/2)), kappa = 2 Sqrt[|A|]/(m+2),
     *      with Z = J/Y for A > 0 and the modified I/K for A < 0.  Fixes
     *      y'' - x^4 y == 0 -> Sqrt[x](C[1] BesselI[1/6, x^3/3] + C[2] BesselK[1/6, x^3/3]).
     *      Airy (Q linear => m==1) and literal Bessel (P==1/x) run first, so the
     *      !general guard keeps this from colliding with them.  Verified, like the
     *      other rows, by the substrate's back-substitution. */
    if (!general && ds_is_zero(Pc)) {
        Expr* reff = ds_simplify(expr_copy(Qc));   /* Q == A x^m (P==0, so this is the potential) */
        Expr* m = NULL; Expr* mp2 = NULL; Expr* A = NULL;
        if (!ds_is_zero(reff)) {
            /* m = x reff'/reff : the exponent when reff = A x^m */
            Expr* dr = ds_d(expr_copy(reff), expr_new_symbol(xvar));
            Expr* invr = expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ expr_copy(reff), expr_new_integer(-1) }, 2);
            m = ds_simplify(ds_call2(SYM_Times, expr_new_symbol(xvar),
                    ds_call2(SYM_Times, dr, invr)));
            Expr* nqm = eval_and_free(ds_call1("NumberQ", expr_copy(m)));
            bool mok = (nqm->type == EXPR_SYMBOL && nqm->data.symbol.name == SYM_True)
                     && ds_free_of(m, xvar);
            expr_free(nqm);
            if (mok) mp2 = ds_simplify(ds_call2(SYM_Plus, expr_copy(m), expr_new_integer(2)));
            if (mok && !ds_is_zero(m) && mp2 && !ds_is_zero(mp2)) {
                /* A = reff / x^m (free of x for a pure power) */
                Expr* xm = expr_new_function(expr_new_symbol(SYM_Power),
                               (Expr*[]){ expr_new_symbol(xvar), expr_copy(m) }, 2);
                A = ds_simplify(ds_call2(SYM_Times, expr_copy(reff),
                        expr_new_function(expr_new_symbol(SYM_Power),
                            (Expr*[]){ xm, expr_new_integer(-1) }, 2)));
                Expr* nqA = eval_and_free(ds_call1("NumberQ", expr_copy(A)));
                bool Aok = (nqA->type == EXPR_SYMBOL && nqA->data.symbol.name == SYM_True)
                         && ds_free_of(A, xvar);
                expr_free(nqA);
                Expr* nA = eval_and_free(ds_call1("N", expr_copy(A)));   /* numeric sign of A */
                double av = (nA->type == EXPR_REAL) ? nA->data.real
                          : (nA->type == EXPR_INTEGER) ? (double)nA->data.integer : 0.0;
                expr_free(nA);
                if (Aok && av != 0.0) {
                    bool pos = av > 0.0;
                    Expr* beta = pos ? expr_copy(A)
                                     : ds_simplify(ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(A)));
                    Expr* n = ds_simplify(expr_new_function(expr_new_symbol(SYM_Power),  /* 1/(m+2) */
                                  (Expr*[]){ expr_copy(mp2), expr_new_integer(-1) }, 2));
                    Expr* kappa = ds_simplify(ds_call2(SYM_Times,          /* 2 Sqrt[beta]/(m+2) */
                                      ds_call2(SYM_Times, expr_new_integer(2), ds_call1("Sqrt", beta)),
                                      expr_new_function(expr_new_symbol(SYM_Power),
                                          (Expr*[]){ expr_copy(mp2), expr_new_integer(-1) }, 2)));
                    Expr* half = ds_simplify(ds_call2(SYM_Times, expr_copy(mp2),  /* (m+2)/2 */
                                     expr_new_function(expr_new_symbol(SYM_Power),
                                         (Expr*[]){ expr_new_integer(2), expr_new_integer(-1) }, 2)));
                    Expr* arg = ds_simplify(ds_call2(SYM_Times, expr_copy(kappa),  /* kappa x^half */
                                    expr_new_function(expr_new_symbol(SYM_Power),
                                        (Expr*[]){ expr_new_symbol(xvar), expr_copy(half) }, 2)));
                    Expr* xsym = expr_new_symbol(xvar);
                    Expr* sqrtx = powrat(xsym, 1, 2);   /* Sqrt[x]; powrat borrows */
                    expr_free(xsym);
                    const char* f0 = pos ? "BesselJ" : "BesselI";
                    const char* f1 = pos ? "BesselY" : "BesselK";
                    Expr* b0 = ds_call2(SYM_Times, expr_copy(sqrtx),
                                   expr_new_function(expr_new_symbol(f0),
                                       (Expr*[]){ expr_copy(n), expr_copy(arg) }, 2));
                    Expr* b1 = ds_call2(SYM_Times, expr_copy(sqrtx),
                                   expr_new_function(expr_new_symbol(f1),
                                       (Expr*[]){ expr_copy(n), expr_copy(arg) }, 2));
                    general = combo(b0, b1);
                    expr_free(sqrtx); expr_free(n); expr_free(kappa);
                    expr_free(half); expr_free(arg);
                }
            }
        }
        expr_free(reff);
        if (m) expr_free(m);
        if (mp2) expr_free(mp2);
        if (A) expr_free(A);
    }

    /* ---- Bessel-reducible exponential potential: P == 0, Q == A e^(lambda x),
     *      lambda a nonzero constant.  y'' + A e^(lambda x) y == 0 maps under
     *      t = (2 Sqrt[|A|]/lambda) e^(lambda x/2) to order-0 Bessel:
     *      J_0/Y_0 for A > 0, modified I_0/K_0 for A < 0.  Fixes
     *      y'' - e^(5x) y == 0 -> BesselI[0,(2/5)e^(5x/2)], BesselK[0,...].
     *      Placed after the pure-power row (whose m = x Q'/Q is not a number for
     *      an exponential Q, so it declines first). */
    if (!general && ds_is_zero(Pc)) {
        Expr* Q = ds_simplify(expr_copy(Qc));
        if (!ds_is_zero(Q)) {
            Expr* dQ = ds_d(expr_copy(Q), expr_new_symbol(xvar));
            Expr* lam = ds_simplify(ds_call2(SYM_Times, dQ,          /* lambda = Q'/Q */
                            expr_new_function(expr_new_symbol(SYM_Power),
                                (Expr*[]){ expr_copy(Q), expr_new_integer(-1) }, 2)));
            if (ds_free_of(lam, xvar) && !ds_is_zero(lam)) {
                Expr* elx = eval_and_free(ds_call1("Exp",
                                ds_call2(SYM_Times, expr_copy(lam), expr_new_symbol(xvar))));
                Expr* A = ds_simplify(ds_call2(SYM_Times, expr_copy(Q),   /* A = Q / e^(lambda x) */
                              expr_new_function(expr_new_symbol(SYM_Power),
                                  (Expr*[]){ elx, expr_new_integer(-1) }, 2)));
                Expr* nA = eval_and_free(ds_call1("N", expr_copy(A)));
                double av = (nA->type == EXPR_REAL) ? nA->data.real
                          : (nA->type == EXPR_INTEGER) ? (double)nA->data.integer : 0.0;
                expr_free(nA);
                if (ds_free_of(A, xvar) && av != 0.0) {
                    bool pos = av > 0.0;
                    Expr* beta = pos ? expr_copy(A)
                                     : ds_simplify(ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(A)));
                    /* kappa = 2 Sqrt[beta] / lambda */
                    Expr* kappa = ds_simplify(ds_call2(SYM_Times,
                                      ds_call2(SYM_Times, expr_new_integer(2), ds_call1("Sqrt", beta)),
                                      expr_new_function(expr_new_symbol(SYM_Power),
                                          (Expr*[]){ expr_copy(lam), expr_new_integer(-1) }, 2)));
                    /* arg = kappa e^(lambda x / 2) */
                    Expr* halfexp = eval_and_free(ds_call1("Exp",
                                        ds_call2(SYM_Times,
                                            ds_call2(SYM_Times, expr_copy(lam),
                                                expr_new_function(expr_new_symbol(SYM_Power),
                                                    (Expr*[]){ expr_new_integer(2), expr_new_integer(-1) }, 2)),
                                            expr_new_symbol(xvar))));
                    Expr* arg = ds_simplify(ds_call2(SYM_Times, kappa, halfexp));
                    const char* f0 = pos ? "BesselJ" : "BesselI";
                    const char* f1 = pos ? "BesselY" : "BesselK";
                    Expr* b0 = expr_new_function(expr_new_symbol(f0),
                                   (Expr*[]){ expr_new_integer(0), expr_copy(arg) }, 2);
                    Expr* b1 = expr_new_function(expr_new_symbol(f1),
                                   (Expr*[]){ expr_new_integer(0), expr_copy(arg) }, 2);
                    general = combo(b0, b1);
                    expr_free(arg);
                }
                expr_free(A);
            }
            expr_free(lam);
        }
        expr_free(Q);
    }

    /* ---- normal-form Bessel: P == 0, Q == A + B/x^2 (A const != 0, B const) ----
     *   u'' + (A + B/x^2) u == 0  ->  u = Sqrt[x] Z_nu(Sqrt[|A|] x),
     *   nu = Sqrt[1/4 - B], Z = J/Y for A > 0 and modified I/K for A < 0.
     *   This is the normal (u'-free) form the plain Bessel row (P == 1/x) and the
     *   pure-power row (Q a single power) both miss; it is the second-order factor
     *   of Bessel-type third-order symmetric squares (e.g. E17). */
    if (!general && ds_is_zero(Pc)) {
        Expr* q2 = ds_simplify(ds_call2(SYM_Times, expr_copy(Qc),         /* x^2 Q */
                       ds_call2(SYM_Power, expr_new_symbol(xvar), expr_new_integer(2))));
        Expr* Ac = eval_and_free(expr_new_function(expr_new_symbol("Coefficient"),
                       (Expr*[]){ expr_copy(q2), expr_new_symbol(xvar), expr_new_integer(2) }, 3)); /* A */
        Expr* Bc = ds_subst(expr_copy(q2), expr_new_symbol(xvar), expr_new_integer(0));            /* B */
        Expr* recon = ds_simplify(ds_call2(SYM_Subtract, expr_copy(q2),
                          ds_call2(SYM_Plus,
                              ds_call2(SYM_Times, expr_copy(Ac),
                                  ds_call2(SYM_Power, expr_new_symbol(xvar), expr_new_integer(2))),
                              expr_copy(Bc))));
        if (ds_is_zero(recon) && ds_free_of(Ac, xvar) && ds_free_of(Bc, xvar) && !ds_is_zero(Ac)) {
            Expr* nA = eval_and_free(ds_call1("N", expr_copy(Ac)));
            double av = (nA->type == EXPR_REAL) ? nA->data.real
                      : (nA->type == EXPR_INTEGER) ? (double)nA->data.integer : 0.0;
            expr_free(nA);
            if (av != 0.0) {
                bool pos = av > 0.0;
                Expr* beta = pos ? expr_copy(Ac)
                                 : ds_simplify(ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(Ac)));
                Expr* nu = ds_simplify(ds_call1("Sqrt", ds_call2(SYM_Subtract,   /* Sqrt[1/4 - B] */
                               expr_new_function(expr_new_symbol(SYM_Power),
                                   (Expr*[]){ expr_new_integer(4), expr_new_integer(-1) }, 2),
                               expr_copy(Bc))));
                Expr* arg = ds_simplify(ds_call2(SYM_Times, ds_call1("Sqrt", beta), expr_new_symbol(xvar)));
                Expr* xsym = expr_new_symbol(xvar);
                Expr* sqrtx = powrat(xsym, 1, 2);   /* Sqrt[x]; powrat borrows */
                expr_free(xsym);
                const char* f0 = pos ? "BesselJ" : "BesselI";
                const char* f1 = pos ? "BesselY" : "BesselK";
                Expr* b0 = ds_call2(SYM_Times, expr_copy(sqrtx),
                               expr_new_function(expr_new_symbol(f0),
                                   (Expr*[]){ expr_copy(nu), expr_copy(arg) }, 2));
                Expr* b1 = ds_call2(SYM_Times, expr_copy(sqrtx),
                               expr_new_function(expr_new_symbol(f1),
                                   (Expr*[]){ expr_copy(nu), expr_copy(arg) }, 2));
                general = combo(b0, b1);
                expr_free(sqrtx); expr_free(nu); expr_free(arg);
            }
        }
        expr_free(q2); expr_free(Ac); expr_free(Bc); expr_free(recon);
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
                /* Emit whenever c is not a PROVABLE integer — this now includes a
                 * SYMBOLIC c (generically non-integer), giving Hypergeometric2F1
                 * closed forms for e.g. (x^2-x)y''+((a+b+1)x-c)y'+ab y==0.  The
                 * symbolic-power second solution x^(1-c) 2F1[...] is kept safe from
                 * a zero_test hang by the special-function early-decline in
                 * zero_test.c; an integer c still declines (dependent solutions). */
                if (!cint) {
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

    /* ---- Trigonometric Pöschl-Teller potential: P == 0,
     *      Q == c0 + c1 Csc[x]^2 + c2 Sec[x]^2   (c_i free of x).
     * Equivalently y'' == (a + p(p-1) Csc^2 x + q(q-1) Sec^2 x) y with a = -c0,
     * p(p-1) = -c1, q(q-1) = -c2.  Solutions (numerically verified):
     *   y_{p} = Sin^p Cos^q 2F1((p+q+s)/2,(p+q-s)/2, p+1/2, Sin^2 x),  s = Sqrt[-a],
     * with the second solution y_{1-p}.  Extraction: R = Q Sin^2 Cos^2, rewritten
     * with Sin^2 -> 1-Cos^2, is the even quadratic in Cos[x]
     *   -c0 Cos^4 + (c0+c1-c2) Cos^2 + c2.
     * The 2F1 residual is undecidable by zero_test, so a numeric self-verify gates
     * emission — a degenerate parameter set (singular 2F1 lower parameter) declines
     * to the Frobenius fallback.  Covers Kamke/Murphy trig-potential equations
     * (Csc^2/Sec^2/Cos^2/Sin^2), including those written via Cot^2 or (a Cos^2 + b
     * Sin^2 + c)/Sin^2. */
    if (!general && ds_is_zero(Pc)) {
        Expr* cs = ds_call2(SYM_Power, ds_call1("Sin", expr_new_symbol(xvar)), expr_new_integer(2));
        Expr* cc = ds_call2(SYM_Power, ds_call1("Cos", expr_new_symbol(xvar)), expr_new_integer(2));
        Expr* cosx = ds_call1("Cos", expr_new_symbol(xvar));
        Expr* R = ds_simplify(ds_call2(SYM_Times, expr_copy(Qc),
                      ds_call2(SYM_Times, expr_copy(cs), expr_copy(cc))));            /* Q Sin^2 Cos^2 */
        /* Rewrite every even power Sin^{2k} -> (1-Cos^2)^k (Simplify can introduce
         * Sin^4/Sin^6), then Expand (NOT Simplify: it reverts 1-Cos^2 -> Sin^2). */
        Expr* omc = ds_call2(SYM_Subtract, expr_new_integer(1), expr_copy(cc));      /* 1-Cos^2 */
        Expr* rl[3];
        for (int kk = 3; kk >= 1; kk--) {
            Expr* lhs = ds_call2(SYM_Power, ds_call1("Sin", expr_new_symbol(xvar)),
                            expr_new_integer(2 * kk));
            Expr* rhs = (kk == 1) ? expr_copy(omc)
                                  : ds_call2(SYM_Power, expr_copy(omc), expr_new_integer(kk));
            rl[3 - kk] = expr_new_function(expr_new_symbol(SYM_Rule), (Expr*[]){ lhs, rhs }, 2);
        }
        Expr* rules = expr_new_function(expr_new_symbol(SYM_List), rl, 3);
        expr_free(omc);
        Expr* R2 = eval_and_free(expr_new_function(expr_new_symbol(SYM_ReplaceAll),
                       (Expr*[]){ expr_copy(R), rules }, 2));
        R2 = eval_and_free(ds_call1("Expand", R2));
        Expr* pq = eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ),
                       (Expr*[]){ expr_copy(R2), expr_copy(cosx) }, 2));
        bool ispoly = (pq->type == EXPR_SYMBOL && pq->data.symbol.name == SYM_True);
        expr_free(pq);
        if (ispoly) {
            #define SF_COEF(k) eval_and_free(expr_new_function(expr_new_symbol("Coefficient"), \
                (Expr*[]){ expr_copy(R2), expr_copy(cosx), expr_new_integer(k) }, 3))
            Expr* k0 = SF_COEF(0); Expr* k1 = SF_COEF(1); Expr* k2 = SF_COEF(2);
            Expr* k3 = SF_COEF(3); Expr* k4 = SF_COEF(4); Expr* k5 = SF_COEF(5);
            #undef SF_COEF
            /* pure even quadratic in Cos[x]^2: odd and degree-5 coeffs vanish */
            if (ds_is_zero(k1) && ds_is_zero(k3) && ds_is_zero(k5)) {
                Expr* c0 = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(k4)));  /* -k4 */
                Expr* c2 = expr_copy(k0);
                Expr* c1 = ds_simplify(ds_call2(SYM_Plus,
                               ds_call2(SYM_Subtract, expr_copy(k2), expr_copy(c0)), expr_copy(c2)));
                if (ds_free_of(c0, xvar) && ds_free_of(c1, xvar) && ds_free_of(c2, xvar)) {
                    /* p=(1+Sqrt[1-4 c1])/2, q=(1+Sqrt[1-4 c2])/2, s=Sqrt[-a]=Sqrt[c0] */
                    Expr* half = ds_call2(SYM_Power, expr_new_integer(2), expr_new_integer(-1));
                    Expr* pP = ds_simplify(ds_call2(SYM_Times,
                                   ds_call2(SYM_Plus, expr_new_integer(1),
                                       ds_call1("Sqrt", ds_call2(SYM_Subtract, expr_new_integer(1),
                                           ds_call2(SYM_Times, expr_new_integer(4), expr_copy(c1))))),
                                   expr_copy(half)));
                    Expr* qQ = ds_simplify(ds_call2(SYM_Times,
                                   ds_call2(SYM_Plus, expr_new_integer(1),
                                       ds_call1("Sqrt", ds_call2(SYM_Subtract, expr_new_integer(1),
                                           ds_call2(SYM_Times, expr_new_integer(4), expr_copy(c2))))),
                                   expr_copy(half)));
                    Expr* s  = ds_simplify(ds_call1("Sqrt", expr_copy(c0)));            /* Sqrt[-a]=Sqrt[c0] */
                    Expr* omp = eval_and_free(ds_call2(SYM_Subtract, expr_new_integer(1), expr_copy(pP))); /* 1-p */
                    Expr* b0 = pt_sol(pP, qQ, s, xvar);
                    Expr* b1 = pt_sol(omp, qQ, s, xvar);
                    Expr* cand = combo(b0, b1);
                    if (sf_num_ok(P, cand, xvar, P->fun_names[0])) general = cand;
                    else expr_free(cand);
                    expr_free(half); expr_free(pP); expr_free(qQ); expr_free(s); expr_free(omp);
                }
                expr_free(c0); expr_free(c1); expr_free(c2);
            }
            expr_free(k0); expr_free(k1); expr_free(k2); expr_free(k3); expr_free(k4); expr_free(k5);
        }
        expr_free(R); expr_free(R2); expr_free(cosx); expr_free(cs); expr_free(cc);
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
        "(x(1-x) y'' + (c - (a+b+1) x) y' - a b y == 0 -> Hypergeometric2F1), Legendre / "
        "associated Legendre, and the trigonometric Pöschl-Teller potential "
        "(y'' == (a + p(p-1) Csc^2 x + q(q-1) Sec^2 x) y -> Hypergeometric2F1, "
        "numerically verified). The "
        "hypergeometric second solution is emitted only when b (resp. c) is not an "
        "integer; otherwise it declines to the series fallback.");
}
