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
 *   Kummer (1F1)    P = b/x-s,  Q = -a s/x      ->  Hypergeometric1F1[a,b,s x],
 *                                                   (s x)^(1-b) 1F1[a-b+1, 2-b, s x]
 *                                                   (integer b: HypergeometricU[a,b,s x]
 *                                                   replaces the dependent one; when the
 *                                                   1F1 is a polynomial, U is too, and the
 *                                                   ExpIntegralEi partner is used instead)
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
#include "dsolve.h"          /* g_dsolve_depth */
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../core.h"
#include <stdlib.h>
#include <math.h>

/* base^-1 (unevaluated); base consumed. */
static Expr* powneg1_sf(Expr* base) {
    return expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ base, expr_new_integer(-1) }, 2);
}

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

/* Cancel[Together[e]] — rational-function reduction in the variable that does NOT
 * attempt to simplify radical/parametric CONSTANTS.  The affine/F-homotopy path
 * feeds coefficients that are rational in the mapped variable but carry radical
 * exponents in several symbolic parameters; full Simplify then blows up (minutes,
 * effectively a hang) trying to canonicalise those radicals, whereas the only
 * reduction actually needed is the rational-in-s one.  e consumed. */
static Expr* sf_ct(Expr* e) {
    return eval_and_free(ds_call1("Cancel", ds_call1("Together", e)));
}

/* True if `e` contains a radical: a Sqrt, or a Power with a non-integer exponent.
 * Used to decline the affine->Gauss row when the regular singular points x1, x2 are
 * radical (roots of an irreducible quadratic leading coefficient): the affine map and
 * the subsequent indicial FactorList/Solve then run on radical-parametric coefficients
 * and take seconds-to-minutes, whereas every genuine target has RATIONAL singular
 * points (Gegenbauer/Jacobi/associated Legendre at +-1, shifted Gauss at rational
 * endpoints).  Such radical-RSP cases are not closed by this recognizer anyway, so
 * declining fast lets the cascade reach the series fallback that used to solve them. */
static bool sf_has_radical(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h && h->type == EXPR_SYMBOL) {
        if (h->data.symbol.name == intern_symbol("Sqrt")) return true;
        if (h->data.symbol.name == SYM_Power && e->data.function.arg_count == 2 &&
            e->data.function.args[1]->type != EXPR_INTEGER) return true;   /* rational/symbolic exponent */
    }
    if (h && sf_has_radical(h)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (sf_has_radical(e->data.function.args[i])) return true;
    return false;
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

/* Exact rational q as an evaluated Expr (Integer or Rational). */
static Expr* sf_mpq_expr(const mpq_t q) {
    Expr* num = expr_bigint_normalize(expr_new_bigint_from_mpz(mpq_numref(q)));
    if (mpz_cmp_ui(mpq_denref(q), 1) == 0) return num;
    Expr* den = expr_bigint_normalize(expr_new_bigint_from_mpz(mpq_denref(q)));
    return eval_and_free(ds_call2(SYM_Times, num, powneg1_sf(den)));
}

/* Upper bound on n + B for the logarithmic Kummer partner below: the coefficient
 * recurrence is O(n + B) exact rationals, but the emitted expression (and the
 * back-substitution verify on it) grows linearly too, so an absurd degree
 * declines to the rest of the cascade instead of building a huge closed form. */
#define SF_KUMMER_LOG_MAX 400
static bool sf_small_int(int64_t v) { return v >= -SF_KUMMER_LOG_MAX && v <= SF_KUMMER_LOG_MAX; }

/* The logarithmic (ExpIntegralEi) second solution of the DEGENERATE Kummer
 * equation  z w'' + (B - z) w' + n w == 0,  n >= 0, B >= 1 integers, i.e.
 * a == -n.  There the polynomial p = 1F1[-n, B, z] (a Laguerre polynomial up to
 * a constant) is one solution and Tricomi's U[-n, B, z] = (-1)^n (B)_n p is the
 * SAME solution (Wronskian -Gamma(B) z^-B e^z / Gamma(a) == 0), so the textbook
 * {1F1, U} basis collapses.  The independent partner is
 *
 *     w2 = p(z) Ei(z) + e^z R(z),   R = sum_{j=1-B}^{n-1} r_j z^j   (Laurent),
 *
 * found by substituting into the ODE: Ei' = e^z/z gives
 *     L[p Ei]   = e^z (2 p' + (B-1) p / z)          (L[p] == 0),
 *     L[e^z R]  = e^z (z R'' + (z + B) R' + (B + n) R),
 * so R solves  z R'' + (z+B) R' + (B+n) R == -(2 p' + (B-1) p/z).  Matching z^j:
 *     (j+B+n) r_j + (j+1)(j+B) r_{j+1} == -(2j+B+1) p_{j+1}     (p_k = 0, k < 0),
 * a two-term recurrence solved top-down from r_n = 0.  j+B+n >= n+1 > 0 on the
 * whole range, so every step is a plain division; the coupling vanishes at
 * j = -B, which is why R stops at z^(1-B).  The Ei term makes w2 non-elementary,
 * hence independent of p (Wronskian == e^z z^-B by Abel).  This is the closed
 * form of e^z U[B+n, B, -z] up to a multiple of p, without its complex branch on
 * real z > 0.  z borrowed; NULL when n + B exceeds SF_KUMMER_LOG_MAX. */
static Expr* specialform_kummer_log_partner(int64_t n, int64_t B, const Expr* z) {
    if (n < 0 || B < 1 || n + B > SF_KUMMER_LOG_MAX) return NULL;
    /* p_k = prod_{i<k} (i - n) / ((B + i)(i + 1)),  k = 0..n */
    mpq_t* p = (mpq_t*)malloc((size_t)(n + 1) * sizeof(mpq_t));
    /* r_j stored at r[j + B - 1], j = 1-B .. n  (r_n = 0) */
    mpq_t* r = (mpq_t*)malloc((size_t)(n + B) * sizeof(mpq_t));
    if (!p || !r) { free(p); free(r); return NULL; }
    mpq_t t, u;
    mpq_init(t); mpq_init(u);
    for (int64_t k = 0; k <= n; k++) {
        mpq_init(p[k]);
        if (k == 0) { mpq_set_ui(p[k], 1, 1); continue; }
        int64_t i = k - 1;
        mpz_set_si(mpq_numref(t), (long)(i - n));
        mpz_set_si(mpq_denref(t), (long)((B + i) * (i + 1)));   /* <= 400^2, fits */
        mpq_canonicalize(t);
        mpq_mul(p[k], p[k - 1], t);
    }
    for (int64_t j = 0; j < n + B; j++) mpq_init(r[j]);   /* all 0, incl. r_n */
    for (int64_t j = n - 1; j >= 1 - B; j--) {
        /* rhs_j = -(2j+B+1) p_{j+1}  (zero when j+1 < 0) */
        if (j + 1 >= 0) {
            mpq_set_si(t, (long)(-(2 * j + B + 1)), 1);
            mpq_mul(t, t, p[j + 1]);
        } else mpq_set_ui(t, 0, 1);
        /* - (j+1)(j+B) r_{j+1} */
        mpq_set_si(u, (long)((j + 1) * (j + B)), 1);
        mpq_mul(u, u, r[j + 1 + B - 1]);
        mpq_sub(t, t, u);
        mpq_set_si(u, (long)(j + B + n), 1);
        mpq_div(r[j + B - 1], t, u);
    }
    /* assemble p(z) Ei(z) + e^z R(z) */
    Expr* poly = expr_new_integer(0);
    for (int64_t k = 0; k <= n; k++)
        poly = ds_call2(SYM_Plus, poly, ds_call2(SYM_Times, sf_mpq_expr(p[k]),
                   ds_call2(SYM_Power, expr_copy((Expr*)z), expr_new_integer(k))));
    Expr* R = expr_new_integer(0);
    for (int64_t j = 1 - B; j <= n - 1; j++) {
        if (mpq_sgn(r[j + B - 1]) == 0) continue;
        R = ds_call2(SYM_Plus, R, ds_call2(SYM_Times, sf_mpq_expr(r[j + B - 1]),
                ds_call2(SYM_Power, expr_copy((Expr*)z), expr_new_integer(j))));
    }
    for (int64_t k = 0; k <= n; k++) mpq_clear(p[k]);
    for (int64_t j = 0; j < n + B; j++) mpq_clear(r[j]);
    mpq_clear(t); mpq_clear(u);
    free(p); free(r);
    Expr* ei = ds_call1(SYM_ExpIntegralEi, expr_copy((Expr*)z));
    Expr* ez = ds_call2(SYM_Power, expr_new_symbol(SYM_E), expr_copy((Expr*)z));
    return eval_and_free(ds_call2(SYM_Plus, ds_call2(SYM_Times, eval_and_free(poly), ei),
                                            ds_call2(SYM_Times, ez, eval_and_free(R))));
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

    /* Complexity backstop: FactorList / Solve on t^2 - S t + Pr can hang for minutes
     * when S, Pr are large radical-parametric expressions (the a+b, a b produced by
     * the affine F-homotopy on an equation with several symbolic parameters, e.g. a
     * generic quadratic potential over (1-x^2)^2).  Every genuine target measured has
     * a quad of <= 29 leaves; the pathological cases run to hundreds.  Above the bound
     * we decline (the caller falls through to the cascade / series fallback) rather
     * than risk the hang.  TimeConstrained is avoided here: this runs inside M14's
     * TimeConstrained wrapper, and nesting it is unreliable. */
    if (leaf_count_internal(quad, true) > 80) { expr_free(quad); return false; }

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

/* ---- reduced-form recognisers: the P == 0 rows (Airy + the three Bessel-
 *      reducible potentials), factored so they can run on either the raw Qc
 *      (when the equation already has P == 0) OR the Liouville normal-form
 *      potential -r (when a y' term was killed).  `Qc` is the coefficient of y
 *      in z'' + Qc z == 0 (borrowed); returns C[1] z0 + C[2] z1 or NULL. */
static Expr* specialform_reduced_basis(Expr* Qc, const char* xvar) {
    Expr* general = NULL;

    /* ---- Airy: P == 0, Q = -(A x + B), A = -dQ/dx constant, B = -Q(0) ---- */
    if (!general) {
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

    /* ---- Bessel-reducible pure-power potential: P == 0, Q == A x^m with m a
     *      number != 0, -2.  The reduced equation y'' + A x^m y == 0 has
     *      y = Sqrt[x] Z_{1/(m+2)}(kappa x^((m+2)/2)), kappa = 2 Sqrt[|A|]/(m+2),
     *      with Z = J/Y for A > 0 and the modified I/K for A < 0.  Fixes
     *      y'' - x^4 y == 0 -> Sqrt[x](C[1] BesselI[1/6, x^3/3] + C[2] BesselK[1/6, x^3/3]).
     *      Airy (Q linear => m==1) and literal Bessel (P==1/x) run first, so the
     *      !general guard keeps this from colliding with them.  Verified, like the
     *      other rows, by the substrate's back-substitution. */
    if (!general) {
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
    if (!general) {
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
    if (!general) {
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

    return general;
}

/* WhittakerM[k, m, z] as the verifiable confluent form
 *   Exp[-z/2] z^(m+1/2) Hypergeometric1F1[m - k + 1/2, 1 + 2 m, z].
 * (WhittakerM/W have no evaluator in Mathilda; the 1F1 form numericizes, so it can
 * be back-substitution-verified.)  k, m, z borrowed. */
static Expr* whittaker_M_1F1(const Expr* k, const Expr* m, const Expr* z) {
    Expr* half = ds_call2(SYM_Power, expr_new_integer(2), expr_new_integer(-1));   /* 1/2 */
    Expr* ez = eval_and_free(ds_call1("Exp",
                   ds_call2(SYM_Times,
                       ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(half)),
                       expr_copy((Expr*)z))));                                      /* Exp[-z/2] */
    Expr* zpow = eval_and_free(ds_call2(SYM_Power, expr_copy((Expr*)z),
                     ds_call2(SYM_Plus, expr_copy((Expr*)m), expr_copy(half))));    /* z^(m+1/2) */
    Expr* a1 = eval_and_free(ds_call2(SYM_Plus,
                   ds_call2(SYM_Subtract, expr_copy((Expr*)m), expr_copy((Expr*)k)),
                   expr_copy(half)));                                              /* m - k + 1/2 */
    Expr* b1 = eval_and_free(ds_call2(SYM_Plus, expr_new_integer(1),
                   ds_call2(SYM_Times, expr_new_integer(2), expr_copy((Expr*)m)))); /* 1 + 2 m */
    Expr* f = expr_new_function(expr_new_symbol("Hypergeometric1F1"),
                  (Expr*[]){ a1, b1, expr_copy((Expr*)z) }, 3);                     /* consumes a1,b1,z-copy */
    Expr* r = eval_and_free(ds_call2(SYM_Times, ez, ds_call2(SYM_Times, zpow, f)));
    expr_free(half);
    return r;
}

/* ---- Generalised power-potential recogniser: P == 0, the reduced potential is
 *      Q == A x^m                        (single power, SYMBOLIC exponent/coeff), or
 *      Q == alpha x^(2c) + beta x^(c-1)   (two arbitrary powers whose exponents
 *      satisfy  P_big - 2 P_sm == 2).
 *
 *  Both are Liouville-solvable.  This generalises the numeric-exponent pure-power
 *  Bessel row (specialform_reduced_basis) to a symbolic exponent, and adds the
 *  two-term family, which is NOT Bessel but confluent (Coulomb / Whittaker):
 *
 *    single power  y'' + A x^m y == 0
 *       -> y = Sqrt[x] Z_{1/(m+2)}(kappa x^((m+2)/2)),  kappa = 2 Sqrt[A]/(m+2),
 *          Z = BesselJ/BesselY.  (A's sign is left inside Sqrt[A]; sf_num_ok verifies.)
 *
 *    two-term  y'' + (alpha x^(2c) + beta x^(c-1)) y == 0,  put d = c+1:
 *       under xi = x^d it becomes the Coulomb equation
 *          w_xixi + (c/d)(1/xi) w_xi + (alpha/d^2 + (beta/d^2)/xi) w == 0,
 *       solved by  y = x^((1-d)/2) WhittakerM[kappa, +-mu, z]  with
 *          mu = 1/(2d),  kappa = beta/(2 d Sqrt[-alpha]),  z = (2 Sqrt[-alpha]/d) x^d.
 *
 *  Every corpus member (2.1.2-434/792/803/804 and kin) carries SYMBOLIC exponents,
 *  so the residual is undecidable by zero_test; the caller gates every emission with
 *  sf_num_ok (numeric back-substitution).  Additive terms are grouped by their
 *  (symbolic) exponent e = x T'/T so an un-collected like-power pair (e.g. the
 *  normal form of 803 carries b x^(n-1) and -a n/2 x^(n-1) separately) is combined;
 *  a potential with three or more distinct powers (Euler + power, 820/822) declines.
 *  Returns C[1] z0 + C[2] z1 or NULL.  Qc borrowed. */
static Expr* specialform_power_potential(Expr* Qc, const char* xvar) {
    if (!Qc || ds_is_zero(Qc)) return NULL;
    Expr* x = expr_new_symbol(xvar);

    /* Expand first: a normal-form potential from dsolve_normal_form arrives factored
     * (e.g. -1/4 (a^2 x^(2n) + (-4b+2an) x^(n-1))), whose single Times/Plus node the
     * per-term extraction cannot read.  Expand distributes it to a sum of monomials
     * (like symbolic powers stay UN-combined -- the exponent grouping below fixes
     * that).  Consumed at the end. */
    Expr* Qe = eval_and_free(ds_call1("Expand", expr_copy(Qc)));

    /* additive power terms (borrowed from Qe) */
    Expr* raw[8]; int nraw = 0;
    if (ds_has_head(Qe, SYM_Plus)) {
        if (Qe->data.function.arg_count > 8) { expr_free(x); expr_free(Qe); return NULL; }
        for (size_t i = 0; i < Qe->data.function.arg_count; i++) raw[nraw++] = Qe->data.function.args[i];
    } else {
        raw[nraw++] = Qe;
    }

    /* group by exponent: exps[g] distinct exponent, coefs[g] accumulated coefficient */
    Expr* exps[3]; Expr* coefs[3]; int ng = 0; bool bad = false;
    for (int i = 0; i < nraw && !bad; i++) {
        Expr* T = raw[i];
        Expr* e = ds_simplify(ds_call2(SYM_Times, expr_copy(x),                     /* e = x T'/T */
                      ds_call2(SYM_Times, ds_d(expr_copy(T), expr_copy(x)),
                          ds_call2(SYM_Power, expr_copy(T), expr_new_integer(-1)))));
        Expr* co = ds_simplify(ds_call2(SYM_Times, expr_copy(T),                    /* coeff = T / x^e */
                       ds_call2(SYM_Power, ds_call2(SYM_Power, expr_copy(x), expr_copy(e)),
                           expr_new_integer(-1))));
        if (!ds_free_of(e, xvar) || !ds_free_of(co, xvar)) { expr_free(e); expr_free(co); bad = true; break; }
        int j = -1;
        for (int g = 0; g < ng; g++) {
            Expr* diff = ds_simplify(ds_call2(SYM_Subtract, expr_copy(e), expr_copy(exps[g])));
            bool same = ds_is_zero(diff); expr_free(diff);
            if (same) { j = g; break; }
        }
        if (j >= 0) {
            coefs[j] = ds_simplify(ds_call2(SYM_Plus, coefs[j], co));   /* consumes old coefs[j], co */
            expr_free(e);
        } else if (ng < 3) {
            exps[ng] = e; coefs[ng] = co; ng++;
        } else {
            expr_free(e); expr_free(co); bad = true; break;
        }
    }

    Expr* general = NULL;

    if (!bad && ng == 1) {
        /* single power: Q = A x^m -> Sqrt[x] Bessel */
        Expr* m = exps[0]; Expr* A = coefs[0];
        Expr* mp2 = ds_simplify(ds_call2(SYM_Plus, expr_copy(m), expr_new_integer(2)));   /* m+2 */
        if (!ds_is_zero(mp2)) {
            Expr* nu = ds_simplify(ds_call2(SYM_Power, expr_copy(mp2), expr_new_integer(-1)));  /* 1/(m+2) */
            Expr* kappa = ds_simplify(ds_call2(SYM_Times,                                  /* 2 Sqrt[A]/(m+2) */
                              ds_call2(SYM_Times, expr_new_integer(2), ds_call1("Sqrt", expr_copy(A))),
                              ds_call2(SYM_Power, expr_copy(mp2), expr_new_integer(-1))));
            Expr* halfp = ds_simplify(ds_call2(SYM_Times, expr_copy(mp2),                   /* (m+2)/2 */
                              ds_call2(SYM_Power, expr_new_integer(2), expr_new_integer(-1))));
            Expr* arg = ds_simplify(ds_call2(SYM_Times, expr_copy(kappa),                   /* kappa x^((m+2)/2) */
                            ds_call2(SYM_Power, expr_copy(x), expr_copy(halfp))));
            Expr* sqrtx = powrat(x, 1, 2);
            Expr* b0 = ds_call2(SYM_Times, expr_copy(sqrtx),
                           expr_new_function(expr_new_symbol("BesselJ"),
                               (Expr*[]){ expr_copy(nu), expr_copy(arg) }, 2));
            Expr* b1 = ds_call2(SYM_Times, expr_copy(sqrtx),
                           expr_new_function(expr_new_symbol("BesselY"),
                               (Expr*[]){ expr_copy(nu), expr_copy(arg) }, 2));
            general = combo(b0, b1);
            expr_free(sqrtx); expr_free(nu); expr_free(kappa); expr_free(halfp); expr_free(arg);
        }
        expr_free(mp2);
    } else if (!bad && ng == 2) {
        /* two-term: pick P_big, P_sm with P_big - 2 P_sm == 2, then c = P_sm + 1 */
        Expr* Ea = exps[0]; Expr* Eb = exps[1];
        Expr* ca = ds_simplify(ds_call2(SYM_Subtract,
                       ds_call2(SYM_Subtract, expr_copy(Ea),
                           ds_call2(SYM_Times, expr_new_integer(2), expr_copy(Eb))),
                       expr_new_integer(2)));                                 /* Ea - 2 Eb - 2 */
        Expr* cb = ds_simplify(ds_call2(SYM_Subtract,
                       ds_call2(SYM_Subtract, expr_copy(Eb),
                           ds_call2(SYM_Times, expr_new_integer(2), expr_copy(Ea))),
                       expr_new_integer(2)));                                 /* Eb - 2 Ea - 2 */
        int big = ds_is_zero(ca) ? 0 : (ds_is_zero(cb) ? 1 : -1);
        expr_free(ca); expr_free(cb);
        if (big >= 0) {
            Expr* Psm   = (big == 0) ? Eb : Ea;
            Expr* alpha = (big == 0) ? coefs[0] : coefs[1];
            Expr* beta  = (big == 0) ? coefs[1] : coefs[0];
            Expr* d = ds_simplify(ds_call2(SYM_Plus, expr_copy(Psm), expr_new_integer(2)));  /* d = c+1 = P_sm+2 */
            if (!ds_is_zero(alpha)) {
                Expr* sqrtna = ds_call1("Sqrt", ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(alpha)));
                Expr* invd = ds_call2(SYM_Power, expr_copy(d), expr_new_integer(-1));         /* 1/d */
                Expr* mu = ds_simplify(ds_call2(SYM_Times,                                    /* 1/(2 d) */
                               ds_call2(SYM_Power, expr_new_integer(2), expr_new_integer(-1)), expr_copy(invd)));
                Expr* kappa = ds_simplify(ds_call2(SYM_Times, expr_copy(beta),                /* beta/(2 d Sqrt[-alpha]) */
                                  ds_call2(SYM_Power,
                                      ds_call2(SYM_Times, expr_new_integer(2),
                                          ds_call2(SYM_Times, expr_copy(d), expr_copy(sqrtna))),
                                      expr_new_integer(-1))));
                Expr* z = ds_simplify(ds_call2(SYM_Times,                                     /* (2 Sqrt[-alpha]/d) x^d */
                              ds_call2(SYM_Times, expr_new_integer(2),
                                  ds_call2(SYM_Times, expr_copy(sqrtna), expr_copy(invd))),
                              ds_call2(SYM_Power, expr_copy(x), expr_copy(d))));
                Expr* pref = ds_simplify(ds_call2(SYM_Power, expr_copy(x),                    /* x^((1-d)/2) */
                                 ds_call2(SYM_Times,
                                     ds_call2(SYM_Subtract, expr_new_integer(1), expr_copy(d)),
                                     ds_call2(SYM_Power, expr_new_integer(2), expr_new_integer(-1)))));
                Expr* negmu = ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(mu));
                Expr* w0 = whittaker_M_1F1(kappa, mu, z);
                Expr* w1 = whittaker_M_1F1(kappa, negmu, z);
                general = combo(ds_call2(SYM_Times, expr_copy(pref), w0),
                                ds_call2(SYM_Times, expr_copy(pref), w1));
                expr_free(sqrtna); expr_free(invd); expr_free(mu); expr_free(kappa);
                expr_free(z); expr_free(pref); expr_free(negmu);
            }
            expr_free(d);
        }
        if (!general) {
            /* General Bessel row:  Q = A x^m + B x^(-2)
             *    ->  y = Sqrt[x] Z_nu(kappa x^((m+2)/2)),
             *        nu = Sqrt[1 - 4B]/(m+2),  kappa = 2 Sqrt[A]/(m+2),  Z = BesselJ/Y.
             * This strictly generalises the single-power row above (B == 0 gives the
             * same nu = 1/(m+2)) and reaches the inverse-power pair {-1, -2} that the
             * Whittaker condition P_big - 2 P_sm == 2 misses -- the shape that the
             * OperatorFactor adjoint peel hands back as its order-2 right factor
             * (2.1.2-294/296: z'' + (1/x + 1/x^2) z == 0 -> Sqrt[x] Z_{I Sqrt[3]}(2 Sqrt[x]),
             * verified to 1e-39).  A complex nu is fine: BesselJ/BesselY numericize
             * there, so the caller's sf_num_ok gate is genuine. */
            int im2 = -1;
            for (int gi = 0; gi < 2; gi++) {
                Expr* t = ds_simplify(ds_call2(SYM_Plus, expr_copy(exps[gi]), expr_new_integer(2)));
                bool is = ds_is_zero(t); expr_free(t);
                if (is) { im2 = gi; break; }
            }
            if (im2 >= 0) {
                Expr* B = coefs[im2], *m = exps[1 - im2], *A = coefs[1 - im2];
                Expr* mp2 = ds_simplify(ds_call2(SYM_Plus, expr_copy(m), expr_new_integer(2)));
                if (!ds_is_zero(mp2) && !ds_is_zero(A)) {
                    Expr* invp = ds_call2(SYM_Power, expr_copy(mp2), expr_new_integer(-1));
                    Expr* nu = ds_simplify(ds_call2(SYM_Times,                    /* Sqrt[1-4B]/(m+2) */
                                   ds_call1("Sqrt", ds_call2(SYM_Subtract, expr_new_integer(1),
                                       ds_call2(SYM_Times, expr_new_integer(4), expr_copy(B)))),
                                   expr_copy(invp)));
                    Expr* kappa = ds_simplify(ds_call2(SYM_Times,                 /* 2 Sqrt[A]/(m+2) */
                                      ds_call2(SYM_Times, expr_new_integer(2),
                                          ds_call1("Sqrt", expr_copy(A))),
                                      expr_copy(invp)));
                    Expr* halfp = ds_simplify(ds_call2(SYM_Times, expr_copy(mp2), /* (m+2)/2 */
                                      ds_call2(SYM_Power, expr_new_integer(2), expr_new_integer(-1))));
                    Expr* arg = ds_simplify(ds_call2(SYM_Times, expr_copy(kappa),
                                    ds_call2(SYM_Power, expr_copy(x), expr_copy(halfp))));
                    Expr* sqrtx = powrat(x, 1, 2);
                    general = combo(
                        ds_call2(SYM_Times, expr_copy(sqrtx),
                            expr_new_function(expr_new_symbol("BesselJ"),
                                (Expr*[]){ expr_copy(nu), expr_copy(arg) }, 2)),
                        ds_call2(SYM_Times, expr_copy(sqrtx),
                            expr_new_function(expr_new_symbol("BesselY"),
                                (Expr*[]){ expr_copy(nu), expr_copy(arg) }, 2)));
                    expr_free(sqrtx); expr_free(arg); expr_free(halfp);
                    expr_free(kappa); expr_free(nu); expr_free(invp);
                }
                expr_free(mp2);
            }
        }
    }

    for (int g = 0; g < ng; g++) { expr_free(exps[g]); expr_free(coefs[g]); }
    expr_free(x); expr_free(Qe);
    return general;
}

/* Numeric self-verify of a Whittaker candidate `base` against the REDUCED equation
 * w'' + Qc w == 0 (Qc = the reduced-form coefficient given to the recogniser).  This
 * is verified BEFORE the recovery factor is applied — crucial, because the composed
 * mu*base stacks the recovery's and Whittaker's same-base symbolic-radical powers,
 * whose re-summation drives the evaluator into $IterationLimit on some multi-parameter
 * confluent forms; the reduced residual (single power per branch, no recovery) never
 * stacks, so it is both cheap and thrash-free.  base solving the reduced equation
 * means mu*base solves the original (the normal-form transformation is exact), so the
 * caller may emit mu*base without the original-equation gate.  base, Qc borrowed. */
static bool whittaker_num_ok(const Expr* base, const Expr* Qc, const char* xv) {
    Expr* b2 = ds_d(ds_d(expr_copy((Expr*)base), expr_new_symbol(xv)), expr_new_symbol(xv));
    Expr* R  = ds_call2(SYM_Plus, b2,
                   ds_call2(SYM_Times, expr_copy((Expr*)Qc), expr_copy((Expr*)base)));  /* w'' + Qc w */
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
    const double xs[] = { 0.2, 0.4, 0.6, 0.8, 1.0 };
    int small = 0, big = 0;
    for (int i = 0; i < 5; i++) {
        double m = sf_abs_at(R, xv, xs[i]);
        if (isnan(m) || !isfinite(m)) continue;
        if (m < 1e-6) small++; else if (m > 1e-3) big++;
    }
    expr_free(R);
    return small >= 2 && big == 0;
}

/* ---- Whittaker / confluent-hypergeometric (1F1) reduced form ----
 *   The reduced equation  w'' + Qc w == 0  has a SINGLE finite regular singular
 *   point x0 (Qc a double pole there) and a rank-1 irregular point at infinity
 *   (Qc -> b2 != 0 as x -> Infinity).  Write Qc = b2 + b1/(x-x0) + b0/(x-x0)^2
 *   and match Whittaker's normal form
 *       W'' + (-1/4 + kappa/z + (1/4 - mu^2)/z^2) W == 0     under  z = c (x-x0):
 *       c = 2 Sqrt[-b2],  mu = Sqrt[1/4 - b0],  kappa = b1/c.
 *   The two independent solutions are emitted as the VERIFIABLE 1F1 form
 *       WhittakerM[kappa, +-mu, z] = Exp[-z/2] z^(1/2+-mu) 1F1[1/2+-mu-kappa, 1+-2mu, z]
 *   -- Hypergeometric1F1 numericizes and carries a z-derivative rule, so the
 *   residual back-substitutes; the inert WhittakerM/W heads are never emitted.
 *   Declines when 2 mu is a numeric integer: the two 1F1 partners are then
 *   dependent (2mu==0) or a lower parameter 1-+2mu is a non-positive integer
 *   (singular 1F1), so the caller falls through to the Frobenius log-second
 *   solution.  The single finite double pole is isolated by the squarefree-part
 *   trick (den / gcd(den,den') of degree 1), so the two-finite-pole Gauss case
 *   (affine->2F1) and the pole-free Airy/polynomial cases decline here.  The
 *   pole-at-0, b1==0 special case (kappa==0) is the modified-Bessel normal form
 *   claimed by the A + B/x^2 row above; this runs after it under the caller's
 *   !general guard, and every emission is gated by the caller's sf_num_ok.
 *   Qc, xvar borrowed; returns C[1] W1 + C[2] W2, or NULL. */
static Expr* specialform_whittaker_basis(const Expr* Qc, const char* xvar, bool degen_ok) {
    Expr* rt = sf_ct(expr_copy((Expr*)Qc));                 /* Cancel[Together[Qc]] */
    if (ds_is_zero(rt)) { expr_free(rt); return NULL; }
    Expr* den  = eval_and_free(ds_call1("Denominator", expr_copy(rt)));
    Expr* dden = ds_d(expr_copy(den), expr_new_symbol(xvar));
    Expr* gg   = eval_and_free(ds_call2("PolynomialGCD", expr_copy(den), dden));   /* consumes dden */
    Expr* rad  = sf_ct(ds_call2(SYM_Times, expr_copy(den),                          /* den/gcd(den,den') */
                     expr_new_function(expr_new_symbol(SYM_Power),
                         (Expr*[]){ expr_copy(gg), expr_new_integer(-1) }, 2)));
    Expr* dDen = eval_and_free(expr_new_function(expr_new_symbol("Exponent"),
                     (Expr*[]){ expr_copy(den), expr_new_symbol(xvar) }, 2));
    Expr* dRad = eval_and_free(expr_new_function(expr_new_symbol("Exponent"),
                     (Expr*[]){ expr_copy(rad), expr_new_symbol(xvar) }, 2));
    bool shape = (dDen->type == EXPR_INTEGER && dDen->data.integer == 2 &&
                  dRad->type == EXPR_INTEGER && dRad->data.integer == 1);
    expr_free(dDen); expr_free(dRad);
    Expr* general = NULL;
    if (shape) {
        /* x0 = -c0/c1 from the squarefree denominator rad = c1 x + c0 */
        Expr* r0 = eval_and_free(expr_new_function(expr_new_symbol("Coefficient"),
                       (Expr*[]){ expr_copy(rad), expr_new_symbol(xvar), expr_new_integer(0) }, 3));
        Expr* r1 = eval_and_free(expr_new_function(expr_new_symbol("Coefficient"),
                       (Expr*[]){ expr_copy(rad), expr_new_symbol(xvar), expr_new_integer(1) }, 3));
        Expr* x0 = sf_ct(ds_call2(SYM_Times, expr_new_integer(-1),
                       ds_call2(SYM_Times, expr_copy(r0),
                           expr_new_function(expr_new_symbol(SYM_Power),
                               (Expr*[]){ expr_copy(r1), expr_new_integer(-1) }, 2))));
        expr_free(r0); expr_free(r1);
        /* g = Cancel[Together[ Qc (x-x0)^2 ]] must be a degree-2 polynomial */
        Expr* g = sf_ct(ds_call2(SYM_Times, expr_copy((Expr*)Qc),
                      ds_call2(SYM_Power,
                          ds_call2(SYM_Subtract, expr_new_symbol(xvar), expr_copy(x0)),
                          expr_new_integer(2))));
        Expr* pq = eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ),
                       (Expr*[]){ expr_copy(g), expr_new_symbol(xvar) }, 2));
        Expr* dG = eval_and_free(expr_new_function(expr_new_symbol("Exponent"),
                       (Expr*[]){ expr_copy(g), expr_new_symbol(xvar) }, 2));
        bool gok = (pq->type == EXPR_SYMBOL && pq->data.symbol.name == SYM_True &&
                    dG->type == EXPR_INTEGER && dG->data.integer == 2);
        expr_free(pq); expr_free(dG);
        if (gok) {
            Expr* b0 = sf_ct(ds_subst(expr_copy(g), expr_new_symbol(xvar), expr_copy(x0)));
            Expr* b1 = sf_ct(ds_subst(ds_d(expr_copy(g), expr_new_symbol(xvar)),
                           expr_new_symbol(xvar), expr_copy(x0)));
            Expr* b2 = eval_and_free(expr_new_function(expr_new_symbol("Coefficient"),
                           (Expr*[]){ expr_copy(g), expr_new_symbol(xvar), expr_new_integer(2) }, 3));
            if (!ds_is_zero(b2)) {
                Expr* quarter = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                                    (Expr*[]){ expr_new_integer(4), expr_new_integer(-1) }, 2));    /* 1/4 */
                Expr* mu = sf_ct(ds_call1("Sqrt", ds_call2(SYM_Subtract, quarter, expr_copy(b0))));
                Expr* c  = sf_ct(ds_call2(SYM_Times, expr_new_integer(2),
                               ds_call1("Sqrt", ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(b2)))));
                Expr* kappa = sf_ct(ds_call2(SYM_Times, expr_copy(b1),
                                  expr_new_function(expr_new_symbol(SYM_Power),
                                      (Expr*[]){ expr_copy(c), expr_new_integer(-1) }, 2)));
                /* decline 2 mu in Z (dependent / singular 1F1 lower parameter) */
                Expr* twomu = sf_ct(ds_call2(SYM_Times, expr_new_integer(2), expr_copy(mu)));
                Expr* iq = eval_and_free(ds_call1("IntegerQ", twomu));   /* consumes twomu */
                bool dep = (iq->type == EXPR_SYMBOL && iq->data.symbol.name == SYM_True);
                expr_free(iq);
                /* Degenerate 2 mu in Z (late SpecialFunctionTransform only, degen_ok):
                 * the +mu 1F1 solution stays, and Tricomi's U[1/2+mu-kappa, 1+2mu, z]
                 * (a solution of the same Kummer equation for every lower parameter)
                 * replaces the dependent/singular -mu partner -- independent unless
                 * the upper parameter 1/2+mu-kappa is a non-positive integer (the 1F1
                 * is then a polynomial onto which U collapses), which still declines. */
                bool use_u = false;
                if (dep && degen_ok) {
                    Expr* up = sf_ct(ds_call2(SYM_Subtract,
                                   ds_call2(SYM_Plus, eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                                       (Expr*[]){ expr_new_integer(2), expr_new_integer(-1) }, 2)), expr_copy(mu)),
                                   expr_copy(kappa)));
                    Expr* iq2 = eval_and_free(ds_call1("IntegerQ", expr_copy(up)));
                    bool upint = (iq2->type == EXPR_SYMBOL && iq2->data.symbol.name == SYM_True);
                    expr_free(iq2);
                    bool polyF = upint && up->type == EXPR_INTEGER && up->data.integer <= 0;
                    expr_free(up);
                    if (!polyF) use_u = true;
                }
                if (!dep || use_u) {
                    Expr* z = sf_ct(ds_call2(SYM_Times, expr_copy(c),
                                  ds_call2(SYM_Subtract, expr_new_symbol(xvar), expr_copy(x0))));
                    Expr* ehalf = eval_and_free(ds_call1("Exp",
                                      ds_call2(SYM_Times,
                                          eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                                              (Expr*[]){ expr_new_integer(-2), expr_new_integer(-1) }, 2)), /* -1/2 */
                                          expr_copy(z))));
                    Expr* half = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                                     (Expr*[]){ expr_new_integer(2), expr_new_integer(-1) }, 2));           /* 1/2 */
                    Expr* branches[2] = { NULL, NULL };
                    for (int s = 0; s < 2; s++) {
                        Expr* smu = (s == 0 || use_u) ? expr_copy(mu)
                                  : ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(mu));   /* +-mu */
                        Expr* exp1  = ds_call2(SYM_Plus, expr_copy(half), expr_copy(smu));        /* 1/2 +- mu */
                        Expr* upper = ds_call2(SYM_Subtract,
                                          ds_call2(SYM_Plus, expr_copy(half), expr_copy(smu)),
                                          expr_copy(kappa));                                       /* 1/2 +- mu - kappa */
                        Expr* lower = ds_call2(SYM_Plus, expr_new_integer(1),
                                          ds_call2(SYM_Times, expr_new_integer(2), expr_copy(smu))); /* 1 +- 2 mu */
                        expr_free(smu);
                        Expr* zpow = eval_and_free(ds_call2(SYM_Power, expr_copy(z), exp1));       /* z^(1/2+-mu) */
                        Expr* onef = expr_new_function(expr_new_symbol(
                                         (use_u && s == 1) ? "HypergeometricU" : "Hypergeometric1F1"),
                                         (Expr*[]){ upper, lower, expr_copy(z) }, 3);
                        branches[s] = eval_and_free(ds_call2(SYM_Times, expr_copy(ehalf),
                                          ds_call2(SYM_Times, zpow, onef)));
                    }
                    general = combo(branches[0], branches[1]);   /* consumes branches[] */
                    expr_free(z); expr_free(ehalf); expr_free(half);
                    /* Self-verify against the reduced equation w'' + Qc w == 0 before
                     * returning (Qc is this recogniser's argument = the reduced-form
                     * coefficient).  Guards a mis-read pole / parameter, and is safe on
                     * the thrash-prone forms because it never applies the recovery. */
                    if (general && !whittaker_num_ok(general, Qc, xvar)) {
                        expr_free(general);
                        general = NULL;
                    }
                }
                expr_free(mu); expr_free(c); expr_free(kappa);
            }
            expr_free(b0); expr_free(b1); expr_free(b2);
        }
        expr_free(x0); expr_free(g);
    }
    expr_free(rt); expr_free(den); expr_free(gg); expr_free(rad);
    return general;
}

/* ---- Gauss (hypergeometric 2F1) recogniser on the CANONICAL interval
 *      W = x(1-x): Q == -a b/W, P == (c-(a+b+1)x)/W.  Factored so the affine
 *      row below can call it on the mapped (P~, Q~, s).  Pc, Qc, xvar borrowed;
 *      returns C[1] 2F1[...] + C[2] x^(1-c) 2F1[...] or NULL. */
static Expr* specialform_gauss_basis(Expr* Pc, Expr* Qc, const char* xvar) {
    Expr* general = NULL;
    Expr* W = eval_and_free(ds_call2(SYM_Times, expr_new_symbol(xvar),
                  ds_call2(SYM_Subtract, expr_new_integer(1), expr_new_symbol(xvar))));  /* x(1-x) */
    Expr* negprod = sf_ct(ds_call2(SYM_Times, expr_copy(W), expr_copy(Qc)));             /* -a b */
    if (ds_free_of(negprod, xvar)) {
        Expr* L = sf_ct(ds_call2(SYM_Times, expr_copy(W), expr_copy(Pc)));                /* c-(a+b+1)x */
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
    return general;
}

/* One F-homotopy attempt on the canonical interval: Y = sv^e0 (1-sv)^e1 F pulls
 * the local exponents e0 at sv=0 and e1 at sv=1 out of Y'' + Pt Y' + Qt Y == 0, so
 * F satisfies the reduced equation  F'' + (Pt + 2 w'/w) F' + (Qt + Pt w'/w + w''/w)
 * F == 0  with w = sv^e0 (1-sv)^e1, w'/w = e0/sv - e1/(1-sv).  If F is in canonical
 * Gauss form, returns sv^e0 (1-sv)^e1 * (Gauss 2F1 combo); else NULL.  This is the
 * reduction of any two-finite-RSP Fuchsian equation to Hypergeometric2F1 — the
 * exponent shift the plain affine map lacks (Gegenbauer/Jacobi/associated Legendre
 * carry nonzero exponents at the finite singular points).  All args borrowed. */
static Expr* specialform_fhomotopy(Expr* Pt, Expr* Qt, const char* sv,
                                   Expr* e0, Expr* e1) {
    Expr* oms = ds_call2(SYM_Subtract, expr_new_integer(1), expr_new_symbol(sv));   /* 1-s */
    Expr* wl  = sf_ct(ds_call2(SYM_Subtract,                                         /* w'/w */
                    ds_call2(SYM_Times, expr_copy(e0),
                        ds_call2(SYM_Power, expr_new_symbol(sv), expr_new_integer(-1))),
                    ds_call2(SYM_Times, expr_copy(e1),
                        ds_call2(SYM_Power, expr_copy(oms), expr_new_integer(-1)))));
    Expr* wll = sf_ct(ds_call2(SYM_Plus,                                             /* w''/w */
                    ds_d(expr_copy(wl), expr_new_symbol(sv)),
                    ds_call2(SYM_Power, expr_copy(wl), expr_new_integer(2))));
    Expr* PtF = sf_ct(ds_call2(SYM_Plus, expr_copy(Pt),
                    ds_call2(SYM_Times, expr_new_integer(2), expr_copy(wl))));
    Expr* QtF = sf_ct(ds_call2(SYM_Plus, expr_copy(Qt),
                    ds_call2(SYM_Plus, ds_call2(SYM_Times, expr_copy(Pt), expr_copy(wl)),
                        expr_copy(wll))));
    Expr* basisF = specialform_gauss_basis(PtF, QtF, sv);
    Expr* out = NULL;
    if (basisF) {
        Expr* pref = eval_and_free(ds_call2(SYM_Times,                               /* sv^e0 (1-sv)^e1 */
                         ds_call2(SYM_Power, expr_new_symbol(sv), expr_copy(e0)),
                         ds_call2(SYM_Power, expr_copy(oms), expr_copy(e1))));
        out = eval_and_free(ds_call2(SYM_Times, pref, basisF));   /* consumes pref, basisF */
    }
    expr_free(oms); expr_free(wl); expr_free(wll); expr_free(PtF); expr_free(QtF);
    return out;
}

/* Numeric self-verify for the affine->Gauss row.  A wrong RSP map yields a 2F1
 * whose residual zero_test cannot prove nonzero (so the substrate would KEEP it),
 * so this gate is the sole guard of the 0-FAIL invariant for the affine row.
 * Mirrors sf_num_ok, but samples x INSIDE the mapped segment [x1, x2] (so the
 * separation s=(x-x1)/h lands in (0,1) and the x^(1-c) factor stays real/finite);
 * x1, h are instantiated at the same generic reals as the residual's params. */
static bool hgc_num_ok(const DSolveProblem* P, const Expr* cand,
                       const char* xv, const char* yname,
                       const Expr* x1, const Expr* h) {
    Expr* R = expr_copy(P->eq_residuals[0]);
    Expr* b0 = expr_copy((Expr*)cand);
    Expr* b1 = ds_d(expr_copy((Expr*)cand), expr_new_symbol(xv));
    Expr* b2 = ds_d(ds_d(expr_copy((Expr*)cand), expr_new_symbol(xv)), expr_new_symbol(xv));
    R = ds_subst(R, ds_make_funcapp(yname, 2, xv), b2);
    R = ds_subst(R, ds_make_funcapp(yname, 1, xv), b1);
    R = ds_subst(R, ds_make_funcapp(yname, 0, xv), b0);
    R = ds_subst(R, ds_const(1), expr_new_real(1.3));
    R = ds_subst(R, ds_const(2), expr_new_real(0.7));
    Expr* x1c = expr_copy((Expr*)x1);
    Expr* hc  = expr_copy((Expr*)h);
    {
        const char* skip[] = { xv, intern_symbol("E"), intern_symbol("Pi"),
            intern_symbol("I"), intern_symbol("C"), intern_symbol("EulerGamma"),
            intern_symbol("Degree"), intern_symbol("GoldenRatio"),
            intern_symbol("Catalan"), intern_symbol("Infinity") };
        const int nskip = (int)(sizeof(skip)/sizeof(skip[0]));
        const char* syms[64]; int ns = 0;
        sf_collect_params(R, syms, &ns, 64);
        sf_collect_params(x1c, syms, &ns, 64);
        sf_collect_params(hc, syms, &ns, 64);
        int pidx = 0;
        for (int i = 0; i < ns; i++) {
            bool sk = false;
            for (int j = 0; j < nskip; j++) if (syms[i] == skip[j]) { sk = true; break; }
            if (sk) continue;
            double v = 0.29 + 0.13 * (double)pidx; pidx++;
            R   = ds_subst(R,   expr_new_symbol(syms[i]), expr_new_real(v));
            x1c = ds_subst(x1c, expr_new_symbol(syms[i]), expr_new_real(v));
            hc  = ds_subst(hc,  expr_new_symbol(syms[i]), expr_new_real(v));
        }
    }
    Expr* x1n = eval_and_free(ds_call1("N", x1c));
    Expr* hn  = eval_and_free(ds_call1("N", hc));
    double x1v = (x1n->type == EXPR_REAL) ? x1n->data.real
               : (x1n->type == EXPR_INTEGER) ? (double)x1n->data.integer : NAN;
    double hv  = (hn->type == EXPR_REAL) ? hn->data.real
               : (hn->type == EXPR_INTEGER) ? (double)hn->data.integer : NAN;
    expr_free(x1n); expr_free(hn);
    bool ok = false;
    if (isfinite(x1v) && isfinite(hv) && hv != 0.0) {
        const double frac[] = { 0.15, 0.3, 0.45, 0.6, 0.85 };   /* s in (0,1) */
        int small = 0, big = 0;
        for (int i = 0; i < 5; i++) {
            double m = sf_abs_at(R, xv, x1v + hv * frac[i]);
            if (isnan(m) || !isfinite(m)) continue;
            if (m < 1e-6) small++; else if (m > 1e-3) big++;
        }
        ok = (small >= 2 && big == 0);
    }
    expr_free(R);
    return ok;
}

/* Kummer row of dsolve_specialform_try (see the comment at its call): the
 * C[1], C[2] basis of y'' + Pc y' + Qc y == 0, or NULL. */
static Expr* specialform_kummer_basis(Expr* Pc, Expr* Qc, const char* xvar, bool sym_b) {
    Expr* out = NULL;
    Expr* xP = ds_simplify(ds_call2(SYM_Times, expr_new_symbol(xvar), expr_copy(Pc)));
    Expr* slope = ds_d(expr_copy(xP), expr_new_symbol(xvar));
    if (ds_free_of(slope, xvar) && !ds_is_zero(slope)) {
        Expr* ks = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(slope)));  /* s */
        Expr* kb = ds_simplify(ds_subst(expr_copy(xP), expr_new_symbol(xvar), expr_new_integer(0)));  /* b */
        Expr* xQ = ds_simplify(ds_call2(SYM_Times, expr_new_symbol(xvar), expr_copy(Qc)));
        bool bnum = specialform_is_number(kb);
        if (ds_free_of(kb, xvar) && ds_free_of(xQ, xvar) && !bnum && sym_b) {
            /* SYMBOLIC b (late SpecialFunctionTransform only): {1F1[a, b, z], U[a, b, z]}
             * -- no symbolic power z^(1-b), so the residual stays a pure pFq sum; U is
             * independent of 1F1 unless a is a non-positive integer (1F1 a polynomial
             * onto which U collapses), which declines.  The caller gates numerically. */
            Expr* ka = sf_ct(ds_call2(SYM_Times, expr_new_integer(-1),
                           ds_call2(SYM_Times, expr_copy(xQ), powneg1_sf(expr_copy(ks)))));
            bool apoly = ka->type == EXPR_INTEGER && ka->data.integer <= 0;
            if (!apoly) {
                Expr* z = eval_and_free(ds_call2(SYM_Times, expr_copy(ks), expr_new_symbol(xvar)));
                Expr* M = eval_and_free(expr_new_function(expr_new_symbol(SYM_Hypergeometric1F1),
                              (Expr*[]){ expr_copy(ka), expr_copy(kb), expr_copy(z) }, 3));
                Expr* U = eval_and_free(expr_new_function(expr_new_symbol(SYM_HypergeometricU),
                              (Expr*[]){ expr_copy(ka), expr_copy(kb), expr_copy(z) }, 3));
                out = combo(M, U);
                expr_free(z);
            }
            expr_free(ka);
        } else if (ds_free_of(kb, xvar) && ds_free_of(xQ, xvar) && bnum) {
            Expr* ka = ds_simplify(ds_call2(SYM_Times, expr_new_integer(-1),
                           ds_call2(SYM_Times, expr_copy(xQ), powneg1_sf(expr_copy(ks)))));  /* a */
            Expr* iq = eval_and_free(ds_call1("IntegerQ", expr_copy(kb)));
            bool bint = (iq->type == EXPR_SYMBOL && iq->data.symbol.name == SYM_True);
            expr_free(iq);
            bool bpos = bint && kb->type == EXPR_INTEGER && kb->data.integer >= 1;
            Expr* z = eval_and_free(ds_call2(SYM_Times, expr_copy(ks), expr_new_symbol(xvar)));
            /* M = 1F1[a, b, z];  Mbar = z^(1-b) 1F1[a-b+1, 2-b, z];  U = U[a, b, z] */
            Expr* M = NULL;
            if (!bint || bpos)
                M = eval_and_free(expr_new_function(expr_new_symbol(SYM_Hypergeometric1F1),
                        (Expr*[]){ expr_copy(ka), expr_copy(kb), expr_copy(z) }, 3));
            Expr* Mbar = NULL;
            if (!bpos) {
                Expr* zpow = eval_and_free(ds_call2(SYM_Power, expr_copy(z),
                                 ds_call2(SYM_Subtract, expr_new_integer(1), expr_copy(kb))));
                Expr* a2 = eval_and_free(ds_call2(SYM_Plus,
                               ds_call2(SYM_Subtract, expr_copy(ka), expr_copy(kb)), expr_new_integer(1)));
                Expr* b2 = eval_and_free(ds_call2(SYM_Subtract, expr_new_integer(2), expr_copy(kb)));
                Expr* h2 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Hypergeometric1F1),
                               (Expr*[]){ a2, b2, expr_copy(z) }, 3));
                Mbar = eval_and_free(ds_call2(SYM_Times, zpow, h2));
            }
            /* Degenerate integer case: U[a, b, z] collapses onto the 1F1
             * solution when that solution is a polynomial (times z^(1-b)):
             *   b >= 1, a == -n:            U = (-1)^n (b)_n 1F1[-n, b, z]
             *   b <= 0, a - b + 1 == -n:    U = z^(1-b) U[-n, 2-b, z] ~ Mbar
             * (Wronskian of {1F1, U} carries 1/Gamma(a) resp. 1/Gamma(a-b+1)).
             * Both reduce to the Kummer equation in (-n, B) with B = b resp.
             * 2-b >= 1 -- the second via w = z^(b-1) y -- whose logarithmic
             * partner specialform_kummer_log_partner supplies. */
            int64_t dn = -1, dB = 0;
            if (bint && kb->type == EXPR_INTEGER && ka->type == EXPR_INTEGER) {
                int64_t ia = ka->data.integer, ib = kb->data.integer;
                if (!sf_small_int(ia) || !sf_small_int(ib)) ;   /* decline below */
                else if (bpos && ia <= 0)              { dn = -ia;         dB = ib; }
                else if (!bpos && ia - ib + 1 <= 0) { dn = -(ia - ib + 1); dB = 2 - ib; }
            }
            if (dn >= 0) {
                Expr* w2 = specialform_kummer_log_partner(dn, dB, z);
                if (w2 && !bpos)
                    w2 = eval_and_free(ds_call2(SYM_Times, eval_and_free(ds_call2(SYM_Power,
                             expr_copy(z), ds_call2(SYM_Subtract, expr_new_integer(1),
                             expr_copy(kb)))), w2));
                if (w2) out = combo(bpos ? M : Mbar, w2);
                else { expr_free(M); expr_free(Mbar); }   /* too large: decline */
            } else if (!bint) out = combo(M, Mbar);
            else {
                Expr* U = eval_and_free(expr_new_function(expr_new_symbol(SYM_HypergeometricU),
                              (Expr*[]){ expr_copy(ka), expr_copy(kb), expr_copy(z) }, 3));
                out = combo(bpos ? M : Mbar, U);
            }
            expr_free(z); expr_free(ka);
        }
        expr_free(ks); expr_free(kb); expr_free(xQ);
    }
    expr_free(slope); expr_free(xP);
    return out;
}

/* ---- affine -> Gauss 2F1 (see the call site in dsolve_specialform_try).
 *      Pc, Qc borrowed; returns the verified general solution or NULL. ---- */
static Expr* specialform_affine_gauss(DSolveProblem* P, Expr* Pc, Expr* Qc,
                                      const char* xvar, bool both_orient) {
    Expr* general = NULL;
    Expr* dP  = eval_and_free(ds_call1("Denominator", ds_call1("Together", expr_copy(Pc))));
    Expr* dQ  = eval_and_free(ds_call1("Denominator", ds_call1("Together", expr_copy(Qc))));
    Expr* u   = eval_and_free(ds_call2("PolynomialLCM", dP, dQ));   /* consumes dP, dQ */
    Expr* du  = ds_d(expr_copy(u), expr_new_symbol(xvar));
    Expr* g   = eval_and_free(ds_call2("PolynomialGCD", expr_copy(u), du));  /* consumes du */
    Expr* rad = eval_and_free(ds_call1("Cancel", ds_call2(SYM_Times, expr_copy(u),
                    expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ expr_copy(g), expr_new_integer(-1) }, 2))));  /* u / gcd(u,u') */
    Expr* pq = eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ),
                   (Expr*[]){ expr_copy(rad), expr_new_symbol(xvar) }, 2));
    bool ratl = (pq->type == EXPR_SYMBOL && pq->data.symbol.name == SYM_True);
    expr_free(pq);
    bool deg2 = false;
    if (ratl) {
        Expr* deg = eval_and_free(expr_new_function(expr_new_symbol("Exponent"),
                        (Expr*[]){ expr_copy(rad), expr_new_symbol(xvar) }, 2));
        deg2 = (deg->type == EXPR_INTEGER && deg->data.integer == 2);
        expr_free(deg);
    }
    if (deg2) {
        #define AF_COEF(k) eval_and_free(expr_new_function(expr_new_symbol("Coefficient"), \
            (Expr*[]){ expr_copy(rad), expr_new_symbol(xvar), expr_new_integer(k) }, 3))
        Expr* c2 = AF_COEF(2); Expr* c1 = AF_COEF(1); Expr* c0 = AF_COEF(0);
        #undef AF_COEF
        /* roots of rad == 0, i.e. of t^2 - S t + Pr, S = -c1/c2, Pr = c0/c2 */
        Expr* S = ds_simplify(ds_call2(SYM_Times, expr_new_integer(-1),
                      ds_call2(SYM_Times, expr_copy(c1),
                          expr_new_function(expr_new_symbol(SYM_Power),
                              (Expr*[]){ expr_copy(c2), expr_new_integer(-1) }, 2))));
        Expr* Pr = ds_simplify(ds_call2(SYM_Times, expr_copy(c0),
                       expr_new_function(expr_new_symbol(SYM_Power),
                           (Expr*[]){ expr_copy(c2), expr_new_integer(-1) }, 2)));
        Expr* x1 = NULL; Expr* x2 = NULL;
        if (specialform_quad_roots(S, Pr, &x1, &x2) &&
            !sf_has_radical(x1) && !sf_has_radical(x2)) {   /* rational RSPs only */
          /* both_orient (late slot): when the map with x1 at s = 0 declines -- e.g.
           * the exponents there are {0, 0}, an integer c the Gauss row refuses --
           * retry with the singular points swapped, so the point with a non-integer
           * exponent difference is the expansion point (4 t^2 (1-t) Y'' +
           * t (2 - 6 t) Y' - a (a+1) Y == 0). */
          for (int orient = 0; orient < (both_orient ? 2 : 1) && !general; orient++) {
            if (orient == 1) { Expr* tx = x1; x1 = x2; x2 = tx; }
            Expr* h = ds_simplify(ds_call2(SYM_Subtract, expr_copy(x2), expr_copy(x1)));
            if (!ds_is_zero(h)) {
                const char* sv = intern_symbol("DSolve`hgs");
                Expr* xsub = eval_and_free(ds_call2(SYM_Plus, expr_copy(x1),
                                 ds_call2(SYM_Times, expr_copy(h), expr_new_symbol(sv))));
                /* Cancel[Together[.]] (not Simplify) for the mapped coefficients:
                 * Simplify can mis-reduce a constant-over-quadratic like
                 * 3/(4 s (s-1)) to -3/(4 s), dropping the second pole so the
                 * mapped equation is no longer canonical Gauss. */
                Expr* Pt = eval_and_free(ds_call1("Cancel", ds_call1("Together",
                               ds_call2(SYM_Times, expr_copy(h),
                                   ds_subst(expr_copy(Pc), expr_new_symbol(xvar), expr_copy(xsub))))));
                Expr* Qt = eval_and_free(ds_call1("Cancel", ds_call1("Together",
                               ds_call2(SYM_Times,
                                   ds_call2(SYM_Power, expr_copy(h), expr_new_integer(2)),
                                   ds_subst(expr_copy(Qc), expr_new_symbol(xvar), expr_copy(xsub))))));
                expr_free(xsub);
                /* local exponents (indicial roots) at s=0 and s=1: rho^2 -
                 * (1-p_i) rho + q_i == 0, p0=(s Pt)|0, q0=(s^2 Qt)|0, and the
                 * s->1 analogue.  Cancel[Together[.]] (not Simplify) reduces the
                 * pole factor to lowest terms so the s^k / s^k cancels BEFORE the
                 * s->0 substitution — else s->0 hits 0*Infinity = Indeterminate. */
                #define SF_POLE(ARG, AT) ds_subst( \
                    eval_and_free(ds_call1("Cancel", ds_call1("Together", (ARG)))), \
                    expr_new_symbol(sv), expr_new_integer(AT))
                Expr* p0 = SF_POLE(ds_call2(SYM_Times, expr_new_symbol(sv), expr_copy(Pt)), 0);
                Expr* q0 = SF_POLE(ds_call2(SYM_Times,
                               ds_call2(SYM_Power, expr_new_symbol(sv), expr_new_integer(2)), expr_copy(Qt)), 0);
                Expr* smo = ds_call2(SYM_Subtract, expr_new_symbol(sv), expr_new_integer(1));   /* s-1 */
                Expr* p1 = SF_POLE(ds_call2(SYM_Times, expr_copy(smo), expr_copy(Pt)), 1);
                Expr* q1 = SF_POLE(ds_call2(SYM_Times,
                               ds_call2(SYM_Power, expr_copy(smo), expr_new_integer(2)), expr_copy(Qt)), 1);
                #undef SF_POLE
                expr_free(smo);
                Expr* e0[2] = { NULL, NULL }; Expr* e1[2] = { NULL, NULL };
                bool have0 = false, have1 = false;
                if (ds_free_of(p0, sv) && ds_free_of(q0, sv)) {
                    Expr* S0 = ds_simplify(ds_call2(SYM_Subtract, expr_new_integer(1), expr_copy(p0)));
                    have0 = specialform_quad_roots(S0, q0, &e0[0], &e0[1]);
                    expr_free(S0);
                }
                if (ds_free_of(p1, sv) && ds_free_of(q1, sv)) {
                    Expr* S1 = ds_simplify(ds_call2(SYM_Subtract, expr_new_integer(1), expr_copy(p1)));
                    have1 = specialform_quad_roots(S1, q1, &e1[0], &e1[1]);
                    expr_free(S1);
                }
                if (have0 && have1) {
                    Expr* sback = eval_and_free(ds_call2(SYM_Times,
                                      ds_call2(SYM_Subtract, expr_new_symbol(xvar), expr_copy(x1)),
                                      expr_new_function(expr_new_symbol(SYM_Power),
                                          (Expr*[]){ expr_copy(h), expr_new_integer(-1) }, 2)));  /* (x-x1)/h */
                    for (int i = 0; i < 2 && !general; i++) {
                        for (int j = 0; j < 2 && !general; j++) {
                            Expr* cand_s = specialform_fhomotopy(Pt, Qt, sv, e0[i], e1[j]);
                            if (!cand_s) continue;
                            /* consumes cand_s */
                            Expr* cand = ds_subst(cand_s, expr_new_symbol(sv), expr_copy(sback));
                            if (hgc_num_ok(P, cand, xvar, P->fun_names[0], x1, h)) general = cand;
                            else expr_free(cand);
                        }
                    }
                    expr_free(sback);
                }
                for (int k = 0; k < 2; k++) {
                    if (e0[k]) expr_free(e0[k]);
                    if (e1[k]) expr_free(e1[k]);
                }
                expr_free(p0); expr_free(q0); expr_free(p1); expr_free(q1);
                expr_free(Pt); expr_free(Qt);
            }
            expr_free(h);
          }
        }
        expr_free(S); expr_free(Pr); expr_free(c2); expr_free(c1); expr_free(c0);
        if (x1) expr_free(x1);
        if (x2) expr_free(x2);
    }
    expr_free(u); expr_free(g); expr_free(rad);
    return general;
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
                /* Ordinary Legendre (mu == 0) -> the 2-arg LegendreP/Q, which
                 * numericize for a symbolic-then-instantiated degree and so verify.
                 * The ASSOCIATED case (mu != 0) is declined HERE so it falls through
                 * to the affine -> Gauss 2F1 row below: the 3-arg LegendreP[nu,mu,x]
                 * numericizes only for integer nu, mu, so its residual verify and the
                 * corpus numeric check both fail, whereas the equivalent
                 * Hypergeometric2F1 form is verifiable at symbolic parameters. */
                if (ds_is_zero(musq)) {
                    Expr* b0 = expr_new_function(expr_new_symbol(SYM_LegendreP),
                                   (Expr*[]){ expr_copy(nu), expr_new_symbol(xvar) }, 2);
                    Expr* b1 = expr_new_function(expr_new_symbol(SYM_LegendreQ),
                                   (Expr*[]){ expr_copy(nu), expr_new_symbol(xvar) }, 2);
                    general = combo(b0, b1);
                }
                expr_free(B); expr_free(nu); expr_free(musq);
            }
            expr_free(qq); expr_free(c0); expr_free(c2); expr_free(recon);
        }
        expr_free(omx2); expr_free(pchk);
    }

    /* ---- P == 0 reduced-form recognisers (Airy + the three Bessel-reducible
     *      potentials), factored into specialform_reduced_basis so the Liouville
     *      normal-form pre-pass below can reuse them on -r when a y' term is
     *      present.  Behaviour on P == 0 equations is unchanged. ---- */
    if (!general && ds_is_zero(Pc)) general = specialform_reduced_basis(Qc, xvar);

    /* ---- Generalised power-potential (symbolic-exponent Bessel + two-term
     *      Coulomb/Whittaker), P == 0.  Numerically gated: the symbolic-parameter
     *      residual is undecidable by zero_test, so a sign/branch error would
     *      otherwise reach the corpus as a numericizing FAIL. ---- */
    if (!general && ds_is_zero(Pc)) {
        Expr* cand = specialform_power_potential(Qc, xvar);
        if (cand) {
            if (sf_num_ok(P, cand, xvar, P->fun_names[0])) general = cand;
            else expr_free(cand);
        }
    }

    /* ---- Whittaker / confluent 1F1 (P == 0): single finite double pole + rank-1
     *      irregular point at infinity, emitted as the verifiable Exp z^mu 1F1 form.
     *      specialform_whittaker_basis self-verifies against the reduced equation
     *      (== the original when P == 0), so it never ships a mis-read pole. ---- */
    if (!general && ds_is_zero(Pc)) general = specialform_whittaker_basis(Qc, xvar, false);

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

    /* ---- Kummer (confluent hypergeometric), scaled: P == b/x - s, Q == -a s/x ----
     * i.e. z y'' + (b - z) y' - a y == 0 in z = s x.  x P == b - s x is affine in x
     * (b its value at x = 0, s = -its slope) and x Q == -a s is free of x.  Basis:
     *   b not an integer:  1F1[a, b, z],              z^(1-b) 1F1[a-b+1, 2-b, z]
     *   b integer >= 1:    1F1[a, b, z],              HypergeometricU[a, b, z]
     *   b integer <= 0:    z^(1-b) 1F1[a-b+1, 2-b, z], HypergeometricU[a, b, z]
     * The integer case is the DEGENERATE one -- the two 1F1 solutions coincide (or
     * one has a singular lower parameter) -- and Tricomi's U is the independent
     * second solution unless the 1F1 solution is itself a polynomial (a == -n for
     * b >= 1, a - b + 1 == -n for b <= 0), where U degenerates onto it and the
     * logarithmic ExpIntegralEi partner below is emitted instead (Laguerre
     * x y'' + (1-x) y' + k y == 0, §2.2.46-4597, is b == 1).  The scale s reaches
     * x y'' - 2x y' - y == 0 (s == 2, §2.2.47-4603).  b must be a NUMBER: a symbolic
     * exponent z^(1-b) makes the verify residual a symbolic-power + pFq sum on which
     * zero_test hangs; a stays free (symbolic). */
    if (!general) general = specialform_kummer_basis(Pc, Qc, xvar, false);

    /* ---- Kummer behind an exponential: P == p0 + p1/x, Q == q0 + q1/x, q0 != 0 ----
     * y == E^(lambda x) w with lambda^2 + p0 lambda + q0 == 0 removes the constant
     * part of Q: w'' + (P + 2 lambda) w' + (Q + lambda P + lambda^2) w == 0, whose
     * coefficients (p0 + 2 lambda) + p1/x and (q1 + lambda p1)/x are the Kummer
     * row's.  The root lambda = (-p0 - Sqrt[p0^2 - 4 q0])/2 makes the scale
     * s = -(p0 + 2 lambda) = Sqrt[disc] positive, so z = s x is real; disc must be
     * a positive constant (disc == 0 leaves no x-term in P, a negative one a
     * complex basis).  x y'' - (x - 1) y' - x y == 0 (§2.2.40-4000) is
     * E^((1 - Sqrt[5]) x/2) (C[1] 1F1[(5 - Sqrt[5])/10, 1, Sqrt[5] x] + C[2] U[..]).
     * The Liouville/Whittaker pass declines it (2 mu == 0, the degenerate case
     * that U resolves here).  Numerically gated like the other composed forms. */
    if (!general) {
        Expr* xP = ds_simplify(ds_call2(SYM_Times, expr_new_symbol(xvar), expr_copy(Pc)));
        Expr* xQ = ds_simplify(ds_call2(SYM_Times, expr_new_symbol(xvar), expr_copy(Qc)));
        Expr* p0 = ds_d(expr_copy(xP), expr_new_symbol(xvar));
        Expr* q0 = ds_d(expr_copy(xQ), expr_new_symbol(xvar));
        if (ds_free_of(p0, xvar) && ds_free_of(q0, xvar) && !ds_is_zero(q0)) {
            Expr* disc = ds_simplify(ds_call2(SYM_Subtract,
                             ds_call2(SYM_Power, expr_copy(p0), expr_new_integer(2)),
                             ds_call2(SYM_Times, expr_new_integer(4), expr_copy(q0))));
            Expr* pos = eval_and_free(ds_call1("Positive", expr_copy(disc)));
            if (pos->type == EXPR_SYMBOL && pos->data.symbol.name == SYM_True) {
                Expr* lam = ds_simplify(ds_call2(SYM_Times,
                                ds_call2(SYM_Plus, ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(p0)),
                                    ds_call2(SYM_Times, expr_new_integer(-1), ds_call1("Sqrt", expr_copy(disc)))),
                                ds_call2(SYM_Power, expr_new_integer(2), expr_new_integer(-1))));
                Expr* P2 = ds_simplify(ds_call2(SYM_Plus, expr_copy(Pc),
                               ds_call2(SYM_Times, expr_new_integer(2), expr_copy(lam))));
                Expr* Q2 = ds_simplify(ds_call2(SYM_Plus, expr_copy(Qc),
                               ds_call2(SYM_Plus, ds_call2(SYM_Times, expr_copy(lam), expr_copy(Pc)),
                                   ds_call2(SYM_Power, expr_copy(lam), expr_new_integer(2)))));
                Expr* base = specialform_kummer_basis(P2, Q2, xvar, false);
                if (base) {
                    Expr* mu = eval_and_free(ds_call2(SYM_Power, expr_new_symbol(SYM_E),
                                   ds_call2(SYM_Times, expr_copy(lam), expr_new_symbol(xvar))));
                    Expr* cand = eval_and_free(ds_call1("Expand", ds_call2(SYM_Times, mu, base)));
                    if (sf_num_ok(P, cand, xvar, P->fun_names[0])) general = cand;
                    else expr_free(cand);
                }
                expr_free(lam); expr_free(P2); expr_free(Q2);
            }
            expr_free(pos); expr_free(disc);
        }
        expr_free(xP); expr_free(xQ); expr_free(p0); expr_free(q0);
    }

    /* ---- Gauss (hypergeometric 2F1) on the canonical interval W = x(1-x) ---- */
    if (!general) general = specialform_gauss_basis(Pc, Qc, xvar);

    /* ---- affine -> Gauss 2F1: a rational-coefficient hypergeometric-class ODE
     *      whose two finite regular singular points x1, x2 are not {0, 1}.  Map
     *      x = x1 + h s (h = x2 - x1) onto the canonical interval, recognise the
     *      transformed (P~, Q~) with the Gauss helper, and compose back
     *      s = (x - x1)/h.  P~ = h P(x1+h s), Q~ = h^2 Q(x1+h s) (the phi''=0
     *      change-of-variable law).  This reaches Gegenbauer / Jacobi /
     *      associated-Legendre at SYMBOLIC (non-integer) degree — where Kovacic
     *      declines (no Liouvillian solution) — as Hypergeometric2F1 closed forms.
     *      The finite RSPs are the poles of P, Q: L = squarefree part of
     *      PolynomialLCM[Denominator P, Denominator Q]; exactly two (deg L == 2) is
     *      the Fuchsian {x1, x2, Infinity} case (deg 1 is confluent = Bessel/Kummer;
     *      deg >= 3 is Heun = future).  A wrong RSP map yields a 2F1 whose residual
     *      zero_test cannot disprove, so hgc_num_ok (numeric back-substitution on
     *      the ORIGINAL equation) is the mandatory guard of the 0-FAIL invariant. */
    if (!general) general = specialform_affine_gauss(P, Pc, Qc, xvar, false);

    /* ---- Liouville normal-form pre-pass: a P != 0 equation whose y'-free normal
     *      form z'' = r z (with y = mu z, mu = Exp[-Integrate[P/2]]) is one of the
     *      reduced-form recognisers above (Airy / Bessel-reducible).  Lets those
     *      rows fire on an equation carrying a first-derivative term.  mu must be
     *      elementary (dsolve_normal_form returns NULL otherwise -> skip).  Gated
     *      by sf_num_ok as a defensive numeric back-substitution: a sign/branch
     *      error in mu would otherwise reach the corpus as a numericizing FAIL. */
    if (!general && !ds_is_zero(Pc) && (g_dsolve_depth <= 1 || dsolve_nh_vop_active || dsolve_cv_active)) {
        /* Top-level only (g_dsolve_depth <= 1): a recursive caller such as
         * OperatorFactor peels a first-order factor and re-solves the order-(n-1)
         * quotient, and the mu = Exp[-Int P/2] recovery factor this pre-pass emits
         * composed back into that reduction can drive the evaluator into an infinite
         * rewrite ($IterationLimit).  The pre-pass only ever helps a genuine top-level
         * equation carrying a y' term, so restricting it to depth 1 loses nothing and
         * keeps OperatorFactor / Riccati recursion clean (cf. dsolve_common.c:642).
         * EXCEPTION: DSolve`VariationOfParameters (dsolve_nh_vop_active) recurses to
         * solve the HOMOGENEOUS part and needs this recogniser for a Bessel/Airy set
         * carrying a y' term (e.g. 4x^2 y'' - 4x y' + (3-16x^2) y == 8 x^(5/2)); it
         * extracts the basis rather than composing mu back into a reduction, so the
         * infinite-rewrite hazard does not apply.  DSolve`ChangeOfVariable (dsolve_cv_active)
         * shares the exemption for the same reason: it composes the transformed basis
         * Y(t) as Y(phi(x)) (t = E^(k x) leaves a Y' term, e.g. Y'' + Y'/t + a/(b^2 t) Y). */
        Expr* mu = NULL;
        Expr* r  = dsolve_normal_form(Pc, Qc, xvar, &mu);
        if (mu) {
            Expr* negr = ds_simplify(ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(r)));  /* -r */
            /* Complexity backstop: the reduced-form recognisers Simplify the potential
             * several times; on a large rational normal form (e.g. a quadratic-over-
             * quadratic-squared coefficient with several symbolic parameters) that runs
             * to seconds and the pre-pass declines anyway.  Every genuine Airy/Bessel-
             * reducible normal form is small (single power / A + B/x^2 / constant). */
            Expr* base = (leaf_count_internal(negr, true) > 50)
                       ? NULL : specialform_reduced_basis(negr, xvar);
            /* Generalised power-potential recogniser (symbolic-exponent Bessel /
             * two-term confluent), on the normal-form potential of a y'-carrying
             * equation whose reduced form is a power potential -- e.g.
             * y'' + a x^n y' + b x^(n-1) y == 0, normal form
             * -a^2/4 x^(2n) + (b - a n/2) x^(n-1).  Unlike the finite-pole Whittaker
             * recogniser the NOTE below excludes, this family's z = kappa x^d is a
             * MONOMIAL (no finite (x - x0) pole) and the recovery factor
             * mu = Exp[-Int P/2] for a power/rational P is itself a monomial-exponent
             * form, so the composed mu*base does not stack finite-pole radicals and
             * the $IterationLimit hazard the NOTE describes does not apply.  Still
             * gated by sf_num_ok below. */
            if (!base && leaf_count_internal(negr, true) <= 50)
                base = specialform_power_potential(negr, xvar);
            /* Confluent (Whittaker / Kummer 1F1) finite-pole normal form, NOW ENABLED
             * in the P != 0 pre-pass.  An earlier revision disabled it here for two
             * reasons, both since resolved:
             *   (1) The recovery factor mu = Exp[-Int P/2] carries (x - x0)^a and the
             *       Whittaker factor carries (c (x - x0))^(1/2 +- mu) -- a DIFFERENT
             *       base (the scale c is inside), so the composed mu*base stacked two
             *       symbolic-radical powers that the evaluator could thrash merging.
             *       PowerExpand[cand] below pulls c out so the two powers share the base
             *       (x - x0) and collapse to a single exponent (4067: (5x)^(5/4)/x^(5/4)
             *       -> 5 5^(1/4) ...; 4009: x^(-1/2) (x)^(1/2+Sqrt7) -> x^Sqrt7).
             *   (2) It was said to "starve the Frobenius series fallback" -- but M68
             *       demoted a SeriesData answer out of the automatic cascade, so that
             *       fallback no longer yields an accepted automatic result; a verified
             *       Whittaker closed form can now only turn an UNEVAL into a PASS.
             * specialform_whittaker_basis already self-verifies (whittaker_num_ok) and
             * declines the degenerate 2 mu in Z case (dependent / singular 1F1), and the
             * composed form is gated again by sf_num_ok -- so a mis-read pole/scale can
             * only decline, never emit a wrong closed form. */
            bool whit = false;
            if (!base && leaf_count_internal(negr, true) <= 50) {
                base = specialform_whittaker_basis(negr, xvar, false);
                if (base) whit = true;
            }
            expr_free(negr);
            if (base) {
                /* mu C[1] z0 + mu C[2] z1: distribute the recovery factor over the
                 * basis terms so it meets each term's own powers of the same base
                 * and cancels against them (4009: E^(x/2) x^(-1/2) times
                 * E^(-x/2) x^(1/2 + Sqrt[7]) is x^Sqrt[7]).  Left as an outer
                 * factor mu (C[1] z0 + C[2] z1), nothing ever combines them. */
                Expr* cand;
                if (base->type == EXPR_FUNCTION &&
                    base->data.function.head->type == EXPR_SYMBOL &&
                    base->data.function.head->data.symbol.name == SYM_Plus) {
                    size_t nb = base->data.function.arg_count;
                    Expr** terms = malloc(nb * sizeof(Expr*));
                    for (size_t k = 0; k < nb; k++)
                        terms[k] = ds_call2(SYM_Times, expr_copy(mu),
                                            expr_copy(base->data.function.args[k]));
                    cand = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), terms, nb));
                    free(terms);
                    expr_free(mu); expr_free(base);
                } else {
                    cand = eval_and_free(ds_call2(SYM_Times, mu, base));
                }
                if (whit)   /* merge the stacked (x - x0) radicals: see (1) above */
                    cand = eval_and_free(ds_call1("PowerExpand", cand));
                if (sf_num_ok(P, cand, xvar, P->fun_names[0])) general = cand;
                else expr_free(cand);
            } else {
                expr_free(mu);
            }
        }
        expr_free(r);
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

/* ===================================================================
 * DSolve`SpecialFunctionTransform -- the LATE special-function slot.
 *
 * SpecialFunctionForm above matches a table of canonical forms in the variable
 * the equation is written in.  A large polynomial-coefficient class reaches those
 * rows only after a change of the dependent and/or independent variable that the
 * table does not try.  This method tries three, each composed back and gated by
 * sf_num_ok (numeric back-substitution into the ORIGINAL equation with every
 * symbolic parameter instantiated), so a mis-read transform can only decline:
 *
 *  (1) Linear coefficients (Laplace's equation, the Kummer class):
 *        (a2 x + b2) y'' + (a1 x + b1) y' + (a0 x + b0) y == 0,  a2 != 0.
 *      Translate the single finite singular point x0 to the origin (t = x - x0),
 *      so  t y'' + (A t + B) y' + (C t + D) y == 0.  y = E^(lam t) w with
 *      lam^2 + A lam + C == 0 removes the t-part of the y-coefficient:
 *        t w'' + ((2 lam + A) t + B) w' + (lam B + D) w == 0,
 *      Kummer's equation in z = s t, s = -(2 lam + A) = Sqrt[A^2 - 4 C] (lam = 0
 *      when C == 0).  A double root (A^2 == 4 C) leaves t w'' + B w' + K w == 0,
 *      K = lam B + D: Bessel, w = t^((1-B)/2) Z_(B-1)(2 Sqrt[K t]).
 *      A symbolic b = B is allowed here (basis {1F1, U}); the early row needs a
 *      numeric b.
 *  (2) The Liouville normal form q = -r (w'' + q w == 0) read by the Whittaker row
 *      with the DEGENERATE case 2 mu in Z admitted (HypergeometricU partner).
 *  (3) Quadratic argument: q even about a centre x0 (the single pole of q, or the
 *      vertex of a polynomial q).  z = (x - x0)^2 with w = (x - x0)^(-1/2) v(z)
 *      gives v'' + Qz v == 0, Qz = (q + 3/(4 (x-x0)^2)) / (4 (x-x0)^2) as a
 *      rational function of z (the Schwarzian term of the Liouville
 *      transformation), read by the Whittaker / reduced / power-potential rows.
 *      Reaches Weber's equation y'' + (c x^2 + b x + a) y == 0, the Hermite type
 *      y'' + (b x + a) y' + c y == 0, and the Bessel operator with an even
 *      quartic potential x^2 y'' + x y' - (c^2 x^4 + b^2 x^2 + a^2) y == 0.
 *
 * Runs after Kovacic and ChangeOfVariable (a Liouvillian/elementary answer, when
 * one exists, is found first and is simpler).  Not recursive.  Killswitch
 * MATHILDA_NO_SF_TRANSFORM.
 * =================================================================== */

/* squarefree part of the common denominator of (a, b) in x; owned. */
static Expr* sft_pole_rad(const Expr* a, const Expr* b, const char* xv) {
    Expr* da = eval_and_free(ds_call1("Denominator", ds_call1("Together", expr_copy((Expr*)a))));
    Expr* db = eval_and_free(ds_call1("Denominator", ds_call1("Together", expr_copy((Expr*)b))));
    Expr* u  = eval_and_free(ds_call2("PolynomialLCM", da, db));
    Expr* du = ds_d(expr_copy(u), expr_new_symbol(xv));
    Expr* g  = eval_and_free(ds_call2("PolynomialGCD", expr_copy(u), du));
    return sf_ct(ds_call2(SYM_Times, u, powneg1_sf(g)));
}
/* degree of e in x, or -1 unless e is a polynomial in x. */
static int sft_degree(const Expr* e, const char* xv) {
    if (e && e->type == EXPR_INTEGER && e->data.integer == 0) return 0;   /* the zero polynomial */
    Expr* pq = eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ),
                   (Expr*[]){ expr_copy((Expr*)e), expr_new_symbol(xv) }, 2));
    bool ok = pq->type == EXPR_SYMBOL && pq->data.symbol.name == SYM_True;
    expr_free(pq);
    if (!ok) return -1;
    Expr* d = eval_and_free(expr_new_function(expr_new_symbol("Exponent"),
                  (Expr*[]){ expr_copy((Expr*)e), expr_new_symbol(xv) }, 2));
    int r = (d->type == EXPR_INTEGER && d->data.integer >= 0 && d->data.integer < 64)
          ? (int)d->data.integer : -1;
    expr_free(d);
    return r;
}
static Expr* sft_coef(const Expr* e, const char* xv, int k) {
    return eval_and_free(expr_new_function(expr_new_symbol("Coefficient"),
               (Expr*[]){ expr_copy((Expr*)e), expr_new_symbol(xv), expr_new_integer(k) }, 3));
}
/* -c0/c1 for the linear polynomial rad = c1 x + c0; owned. */
static Expr* sft_root1(const Expr* rad, const char* xv) {
    Expr* r0 = sft_coef(rad, xv, 0);
    Expr* r1 = sft_coef(rad, xv, 1);
    return sf_ct(ds_call2(SYM_Times, expr_new_integer(-1), ds_call2(SYM_Times, r0, powneg1_sf(r1))));
}
static Expr* sft_half(int num) {   /* num/2 */
    return eval_and_free(ds_call2(SYM_Times, expr_new_integer(num),
               powneg1_sf(expr_new_integer(2))));
}
/* f * (C[1] b0 + C[2] b1), distributed over the basis terms so the factor meets
 * each term's own powers (as the pre-pass above does).  f, base consumed. */
static Expr* sft_distribute(Expr* f, Expr* base) {
    if (base->type == EXPR_FUNCTION && base->data.function.head->type == EXPR_SYMBOL &&
        base->data.function.head->data.symbol.name == SYM_Plus) {
        size_t nb = base->data.function.arg_count;
        Expr** terms = malloc(nb * sizeof(Expr*));
        for (size_t k = 0; k < nb; k++)
            terms[k] = ds_call2(SYM_Times, expr_copy(f), expr_copy(base->data.function.args[k]));
        Expr* r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), terms, nb));
        free(terms);
        expr_free(f); expr_free(base);
        return r;
    }
    return eval_and_free(ds_call2(SYM_Times, f, base));
}

/* (1) linear coefficients.  Pc, Qc borrowed. */
static Expr* sft_laplace(DSolveProblem* P, const Expr* Pc, const Expr* Qc, const char* xv) {
    Expr* rad = sft_pole_rad(Pc, Qc, xv);
    int drad = sft_degree(rad, xv);
    Expr* out = NULL;
    if (drad != 1) { expr_free(rad); return NULL; }
    Expr* x0 = sft_root1(rad, xv);
    expr_free(rad);
    Expr* xs = eval_and_free(ds_call2(SYM_Plus, expr_new_symbol(xv), expr_copy(x0)));
    Expr* Pt = sf_ct(ds_call2(SYM_Times, expr_new_symbol(xv),
                   ds_subst(expr_copy((Expr*)Pc), expr_new_symbol(xv), expr_copy(xs))));
    Expr* Qt = sf_ct(ds_call2(SYM_Times, expr_new_symbol(xv),
                   ds_subst(expr_copy((Expr*)Qc), expr_new_symbol(xv), expr_copy(xs))));
    expr_free(xs);
    int dp = sft_degree(Pt, xv), dq = sft_degree(Qt, xv);
    if (dp >= 0 && dp <= 1 && dq >= 0 && dq <= 1) {
        Expr* A = sft_coef(Pt, xv, 1); Expr* B = sft_coef(Pt, xv, 0);
        Expr* C = sft_coef(Qt, xv, 1); Expr* D = sft_coef(Qt, xv, 0);
        Expr* disc = sf_ct(ds_call2(SYM_Subtract,
                         ds_call2(SYM_Power, expr_copy(A), expr_new_integer(2)),
                         ds_call2(SYM_Times, expr_new_integer(4), expr_copy(C))));
        Expr* lam = NULL; Expr* base = NULL;
        if (ds_is_zero(disc)) {
            /* double root: t w'' + B w' + K w == 0 -> t^((1-B)/2) Z_(B-1)(2 Sqrt[K t]) */
            lam = sf_ct(ds_call2(SYM_Times, sft_half(-1), expr_copy(A)));
            Expr* K = sf_ct(ds_call2(SYM_Plus, ds_call2(SYM_Times, expr_copy(lam), expr_copy(B)),
                                     expr_copy(D)));
            if (!ds_is_zero(K)) {
                Expr* nu  = eval_and_free(ds_call2(SYM_Subtract, expr_copy(B), expr_new_integer(1)));
                Expr* arg = eval_and_free(ds_call2(SYM_Times, expr_new_integer(2),
                                ds_call1("Sqrt", ds_call2(SYM_Times, expr_copy(K), expr_new_symbol(xv)))));
                Expr* pre = eval_and_free(ds_call2(SYM_Power, expr_new_symbol(xv),
                                ds_call2(SYM_Times, sft_half(1),
                                    ds_call2(SYM_Subtract, expr_new_integer(1), expr_copy(B)))));
                Expr* b0 = ds_call2(SYM_Times, expr_copy(pre), expr_new_function(expr_new_symbol("BesselJ"),
                               (Expr*[]){ expr_copy(nu), expr_copy(arg) }, 2));
                Expr* b1 = ds_call2(SYM_Times, expr_copy(pre), expr_new_function(expr_new_symbol("BesselY"),
                               (Expr*[]){ expr_copy(nu), expr_copy(arg) }, 2));
                base = combo(b0, b1);
                expr_free(nu); expr_free(arg); expr_free(pre);
            }
            expr_free(K);
        } else {
            lam = ds_is_zero(C) ? expr_new_integer(0)
                : sf_ct(ds_call2(SYM_Times, sft_half(1),
                      ds_call2(SYM_Subtract, ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(A)),
                          ds_call1("Sqrt", expr_copy(disc)))));
            /* P2 = (2 lam + A) + B/t,  Q2 = (lam B + D)/t */
            Expr* P2 = eval_and_free(ds_call2(SYM_Plus,
                           ds_call2(SYM_Plus, ds_call2(SYM_Times, expr_new_integer(2), expr_copy(lam)),
                                    expr_copy(A)),
                           ds_call2(SYM_Times, expr_copy(B), powneg1_sf(expr_new_symbol(xv)))));
            Expr* Q2 = eval_and_free(ds_call2(SYM_Times,
                           ds_call2(SYM_Plus, ds_call2(SYM_Times, expr_copy(lam), expr_copy(B)),
                                    expr_copy(D)),
                           powneg1_sf(expr_new_symbol(xv))));
            base = specialform_kummer_basis(P2, Q2, xv, true);
            expr_free(P2); expr_free(Q2);
        }
        if (base) {
            Expr* ef = eval_and_free(ds_call1("Exp", ds_call2(SYM_Times, expr_copy(lam),
                                                               expr_new_symbol(xv))));
            Expr* cand = sft_distribute(ef, base);
            Expr* back = eval_and_free(ds_call2(SYM_Subtract, expr_new_symbol(xv), expr_copy(x0)));
            cand = ds_subst(cand, expr_new_symbol(xv), back);
            if (sf_num_ok(P, cand, xv, P->fun_names[0])) out = cand;
            else expr_free(cand);
        }
        if (lam) expr_free(lam);
        expr_free(A); expr_free(B); expr_free(C); expr_free(D); expr_free(disc);
    }
    expr_free(Pt); expr_free(Qt); expr_free(x0);
    return out;
}

/* y = f * base, PowerExpand-merged (stacked radicals of one base), gated.
 * f, base consumed. */
static Expr* sft_compose(DSolveProblem* P, Expr* f, Expr* base, const char* xv) {
    Expr* cand = sft_distribute(f, base);
    cand = eval_and_free(ds_call1("PowerExpand", cand));
    if (sf_num_ok(P, cand, xv, P->fun_names[0])) return cand;
    expr_free(cand);
    return NULL;
}

/* (3) quadratic argument on the normal form q (w'' + q w == 0).  q, mu borrowed
 * (mu the recovery factor, y = mu w). */
static Expr* sft_quadarg(DSolveProblem* P, const Expr* q, const Expr* mu, const char* xv) {
    Expr* zero = expr_new_integer(0);
    Expr* rad = sft_pole_rad(q, zero, xv);
    expr_free(zero);
    int drad = sft_degree(rad, xv);
    Expr* x0 = NULL;
    if (drad == 1) x0 = sft_root1(rad, xv);
    else if (drad == 0) {
        int d = sft_degree(q, xv);
        if (d >= 2 && d % 2 == 0 && d <= 8) {    /* centre = vertex of the even-degree polynomial */
            Expr* cd = sft_coef(q, xv, d); Expr* cd1 = sft_coef(q, xv, d - 1);
            x0 = sf_ct(ds_call2(SYM_Times, expr_new_integer(-1),
                     ds_call2(SYM_Times, cd1, powneg1_sf(ds_call2(SYM_Times, expr_new_integer(d), cd)))));
        }
    }
    expr_free(rad);
    if (!x0) return NULL;
    Expr* out = NULL;
    Expr* xs = eval_and_free(ds_call2(SYM_Plus, expr_new_symbol(xv), expr_copy(x0)));
    Expr* qs = sf_ct(ds_subst(expr_copy((Expr*)q), expr_new_symbol(xv), xs));     /* q(x0 + u) */
    Expr* qm = ds_subst(expr_copy(qs), expr_new_symbol(xv),
                   eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), expr_new_symbol(xv))));
    Expr* odd = sf_ct(ds_call2(SYM_Subtract, expr_copy(qs), qm));
    bool even = ds_is_zero(odd);
    expr_free(odd);
    if (even && !ds_free_of(qs, xv)) {
        /* Qz = (q(Sqrt z) + 3/(4 z)) / (4 z) */
        Expr* qz = ds_subst(expr_copy(qs), expr_new_symbol(xv), ds_call1("Sqrt", expr_new_symbol(xv)));
        Expr* three4 = eval_and_free(ds_call2(SYM_Times, expr_new_integer(3),
                           powneg1_sf(ds_call2(SYM_Times, expr_new_integer(4), expr_new_symbol(xv)))));
        Expr* Qz = sf_ct(ds_call2(SYM_Times, ds_call2(SYM_Plus, qz, three4),
                         powneg1_sf(ds_call2(SYM_Times, expr_new_integer(4), expr_new_symbol(xv)))));
        Expr* num = eval_and_free(ds_call1("Numerator", expr_copy(Qz)));
        Expr* den = eval_and_free(ds_call1("Denominator", expr_copy(Qz)));
        bool rat = sft_degree(num, xv) >= 0 && sft_degree(den, xv) >= 0;
        expr_free(num); expr_free(den);
        Expr* base = NULL;
        if (rat && leaf_count_internal(Qz, true) <= 80) {
            base = specialform_whittaker_basis(Qz, xv, true);
            if (!base) base = specialform_reduced_basis(Qz, xv);
            if (!base) base = specialform_power_potential(Qz, xv);
        }
        expr_free(Qz);
        if (base) {
            /* v(z) -> v((x - x0)^2);  y = mu (x - x0)^(-1/2) v */
            Expr* t = eval_and_free(ds_call2(SYM_Subtract, expr_new_symbol(xv), expr_copy(x0)));
            base = ds_subst(base, expr_new_symbol(xv),
                       eval_and_free(ds_call2(SYM_Power, expr_copy(t), expr_new_integer(2))));
            Expr* f = eval_and_free(ds_call2(SYM_Times, expr_copy((Expr*)mu),
                          ds_call2(SYM_Power, t, sft_half(-1))));
            out = sft_compose(P, f, base, xv);
        }
    }
    expr_free(qs); expr_free(x0);
    return out;
}

Expr** dsolve_specialform_transform_try(DSolveProblem* P, size_t* nbranch) {
    if (getenv("MATHILDA_NO_SF_TRANSFORM")) return NULL;
    Expr* Pc; Expr* Qc;
    if (!dsolve_second_order_PQ(P, &Pc, &Qc)) return NULL;
    const char* xv = P->ind_names[0];
    Expr* general = NULL;
    if (ds_is_rational_in(Pc, xv) && ds_is_rational_in(Qc, xv) &&
        leaf_count_internal(Pc, true) + leaf_count_internal(Qc, true) <= 120) {
        general = sft_laplace(P, Pc, Qc, xv);
        if (!general) general = specialform_affine_gauss(P, Pc, Qc, xv, true);
        if (!general) {
            Expr* mu = NULL;
            Expr* r = dsolve_normal_form(Pc, Qc, xv, &mu);
            if (mu && r) {
                Expr* q = sf_ct(ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(r)));
                if (leaf_count_internal(q, true) <= 80) {
                    Expr* base = specialform_whittaker_basis(q, xv, true);
                    if (base) general = sft_compose(P, expr_copy(mu), base, xv);
                    if (!general) general = sft_quadarg(P, q, mu, xv);
                }
                expr_free(q);
            }
            if (mu) expr_free(mu);
            if (r) expr_free(r);
        }
    }
    expr_free(Pc); expr_free(Qc);
    if (!general) return NULL;
    Expr** out = malloc(sizeof(Expr*));
    out[0] = general;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_specialform_transform(Expr* res) {
    return dsolve_method_builtin(res, dsolve_specialform_transform_try);
}

static Expr* builtin_dsolve_specialform(Expr* res) {
    return dsolve_method_builtin(res, dsolve_specialform_try);
}

void dsolve_specialform_init(void) {
    symtab_add_builtin("DSolve`SpecialFunctionTransform", builtin_dsolve_specialform_transform);
    symtab_get_def("DSolve`SpecialFunctionTransform")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`SpecialFunctionTransform",
        "DSolve`SpecialFunctionTransform[eqn, y, x] solves a homogeneous second-order linear "
        "ODE with rational coefficients that reaches a confluent special-function row only "
        "after a change of variables: linear coefficients (a2 x + b2) y'' + (a1 x + b1) y' + "
        "(a0 x + b0) y == 0 (shift to the singular point, y = E^(lam x) w -> Kummer "
        "Hypergeometric1F1/HypergeometricU, or Bessel at a double root), the degenerate "
        "Whittaker normal form (HypergeometricU partner), and a normal form even about a "
        "centre x0 read in z = (x - x0)^2 (Weber / Hermite / Whittaker in x^2). Every answer "
        "is numerically back-substituted with the parameters sampled.");
    symtab_add_builtin("DSolve`SpecialFunctionForm", builtin_dsolve_specialform);
    symtab_get_def("DSolve`SpecialFunctionForm")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`SpecialFunctionForm",
        "DSolve`SpecialFunctionForm[eqn, y, x] recognises second-order linear ODEs "
        "whose solutions are named special functions: Airy (y'' == (A x + B) y), "
        "the power potential y'' + A x^m y == 0 (symbolic m -> Sqrt[x] Bessel) and the "
        "two-term y'' + (a x^(2c) + b x^(c-1)) y == 0 (-> Coulomb/Whittaker -> "
        "Hypergeometric1F1), "
        "Bessel / modified Bessel (x^2 y'' + x y' +- (x^2 -+ v^2) y == 0), Kummer "
        "confluent hypergeometric (x y'' + (b - s x) y' - a s y == 0 -> "
        "Hypergeometric1F1[a, b, s x], with HypergeometricU for integer b), and Gauss hypergeometric "
        "(x(1-x) y'' + (c - (a+b+1) x) y' - a b y == 0 -> Hypergeometric2F1), ordinary "
        "Legendre, and the trigonometric Pöschl-Teller potential "
        "(y'' == (a + p(p-1) Csc^2 x + q(q-1) Sec^2 x) y -> Hypergeometric2F1, "
        "numerically verified). Any rational-coefficient hypergeometric-class equation "
        "with two finite regular singular points {x1, x2} (Gegenbauer, Jacobi, "
        "associated Legendre at symbolic degree) is mapped affinely onto x(1-x) and the "
        "local exponents pulled out (Y = s^r0 (1-s)^r1 F) to reach Hypergeometric2F1, "
        "gated by a numeric self-verify. Equations carrying a y' term are also tried "
        "through the Liouville normal form (y = z Exp[-Int P/2]) against the y'-free "
        "Airy/Bessel recognisers. A normal form with a single finite double pole and a "
        "rank-1 irregular point at infinity is Whittaker's equation, emitted as the "
        "confluent form Exp[-z/2] z^(1/2+-mu) Hypergeometric1F1[1/2+-mu-kappa, 1+-2mu, z] "
        "(declining when 2 mu is an integer). The hypergeometric second solution is "
        "emitted only when b (resp. c) is not an integer; otherwise it declines to the "
        "series fallback.");
}
