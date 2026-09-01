/*
 * dsolve_kovacic.c — DSolve`Kovacic.
 *
 * Liouvillian solutions of a homogeneous second-order linear ODE
 *     y'' + P(x) y' + Q(x) y == 0
 * via the Kovacic reduction to the normal form  z'' == r z  (r rational,
 * y == w z with w == Exp[-Integrate[P/2]]) and a search for the logarithmic
 * derivative ω of a solution z == Exp[Integrate[ω]].
 *
 * Case 1 (ω ∈ C(x)):  ω satisfies the Riccati equation ω' + ω² == r.  We posit ω
 * with undetermined coefficients matching the pole structure of r (a polynomial
 * part whose degree is fixed by the order of r at ∞, plus a principal part
 * Σ N_{i,k}(x)/f_i(x)^k over the irreducible factors f_i of the denominator of r,
 * with k up to ceil(o_i/2)), clear denominators, and solve ω'+ω²-r == 0 for the
 * coefficients (Solve).  A solution gives z1 = Exp[Integrate[ω]] and the second
 * solution z2 = z1 Integrate[1/z1^2] (the Wronskian is constant for z''==r z).
 * These recover elementary y-solutions, which the substrate verifies by
 * back-substitution.
 *
 * Case 2 (ω algebraic of degree 2 over C(x)):  the pair {ω1,ω2} of logarithmic
 * derivatives has σ == ω1+ω2 ∈ C(x) satisfying D' + 2σD == 0 with
 * D == 4r - 2σ' - σ²; a rational σ (same undetermined-coefficient search) yields
 * ω == (σ ± Sqrt[D])/2.  Case-2 solutions are algebraic, so their residual is not
 * decidable symbolically; the method numerically back-substitutes each candidate
 * before returning it, so a wrong candidate is never shipped.
 *
 * Case 3 (degree 4/6/12) is not implemented: the method declines.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include "../internal.h"
#include "../parse.h"
#include "integrate.h"          /* g_integrate_quiet: silence speculative nonelem */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* ---- small evaluated builders (args consumed, result owned) ---- */
static Expr* T2(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Times, a, b)); }
static Expr* A2(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Plus,  a, b)); }
static Expr* Sub(Expr* a, Expr* b){ return eval_and_free(ds_call2(SYM_Subtract, a, b)); }
static Expr* Neg(Expr* a)         { return eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), a)); }
static Expr* Powi(Expr* b, int e) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ b, expr_new_integer(e) }, 2));
}
static Expr* fn1(const char* h, Expr* a)          { return eval_and_free(ds_call1(h, a)); }
static Expr* fn2(const char* h, Expr* a, Expr* b) { return eval_and_free(ds_call2(h, a, b)); }

/* integer degree of poly in x (Exponent[poly, x]); -inf reported as -1 */
static int degree_in(const Expr* poly, const char* x) {
    Expr* e = fn2("Exponent", expr_copy((Expr*)poly), expr_new_symbol(x));
    int d = (e->type == EXPR_INTEGER) ? (int)e->data.integer : -1;
    expr_free(e);
    return d;
}

/* Coefficient[e, x, k] (owned). */
static Expr* coeff_k(const Expr* e, const char* x, int k) {
    return eval_and_free(expr_new_function(expr_new_symbol("Coefficient"),
               (Expr*[]){ expr_copy((Expr*)e), expr_new_symbol(x), expr_new_integer(k) }, 3));
}

/* Unevaluated x^k. */
static Expr* xpow(const char* x, int k) {
    return expr_new_function(expr_new_symbol(SYM_Power),
               (Expr*[]){ expr_new_symbol(x), expr_new_integer(k) }, 2);
}

/* Try one Riccati/σ ansatz `w` (in x and the unknown symbols `unk[0..nu-1]`) that
 * must satisfy `eqexpr == 0`.  Returns the solved, fully-determined expression
 * (unknowns fixed by Solve, any free ones set to 0), or NULL if no solution. */
static Expr* solve_ansatz(Expr* eqexpr, const char* x, Expr** unk, size_t nu, Expr* w) {
    Expr* num = fn1("Numerator", fn1("Together", eqexpr));   /* consumes eqexpr */
    Expr* clist = fn2("CoefficientList", num, expr_new_symbol(x));
    if (!clist || !head_is(clist, SYM_List)) { if (clist) expr_free(clist); expr_free(w); return NULL; }
    size_t nc = clist->data.function.arg_count;
    Expr** eqs = malloc((nc ? nc : 1) * sizeof(Expr*));
    size_t ne = 0;
    for (size_t i = 0; i < nc; i++) {
        Expr* co = clist->data.function.args[i];
        if (ds_is_zero(co)) continue;   /* a 0 coefficient is a trivially-True equation */
        eqs[ne++] = expr_new_function(expr_new_symbol(SYM_Equal),
                        (Expr*[]){ expr_copy(co), expr_new_integer(0) }, 2);
    }
    expr_free(clist);
    Expr* eqlist = expr_new_function(expr_new_symbol(SYM_List), eqs, ne);
    free(eqs);

    /* No unknowns: the ansatz is fully fixed, so it is a solution iff every
     * coefficient equation was trivially satisfied (none survived).  Handle this
     * directly rather than calling Solve[eqs, {}] (which errors "ivar: {} is not
     * a valid variable"). */
    if (nu == 0) {
        bool consistent = (ne == 0);
        expr_free(eqlist);
        if (consistent) return ds_simplify(w);
        expr_free(w); return NULL;
    }

    Expr** vs = malloc(nu * sizeof(Expr*));
    for (size_t i = 0; i < nu; i++) vs[i] = expr_copy(unk[i]);
    Expr* varlist = expr_new_function(expr_new_symbol(SYM_List), vs, nu);
    free(vs);

    Expr* sol = ds_solve(eqlist, varlist);
    Expr* out = NULL;
    if (sol && head_is(sol, SYM_List) && sol->data.function.arg_count >= 1) {
        Expr* branch = sol->data.function.args[0];
        if (head_is(branch, SYM_List)) {
            out = eval_and_free(internal_replace_all(
                      (Expr*[]){ expr_copy(w), expr_copy(branch) }, 2));
            /* any unknown not fixed by Solve: set it to 0 */
            for (size_t i = 0; i < nu; i++)
                if (ds_contains(out, unk[i]->data.symbol.name))
                    out = ds_subst(out, expr_copy(unk[i]), expr_new_integer(0));
            out = ds_simplify(out);
        }
    }
    if (sol) expr_free(sol);
    expr_free(w);
    return out;
}

/* Solve the coefficients of x^k (k = kmin..kmax) of `g` to zero for the unknowns
 * `unk[0..nu-1]`, and return `w` with the solution substituted (free unknowns set
 * to 0), or NULL if inconsistent.  `g` and `w` are consumed.  Used to fit the
 * high-order part of ω_base (partial matching, unlike solve_ansatz's full one). */
static Expr* solve_high_coeffs(Expr* g, const char* x, int kmin, int kmax,
                               Expr** unk, size_t nu, Expr* w) {
    Expr** eqs = malloc((size_t)(kmax - kmin + 1) * sizeof(Expr*));
    size_t ne = 0;
    for (int k = kmax; k >= kmin; k--) {
        Expr* co = coeff_k(g, x, k);
        if (ds_is_zero(co)) { expr_free(co); continue; }
        eqs[ne++] = expr_new_function(expr_new_symbol(SYM_Equal),
                        (Expr*[]){ co, expr_new_integer(0) }, 2);
    }
    expr_free(g);
    Expr* eqlist = expr_new_function(expr_new_symbol(SYM_List), eqs, ne);
    free(eqs);
    Expr** vs = malloc((nu ? nu : 1) * sizeof(Expr*));
    for (size_t i = 0; i < nu; i++) vs[i] = expr_copy(unk[i]);
    Expr* varlist = expr_new_function(expr_new_symbol(SYM_List), vs, nu);
    free(vs);
    Expr* sol = ds_solve(eqlist, varlist);
    Expr* out = NULL;
    if (sol && head_is(sol, SYM_List) && sol->data.function.arg_count >= 1) {
        Expr* branch = sol->data.function.args[0];
        if (head_is(branch, SYM_List)) {
            out = eval_and_free(internal_replace_all(
                      (Expr*[]){ expr_copy(w), expr_copy(branch) }, 2));
            for (size_t i = 0; i < nu; i++)
                if (ds_contains(out, unk[i]->data.symbol.name))
                    out = ds_subst(out, expr_copy(unk[i]), expr_new_integer(0));
            out = ds_simplify(out);
        }
    }
    if (sol) expr_free(sol);
    expr_free(w);
    return out;
}

/* Build ω = Σ_{j=0}^d k_j x^j  +  Σ_i Σ_{k=1}^{m_i} (Σ_{l<deg f_i} k x^l)/f_i^k,
 * appending each fresh unknown symbol to `unk` (grown in place).  `factors` is a
 * FactorList result; `poly_deg` is d.  Returns the ansatz expression (unevaluated
 * sum) and fills *nu / *unk. */
static Expr* build_riccati_ansatz(const char* x, int poly_deg,
                                  Expr* factors, int* counter,
                                  Expr*** unk, size_t* nu) {
    size_t cap = 8; Expr** U = malloc(cap * sizeof(Expr*)); size_t n = 0;
    size_t tcap = 8; Expr** terms = malloc(tcap * sizeof(Expr*)); size_t nt = 0;
    #define PUSH_U(sym) do { if (n==cap){cap*=2;U=realloc(U,cap*sizeof(Expr*));} U[n++]=(sym); } while(0)
    #define PUSH_T(t)  do { if (nt==tcap){tcap*=2;terms=realloc(terms,tcap*sizeof(Expr*));} terms[nt++]=(t); } while(0)

    /* polynomial part */
    for (int j = 0; j <= poly_deg; j++) {
        char buf[32]; snprintf(buf, sizeof(buf), "DSolve`kv%d", (*counter)++);
        const char* sn = intern_symbol(buf);
        PUSH_U(expr_new_symbol(sn));
        PUSH_T(ds_call2(SYM_Times, expr_new_symbol(sn),
                        expr_new_function(expr_new_symbol(SYM_Power),
                            (Expr*[]){ expr_new_symbol(x), expr_new_integer(j) }, 2)));
    }
    /* pole parts, one per irreducible factor (skip the constant content) */
    if (factors && head_is(factors, SYM_List)) {
        for (size_t fi = 0; fi < factors->data.function.arg_count; fi++) {
            Expr* pair = factors->data.function.args[fi];
            if (!head_is(pair, SYM_List) || pair->data.function.arg_count != 2) continue;
            Expr* f = pair->data.function.args[0];
            Expr* om = pair->data.function.args[1];
            int df = degree_in(f, x);
            if (df < 1) continue;                        /* constant content */
            int o = (om->type == EXPR_INTEGER) ? (int)om->data.integer : 1;
            int m = (o + 1) / 2;                          /* ceil(o/2) */
            for (int k = 1; k <= m; k++) {
                /* numerator poly of degree < df */
                for (int l = 0; l < df; l++) {
                    char buf[32]; snprintf(buf, sizeof(buf), "DSolve`kv%d", (*counter)++);
                    const char* sn = intern_symbol(buf);
                    PUSH_U(expr_new_symbol(sn));
                    Expr* term = ds_call2(SYM_Times, expr_new_symbol(sn),
                                     ds_call2(SYM_Times,
                                         expr_new_function(expr_new_symbol(SYM_Power),
                                             (Expr*[]){ expr_new_symbol(x), expr_new_integer(l) }, 2),
                                         expr_new_function(expr_new_symbol(SYM_Power),
                                             (Expr*[]){ expr_copy(f), expr_new_integer(-k) }, 2)));
                    PUSH_T(term);
                }
            }
        }
    }
    Expr* w = expr_new_function(expr_new_symbol(SYM_Plus), terms, nt);
    free(terms);
    *unk = U; *nu = n;
    #undef PUSH_U
    #undef PUSH_T
    return w;
}

/* Fold the x-free multiplicative content of `e` into the arbitrary constant that
 * will multiply it: divide `e` by the product of its top-level Times factors that
 * are free of x.  A reduction-of-order second solution z2 = z1 ∫1/z1² carries a
 * scalar like 1/(1-2r1) = -2/Sqrt[7]; since it is only ever seen as C[2] z2, that
 * scalar belongs in C[2].  `e` is consumed; result owned. */
static Expr* fold_const_factor(Expr* e, const char* x) {
    if (!head_is(e, SYM_Times)) return e;
    Expr* content = expr_new_integer(1);
    bool any = false;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        Expr* a = e->data.function.args[i];
        Expr* fq = fn2("FreeQ", expr_copy(a), expr_new_symbol(x));
        bool freex = (fq->type == EXPR_SYMBOL && fq->data.symbol.name == intern_symbol("True"));
        expr_free(fq);
        if (freex) { content = T2(content, expr_copy(a)); any = true; }
    }
    if (!any) { expr_free(content); return e; }
    return ds_simplify(T2(e, Powi(content, -1)));   /* e / content */
}

/* Second solution of z'' == r z from a known z1: z2 = z1 ∫ 1/z1^2 dx.  NULL if the
 * integral is not elementary. */
static Expr* second_solution(const Expr* z1, const char* x) {
    Expr* inv = Powi(expr_copy((Expr*)z1), -2);
    g_integrate_quiet++;   /* speculative: a non-elementary ∫ means "decline" */
    Expr* integ = ds_integrate(inv, expr_new_symbol(x));
    g_integrate_quiet--;
    if (ds_has_head(integ, SYM_Integrate)) { expr_free(integ); return NULL; }
    return fold_const_factor(ds_simplify(T2(expr_copy((Expr*)z1), integ)), x);
}

/* z = Exp[∫ω] with the antiderivative guarded (D[∫ω] == ω).  NULL on failure. */
static Expr* exp_integral(const Expr* w, const char* x) {
    g_integrate_quiet++;   /* speculative: a non-elementary ∫ means "decline" */
    Expr* integ = ds_integrate(expr_copy((Expr*)w), expr_new_symbol(x));
    g_integrate_quiet--;
    if (ds_has_head(integ, SYM_Integrate)) { expr_free(integ); return NULL; }
    Expr* back = ds_d(expr_copy(integ), expr_new_symbol(x));
    Expr* chk = Sub(back, expr_copy((Expr*)w));
    bool ok = ds_is_zero(chk);
    expr_free(chk);
    if (!ok) { expr_free(integ); return NULL; }
    return fn1("Exp", integ);
}

/* Apply a `//.` rewrite whose rule(s) are given as a source string (parsed once). */
static Expr* apply_rules(Expr* body, const char* rules_src) {
    Expr* rules = parse_expression(rules_src);
    if (!rules) return body;
    return eval_and_free(expr_new_function(expr_new_symbol("ReplaceRepeated"),
                             (Expr*[]){ body, rules }, 2));
}

/* Does `e` still carry the imaginary unit? (FreeQ[e, Complex] === False) */
static bool has_imaginary(const Expr* e) {
    Expr* fq = fn2("FreeQ", expr_copy((Expr*)e), expr_new_symbol("Complex"));
    bool has = (fq->type == EXPR_SYMBOL && fq->data.symbol.name == intern_symbol("False"));
    expr_free(fq);
    return has;
}

/* Realify a body that may carry the imaginary unit.  First fold the complex forms
 * of Erf/Erfi that the z2 integral of an Exp[poly] solution produces
 * (Erfi[I z] == I Erf[z]) — zero_test cannot cheaply sample Erfi at a complex
 * irrational argument, so leaving them makes verification hang.  Then, if the body
 * is still complex (a genuinely oscillatory solution), ComplexExpand it and fold
 * the Simplify-introduced hyperbolic exponentials back to E (the M8 realifier). */
static Expr* realify(Expr* body) {
    body = apply_rules(body, "{Erfi[I z_] :> I Erf[z], Erf[I z_] :> I Erfi[z]}");
    body = ds_simplify(body);
    if (!has_imaginary(body)) return body;
    Expr* ce = ds_simplify(fn1("ComplexExpand", body));
    return apply_rules(ce, "Cosh[a_] + Sinh[a_] :> E^a");
}

/* Numeric back-substitution check for an algebraic (Case-2) candidate body: the
 * residual of the original equation must be ~0 at several sample points.  The
 * residual is linear in the generated constants C[1],C[2], so we test the two
 * basis solutions independently ((C[1],C[2]) = (1,0) and (0,1)); otherwise
 * N[Abs[...]] would stay symbolic in C[1],C[2] and never evaluate to a number. */
static bool numeric_verify(const DSolveProblem* P, const Expr* body) {
    const char* yname = P->fun_names[0];
    const char* x = P->ind_names[0];
    int maxord = P->max_order[0];
    const long pts[] = { 3, 5, 7, 11 };            /* x = 0.3, 0.5, 0.7, 1.1 */
    const int cval[2][2] = { {1, 0}, {0, 1} };     /* independent basis choices */
    for (size_t e = 0; e < P->neq; e++) {
        Expr* sub = expr_copy(P->eq_residuals[e]);
        for (int k = maxord; k >= 1; k--) {
            Expr* dk = expr_copy((Expr*)body);
            for (int i = 0; i < k; i++) dk = ds_d(dk, expr_new_symbol(x));
            sub = ds_subst(sub, ds_make_funcapp(yname, k, x), dk);
        }
        sub = ds_subst(sub, ds_make_funcapp(yname, 0, x), expr_copy((Expr*)body));
        for (size_t cv = 0; cv < 2; cv++) {
            Expr* subc = ds_subst(expr_copy(sub), ds_const(1), expr_new_integer(cval[cv][0]));
            subc = ds_subst(subc, ds_const(2), expr_new_integer(cval[cv][1]));
            bool all_small = true;
            for (size_t i = 0; i < sizeof(pts)/sizeof(pts[0]); i++) {
                Expr* at = ds_subst(expr_copy(subc), expr_new_symbol(x),
                                    expr_new_real((double)pts[i] / 10.0));
                Expr* val = fn1("N", fn1("Abs", at));
                double v = (val->type == EXPR_REAL) ? val->data.real
                         : (val->type == EXPR_INTEGER) ? (double)val->data.integer : 1e9;
                expr_free(val);
                if (!(v < 1e-6)) { all_small = false; break; }
            }
            expr_free(subc);
            if (!all_small) { expr_free(sub); return false; }
        }
        expr_free(sub);
    }
    return true;
}

/* A reduction-of-order integral in terms of these special functions is a form
 * whose Simplify (and downstream zero_test) hangs on the complex arguments that a
 * growth-at-∞ Kovacic solution produces (Erfi[I z] etc.).  When the second-
 * solution integral lands here we decline the Case-1c candidate and let Case 2 /
 * the series fallback take over, rather than hang. */
static bool has_hang_special(const Expr* e) {
    static const char* fns[] = { "Erf", "Erfi", "Erfc", "ExpIntegralEi",
                                 "ExpIntegralE", "FresnelC", "FresnelS",
                                 "SinIntegral", "CosIntegral", "LogIntegral" };
    for (size_t i = 0; i < sizeof(fns)/sizeof(fns[0]); i++)
        if (ds_contains(e, intern_symbol(fns[i]))) return true;
    return false;
}

/* Numeric linear-independence of two solution candidates: ya/yb is non-constant,
 * i.e. d/dx(ya/yb) is nonzero at some real sample point.  Done numerically (not by
 * Simplify) — a Simplify of a ratio of radical/exp candidates is a multi-second
 * blow-up that the reduction-of-order path deliberately avoids. */
static bool candidates_independent(const Expr* ya, const Expr* yb, const char* x) {
    Expr* ratio = T2(expr_copy((Expr*)ya), Powi(expr_copy((Expr*)yb), -1));
    Expr* dr = ds_d(ratio, expr_new_symbol(x));         /* consumes ratio */
    const double pts[3] = { 0.7, 1.3, 1.9 };
    bool indep = false;
    for (int i = 0; i < 3 && !indep; i++) {
        Expr* at = ds_subst(expr_copy(dr), expr_new_symbol(x), expr_new_real(pts[i]));
        Expr* nv = fn1("N", fn1("Abs", at));
        double v = (nv->type == EXPR_REAL) ? nv->data.real
                 : (nv->type == EXPR_INTEGER) ? (double)nv->data.integer : 0.0;
        expr_free(nv);
        if (v > 1e-9) indep = true;
    }
    expr_free(dr);
    return indep;
}

/* Assemble the general solution C[1] y1 + C[2] y2 from z1,z2 and recovery w. */
static Expr* assemble_general(const Expr* z1, const Expr* z2, const Expr* recovery) {
    Expr* y1 = T2(expr_copy((Expr*)recovery), expr_copy((Expr*)z1));
    Expr* y2 = T2(expr_copy((Expr*)recovery), expr_copy((Expr*)z2));
    Expr* gen = A2(T2(ds_const(1), y1), T2(ds_const(2), y2));
    return realify(gen);
}

/* Case 1 for a polynomial r of even degree 2m: z = P(x) Exp[∫ω_base], where
 * ω_base is the polynomial part of ±√r (degree m) and P is a polynomial of the
 * determined degree d = -[R2]_{m-1}/(2 c_m) solving P'' + 2ω_base P' + R2 P == 0
 * with R2 = ω_base' + ω_base² - r.  Captures apparent singularities (roots of P
 * that are not poles of r), e.g. z''=(x²+3)z -> x Exp[x²/2]. */
static Expr* kovacic_case1_poly_r(const Expr* r, const Expr* recovery,
                                  const char* x, int* counter) {
    int D = degree_in(r, x);
    if (D < 2 || (D % 2) != 0) return NULL;
    int m = D / 2;
    Expr* lead = coeff_k(r, x, D);
    Expr* body = NULL;

    for (int sgn = 1; sgn >= -1 && !body; sgn -= 2) {
        Expr* cm = ds_simplify(T2(expr_new_integer(sgn), fn1("Sqrt", expr_copy(lead))));
        /* ω_base = cm x^m + Σ_{j=0}^{m-1} u_j x^j */
        Expr** U = malloc((size_t)m * sizeof(Expr*));
        Expr** terms = malloc((size_t)(m + 1) * sizeof(Expr*));
        terms[0] = ds_call2(SYM_Times, expr_copy(cm), xpow(x, m));
        for (int j = 0; j < m; j++) {
            char buf[32]; snprintf(buf, sizeof(buf), "DSolve`kv%d", (*counter)++);
            const char* sn = intern_symbol(buf);
            U[j] = expr_new_symbol(sn);
            terms[j + 1] = ds_call2(SYM_Times, expr_new_symbol(sn), xpow(x, j));
        }
        Expr* wbase = expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)(m + 1));
        free(terms);
        Expr* g = Sub(A2(ds_d(expr_copy(wbase), expr_new_symbol(x)), Powi(expr_copy(wbase), 2)),
                      expr_copy((Expr*)r));
        Expr* wb = solve_high_coeffs(g, x, m, 2 * m, U, (size_t)m, wbase);  /* consumes g, wbase */
        for (int j = 0; j < m; j++) expr_free(U[j]);
        free(U);
        if (!wb) { expr_free(cm); continue; }

        Expr* R2 = ds_simplify(Sub(A2(ds_d(expr_copy(wb), expr_new_symbol(x)), Powi(expr_copy(wb), 2)),
                                   expr_copy((Expr*)r)));
        Expr* dcoef = coeff_k(R2, x, m - 1);
        Expr* dexpr = ds_simplify(T2(Neg(dcoef), Powi(T2(expr_new_integer(2), expr_copy(cm)), -1)));
        expr_free(cm);
        int dd = -1;
        if (dexpr->type == EXPR_INTEGER && dexpr->data.integer >= 0
            && dexpr->data.integer < 64)
            dd = (int)dexpr->data.integer;
        expr_free(dexpr);
        if (dd < 0) { expr_free(wb); expr_free(R2); continue; }

        /* P = x^dd + Σ_{i=0}^{dd-1} b_i x^i (monic) */
        Expr** B = malloc((size_t)(dd ? dd : 1) * sizeof(Expr*));
        Expr** pterms = malloc((size_t)(dd + 1) * sizeof(Expr*));
        pterms[0] = xpow(x, dd);
        for (int i = 0; i < dd; i++) {
            char buf[32]; snprintf(buf, sizeof(buf), "DSolve`kv%d", (*counter)++);
            const char* sn = intern_symbol(buf);
            B[i] = expr_new_symbol(sn);
            pterms[i + 1] = ds_call2(SYM_Times, expr_new_symbol(sn), xpow(x, i));
        }
        Expr* Pp = expr_new_function(expr_new_symbol(SYM_Plus), pterms, (size_t)(dd + 1));
        free(pterms);
        /* P'' + 2 ω_base P' + R2 P == 0 */
        Expr* Pode = A2(A2(ds_d(ds_d(expr_copy(Pp), expr_new_symbol(x)), expr_new_symbol(x)),
                           T2(T2(expr_new_integer(2), expr_copy(wb)),
                              ds_d(expr_copy(Pp), expr_new_symbol(x)))),
                        T2(expr_copy(R2), expr_copy(Pp)));
        Expr* Psol = solve_ansatz(Pode, x, B, (size_t)dd, Pp);   /* consumes Pode, Pp */
        for (int i = 0; i < dd; i++) expr_free(B[i]);
        free(B);
        expr_free(R2);
        if (Psol) {
            Expr* e = exp_integral(wb, x);
            if (e) {
                Expr* z1 = ds_simplify(T2(expr_copy(Psol), e));
                Expr* z2 = second_solution(z1, x);
                if (z2) body = assemble_general(z1, z2, recovery);
                if (z2) expr_free(z2);
                expr_free(z1);
            }
            expr_free(Psol);
        }
        expr_free(wb);
    }
    expr_free(lead);
    return body;
}

/* True iff N[e] is (numerically) a nonnegative integer ≤ cap, in which case *out
 * receives it.  Used to test the Kovacic degree bound d = α_∞ - Σα_c cheaply:
 * the α at complex poles are complex algebraic numbers (e.g. built on (-1)^(1/3))
 * whose symbolic Simplify/IntegerQ is minutes-slow across every sign mask, but a
 * numeric evaluation is instant — a wrong candidate is caught later when
 * solve_monic_P finds no polynomial. */
static bool numeric_nonneg_int(const Expr* e, int cap, int* out) {
    Expr* v = fn1("N", expr_copy((Expr*)e));
    double re = 0.0, im = 0.0; bool okr = false;
    if (v->type == EXPR_REAL) { re = v->data.real; okr = true; }
    else if (v->type == EXPR_INTEGER) { re = (double)v->data.integer; okr = true; }
    else if (head_is(v, SYM_Complex) && v->data.function.arg_count == 2) {
        Expr* a = v->data.function.args[0]; Expr* b = v->data.function.args[1];
        if ((a->type == EXPR_REAL || a->type == EXPR_INTEGER) &&
            (b->type == EXPR_REAL || b->type == EXPR_INTEGER)) {
            re = (a->type == EXPR_REAL) ? a->data.real : (double)a->data.integer;
            im = (b->type == EXPR_REAL) ? b->data.real : (double)b->data.integer;
            okr = true;
        }
    }
    expr_free(v);
    if (!okr || fabs(im) > 1e-7) return false;
    long n = (long)(re + (re < 0 ? -0.5 : 0.5));
    if (n < 0 || n > cap || fabs(re - (double)n) > 1e-7) return false;
    *out = (int)n;
    return true;
}

/* Solve for a monic polynomial P of degree dd satisfying the Kovacic Case-1
 * P-equation  P'' + 2 θ P' + R2 P == 0,  R2 = θ' + θ² - r  (undetermined
 * coefficients on P = x^dd + Σ b_i x^i).  R2 is precomputed by the caller (it
 * does not depend on dd).  Returns P (owned, fully determined) or NULL if no such
 * P exists.  `theta`, `R2` borrowed.  The roots of P are the apparent
 * singularities of the solution z1 = P Exp[∫θ] — the zeros the pole-only Riccati
 * ansatz cannot represent. */
static Expr* solve_monic_P(const Expr* theta, const Expr* R2, int dd,
                           const char* x, int* counter) {
    if (dd < 0) return NULL;
    /* P = x^dd + Σ_{i=0}^{dd-1} b_i x^i (monic) */
    Expr** B = malloc((size_t)(dd ? dd : 1) * sizeof(Expr*));
    Expr** pterms = malloc((size_t)(dd + 1) * sizeof(Expr*));
    pterms[0] = xpow(x, dd);
    for (int i = 0; i < dd; i++) {
        char buf[32]; snprintf(buf, sizeof(buf), "DSolve`kv%d", (*counter)++);
        const char* sn = intern_symbol(buf);
        B[i] = expr_new_symbol(sn);
        pterms[i + 1] = ds_call2(SYM_Times, expr_new_symbol(sn), xpow(x, i));
    }
    Expr* Pp = expr_new_function(expr_new_symbol(SYM_Plus), pterms, (size_t)(dd + 1));
    free(pterms);
    /* P'' + 2 θ P' + R2 P */
    Expr* Pode = A2(A2(ds_d(ds_d(expr_copy(Pp), expr_new_symbol(x)), expr_new_symbol(x)),
                       T2(T2(expr_new_integer(2), expr_copy((Expr*)theta)),
                          ds_d(expr_copy(Pp), expr_new_symbol(x)))),
                    T2(expr_copy((Expr*)R2), expr_copy(Pp)));
    Expr* Psol = solve_ansatz(Pode, x, B, (size_t)dd, Pp);   /* consumes Pode, Pp */
    for (int i = 0; i < dd; i++) expr_free(B[i]);
    free(B);
    return Psol;
}

/* SeriesCoefficient[f, {var, c, n}], simplified.  f, c borrowed. */
static Expr* series_coeff(const Expr* f, const char* var, const Expr* c, int n) {
    Expr* spec = expr_new_function(expr_new_symbol(SYM_List),
                     (Expr*[]){ expr_new_symbol(var), expr_copy((Expr*)c), expr_new_integer(n) }, 3);
    Expr* sc = expr_new_function(expr_new_symbol(intern_symbol("SeriesCoefficient")),
                   (Expr*[]){ expr_copy((Expr*)f), spec }, 2);
    return ds_simplify(sc);
}

/* A Laurent coefficient we can use downstream: not an unevaluated SeriesCoefficient
 * / Series and not an infinity/indeterminate. */
static bool kv_clean(const Expr* e) {
    static const char* bad[] = { "SeriesCoefficient", "Series", "SeriesData",
                                 "Infinity", "ComplexInfinity", "DirectedInfinity",
                                 "Indeterminate" };
    for (size_t i = 0; i < sizeof(bad)/sizeof(bad[0]); i++)
        if (ds_contains(e, intern_symbol(bad[i]))) return false;
    return true;
}

/* Classical Kovacic Case-1 data at a pole c of r of order `ord`:
 *   *sqc = [√r]_c, the negative part of the Laurent series of √r keeping the
 *          powers (x-c)^{-i} for 2 ≤ i ≤ ν (ν = ord/2); 0 for ord ∈ {1,2}.
 *   *ap,*am = α_c^± ; *degen marks α^+==α^- with sqc==0 (a redundant sign bit).
 * Returns false when the pole order is incompatible with Case 1 (odd ≥ 3) or a
 * Laurent coefficient did not come out clean. */
static bool kovacic_pole_data(const Expr* r, const char* x, const Expr* c, int ord,
                              Expr** sqc, Expr** ap, Expr** am, bool* degen) {
    *sqc = NULL; *ap = NULL; *am = NULL; *degen = false;
    if (ord <= 0) return false;
    if (ord == 1) {                          /* [√r]_c = 0, α^+ = α^- = 1 */
        *sqc = expr_new_integer(0);
        *ap = expr_new_integer(1); *am = expr_new_integer(1);
        *degen = true;
        return true;
    }
    if (ord == 2) {                          /* [√r]_c = 0, α^± = (1 ± √(1+4b))/2 */
        *sqc = expr_new_integer(0);
        Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
                         (Expr*[]){ expr_new_symbol(x), expr_copy((Expr*)c) }, 2);
        Expr* b = ds_simplify(fn2("Limit",
                      T2(Powi(Sub(expr_new_symbol(x), expr_copy((Expr*)c)), 2), expr_copy((Expr*)r)),
                      rule));
        Expr* disc = fn1("Sqrt", A2(expr_new_integer(1), T2(expr_new_integer(4), b)));
        { int dummy; *degen = numeric_nonneg_int(disc, 0, &dummy); }      /* disc ≈ 0 */
        *ap = T2(A2(expr_new_integer(1), expr_copy(disc)), Powi(expr_new_integer(2), -1));
        *am = T2(Sub(expr_new_integer(1), disc), Powi(expr_new_integer(2), -1));
        return true;
    }
    if (ord % 2 == 1) return false;          /* odd order ≥ 3: only Case 2 */
    int nu = ord / 2;                        /* even order 2ν, ν ≥ 2 */
    Expr* sqr = fn1("Sqrt", expr_copy((Expr*)r));
    Expr* sq = expr_new_integer(0);
    Expr* alead = NULL; bool clean = true;
    for (int j = nu; j >= 2; j--) {          /* keep powers -ν .. -2 (drop the -1 residue) */
        Expr* aj = series_coeff(sqr, x, c, -j);
        if (!kv_clean(aj)) clean = false;
        if (j == nu) alead = expr_copy(aj);
        sq = A2(sq, T2(aj, Powi(Sub(expr_new_symbol(x), expr_copy((Expr*)c)), -j)));
    }
    sq = ds_simplify(sq);
    /* b = coeff of (x-c)^{-(ν+1)} in r - ([√r]_c)² (uses ([√r]_c)², NOT the raw √r
     * Laurent series — they diverge for ν ≥ 2). */
    Expr* rmsq = ds_simplify(Sub(expr_copy((Expr*)r), Powi(expr_copy(sq), 2)));
    Expr* b = series_coeff(rmsq, x, c, -(nu + 1));
    expr_free(rmsq);
    bool okc = clean && kv_clean(b) && alead && kv_clean(alead) && !ds_is_zero(alead);
    if (!okc) { expr_free(sqr); expr_free(sq); if (alead) expr_free(alead); expr_free(b); return false; }
    Expr* ba = T2(expr_copy(b), Powi(expr_copy(alead), -1));            /* b/a */
    *ap = ds_simplify(T2(A2(expr_copy(ba), expr_new_integer(nu)), Powi(expr_new_integer(2), -1)));
    *am = ds_simplify(T2(A2(Neg(ba), expr_new_integer(nu)), Powi(expr_new_integer(2), -1)));
    *sqc = sq;
    expr_free(sqr); expr_free(alead); expr_free(b);
    return true;
}

/* Classical Kovacic Case-1 data at infinity for O(∞) = degd - degn:
 *   O(∞) ≥ 3  : [√r]_∞ = 0, α^+ = 0, α^- = 1.
 *   O(∞) == 2 : [√r]_∞ = 0, α^± = (1 ± √(1+4b))/2, b = lim x²r.
 *   O(∞) ≤ 0 even (growth): [√r]_∞ = the degree-ν′ polynomial part of √r at ∞
 *       (ν′ = -O(∞)/2), α^± = (±b/a − ν′)/2 with b = coeff of x^{ν′-1} in r-([√r]_∞)².
 * The polynomial part and the "b via ([√r]_∞)²" are computed at y=0 after x→1/y,
 * since SeriesCoefficient does not evaluate directly at Infinity.  Returns false
 * for O(∞)==1 or negative-odd (incompatible with Case 1), or an unclean coeff. */
static bool kovacic_inf_data(const Expr* r, const char* x, int oinf,
                             Expr** sqinf, Expr** ap, Expr** am) {
    *sqinf = NULL; *ap = NULL; *am = NULL;
    if (oinf >= 3) {
        *sqinf = expr_new_integer(0);
        *ap = expr_new_integer(0); *am = expr_new_integer(1);
        return true;
    }
    if (oinf == 2) {
        *sqinf = expr_new_integer(0);
        Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
                        (Expr*[]){ expr_new_symbol(x),
                                   expr_new_symbol(intern_symbol("Infinity")) }, 2);
        Expr* binf = ds_simplify(fn2("Limit",
                         T2(Powi(expr_new_symbol(x), 2), expr_copy((Expr*)r)), rule));
        Expr* di = ds_simplify(fn1("Sqrt", A2(expr_new_integer(1), T2(expr_new_integer(4), binf))));
        *ap = ds_simplify(T2(A2(expr_new_integer(1), expr_copy(di)), Powi(expr_new_integer(2), -1)));
        *am = ds_simplify(T2(Sub(expr_new_integer(1), di), Powi(expr_new_integer(2), -1)));
        return true;
    }
    if (oinf == 1 || (oinf % 2 != 0)) return false;   /* odd order at ∞: not Case 1 */
    int nup = -oinf / 2;                              /* ν′ ≥ 0 */
    const char* y = intern_symbol("DSolve`ky");
    Expr* zero = expr_new_integer(0);
    Expr* g = ds_simplify(ds_subst(expr_copy((Expr*)r), expr_new_symbol(x),
                  Powi(expr_new_symbol(y), -1)));                     /* r(1/y) */
    Expr* sqg = fn1("Sqrt", expr_copy(g));
    Expr* sq = expr_new_integer(0);
    Expr* alead = NULL; bool clean = true;
    for (int j = nup; j >= 0; j--) {                                  /* powers x^ν′..x^0 */
        Expr* aj = series_coeff(sqg, y, zero, -j);                    /* coeff of y^{-j} = x^j */
        if (!kv_clean(aj)) clean = false;
        if (j == nup) alead = expr_copy(aj);
        sq = A2(sq, T2(aj, Powi(expr_new_symbol(x), j)));
    }
    sq = ds_simplify(sq);
    Expr* rmsq = ds_simplify(Sub(expr_copy((Expr*)r), Powi(expr_copy(sq), 2)));
    Expr* rmsqy = ds_simplify(ds_subst(rmsq, expr_new_symbol(x), Powi(expr_new_symbol(y), -1)));
    Expr* b = series_coeff(rmsqy, y, zero, -(nup - 1));               /* coeff of x^{ν′-1} */
    expr_free(rmsqy);
    bool okc = clean && kv_clean(b) && alead && kv_clean(alead) && !ds_is_zero(alead);
    if (!okc) { expr_free(g); expr_free(sqg); expr_free(sq); if (alead) expr_free(alead); expr_free(b); expr_free(zero); return false; }
    Expr* ba = T2(expr_copy(b), Powi(expr_copy(alead), -1));          /* b/a */
    *ap = ds_simplify(T2(Sub(expr_copy(ba), expr_new_integer(nup)), Powi(expr_new_integer(2), -1)));
    *am = ds_simplify(T2(Sub(Neg(ba), expr_new_integer(nup)), Powi(expr_new_integer(2), -1)));
    *sqinf = sq;
    expr_free(g); expr_free(sqg); expr_free(alead); expr_free(b); expr_free(zero);
    return true;
}

/* General classical Kovacic Case 1 for a rational r = rn/rd with genuine poles
 * (degd >= 1).  The pole-only Riccati ansatz (Case 1 above) misses solutions
 * whose z1 has zeros away from the poles of r — an *apparent singularity*, the
 * P'/P term of ω.  Here we build ω = θ + P'/P with:
 *   ω  = s(∞)[√r]_∞ + Σ_c ( s(c)[√r]_c + α_c^{s(c)}/(x - c) )    and
 *   P  = a monic polynomial of degree d (its roots are the apparent
 *        singularities), found by solve_monic_P.
 * The per-pole [√r]_c / α_c^± come from kovacic_pole_data and the α at infinity
 * plus the growing-r polynomial part [√r]_∞ from kovacic_inf_data — the FULL
 * classical Case 1 (Kovacic 1986), so this now covers even-order poles ≥ 4 and r
 * that grows at ∞ (O(∞) = degd-degn ≤ 0), not just the δ ≥ 2 orthogonal-polynomial
 * slice.  We enumerate the ∞ sign and the ± sign per pole (complex poles included),
 * and take the degree bound d = α_∞^{s(∞)} − Σ_c α_c^{s(c)} where it is a
 * nonnegative integer.  A valid (ω,d) gives z1 = P Exp[∫ω]: the algebraic identity
 * P'' + 2ω P' + (ω'+ω²-r)P == 0 makes z1'' == r z1 hold *exactly*; dsolve_run
 * verifies the assembled body symbolically as the backstop.  Examples: Legendre /
 * Chebyshev / Gegenbauer / Jacobi (order-2 poles), and Kovacic's own paper example
 * y'' = ((4x⁶-8x⁵+12x⁴+4x³+7x²-20x+4)/(4x⁴)) y (order-4 pole + x² growth). */
#define KOV_C1G_DMAX 10     /* max apparent-singularity polynomial degree searched */
#define KOV_C1G_MAXPOLES 6  /* bound the 2^k sign enumeration                      */
#define KOV_DENOM_DEGMAX 8  /* bound deg(denominator of r) before any solve/roots  */
static Expr* kovacic_case1_general(const Expr* r, const Expr* rd,
                                   int degn, int degd,
                                   const Expr* recovery, const char* x, int* counter) {
    if (degd < 1) return NULL;              /* polynomial r: kovacic_case1_poly_r */
    int oinf = degd - degn;                 /* order of r at infinity */

    /* [√r]_∞ and α_∞^± (declines an infinity order incompatible with Case 1). */
    Expr* sqinf; Expr* ainf_p; Expr* ainf_m;
    if (!kovacic_inf_data(r, x, oinf, &sqinf, &ainf_p, &ainf_m)) return NULL;

    DSolveRoots pr;
    if (!dsolve_analyze_roots(rd, x, degd, &pr)) {
        expr_free(sqinf); expr_free(ainf_p); expr_free(ainf_m); return NULL;
    }
    size_t k = pr.ndist;
    if (k == 0 || k > KOV_C1G_MAXPOLES) {
        dsolve_roots_free(&pr); expr_free(sqinf); expr_free(ainf_p); expr_free(ainf_m); return NULL;
    }

    /* per-pole [√r]_c and α_c^± (kovacic_pole_data); `degen[i]` marks a redundant
     * sign bit (α^+==α^- with sqc==0).  A pole order incompatible with Case 1
     * (odd ≥ 3) or an unclean Laurent coefficient makes the whole case decline. */
    Expr** sqc = malloc(k * sizeof(Expr*));
    Expr** ap  = malloc(k * sizeof(Expr*));
    Expr** am  = malloc(k * sizeof(Expr*));
    bool*  degen = malloc(k * sizeof(bool));
    for (size_t i = 0; i < k; i++) { sqc[i] = NULL; ap[i] = NULL; am[i] = NULL; degen[i] = false; }
    bool ok = true;
    for (size_t i = 0; i < k && ok; i++)
        if (!kovacic_pole_data(r, x, pr.roots[i], pr.mult[i], &sqc[i], &ap[i], &am[i], &degen[i]))
            ok = false;

    /* w² = Exp[-∫P] = recovery², the reduction-of-order weight (independent of ω). */
    Expr* w2 = ds_simplify(Powi(expr_copy((Expr*)recovery), 2));

    Expr* body = NULL;
    if (ok) {
        size_t ncomb = (size_t)1 << k;
        /* Pass 1: collect a first-solution candidate y1 = w·P·Exp[∫ω] for each
         * (∞ sign, pole sign mask).  Different sign families give different (but
         * equally valid) members of the solution space; the second-solution
         * integral below is the expensive step, done only for the *cleanest* y1
         * (smallest LeafCount — the orthogonal-polynomial member, vs a
         * fractional-power sibling). */
        Expr** cand = malloc(2 * (ncomb ? ncomb : 1) * sizeof(Expr*));
        long*  score = malloc(2 * (ncomb ? ncomb : 1) * sizeof(long));
        size_t ncand = 0;
        for (int sinf = 0; sinf < 2; sinf++) {          /* ∞ sign: 1 = '+', 0 = '-' */
            Expr* ainf = sinf ? ainf_p : ainf_m;
            for (size_t mask = 0; mask < ncomb; mask++) {
                bool redundant = false;                 /* skip fixed bits of degenerate poles */
                for (size_t i = 0; i < k; i++)
                    if (degen[i] && ((mask >> i) & 1)) { redundant = true; break; }
                if (redundant) continue;
                /* d = α_∞^{s(∞)} - Σ_c α_c^{s(c)} — tested *numerically* for
                 * nonnegative-integrality (numeric_nonneg_int) before the expensive
                 * ω/exp_integral, so a family with no integer degree bound is skipped
                 * instantly. */
                Expr* suma = expr_new_integer(0);
                for (size_t i = 0; i < k; i++)
                    suma = A2(suma, expr_copy(((mask >> i) & 1) ? ap[i] : am[i]));
                Expr* de = Sub(expr_copy(ainf), suma);
                int dv;
                bool okd = numeric_nonneg_int(de, KOV_C1G_DMAX, &dv);
                expr_free(de);
                if (!okd) continue;                     /* no valid degree: skip */

                /* ω = s(∞)[√r]_∞ + Σ_i ( s_i[√r]_c + α_i^{s_i}/(x - c_i) ) */
                Expr* theta = sinf ? expr_copy(sqinf) : Neg(expr_copy(sqinf));
                for (size_t i = 0; i < k; i++) {
                    bool si = (mask >> i) & 1;
                    Expr* sct = si ? expr_copy(sqc[i]) : Neg(expr_copy(sqc[i]));
                    Expr* a   = si ? ap[i] : am[i];
                    theta = A2(theta, A2(sct, T2(expr_copy(a),
                                 Powi(Sub(expr_new_symbol(x), expr_copy(pr.roots[i])), -1))));
                }
                theta = ds_simplify(theta);
                Expr* e = exp_integral(theta, x);       /* Exp[∫ω]; borrows theta */
                if (!e) { expr_free(theta); continue; }
                /* R2 = ω' + ω² - r */
                Expr* R2 = ds_simplify(Sub(A2(ds_d(expr_copy(theta), expr_new_symbol(x)),
                                              Powi(expr_copy(theta), 2)),
                                           expr_copy((Expr*)r)));
                Expr* Psol = solve_monic_P(theta, R2, dv, x, counter);
                if (Psol) {
                    Expr* y1 = ds_simplify(T2(T2(expr_copy((Expr*)recovery), expr_copy(Psol)),
                                              expr_copy(e)));
                    expr_free(Psol);
                    Expr* lc = fn1("LeafCount", expr_copy(y1));
                    cand[ncand] = y1;
                    score[ncand] = (lc->type == EXPR_INTEGER) ? lc->data.integer : 1L << 30;
                    expr_free(lc);
                    ncand++;
                }
                expr_free(R2);
                expr_free(e);
                expr_free(theta);
            }
        }

        /* Pass 2: sort candidates cleanest-first.  Prefer TWO INDEPENDENT
         * Liouvillian solutions coming from distinct sign families — no integral
         * needed — which is what makes Kovacic's own example solve in closed form
         * (its d=2 and d=0 families give independent y1,y2, while its reduction-of-
         * order second solution is non-elementary).  Only if no independent pair
         * exists do we fall back to reduction of order y2 = y1 ∫ w²/y1² from the
         * cleanest candidate (e.g. Legendre/Chebyshev, whose second solution is the
         * elementary log integral). */
        for (size_t i = 0; i < ncand; i++)              /* selection sort by LeafCount */
            for (size_t j = i + 1; j < ncand; j++)
                if (score[j] < score[i]) {
                    Expr* te = cand[i]; cand[i] = cand[j]; cand[j] = te;
                    long ts = score[i]; score[i] = score[j]; score[j] = ts;
                }
        for (size_t a = 0; a < ncand && !body; a++) {
            for (size_t b = a + 1; b < ncand && !body; b++) {
                if (candidates_independent(cand[a], cand[b], x))
                    body = A2(T2(ds_const(1), expr_copy(cand[a])),
                              T2(ds_const(2), expr_copy(cand[b])));
            }
        }
        for (size_t pick = 0; pick < ncand && !body; pick++) {   /* reduction of order */
            Expr* y1 = cand[pick];
            Expr* integrand = ds_simplify(T2(expr_copy(w2), Powi(expr_copy(y1), -2)));
            g_integrate_quiet++;   /* speculative: non-elementary ∫ => this y1 declines */
            Expr* integ = ds_integrate(integrand, expr_new_symbol(x));  /* consumes integrand */
            g_integrate_quiet--;
            if (!ds_has_head(integ, SYM_Integrate) && !has_hang_special(integ)) {   /* ∫ elementary & clean */
                Expr* y2 = ds_simplify(T2(expr_copy(y1), integ));       /* consumes integ */
                /* Assemble C[1] y1 + C[2] y2 with y1, y2 already individually
                 * Simplify'd.  Deliberately NOT realify()'d: a final Simplify of
                 * the whole sum tries to cross-factor over a radical second
                 * solution (e.g. the (√(x²-1)-x)^6 form the integrator emits for
                 * Gegenbauer), a multi-second blow-up for no gain. */
                body = A2(T2(ds_const(1), expr_copy(y1)), T2(ds_const(2), y2));
            } else {
                expr_free(integ);
            }
        }
        for (size_t i = 0; i < ncand; i++) expr_free(cand[i]);
        free(cand); free(score);
    }
    expr_free(w2); expr_free(sqinf); expr_free(ainf_p); expr_free(ainf_m);

    for (size_t i = 0; i < k; i++) {
        if (sqc[i]) expr_free(sqc[i]);
        if (ap[i]) expr_free(ap[i]);
        if (am[i]) expr_free(am[i]);
    }
    free(sqc); free(ap); free(am); free(degen);
    dsolve_roots_free(&pr);
    return body;
}

Expr** dsolve_kovacic_try(DSolveProblem* P, size_t* nbranch) {
    Expr* Pc; Expr* Qc;
    if (!dsolve_second_order_PQ(P, &Pc, &Qc)) return NULL;
    const char* x = P->ind_names[0];

    Expr* recovery = NULL;
    Expr* r = dsolve_normal_form(Pc, Qc, x, &recovery);
    expr_free(Pc); expr_free(Qc);
    if (!recovery) { expr_free(r); return NULL; }   /* cannot recover y = w z */

    /* r = rn/rd, factor the denominator for the pole structure */
    Expr* rt = fn1("Together", expr_copy(r));
    Expr* rn = fn1("Numerator", expr_copy(rt));
    Expr* rd = fn1("Denominator", rt);          /* consumes rt */
    /* Kovacic requires r ∈ C(x): decline non-rational r (Frobenius handles it). */
    bool rational;
    {
        const char* T = intern_symbol("True");
        Expr* q1 = fn2("PolynomialQ", expr_copy(rn), expr_new_symbol(x));
        Expr* q2 = fn2("PolynomialQ", expr_copy(rd), expr_new_symbol(x));
        rational = (q1->type == EXPR_SYMBOL && q1->data.symbol.name == T)
                && (q2->type == EXPR_SYMBOL && q2->data.symbol.name == T);
        expr_free(q1); expr_free(q2);
    }
    if (!rational) { expr_free(rn); expr_free(rd); expr_free(r); expr_free(recovery); return NULL; }
    int degn = degree_in(rn, x);
    int degd = degree_in(rd, x);
    Expr* factors = fn1("FactorList", expr_copy(rd));
    expr_free(rn);   /* rd kept alive for the general Case-1 pole enumeration */

    /* ---- Early complexity gate (before any ds_solve / analyze_roots) ----
     * Both the Case-1 Riccati ansatz (build_riccati_ansatz + solve_ansatz) and
     * the Case-1c pole enumeration (dsolve_analyze_roots) are driven by the pole
     * structure of r.  On a high-degree denominator — e.g. the normal form of
     * y'' + ((x^10-1)/(x^12+1)) y == 0, with rd = x^12+1 — clearing denominators
     * yields a ~degree-24 nonlinear system and dsolve_analyze_roots must solve a
     * degree-12 polynomial symbolically; either can run effectively unbounded.
     * KOV_C1G_MAXPOLES already bounds the sign enumeration, but only *inside*
     * Case 1c, after that expensive work.  Apply the same philosophy up front,
     * cheaply, from the (polynomial-time) FactorList: if r has too many distinct
     * poles or too high a denominator degree, decline to the Frobenius fallback
     * (a fast Taylor series about the ordinary point) rather than hang. */
    {
        int npole = 0;
        if (factors && head_is(factors, SYM_List)) {
            for (size_t fi = 0; fi < factors->data.function.arg_count; fi++) {
                Expr* pair = factors->data.function.args[fi];
                if (head_is(pair, SYM_List) && pair->data.function.arg_count == 2 &&
                    degree_in(pair->data.function.args[0], x) >= 1) npole++;
            }
        }
        if (npole > KOV_C1G_MAXPOLES || degd > KOV_DENOM_DEGMAX) {
            expr_free(r); expr_free(rd); expr_free(factors); expr_free(recovery);
            return NULL;
        }
    }

    int ddiff = degn - degd;
    int poly_deg = (ddiff > 0) ? (ddiff + 1) / 2 : 0;
    if (poly_deg > 8) poly_deg = 8;

    Expr* body = NULL;

    /* ---- Case 1: ω' + ω² == r, ω ∈ C(x) ---- */
    {
        int counter = 1;
        Expr** unk; size_t nu;
        Expr* w = build_riccati_ansatz(x, poly_deg, factors, &counter, &unk, &nu);
        Expr* eq = Sub(A2(ds_d(expr_copy(w), expr_new_symbol(x)), Powi(expr_copy(w), 2)),
                       expr_copy(r));
        Expr* ws = solve_ansatz(eq, x, unk, nu, w);   /* consumes eq and w */
        for (size_t i = 0; i < nu; i++) expr_free(unk[i]);
        free(unk);
        if (ws) {
            Expr* z1 = exp_integral(ws, x);
            if (z1) {
                Expr* z2 = second_solution(z1, x);
                if (z2) body = assemble_general(z1, z2, recovery);
                if (z2) expr_free(z2);
                expr_free(z1);
            }
            expr_free(ws);
        }
    }

    /* ---- Case 1b: polynomial r, z = P Exp[∫ω_base] (apparent singularities) ---- */
    if (!body && degd == 0) {
        int counter = 1;
        body = kovacic_case1_poly_r(r, recovery, x, &counter);
    }

    /* ---- Case 1c: rational r with genuine poles, apparent singularities via a
     *      monic polynomial P over the local pole exponents (the classical Case-1
     *      completion the pole-only Riccati ansatz misses — e.g. Legendre,
     *      Chebyshev, Gegenbauer, Jacobi).  Runs before the heavier Case 2. ---- */
    if (!body && degd >= 1) {
        int counter = 1;
        body = kovacic_case1_general(r, rd, degn, degd, recovery, x, &counter);
    }

    /* Guard: Case 2's rational-σ search builds an undetermined-coefficient system
     * over the irreducible factors of the denominator and hands it to Solve.  When
     * a factor has degree ≥ 2 (a complex-conjugate pole pair), that system is
     * large and coupled and Solve can run for minutes — while the complex-pole
     * Liouvillian solutions it might find are already covered by Case 1c.  So only
     * attempt Case 2 when every pole is simple (all factors linear); otherwise
     * decline to the series fallback rather than hang. */
    bool case2_ok = true;
    if (!body && factors && head_is(factors, SYM_List)) {
        for (size_t fi = 0; fi < factors->data.function.arg_count; fi++) {
            Expr* pair = factors->data.function.args[fi];
            if (head_is(pair, SYM_List) && pair->data.function.arg_count == 2 &&
                degree_in(pair->data.function.args[0], x) >= 2) { case2_ok = false; break; }
        }
    }

    /* ---- Case 2: σ ∈ C(x), D' + 2σD == 0, D = 4r - 2σ' - σ² ---- */
    if (!body && case2_ok) {
        int counter = 1;
        Expr** unk; size_t nu;
        Expr* sig = build_riccati_ansatz(x, poly_deg, factors, &counter, &unk, &nu);
        /* D = 4r - 2 σ' - σ² */
        Expr* Dd = Sub(Sub(T2(expr_new_integer(4), expr_copy(r)),
                           T2(expr_new_integer(2), ds_d(expr_copy(sig), expr_new_symbol(x)))),
                       Powi(expr_copy(sig), 2));
        /* eq: D' + 2 σ D == 0 */
        Expr* eq = A2(ds_d(expr_copy(Dd), expr_new_symbol(x)),
                      T2(T2(expr_new_integer(2), expr_copy(sig)), expr_copy(Dd)));
        Expr* sigs = solve_ansatz(eq, x, unk, nu, sig);   /* consumes eq and sig */
        for (size_t i = 0; i < nu; i++) expr_free(unk[i]);
        free(unk);
        if (sigs) {
            /* ω = (σ + Sqrt[D]) / 2 with D re-evaluated at the solved σ */
            Expr* Ds = Sub(Sub(T2(expr_new_integer(4), expr_copy(r)),
                               T2(expr_new_integer(2), ds_d(expr_copy(sigs), expr_new_symbol(x)))),
                           Powi(expr_copy(sigs), 2));
            Ds = ds_simplify(Ds);
            Expr* rootD = fn1("Sqrt", Ds);
            Expr* w1 = T2(A2(expr_copy(sigs), expr_copy(rootD)), Powi(expr_new_integer(2), -1));
            Expr* w2 = T2(Sub(expr_copy(sigs), expr_copy(rootD)), Powi(expr_new_integer(2), -1));
            expr_free(rootD);
            Expr* z1 = exp_integral(w1, x);
            Expr* z2 = exp_integral(w2, x);
            if (z1 && z2) {
                Expr* cand = assemble_general(z1, z2, recovery);
                if (numeric_verify(P, cand)) body = cand; else expr_free(cand);
            }
            if (z1) expr_free(z1);
            if (z2) expr_free(z2);
            expr_free(w1); expr_free(w2);
            expr_free(sigs);
        }
    }

    expr_free(r); expr_free(rd); expr_free(factors); expr_free(recovery);   /* rt was consumed by Denominator */
    if (!body) return NULL;
    Expr** out = malloc(sizeof(Expr*));
    out[0] = body;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_kovacic(Expr* res) {
    return dsolve_method_builtin(res, dsolve_kovacic_try);
}

void dsolve_kovacic_init(void) {
    symtab_add_builtin("DSolve`Kovacic", builtin_dsolve_kovacic);
    symtab_get_def("DSolve`Kovacic")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`Kovacic",
        "DSolve`Kovacic[eqn, y, x] finds Liouvillian solutions of a second-order "
        "linear ODE y'' + P y' + Q y == 0 by reducing to z'' == r z and searching "
        "for the logarithmic derivative of a solution: Case 1 (rational omega, z = "
        "P Exp[Integrate[theta]], including the apparent-singularity monic P over the "
        "local pole exponents that yields the elementary Legendre/Chebyshev/Gegenbauer "
        "family) and Case 2 (degree-2 algebraic). Declines otherwise.");
}
