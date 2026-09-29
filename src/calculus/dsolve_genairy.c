/*
 * dsolve_genairy.c — DSolve`GeneralizedAiry.
 *
 * The n-th order PURE-POWER potential, the higher-order analogue of the Airy
 * equation:
 *
 *      u^(n) == A x^m u          (n >= 3,  p := m + n != 0)
 *
 * has the fundamental set, for j = 0 .. n-1,
 *
 *      u_j = x^j 0F_{n-1}( ; { 1 + (j-i)/p : i = 0..n-1, i != j } ; A x^p / p^n ).
 *
 * Derivation: substituting u = x^j Sum_k c_k x^(p k) balances x^(j+pk-n) against
 * the right-hand side (p - m = n), giving
 *      c_k / c_{k-1} = (A/p^n) / Prod_{i=0}^{n-1} (k + (j-i)/p),
 * whose i == j factor is exactly k — the 0F_{n-1} denominator's k! — while the
 * other n-1 factors are the lower parameters (shifted by one).  For n = 2, m = 1
 * this is Airy; for n = 3, m = 1 it is 0F2(;{1/2,3/4}; x^4/64) and its partners.
 *
 * The equation reaches this shape through a DEPRESSION (gauge) pre-pass.  For a
 * monic linear operator  L = D^n + c_{n-1} D^(n-1) + ... + c_0,  the substitution
 * y = w u with the logarithmic derivative v := -c_{n-1}/n,  w = Exp[Integrate[v]],
 * gives  L[w u]/w = Sum_j d_j u^(j)  with
 *
 *      d_j = Sum_{k=j}^{n} c_k Binomial[k,j] W_{k-j},   W_0 = 1, W_{i+1} = W_i' + v W_i
 *
 * (the W_i = w^(i)/w are the Bell polynomials in v, so no Exp ever has to cancel).
 * By construction d_n == 1 and d_{n-1} == 0; when every remaining intermediate
 * coefficient vanishes the depressed equation is u^(n) == -d_0 u, and the pure
 * power is read off -d_0.  That gauge is what maps the corpus operator
 * x y''' + 3 y'' - A x^2 y == 0  (w = 1/x)  onto  u''' == A x u.
 *
 * Correctness.  HypergeometricPFQ numericizes and carries a z-derivative rule, so
 * every emitted branch is checked by an in-method numeric back-substitution on the
 * ORIGINAL residual (ga_num_ok) BEFORE it is returned, with free symbolic
 * parameters instantiated at distinct generic reals — the corpus shape is
 * y''' == A x^b y with symbolic A, b, whose residual zero_test cannot decide.  A
 * wrong read is therefore a clean decline (to Frobenius), never a wrong answer.
 *
 * Declines (all fall through to the series fallback): p == 0 (the Euler case);
 * a lower parameter that is a non-positive integer (exponents differing by a
 * multiple of p — the logarithmic Frobenius case, e.g. x w''' + w == 0 and
 * x^2 y'''' == A y); a non-power or undefined potential (y''' == Sin[x] y); a
 * non-elementary depression integral; forcing whose variation-of-parameters
 * integrals do not close.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>
#include <math.h>

#define GA_MAX_ORDER 6     /* beyond this the depression/Bell expansion is not worth it */

/* ---- small builders (args consumed, result owned) ---- */
static Expr* ga_mul(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Times, a, b)); }
static Expr* ga_add(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Plus,  a, b)); }
static Expr* ga_sub(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Subtract, a, b)); }
static Expr* ga_pow(Expr* b, Expr* e) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ b, e }, 2));
}
static Expr* ga_powi(Expr* b, int e) { return ga_pow(b, expr_new_integer(e)); }
static Expr* ga_div(Expr* a, Expr* b) { return ds_simplify(ga_mul(a, ga_powi(b, -1))); }
static Expr* ga_fn1(const char* h, Expr* a) { return eval_and_free(ds_call1(h, a)); }
static Expr* ga_fn2(const char* h, Expr* a, Expr* b) { return eval_and_free(ds_call2(h, a, b)); }

static long ga_binom(int nn, int kk) {
    if (kk < 0 || kk > nn) return 0;
    long r = 1;
    for (int i = 0; i < kk; i++) r = r * (nn - i) / (i + 1);
    return r;
}

/* ---- numeric self-verify ------------------------------------------------- *
 * Back-substitute `general` (with all C[k] and any free symbolic parameters
 * instantiated at distinct generic reals) into the ORIGINAL residual and sample.
 * Modelled on dsolve_specialform.c's sf_num_ok, lifted to arbitrary order.     */

static void ga_collect_params(const Expr* e, const char** names, int* n, int cap) {
    if (!e || *n >= cap) return;
    if (e->type == EXPR_SYMBOL) {
        const char* nm = e->data.symbol.name;
        for (int i = 0; i < *n; i++) if (names[i] == nm) return;
        names[(*n)++] = nm; return;
    }
    if (e->type == EXPR_FUNCTION) {
        /* skip a function HEAD (Sin, C, ...): only argument-position symbols are
         * parameters, as in sf_num_ok. */
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            ga_collect_params(e->data.function.args[i], names, n, cap);
    }
}

static double ga_abs_at(const Expr* R, const char* xv, double xval) {
    Expr* e = ds_subst(expr_copy((Expr*)R), expr_new_symbol(xv), expr_new_real(xval));
    e = eval_and_free(ds_call1("Abs", ga_fn2("N", e, expr_new_integer(30))));
    double m = (e && e->type == EXPR_REAL) ? e->data.real
             : (e && e->type == EXPR_INTEGER) ? (double)e->data.integer : NAN;
    expr_free(e);
    return m;
}

static bool ga_num_ok(const DSolveProblem* P, const Expr* general,
                      const char* xv, const char* yname, int order) {
    Expr* R = expr_copy(P->eq_residuals[0]);
    for (int k = order; k >= 0; k--) {
        Expr* d = expr_copy((Expr*)general);
        for (int i = 0; i < k; i++) d = ds_d(d, expr_new_symbol(xv));
        R = ds_subst(R, ds_make_funcapp(yname, k, xv), d);
    }
    for (int k = 1; k <= order + 1; k++)
        R = ds_subst(R, ds_const(k), expr_new_real(0.31 + 0.17 * (double)k));
    {
        const char* skip[] = { xv, intern_symbol("E"), intern_symbol("Pi"),
            intern_symbol("I"), intern_symbol("C"), intern_symbol("EulerGamma"),
            intern_symbol("Degree"), intern_symbol("GoldenRatio"),
            intern_symbol("Catalan"), intern_symbol("Infinity") };
        const int nskip = (int)(sizeof(skip) / sizeof(skip[0]));
        const char* syms[64]; int ns = 0;
        ga_collect_params(R, syms, &ns, 64);
        int pi = 0;
        for (int i = 0; i < ns; i++) {
            bool sk = false;
            for (int j = 0; j < nskip; j++) if (syms[i] == skip[j]) { sk = true; break; }
            if (sk) continue;
            double v = 0.37 + 0.11 * (double)pi; pi++;
            R = ds_subst(R, expr_new_symbol(syms[i]), expr_new_real(v));
        }
    }
    /* Sample away from x = 0 (a regular singular point of every member with
     * m < 0) and below the radius where the pFq argument gets large. */
    const double xs[] = { 0.6, 0.9, 1.2, 1.5, 1.8 };
    int small = 0, big = 0;
    for (int i = 0; i < 5; i++) {
        double m = ga_abs_at(R, xv, xs[i]);
        if (isnan(m) || !isfinite(m)) continue;
        if (m < 1e-8) small++; else if (m > 1e-5) big++;
    }
    expr_free(R);
    return small >= 3 && big == 0;
}

/* ---- the recogniser ------------------------------------------------------ */

/* True when `e` is a decidable non-positive integer (0, -1, -2, ...), which makes
 * a 0F_{n-1} lower parameter singular. */
static bool ga_bad_lower(const Expr* e) {
    Expr* v = ds_simplify(expr_copy((Expr*)e));
    /* A lower parameter is 1 + (j-i)/p with |j-i| < n, so it is a small rational; an
     * int64 test covers every representable case (a BIGINT would need |p| < 1/2^63). */
    bool bad = (v->type == EXPR_INTEGER) && (v->data.integer <= 0);
    expr_free(v);
    return bad;
}

/* x^j 0F_{n-1}( ; {1 + (j-i)/p : i != j} ; A x^p / p^n ).  p, A borrowed. */
static Expr* ga_basis_member(int j, int n, const Expr* p, const Expr* A, const char* xv,
                             bool* degenerate) {
    Expr** lower = malloc((size_t)(n - 1) * sizeof(Expr*));
    int nl = 0;
    for (int i = 0; i < n; i++) {
        if (i == j) continue;
        Expr* b = ds_simplify(ga_add(expr_new_integer(1),
                      ga_div(expr_new_integer(j - i), expr_copy((Expr*)p))));
        if (ga_bad_lower(b)) { *degenerate = true; expr_free(b); break; }
        lower[nl++] = b;
    }
    if (*degenerate) {
        for (int i = 0; i < nl; i++) expr_free(lower[i]);
        free(lower);
        return NULL;
    }
    Expr* lowlist = expr_new_function(expr_new_symbol(SYM_List), lower, (size_t)nl);
    free(lower);
    Expr* z = ds_simplify(ga_div(ga_mul(expr_copy((Expr*)A),
                                        ga_pow(expr_new_symbol(xv), expr_copy((Expr*)p))),
                                 ga_pow(expr_copy((Expr*)p), expr_new_integer(n))));
    Expr* pfq = eval_and_free(expr_new_function(expr_new_symbol("HypergeometricPFQ"),
                    (Expr*[]){ expr_new_function(expr_new_symbol(SYM_List), NULL, 0),
                               lowlist, z }, 3));
    if (j == 0) return pfq;
    return ga_mul(ga_powi(expr_new_symbol(xv), j), pfq);
}

Expr** dsolve_genairy_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1 || P->is_pde) return NULL;
    if (P->max_order[0] < 3 || P->max_order[0] > GA_MAX_ORDER) return NULL;
    const char* xv = P->ind_names[0];

    Expr** a = NULL; Expr* g = NULL; int n = 0;
    if (!dsolve_linear_coeffs(P, &a, &g, &n)) return NULL;
    if (n < 3 || n > GA_MAX_ORDER || ds_is_zero(a[n])) {
        for (int k = 0; k <= n; k++) expr_free(a[k]);
        free(a); if (g) expr_free(g);
        return NULL;
    }

    /* monic coefficients c_k = a_k / a_n, forcing gg = g / a_n */
    Expr** c = malloc((size_t)(n + 1) * sizeof(Expr*));
    for (int k = 0; k <= n; k++) c[k] = ds_simplify(ga_div(expr_copy(a[k]), expr_copy(a[n])));
    Expr* gg = ds_is_zero(g) ? NULL : ds_simplify(ga_div(expr_copy(g), expr_copy(a[n])));
    for (int k = 0; k <= n; k++) expr_free(a[k]);
    free(a); expr_free(g);

    /* ---- depression: v = -c_{n-1}/n, W_i = w^(i)/w by the Bell recurrence ---- */
    Expr* v = ds_simplify(ga_div(ga_mul(expr_new_integer(-1), expr_copy(c[n - 1])),
                                 expr_new_integer(n)));
    Expr** W = malloc((size_t)(n + 1) * sizeof(Expr*));
    W[0] = expr_new_integer(1);
    for (int i = 0; i < n; i++)
        W[i + 1] = ds_simplify(ga_add(ds_d(expr_copy(W[i]), expr_new_symbol(xv)),
                                      ga_mul(expr_copy(v), expr_copy(W[i]))));

    /* depressed coefficients d_j = Sum_{k>=j} c_k Binomial[k,j] W_{k-j} */
    Expr** d = malloc((size_t)(n + 1) * sizeof(Expr*));
    for (int j = 0; j <= n; j++) {
        Expr* s = expr_new_integer(0);
        for (int k = j; k <= n; k++)
            s = ga_add(s, ga_mul(ga_mul(expr_new_integer((long)ga_binom(k, j)), expr_copy(c[k])),
                                 expr_copy(W[k - j])));
        d[j] = ds_simplify(s);
    }
    for (int i = 0; i <= n; i++) expr_free(W[i]);
    free(W);

    /* every intermediate coefficient must vanish: u^(n) == -d_0 u */
    bool pure = true;
    for (int j = 1; j <= n - 1 && pure; j++) pure = ds_is_zero(d[j]);

    Expr* body = NULL;
    if (pure && !ds_is_zero(d[0])) {
        /* T = -d_0 = A x^m */
        Expr* T = ds_simplify(ga_mul(expr_new_integer(-1), expr_copy(d[0])));
        Expr* x = expr_new_symbol(xv);
        Expr* m = ds_simplify(ga_mul(expr_copy(x),                        /* m = x T'/T */
                      ga_mul(ds_d(expr_copy(T), expr_copy(x)), ga_powi(expr_copy(T), -1))));
        Expr* A = ds_simplify(ga_mul(expr_copy(T),                        /* A = T / x^m */
                      ga_powi(ga_pow(expr_copy(x), expr_copy(m)), -1)));
        if (ds_free_of(m, xv) && ds_free_of(A, xv) && !ds_is_zero(A)) {
            Expr* p = ds_simplify(ga_add(expr_copy(m), expr_new_integer(n)));
            if (!ds_is_zero(p)) {
                /* the gauge factor w = Exp[Integrate[v]] (1 when v == 0) */
                Expr* w = NULL;
                if (ds_is_zero(v)) {
                    w = expr_new_integer(1);
                } else {
                    Expr* iv = ds_integrate(expr_copy(v), expr_new_symbol(xv));
                    bool okv = !ds_has_head(iv, SYM_Integrate);
                    if (okv) {
                        Expr* chk = ga_sub(ds_d(expr_copy(iv), expr_new_symbol(xv)), expr_copy(v));
                        okv = ds_is_zero(chk);
                        expr_free(chk);
                    }
                    if (okv) w = ds_simplify(ga_fn1("Exp", iv)); else expr_free(iv);
                }
                if (w) {
                    bool degenerate = false;
                    Expr** basis = malloc((size_t)n * sizeof(Expr*));
                    int nb = 0;
                    for (int j = 0; j < n && !degenerate; j++) {
                        Expr* uj = ga_basis_member(j, n, p, A, xv, &degenerate);
                        if (uj) basis[nb++] = ds_simplify(ga_mul(expr_copy(w), uj));
                    }
                    if (!degenerate && nb == n) {
                        Expr* gen = expr_new_integer(0);
                        for (int j = 0; j < n; j++)
                            gen = ga_add(gen, ga_mul(ds_const(j + 1), expr_copy(basis[j])));
                        if (gg) {                 /* forcing: variation of parameters */
                            Expr* one = expr_new_integer(1);   /* monic leading coefficient */
                            Expr* part = dsolve_variation_of_parameters(basis, (size_t)n, gg,
                                             one, xv);
                            expr_free(one);
                            if (part && !ds_has_head(part, SYM_Integrate)) {
                                gen = ga_add(gen, part);
                            } else {
                                if (part) expr_free(part);
                                expr_free(gen); gen = NULL;
                            }
                        }
                        if (gen && ga_num_ok(P, gen, xv, P->fun_names[0], n)) body = gen;
                        else if (gen) expr_free(gen);
                    }
                    for (int j = 0; j < nb; j++) expr_free(basis[j]);
                    free(basis);
                    expr_free(w);
                }
            }
            expr_free(p);
        }
        expr_free(m); expr_free(A); expr_free(T); expr_free(x);
    }

    for (int j = 0; j <= n; j++) expr_free(d[j]);
    free(d);
    for (int k = 0; k <= n; k++) expr_free(c[k]);
    free(c);
    expr_free(v);
    if (gg) expr_free(gg);

    if (!body) return NULL;
    Expr** out = malloc(sizeof(Expr*));
    out[0] = body;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_genairy(Expr* res) {
    return dsolve_method_builtin(res, dsolve_genairy_try);
}

void dsolve_genairy_init(void) {
    symtab_add_builtin("DSolve`GeneralizedAiry", builtin_dsolve_genairy);
    symtab_get_def("DSolve`GeneralizedAiry")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`GeneralizedAiry",
        "DSolve`GeneralizedAiry[eqn, y, x] solves a linear ODE of order n >= 3 that "
        "reduces, after the depression y = w u (w = Exp[-Integrate[c[n-1]/n]]), to the "
        "pure-power potential u^(n) == A x^m u -- the higher-order analogue of Airy's "
        "equation. The fundamental set is x^j HypergeometricPFQ[{}, lower_j, A x^p/p^n] "
        "with p = m + n, for j = 0..n-1. Symbolic A and m are supported; forcing is "
        "added by variation of parameters. Declines p == 0, a degenerate (logarithmic) "
        "exponent pattern, and any non-power potential.");
}
