/*
 * dsolve_frobenius.c — DSolve`PowerSeries / DSolve`FrobeniusSeries.
 *
 * The always-available series fallback for a homogeneous second-order linear
 * ODE  y'' + P(x) y' + Q(x) y == 0, expanded about x0 = 0.  It is the last method
 * in the scalar cascade: it fires only when every closed-form method has
 * declined, and returns a truncated SeriesData solution.
 *
 * The point x0 = 0 is classified by the pole orders of P and Q there:
 *   - ordinary        (P, Q analytic)          -> two power series, y = C1 s1 + C2 s2
 *   - regular singular (xP, x^2 Q analytic)     -> Frobenius y = x^s Sum a_n x^n:
 *        * indicial quadratic  s(s-1) + P0 s + Q0 == 0  (P0 = lim xP, Q0 = lim x^2 Q)
 *        * root difference not a non-negative integer -> two independent series
 *        * equal roots -> second solution carries a Log (via d/ds of the series)
 *        * positive-integer difference -> second series when unobstructed, else
 *          (a genuine Log is required) the method declines
 *   - irregular singular                          -> decline
 *
 * The recurrence coefficients are exact (each a_n is solved to annihilate the
 * corresponding residual coefficient), so the truncated residual is O[x]^k and
 * back-substitution verification (dsolve_run) keeps the solution.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

#define FROB_ORDER 6   /* number of series terms past the leading one (a_0..a_N) */

/* ---- small evaluated builders (args consumed, result owned) ---- */
static Expr* T2(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Times, a, b)); }
static Expr* A2(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Plus,  a, b)); }
static Expr* Neg(Expr* a)         { return eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), a)); }
static Expr* Inv(Expr* a) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ a, expr_new_integer(-1) }, 2));
}
static Expr* PowE(Expr* base, Expr* e) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ base, e }, 2));
}

/* Taylor coefficient [x^k] e about x=0 = (1/k!) (D[e,{x,k}] /. x->0).  e borrowed. */
static Expr* taylor_coeff(const Expr* e, const char* x, int k) {
    Expr* d = expr_copy((Expr*)e);
    for (int i = 0; i < k; i++) d = ds_d(d, expr_new_symbol(x));
    Expr* at0 = ds_subst(d, expr_new_symbol(x), expr_new_integer(0));
    Expr* kfac = eval_and_free(ds_call1("Factorial", expr_new_integer(k)));
    return ds_simplify(T2(at0, Inv(kfac)));
}

/* Is `v` a finite value (no infinity / indeterminate / undecided Limit)? */
static bool is_finite_value(const Expr* v) {
    static const char* bad[] = { "ComplexInfinity", "DirectedInfinity",
                                 "Infinity", "Indeterminate", "Limit" };
    for (size_t i = 0; i < 5; i++)
        if (ds_contains(v, intern_symbol(bad[i]))) return false;
    return true;
}

/* Pole order of `e` at x=0 (0 = analytic); -1 if worse than order `maxm`. */
static int pole_order(const Expr* e, const char* x, int maxm) {
    for (int m = 0; m <= maxm; m++) {
        Expr* t = (m == 0) ? expr_copy((Expr*)e)
                           : T2(PowE(expr_new_symbol(x), expr_new_integer(m)), expr_copy((Expr*)e));
        Expr* lim = eval_and_free(ds_call2("Limit", t,
                        expr_new_function(expr_new_symbol(SYM_Rule),
                            (Expr*[]){ expr_new_symbol(x), expr_new_integer(0) }, 2)));
        bool fin = is_finite_value(lim);
        expr_free(lim);
        if (fin) return m;
    }
    return -1;
}

/* SeriesData[x, 0, {coeffs}, nmin, nmax, den]; coeffs array elements adopted. */
static Expr* mk_seriesdata(const char* x, Expr** coeffs, size_t ncoef,
                           int nmin, int nmax, int den) {
    Expr* clist = expr_new_function(expr_new_symbol(SYM_List), coeffs, ncoef);
    Expr* args[6] = { expr_new_symbol(x), expr_new_integer(0), clist,
                      expr_new_integer(nmin), expr_new_integer(nmax), expr_new_integer(den) };
    return expr_new_function(expr_new_symbol("SeriesData"), args, 6);
}

/* x^r * Sum_{i=0}^N a[i] x^i, as Power[x,r] * SeriesData (a borrowed). */
static Expr* xr_series(const char* x, const Expr* r, Expr** a, int N) {
    Expr** cc = malloc((size_t)(N + 1) * sizeof(Expr*));
    for (int i = 0; i <= N; i++) cc[i] = expr_copy(a[i]);
    Expr* sd = mk_seriesdata(x, cc, (size_t)(N + 1), 0, N + 1, 1);
    free(cc);
    return T2(PowE(expr_new_symbol(x), expr_copy((Expr*)r)), sd);
}

/* Indicial polynomial value m(m-1) + P0 m + Q0; m consumed, P0/Q0 borrowed. */
static Expr* indicial_poly(const Expr* P0, const Expr* Q0, Expr* m) {
    Expr* mm1 = T2(expr_copy(m), A2(expr_copy(m), expr_new_integer(-1)));
    Expr* pm  = T2(expr_copy((Expr*)P0), expr_copy(m));
    expr_free(m);
    return A2(A2(mm1, pm), expr_copy((Expr*)Q0));
}

/* ------------------------------------------------------------------ *
 *  Ordinary point: two power series folded into one SeriesData        *
 * ------------------------------------------------------------------ */
static Expr* frobenius_ordinary(const Expr* Pc, const Expr* Qc, const char* x, int N) {
    Expr** p = malloc((size_t)(N + 1) * sizeof(Expr*));
    Expr** q = malloc((size_t)(N + 1) * sizeof(Expr*));
    for (int k = 0; k <= N; k++) { p[k] = taylor_coeff(Pc, x, k); q[k] = taylor_coeff(Qc, x, k); }

    Expr** a = malloc((size_t)(N + 1) * sizeof(Expr*));
    a[0] = ds_const(1);   /* a_0 = C[1] */
    a[1] = ds_const(2);   /* a_1 = C[2] */
    for (int n = 0; n <= N - 2; n++) {
        /* (n+2)(n+1) a_{n+2} + Sum_{k=0}^n [p_k (n-k+1) a_{n-k+1} + q_k a_{n-k}] = 0 */
        Expr* S = expr_new_integer(0);
        for (int k = 0; k <= n; k++) {
            Expr* t1 = T2(T2(expr_copy(p[k]), expr_new_integer(n - k + 1)), expr_copy(a[n - k + 1]));
            Expr* t2 = T2(expr_copy(q[k]), expr_copy(a[n - k]));
            S = A2(S, A2(t1, t2));
        }
        a[n + 2] = ds_simplify(T2(Neg(S), Inv(expr_new_integer((n + 2) * (n + 1)))));
    }

    Expr** cc = malloc((size_t)(N + 1) * sizeof(Expr*));
    for (int i = 0; i <= N; i++) cc[i] = expr_copy(a[i]);
    Expr* body = mk_seriesdata(x, cc, (size_t)(N + 1), 0, N + 1, 1);
    free(cc);

    for (int k = 0; k <= N; k++) { expr_free(p[k]); expr_free(q[k]); expr_free(a[k]); }
    free(p); free(q); free(a);
    return body;
}

/* ------------------------------------------------------------------ *
 *  Frobenius recurrence for one root r                                *
 * ------------------------------------------------------------------ */
/* Fill a[0..N] (a[0]=1) for the exponent `r`.  Returns 0 on success, 1 when a
 * genuine Log is required (the recurrence is obstructed with nonzero numerator),
 * in which case a[] is filled with placeholders for the caller to free. */
static int build_coeffs(Expr** Pk, Expr** Qk, const Expr* P0, const Expr* Q0,
                        int N, const Expr* r, Expr** a) {
    a[0] = expr_new_integer(1);
    for (int n = 1; n <= N; n++) {
        Expr* num = expr_new_integer(0);
        for (int k = 1; k <= n; k++) {
            Expr* rnk  = A2(expr_copy((Expr*)r), expr_new_integer(n - k));   /* r+n-k */
            Expr* coef = A2(T2(expr_copy(Pk[k]), rnk), expr_copy(Qk[k]));    /* P_k(r+n-k)+Q_k */
            num = A2(num, T2(coef, expr_copy(a[n - k])));
        }
        num = ds_simplify(num);
        Expr* Frn = ds_simplify(indicial_poly(P0, Q0, A2(expr_copy((Expr*)r), expr_new_integer(n))));
        if (ds_is_zero(Frn)) {
            expr_free(Frn);
            if (ds_is_zero(num)) { expr_free(num); a[n] = expr_new_integer(0); continue; }
            expr_free(num);
            for (int j = n; j <= N; j++) a[j] = expr_new_integer(0);
            return 1;   /* Log required */
        }
        a[n] = ds_simplify(T2(Neg(num), Inv(Frn)));
    }
    return 0;
}

/* ------------------------------------------------------------------ *
 *  Regular singular point                                             *
 * ------------------------------------------------------------------ */
static Expr* frobenius_regsing(const Expr* Pc, const Expr* Qc, const char* x, int N) {
    /* xP and x^2 Q are analytic at 0 */
    Expr* xP  = ds_simplify(T2(PowE(expr_new_symbol(x), expr_new_integer(1)), expr_copy((Expr*)Pc)));
    Expr* x2Q = ds_simplify(T2(PowE(expr_new_symbol(x), expr_new_integer(2)), expr_copy((Expr*)Qc)));
    Expr** Pk = malloc((size_t)(N + 1) * sizeof(Expr*));
    Expr** Qk = malloc((size_t)(N + 1) * sizeof(Expr*));
    for (int k = 0; k <= N; k++) { Pk[k] = taylor_coeff(xP, x, k); Qk[k] = taylor_coeff(x2Q, x, k); }
    expr_free(xP); expr_free(x2Q);
    Expr* P0 = Pk[0]; Expr* Q0 = Qk[0];

    /* indicial roots */
    const char* M = intern_symbol("DSolve`fs");
    Expr* Fpoly = indicial_poly(P0, Q0, expr_new_symbol(M));
    DSolveRoots R;
    bool haveroots = dsolve_analyze_roots(Fpoly, M, 2, &R);
    expr_free(Fpoly);

    Expr* body = NULL;
    if (haveroots && R.total == 2) {
        Expr* r1 = NULL; Expr* r2 = NULL; bool equal = false;
        if (R.ndist == 1) { r1 = expr_copy(R.roots[0]); r2 = expr_copy(R.roots[0]); equal = true; }
        else {
            /* order by real part so r1 is the larger root */
            Expr* d = eval_and_free(ds_call1("Re",
                          eval_and_free(ds_call2(SYM_Subtract, expr_copy(R.roots[0]), expr_copy(R.roots[1])))));
            Expr* sgn = eval_and_free(ds_call1("Sign", d));
            bool first_bigger = !(sgn->type == EXPR_INTEGER && sgn->data.integer < 0);
            expr_free(sgn);
            r1 = expr_copy(R.roots[first_bigger ? 0 : 1]);
            r2 = expr_copy(R.roots[first_bigger ? 1 : 0]);
        }

        Expr** a1 = malloc((size_t)(N + 1) * sizeof(Expr*));
        if (equal) {
            /* symbolic recurrence in RR so we can differentiate w.r.t. the exponent */
            const char* RR = intern_symbol("DSolve`fr");
            Expr* rsym = expr_new_symbol(RR);
            Expr** aR = malloc((size_t)(N + 1) * sizeof(Expr*));
            (void)build_coeffs(Pk, Qk, P0, Q0, N, rsym, aR);   /* never obstructed for RR symbolic */
            Expr** a2 = malloc((size_t)(N + 1) * sizeof(Expr*));
            for (int i = 0; i <= N; i++) {
                a1[i] = ds_simplify(ds_subst(expr_copy(aR[i]), expr_new_symbol(RR), expr_copy(r1)));
                Expr* dai = ds_d(expr_copy(aR[i]), expr_new_symbol(RR));
                a2[i] = ds_simplify(ds_subst(dai, expr_new_symbol(RR), expr_copy(r1)));
            }
            Expr* y1 = xr_series(x, r1, a1, N);
            Expr* y2corr = xr_series(x, r1, a2, N);
            /* y2 = y1 Log[x] + x^s Sum a_n'(s) x^n */
            Expr* y2 = A2(T2(expr_copy(y1), ds_call1("Log", expr_new_symbol(x))), y2corr);
            body = A2(T2(ds_const(1), y1), T2(ds_const(2), y2));
            for (int i = 0; i <= N; i++) { expr_free(aR[i]); expr_free(a2[i]); }
            free(aR); free(a2);
            expr_free(rsym);
        } else {
            int obstr1 = build_coeffs(Pk, Qk, P0, Q0, N, r1, a1);
            Expr** a2 = malloc((size_t)(N + 1) * sizeof(Expr*));
            int obstr2 = build_coeffs(Pk, Qk, P0, Q0, N, r2, a2);
            if (!obstr1 && !obstr2) {
                Expr* s1 = xr_series(x, r1, a1, N);
                Expr* s2 = xr_series(x, r2, a2, N);
                body = A2(T2(ds_const(1), s1), T2(ds_const(2), s2));
            }
            for (int i = 0; i <= N; i++) expr_free(a2[i]);
            free(a2);
        }
        for (int i = 0; i <= N; i++) expr_free(a1[i]);
        free(a1);
        expr_free(r1); expr_free(r2);
    }
    if (haveroots) dsolve_roots_free(&R);

    for (int k = 0; k <= N; k++) { expr_free(Pk[k]); expr_free(Qk[k]); }
    free(Pk); free(Qk);
    return body;
}

/* ------------------------------------------------------------------ *
 *  Cascade entry                                                      *
 * ------------------------------------------------------------------ */
Expr** dsolve_frobenius_try(DSolveProblem* P, size_t* nbranch) {
    Expr* Pc; Expr* Qc;
    if (!dsolve_second_order_PQ(P, &Pc, &Qc)) return NULL;
    const char* x = P->ind_names[0];
    int N = FROB_ORDER;

    int pord = pole_order(Pc, x, 2);
    int qord = pole_order(Qc, x, 2);

    Expr* body = NULL;
    if (pord == 0 && qord == 0)
        body = frobenius_ordinary(Pc, Qc, x, N);
    else if (pord >= 0 && pord <= 1 && qord >= 0 && qord <= 2)
        body = frobenius_regsing(Pc, Qc, x, N);   /* else irregular -> decline */

    expr_free(Pc); expr_free(Qc);
    if (!body) return NULL;
    Expr** out = malloc(sizeof(Expr*));
    out[0] = body;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_frobenius(Expr* res) {
    return dsolve_method_builtin(res, dsolve_frobenius_try);
}

void dsolve_frobenius_init(void) {
    symtab_add_builtin("DSolve`PowerSeries", builtin_dsolve_frobenius);
    symtab_get_def("DSolve`PowerSeries")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`PowerSeries",
        "DSolve`PowerSeries[eqn, y, x] gives a truncated power-series solution of a "
        "second-order linear ODE about an ordinary point x == 0, as a SeriesData.");

    symtab_add_builtin("DSolve`FrobeniusSeries", builtin_dsolve_frobenius);
    symtab_get_def("DSolve`FrobeniusSeries")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`FrobeniusSeries",
        "DSolve`FrobeniusSeries[eqn, y, x] gives a truncated Frobenius series solution "
        "y == x^s Sum a_n x^n of a second-order linear ODE about a regular singular "
        "point x == 0 (indicial quadratic s(s-1)+P0 s+Q0 == 0; equal or "
        "non-integer-difference roots handled, with a Log term for equal roots).");
}
