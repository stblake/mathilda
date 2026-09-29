/*
 * dsolve_operator_factor.c — DSolve`OperatorFactor and DSolve`DFactor.
 *
 * Factors a linear differential operator  L = Sum_{k=0}^{n} a_k(x) D^k  by finding
 * a FIRST-ORDER RIGHT FACTOR  (D - r),  r in C(x).  (D - r) right-divides L iff
 * y = Exp[Integrate[r]] solves L[y] == 0 (a hyperexponential solution), iff the
 * "Riccati" residual  R(r) = Sum a_k P_k(r)  vanishes, where the P_k are the Bell
 * polynomials  P_0 = 1,  P_{k+1} = P_k' + r P_k  (so P_k(r) = D^k(Exp[Int r])/Exp[Int r]).
 *
 * DSolve`OperatorFactor (a scalar cascade method, order >= 2) peels TWO ways:
 *
 *   RIGHT factor.  Find one rational r, peel via operator right-division to the
 *   order-(n-1) quotient Q = L/(D-r), recurse DSolve on Q[z] == g, then close with the
 *   first-order linear solve (D - r)y == z.  This reuses the whole cascade (the
 *   quotient may be solved by const-coeff / Euler / Kovacic / a further peel), exactly
 *   like ExactODE / ReductionOfOrder.  Order 2 is admitted too: Kovacic runs FIRST and
 *   owns the tidy answers there, so what reaches here is the rational-Riccati Case 1 it
 *   declined -- and a closed form from this search beats the truncated Frobenius series
 *   that would otherwise win (it is also what makes the order-3 peel's recursion on its
 *   order-2 quotient terminate in closed form).
 *
 *   LEFT factor (Beke, order-(n-1) right factors).  When no first-order right factor
 *   exists, run the SAME search on the adjoint L* = Sum (-1)^k D^k o a_k: since
 *   (A o B)* = B* o A* and (D - s)* = -(D + s), a first-order right factor of L* is a
 *   first-order LEFT factor (D + s) of L, i.e. an order-(n-1) RIGHT factor Q.  Solve
 *   Q[z] == 0 by recursion, then L[y] == g reduces to Q[y] == W with
 *   W = Exp[-Int s](Int g Exp[Int s] dx + C[n]), closed by variation of parameters over
 *   Q's fundamental set.  For n = 3 this is exactly the classical "2nd-order right
 *   factor" case, reached without exterior powers.
 *
 * Forcing is carried through both peels (the quotient equation keeps the monic-
 * normalised right-hand side), so an inhomogeneous reducible operator solves.
 *
 * DSolve`DFactor[eqn, y, x]: the standalone factoriser.  Returns the ordered factor
 * list  {Dx - r1, Dx - r2, ...}  (Dx an inert d/dx operator symbol), INNERMOST FIRST:
 * L = (Dx - rk) o ... o (Dx - r2) o (Dx - r1), so r1 is the right factor found first
 * and any left factor peeled through the adjoint is emitted LAST.  A not-fully-
 * reducible operator returns its peeled factors plus the inert remainder operator
 * Sum q_j Dx^j.
 *
 * Right-division recurrence (monic a_n = 1):  q_{n-1} = a_n;
 *   q_{m-1} = a_m + Sum_{j=m}^{n-1} C(j,m) D^(j-m)[r] q_j   (m = n-1..1);
 * remainder rho = a_0 + Sum_{j} q_j r^(j)  ==  R(r)  (the factor divides iff rho == 0).
 * Left-division recurrence:  q_{n-1} = a_n;  q_{m-1} = a_m - q_m' - s q_m;
 * remainder rho = a_0 - q_0' - s q_0.
 *
 * Self-contained: reuses only the shared ds_* substrate and front-end builtins; makes
 * NO changes to dsolve_kovacic.c (the small ansatz/solve overlap is duplicated to keep
 * that engine at zero regression risk).
 *
 * Scope: linear, order >= 2 (OperatorFactor) / >= 1 (DFactor), homogeneous or forced;
 * first-order right and left factors with r, s in C(x) (simple poles + low-degree
 * polynomial part), rational coefficients with numeric constants.  Deferred:
 * irregular-singular (double-pole) r, irreducible-quadratic-denominator r, symbolic
 * parameters in the coefficients, and right factors of order m with 2 <= m <= n-2
 * (which genuinely need the m-th exterior power).
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include "../internal.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* ---- small evaluated builders (args consumed, result owned) ---- */
static Expr* T2(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Times, a, b)); }
static Expr* A2(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Plus,  a, b)); }
static Expr* Sub(Expr* a, Expr* b){ return eval_and_free(ds_call2(SYM_Subtract, a, b)); }
static Expr* Powi(Expr* b, int e) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ b, expr_new_integer(e) }, 2));
}
static Expr* fn1(const char* h, Expr* a)          { return eval_and_free(ds_call1(h, a)); }
static Expr* fn2(const char* h, Expr* a, Expr* b) { return eval_and_free(ds_call2(h, a, b)); }

static int of_degree_in(const Expr* poly, const char* x) {
    Expr* e = fn2("Exponent", expr_copy((Expr*)poly), expr_new_symbol(x));
    int d = (e->type == EXPR_INTEGER) ? (int)e->data.integer : -1;
    expr_free(e);
    return d;
}

static long of_binom(int nn, int kk) {
    if (kk < 0 || kk > nn) return 0;
    long r = 1;
    for (int i = 0; i < kk; i++) r = r * (nn - i) / (i + 1);
    return r;
}

/* Extract the RHS body of {{z[x] -> expr}} (applied form) from a DSolve result. */
static Expr* extract_applied(Expr* r, const char* zfun) {
    if (!r || !head_is(r, SYM_List) || r->data.function.arg_count < 1) return NULL;
    Expr* inner = r->data.function.args[0];
    if (!head_is(inner, SYM_List)) return NULL;
    for (size_t k = 0; k < inner->data.function.arg_count; k++) {
        Expr* rule = inner->data.function.args[k];
        if (head_is(rule, SYM_Rule) && rule->data.function.arg_count == 2) {
            Expr* lhs = rule->data.function.args[0];
            if (lhs->type == EXPR_FUNCTION && lhs->data.function.head->type == EXPR_SYMBOL
                && lhs->data.function.head->data.symbol.name == zfun)
                return expr_copy(rule->data.function.args[1]);
        }
    }
    return NULL;
}

/* Largest index k of any C[k] occurring in e (0 if none). */
static void of_scan_const(const Expr* e, int* mx) {
    if (!e || e->type != EXPR_FUNCTION) return;
    Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && strcmp(h->data.symbol.name, "C") == 0
        && e->data.function.arg_count == 1
        && e->data.function.args[0]->type == EXPR_INTEGER) {
        int idx = (int)e->data.function.args[0]->data.integer;
        if (idx > *mx) *mx = idx;
    }
    of_scan_const(h, mx);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        of_scan_const(e->data.function.args[i], mx);
}

/* Require every coefficient to be a rational function of x with numeric constants
 * (no symbolic parameters) — probe x -> 17/13 and demand NumberQ. */
static bool of_coeffs_numeric(Expr** c, int n, const char* x) {
    const char* T = intern_symbol("True");
    Expr* probe = T2(expr_new_integer(17), Powi(expr_new_integer(13), -1));  /* 17/13 */
    bool ok = true;
    for (int k = 0; k <= n && ok; k++) {
        Expr* v = ds_subst(expr_copy(c[k]), expr_new_symbol(x), expr_copy(probe));
        Expr* q = fn1("NumberQ", v);
        ok = (q->type == EXPR_SYMBOL && q->data.symbol.name == T);
        expr_free(q);
    }
    expr_free(probe);
    return ok;
}

/* Build r-ansatz: polynomial part (degree poly_deg) + principal parts at each
 * denominator factor (numerator degree < deg f, powers 1..pole_order).  Unknown
 * symbols "DSolve`of<counter>" are pushed into *unk. */
static Expr* of_build_ansatz(const char* x, int poly_deg, Expr* factors,
                             int pole_order, int* counter, Expr*** unk, size_t* nu) {
    /* count exactly (poly part + per-factor principal parts) to size the arrays */
    size_t total = (size_t)(poly_deg + 1);
    if (factors && head_is(factors, SYM_List)) {
        for (size_t fi = 0; fi < factors->data.function.arg_count; fi++) {
            Expr* pair = factors->data.function.args[fi];
            if (!head_is(pair, SYM_List) || pair->data.function.arg_count != 2) continue;
            int df = of_degree_in(pair->data.function.args[0], x);
            if (df >= 1) total += (size_t)(pole_order * df);
        }
    }
    Expr** U = malloc(total * sizeof(Expr*));
    Expr** terms = malloc(total * sizeof(Expr*));
    size_t nU = 0, nt = 0;

    for (int j = 0; j <= poly_deg; j++) {
        char buf[32]; snprintf(buf, sizeof(buf), "DSolve`of%d", (*counter)++);
        const char* sn = intern_symbol(buf);
        U[nU++] = expr_new_symbol(sn);
        terms[nt++] = ds_call2(SYM_Times, expr_new_symbol(sn),
                        expr_new_function(expr_new_symbol(SYM_Power),
                            (Expr*[]){ expr_new_symbol(x), expr_new_integer(j) }, 2));
    }
    if (factors && head_is(factors, SYM_List)) {
        for (size_t fi = 0; fi < factors->data.function.arg_count; fi++) {
            Expr* pair = factors->data.function.args[fi];
            if (!head_is(pair, SYM_List) || pair->data.function.arg_count != 2) continue;
            Expr* f = pair->data.function.args[0];
            int df = of_degree_in(f, x);
            if (df < 1) continue;                     /* constant content */
            for (int k = 1; k <= pole_order; k++)
                for (int l = 0; l < df; l++) {
                    char buf[32]; snprintf(buf, sizeof(buf), "DSolve`of%d", (*counter)++);
                    const char* sn = intern_symbol(buf);
                    U[nU++] = expr_new_symbol(sn);
                    terms[nt++] = ds_call2(SYM_Times, expr_new_symbol(sn),
                                    ds_call2(SYM_Times,
                                        expr_new_function(expr_new_symbol(SYM_Power),
                                            (Expr*[]){ expr_new_symbol(x), expr_new_integer(l) }, 2),
                                        expr_new_function(expr_new_symbol(SYM_Power),
                                            (Expr*[]){ expr_copy(f), expr_new_integer(-k) }, 2)));
                }
        }
    }
    Expr* w = expr_new_function(expr_new_symbol(SYM_Plus), terms, nt);
    free(terms);
    *unk = U; *nu = nU;
    return w;
}

/* Riccati residual  R(r) = Sum_{k=0}^n a[k] P_k(r),  P_0=1, P_{k+1}=P_k'+r P_k. */
static Expr* of_riccati_residual(Expr** a, int n, const Expr* r, const char* x) {
    Expr* P = expr_new_integer(1);          /* P_0 */
    Expr* R = expr_copy(a[0]);              /* a[0] P_0 */
    for (int k = 1; k <= n; k++) {
        Expr* dP = ds_d(expr_copy(P), expr_new_symbol(x));
        Expr* rP = T2(expr_copy((Expr*)r), P);     /* consumes P */
        P = A2(dP, rP);                            /* P_k */
        R = A2(R, T2(expr_copy(a[k]), expr_copy(P)));
    }
    expr_free(P);
    return R;
}

/* Exact right-division of monic L (a[0..n], a[n]=1) by (D - r).  Returns q[0..n-1]
 * (malloc'd) and, via *rho_out, the remainder (== the Riccati residual). */
static Expr** of_divide(Expr** a, int n, const Expr* r, const char* x, Expr** rho_out) {
    Expr** rd = malloc((size_t)n * sizeof(Expr*));
    rd[0] = expr_copy((Expr*)r);
    for (int i = 1; i < n; i++) rd[i] = ds_d(expr_copy(rd[i-1]), expr_new_symbol(x));

    Expr** q = malloc((size_t)n * sizeof(Expr*));
    q[n-1] = expr_copy(a[n]);
    for (int m = n-1; m >= 1; m--) {
        Expr* s = expr_copy(a[m]);
        for (int j = m; j <= n-1; j++) {
            long b = of_binom(j, m);
            Expr* term = T2(T2(expr_new_integer(b), expr_copy(rd[j-m])), expr_copy(q[j]));
            s = A2(s, term);
        }
        q[m-1] = ds_simplify(s);
    }
    Expr* rho = expr_copy(a[0]);
    for (int j = 0; j < n; j++) rho = A2(rho, T2(expr_copy(q[j]), expr_copy(rd[j])));
    *rho_out = ds_simplify(rho);

    for (int i = 0; i < n; i++) expr_free(rd[i]);
    free(rd);
    return q;
}

/* Search for a first-order right factor (D - r), r in C(x), of the monic operator
 * a[0..n].  Returns r (owned) or NULL.  Iterates cheapest-first over polynomial
 * degree and pole order; for each Solve branch tests the exact division remainder. */
#define OF_MAX_POLES 3   /* bound the cubic determining system before any ds_solve */
/* ...and bound its WIDTH as well.  The pole COUNT alone does not bound the ansatz:
 * two poles at multiplicity 2 with a degree-2 polynomial part already give 11
 * unknowns, and ds_solve on a cubic system that wide does not return (an order-3
 * operator with x^2 (1+x) denominators hangs there).  Every factor this search has
 * ever found needed a narrow ansatz -- 2.1.2-250 uses 4 unknowns, -253 uses 3 -- so
 * the wide combinations are cost with no recorded benefit. */
#define OF_MAX_UNKNOWNS 6
/* Pole order of the ansatz.  A double pole in r is an IRREGULAR singularity, which
 * the file header lists as deferred scope anyway; searching for one doubles the
 * width of every combination and is what makes a FAILING search expensive. */
#define OF_MAX_POLE_ORDER 1

static Expr* of_find_factor(Expr** a, int n, const char* x) {
    Expr* sum = expr_copy(a[0]);
    for (int k = 1; k <= n; k++) sum = A2(sum, expr_copy(a[k]));
    Expr* factors = fn1("FactorList", fn1("Denominator", fn1("Together", sum)));

    /* Complexity gate: the determining system is CUBIC (Bell polynomials
     * P_0..P_3), so a many-pole ansatz sends ds_solve into unbounded parametric
     * elimination — e.g. the rational-coefficient E18-class ODE with four
     * distinct poles.  Bound the distinct pole count up front from the
     * (polynomial-time) FactorList and decline to the series fallback rather than
     * hang.  (The rational/polynomial-solution route that would close those is a
     * later, linear-system method.) */
    {
        int npoles = 0;
        if (factors && head_is(factors, SYM_List))
            for (size_t fi = 0; fi < factors->data.function.arg_count; fi++) {
                Expr* pair = factors->data.function.args[fi];
                if (head_is(pair, SYM_List) && pair->data.function.arg_count == 2
                    && !ds_free_of(pair->data.function.args[0], x)) npoles++;
            }
        if (npoles > OF_MAX_POLES) { expr_free(factors); return NULL; }
    }

    /* total principal-part width: Sum over non-constant factors of deg(f), used with
     * the pole order and polynomial degree to size the ansatz below */
    int pole_width = 0;
    if (factors && head_is(factors, SYM_List))
        for (size_t fi = 0; fi < factors->data.function.arg_count; fi++) {
            Expr* pair = factors->data.function.args[fi];
            if (!head_is(pair, SYM_List) || pair->data.function.arg_count != 2) continue;
            int df = of_degree_in(pair->data.function.args[0], x);
            if (df >= 1) pole_width += df;
        }

    Expr* found = NULL;
    for (int pd = 0; pd <= 2 && !found; pd++) {
        for (int po = 1; po <= OF_MAX_POLE_ORDER && !found; po++) {
            if (pd + 1 + po * pole_width > OF_MAX_UNKNOWNS) continue;
            int counter = 1;
            Expr** unk; size_t nu;
            Expr* r = of_build_ansatz(x, pd, factors, po, &counter, &unk, &nu);
            Expr* R = of_riccati_residual(a, n, r, x);
            Expr* num = fn1("Numerator", fn1("Together", R));
            Expr* clist = fn2("CoefficientList", num, expr_new_symbol(x));
            if (clist && head_is(clist, SYM_List)) {
                size_t nc = clist->data.function.arg_count;
                Expr** eqs = malloc((nc ? nc : 1) * sizeof(Expr*));
                size_t ne = 0;
                for (size_t i = 0; i < nc; i++) {
                    Expr* co = clist->data.function.args[i];
                    if (ds_is_zero(co)) continue;
                    eqs[ne++] = expr_new_function(expr_new_symbol(SYM_Equal),
                                    (Expr*[]){ expr_copy(co), expr_new_integer(0) }, 2);
                }
                Expr* eqlist = expr_new_function(expr_new_symbol(SYM_List), eqs, ne);
                free(eqs);
                Expr** vs = malloc((nu ? nu : 1) * sizeof(Expr*));
                for (size_t i = 0; i < nu; i++) vs[i] = expr_copy(unk[i]);
                Expr* varlist = expr_new_function(expr_new_symbol(SYM_List), vs, nu);
                free(vs);
                Expr* sol = ds_solve(eqlist, varlist);
                if (sol && head_is(sol, SYM_List)) {
                    for (size_t bi = 0; bi < sol->data.function.arg_count && !found; bi++) {
                        Expr* branch = sol->data.function.args[bi];
                        if (!head_is(branch, SYM_List)) continue;
                        Expr* cand = eval_and_free(internal_replace_all(
                                        (Expr*[]){ expr_copy(r), expr_copy(branch) }, 2));
                        for (size_t i = 0; i < nu; i++)
                            if (ds_contains(cand, unk[i]->data.symbol.name))
                                cand = ds_subst(cand, expr_copy(unk[i]), expr_new_integer(0));
                        cand = ds_simplify(cand);
                        Expr* rho; Expr** q = of_divide(a, n, cand, x, &rho);
                        bool ok = ds_is_zero(rho);
                        expr_free(rho);
                        for (int i = 0; i < n; i++) expr_free(q[i]);
                        free(q);
                        if (ok) found = expr_copy(cand);
                        expr_free(cand);
                    }
                }
                if (sol) expr_free(sol);
            }
            if (clist) expr_free(clist);
            for (size_t i = 0; i < nu; i++) expr_free(unk[i]);
            free(unk);
            expr_free(r);
        }
    }
    expr_free(factors);
    return found;
}

/* Adjoint coefficients  b_j = Sum_{k>=j} (-1)^k Binomial[k,j] D^(k-j)[a_k].
 *
 * (D - s) right-divides L* exactly when (D + s) LEFT-divides L, because
 * (A o B)* = B* o A* and (D - s)* = -(D + s).  So one run of the existing
 * first-order right-factor search on L* hands back an order-(n-1) RIGHT factor of
 * L — the "2nd-order right factor (Beke)" case for n = 3 — with no exterior
 * powers.  The candidate is never trusted: of_left_divide's remainder is the
 * exact acceptance test. */
static Expr** of_adjoint(Expr** a, int n, const char* x) {
    Expr** b = malloc((size_t)(n + 1) * sizeof(Expr*));
    for (int j = 0; j <= n; j++) {
        Expr* s = expr_new_integer(0);
        for (int k = j; k <= n; k++) {
            Expr* d = expr_copy(a[k]);
            for (int i = 0; i < k - j; i++) d = ds_d(d, expr_new_symbol(x));
            long co = of_binom(k, j) * ((k % 2) ? -1 : 1);
            s = A2(s, T2(expr_new_integer(co), d));
        }
        b[j] = ds_simplify(s);
    }
    return b;
}

/* Exact LEFT division of L (a[0..n]) by (D + s):  matching the coefficient of
 * y^(m) in (D + s)(Sum q_j D^j) = Sum (q_{m-1} + q_m' + s q_m) y^(m) gives
 *   q_{n-1} = a_n,   q_{m-1} = a_m - q_m' - s q_m   (m = n-1 .. 1),
 * with remainder rho = a_0 - q_0' - s q_0.  Returns q[0..n-1] (malloc'd). */
static Expr** of_left_divide(Expr** a, int n, const Expr* s, const char* x, Expr** rho_out) {
    Expr** q = malloc((size_t)n * sizeof(Expr*));
    q[n-1] = expr_copy(a[n]);
    for (int m = n-1; m >= 1; m--)
        q[m-1] = ds_simplify(Sub(Sub(expr_copy(a[m]),
                                     ds_d(expr_copy(q[m]), expr_new_symbol(x))),
                                 T2(expr_copy((Expr*)s), expr_copy(q[m]))));
    *rho_out = ds_simplify(Sub(Sub(expr_copy(a[0]),
                                   ds_d(expr_copy(q[0]), expr_new_symbol(x))),
                               T2(expr_copy((Expr*)s), expr_copy(q[0]))));
    return q;
}

/* Monic coefficient array a[0..n] from the parsed problem, or NULL.  Sets *order.
 * When `forcing_out` is non-NULL it receives the monic-normalised forcing g/a_n
 * (owned), or NULL when the equation is homogeneous. */
static Expr** of_monic_coeffs(DSolveProblem* P, int* order, Expr** forcing_out) {
    Expr** c; Expr* g; int n;
    if (!dsolve_linear_coeffs(P, &c, &g, &n)) return NULL;
    bool homog = ds_is_zero(g);
    const char* x = P->ind_names[0];
    if (n < 1 || ds_is_zero(c[n]) || !of_coeffs_numeric(c, n, x)) {
        for (int k = 0; k <= n; k++) expr_free(c[k]);
        free(c);
        expr_free(g);
        return NULL;
    }
    Expr** a = malloc((size_t)(n + 1) * sizeof(Expr*));
    for (int k = 0; k <= n; k++)
        a[k] = ds_simplify(T2(expr_copy(c[k]), Powi(expr_copy(c[n]), -1)));
    if (forcing_out)
        *forcing_out = homog ? NULL : ds_simplify(T2(expr_copy(g), Powi(expr_copy(c[n]), -1)));
    expr_free(g);
    for (int k = 0; k <= n; k++) expr_free(c[k]);
    free(c);
    *order = n;
    return a;
}

/* Build and solve  Sum_{j} q_j z^(j) == rhs  by recursing DSolve; returns the
 * solution body in a fresh function, or NULL.  `rhs` borrowed (NULL = 0). */
static Expr* of_solve_quotient(Expr** q, int m, const Expr* rhs, const char* xvar,
                               const char* zf) {
    Expr** terms = malloc((size_t)(m + 1) * sizeof(Expr*));
    for (int j = 0; j <= m; j++)
        terms[j] = ds_call2(SYM_Times, expr_copy(q[j]), ds_make_funcapp(zf, j, xvar));
    Expr* lhs = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)(m + 1)));
    free(terms);
    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                   (Expr*[]){ lhs, rhs ? expr_copy((Expr*)rhs) : expr_new_integer(0) }, 2);
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ eq, ds_make_funcapp(zf, 0, xvar), expr_new_symbol(xvar) }, 3);
    Expr* res = eval_and_free(call);
    Expr* z = extract_applied(res, zf);
    expr_free(res);
    return z;
}

/* Integrate `e` and require the result to be elementary AND to differentiate back
 * to `e` (the guard dsolve_operfactor_try has always applied to Integrate[r]).
 * Returns the antiderivative (owned) or NULL.  `e` borrowed. */
static Expr* of_safe_integral(const Expr* e, const char* xvar) {
    Expr* I = ds_integrate(expr_copy((Expr*)e), expr_new_symbol(xvar));
    bool ok = !ds_has_head(I, SYM_Integrate);
    if (ok) {
        Expr* chk = Sub(ds_d(expr_copy(I), expr_new_symbol(xvar)), expr_copy((Expr*)e));
        ok = ds_is_zero(chk);
        expr_free(chk);
    }
    if (!ok) { expr_free(I); return NULL; }
    return I;
}

/* ---------- DSolve`OperatorFactor ---------- */

/* Integrate  z * winv  by SPLITTING z along its generated constants first:
 * z = Sum_k C[k] z_k + z_0, so the integral is Sum_k C[k] Int[z_k winv] + Int[z_0 winv].
 * Linearity makes the split exact, and it matters: the algebraic integrator is two
 * orders of magnitude slower on an integrand that carries symbolic C[k] coefficients
 * than on the same integrand with numeric ones (2.06 s vs 0.011 s for the trailing
 * integrand (C1 + C2 x)/(x^3 (1+x^2)^(3/2)) of 2.1.2-250).  Returns NULL when any
 * piece fails to close in elementary form.  z, winv borrowed. */
static Expr* of_linear_integral(const Expr* z, const Expr* winv, const char* xvar) {
    int mx = 0; of_scan_const(z, &mx);
    Expr* acc  = expr_new_integer(0);
    Expr* rest = expr_copy((Expr*)z);
    bool ok = true;
    for (int k = 1; k <= mx && ok; k++) {
        Expr* zk = fn2("Coefficient", expr_copy((Expr*)z), ds_const(k));
        if (ds_is_structural_zero(zk)) { expr_free(zk); continue; }
        rest = Sub(rest, T2(ds_const(k), expr_copy(zk)));
        Expr* Ik = ds_integrate(T2(zk, expr_copy((Expr*)winv)), expr_new_symbol(xvar));
        if (ds_has_head(Ik, SYM_Integrate)) { expr_free(Ik); ok = false; break; }
        acc = A2(acc, T2(ds_const(k), Ik));
    }
    if (ok) {
        rest = ds_simplify(rest);
        if (!ds_is_structural_zero(rest)) {
            Expr* I0 = ds_integrate(T2(expr_copy(rest), expr_copy((Expr*)winv)),
                                    expr_new_symbol(xvar));
            if (ds_has_head(I0, SYM_Integrate)) { expr_free(I0); ok = false; }
            else acc = A2(acc, I0);
        }
    }
    expr_free(rest);
    if (!ok) { expr_free(acc); return NULL; }
    return acc;
}

/* True when `e` is built only from elementary functions — no truncated series and
 * no special-function head whose variation-of-parameters integrals do not close. */
static bool of_basis_elementary(const Expr* e) {
    static const char* bad[] = { "SeriesData", "BesselJ", "BesselY", "BesselI", "BesselK",
        "AiryAi", "AiryBi", "HypergeometricPFQ", "Hypergeometric1F1", "Hypergeometric2F1",
        "LegendreP", "LegendreQ", "WhittakerM", "WhittakerW", "MathieuC", "MathieuS" };
    for (size_t i = 0; i < sizeof(bad) / sizeof(bad[0]); i++)
        if (ds_contains(e, intern_symbol(bad[i]))) return false;
    return true;
}

/* RIGHT-factor peel: L = Q o (D - r).  L[y] == g  <=>  Q[z] == g with z = (D-r)y,
 * then the trailing first-order solve  y = Exp[Int r](Int z Exp[-Int r] + C[k+1]).
 * `gg` borrowed (NULL = homogeneous). */
static Expr* of_peel_right(Expr** a, int n, const Expr* r, const Expr* gg, const char* xvar) {
    Expr* rho; Expr** q = of_divide(a, n, r, xvar, &rho); expr_free(rho);
    Expr* z = of_solve_quotient(q, n - 1, gg, xvar, intern_symbol("DSolve`ofz"));
    for (int j = 0; j < n; j++) expr_free(q[j]);
    free(q);
    if (!z) return NULL;

    Expr* body = NULL;
    Expr* intR = of_safe_integral(r, xvar);
    if (intR) {
        /* Simplify the hyperexponential factor before it reaches the integrator.
         * Int r is typically a sum of logs, and the raw Exp[-(Log[x] + Log[1+x^2]/2)]
         * sends Integrate down the transcendental-tower path, where the trailing
         * integrand of 2.1.2-250 costs 3.98 s and still FAILS; the same integrand as
         * 1/(x Sqrt[1+x^2]) is algebraic, closes, and costs 0.015 s.  The branch the
         * simplification picks is absorbed into the arbitrary constants, and
         * dsolve_run back-substitutes the assembled body regardless. */
        Expr* w    = ds_simplify(fn1("Exp", expr_copy(intR)));
        Expr* winv = ds_simplify(fn1("Exp", Sub(expr_new_integer(0), expr_copy(intR))));
        /* The whole construction rests on w' == r w, so check it rather than trust the
         * simplifier's branch choice; on failure fall back to the raw Exp form, which
         * satisfies it by construction. */
        {
            Expr* chk = Sub(ds_d(expr_copy(w), expr_new_symbol(xvar)),
                            T2(expr_copy((Expr*)r), expr_copy(w)));
            bool ok = ds_is_zero(chk);
            expr_free(chk);
            if (!ok) {
                expr_free(w); expr_free(winv);
                w    = fn1("Exp", expr_copy(intR));
                winv = fn1("Exp", Sub(expr_new_integer(0), expr_copy(intR)));
            }
        }
        Expr* part = of_linear_integral(z, winv, xvar);
        if (part) {
            int mx = 0; of_scan_const(z, &mx);
            body = T2(w, A2(part, ds_const(mx + 1)));
        } else {
            expr_free(w);
        }
        expr_free(winv);
        expr_free(intR);
    }
    expr_free(z);
    return body;
}

/* LEFT-factor peel: L = (D + s) o Q, Q of order n-1 (so Q is an order-(n-1) RIGHT
 * factor of L -- the Beke case).  L[y] == g  <=>  Q[y] == W with
 * W' + s W == g, i.e. W = Exp[-Int s](Int g Exp[Int s] dx + C[n]); the general
 * solution is Q's own general solution plus a variation-of-parameters particular
 * for the forcing W.  `gg` borrowed (NULL = homogeneous). */
static Expr* of_peel_left(Expr** a, int n, const Expr* s, const Expr* gg, const char* xvar) {
    Expr* intS = of_safe_integral(s, xvar);
    if (!intS) return NULL;

    Expr* rho; Expr** q = of_left_divide(a, n, s, xvar, &rho);
    bool divides = ds_is_zero(rho);
    expr_free(rho);
    if (!divides) {
        for (int j = 0; j < n; j++) expr_free(q[j]);
        free(q); expr_free(intS);
        return NULL;
    }

    Expr* z = of_solve_quotient(q, n - 1, NULL, xvar, intern_symbol("DSolve`ofw"));
    Expr* body = NULL;
    /* The quotient's fundamental set must be ELEMENTARY before variation of
     * parameters is attempted.  A series (Frobenius) basis cannot carry a VoP
     * particular at all, and a special-function basis essentially never yields an
     * elementary Wronskian integral -- but Integrate churns for many seconds finding
     * that out (11.9 s on 2.1.2-294, whose quotient basis is
     * Sqrt[x] BesselJ[I Sqrt[3], 2 Sqrt[x]]).  A syntactic scan is the cheap, honest
     * gate: a missed elementary case costs one decline, a missing gate costs the
     * whole solve budget. */
    if (z && !of_basis_elementary(z)) { expr_free(z); z = NULL; }
    if (z) {
        /* fundamental set of Q from the C[k] coefficients of its general solution */
        Expr** basis = malloc((size_t)(n - 1) * sizeof(Expr*));
        int nb = 0; bool ok = true;
        for (int i = 1; i <= n - 1 && ok; i++) {
            Expr* bi = eval_and_free(ds_call2("Coefficient", expr_copy(z), ds_const(i)));
            if (ds_is_structural_zero(bi)) { expr_free(bi); ok = false; }
            else basis[nb++] = bi;
        }
        if (ok) {
            int mx = 0; of_scan_const(z, &mx);
            /* W = Exp[-Int s] (C[mx+1] + Int gg Exp[Int s] dx).  Run variation of
             * parameters on the two pieces SEPARATELY (same linearity argument as
             * of_linear_integral: a forcing carrying a symbolic C[k] makes every
             * Wronskian integral far more expensive than the numeric one). */
            Expr* winv = ds_simplify(fn1("Exp", Sub(expr_new_integer(0), expr_copy(intS))));
            Expr* ph = dsolve_variation_of_parameters(basis, (size_t)(n - 1), winv,
                                                      a[n], xvar);
            bool good = (ph && !ds_has_head(ph, SYM_Integrate));
            Expr* total = good ? T2(ds_const(mx + 1), expr_copy(ph)) : NULL;
            if (ph) expr_free(ph);
            if (good && gg) {
                Expr* ig = ds_integrate(T2(expr_copy((Expr*)gg),
                                            ds_simplify(fn1("Exp", expr_copy(intS)))),
                                        expr_new_symbol(xvar));
                if (ds_has_head(ig, SYM_Integrate)) { expr_free(ig); good = false; }
                else {
                    Expr* f2 = T2(expr_copy(winv), ig);
                    Expr* pp = dsolve_variation_of_parameters(basis, (size_t)(n - 1), f2,
                                                              a[n], xvar);
                    expr_free(f2);
                    if (pp && !ds_has_head(pp, SYM_Integrate)) total = A2(total, pp);
                    else { if (pp) expr_free(pp); good = false; }
                }
            }
            if (good && total) body = A2(expr_copy(z), total);
            else if (total) expr_free(total);
            expr_free(winv);
        }
        for (int i = 0; i < nb; i++) expr_free(basis[i]);
        free(basis);
    }
    if (z) expr_free(z);
    for (int j = 0; j < n; j++) expr_free(q[j]);
    free(q);
    expr_free(intS);
    return body;
}

Expr** dsolve_operfactor_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1 || P->is_pde) return NULL;
    const char* xvar = P->ind_names[0];

    int n; Expr* gg = NULL;
    Expr** a = of_monic_coeffs(P, &n, &gg);
    if (!a) return NULL;
    /* Order 2 is admitted (Kovacic runs first and owns the tidy answers there);
     * this is the rational-Riccati Case-1 search it declines, and a closed form
     * from it beats the Frobenius series that would otherwise win. */
    if (n < 2) {
        for (int k = 0; k <= n; k++) expr_free(a[k]);
        free(a); if (gg) expr_free(gg);
        return NULL;
    }

    Expr* body = NULL;
    Expr* r = of_find_factor(a, n, xvar);
    if (r) { body = of_peel_right(a, n, r, gg, xvar); expr_free(r); }

    if (!body) {
        /* No first-order RIGHT factor: look for a first-order LEFT factor via the
         * adjoint, i.e. an order-(n-1) right factor (Beke). */
        Expr** b = of_adjoint(a, n, xvar);
        Expr* s = of_find_factor(b, n, xvar);
        for (int k = 0; k <= n; k++) expr_free(b[k]);
        free(b);
        if (s) { body = of_peel_left(a, n, s, gg, xvar); expr_free(s); }
    }

    for (int k = 0; k <= n; k++) expr_free(a[k]);
    free(a);
    if (gg) expr_free(gg);
    if (!body) return NULL;

    Expr** out = malloc(sizeof(Expr*));
    out[0] = body;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_operfactor(Expr* res) {
    return dsolve_method_builtin(res, dsolve_operfactor_try);
}

/* ---------- DSolve`DFactor ---------- */

/* Inert operator (Dx - r). */
static Expr* dfactor_first_order(const Expr* r) {
    return Sub(expr_new_symbol("Dx"), expr_copy((Expr*)r));
}

/* Inert operator  Sum_{j=0}^{m} a[j] Dx^j  (the irreducible remainder of order m). */
static Expr* dfactor_remainder(Expr** a, int m) {
    Expr** terms = malloc((size_t)(m + 1) * sizeof(Expr*));
    for (int j = 0; j <= m; j++)
        terms[j] = ds_call2(SYM_Times, expr_copy(a[j]),
                       expr_new_function(expr_new_symbol(SYM_Power),
                           (Expr*[]){ expr_new_symbol("Dx"), expr_new_integer(j) }, 2));
    Expr* op = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)(m + 1)));
    free(terms);
    return op;
}

static Expr* builtin_dsolve_dfactor(Expr* res) {
    DSolveProblem P;
    if (!dsolve_parse(res, &P)) return NULL;
    if (P.is_pde || P.nfun != 1) { dsolve_problem_free(&P); return NULL; }
    const char* xvar = P.ind_names[0];

    int n;
    Expr** a = of_monic_coeffs(&P, &n, NULL);
    dsolve_problem_free(&P);
    if (!a) return NULL;

    size_t cap = 8, nf = 0, nl = 0;
    Expr** facs = malloc(cap * sizeof(Expr*));
    Expr** lefts = malloc((size_t)(n + 1) * sizeof(Expr*));   /* outermost-first */
    int cur = n;
    while (cur >= 1) {
        if (cur == 1) {                                  /* Dx + a[0] = Dx - (-a[0]) */
            Expr* r = Sub(expr_new_integer(0), expr_copy(a[0]));
            if (nf == cap) { cap *= 2; facs = realloc(facs, cap * sizeof(Expr*)); }
            facs[nf++] = dfactor_first_order(r);
            expr_free(r);
            for (int k = 0; k <= cur; k++) expr_free(a[k]);
            free(a); a = NULL; cur = 0;
            break;
        }
        Expr* r = of_find_factor(a, cur, xvar);
        if (r) {
            if (nf == cap) { cap *= 2; facs = realloc(facs, cap * sizeof(Expr*)); }
            facs[nf++] = dfactor_first_order(r);
            Expr* rho; Expr** q = of_divide(a, cur, r, xvar, &rho); expr_free(rho);
            expr_free(r);
            for (int k = 0; k <= cur; k++) expr_free(a[k]);
            free(a);
            a = q;                                       /* order cur-1, monic q[cur-1]=1 */
            cur -= 1;
            continue;
        }
        /* No first-order right factor: try a first-order LEFT factor (D + s) via
         * the adjoint, which peels an order-(cur-1) RIGHT factor (Beke).  A left
         * factor is the OUTERMOST operator, so it is emitted last (the list is
         * innermost-first). */
        Expr** b = of_adjoint(a, cur, xvar);
        Expr* s = of_find_factor(b, cur, xvar);
        for (int k = 0; k <= cur; k++) expr_free(b[k]);
        free(b);
        if (!s) break;                                   /* irreducible tail of order cur */
        Expr* lrho; Expr** lq = of_left_divide(a, cur, s, xvar, &lrho);
        bool ok = ds_is_zero(lrho);
        expr_free(lrho);
        if (!ok) {
            for (int j = 0; j < cur; j++) expr_free(lq[j]);
            free(lq); expr_free(s);
            break;
        }
        lefts[nl++] = A2(expr_new_symbol("Dx"), expr_copy(s));   /* Dx + s */
        expr_free(s);
        for (int k = 0; k <= cur; k++) expr_free(a[k]);
        free(a);
        a = lq;
        cur -= 1;
    }
    if (a && cur >= 2) {                                 /* append irreducible remainder */
        if (nf == cap) { cap *= 2; facs = realloc(facs, cap * sizeof(Expr*)); }
        facs[nf++] = dfactor_remainder(a, cur);
    }
    if (a) { for (int k = 0; k <= cur; k++) expr_free(a[k]); free(a); }
    while (nl > 0) {                                     /* outermost last */
        if (nf == cap) { cap *= 2; facs = realloc(facs, cap * sizeof(Expr*)); }
        facs[nf++] = lefts[--nl];
    }
    free(lefts);

    Expr* out = expr_new_function(expr_new_symbol(SYM_List), facs, nf);
    free(facs);
    return out;
}

void dsolve_operfactor_init(void) {
    symtab_add_builtin("DSolve`OperatorFactor", builtin_dsolve_operfactor);
    symtab_get_def("DSolve`OperatorFactor")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`OperatorFactor",
        "DSolve`OperatorFactor[eqn, y, x] solves a linear ODE of order >= 2 (forced or "
        "homogeneous) by factoring its operator. It first seeks a first-order right "
        "factor (D - r) with r rational (a hyperexponential solution Exp[Integrate[r]]), "
        "peels it off by operator right-division, recurses on the lower-order quotient, "
        "and closes with one first-order linear solve. Failing that it seeks a "
        "first-order LEFT factor (D + s) through the adjoint, which yields an "
        "order-(n-1) right factor (Beke); that quotient is solved by recursion and "
        "closed by variation of parameters.");

    symtab_add_builtin("DSolve`DFactor", builtin_dsolve_dfactor);
    symtab_get_def("DSolve`DFactor")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`DFactor",
        "DSolve`DFactor[eqn, y, x] factors the linear differential operator of eqn into "
        "first-order factors, returning {Dx - r1, Dx - r2, ...} INNERMOST FIRST, so "
        "L = (Dx - rk) ... (Dx - r1) (Dx an inert d/dx operator). Right factors are "
        "found directly; a left factor, found through the adjoint, is emitted last. A "
        "not-fully-reducible operator returns its factors plus an inert remainder "
        "operator Sum q_j Dx^j.");
}
