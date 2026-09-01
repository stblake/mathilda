/*
 * dsolve_operator_factor.c — DSolve`OperatorFactor and DSolve`DFactor.
 *
 * Factors a linear differential operator  L = Sum_{k=0}^{n} a_k(x) D^k  by finding
 * a FIRST-ORDER RIGHT FACTOR  (D - r),  r in C(x).  (D - r) right-divides L iff
 * y = Exp[Integrate[r]] solves L[y] == 0 (a hyperexponential solution), iff the
 * "Riccati" residual  R(r) = Sum a_k P_k(r)  vanishes, where the P_k are the Bell
 * polynomials  P_0 = 1,  P_{k+1} = P_k' + r P_k  (so P_k(r) = D^k(Exp[Int r])/Exp[Int r]).
 *
 * DSolve`OperatorFactor (a scalar cascade method, homogeneous order >= 3): find one
 * rational r, peel it via operator right-division to the order-(n-1) quotient
 * L1 = L/(D-r), recurse DSolve on L1[z]==0, then close with the first-order linear
 * solve (D - r)y == z.  This reuses the whole cascade (the quotient may be solved by
 * const-coeff / Euler / Kovacic / a further OperatorFactor peel), exactly like
 * ExactODE / ReductionOfOrder.  Order 2 is left to Kovacic (which runs first and whose
 * Case 1 is the same rational first-order-factor search).
 *
 * DSolve`DFactor[eqn, y, x]: the standalone factoriser.  Returns the ordered factor
 * list  {Dx - r1, Dx - r2, ...}  (Dx an inert d/dx operator symbol; the factors
 * compose left-to-right, L = (Dx-r1)(Dx-r2)...).  A not-fully-reducible operator
 * returns its peeled first-order factors plus the inert remainder operator Sum q_j Dx^j.
 *
 * Right-division recurrence (monic a_n = 1):  q_{n-1} = a_n;
 *   q_{m-1} = a_m + Sum_{j=m}^{n-1} C(j,m) D^(j-m)[r] q_j   (m = n-1..1);
 * remainder rho = a_0 + Sum_{j} q_j r^(j)  ==  R(r)  (the factor divides iff rho == 0).
 *
 * Self-contained: reuses only the shared ds_* substrate and front-end builtins; makes
 * NO changes to dsolve_kovacic.c (the small ansatz/solve overlap is duplicated to keep
 * that engine at zero regression risk).
 *
 * Scope (first cut): linear, homogeneous, order >= 3 (OperatorFactor) / >= 1 (DFactor);
 * first-order right factors with r in C(x) (simple poles + low-degree polynomial part),
 * rational coefficients with numeric constants.  Deferred: inhomogeneous forcing,
 * irregular-singular (double-pole) r, irreducible-quadratic-denominator r, and
 * 2nd-order right factors (Beke).
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
static Expr* of_find_factor(Expr** a, int n, const char* x) {
    Expr* sum = expr_copy(a[0]);
    for (int k = 1; k <= n; k++) sum = A2(sum, expr_copy(a[k]));
    Expr* factors = fn1("FactorList", fn1("Denominator", fn1("Together", sum)));

    Expr* found = NULL;
    for (int pd = 0; pd <= 2 && !found; pd++) {
        for (int po = 1; po <= 2 && !found; po++) {
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

/* Monic coefficient array a[0..n] from the parsed problem, or NULL.  Sets *order. */
static Expr** of_monic_coeffs(DSolveProblem* P, int* order, bool* homog) {
    Expr** c; Expr* g; int n;
    if (!dsolve_linear_coeffs(P, &c, &g, &n)) return NULL;
    *homog = ds_is_zero(g);
    expr_free(g);
    const char* x = P->ind_names[0];
    if (n < 1 || ds_is_zero(c[n]) || !of_coeffs_numeric(c, n, x)) {
        for (int k = 0; k <= n; k++) expr_free(c[k]);
        free(c);
        return NULL;
    }
    Expr** a = malloc((size_t)(n + 1) * sizeof(Expr*));
    for (int k = 0; k <= n; k++)
        a[k] = ds_simplify(T2(expr_copy(c[k]), Powi(expr_copy(c[n]), -1)));
    for (int k = 0; k <= n; k++) expr_free(c[k]);
    free(c);
    *order = n;
    return a;
}

/* ---------- DSolve`OperatorFactor ---------- */

Expr** dsolve_operfactor_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    const char* xvar = P->ind_names[0];

    int n; bool homog;
    Expr** a = of_monic_coeffs(P, &n, &homog);
    if (!a) return NULL;
    if (n < 3 || !homog) { for (int k = 0; k <= n; k++) expr_free(a[k]); free(a); return NULL; }

    Expr* r = of_find_factor(a, n, xvar);
    Expr* body = NULL;
    if (r) {
        Expr* rho; Expr** q = of_divide(a, n, r, xvar, &rho); expr_free(rho);

        /* reduced ODE  Sum q[j] z^(j) == 0  in a fresh function ofz */
        const char* zf = intern_symbol("DSolve`ofz");
        Expr** terms = malloc((size_t)n * sizeof(Expr*));
        for (int j = 0; j < n; j++)
            terms[j] = ds_call2(SYM_Times, expr_copy(q[j]), ds_make_funcapp(zf, j, xvar));
        Expr* lhs = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)n));
        free(terms);
        for (int j = 0; j < n; j++) expr_free(q[j]);
        free(q);
        Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                       (Expr*[]){ lhs, expr_new_integer(0) }, 2);
        Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                         (Expr*[]){ eq, ds_make_funcapp(zf, 0, xvar),
                                    expr_new_symbol(xvar) }, 3);
        Expr* res = eval_and_free(call);
        Expr* z = extract_applied(res, zf);
        expr_free(res);

        if (z) {
            /* trailing first-order solve  (D - r)y == z:
             *   y = Exp[Int r] ( Int z Exp[-Int r] dx + C[kk] ),  kk = max C-index(z)+1. */
            Expr* intR = ds_integrate(expr_copy(r), expr_new_symbol(xvar));
            bool okR = !ds_has_head(intR, SYM_Integrate);
            if (okR) {
                Expr* chk = Sub(ds_d(expr_copy(intR), expr_new_symbol(xvar)), expr_copy(r));
                okR = ds_is_zero(chk);
                expr_free(chk);
            }
            if (okR) {
                Expr* w    = fn1("Exp", expr_copy(intR));
                Expr* winv = fn1("Exp", Sub(expr_new_integer(0), expr_copy(intR)));
                Expr* part = ds_integrate(T2(expr_copy(z), winv), expr_new_symbol(xvar));
                if (!ds_has_head(part, SYM_Integrate)) {
                    int mx = 0; of_scan_const(z, &mx);
                    body = T2(w, A2(part, ds_const(mx + 1)));
                } else {
                    expr_free(part); expr_free(w);
                }
            }
            expr_free(intR);
            expr_free(z);
        }
    }

    for (int k = 0; k <= n; k++) expr_free(a[k]);
    free(a);
    if (r) expr_free(r);
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

    int n; bool homog; (void)homog;
    Expr** a = of_monic_coeffs(&P, &n, &homog);
    dsolve_problem_free(&P);
    if (!a) return NULL;

    size_t cap = 8, nf = 0;
    Expr** facs = malloc(cap * sizeof(Expr*));
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
        if (!r) break;                                   /* irreducible tail of order cur */
        if (nf == cap) { cap *= 2; facs = realloc(facs, cap * sizeof(Expr*)); }
        facs[nf++] = dfactor_first_order(r);
        Expr* rho; Expr** q = of_divide(a, cur, r, xvar, &rho); expr_free(rho);
        expr_free(r);
        for (int k = 0; k <= cur; k++) expr_free(a[k]);
        free(a);
        a = q;                                           /* order cur-1, monic q[cur-1]=1 */
        cur -= 1;
    }
    if (a && cur >= 2) {                                 /* append irreducible remainder */
        if (nf == cap) { cap *= 2; facs = realloc(facs, cap * sizeof(Expr*)); }
        facs[nf++] = dfactor_remainder(a, cur);
    }
    if (a) { for (int k = 0; k <= cur; k++) expr_free(a[k]); free(a); }

    Expr* out = expr_new_function(expr_new_symbol(SYM_List), facs, nf);
    free(facs);
    return out;
}

void dsolve_operfactor_init(void) {
    symtab_add_builtin("DSolve`OperatorFactor", builtin_dsolve_operfactor);
    symtab_get_def("DSolve`OperatorFactor")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`OperatorFactor",
        "DSolve`OperatorFactor[eqn, y, x] solves a homogeneous linear ODE of order "
        ">= 3 by factoring its operator: it finds a first-order right factor (D - r) "
        "with r rational (a hyperexponential solution Exp[Integrate[r]]), peels it off "
        "(operator right-division), recurses on the lower-order quotient, and closes "
        "with one first-order linear solve.");

    symtab_add_builtin("DSolve`DFactor", builtin_dsolve_dfactor);
    symtab_get_def("DSolve`DFactor")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`DFactor",
        "DSolve`DFactor[eqn, y, x] factors the linear differential operator of eqn into "
        "first-order right factors, returning {Dx - r1, Dx - r2, ...} (Dx an inert d/dx "
        "operator; the factors compose left-to-right).  A not-fully-reducible operator "
        "returns its first-order factors plus an inert remainder operator Sum q_j Dx^j.");
}
