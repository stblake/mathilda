/*
 * dsolve_knownsol.c — DSolve`KnownSolutionReduction.
 *
 * Reduction of order from an ELEMENTARY KNOWN SOLUTION, with the quadratures
 * held inert when they do not close.  For a second-order linear equation
 *
 *     c2(x) y'' + c1(x) y' + c0(x) y == g(x),        P = c1/c2,
 *
 * a solution y1 of the homogeneous part turns y = y1 v into a FIRST-order linear
 * equation in w = v':
 *
 *     w' + (2 y1'/y1 + P) w == g/(c2 y1),   integrating factor mu = y1^2 E^Int[P],
 *
 * so the general solution is
 *
 *     y = y1 (C[2] + C[1] Int[1/mu] + Int[(1/mu) Int[mu g/(c2 y1)]]).
 *
 * The method is the general fallback for the many corpus equations whose second
 * solution (or particular) is non-elementary but one solution is a simple closed
 * form the special-function rows and Kovacic do not reach -- Kovacic needs
 * rational coefficients with numeric data, so it misses
 *   x y'' + (2 + f[x] x) y' + f[x] y == 0     (y1 = 1/x, arbitrary f; §2.2.60-5933)
 *   x^2 y'' + a y' - b (b x^2 + a) y == 0      (y1 = E^(b x), symbolic b; §2.2.60-5968)
 * and its forcing closure; and the forced y'' - x^2 y' + x y == x (y1 = x,
 * §2.2.59-5836).
 *
 * Candidates for y1, each certified by L[y1] == 0 symbolically (Together, then
 * the zero test):
 *   - the monomials x^m, m = -3..3 (m = 0 is the constant solution of a missing-y
 *     operator, which DSolve`ReductionOfOrder already owns -- kept for c0 == 0
 *     equations it declines);
 *   - the exponentials E^(lambda x), lambda a root of c2 lambda^2 + c1 lambda + c0
 *     that vanishes IDENTICALLY in x (the coefficients of its numerator, as a
 *     polynomial in x, must vanish -- Solve over those).
 * Each quadrature is attempted bounded (ds_integrate_or_inert, 3 s) and held as
 * Inactive[Integrate] when it does not close -- the reference's shape for this
 * family.
 *
 * Correctness gate: an inert body is accepted UNCONDITIONALLY by the substrate's
 * verifiers (see ds_inert_vop_verified), so this method verifies its own answer:
 * the residual, with every inert integral replaced by an opaque symbol after
 * differentiation (the derivative of Inactive[Integrate][f, x] is f), must vanish
 * symbolically (Together + zero test) or pass ds_inert_vop_verified's
 * exact-rational sampling with the parameters instantiated.  A closed body goes
 * through dsolve_run's ordinary back-substitution as well.
 *
 * Bounds (BRIEF rules for a recursive / late method): top-level only
 * (g_dsolve_depth <= 1), re-entry guard registered with
 * tc_register_unwind_counter, every quadrature time-bounded, decline memo per
 * top-level evaluation, MATHILDA_NO_KNOWNSOL killswitch.  It runs LATE, just
 * before the M61 inert variation-of-parameters slot, so it can only turn an
 * UNEVAL into an answer.  General solutions with an inert quadrature decline for
 * an IVP (an indefinite inert integral cannot be fitted at a point).
 */
#include "dsolve_common.h"
#include "dsolve.h"               /* g_dsolve_depth */
#include "../core.h"              /* tc_register_unwind_counter */
#include "../sym_names.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../eval.h"
#include "../common.h"            /* head_is */
#include <stdlib.h>
#include <stdio.h>

static int g_ks_active = 0;       /* re-entry guard */

#define KS_MEMO_SLOTS 32
static uint64_t ks_epoch = 0;
static int ks_memo_n = 0;
static uint64_t ks_memo[KS_MEMO_SLOTS];
static void ks_memo_sync(uint64_t tid){ if (tid != ks_epoch){ ks_epoch = tid; ks_memo_n = 0; } }
static bool ks_memo_seen(uint64_t h){ for (int i = 0; i < ks_memo_n; i++) if (ks_memo[i] == h) return true; return false; }
static void ks_memo_add(uint64_t h){ if (ks_memo_n < KS_MEMO_SLOTS && !ks_memo_seen(h)) ks_memo[ks_memo_n++] = h; }

#define KS_INT_BUDGET 3

/* L[u] = c2 u'' + c1 u' + c0 u (u borrowed), evaluated. */
static Expr* ks_apply(Expr** c, const Expr* u, const char* xv) {
    Expr* d1 = ds_d(expr_copy((Expr*)u), expr_new_symbol(xv));
    Expr* d2 = ds_d(expr_copy(d1), expr_new_symbol(xv));
    return eval_and_free(ds_call2(SYM_Plus,
               ds_call2(SYM_Times, expr_copy(c[2]), d2),
               ds_call2(SYM_Plus, ds_call2(SYM_Times, expr_copy(c[1]), d1),
                                  ds_call2(SYM_Times, expr_copy(c[0]), expr_copy((Expr*)u)))));
}

/* Is L[u] identically zero?  Together first (cheap and exact on rational / arbitrary-
 * function coefficients), then the zero test. */
static bool ks_annihilates(Expr** c, const Expr* u, const char* xv) {
    Expr* r = eval_and_free(ds_call1(SYM_Together, ks_apply(c, u, xv)));
    bool z = ds_is_structural_zero(r) || ds_is_zero(r);
    expr_free(r);
    return z;
}

/* Exponential candidates: lambda with c2 l^2 + c1 l + c0 == 0 identically in x.
 * Returns up to 2 owned roots in out[]. */
static int ks_exp_rates(Expr** c, const char* xv, Expr** out) {
    const char* lv = intern_symbol("DSolve`ksL");
    Expr* E = eval_and_free(ds_call1(SYM_Together, ds_call2(SYM_Plus,
                  ds_call2(SYM_Times, expr_copy(c[2]),
                      ds_call2(SYM_Power, expr_new_symbol(lv), expr_new_integer(2))),
                  ds_call2(SYM_Plus, ds_call2(SYM_Times, expr_copy(c[1]), expr_new_symbol(lv)),
                                     expr_copy(c[0])))));
    Expr* num = eval_and_free(ds_call1(SYM_Numerator, E));
    num = eval_and_free(ds_call1("Expand", num));
    int n = 0;
    Expr* pq = eval_and_free(ds_call2(SYM_PolynomialQ, expr_copy(num), expr_new_symbol(xv)));
    bool poly = (pq->type == EXPR_SYMBOL && pq->data.symbol.name == SYM_True);
    expr_free(pq);
    if (poly && ds_contains(num, lv)) {
        Expr* cl = eval_and_free(ds_call2("CoefficientList", expr_copy(num), expr_new_symbol(xv)));
        if (head_is(cl, SYM_List)) {
            size_t m = cl->data.function.arg_count;
            Expr** eqs = malloc((m ? m : 1) * sizeof(Expr*));
            size_t ne = 0;
            bool contra = false;
            for (size_t i = 0; i < m; i++) {
                Expr* ci = cl->data.function.args[i];
                if (ds_is_structural_zero(ci)) continue;              /* 0 == 0 */
                if (!ds_contains(ci, lv)) { contra = true; continue; } /* k == 0, k != 0 */
                eqs[ne++] = ds_call2(SYM_Equal, expr_copy(ci), expr_new_integer(0));
            }
            if (contra || ne == 0) {
                for (size_t i = 0; i < ne; i++) expr_free(eqs[i]);
                free(eqs); expr_free(cl); expr_free(num);
                return 0;
            }
            Expr* sys = expr_new_function(expr_new_symbol(SYM_List), eqs, ne);
            free(eqs);
            Expr* sol = ds_solve(sys, expr_new_symbol(lv));
            size_t k = 0;
            Expr** rs = sol ? dsolve_extract_solutions(sol, lv, &k) : NULL;
            if (sol) expr_free(sol);
            for (size_t j = 0; j < k; j++) {
                if (n < 2 && ds_free_of(rs[j], xv) && !ds_is_zero(rs[j])) out[n++] = rs[j];
                else expr_free(rs[j]);
            }
            free(rs);
        }
        expr_free(cl);
    }
    expr_free(num);
    return n;
}

/* Integrate each additive term of f (bounded, inert on refusal), pulling the
 * x-free factors out of each term.  f consumed. */
/* Elementary closed form: only arithmetic, Exp/Log, the trigonometric and
 * hyperbolic functions and their inverses.  (An incomplete Gamma / ExpIntegral
 * antiderivative is exact but may not numericize at the sample points -- Gamma[2/3,
 * -x^3/3] does not -- leaving the answer uncertifiable.) */
static bool ks_elementary(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return true;
    const Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL) return false;
    static const char* ok[] = { "Plus", "Times", "Power", "Log", "Exp", "Sin", "Cos", "Tan",
        "Cot", "Sec", "Csc", "Sinh", "Cosh", "Tanh", "Coth", "Sech", "Csch", "ArcSin",
        "ArcCos", "ArcTan", "ArcCot", "ArcSinh", "ArcCosh", "ArcTanh", "ArcCoth",
        "Rational", "Complex", "Sqrt" };
    bool good = false;
    for (size_t i = 0; i < sizeof(ok)/sizeof(ok[0]) && !good; i++)
        good = (h->data.symbol.name == intern_symbol(ok[i]));
    if (!good) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!ks_elementary(e->data.function.args[i])) return false;
    return true;
}

static bool g_ks_elem_only = false;   /* second pass: hold non-elementary antiderivatives inert */

static Expr* ks_integrate(Expr* f, const char* xv) {
    Expr* fe = eval_and_free(ds_call1("Expand", f));
    size_t nt = head_is(fe, SYM_Plus) ? fe->data.function.arg_count : 1;
    Expr** parts = malloc(nt * sizeof(Expr*));
    for (size_t i = 0; i < nt; i++) {
        Expr* t = head_is(fe, SYM_Plus) ? expr_copy(fe->data.function.args[i]) : expr_copy(fe);
        Expr* k = expr_new_integer(1);
        Expr* g = t;
        if (head_is(t, SYM_Times)) {
            size_t n = t->data.function.arg_count, nk = 0, ng = 0;
            Expr** ka = malloc(n * sizeof(Expr*));
            Expr** ga = malloc(n * sizeof(Expr*));
            for (size_t j = 0; j < n; j++) {
                Expr* fj = t->data.function.args[j];
                if (!ds_contains(fj, intern_symbol(xv))) ka[nk++] = expr_copy(fj);
                else ga[ng++] = expr_copy(fj);
            }
            expr_free(k);
            k = nk ? eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), ka, nk))
                   : expr_new_integer(1);
            g = ng ? eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), ga, ng))
                   : expr_new_integer(1);
            free(ka); free(ga);
            expr_free(t);
        }
        Expr* gi = ds_integrate_or_inert(expr_copy(g), xv, KS_INT_BUDGET);
        if (g_ks_elem_only && !ds_has_inactive_integrate(gi) && !ks_elementary(gi)) {
            expr_free(gi);
            gi = ds_inactivate_integrate(expr_new_function(expr_new_symbol(SYM_Integrate),
                     (Expr*[]){ expr_copy(g), expr_new_symbol(xv) }, 2));
        }
        expr_free(g);
        parts[i] = ds_call2(SYM_Times, k, gi);
    }
    expr_free(fe);
    Expr* r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), parts, nt));
    free(parts);
    return r;
}

/* Self-verification of a body carrying inert integrals (see header).  Borrowed. */
static bool ks_verified(const DSolveProblem* P, Expr** c, const Expr* g,
                        const Expr* body, const char* xv) {
    Expr* r = eval_and_free(ds_call2(SYM_Subtract, ks_apply(c, body, xv), expr_copy((Expr*)g)));
    /* opaque symbols for the inert integrals, then an exact rational identity */
    if (ds_has_inactive_integrate(r)) {
        for (int i = 0; i < 12 && ds_has_inactive_integrate(r); i++) {
            /* find one innermost-first occurrence via the substrate's collector:
             * substituting the outermost first keeps inner ones inside it, so a
             * fixed number of passes over the distinct integrals suffices. */
            const Expr* e = r;
            /* depth-first search for an Inactive[Integrate][..] node */
            const Expr* found = NULL;
            const Expr* todo[256]; int nt = 0;
            todo[nt++] = e;
            while (nt && !found) {
                const Expr* q = todo[--nt];
                if (!q || q->type != EXPR_FUNCTION) continue;
                const Expr* h = q->data.function.head;
                if (h->type == EXPR_FUNCTION && head_is((Expr*)h, SYM_Inactive)
                    && q->data.function.arg_count == 2) { found = q; break; }
                for (size_t j = 0; j < q->data.function.arg_count && nt < 256; j++)
                    todo[nt++] = q->data.function.args[j];
            }
            if (!found) break;
            char buf[32];
            snprintf(buf, sizeof(buf), "DSolve`ksZ%d", i);
            r = ds_subst(r, expr_copy((Expr*)found), expr_new_symbol(intern_symbol(buf)));
        }
        if (ds_has_inactive_integrate(r)) { expr_free(r); return false; }
    }
    Expr* t = eval_and_free(ds_call1(SYM_Together, r));
    bool ok = ds_is_structural_zero(t) || ds_is_zero(t);
    expr_free(t);
    if (!ok) ok = ds_has_inactive_integrate(body) ? ds_inert_vop_verified(P, body)
                                                  : ds_residual_rel_verified(P, body);
    return ok;
}

/* Build the general solution from y1 (borrowed); NULL when a step fails. */
static Expr* ks_build(const DSolveProblem* P, Expr** c, const Expr* g, const Expr* y1,
                      const char* xv) {
    Expr* Pc = eval_and_free(ds_call1(SYM_Together,
                   ds_call2(SYM_Times, expr_copy(c[1]),
                       ds_call2(SYM_Power, expr_copy(c[2]), expr_new_integer(-1)))));
    Expr* IP = ks_integrate(Pc, xv);                                  /* Int[P] */
    /* mu = y1^2 E^Int[P]: simplify the CLOSED part of the exponent with y1^2 and
     * keep the inert part as a separate factor (x^2 E^(-2 Log[x] - Z) -> E^-Z). */
    Expr* closed_terms[64]; Expr* inert_terms[64]; size_t nc = 0, ni = 0;
    if (head_is(IP, SYM_Plus) && IP->data.function.arg_count <= 64) {
        for (size_t i = 0; i < IP->data.function.arg_count; i++) {
            Expr* t = expr_copy(IP->data.function.args[i]);
            if (ds_has_inactive_integrate(t)) inert_terms[ni++] = t; else closed_terms[nc++] = t;
        }
    } else if (ds_has_inactive_integrate(IP)) inert_terms[ni++] = expr_copy(IP);
    else closed_terms[nc++] = expr_copy(IP);
    expr_free(IP);
    Expr* ec = nc ? eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), closed_terms, nc))
                  : expr_new_integer(0);
    Expr* mu = eval_and_free(ds_call2(SYM_Times,
                   ds_call2(SYM_Power, expr_copy((Expr*)y1), expr_new_integer(2)),
                   ds_call1("Exp", ec)));
    mu = eval_and_free(ds_call1(SYM_PowerExpand, ds_simplify(mu)));
    if (ni) {
        Expr* ei = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), inert_terms, ni));
        mu = eval_and_free(ds_call2(SYM_Times, mu, ds_call1("Exp", ei)));
    }
    Expr* wh = eval_and_free(ds_call2(SYM_Power, expr_copy(mu), expr_new_integer(-1)));
    if (!ds_has_inactive_integrate(wh)) wh = ds_simplify(wh);
    Expr* v = ds_call2(SYM_Plus, ds_const(2), ds_call2(SYM_Times, ds_const(1), ks_integrate(wh, xv)));
    if (!ds_is_structural_zero(g)) {
        Expr* src = eval_and_free(ds_call2(SYM_Times, expr_copy(mu),
                        ds_call2(SYM_Times, expr_copy((Expr*)g),
                            ds_call2(SYM_Power,
                                ds_call2(SYM_Times, expr_copy(c[2]), expr_copy((Expr*)y1)),
                                expr_new_integer(-1)))));
        if (!ds_has_inactive_integrate(src)) src = ds_simplify(src);
        Expr* inner = ks_integrate(src, xv);
        Expr* wp = eval_and_free(ds_call2(SYM_Times, inner,
                       ds_call2(SYM_Power, expr_copy(mu), expr_new_integer(-1))));
        if (!ds_has_inactive_integrate(wp)) wp = ds_simplify(wp);
        v = ds_call2(SYM_Plus, v, ks_integrate(wp, xv));
    }
    expr_free(mu);
    Expr* body = eval_and_free(ds_call2(SYM_Times, expr_copy((Expr*)y1), v));
    body = eval_and_free(ds_call1("Expand", body));
    if (ds_has_active_integrate(body)) body = ds_inactivate_integrate(body);
    if (ds_has_inactive_integrate(body) && P->ncond > 0) { expr_free(body); return NULL; }
    if (!ks_verified(P, c, g, body, xv)) { expr_free(body); return NULL; }
    return body;
}

/* ks_build, then (if the closed form could not be certified) again with every
 * non-elementary antiderivative held inert. */
static Expr* ks_build2(const DSolveProblem* P, Expr** c, const Expr* g, const Expr* y1,
                       const char* xv) {
    Expr* b = ks_build(P, c, g, y1, xv);
    if (b) return b;
    g_ks_elem_only = true;
    b = ks_build(P, c, g, y1, xv);
    g_ks_elem_only = false;
    return b;
}

Expr** dsolve_knownsol_try(DSolveProblem* P, size_t* nbranch) {
    if (getenv("MATHILDA_NO_KNOWNSOL")) return NULL;
    if (g_dsolve_depth > 1 || g_ks_active) return NULL;
    if (P->nfun != 1 || P->neq != 1 || P->max_order[0] != 2) return NULL;
    const char* xv = P->ind_names[0];

    ks_memo_sync(eval_toplevel_id());
    uint64_t h = expr_hash(P->eq_residuals[0]);
    if (ks_memo_seen(h)) return NULL;

    Expr** c = NULL; Expr* g = NULL; int n = 0;
    if (!dsolve_linear_coeffs(P, &c, &g, &n)) { ks_memo_add(h); return NULL; }
    if (n != 2 || ds_is_structural_zero(c[2])) {
        for (int k = 0; k <= n; k++) expr_free(c[k]);
        free(c); expr_free(g); ks_memo_add(h); return NULL;
    }

    g_ks_active++;
    Expr* body = NULL;
    /* monomial candidates, small |m| first */
    static const int ms[] = { 1, -1, 2, -2, 3, -3, 0 };
    for (size_t i = 0; i < sizeof(ms)/sizeof(ms[0]) && !body; i++) {
        Expr* y1 = (ms[i] == 0) ? expr_new_integer(1)
                 : eval_and_free(ds_call2(SYM_Power, expr_new_symbol(xv), expr_new_integer(ms[i])));
        if (ks_annihilates(c, y1, xv)) body = ks_build2(P, c, g, y1, xv);
        expr_free(y1);
    }
    if (!body) {
        Expr* lam[2] = { NULL, NULL };
        int nl = ks_exp_rates(c, xv, lam);
        for (int i = 0; i < nl; i++) {
            if (!body) {
                Expr* y1 = eval_and_free(ds_call1("Exp",
                               ds_call2(SYM_Times, expr_copy(lam[i]), expr_new_symbol(xv))));
                if (ks_annihilates(c, y1, xv)) body = ks_build2(P, c, g, y1, xv);
                expr_free(y1);
            }
            expr_free(lam[i]);
        }
    }
    g_ks_active--;

    for (int k = 0; k <= n; k++) expr_free(c[k]);
    free(c); expr_free(g);
    if (!body) { ks_memo_add(h); return NULL; }
    Expr** out = malloc(sizeof(Expr*));
    out[0] = body;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_knownsol(Expr* res) {
    return dsolve_method_builtin(res, dsolve_knownsol_try);
}

void dsolve_knownsol_init(void) {
    tc_register_unwind_counter(&g_ks_active);   /* restored on a timeout unwind */
    symtab_add_builtin("DSolve`KnownSolutionReduction", builtin_dsolve_knownsol);
    symtab_get_def("DSolve`KnownSolutionReduction")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`KnownSolutionReduction",
        "DSolve`KnownSolutionReduction[eqn, y, x] solves a second-order linear ODE "
        "c2 y'' + c1 y' + c0 y == g with an elementary solution y1 of its homogeneous "
        "part (a monomial x^m, |m| <= 3, or an exponential E^(lambda x)) by reduction "
        "of order: y = y1 (C[2] + C[1] Int[1/mu] + Int[Int[mu g/(c2 y1)]/mu]), "
        "mu = y1^2 E^Int[c1/c2]; quadratures that do not close are kept as "
        "Inactive[Integrate] and the answer is self-verified.");
}
