/* risch_field_integrate.c — recursive transcendental-tower integrator.
 *
 * The Bronstein field integrator over the differential tower (risch_tower.c):
 * canonical split + Hermite reduction + LRT log part + coupled bounded ansatz +
 * field Risch DE, mutually recursive with the RDE via rt_field_integrate.  Also
 * the flat Phase-B tower cases (rt_log_tower_case / rt_exp_tower_case), the
 * top-level recursive assembly (rt_recursive_tower_case), and the shared
 * elementary-integrability decision state.  See risch_field_integrate.h.
 */

#include "risch_field_integrate.h"
#include "integrate_risch_rde.h"
#include "risch_tower.h"
#include "risch_util.h"

#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "attr.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "arithmetic.h"
#include "flint_bridge.h"
#include "risch_field.h"
#include "risch_hermite.h"
#include "risch_canonical.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/* Multiplicity of the variable v at 0 in the polynomial p: the lowest power of v
 * that appears with a nonzero coefficient (0 if v does not divide p).  This is the
 * exact negative Laurent extent for an EXPONENTIAL kernel v — the order of the pole
 * that v contributes to a proper fraction — derived from the input, no cap. */
static long rt_var_mult_at_zero(Expr* p, Expr* v) {
    long a = 0;
    Expr* cl = rt_eval2("CoefficientList", expr_copy(p), expr_copy(v));
    if (cl && cl->type == EXPR_FUNCTION
        && cl->data.function.head->type == EXPR_SYMBOL
        && cl->data.function.head->data.symbol.name == intern_symbol("List")) {
        for (size_t i = 0; i < cl->data.function.arg_count; i++)
            if (!rt_is_zero(cl->data.function.args[i])) { a = (long)i; break; }
    }
    if (cl) expr_free(cl);
    return a;
}
/* ================================================================== */
/* The tower cases (rt_log_tower_case / rt_exp_tower_case) are bounded-ansatz
 * SEARCHES over a multi-variable SolveAlways, not pure decision procedures: the
 * certificate "D_tower[Q] = F as a polynomial identity in the tower variables"
 * is sound only when SolveAlways returns a genuine solution, and a non-elementary
 * integrand can yield a spurious one (e.g. E^(E^x)/(1+E^(E^x)) is non-elementary
 * yet the ansatz once returned Log[1+E^(E^x)]/E^x).  As with the other
 * search-based integrators, the tower results are therefore diff-back verified:
 * Simplify[D[result,x] - f] must be exactly 0, else the case declines.  This
 * keeps the genuine closures (E^x E^(E^x) -> E^(E^x),
 * E^x E^(E^x)/(1+E^(E^x)) -> Log[1+E^(E^x)]) and rejects the spurious ones. */
/* ================================================================== */
/* Case: nested logarithmic tower (Phase B, second increment).         */
/* ================================================================== */

/* Nested logarithmic tower.  When the integrand is a rational function of a
 * chain of NESTED logarithms  t_1 = Log[u_1](x), t_2 = Log[u_2](x, t_1), ...,
 * t_n = Log[u_n](x, t_1..t_{n-1})  (n >= 2), the single-extension cases (which
 * model one kernel over C(x)) do not apply.  Generalize the single-kernel
 * derivation  D = d/dx + Dt d/dt  to the full tower
 *   D_tower = d/dx + sum_i Dt_i d/dt_i,   Dt_i = Cancel[D[u_i,x]/u_i] |_{ker->t}
 * (triangular: Dt_i lies in C(x, t_1..t_{i-1})), and solve ONE unified ansatz
 *   Q = sum_{k=0}^{Ntop} P_k(x, t_1..t_{n-1}) t_n^k  +  sum_j c_j Log(g_j)
 * (P_k bounded-degree polynomials with unknown constant coefficients, g_j the
 * squarefree t_n-denominator factors, c_j constants) by requiring
 * D_tower[Q] = f for all {t_n,...,t_1,x} via SolveAlways.  Correct by
 * construction (an exact solution certifies the identity).  Closes e.g.
 * Integrate[1/(x Log[x] Log[Log[x]]), x] = Log[Log[Log[x]]] and
 * Integrate[Log[Log[x]]/(x Log[x]), x] = Log[Log[x]]^2/2.  Declines cleanly
 * when the tower is not a valid nested chain, the top denominator has a
 * repeated pole (tower Hermite — later), or a lower-field coefficient must be
 * rational (the full recursion — later). */
Expr* rt_log_tower_case(Expr* f, Expr* x) {
    /* Normalize composite logs (Log[a b] -> Log a + Log b, Log[b^p] -> p Log b)
     * so the tower is built over the MINIMAL independent generator set.  `fn`
     * drives kernel collection and the ansatz; correctness is still gated by
     * diff-back against the ORIGINAL `f` (rt_verify_antideriv, below). */
    Expr* fn = rt_eval_own(rt_expand_logs(f));
    if (!fn) fn = expr_copy(f);

    Expr** logs = NULL; size_t nl = 0, lcap = 0;
    rt_collect_logs(fn, x, &logs, &nl, &lcap);
    if (nl < 2 || nl > 4) { for (size_t i = 0; i < nl; i++) expr_free(logs[i]);
                            free(logs); expr_free(fn); return NULL; }

    /* Order innermost-first: if logs[i] contains logs[k], logs[k] is deeper. */
    for (size_t i = 0; i < nl; i++)
        for (size_t k = i + 1; k < nl; k++)
            if (rt_contains(logs[i], logs[k])) {
                Expr* t = logs[i]; logs[i] = logs[k]; logs[k] = t;
            }

    /* t-symbols and substitution rules kernel_i -> t_i. */
    Expr** ts = malloc(nl * sizeof(Expr*));
    Expr** rules = malloc(nl * sizeof(Expr*));
    for (size_t i = 0; i < nl; i++) {
        char nm[32]; snprintf(nm, sizeof(nm), "rmtw%zu", i);
        ts[i] = expr_new_symbol(nm);
        rules[i] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_copy(logs[i]), expr_copy(ts[i]) }, 2);
    }
    Expr* rl = expr_new_function(expr_new_symbol("List"), rules, nl);
    free(rules);

    Expr* F = rt_eval1("Together", rt_eval_own(expr_new_function(
        expr_new_symbol("ReplaceAll"), (Expr*[]){ expr_copy(fn), expr_copy(rl) }, 2)));
    Expr* top = ts[nl - 1];
    Expr* num = NULL; Expr* den = NULL;
    Expr** Dt = NULL;
    Expr* result = NULL;

    bool ok = F && rt_find_log_of_x(F, x) == NULL && rt_find_exp_of_x(F, x) == NULL;
    if (ok) {
        num = rt_eval1("Numerator", expr_copy(F));
        den = rt_eval1("Denominator", expr_copy(F));
        /* F must be a genuine RATIONAL function of the whole tower {x, t_1..t_n}
         * — i.e. num and den are multivariate polynomials in those variables and
         * nothing else.  Without this, a residual non-rational kernel of an inner
         * argument (e.g. Sin[Log[x]] -> Sin[t_1]) would slip through the per-top-
         * variable test and the ansatz could certify a WRONG closed form. */
        Expr** vv = malloc((nl + 1) * sizeof(Expr*));
        vv[0] = expr_copy(x);
        for (size_t i = 0; i < nl; i++) vv[i + 1] = expr_copy(ts[i]);
        Expr* vlist = expr_new_function(expr_new_symbol("List"), vv, nl + 1);
        free(vv);
        Expr* pqn = num ? rt_eval2("PolynomialQ", expr_copy(num), expr_copy(vlist)) : NULL;
        Expr* pqd = den ? rt_eval2("PolynomialQ", expr_copy(den), expr_copy(vlist)) : NULL;
        ok = num && den && rt_is_true(pqn) && rt_is_true(pqd)
             && !rt_free_of_x(F, top);
        if (pqn) expr_free(pqn);
        if (pqd) expr_free(pqd);
        expr_free(vlist);
    }

    /* Tower derivation coefficients Dt_i, checked triangular. */
    if (ok) {
        Dt = calloc(nl, sizeof(Expr*));
        for (size_t i = 0; i < nl && ok; i++) {
            Expr* arg = logs[i]->data.function.args[0];
            Expr* darg = rt_eval2("D", expr_copy(arg), expr_copy(x));
            Expr* q = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ darg, expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(arg), expr_new_integer(-1) }, 2) }, 2);
            Expr* dc = rt_eval1("Cancel", q);
            Dt[i] = dc ? rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                (Expr*[]){ dc, expr_copy(rl) }, 2)) : NULL;
            if (!Dt[i]) { ok = false; break; }
            for (size_t j = i; j < nl; j++)
                if (!rt_free_of_x(Dt[i], ts[j])) { ok = false; break; }
        }
    }

    if (ok) {
        long dtn_num = rt_degree(num, top), dtn_den = rt_degree(den, top);
        /* Exact top-degree bound: t_n is a LOGARITHM (D lowers deg_{t_n} by 1), so
         * deg_{t_n}(Q) = deg_{t_n}(f) + 1.  No cap. */
        long Ntop = dtn_num - dtn_den + 1;
        if (Ntop < 0) Ntop = 0;

        /* Lower field variables: x, t_1..t_{n-1}. */
        size_t nlv = nl;                 /* x plus (nl-1) inner kernels */
        Expr** lv = malloc(nlv * sizeof(Expr*));
        long* bd = malloc(nlv * sizeof(long));
        lv[0] = x;
        for (size_t i = 0; i + 1 < nl; i++) lv[i + 1] = ts[i];
        for (size_t j = 0; j < nlv; j++) {
            long a = rt_degree(num, lv[j]); if (a < 0) a = 0;
            long b = rt_degree(den, lv[j]); if (b < 0) b = 0;
            long d = a + b + 1;              /* derived lower-field proxy, no cap */
            bd[j] = d;
        }
        long nmono = 1;
        for (size_t j = 0; j < nlv; j++) nmono *= (bd[j] + 1);

        /* Squarefree t_n-dependent factors of den -> log terms. */
        Expr* g[16]; size_t ng = 0; bool bad = false;
        Expr* factored = rt_eval1("Factor", expr_copy(den));
        if (!factored) bad = true;
        else {
            Expr** fa; size_t nf; Expr* single[1];
            if (factored->type == EXPR_FUNCTION
                && factored->data.function.head->type == EXPR_SYMBOL
                && factored->data.function.head->data.symbol.name == intern_symbol("Times")) {
                fa = factored->data.function.args; nf = factored->data.function.arg_count;
            } else { single[0] = factored; fa = single; nf = 1; }
            for (size_t i = 0; i < nf && !bad; i++) {
                Expr* term = fa[i]; Expr* base = term; long e = 1;
                if (term->type == EXPR_FUNCTION
                    && term->data.function.head->type == EXPR_SYMBOL
                    && term->data.function.head->data.symbol.name == intern_symbol("Power")
                    && term->data.function.arg_count == 2
                    && term->data.function.args[1]->type == EXPR_INTEGER) {
                    base = term->data.function.args[0];
                    e = (long)term->data.function.args[1]->data.integer;
                }
                if (rt_free_of_x(base, top)) continue;     /* not a t_n factor */
                if (e != 1 || ng >= 16) { bad = true; break; }  /* repeated: later */
                g[ng++] = expr_copy(base);
            }
        }

        long nunk = (Ntop + 1) * nmono + (long)ng;
        if (!bad && nunk > 0) {
            /* Build Q = sum_{k,mono} rmLp{k}_{m} (mono) t_n^k + sum_j rmLc{j} Log(g_j). */
            size_t nq = (size_t)((Ntop + 1) * nmono + (long)ng);
            Expr** qterms = malloc(nq * sizeof(Expr*));
            size_t ntq = 0;
            long* ev = malloc(nlv * sizeof(long));
            for (long k = 0; k <= Ntop; k++)
                for (long m = 0; m < nmono; m++) {
                    char nm[64]; snprintf(nm, sizeof(nm), "rmLp%ld_%ld", k, m);
                    rt_decode_mono(m, bd, nlv, ev);
                    Expr* mono = rt_build_monomial(lv, ev, nlv);
                    qterms[ntq++] = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_symbol(nm), mono,
                          expr_new_function(expr_new_symbol("Power"),
                            (Expr*[]){ expr_copy(top), expr_new_integer(k) }, 2) }, 3);
                }
            for (size_t j = 0; j < ng; j++) {
                char nm[32]; snprintf(nm, sizeof(nm), "rmLc%zu", j);
                qterms[ntq++] = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_symbol(nm),
                      expr_new_function(expr_new_symbol("Log"),
                        (Expr*[]){ expr_copy(g[j]) }, 1) }, 2);
            }
            Expr* Q = expr_new_function(expr_new_symbol("Plus"), qterms, ntq);
            free(qterms);

            /* D_tower[Q] = D[Q,x] + sum_i Dt_i D[Q,t_i]. */
            Expr** dterms = malloc((nl + 1) * sizeof(Expr*));
            dterms[0] = rt_eval2("D", expr_copy(Q), expr_copy(x));
            for (size_t i = 0; i < nl; i++) {
                Expr* dqi = rt_eval2("D", expr_copy(Q), expr_copy(ts[i]));
                dterms[i + 1] = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_copy(Dt[i]), dqi }, 2);
            }
            Expr* Qder = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                dterms, nl + 1));
            free(dterms);

            Expr* diff = expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ Qder, expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(F) }, 2) }, 2);
            Expr* tog = rt_eval1("Together", diff);
            Expr* rnum = tog ? rt_eval1("Numerator", tog) : NULL;

            /* SolveAlways over {t_n, ..., t_1, x}. */
            Expr* sol = NULL;
            if (rnum) {
                Expr** vl = malloc((nl + 1) * sizeof(Expr*));
                for (size_t i = 0; i < nl; i++) vl[i] = expr_copy(ts[nl - 1 - i]);
                vl[nl] = expr_copy(x);
                Expr* varlist = expr_new_function(expr_new_symbol("List"), vl, nl + 1);
                free(vl);
                Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
                    (Expr*[]){ rnum, expr_new_integer(0) }, 2);
                sol = rt_eval2("SolveAlways", eqn, varlist);
            }
            if (sol && sol->type == EXPR_FUNCTION
                && sol->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.head->data.symbol.name == intern_symbol("List")
                && sol->data.function.arg_count >= 1
                && sol->data.function.args[0]->type == EXPR_FUNCTION
                && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.args[0]->data.function.head->data.symbol.name
                     == intern_symbol("List")) {
                Expr* Qs = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                    (Expr*[]){ expr_copy(Q), expr_copy(sol->data.function.args[0]) }, 2));
                if (Qs) {
                    /* Pin any free ansatz parameters to 0. */
                    Expr** zero = malloc(nq * sizeof(Expr*));
                    size_t zi = 0;
                    for (long k = 0; k <= Ntop; k++)
                        for (long m = 0; m < nmono; m++) {
                            char nm[64]; snprintf(nm, sizeof(nm), "rmLp%ld_%ld", k, m);
                            zero[zi++] = expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_new_symbol(nm), expr_new_integer(0) }, 2);
                        }
                    for (size_t j = 0; j < ng; j++) {
                        char nm[32]; snprintf(nm, sizeof(nm), "rmLc%zu", j);
                        zero[zi++] = expr_new_function(expr_new_symbol("Rule"),
                            (Expr*[]){ expr_new_symbol(nm), expr_new_integer(0) }, 2);
                    }
                    Expr* zl = expr_new_function(expr_new_symbol("List"), zero, zi);
                    free(zero);
                    Qs = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ Qs, zl }, 2));
                    if (Qs) {
                        /* Back-substitute t_i -> Log kernels (innermost first so a
                         * t appearing inside a Log argument is restored too). */
                        Expr** back = malloc(nl * sizeof(Expr*));
                        for (size_t i = 0; i < nl; i++)
                            back[i] = expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_copy(ts[i]), expr_copy(logs[i]) }, 2);
                        Expr* bl = expr_new_function(expr_new_symbol("List"), back, nl);
                        free(back);
                        result = rt_eval_own(expr_new_function(
                            expr_new_symbol("ReplaceAll"), (Expr*[]){ Qs, bl }, 2));
                    }
                }
            }
            if (sol) expr_free(sol);
            free(ev);
            expr_free(Q);
        }
        for (size_t j = 0; j < ng; j++) expr_free(g[j]);
        if (factored) expr_free(factored);
        free(lv); free(bd);
    }

    /* Diff-back safety gate (see rt_verify_antideriv). */
    if (result && !rt_verify_antideriv(result, f, x)) { expr_free(result); result = NULL; }

    if (Dt) { for (size_t i = 0; i < nl; i++) if (Dt[i]) expr_free(Dt[i]); free(Dt); }
    if (num) expr_free(num);
    if (den) expr_free(den);
    if (F) expr_free(F);
    expr_free(fn);
    expr_free(rl);
    for (size_t i = 0; i < nl; i++) { expr_free(ts[i]); expr_free(logs[i]); }
    free(ts); free(logs);
    return result;
}

/* Shared tower solve tail: given the ansatz Q (in {x, t_1..t_n}), its unknown
 * coefficient symbols `syms`, the tower derivation coefficients Dt, the tower
 * variables ts, the target F, and the back-substitution rule list `backrules`
 * (t_i -> kernel_i), form D_tower[Q] - F, solve the polynomial identity over
 * {t_n,...,t_1,x} with SolveAlways, pin free parameters to 0, and return the
 * back-substituted antiderivative (owned) or NULL.  Borrows all arguments.
 *
 * EXTENDED-LIOUVILLE EXTENSION POINT (Cherry substrate — CHERRY_DESIGN.md §2.2/§3.1).
 * This IS the "reduce to a linear system over the constants" solver both Cherry
 * algorithms need (Risch69 Main Thm part (b) in the 1986 paper; the undetermined-
 * coefficient E1/E2/E3 solve in the 1989 paper).  It already accepts special-
 * function basis terms with NO change: put  Q = poly_ansatz + Sum_i k_i SF(a_i)
 * (SF in {LogIntegral, ExpIntegralEi, Erf, PolyLog, ...}, a_i a generated
 * argument) with the k_i listed in `syms`.  The tower derivative it forms,
 * D[Q,x] + Sum_j Dt_j D[Q,t_j], differentiates the SF terms correctly because
 * Mathilda's D[] knows their derivatives (D[LogIntegral[u]] = u'/Log[u], etc.), so
 * for a kernel argument u the SF term contributes exactly D_tower[u]/Log[u] to the
 * linear system.  The Cherry driver (future) generates the a_i, appends the SF
 * basis terms to Q and the k_i to `syms`, and calls this unchanged. */
static Expr* rt_tower_solve(Expr* Q, Expr** syms, size_t nsym,
                            Expr** Dt, Expr** ts, size_t nl,
                            Expr* F, Expr* x, Expr* backrules) {
    Expr** dterms = malloc((nl + 1) * sizeof(Expr*));
    dterms[0] = rt_eval2("D", expr_copy(Q), expr_copy(x));
    for (size_t i = 0; i < nl; i++) {
        Expr* dqi = rt_eval2("D", expr_copy(Q), expr_copy(ts[i]));
        dterms[i + 1] = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(Dt[i]), dqi }, 2);
    }
    Expr* Qder = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
        dterms, nl + 1));
    free(dterms);
    Expr* diff = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ Qder, expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(F) }, 2) }, 2);
    Expr* tog = rt_eval1("Together", diff);
    Expr* rnum = tog ? rt_eval1("Numerator", tog) : NULL;

    Expr* result = NULL; Expr* sol = NULL;
    if (rnum) {
        Expr** vl = malloc((nl + 1) * sizeof(Expr*));
        for (size_t i = 0; i < nl; i++) vl[i] = expr_copy(ts[nl - 1 - i]);
        vl[nl] = expr_copy(x);
        Expr* varlist = expr_new_function(expr_new_symbol("List"), vl, nl + 1);
        free(vl);
        Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
            (Expr*[]){ rnum, expr_new_integer(0) }, 2);
        sol = rt_eval2("SolveAlways", eqn, varlist);
    }
    if (sol && sol->type == EXPR_FUNCTION
        && sol->data.function.head->type == EXPR_SYMBOL
        && sol->data.function.head->data.symbol.name == intern_symbol("List")
        && sol->data.function.arg_count >= 1
        && sol->data.function.args[0]->type == EXPR_FUNCTION
        && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
        && sol->data.function.args[0]->data.function.head->data.symbol.name
             == intern_symbol("List")) {
        Expr* Qs = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
            (Expr*[]){ expr_copy(Q), expr_copy(sol->data.function.args[0]) }, 2));
        if (Qs) {
            Expr** zero = malloc((nsym ? nsym : 1) * sizeof(Expr*));
            for (size_t j = 0; j < nsym; j++)
                zero[j] = expr_new_function(expr_new_symbol("Rule"),
                    (Expr*[]){ expr_copy(syms[j]), expr_new_integer(0) }, 2);
            Expr* zl = expr_new_function(expr_new_symbol("List"), zero, nsym);
            free(zero);
            Qs = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                (Expr*[]){ Qs, zl }, 2));
            if (Qs)
                result = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                    (Expr*[]){ Qs, expr_copy(backrules) }, 2));
        }
    }
    if (sol) expr_free(sol);
    return result;
}

/* Nested EXPONENTIAL tower (Phase B, third increment) — the dual of
 * rt_log_tower_case for a chain of nested exponentials t_i = E^(u_i) with u_i in
 * C(x, t_1..t_{i-1}) (e.g. t_1 = E^x, t_2 = E^(E^x)).  An exponential kernel's
 * derivative is a multiple of itself, so the tower derivation is
 *   D_tower = d/dx + sum_i Dt_i d/dt_i,   Dt_i = (D[u_i,x] |_{ker->t}) t_i,
 * and a unified Laurent ansatz
 *   Q = sum_{i=ilo}^{ihi} P_i(x, t_1..t_{n-1}) t_n^i  +  sum_j c_j Log(g_j)
 * (P_i bounded-degree lower-field polynomials, ilo = -(t_n multiplicity in den),
 * g_j the squarefree t_n-coprime denominator factors) is solved by SolveAlways
 * over {t_n,...,t_1,x}.  Same whole-tower rationality gate and correct-by-
 * construction certificate as the log tower.  Closes e.g.
 * Integrate[E^x E^(E^x), x] = E^(E^x).  Repeated t_n poles (tower Hermite) and
 * rational lower-field coefficients are out of scope and decline. */
Expr* rt_exp_tower_case(Expr* f, Expr* x) {
    if (!rt_find_exp_of_x(f, x)) return NULL;
    Expr** us = NULL; size_t nl = 0, ucap = 0;
    rt_collect_exp_exponents(f, x, &us, &nl, &ucap);   /* us[i] = exp exponents */
    if (nl < 2 || nl > 4) { for (size_t i = 0; i < nl; i++) expr_free(us[i]);
                            free(us); return NULL; }

    /* Kernel forms E^(u_i). */
    Expr** kf = malloc(nl * sizeof(Expr*));
    for (size_t i = 0; i < nl; i++)
        kf[i] = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_new_symbol("E"), expr_copy(us[i]) }, 2);

    /* Order innermost-first: if kf[i] contains kf[k], kf[k] is deeper. */
    for (size_t i = 0; i < nl; i++)
        for (size_t k = i + 1; k < nl; k++)
            if (rt_contains(kf[i], kf[k])) {
                Expr* t = kf[i]; kf[i] = kf[k]; kf[k] = t;
                t = us[i]; us[i] = us[k]; us[k] = t;
            }

    /* t-symbols and substitution rules (Exp[u]->t and Power[E,u]->t). */
    Expr** ts = malloc(nl * sizeof(Expr*));
    Expr** rules = malloc(2 * nl * sizeof(Expr*));
    for (size_t i = 0; i < nl; i++) {
        char nm[32]; snprintf(nm, sizeof(nm), "rmte%zu", i);
        ts[i] = expr_new_symbol(nm);
        rules[2 * i] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Exp"),
                (Expr*[]){ expr_copy(us[i]) }, 1), expr_copy(ts[i]) }, 2);
        rules[2 * i + 1] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_copy(kf[i]), expr_copy(ts[i]) }, 2);
    }
    Expr* rl = expr_new_function(expr_new_symbol("List"), rules, 2 * nl);
    free(rules);

    Expr* F = rt_eval1("Together", rt_eval_own(expr_new_function(
        expr_new_symbol("ReplaceAll"), (Expr*[]){ expr_copy(f), expr_copy(rl) }, 2)));
    Expr* top = ts[nl - 1];
    Expr* num = NULL; Expr* den = NULL; Expr** Dt = NULL; Expr* result = NULL;

    bool ok = F && rt_find_log_of_x(F, x) == NULL && rt_find_exp_of_x(F, x) == NULL;
    if (ok) {
        num = rt_eval1("Numerator", expr_copy(F));
        den = rt_eval1("Denominator", expr_copy(F));
        Expr** vv = malloc((nl + 1) * sizeof(Expr*));
        vv[0] = expr_copy(x);
        for (size_t i = 0; i < nl; i++) vv[i + 1] = expr_copy(ts[i]);
        Expr* vlist = expr_new_function(expr_new_symbol("List"), vv, nl + 1);
        free(vv);
        Expr* pqn = num ? rt_eval2("PolynomialQ", expr_copy(num), expr_copy(vlist)) : NULL;
        Expr* pqd = den ? rt_eval2("PolynomialQ", expr_copy(den), expr_copy(vlist)) : NULL;
        ok = num && den && rt_is_true(pqn) && rt_is_true(pqd)
             && !rt_free_of_x(F, top);
        if (pqn) expr_free(pqn);
        if (pqd) expr_free(pqd);
        expr_free(vlist);
    }

    /* Dt_i = Cancel[D[u_i,x]] |_{ker->t} * t_i. */
    if (ok) {
        Dt = calloc(nl, sizeof(Expr*));
        for (size_t i = 0; i < nl && ok; i++) {
            Expr* dui = rt_eval2("D", expr_copy(us[i]), expr_copy(x));
            Expr* dc = rt_eval1("Cancel", dui);
            Expr* dcs = dc ? rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                (Expr*[]){ dc, expr_copy(rl) }, 2)) : NULL;
            Dt[i] = dcs ? rt_eval_own(expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ dcs, expr_copy(ts[i]) }, 2)) : NULL;
            if (!Dt[i]) { ok = false; break; }
        }
    }

    if (ok) {
        /* a = multiplicity of t_n at 0 in den (Laurent negative extent). */
        long a = rt_var_mult_at_zero(den, top);
        Expr* Dtil = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(top), expr_new_integer(-a) }, 2) }, 2));

        /* Squarefree t_n-dependent factors of Dtil -> log terms. */
        Expr* g[16]; size_t ng = 0; bool bad = false;
        Expr* factored = Dtil ? rt_eval1("Factor", expr_copy(Dtil)) : NULL;
        if (!factored) bad = true;
        else {
            Expr** fa; size_t nf; Expr* single[1];
            if (factored->type == EXPR_FUNCTION
                && factored->data.function.head->type == EXPR_SYMBOL
                && factored->data.function.head->data.symbol.name == intern_symbol("Times")) {
                fa = factored->data.function.args; nf = factored->data.function.arg_count;
            } else { single[0] = factored; fa = single; nf = 1; }
            for (size_t i = 0; i < nf && !bad; i++) {
                Expr* term = fa[i]; Expr* base = term; long e = 1;
                if (term->type == EXPR_FUNCTION
                    && term->data.function.head->type == EXPR_SYMBOL
                    && term->data.function.head->data.symbol.name == intern_symbol("Power")
                    && term->data.function.arg_count == 2
                    && term->data.function.args[1]->type == EXPR_INTEGER) {
                    base = term->data.function.args[0];
                    e = (long)term->data.function.args[1]->data.integer;
                }
                if (rt_free_of_x(base, top)) continue;    /* t_n-coprime constant part */
                if (e != 1 || ng >= 16) { bad = true; break; }  /* repeated: later */
                g[ng++] = expr_copy(base);
            }
        }

        long dtn_num = rt_degree(num, top), dtn_den = rt_degree(den, top);
        /* Top kernel t_n is an EXPONENTIAL (D preserves deg_{t_n}), so Q's Laurent
         * range in t_n matches f's: high extent deg_{t_n}(num) - deg_{t_n}(den),
         * low extent -(multiplicity of t_n at 0 in den).  Exact, no cap. */
        long ihi = dtn_num - dtn_den;
        long ilo = -a;

        size_t nlv = nl;                 /* x plus (nl-1) inner kernels */
        Expr** lv = malloc(nlv * sizeof(Expr*));
        long* bd = malloc(nlv * sizeof(long));   /* per-var odometer range size-1 */
        long* lo = malloc(nlv * sizeof(long));   /* per-var lowest exponent       */
        lv[0] = x; lo[0] = 0;                    /* x: polynomial, degree >= 0     */
        for (size_t i = 0; i + 1 < nl; i++) lv[i + 1] = ts[i];
        {
            /* x is polynomial; degree bound derived from num/den, no cap. */
            long p = rt_degree(num, x); if (p < 0) p = 0;
            long q = rt_degree(den, x); if (q < 0) q = 0;
            long d = p + q + 1;
            bd[0] = d;
        }
        /* Inner kernels t_1..t_{n-1} are EXPONENTIALS and invertible, so their
         * coefficient exponents are LAURENT — the antiderivative can carry a power
         * the integrand lacks (e.g. INT E^(x+E^x) dx = E^(E^x)).  Each inner window
         * is f's own per-kernel Laurent extent WIDENED by the reach of the top
         * derivation coefficient w' = D[u_n] (Dt[top] = w' t_n): one factor of w'
         * enters via t_n^i in D_tower[Q], shifting inner-kernel exponents by up to
         * deg_{t_j}(w') upward and w''s t_j-pole order downward.  Fully derived from
         * the input — no hardcoded window. */
        Expr* wtop = rt_eval1("Together", expr_copy(Dt[nl - 1]));
        Expr* wnum = wtop ? rt_eval1("Numerator", expr_copy(wtop)) : NULL;
        Expr* wden = wtop ? rt_eval1("Denominator", expr_copy(wtop)) : NULL;
        for (size_t j = 1; j < nlv; j++) {
            long hij = rt_degree(num, ts[j - 1]) - rt_degree(den, ts[j - 1]);
            if (hij < 0) hij = 0;
            long loj = -rt_var_mult_at_zero(den, ts[j - 1]);
            long wr = 0;
            if (wnum && wden) {
                long wu = rt_degree(wnum, ts[j - 1]); if (wu < 0) wu = 0;
                wr = wu + rt_var_mult_at_zero(wden, ts[j - 1]);
            }
            lo[j] = loj - wr;
            bd[j] = (hij + wr) - (loj - wr);       /* odometer count-1 (>= 0) */
        }
        if (wtop) expr_free(wtop);
        if (wnum) expr_free(wnum);
        if (wden) expr_free(wden);
        long nmono = 1;
        for (size_t j = 0; j < nlv; j++) nmono *= (bd[j] + 1);
        long nwi = (ihi >= ilo) ? (ihi - ilo + 1) : 0;
        long nunk = nwi * nmono + (long)ng;

        if (!bad && nwi > 0 && nunk > 0) {
            size_t nq = (size_t)nunk;
            Expr** qterms = malloc(nq * sizeof(Expr*));
            Expr** syms = malloc(nq * sizeof(Expr*));
            size_t ntq = 0, nsym = 0;
            long* ev = malloc(nlv * sizeof(long));
            for (long i = ilo; i <= ihi; i++)
                for (long m = 0; m < nmono; m++) {
                    char nm[40]; snprintf(nm, sizeof(nm), "rmEp%ld_%ld", i - ilo, m);
                    rt_decode_mono(m, bd, nlv, ev);
                    for (size_t j = 0; j < nlv; j++) ev[j] += lo[j];  /* Laurent shift */
                    Expr* mono = rt_build_monomial(lv, ev, nlv);
                    syms[nsym++] = expr_new_symbol(nm);
                    qterms[ntq++] = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_symbol(nm), mono,
                          expr_new_function(expr_new_symbol("Power"),
                            (Expr*[]){ expr_copy(top), expr_new_integer(i) }, 2) }, 3);
                }
            for (size_t j = 0; j < ng; j++) {
                char nm[40]; snprintf(nm, sizeof(nm), "rmEc%zu", j);
                syms[nsym++] = expr_new_symbol(nm);
                qterms[ntq++] = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_symbol(nm),
                      expr_new_function(expr_new_symbol("Log"),
                        (Expr*[]){ expr_copy(g[j]) }, 1) }, 2);
            }
            free(ev);
            Expr* Q = expr_new_function(expr_new_symbol("Plus"), qterms, ntq);
            free(qterms);

            Expr** back = malloc(nl * sizeof(Expr*));
            for (size_t i = 0; i < nl; i++)
                back[i] = expr_new_function(expr_new_symbol("Rule"),
                    (Expr*[]){ expr_copy(ts[i]), expr_copy(kf[i]) }, 2);
            Expr* bl = expr_new_function(expr_new_symbol("List"), back, nl);
            free(back);

            result = rt_tower_solve(Q, syms, nsym, Dt, ts, nl, F, x, bl);
            expr_free(bl);
            for (size_t j = 0; j < nsym; j++) expr_free(syms[j]);
            free(syms);
            expr_free(Q);
        }
        for (size_t j = 0; j < ng; j++) expr_free(g[j]);
        if (factored) expr_free(factored);
        if (Dtil) expr_free(Dtil);
        free(lv); free(bd); free(lo);
    }

    /* Diff-back safety gate (see rt_verify_antideriv). */
    if (result && !rt_verify_antideriv(result, f, x)) { expr_free(result); result = NULL; }

    if (Dt) { for (size_t i = 0; i < nl; i++) if (Dt[i]) expr_free(Dt[i]); free(Dt); }
    if (num) expr_free(num);
    if (den) expr_free(den);
    if (F) expr_free(F);
    expr_free(rl);
    for (size_t i = 0; i < nl; i++) { expr_free(ts[i]); expr_free(kf[i]); expr_free(us[i]); }
    free(ts); free(kf); free(us);
    return result;
}

/* ================================================================== */
/* Genuine one-extension-at-a-time recursive Risch engine.             */
/* ================================================================== */
/* The tower cases above use a FLAT ansatz — one SolveAlways over every tower
 * variable at once, with bounded-degree POLYNOMIAL lower-field coefficients.
 * That misses two whole families the recursive algorithm (Bronstein ch. 5)
 * handles: (1) MIXED exp/log towers (the flat cases are each
 * single-kind), and (2) RATIONAL lower-field coefficients (a t_n coefficient
 * that is 1/x, not a polynomial — a nonlinear unknown for the flat ansatz).
 *
 * This engine peels ONE extension at a time.  To integrate F in the top field
 * K_n = K_{n-1}(t_n): split F into its polynomial/Laurent part in t_n and a
 * proper rational part; integrate the polynomial part coefficient-by-coefficient,
 * where each coefficient integral is itself an integration in the LOWER field
 * K_{n-1} (the recursion), bottoming out at the rational base case C(x)
 * (Integrate`BronsteinRational).  Because the coefficients are integrated in
 * their own field they may be rational, and because each level dispatches on its
 * own kernel kind the tower may mix logs and exponentials.
 *
 * The genuine algorithm's proper-rational part (tower Hermite + Rothstein-Trager
 * over K_{n-1}) and its general field Risch-DE are deferred to a later increment;
 * a nonzero proper part or a non-base RDE declines cleanly here.  As with the
 * flat tower cases the whole result is a bounded search, so it is diff-back
 * VERIFIED (rt_verify_antideriv) — a spurious certificate cannot ship. */



/* ================================================================== */
/* Elementary-integrability decision (P3: Bronstein §5.6 residue        */
/* criterion + Ch.6 Risch-DE no-solution certificates).                 */
/* ================================================================== */
/* RtDecision + the g_rt_decide_mode / g_rt_decision externs are declared in
 * risch_field_integrate.h; the storage is defined below. */

/* Decision context.  The field integrator (rt_field_integrate and its
 * callees) is a genuine decision procedure — every NULL it returns is either a
 * theorem-backed "no elementary integral" or a dispatch-to-sibling decline.  When
 * a caller wants the VERDICT (not just the antiderivative) it sets
 * g_rt_decide_mode and runs the same path; the sub-algorithms raise
 * g_rt_decision to RT_DEC_NONELEMENTARY at their AUTHORITATIVE decline points
 * ONLY — a non-constant residue (Thm 5.6.1(ii)), a Risch DE with no rational
 * solution (Ch.6, via the complete rde_base / exact-degree-bound ansatz), or the
 * hypertangent Dc≠0 certificate — and NEVER on a routing decline.  The flag is
 * write-once/sticky within a run; the driver (rt_decide) resets it, then reads
 * ELEMENTARY iff a full antiderivative with no unintegrated remainder is produced
 * (an elementary result always wins over a stray flag from an abandoned path). */
bool       g_rt_decide_mode = false;
RtDecision g_rt_decision     = RT_DEC_UNKNOWN;

/* Raise the decision to NONELEMENTARY at an authoritative decline (no-op unless a
 * decision is being requested). */
static void rt_dec_nonelem(void) {
    if (g_rt_decide_mode) g_rt_decision = RT_DEC_NONELEMENTARY;
}






/* Mutually recursive field-integration primitives (forward declarations).
 * `rem_out` (when non-NULL) receives the unintegrated non-constant-residue
 * remainder of a partial log part (Bronstein Thm 5.6.1); it is threaded only
 * from the top-level recursive-tower assembly point — inner recursions pass
 * NULL, which makes the residue criterion decline a partial rather than surface
 * one mid-recursion. */
Expr* rt_field_integrate(Expr* F, RtTower* T, long L, Expr* x, Expr** rem_out);
static int   rt_limited_field_integrate(Expr* r, RtTower* T, long L, Expr* x, int is_bottom,
                                        Expr** s_out, Expr** c_out);
static Expr* rt_int_primitive_poly(Expr* num, Expr* den, RtTower* T, long L, Expr* x);
static Expr* rt_int_hyperexp_poly(Expr* num, Expr* den, RtTower* T, long L, Expr* x);
static Expr* rt_field_rde(Expr* p, long i, RtTower* T, long L, Expr* x);
static Expr* rt_field_ratint(Expr* num, Expr* den, RtTower* T, long L, Expr* x, Expr** rem_out);
static Expr* rt_field_hyperexp_coupled(Expr* num, Expr* den, RtTower* T, long L, Expr* x,
                                       Expr** rem_out);



/* Pure resultant Lazard-Rioboo-Trager log part for a TOWER proper rational
 * part Rr/den in t = t_L (deg_t Rr < deg_t den), lifting rt_frac_lrt's single-
 * kernel algebraic-residue closure to the recursive-Risch tower.  The bounded
 * SolveAlways ansatz in rt_field_ratint / rt_field_hyperexp_coupled expresses
 * exactly one CONSTANT residue per squarefree t-factor, so it cannot integrate
 * a factor whose Rothstein-Trager residues are ALGEBRAIC (e.g. the +-I/2
 * residues that split t^2+1 into a conjugate log pair = ArcTan) — this closes
 * exactly that class at the tower level, e.g.
 *   1/(x Log[x] (Log[Log[x]]^2+1)) -> ArcTan[Log[Log[x]]].
 *
 * Only the squarefree case (den has no repeated t-factor) is handled: den
 * factors as cont(x, t_0..t_{L-1}) * rad with rad squarefree in t, and the log
 * part is
 *   Integrate`TranscendentalLogPart[Rr/cont, rad, t, z, D_tower[rad],
 *                                   {x, t_0, ..., t_{L-1}}],
 * whose gate list requires the residues to be constants of the tower derivation
 * (free of every lower-field variable).  A genuine Hermite (repeated-pole) part
 * is left to the coupled ansatz.  Returns the log part in tower-variable form
 * (owned, in terms of the t_L symbol), or NULL.  The top-level caller diff-back
 * verifies, so an incomplete reduction (e.g. an exponential Laurent part not
 * captured here) declines rather than shipping a wrong form. */
static Expr* rt_field_lrt_logpart(Expr* Rr, Expr* den, RtTower* T, long L,
                                  Expr* x, Expr** rem_out) {
    if (rem_out) *rem_out = NULL;
    Expr* t = T->t[L];
    long dn = rt_degree(Rr, t), dd = rt_degree(den, t);
    if (dn < 0 || dd <= 0 || dn >= dd) return NULL;   /* proper, den depends on t */

    /* Strip the t-free content of den: Hden = gcd(den, dden/dt); when den is
     * squarefree in t this is exactly the t-free factor (dH == 0). */
    Expr* dden = rt_eval2("D", expr_copy(den), expr_copy(t));
    if (!dden) return NULL;
    Expr* Hden = rt_eval_call("PolynomialGCD",
        (Expr*[]){ expr_copy(den), dden }, 2);
    if (!Hden) return NULL;
    if (rt_degree(Hden, t) != 0) { expr_free(Hden); return NULL; }  /* repeated pole */

    /* rad = den / Hden (squarefree in t); a = Rr / Hden. */
    Expr* rad = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(Hden), expr_new_integer(-1) }, 2) }, 2));
    Expr* a = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(Rr), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(Hden), expr_new_integer(-1) }, 2) }, 2));
    expr_free(Hden);
    if (!rad || !a) { if (rad) expr_free(rad); if (a) expr_free(a); return NULL; }

    /* Monomial (tower) derivation of the squarefree radical. */
    Expr* Dd = rt_tower_deriv(rad, T, x);
    if (!Dd) { expr_free(rad); expr_free(a); return NULL; }

    /* Clear lower-field (x, t_0..t_{L-1}) denominators of a and Dd so the resultant
     * machinery sees polynomial coefficients: scaling BOTH by the same t_L-free
     * factor leaves the residues z = a/Dd and the log arguments gcd(a - z Dd, rad)
     * unchanged.  This is essential when a lower monomial is NONLINEAR — a tangent
     * t_j has Dt_j = t_j^2 + 1, so its logarithmic-derivative coefficient
     * Dcoef = (t_j^2+1)/t_j is rational, giving a, Dd a t_j denominator. */
    Expr* ad = rt_eval1("Denominator", rt_eval1("Together", expr_copy(a)));
    Expr* Ddd = rt_eval1("Denominator", rt_eval1("Together", expr_copy(Dd)));
    Expr* lcd = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ ad, Ddd }, 2));                       /* adopts ad, Ddd */
    if (lcd && !rt_is_zero(lcd) && rt_free_of_x(lcd, t)) {
        Expr* a2 = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ a, expr_copy(lcd) }, 2));          /* adopts a */
        Expr* D2 = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ Dd, lcd }, 2));                    /* adopts Dd, lcd */
        a = a2; Dd = D2;
    } else if (lcd) expr_free(lcd);

    /* Gate list {x, t_0, ..., t_{L-1}}: residues must be constants of D_tower. */
    size_t ng = (size_t)L + 1;
    Expr** gv = malloc(ng * sizeof(Expr*));
    gv[0] = expr_copy(x);
    for (long i = 0; i < L; i++) gv[i + 1] = expr_copy(T->t[i]);
    Expr* gate = expr_new_function(expr_new_symbol("List"), gv, ng);
    free(gv);

    /* In decision mode, append a 7th "decide" marker so a non-constant-residue
     * resultant returns the inert marker Integrate`$NonConstantResidue (the simple
     * part is non-elementary, Thm 5.6.1(ii)) instead of declining silently. */
    Expr* logpart = g_rt_decide_mode
        ? rt_eval_call("Integrate`TranscendentalLogPart",
              (Expr*[]){ a, rad, expr_copy(t), expr_new_symbol("rmZ"), Dd, gate,
                         expr_new_symbol("True") }, 7)
        : rt_eval_call("Integrate`TranscendentalLogPart",
              (Expr*[]){ a, rad, expr_copy(t), expr_new_symbol("rmZ"), Dd, gate }, 6);
    /* Non-constant-residue verdict marker (decision mode): the simple part has no
     * elementary integral — raise the decision and decline the (empty) log part. */
    if (logpart && logpart->type == EXPR_SYMBOL
        && logpart->data.symbol.name == intern_symbol("Integrate`$NonConstantResidue")) {
        expr_free(logpart);
        rt_dec_nonelem();
        return NULL;
    }
    /* Declined iff the head is unchanged (a, rad, Dd, gate adopted by the call). */
    if (logpart && logpart->type == EXPR_FUNCTION
        && logpart->data.function.head->type == EXPR_SYMBOL
        && logpart->data.function.head->data.symbol.name
             == intern_symbol("Integrate`TranscendentalLogPart")) {
        expr_free(logpart);
        return NULL;
    }
    /* Partial log part (Bronstein Thm 5.6.1): the builtin returned
     * Integrate`PartialLogPart[logs, remainder] — the residue resultant mixed
     * constant and non-constant residues.  Only a caller that can carry the
     * remainder (rem_out != NULL) may accept it; otherwise decline (the
     * exponential coupled path cannot absorb a proper simple remainder). */
    if (logpart && logpart->type == EXPR_FUNCTION
        && logpart->data.function.head->type == EXPR_SYMBOL
        && logpart->data.function.head->data.symbol.name
             == intern_symbol("Integrate`PartialLogPart")
        && logpart->data.function.arg_count == 2) {
        if (!rem_out) { expr_free(logpart); return NULL; }
        Expr* logs = expr_copy(logpart->data.function.args[0]);
        *rem_out = expr_copy(logpart->data.function.args[1]);
        expr_free(logpart);
        return logs;   /* tower-variable elementary logs; remainder via rem_out */
    }
    return logpart;   /* tower-variable form (in terms of the t_L symbol) */
}

/* Proper rational part in t = t_L (deg_t num < deg_t den, coefficients in the
 * lower field K_{L-1}): integrate to  Q = H(t)/Hden(t) + sum_j c_j Log(g_j),
 * where Hden = gcd(den, d den/dt) is the repeated part, g_j the distinct
 * t-dependent factors of the squarefree radical rad = den/Hden (Rothstein-Trager
 * log arguments), H a polynomial in t of degree < deg(Hden) with bounded-degree
 * lower-field (x, t_1..t_{L-1}) polynomial coefficients, and c_j constants.  All
 * unknowns are solved at once by requiring D_tower[Q] = num/den for every tower
 * variable via SolveAlways — the constant-residue and denominator-structure
 * certificates are exactly Hermite reduction + the residue criterion, lifted to
 * the tower derivation.  (A genuinely rational Hermite numerator coefficient and
 * the pure resultant LRT remain a later refinement; a bounded ansatz that misses
 * them declines, and the caller diff-back verifies.)  Returns Q (tower-variable
 * form, owned) or NULL. */

/* Bronstein canonical representation (§3.5) at tower level L: split
 * F = f_p + f_s + f_n into its polynomial part f_p in t = t_L, the special part
 * f_s = b/d_s (special denominator — poles fixed by the derivation, e.g. t^k for
 * an exponential, t^2+1 for a hypertangent), and the normal part f_n = c/d_n
 * (proper, denominator coprime to the special part).  Reuses the tested
 * differential-field C API (risch_canonical_representation) over the tower's
 * monomial derivation.  Writes owned *f_p, *f_s, *f_n; returns false (leaving
 * them untouched) only if the derivation rules are malformed. */
static bool rt_canonical_split(Expr* F, RtTower* T, long L, Expr* x,
                               Expr** f_p, Expr** f_s, Expr** f_n) {
    Expr* rules = rt_build_deriv_rules(T, x);
    RischDeriv d;
    if (!risch_deriv_from_rules(rules, &d)) { expr_free(rules); return false; }
    risch_canonical_representation(F, T->t[L], &d, f_p, f_s, f_n);
    risch_deriv_free(&d);
    expr_free(rules);
    return true;
}

/* Literal Hermite reduction of the proper part Rr/den in t = t_L (Bronstein
 * §5.3) — the exact algorithm, with no undetermined-coefficient ansatz.  For a
 * primitive (log) monomial t_L is normal and Rr/den is proper, so the reduced
 * part r is 0: Rr/den = D_tower[g] + h with h simple, and the antiderivative is
 * g + (log part of h), the log part coming from the residue criterion (§5.6).
 * The literal reduction handles arbitrary rational K_{L-1} = C(x, t_0..t_{L-1})
 * numerator coefficients.  Returns g + logpart, self-verified by a tower
 * diff-back, or NULL when there is no elementary integral. */
static Expr* rt_field_ratint_hermite(Expr* num, Expr* den, RtTower* T, long L, Expr* x,
                                     Expr** rem_out) {
    if (rem_out) *rem_out = NULL;
    Expr* t = T->t[L];
    Expr* f = rt_eval1("Cancel", rt_eval1("Together",
        expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(num), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(den), expr_new_integer(-1) }, 2) }, 2)));
    Expr* rules = rt_build_deriv_rules(T, x);
    RischDeriv d;
    if (!risch_deriv_from_rules(rules, &d)) { expr_free(rules); expr_free(f); return NULL; }
    Expr* g = NULL; Expr* h = NULL; Expr* r = NULL;
    bool ok = risch_hermite_reduce(f, t, &d, &g, &h, &r);
    risch_deriv_free(&d);
    expr_free(rules); expr_free(f);
    if (!ok) return NULL;

    /* Proper primitive input: the reduced part must vanish (no polynomial/special
     * leftover).  If not (unexpected), decline. */
    if (!rt_is_zero(r)) { expr_free(g); expr_free(h); expr_free(r); return NULL; }
    expr_free(r);

    Expr* logpart = NULL;
    Expr* rem_tower = NULL;   /* partial log part: unintegrated r_n remainder */
    if (!rt_is_zero(h)) {
        /* h is simple (squarefree normal denominator): its integral is the
         * residue-criterion / Rothstein-Trager log part.  With a partial log
         * part (Thm 5.6.1) it returns the constant-residue logs plus rem_tower. */
        Expr* hn = rt_eval1("Numerator", rt_eval1("Together", expr_copy(h)));
        Expr* hd = rt_eval1("Denominator", rt_eval1("Together", expr_copy(h)));
        logpart = (hn && hd) ? rt_field_lrt_logpart(hn, hd, T, L, x, rem_out ? &rem_tower : NULL) : NULL;
        if (hn) expr_free(hn);
        if (hd) expr_free(hd);
        if (!logpart) { expr_free(g); expr_free(h); return NULL; }  /* not elementary */
    }
    expr_free(h);

    Expr* result = logpart
        ? rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
              (Expr*[]){ g, logpart }, 2))
        : g;

    /* Self-verify the (possibly partial) identity D_tower[result] == num/den −
     * rem_tower.  Exact and decidable in the tower field; identical to the full
     * check when rem_tower is NULL.  A mismatch declines rather than shipping a
     * wrong form. */
    Expr* Dres = rt_tower_deriv(result, T, x);
    Expr* tgt = expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(num), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(den), expr_new_integer(-1) }, 2) }, 2);
    Expr* diff = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ Dres, expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), tgt }, 2) }, 2);
    if (rem_tower)   /* D[result] − num/den + rem_tower == 0 */
        diff = expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ diff, expr_copy(rem_tower) }, 2);
    Expr* chk = rt_eval1("Together", diff);
    bool good = chk && rt_is_zero(chk);
    if (chk) expr_free(chk);
    if (!good) { expr_free(result); if (rem_tower) expr_free(rem_tower); return NULL; }
    if (rem_out) *rem_out = rem_tower;
    else if (rem_tower) expr_free(rem_tower);
    return result;
}

static Expr* rt_field_ratint(Expr* num, Expr* den, RtTower* T, long L, Expr* x,
                             Expr** rem_out) {
    /* The proper rational part in t = t_L, integrated by the LITERAL Bronstein
     * pipeline: Hermite reduction (Sect. 5.3) peels the repeated normal poles
     * into an exact derivative D_tower[g], leaving a simple part h whose integral
     * is the residue-criterion / Lazard-Rioboo-Trager log part (Sect. 5.6).  This
     * is the exact decision procedure; the former bounded undetermined-
     * coefficient ansatz has been removed.  Returns g + logpart (self-verified by
     * a tower diff-back), or NULL when there is no elementary integral.  A mixed
     * residue resultant returns the constant-residue logs and the r_n remainder
     * via rem_out (Thm 5.6.1) when rem_out != NULL. */
    return rt_field_ratint_hermite(num, den, T, L, x, rem_out);
}

/* Coupled hyperexponential proper part (exponential top level t = t_L = E^(w_L)),
 * integrated by the LITERAL Bronstein pipeline (Sect. 5.3 + 5.6 + 5.9):
 *
 *   1. HermiteReduce(f)  ->  (g, h, r)              [Sect. 5.3, over D_tower]
 *      g = exact rational part (repeated NORMAL poles peeled, with arbitrary
 *      rational lower-field coefficients); h = simple normal part (squarefree
 *      denominator, coprime to t); r = reduced Laurent polynomial in t (the
 *      special part f_s at t = 0 plus the polynomial part f_p — CanonicalRepresen-
 *      tation puts the t-power denominator into f_s for an exponential monomial).
 *   2. residue criterion on h  ->  L = sum_j c_j Log(g_j)          [Sect. 5.6]
 *      constant residues c_j (rt_field_lrt_logpart gate); the g_j are coprime to t.
 *   3. leftover  P = h + r - D_tower[L].  For an exponential monomial
 *      D_tower[Log g] = D_tower[g]/g is NOT proper (deg_t D[g] = deg_t g since
 *      D[t] = w' t), so the logs of step 2 spill a t-polynomial into the Laurent
 *      part — Bronstein's coupling.  Subtracting D_tower[L] reconciles it exactly,
 *      leaving P a genuine Laurent polynomial in t (denominator a power of t).
 *   4. IntegrateHyperexponentialPolynomial(P)  ->  Q               [Sect. 5.9]
 *      per-coefficient Risch DE (rt_int_hyperexp_poly): each t^i coefficient (i!=0)
 *      solves q_i' + i w' q_i = p_i, the t^0 coefficient recurses to the lower field.
 *   return g + L + Q, self-verified by a tower diff-back.
 *
 * This is the exact decision procedure — Hermite reduction + residue criterion +
 * hyperexponential-polynomial integration — and is strictly more complete than the
 * former unified SolveAlways ansatz (which used bounded-degree POLYNOMIAL Hermite,
 * Laurent and log coefficients, declining rational lower-field coefficients).  The
 * ansatz has been removed.  Returns g + L + Q (tower-variable form, owned), or NULL
 * when there is no elementary integral. */
static Expr* rt_field_hyperexp_hermite(Expr* num, Expr* den, RtTower* T, long L, Expr* x,
                                       Expr** rem_out) {
    if (rem_out) *rem_out = NULL;
    Expr* t = T->t[L];
    Expr* f = rt_eval1("Cancel", rt_eval1("Together",
        expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(num), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(den), expr_new_integer(-1) }, 2) }, 2)));
    if (!f) return NULL;
    Expr* rules = rt_build_deriv_rules(T, x);
    RischDeriv d;
    if (!risch_deriv_from_rules(rules, &d)) { expr_free(rules); expr_free(f); return NULL; }
    Expr* g = NULL; Expr* h = NULL; Expr* r = NULL;
    bool ok = risch_hermite_reduce(f, t, &d, &g, &h, &r);
    risch_deriv_free(&d);
    expr_free(rules); expr_free(f);
    if (!ok) return NULL;

    /* Step 2: residue-criterion log part of the simple h.  When the caller can carry
     * a remainder (rem_out != NULL) we ACCEPT a partial log part (Bronstein Thm
     * 5.6.1): a mixed residue resultant yields the constant-residue logs L_s plus the
     * non-constant-residue remainder h_n (a proper simple fraction in t).  The
     * elementary part is g + L_s + Q_s, with Q_s reconciled from h_s = h - h_n, and
     * h_n reported unintegrated — mirroring the primitive/log path.  With rem_out ==
     * NULL a mixed resultant declines (rt_field_lrt_logpart returns NULL). */
    Expr* Llog = NULL;
    Expr* h_n = NULL;   /* non-constant-residue remainder (owned), NULL if none */
    if (!rt_is_zero(h)) {
        Expr* hn = rt_eval1("Numerator", rt_eval1("Together", expr_copy(h)));
        Expr* hd = rt_eval1("Denominator", rt_eval1("Together", expr_copy(h)));
        Llog = (hn && hd) ? rt_field_lrt_logpart(hn, hd, T, L, x, rem_out ? &h_n : NULL) : NULL;
        if (hn) expr_free(hn);
        if (hd) expr_free(hd);
        if (!Llog) { expr_free(g); expr_free(h); expr_free(r); return NULL; }  /* non-elementary */
    } else {
        Llog = expr_new_integer(0);
    }

    /* h_s = h - h_n (the constant-residue part that reconciles into the Laurent
     * polynomial); h_s == h when there is no partial remainder.  h is consumed. */
    Expr* h_s;
    if (h_n) {
        h_s = rt_eval1("Together", expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ h, expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), expr_copy(h_n) }, 2) }, 2));   /* adopts h */
        if (!h_s) { expr_free(g); expr_free(r); expr_free(Llog); expr_free(h_n); return NULL; }
    } else {
        h_s = h;
    }

    /* Step 3: Laurent leftover  P = h_s + r - D_tower[L]  (a Laurent polynomial in t).
     * (h_s, r are consumed by the Plus below regardless of success.) */
    Expr* DL = rt_tower_deriv(Llog, T, x);
    if (!DL) {
        expr_free(g); expr_free(h_s); expr_free(r); expr_free(Llog);
        if (h_n) expr_free(h_n);
        return NULL;
    }
    Expr* P = rt_eval1("Cancel", rt_eval1("Together",
        expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ h_s, r,
              expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), DL }, 2) }, 3)));
    if (!P) { expr_free(g); expr_free(Llog); if (h_n) expr_free(h_n); return NULL; }

    /* Step 4: integrate the Laurent polynomial P (per-coefficient Risch DE).  If P
     * has any t-coprime denominator factor (the logs did not reconcile), the
     * hyperexponential-polynomial gate declines and so do we. */
    Expr* Pn = rt_eval1("Numerator", expr_copy(P));
    Expr* Pd = rt_eval1("Denominator", expr_copy(P));
    expr_free(P);
    Expr* Q = (Pn && Pd) ? rt_int_hyperexp_poly(Pn, Pd, T, L, x) : NULL;
    if (Pn) expr_free(Pn);
    if (Pd) expr_free(Pd);
    if (!Q) { expr_free(g); expr_free(Llog); if (h_n) expr_free(h_n); return NULL; }

    Expr* result = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ g, Llog, Q }, 3));

    /* Self-verify the (possibly partial) identity D_tower[result] + h_n == num/den.
     * Exact and decidable in the tower field; with h_n == NULL this is the usual
     * D_tower[result] == num/den.  A mismatch declines (a wrong h_n reconstruction —
     * the κ_D gcd reconstruction is not exact over a transcendental monomial — can
     * therefore only decline, never ship a wrong partial). */
    Expr* Dres = rt_tower_deriv(result, T, x);
    Expr* tgt = expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(num), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(den), expr_new_integer(-1) }, 2) }, 2);
    Expr* diff = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ Dres, expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), tgt }, 2) }, 2);
    if (h_n) {
        diff = expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ diff, expr_copy(h_n) }, 2);
    }
    Expr* chk = rt_eval1("Together", diff);
    bool good = chk && rt_is_zero(chk);
    if (chk) expr_free(chk);
    if (!good) { expr_free(result); if (h_n) expr_free(h_n); return NULL; }
    if (rem_out) *rem_out = h_n;
    else if (h_n) expr_free(h_n);
    return result;
}

static Expr* rt_field_hyperexp_coupled(Expr* num, Expr* den, RtTower* T, long L, Expr* x,
                                       Expr** rem_out) {
    /* The exponential-top proper part, integrated by the literal Bronstein pipeline
     * (Hermite reduction + residue criterion + hyperexponential-polynomial), with the
     * former undetermined-coefficient ansatz removed.  See rt_field_hyperexp_hermite.
     * A mixed residue resultant surfaces its non-constant-residue part via rem_out. */
    return rt_field_hyperexp_hermite(num, den, T, L, x, rem_out);
}

/* Hypertangent top monomial (Bronstein §5.10): Fg is rational in t = t_L over the
 * lower tower field K_{L-1}, with t_L a hypertangent (Dt = Dcoef*(t^2+sigma),
 * sigma = +1 Tan/Cot, -1 Tanh/Coth).  Build the full tower derivation as a
 * Risch`Derivation rule-list {x->1, t_0->Dt_0, ..., t_L->Dt_L}, dispatch to the
 * §5.10 driver (Risch`IntegrateHypertangent for sigma=+1, Risch`IntegrateHypertanh
 * for sigma=-1), then integrate the t_L-free base-field remainder base = F - D[g]
 * recursively in K_{L-1}.  Returns the antiderivative in tower-variable form (g in
 * terms of the t_i, with the real special log Log[t_L^2+sigma]) or NULL.
 *
 * The driver's sub-steps (HermiteReduce, ResidueReduce, IntegrateHypertangent-
 * Polynomial) are tower-general — they treat every deriv variable except t_L
 * symbolically — so a polynomial-in-t_L integrand (no normal poles) is handled
 * fully over the tower.  The reduced/pole-peeling sub-step solves its base Risch
 * DE only over C(x) (Risch`RischDE with the single base var); when the true base
 * field is a deeper tower a wrong g may result, but the caller's exact
 * D_tower[Q] == F verification (rt_recursive_tower_case) rejects it — so this path
 * is SOUND (declines) even where it is not yet complete. */
static Expr* rt_int_hypertangent_field(Expr* Fg, RtTower* T, long L, Expr* x) {
    /* deriv = {x -> 1, t_0 -> Dt_0, ..., t_L -> Dt_L}, all in tower-variable form. */
    Expr** rules = malloc((size_t)(L + 2) * sizeof(Expr*));
    rules[0] = expr_new_function(expr_new_symbol("Rule"),
        (Expr*[]){ expr_copy(x), expr_new_integer(1) }, 2);
    for (long i = 0; i <= L; i++)
        rules[i + 1] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_copy(T->t[i]), rt_dt_i(T, (size_t)i) }, 2);
    Expr* deriv = expr_new_function(expr_new_symbol("List"), rules, (size_t)(L + 2));
    free(rules);

    const char* driver = (T->tsg[L] > 0) ? "Risch`IntegrateHypertangent"
                                         : "Risch`IntegrateHypertanh";
    Expr* res = rt_eval_call(driver,
        (Expr*[]){ expr_copy(Fg), expr_copy(T->t[L]), expr_copy(deriv) }, 3);
    Expr* result = NULL;
    if (res && rt_head_is(res, "List") && res->data.function.arg_count == 2
        && rt_is_true(res->data.function.args[1])) {
        Expr* g = res->data.function.args[0];    /* borrowed, tower-var form */
        /* base = F - D_tower[g]  (element of K_{L-1}: must be free of t_L).  Expand
         * after Together: the derivation returns unexpanded products whose t-terms
         * only cancel once expanded (cf. rt_hypertan_family). */
        Expr* Dg = rt_eval2("Risch`Derivation", expr_copy(g), expr_copy(deriv));
        Expr* base = Dg ? rt_eval1("Expand", rt_eval1("Together",
            expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ expr_copy(Fg), expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), Dg }, 2) }, 2))) : NULL;
        if (base && rt_free_of_x(base, T->t[L])) {
            Expr* ib = rt_is_zero(base) ? expr_new_integer(0)
                                        : rt_field_integrate(base, T, L - 1, x, NULL);
            if (ib)
                result = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                    (Expr*[]){ expr_copy(g), ib }, 2));
        }
        if (base) expr_free(base);
    }
    if (res) expr_free(res);
    expr_free(deriv);
    return result;
}

/* Integrate F (a rational function of x, t_1..t_L in tower-variable form) with
 * respect to the tower derivation, returning an antiderivative in tower-variable
 * form (owned) or NULL.  L < 0 is the rational base case C(x). */
Expr* rt_field_integrate(Expr* F, RtTower* T, long L, Expr* x, Expr** rem_out) {
    if (rem_out) *rem_out = NULL;
    if (!F) return NULL;
    if (L < 0) {
        Expr* r = rt_eval_call("Integrate`BronsteinRational",
            (Expr*[]){ expr_copy(F), expr_copy(x) }, 2);
        if (!r) return NULL;
        if (rt_head_is(r, "Integrate`BronsteinRational")) { expr_free(r); return NULL; }
        return r;
    }
    Expr* Fg = rt_eval1("Together", expr_copy(F));
    if (!Fg) return NULL;
    Expr* num = rt_eval1("Numerator", expr_copy(Fg));
    Expr* den = rt_eval1("Denominator", expr_copy(Fg));
    Expr* result = NULL;
    if (num && den) {
        if (T->kind[L] == RT_LOG) {
            /* Bronstein canonical representation (§3.5): F = f_p + f_s + f_n over
             * the primitive (log) monomial t_L.  A primitive monomial has no
             * non-trivial special polynomial (Dt ∈ k ⇒ every special factor is a
             * unit), so f_s ≡ 0; integrate the polynomial part f_p by the
             * primitive-polynomial recursion (§5.8) and the normal proper part
             * f_n by Hermite reduction + the Rothstein-Trager log part
             * (§5.3 + §5.6).  This unifies the former ad-hoc PolynomialQuotient/
             * Remainder denominator split (poly-quotient = f_p, remainder/den =
             * f_n for a primitive) behind the canonical decomposition. */
            Expr* f_p = NULL; Expr* f_s = NULL; Expr* f_n = NULL;
            if (rt_canonical_split(Fg, T, L, x, &f_p, &f_s, &f_n)) {
                if (rt_is_zero(f_s)) {   /* primitive: special part must vanish */
                    Expr* one = expr_new_integer(1);
                    Expr* poly_int = rt_int_primitive_poly(f_p, one, T, L, x);
                    expr_free(one);
                    Expr* prop_int;
                    if (rt_is_zero(f_n)) {
                        prop_int = expr_new_integer(0);
                    } else {
                        Expr* nn = rt_eval1("Numerator", expr_copy(f_n));
                        Expr* dd = rt_eval1("Denominator", expr_copy(f_n));
                        prop_int = (nn && dd) ? rt_field_ratint(nn, dd, T, L, x, rem_out) : NULL;
                        if (nn) expr_free(nn);
                        if (dd) expr_free(dd);
                    }
                    if (poly_int && prop_int)
                        result = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                            (Expr*[]){ poly_int, prop_int }, 2));
                    else {
                        if (poly_int) expr_free(poly_int);
                        if (prop_int) expr_free(prop_int);
                    }
                }
                expr_free(f_p); expr_free(f_s); expr_free(f_n);
            }
        } else if (T->kind[L] == RT_TAN) {
            /* Hypertangent top (§5.10): dispatch to the IntegrateHypertangent /
             * -Hypertanh driver over the tower, then integrate the base remainder
             * in K_{L-1}.  Correct-by-construction / decline (sound) — see helper. */
            result = rt_int_hypertangent_field(Fg, T, L, x);
        } else {
            /* Exponential top: the pure-Laurent recursion first (genuine RDE,
             * rational lower-field coefficients); if it declines because a proper
             * part is present, the coupled hyperexponential pipeline.  The latter
             * (rt_field_hyperexp_hermite) already performs the canonical split
             * implicitly: its HermiteReduce -> (g, h, r) yields the f_n-derivative
             * part g, the simple normal part h (= f_n), and the reduced Laurent
             * polynomial r (= f_p + f_s, the special part living at t = 0), which
             * it must keep COUPLED through the §5.9 residue reconciliation — so an
             * independent f_p/f_s/f_n split is deliberately not applied here. */
            result = rt_int_hyperexp_poly(num, den, T, L, x);
            if (!result)
                result = rt_field_hyperexp_coupled(num, den, T, L, x, rem_out);
        }
    }
    if (num) expr_free(num);
    if (den) expr_free(den);
    expr_free(Fg);
    /* A partial remainder is meaningful only alongside a successful elementary
     * part; if the overall integration failed, drop it. */
    if (!result && rem_out && *rem_out) { expr_free(*rem_out); *rem_out = NULL; }
    return result;
}

/* Primitive (logarithmic top) polynomial part.  num/den is a polynomial in
 * t = t_L (den free of t) with lower-field coefficients; Dt = Dcoef_L = u_L'/u_L.
 * Integrate P(t) = sum_i p_i t^i to Q = sum_i q_i t^i by the recursive coefficient
 * matching q_i' + (i+1) q_{i+1} Dt = p_i (top-down): each residual r_i is
 * integrated in the LOWER field K_L (recursion), and a would-be new logarithm
 * equal to t is folded back into q_{i+1}. */
static Expr* rt_int_primitive_poly(Expr* num, Expr* den, RtTower* T, long L, Expr* x) {
    Expr* t = T->t[L];
    Expr* Dt = T->Dcoef[L];
    long m = rt_degree(num, t);
    if (m < 0) return expr_new_integer(0);    /* num == 0: integral is 0 */
    Expr** p = malloc((size_t)(m + 1) * sizeof(Expr*));
    for (long i = 0; i <= m; i++)
        p[i] = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ rt_coeff(num, t, i), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(den), expr_new_integer(-1) }, 2) }, 2));
    Expr** q = calloc((size_t)(m + 2), sizeof(Expr*));
    bool fail = false;
    for (long i = m; i >= 0 && !fail; i--) {
        Expr* r_i;
        if (!q[i + 1]) {
            r_i = expr_copy(p[i]);
        } else {
            Expr* term = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(i + 1), expr_copy(q[i + 1]),
                          expr_copy(Dt) }, 3);
            r_i = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ expr_copy(p[i]),
                    expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_integer(-1), term }, 2) }, 2));
        }
        Expr* s = NULL; Expr* c = NULL;
        int rc = r_i ? rt_limited_field_integrate(r_i, T, L, x, i == 0, &s, &c) : -1;
        if (r_i) expr_free(r_i);
        if (rc != 0) { fail = true; break; }
        if (c) {
            Expr* bump = rt_eval_own(expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ c, expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_new_integer(i + 1), expr_new_integer(-1) }, 2) }, 2));
            if (!q[i + 1]) q[i + 1] = bump;
            else q[i + 1] = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ q[i + 1], bump }, 2));
        }
        q[i] = s;
    }
    Expr* result = NULL;
    if (!fail) {
        Expr** terms = malloc((size_t)(m + 2) * sizeof(Expr*));
        size_t nt = 0;
        for (long i = 0; i <= m + 1; i++) {
            if (!q[i] || rt_is_zero(q[i])) continue;
            Expr* pw = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(t), expr_new_integer(i) }, 2);
            terms[nt++] = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(q[i]), pw }, 2);
        }
        result = rt_eval_own(expr_new_function(expr_new_symbol("Plus"), terms, nt));
        free(terms);
    }
    for (long i = 0; i <= m; i++) if (p[i]) expr_free(p[i]);
    for (long i = 0; i <= m + 1; i++) if (q[i]) expr_free(q[i]);
    free(p); free(q);
    return result;
}

/* True if `e` contains an inverse-tangent / -cotangent (circular or hyperbolic)
 * kernel with an x-dependent argument.  Such a term is what a rational integral
 * over an IRREDUCIBLE-QUADRATIC denominator produces (LogToReal of the complex
 * conjugate log pair).  Appearing at a theta^i coefficient with i >= 1 it is the
 * dilogarithm obstruction: e.g. Integrate[Log[1+x^2]/(1+x^2), x] is not
 * elementary, so a primitive polynomial whose (i>=1)-level residual integrates to
 * such an ArcTan has NO elementary integral (Bronstein Thm 5.6.1(ii)). */
static bool rt_has_arctan_of_x(Expr* e, Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if ((h == intern_symbol("ArcTan") || h == intern_symbol("ArcCot")
             || h == intern_symbol("ArcTanh") || h == intern_symbol("ArcCoth"))) {
            for (size_t i = 0; i < e->data.function.arg_count; i++)
                if (!rt_free_of_x(e->data.function.args[i], x)) return true;
        }
    }
    if (rt_has_arctan_of_x(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (rt_has_arctan_of_x(e->data.function.args[i], x)) return true;
    return false;
}

/* Limited integration in the lower field: integrate r (free of t_L) in K_L and
 * split the result as s + c*t_L with s in K_L and c a constant (the single new
 * logarithm equal to t_L = Log[u_L] is folded into c; any other new logarithm or
 * non-constant t_L-coefficient means r is not integrable within this tower here,
 * so decline).  On success sets *s_out (owned, may be NULL for 0) and *c_out
 * (owned, NULL when 0) and returns 0; else -1.
 *
 * `is_bottom` selects the theta^0 coefficient, which may be ANY elementary
 * function of x (a new logarithm or an ArcTan from an irreducible-quadratic log
 * argument): the whole lower-field antiderivative is accepted wholesale.  At a
 * theta^{>=1} coefficient an x-dependent ArcTan (or other non-theta transcendental
 * the fold cannot absorb) is instead the dilogarithm obstruction and yields an
 * authoritative NON-elementary verdict. */
static int rt_limited_field_integrate(Expr* r, RtTower* T, long L, Expr* x, int is_bottom,
                                      Expr** s_out, Expr** c_out) {
    *s_out = NULL; *c_out = NULL;
    Expr* R = rt_field_integrate(r, T, L - 1, x, NULL);  /* lower field: no partial */
    if (!R) return -1;
    /* LimitedIntegrate (§7.2, m=1) recognition of the single new logarithm.
     * rt_field_integrate returns ∫r in TOWER-VARIABLE form, so a new logarithm it
     * introduces appears as Log[<u_L in t-vars>] — NOT the kernel-form Log[u_L] that
     * T->subrules maps to t_L.  When t_L is logarithmic, fold that tower-variable
     * new-log form to t_L first, so the primitive-polynomial fold-back captures it
     * as the c·t_L term (closes e.g. ∫ Log[Log[x]]/(x Log[x]) = Log[Log[x]]^2/2). */
    if (T->kind[L] == RT_LOG) {
        Expr* uL_tv = rt_subst_kernels(T->arg[L], T);          /* u_L in t-vars */
        Expr* newlog = expr_new_function(expr_new_symbol("Log"), (Expr*[]){ uL_tv }, 1);
        Expr* nlrule = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ newlog, expr_copy(T->t[L]) }, 2);
        Expr* nlist = expr_new_function(expr_new_symbol("List"), (Expr*[]){ nlrule }, 1);
        R = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
            (Expr*[]){ R, nlist }, 2));                        /* adopts R, nlist */
        if (!R) return -1;
    }
    Expr* Rs = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ R, expr_copy(T->subrules) }, 2));               /* adopts R */
    if (!Rs) return -1;
    Expr* tL = T->t[L];
    int rc = -1;
    if (is_bottom) {
        /* theta^0 coefficient: accept the complete elementary lower-field
         * antiderivative (new logs / ArcTan permitted — see the header). */
        *s_out = expr_copy(Rs); *c_out = NULL; rc = 0;
        expr_free(Rs);
        return rc;
    }
    bool foreign = rt_find_exp_of_x(Rs, x) != NULL || rt_find_log_of_x(Rs, x) != NULL
                   || rt_has_arctan_of_x(Rs, x);
    if (!foreign && rt_is_poly(Rs, tL) && rt_degree(Rs, tL) <= 1) {
        Expr* c = rt_coeff(Rs, tL, 1);
        Expr* s = rt_coeff(Rs, tL, 0);
        bool cconst = rt_free_of_x(c, x);
        for (size_t j = 0; j < T->n && cconst; j++)
            if (!rt_free_of_x(c, T->t[j])) cconst = false;
        if (cconst) {
            if (rt_is_zero(c)) { expr_free(c); *c_out = NULL; }
            else *c_out = c;
            *s_out = s;
            rc = 0;
        } else {
            /* The would-be new-logarithm coefficient c is NON-constant (Dc != 0):
             * by Bronstein's primitive-polynomial criterion (§5.8, p.158) the term
             * has no elementary integral.  Authoritative non-elementary verdict. */
            rt_dec_nonelem();
            expr_free(c); expr_free(s);
        }
    } else if (rt_has_arctan_of_x(Rs, x)) {
        /* A theta^{>=1} coefficient would have to carry an x-dependent ArcTan — the
         * irreducible-quadratic dilogarithm obstruction (Thm 5.6.1(ii)).  The lower
         * integral is complete (rational RHS), so this is an authoritative
         * non-elementary verdict, e.g. Integrate[(3x+1) Log[x^2+1]^5, x].  (A merely
         * unfolded x-dependent LOG is left as a plain decline — it can be a content /
         * combination artifact of a reducible argument, not proven non-elementary.) */
        rt_dec_nonelem();
    }
    expr_free(Rs);
    return rc;
}

/* Hyperexponential (exponential top) Laurent part.  num/den is a Laurent
 * polynomial in t = t_L (den = c t^a, a monomial in t) with lower-field
 * coefficients; Dcoef_L = w_L'.  Integrate sum_i p_i t^i to sum_i q_i t^i: the
 * i = 0 term integrates in the lower field K_L (recursion); each i != 0 term
 * solves the Risch differential equation q_i' + i w_L' q_i = p_i.  A genuine
 * (non-monomial) t-denominator is a proper fraction — declined here. */
static Expr* rt_int_hyperexp_poly(Expr* num, Expr* den, RtTower* T, long L, Expr* x) {
    Expr* t = T->t[L];
    long a = 0;
    Expr* cl = rt_eval2("CoefficientList", expr_copy(den), expr_copy(t));
    if (cl && cl->type == EXPR_FUNCTION
        && cl->data.function.head->type == EXPR_SYMBOL
        && cl->data.function.head->data.symbol.name == intern_symbol("List")) {
        for (size_t i = 0; i < cl->data.function.arg_count; i++)
            if (!rt_is_zero(cl->data.function.args[i])) { a = (long)i; break; }
    }
    if (cl) expr_free(cl);
    Expr* Dtil = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(t), expr_new_integer(-a) }, 2) }, 2));
    if (!Dtil || !rt_free_of_x(Dtil, t)) { if (Dtil) expr_free(Dtil); return NULL; }
    Expr* nnum = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(num), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ Dtil, expr_new_integer(-1) }, 2) }, 2));     /* adopts Dtil */
    if (!nnum) return NULL;
    long dnn = rt_degree(nnum, t);
    if (dnn < 0) { expr_free(nnum); return expr_new_integer(0); }  /* num == 0 */

    Expr** q = calloc((size_t)(dnn + 1), sizeof(Expr*));
    long* pw = malloc((size_t)(dnn + 1) * sizeof(long));
    bool fail = false;
    for (long j = 0; j <= dnn && !fail; j++) {
        long ip = j - a;
        pw[j] = ip;
        Expr* pj = rt_coeff(nnum, t, j);
        if (rt_is_zero(pj)) { expr_free(pj); q[j] = NULL; continue; }
        Expr* qj = (ip == 0) ? rt_field_integrate(pj, T, L - 1, x, NULL)
                             : rt_field_rde(pj, ip, T, L, x);
        expr_free(pj);
        if (!qj) { fail = true; break; }
        q[j] = qj;
    }
    Expr* result = NULL;
    if (!fail) {
        Expr** terms = malloc((size_t)(dnn + 1) * sizeof(Expr*));
        size_t nt = 0;
        for (long j = 0; j <= dnn; j++) {
            if (!q[j]) continue;
            Expr* term;
            if (pw[j] == 0) term = expr_copy(q[j]);
            else term = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(q[j]), expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(t), expr_new_integer(pw[j]) }, 2) }, 2);
            terms[nt++] = term;
        }
        result = rt_eval_own(expr_new_function(expr_new_symbol("Plus"), terms, nt));
        free(terms);
    }
    for (long j = 0; j <= dnn; j++) if (q[j]) expr_free(q[j]);
    free(q); free(pw);
    expr_free(nnum);
    return result;
}

/* ==================================================================== */

/* Risch differential equation q' + i w_L' q = p for q in the lower field K_L.
 * Base case (w_L and p both in C(x)): the bounded polynomial-in-x ansatz of
 * rt_solve_rde.  General case: q lives in K_L = C(x, t_0..t_{L-1}) and may be
 * rational there (a monomial denominator, e.g. 1/Log[x]), which the coupled-
 * hyperexponential polynomial coefficients cannot express — solve by a bounded
 * LAURENT ansatz over {x, t_0..t_{L-1}} (each variable ranging over negative powers
 * too), requiring D_tower[q] + i Dcoef_L q = p for every lower tower variable via
 * SolveAlways.  A non-monomial denominator (needing the full Bronstein SPDE) is out
 * of a monomial-Laurent ansatz's reach and declines. */
static Expr* rt_field_rde(Expr* p, long i, RtTower* T, long L, Expr* x) {
    Expr* w = T->arg[L];
    /* Base field C(x): the RDE q' + i w' q = p has both coefficient (i w', with
     * w = the exponent of the current monomial) and RHS p in C(x) — i.e. free of
     * every lower tower variable.  Route it to the COMPLETE base-field solver
     * rde_base (via rt_solve_rde), whose every "no solution" is authoritative
     * (Bronstein Ch. 6, A1/F5 audit).  This deliberately does NOT require w or p to
     * be POLYNOMIAL in x — rde_base handles rational coefficients (weak
     * normalization + normal-denominator reduction), so a rational RHS such as the
     * 1/x of E^x/x (∫ non-elementary, Ei) or E^(x^2) (Erf) is decided here rather
     * than diverted to the bounded-ansatz general path below (the former
     * rt_is_poly(p,x) guard was a spurious scope decline: it sent every rational-p
     * base RDE to SolveAlways). */
    bool base = true;
    for (size_t j = 0; j < T->n && base; j++)
        if (!rt_free_of_x(w, T->t[j]) || !rt_free_of_x(p, T->t[j])) base = false;
    if (base) {
        Expr* q = rt_solve_rde(p, i, w, x);
        if (!q) rt_dec_nonelem();   /* rde_base "no solution" is authoritative (Ch.6) */
        return q;
    }

    /* Gap 1/2: recursive Risch DE over the tower.  The equation is
     * D_tower[q] + (i Dcoef_L) q = p over K_L = C(x, t_0..t_{L-1}); when the top of
     * K_L (t_{L-1}) is a PRIMITIVE (RT_LOG) or EXPONENTIAL (RT_EXP) monomial, solve
     * it with the literal Bronstein Ch.6 stack (rde_tower: normal/special denominator,
     * degree bound, SPDE, non-cancellation / cancellation, antidifferentiation) rather
     * than the bounded ansatz.  The result is exact-identity certified inside
     * rde_tower, so on success it is authoritative; on NULL we fall through to the
     * ansatz (non-regressing) for the residual cases rde_tower does not yet cover
     * (the b=Dz/z limited-integration branch, tangent tops). */
    if (L >= 1 && (T->kind[L - 1] == RT_LOG || T->kind[L - 1] == RT_EXP)) {
        Expr* fco = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(i), expr_copy(T->Dcoef[L]) }, 2));
        RdeCtx C = { T, L - 1, x, T->t[L - 1] };
        Expr* q = rde_tower(fco, p, &C);
        expr_free(fco);
        if (q) return q;
        /* Gap 1e: the recursive Bronstein Chapter-6 solver (rde_tower) IS the
         * decision procedure for the field Risch DE over a log/exp tower —
         * normal/special denominator, exact degree bound, SPDE, non-cancellation,
         * cancellation (recursive), antidifferentiation.  With a log/exp lower top
         * it was actually run to completion, so its NULL is an authoritative
         * "no rational solution in K_L". */
        rt_dec_nonelem();
        return NULL;
    }

    /* The lower field's top monomial is a hypertangent (RT_TAN) or otherwise
     * outside the recursive RDE's current scope (Gap 3): the field Risch DE was
     * NOT solved here at all.  This is a SCOPE decline, NOT a theorem — so return
     * NULL WITHOUT raising the non-elementary certificate, leaving
     * ElementaryIntegralQ `undec` rather than a false `False`.  Previously an
     * unconditional rt_dec_nonelem() here declared elementary integrands such as
     * Sec[x]^2 E^Tan[x] (antiderivative E^Tan[x]) "provably non-elementary".
     * See RISCH_AUDIT_A4.md Finding 2. */
    return NULL;
}

/* Recursive-tower case: build the differential tower of f, substitute all kernels
 * to tower variables, gate on a genuine rational function of the whole tower,
 * integrate by the one-extension recursion, back-substitute, and diff-back verify.
 * Closes mixed exp/log towers and rational lower-field coefficients that the flat
 * tower cases decline; runs after them, before the trig front-end. */
/* Structural pre-pass for the recursive tower: expand a MERGED exponential
 * monomial  E^(a + b + ...)  into the product  E^a E^b ...,  built directly
 * (NOT evaluated) so the evaluator cannot re-merge it.  The evaluator
 * automatically combines a product of exponentials into one power with a
 * summed exponent (E^x E^(E^x) -> E^(x + E^x)); that merged exponent is not a
 * valid tower monomial — its argument x + E^x contains the foreign kernel E^x,
 * so rt_tower_build's structure-theorem check rejects the tower and the case
 * declines.  Splitting E^(x + E^x) back to E^x E^(E^x) restores the
 * independent tower basis {E^x, E^(E^x)}, so integrands the evaluator merged
 * (e.g. E^x E^(E^x)/(1+E^(E^x)) = E^(x+E^x)/(1+E^(E^x))) close instead.  The
 * rewrite is exact (E^(a+b) = E^a E^b) and the whole recursive case is
 * diff-back verified, so it can never ship a wrong form.  Returns a
 * freshly-owned tree (caller frees). */



Expr* rt_recursive_tower_case(Expr* f, Expr* x) {
    /* Split any evaluator-merged exponential monomial back into an independent
     * tower basis before building the tower (see rt_expand_exp_sums).  Circular/
     * hyperbolic trig of a tangent argument is rationalised to the tower variable
     * during substitution (rt_subst_kernels), where the fresh symbol prevents the
     * evaluator from canonicalising the Tan-rational form back to Csc/Sec.
     *
     * rt_powers_to_exp first re-exposes any transcendental general power b^e (e.g.
     * x^(I x) = E^(I x Log x), hidden by the E^(c Log b) -> b^c collapse) as a raw
     * base-e kernel, so a mixed log/exp tower over such a power is recognised. */
    Expr* fp = rt_powers_to_exp(f, x);
    Expr* fx = rt_expand_exp_sums(fp);
    expr_free(fp);
    RtTower T;
    /* min_n = 2: single-kernel integrands are handled by the dedicated flat-tower
     * cases above; the general recursion only takes depth->=2 towers. */
    if (!rt_tower_build_min(fx, x, &T, 2)) { rt_tower_free(&T); expr_free(fx); return NULL; }

    /* Alias the kernels to tower variables structurally (NOT via an evaluated
     * ReplaceAll, which would re-merge a split exponential product before
     * substitution — see rt_subst_kernels), then normalise with Together. */
    Expr* F = rt_eval1("Together", rt_subst_kernels(fx, &T));
    Expr* result = NULL;
    bool tower_verified = false;
    if (F && rt_find_exp_of_x(F, x) == NULL && rt_find_log_of_x(F, x) == NULL) {
        Expr* num = rt_eval1("Numerator", expr_copy(F));
        Expr* den = rt_eval1("Denominator", expr_copy(F));
        Expr** vv = malloc((T.n + 1) * sizeof(Expr*));
        vv[0] = expr_copy(x);
        for (size_t i = 0; i < T.n; i++) vv[i + 1] = expr_copy(T.t[i]);
        Expr* vlist = expr_new_function(expr_new_symbol("List"), vv, T.n + 1);
        free(vv);
        Expr* pqn = num ? rt_eval2("PolynomialQ", expr_copy(num), expr_copy(vlist)) : NULL;
        Expr* pqd = den ? rt_eval2("PolynomialQ", expr_copy(den), expr_copy(vlist)) : NULL;
        bool gate = num && den && rt_is_true(pqn) && rt_is_true(pqd);
        if (pqn) expr_free(pqn);
        if (pqd) expr_free(pqd);
        expr_free(vlist);
        if (num) expr_free(num);
        if (den) expr_free(den);
        if (gate) {
            Expr* rem_tower = NULL;   /* partial log part: unintegrated r_n */
            Expr* Q = rt_field_integrate(F, &T, (long)T.n - 1, x, &rem_tower);
            if (Q) {
                /* Reduce Q to lowest terms in the TOWER VARIABLES (each kernel
                 * t_i an independent symbol) BEFORE back-substituting.  A common
                 * factor like (1 + t_1) in Q = t_3 (1+t_1) / ((1+t_1) t_2)
                 * cancels cleanly here, but not afterwards when t_1 = E^x
                 * reappears inside t_2 = Log[1 + E^x] and blocks the gcd —
                 * otherwise a correct answer (Log[x]/Log[1 + E^x]) is left in an
                 * unrecognisable form.  Together first: Q is a SUM of Laurent/log
                 * terms, which Cancel alone does not combine before reducing. */
                Q = rt_eval1("Cancel", rt_eval1("Together", Q));
                /* EXACT verification in the tower variables: with a partial log
                 * part (Bronstein Thm 5.6.1) the identity is D_tower[Q] + rem =
                 * F, i.e. the elementary Q integrates all of F except the
                 * unintegrated non-constant-residue remainder rem_tower.  When
                 * rem_tower is NULL this is the usual exact D_tower[Q] == F. */
                Expr* DQ = rt_tower_deriv(Q, &T, x);
                Expr* diff = DQ ? expr_new_function(expr_new_symbol("Plus"),
                    (Expr*[]){ expr_copy(DQ), expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_integer(-1), expr_copy(F) }, 2) }, 2) : NULL;
                if (diff && rem_tower)
                    diff = expr_new_function(expr_new_symbol("Plus"),
                        (Expr*[]){ diff, expr_copy(rem_tower) }, 2);
                /* Expand before Together: a summed derivation coefficient (e.g.
                 * Dcoef = (k - t_tan x)/x for Log[x^k Cos x]) leaves an
                 * undistributed product like x^5 (k/x - t_tan) in D_tower[Q], which
                 * Together alone does not multiply out — so a genuinely-zero
                 * residual would read as nonzero and a correct antiderivative would
                 * be spuriously rejected.  Expand distributes the products; Together
                 * then clears any lower-field fractions so rt_is_zero is exact. */
                Expr* chk = diff ? rt_eval1("Together", rt_eval1("Expand", diff)) : NULL;
                bool tower_ok = chk && rt_is_zero(chk);
                if (DQ) expr_free(DQ);
                if (chk) expr_free(chk);
                if (tower_ok) {
                    Expr** back = malloc(T.n * sizeof(Expr*));
                    for (size_t i = 0; i < T.n; i++)
                        back[i] = expr_new_function(expr_new_symbol("Rule"),
                            (Expr*[]){ expr_copy(T.t[i]), expr_copy(T.kernel[i]) }, 2);
                    Expr* bl = expr_new_function(expr_new_symbol("List"), back, T.n);
                    free(back);
                    result = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ Q, bl }, 2));                   /* adopts Q, bl */
                    if (result) result = rt_eval1("Cancel", result);
                    /* Surface the unintegrated non-constant-residue part as a
                     * genuine Integrate[rem, x] in the original kernels (a
                     * separate back-rule list — ReplaceAll consumes its list).
                     * rem's residues are all non-constant, so re-integrating it
                     * declines (r_s empty ⇒ NULL): a fixed point, no loop. */
                    if (result && rem_tower) {
                        Expr** back2 = malloc(T.n * sizeof(Expr*));
                        for (size_t i = 0; i < T.n; i++)
                            back2[i] = expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_copy(T.t[i]), expr_copy(T.kernel[i]) }, 2);
                        Expr* bl2 = expr_new_function(expr_new_symbol("List"), back2, T.n);
                        free(back2);
                        Expr* rem_k = rt_eval1("Cancel", rt_eval_own(
                            expr_new_function(expr_new_symbol("ReplaceAll"),
                                (Expr*[]){ expr_copy(rem_tower), bl2 }, 2)));
                        Expr* rem_int = expr_new_function(expr_new_symbol("Integrate"),
                            (Expr*[]){ rem_k, expr_copy(x) }, 2);
                        result = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                            (Expr*[]){ result, rem_int }, 2));
                    }
                    tower_verified = true;
                } else {
                    expr_free(Q);
                }
            }
            if (rem_tower) expr_free(rem_tower);
        }
    }
    if (F) expr_free(F);

    /* The tower-variable check above is exact; fall back to the (bounded)
     * Simplify diff-back only if it was not reached (defensive). */
    if (result && !tower_verified && !rt_verify_antideriv(result, f, x)) {
        expr_free(result); result = NULL;
    }

    rt_tower_free(&T);
    expr_free(fx);
    return result;
}

/* Solve the Risch differential equation D[y] + f y = g over the transcendental
 * TOWER inferred from the Log/Exp kernels of f and g (Gap 1, rde_tower).  Builds
 * the tower over C(x), aliases the kernels to tower variables, solves in the full
 * field K_{n-1}, and back-substitutes.  Returns y (owned) or NULL (no rational
 * solution in the tower field / out of the 1a scope).  Exposed to the REPL via
 * Risch`RischDE for direct unit testing at depth. */
Expr* rde_solve_tower(Expr* f, Expr* g, Expr* x) {
    Expr* fg = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ expr_copy(f), expr_copy(g) }, 2);
    Expr* fgx = rt_expand_exp_sums(fg);
    expr_free(fg);
    RtTower T;
    if (!rt_tower_build_min(fgx, x, &T, 1)) { rt_tower_free(&T); expr_free(fgx); return NULL; }
    expr_free(fgx);
    Expr* fx = rt_expand_exp_sums(f);
    Expr* gx = rt_expand_exp_sums(g);
    Expr* fs = rt_eval1("Together", rt_subst_kernels(fx, &T));
    Expr* gs = rt_eval1("Together", rt_subst_kernels(gx, &T));
    expr_free(fx); expr_free(gx);
    Expr* q = NULL;
    if (fs && gs
        && rt_find_exp_of_x(fs, x) == NULL && rt_find_log_of_x(fs, x) == NULL
        && rt_find_exp_of_x(gs, x) == NULL && rt_find_log_of_x(gs, x) == NULL) {
        RdeCtx C = { &T, (long)T.n - 1, x, T.t[T.n - 1] };
        q = rde_tower(fs, gs, &C);
    }
    if (fs) expr_free(fs);
    if (gs) expr_free(gs);
    Expr* y = NULL;
    if (q) {
        Expr** back = malloc(T.n * sizeof(Expr*));
        for (size_t i = 0; i < T.n; i++)
            back[i] = expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ expr_copy(T.t[i]), expr_copy(T.kernel[i]) }, 2);
        Expr* bl = expr_new_function(expr_new_symbol("List"), back, T.n);
        free(back);
        y = rt_eval1("Cancel", rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
            (Expr*[]){ q, bl }, 2)));
    }
    rt_tower_free(&T);
    return y;
}

/* ---- real reconstruction of I-laden trigonometric antiderivatives --------- */
/* The complex-exponential route integrates over the kernel E^(I x); its log part
 * comes back as complex logs / ArcTan of E^(±I x), a correct but I-laden closed
 * form (e.g. ∫Csc x -> Log[E^(I x)-1] - Log[E^(I x)+1], which equals the real
 * Log[Tan[x/2]] up to a constant).  cx_reim decomposes such an expression of the
 * REAL variable x into real + I·imaginary parts (treating x, Pi, E and any
 * trig/Sqrt of a real argument as real); rt_realify keeps the real part — a valid
 * antiderivative since the discarded imaginary part is a constant for a real
 * integrand — and returns it only behind an exact diff-back check.  This is the
 * LogToReal real reconstruction for the exponential tower, not a rewrite of the
 * integrand. */
