/* risch_singleext.c — flat single-extension transcendental Risch cases.
 *
 * One-kernel integrators over C(x)(theta) for a single monomial theta = Log[u]
 * or E^u: the exponential Laurent case, the fractional Hermite / LRT cases, and
 * the non-commensurate exponential-sum case, plus the Bronstein RdeBoundDegree
 * helper.  Dispatched before the recursive tower engine.  See risch_singleext.h.
 */

#include "risch_singleext.h"
#include "integrate_risch_transcendental.h"   /* rt_rde_var_bound decl */
#include "integrate_risch_rde.h"              /* rt_solve_rde */
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

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/* Integrate g in K = C(x) allowing genuine logarithms in the result (the
 * i = 0 term of the exponential case is an ordinary base-field integral). */
static Expr* rt_integrate_in_K_with_logs(Expr* g, Expr* x) {
    Expr* r = rt_eval_call("Integrate`BronsteinRational",
        (Expr*[]){ expr_copy(g), expr_copy(x) }, 2);
    if (!r) return NULL;
    if (rt_head_is(r, "Integrate`BronsteinRational")) { expr_free(r); return NULL; }
    return r;
}

/* Exact degree bound for the solution q of a Risch DE  D[q] + f q = p, in a single
 * monomial variable v (Bronstein RdeBoundDegree).  Returns an UPPER bound on deg_v(q)
 * with NO arbitrary cap — a function of the equation's degrees alone (plus the exact
 * cancellation/resonance sub-case, below).  `deriv_lowers` is true when D lowers deg_v
 * by one (v = x under d/dx, or a LOGARITHMIC monomial), false when D preserves deg_v (an
 * EXPONENTIAL monomial, whose self-derivative D[t^k] = k w' t^k keeps the degree).
 * dpv = deg_v(p), dfv = deg_v(f) (each as deg(numerator) - deg(denominator)).
 *
 * Leading-degree balance.  Match the top v-degree of D[q] + f q = p.  For a
 * deriv-lowering v, deg_v(D[q]) = deg_v(q) - 1 and deg_v(f q) = dfv + deg_v(q); when
 * dfv >= 0 the f q term strictly dominates (distinct degrees, no cancellation) so
 * deg_v(q) = dpv - dfv exactly, and when dfv < 0 the balance is between D[q]
 * (integration raises degree by one: dpv + 1) and a pole in f (dpv - dfv), so the bound
 * is their max.  For a deriv-preserving (exponential) v, both terms sit at
 * deg_v(q) + max(0, dfv), giving deg_v(q) = dpv - dfv.
 *
 * Cancellation / resonance (Bronstein's recursive degree reduction).  The leading-degree
 * balance is exact EXCEPT where the two leading coefficients cancel, so the top term of
 * D[q] + f q vanishes and deg_v(q) can exceed the naive value.  This happens in exactly
 * two configurations, each keyed by an integer `m_res` (the Bronstein resonance integer,
 * or -1 when none):
 *   - deriv-PRESERVING v (exponential), dfv == 0: D[c v^n] + f(c v^n) has leading
 *     coefficient (n D[v]/v + f) c, which vanishes when n = -f/(D[v]/v) is a
 *     nonnegative integer m_res; then deg_v(q) can reach m_res, so bound = max(naive,
 *     m_res).  m_res is supplied by the caller (rt_resonance_int), which alone has the
 *     coefficients.
 *   - deriv-LOWERING v (base x / logarithmic), dfv == -1: the primitive
 *     leading-coefficient cancellation, whose homogeneous-solution degree can likewise
 *     reach m_res, so bound = max(naive, m_res).
 * The widening is MONOTONE (only ever raises the bound), so a spurious m_res can never
 * ship a wrong result (the solution stays SolveAlways-certified and the enclosing tower
 * case diff-back verifies); a missed one merely declines.  Pass m_res = -1 to disable.
 *
 * Reachability note.  In the current tower architecture BOTH cancellation configurations
 * are pre-empted, so m_res is -1 on every reachable call: (a) an exponential resonance
 * n = -(i w_L')/w_j' being an integer means the top and lower exp exponents are
 * commensurate, and the commensurate-exponent reduction in rt_tower_build collapses such
 * kernels to one primitive before any RDE solve; (b) dfv == -1 requires a simple pole in
 * f = i Dcoef, which a rational tower element's derivative cannot have (it would
 * integrate to a Log), and the only kernel that could — a log top with Dcoef = u'/u —
 * never routes through the field RDE (it uses the primitive-polynomial recursion).  The
 * detection is nonetheless computed live and folded in here so the degree bound is exact
 * per Bronstein should a future kernel type expose either configuration. */
long rt_rde_var_bound(long dpv, long dfv, bool deriv_lowers, long m_res) {
    long bq;
    if (deriv_lowers) {
        if (dfv >= 0) {
            bq = dpv - dfv;                         /* f q dominates: exact */
        } else {                                   /* integration rise vs pole in f */
            long a = dpv + 1, b = dpv - dfv;
            bq = (a > b) ? a : b;
            /* primitive leading-coefficient cancellation (dfv == -1) */
            if (dfv == -1 && m_res >= 0 && m_res > bq) bq = m_res;
        }
    } else {                                       /* exponential: derivation preserves deg */
        bq = dpv - dfv;
        /* exponential integer resonance (dfv == 0) */
        if (dfv == 0 && m_res >= 0 && m_res > bq) bq = m_res;
    }
    return (bq < 0) ? 0 : bq;
}





/* Exponential (hyperexponential) Laurent-polynomial case.  For a single
 * kernel theta = E^u (u a polynomial in x, theta' = u' theta), an integrand
 * that is a Laurent polynomial  f = sum_i p_i theta^i  (i possibly negative)
 * integrates to Q = sum_i q_i theta^i: the powers DECOUPLE, the i = 0 term is
 * an ordinary integral q_0 = INT p_0 in K, and each i != 0 term solves the
 * Risch differential equation q_i' + i u' q_i = p_i.  All E^(k u) present must
 * share the one primitive exponent u (k integer); a genuine theta-denominator
 * beyond a pure power of theta (a real fractional part) is declined here and
 * left to the fractional case.  Closes e.g. Cosh/Sinh and mixed +/- powers. */
/* Kernelize a single-primitive exponential extension: find a primitive
 * exponent u so every E^w kernel of `f` (w depending on x) has w = k u with
 * k a nonzero integer, substitute every E^(k u) -> rmT^k, and return the
 * resulting expression in rmT (sets *u_out to an owned copy of u).  Returns
 * NULL when there is no single-primitive exponential structure. */
Expr* rt_exp_kernelize(Expr* f, Expr* x, Expr** u_out) {
    *u_out = NULL;
    if (!rt_find_exp_of_x(f, x)) return NULL;
    Expr** ws = NULL; size_t nw = 0, cap = 0;
    rt_collect_exp_exponents(f, x, &ws, &nw, &cap);
    if (nw == 0) { free(ws); return NULL; }

    /* Synthesize the single primitive exponent u so every E^(w_j) = (E^u)^k_j
     * (k_j integer): u = ws[0]/lcm(ratio denominators).  Handles integer ratios
     * (E^(2u) = (E^u)^2) and genuine rational ratios (E^(x/2), E^(x/3) → E^(x/6))
     * alike; declines (NULL) when the exponents are not all a rational multiple
     * of one another — genuinely independent kernels (a sum / tower). */
    long* kof = malloc(nw * sizeof(long));
    Expr* u = rt_class_primitive(ws, nw, kof);   /* owned primitive, or NULL */
    /* The primitive exponent must be rational in x alone; a nested exponent
     * (e.g. u = E^x in E^(E^x)) is a two-extension tower left to rt_exp_tower_case
     * (else the single-kernel RDE would carry the inner kernel as a free param). */
    if (u && !rt_kernel_simple(u, x)) { expr_free(u); u = NULL; }
    if (!u) { for (size_t i = 0; i < nw; i++) expr_free(ws[i]); free(ws); free(kof); return NULL; }

    Expr** rules = malloc(2 * nw * sizeof(Expr*));
    for (size_t j = 0; j < nw; j++) {
        Expr* tk = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_new_symbol("rmT"), expr_new_integer(kof[j]) }, 2);
        rules[2 * j] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Exp"),
                (Expr*[]){ expr_copy(ws[j]) }, 1), expr_copy(tk) }, 2);
        rules[2 * j + 1] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_symbol("E"), expr_copy(ws[j]) }, 2), tk }, 2);
    }
    Expr* rl = expr_new_function(expr_new_symbol("List"), rules, 2 * nw);
    free(rules);
    Expr* uexp = u;                       /* adopt the synthesized primitive */
    for (size_t i = 0; i < nw; i++) expr_free(ws[i]);
    free(ws); free(kof);

    Expr* F = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ expr_copy(f), rl }, 2));
    if (!F) { expr_free(uexp); return NULL; }
    *u_out = uexp;
    return F;
}

Expr* rt_exp_poly_case(Expr* f, Expr* x) {
    Expr* uexp = NULL;
    Expr* F = rt_exp_kernelize(f, x, &uexp);
    if (!F) return NULL;

    Expr* tsym = expr_new_symbol("rmT");
    /* F must be a Laurent polynomial in t: num/den with den a pure power of t.
     * Split into num and the offset M (den = c t^M). */
    Expr* G = rt_eval1("Together", expr_copy(F));
    Expr* num = G ? rt_eval1("Numerator", expr_copy(G)) : NULL;
    Expr* den = G ? rt_eval1("Denominator", expr_copy(G)) : NULL;
    Expr* result = NULL;
    long M = 0;
    bool ok = num && den && rt_is_poly(num, tsym) && rt_is_poly(den, tsym)
        && !rt_free_of_x(F, tsym)
        && rt_find_exp_of_x(F, x) == NULL && rt_find_log_of_x(F, x) == NULL;
    if (ok) {
        M = rt_degree(den, tsym);
        /* den must be a monomial c t^M: den / t^M free of t. */
        Expr* dq = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(tsym), expr_new_integer(-M) }, 2) }, 2));
        if (!dq || !rt_free_of_x(dq, tsym)) ok = false;
        if (dq) expr_free(dq);
        /* Fold the constant den/t^M into num so p_i = Coefficient[num/c, ...]. */
        if (ok) {
            Expr* nnew = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(num), expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(den), expr_new_integer(-1) }, 2),
                  expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(tsym), expr_new_integer(M) }, 2) }, 3));
            if (nnew) { expr_free(num); num = nnew; }
            else ok = false;
        }
    }
    long dnum = ok ? rt_degree(num, tsym) : -1;
    if (ok && dnum >= 0) {
        Expr** q = calloc((size_t)(dnum + 1), sizeof(Expr*));
        long* qi_pow = malloc((size_t)(dnum + 1) * sizeof(long));
        bool fail = false;
        for (long j = 0; j <= dnum && !fail; j++) {
            long i = j - M;                       /* Laurent power */
            Expr* pi = rt_coeff(num, tsym, j);
            if (rt_is_zero(pi)) { expr_free(pi); q[j] = NULL; qi_pow[j] = i; continue; }
            Expr* qi = (i == 0) ? rt_integrate_in_K_with_logs(pi, x)
                                : rt_solve_rde(pi, i, uexp, x);
            expr_free(pi);
            qi_pow[j] = i;
            if (!qi) { fail = true; break; }
            q[j] = qi;
        }
        if (!fail) {
            Expr** terms = malloc((size_t)(dnum + 1) * sizeof(Expr*));
            size_t nt = 0;
            for (long j = 0; j <= dnum; j++) {
                if (!q[j]) continue;
                long i = qi_pow[j];
                Expr* term;
                if (i == 0) {
                    term = expr_copy(q[j]);
                } else {
                    Expr* iu = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_integer(i), expr_copy(uexp) }, 2);
                    Expr* eiu = expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_new_symbol("E"), iu }, 2);
                    term = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_copy(q[j]), eiu }, 2);
                }
                terms[nt++] = term;
            }
            result = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                terms, nt));
            free(terms);
        }
        for (long j = 0; j <= dnum; j++) if (q[j]) expr_free(q[j]);
        free(q); free(qi_pow);
    }

    if (G) expr_free(G);
    if (num) expr_free(num);
    if (den) expr_free(den);
    expr_free(F);
    expr_free(tsym);
    expr_free(uexp);
    return result;
}

/* Fractional (Rothstein-Trager) log-part for a single monomial extension
 * theta (theta = Log[u], D theta = u'/u; or theta = E^u, D theta = u' theta).
 * With theta -> t, an integrand whose t-denominator d = prod g_i is
 * squarefree integrates (when elementary and free of a polynomial/Hermite
 * remainder) to  sum_i c_i Log(g_i)  with each c_i a CONSTANT.  Those
 * constants are found by solving the exact polynomial identity
 *     num  =  sum_i c_i D(g_i) (d / g_i)          (D = d/dx + (D t) d/dt)
 * for all t and x with SolveAlways[..., {t, x}]: the solver returns the
 * constant residues, or no solution (which declines) when the integrand has
 * a repeated pole, an irreducible factor with non-constant residue, or a
 * polynomial part.  Correct by construction — a solution certifies
 * D(sum c_i Log g_i) = num/d exactly. */
static Expr* rt_frac_try(Expr* f, Expr* x, Expr* u, bool is_log) {
    Expr* tsym = expr_new_symbol("rmT");
    Expr* rules; Expr* Dt; Expr* kernel_back;
    if (is_log) {
        Expr* logu = expr_new_function(expr_new_symbol("Log"),
            (Expr*[]){ expr_copy(u) }, 1);
        rules = expr_new_function(expr_new_symbol("List"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ logu, expr_new_symbol("rmT") }, 2) }, 1);
        Expr* du = rt_eval2("D", expr_copy(u), expr_copy(x));
        Expr* invu = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(u), expr_new_integer(-1) }, 2);
        Dt = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ du, invu }, 2));
        kernel_back = expr_new_function(expr_new_symbol("Log"),
            (Expr*[]){ expr_copy(u) }, 1);
    } else {
        Expr* eu = expr_new_function(expr_new_symbol("Exp"),
            (Expr*[]){ expr_copy(u) }, 1);
        Expr* peu = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_new_symbol("E"), expr_copy(u) }, 2);
        rules = expr_new_function(expr_new_symbol("List"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ eu, expr_new_symbol("rmT") }, 2),
              expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ peu, expr_new_symbol("rmT") }, 2) }, 2);
        Expr* up = rt_eval2("D", expr_copy(u), expr_copy(x));
        Dt = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ up, expr_new_symbol("rmT") }, 2);
        kernel_back = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_new_symbol("E"), expr_copy(u) }, 2);
    }
    if (!Dt) { expr_free(tsym); expr_free(rules); expr_free(kernel_back); return NULL; }

    Expr* F0 = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ expr_copy(f), rules }, 2));
    Expr* F = F0 ? rt_eval1("Together", F0) : NULL;
    if (!F) { expr_free(tsym); expr_free(Dt); expr_free(kernel_back); return NULL; }

    Expr* num = rt_eval1("Numerator", expr_copy(F));
    Expr* den = rt_eval1("Denominator", expr_copy(F));
    Expr* result = NULL;

    bool ok = num && den && rt_is_poly(num, tsym) && rt_is_poly(den, tsym)
        && !rt_free_of_x(den, tsym)
        && rt_find_exp_of_x(F, x) == NULL && rt_find_log_of_x(F, x) == NULL;

    if (ok) {
        Expr* factored = rt_eval1("Factor", expr_copy(den));
        Expr* g[16]; size_t ng = 0; bool bad = false;
        if (factored) {
            Expr** fa; size_t nf; Expr* single[1];
            if (factored->type == EXPR_FUNCTION
                && factored->data.function.head->type == EXPR_SYMBOL
                && factored->data.function.head->data.symbol.name == intern_symbol("Times")) {
                fa = factored->data.function.args;
                nf = factored->data.function.arg_count;
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
                if (rt_free_of_x(base, tsym)) continue;   /* FreeQ[base, t] */
                if (e != 1 || ng >= 16) { bad = true; break; }
                g[ng++] = expr_copy(base);
            }
        }
        if (!bad && ng >= 1) {
            /* residual = num - sum_i c_i D(g_i) (den/g_i). */
            Expr** negterms = malloc((ng + 1) * sizeof(Expr*));
            negterms[0] = expr_copy(num);
            for (size_t i = 0; i < ng; i++) {
                char nm[24]; snprintf(nm, sizeof(nm), "rmK%zu", i);
                Expr* dgx = rt_eval2("D", expr_copy(g[i]), expr_copy(x));
                Expr* dgt = rt_eval2("D", expr_copy(g[i]), expr_copy(tsym));
                Expr* dgi = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                    (Expr*[]){ dgx, expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_copy(Dt), dgt }, 2) }, 2));
                Expr* cof = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(g[i]), expr_new_integer(-1) }, 2) }, 2));
                negterms[i + 1] = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), expr_new_symbol(nm),
                              dgi ? dgi : expr_new_integer(0),
                              cof ? cof : expr_new_integer(0) }, 4);
            }
            Expr* residual = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                negterms, ng + 1));
            free(negterms);

            Expr* varlist = expr_new_function(expr_new_symbol("List"),
                (Expr*[]){ expr_copy(tsym), expr_copy(x) }, 2);
            Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
                (Expr*[]){ residual, expr_new_integer(0) }, 2);
            Expr* sol = rt_eval2("SolveAlways", eqn, varlist);
            if (sol && sol->type == EXPR_FUNCTION
                && sol->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.head->data.symbol.name == intern_symbol("List")
                && sol->data.function.arg_count >= 1
                && sol->data.function.args[0]->type == EXPR_FUNCTION
                && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.args[0]->data.function.head->data.symbol.name
                     == intern_symbol("List")) {
                /* Rothstein-Trager residues MUST be constants.  SolveAlways can
                 * return an x-dependent pseudo-solution for a residual with a
                 * constant-term structure (e.g. k -> x + c for the residual
                 * 1 - k/(x + c) of 1/Log[x + c]); accepting it would certify a
                 * WRONG polynomial-coefficient logarithm ((x+c) Log[Log[x+c]]).
                 * Require every residue free of x and t, else decline. */
                Expr* srules = sol->data.function.args[0];
                bool const_res = true;
                for (size_t si = 0;
                     si < srules->data.function.arg_count && const_res; si++) {
                    Expr* rule = srules->data.function.args[si];
                    if (rule->type == EXPR_FUNCTION
                        && rule->data.function.arg_count == 2
                        && (!rt_free_of_x(rule->data.function.args[1], x)
                            || !rt_free_of_x(rule->data.function.args[1], tsym)))
                        const_res = false;
                }
                if (const_res) {
                Expr** rterms = malloc(ng * sizeof(Expr*));
                for (size_t i = 0; i < ng; i++) {
                    char nm[24]; snprintf(nm, sizeof(nm), "rmK%zu", i);
                    Expr* gib = rt_eval_own(expr_new_function(
                        expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ expr_copy(g[i]),
                          expr_new_function(expr_new_symbol("List"),
                            (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_copy(tsym), expr_copy(kernel_back) },
                                2) }, 1) }, 2));
                    Expr* logg = expr_new_function(expr_new_symbol("Log"),
                        (Expr*[]){ gib }, 1);
                    rterms[i] = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_symbol(nm), logg }, 2);
                }
                Expr* R = expr_new_function(expr_new_symbol("Plus"), rterms, ng);
                free(rterms);
                R = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                    (Expr*[]){ R, expr_copy(sol->data.function.args[0]) }, 2));
                if (R) {
                    Expr** zero = malloc(ng * sizeof(Expr*));
                    for (size_t i = 0; i < ng; i++) {
                        char nm[24]; snprintf(nm, sizeof(nm), "rmK%zu", i);
                        zero[i] = expr_new_function(expr_new_symbol("Rule"),
                            (Expr*[]){ expr_new_symbol(nm), expr_new_integer(0) }, 2);
                    }
                    Expr* zl = expr_new_function(expr_new_symbol("List"), zero, ng);
                    free(zero);
                    result = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ R, zl }, 2));
                }
                }
            }
            if (sol) expr_free(sol);
        }
        for (size_t i = 0; i < ng; i++) expr_free(g[i]);
        if (factored) expr_free(factored);
    }
    if (num) expr_free(num);
    if (den) expr_free(den);
    expr_free(F);
    expr_free(tsym);
    expr_free(Dt);
    expr_free(kernel_back);
    return result;
}

/* Pure resultant Lazard-Rioboo-Trager log part for a single monomial
 * extension theta (theta = Log[u] or E^u), for the residue class the
 * SolveAlways single-constant-per-factor ansatz of rt_frac_try cannot
 * express: an irreducible-over-Q denominator factor in theta whose
 * Rothstein-Trager residues are ALGEBRAIC (non-rational) constants — e.g.
 * the +-I/2 residues of 1/(x (Log[x]^2+1)) that split the irreducible
 * theta^2+1 into a conjugate log pair (= ArcTan[Log[x]]).
 *
 * With theta -> t (rmT) and D t = u'/u (log) or u' t (exp), a PROPER,
 * SQUAREFREE-denominator fraction a(t)/d(t) has log part
 *   sum over roots c of Res_t(a - z D(d), d):  c Log(gcd_t(a - c D(d), d)),
 * where D(d) = d/dx(d) + (D t) d/dt(d) is the monomial derivation.  The
 * exact resultant + Rioboo LogToReal real-form collapse are delegated to
 * Integrate`TranscendentalLogPart (intrat.c), which reuses the tested
 * rational-LRT subresultant machinery; this routine just builds the
 * monomial derivation, hands off, substitutes t -> kernel, and DIFF-BACK
 * VERIFIES (the reuse crosses content/denominator boundaries, so — like the
 * tower cases — a mis-reduction must decline, never ship a wrong form).
 * Runs only after rt_frac_try declines. */
static Expr* rt_frac_lrt(Expr* f, Expr* x, Expr* u, bool is_log) {
    Expr* tsym = expr_new_symbol("rmT");
    Expr* Dt = NULL; Expr* kernel_back = NULL; Expr* F = NULL;
    if (is_log) {
        /* Log[u] -> t; Log[u]^k is Power[Log[u], k] so the inner rewrite
         * lifts to t^k automatically. */
        Expr* rules = expr_new_function(expr_new_symbol("List"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ expr_new_function(expr_new_symbol("Log"),
                    (Expr*[]){ expr_copy(u) }, 1), expr_new_symbol("rmT") }, 2) }, 1);
        Expr* du = rt_eval2("D", expr_copy(u), expr_copy(x));
        Expr* invu = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(u), expr_new_integer(-1) }, 2);
        Dt = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ du, invu }, 2));                 /* eta = u'/u */
        kernel_back = expr_new_function(expr_new_symbol("Log"),
            (Expr*[]){ expr_copy(u) }, 1);
        Expr* F0 = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
            (Expr*[]){ expr_copy(f), rules }, 2));
        F = F0 ? rt_eval1("Together", F0) : NULL;
    } else {
        /* Kernelize so multiplicatively commensurate exponents collapse onto
         * one primitive t = E^uexp (E^(2x) -> t^2), which a plain E^x -> t
         * substitution cannot do. */
        Expr* uexp = NULL;
        Expr* Fk = rt_exp_kernelize(f, x, &uexp);
        if (Fk && uexp) {
            Expr* up = rt_eval2("D", expr_copy(uexp), expr_copy(x));
            Dt = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ up, expr_new_symbol("rmT") }, 2);   /* u' t */
            kernel_back = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_symbol("E"), expr_copy(uexp) }, 2);
            F = rt_eval1("Together", Fk);
        } else if (Fk) { expr_free(Fk); }
        if (uexp) expr_free(uexp);
    }
    (void)u;   /* exp branch derives its own primitive via rt_exp_kernelize */
    if (!Dt || !F) {
        expr_free(tsym); if (Dt) expr_free(Dt);
        if (kernel_back) { expr_free(kernel_back); } if (F) { expr_free(F); }
        return NULL;
    }

    Expr* num = rt_eval1("Numerator", expr_copy(F));
    Expr* den = rt_eval1("Denominator", expr_copy(F));
    Expr* result = NULL;

    /* Gates: proper squarefree rational function of the single kernel t,
     * with the kernel fully substituted (no residual nested exp/log of x —
     * that would let the resultant treat it as a free parameter). */
    bool ok = num && den && rt_is_poly(num, tsym) && rt_is_poly(den, tsym)
        && !rt_free_of_x(den, tsym)
        && rt_find_exp_of_x(F, x) == NULL && rt_find_log_of_x(F, x) == NULL
        && rt_degree(num, tsym) < rt_degree(den, tsym);

    if (ok) {
        /* Squarefree denominator: gcd(d, dd/dt) constant. */
        Expr* ddt0 = rt_eval2("D", expr_copy(den), expr_copy(tsym));
        Expr* g = ddt0 ? rt_eval_call("PolynomialGCD",
            (Expr*[]){ expr_copy(den), ddt0 }, 2) : NULL;
        bool squarefree = g && rt_degree(g, tsym) <= 0;
        if (g) expr_free(g);

        if (squarefree) {
            /* Monomial derivation D(d) = d/dx(d) + (D t) d/dt(d). */
            Expr* ddx = rt_eval2("D", expr_copy(den), expr_copy(x));
            Expr* ddt = rt_eval2("D", expr_copy(den), expr_copy(tsym));
            Expr* Dd = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ ddx, expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_copy(Dt), ddt }, 2) }, 2));
            if (Dd) {
                Expr* logpart = rt_eval_call("Integrate`TranscendentalLogPart",
                    (Expr*[]){ expr_copy(num), expr_copy(den), expr_copy(tsym),
                              expr_new_symbol("rmZ"), Dd, expr_copy(x) }, 6);
                /* Unevaluated (declined) iff the head is unchanged. */
                bool declined = !logpart
                    || (logpart->type == EXPR_FUNCTION
                        && logpart->data.function.head->type == EXPR_SYMBOL
                        && logpart->data.function.head->data.symbol.name
                             == intern_symbol("Integrate`TranscendentalLogPart"));
                /* Partial log part (Bronstein Thm 5.6.1): the residue resultant
                 * mixed constant and non-constant residues.  Surface the
                 * elementary constant-residue logs plus an unevaluated
                 * Integrate[remainder, x] for the non-constant part.  The FTC
                 * rule D[Integrate[f,x],x] = f makes rt_verify_antideriv close
                 * the diff-back on the assembled partial form. */
                bool partial = logpart && logpart->type == EXPR_FUNCTION
                    && logpart->data.function.head->type == EXPR_SYMBOL
                    && logpart->data.function.head->data.symbol.name
                         == intern_symbol("Integrate`PartialLogPart")
                    && logpart->data.function.arg_count == 2;
                if (partial) {
                    Expr* krule = expr_new_function(expr_new_symbol("List"),
                        (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                            (Expr*[]){ expr_copy(tsym), expr_copy(kernel_back) }, 2) }, 1);
                    Expr* logs_k = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ expr_copy(logpart->data.function.args[0]), expr_copy(krule) }, 2));
                    Expr* rem_k = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ expr_copy(logpart->data.function.args[1]), krule }, 2));
                    Expr* Q = (logs_k && rem_k) ? rt_eval_own(expr_new_function(
                        expr_new_symbol("Plus"),
                        (Expr*[]){ logs_k,
                            expr_new_function(expr_new_symbol("Integrate"),
                                (Expr*[]){ rem_k, expr_copy(x) }, 2) }, 2)) : NULL;
                    if (!Q) { if (logs_k) expr_free(logs_k); if (rem_k) expr_free(rem_k); }
                    if (Q && rt_verify_antideriv(Q, f, x)) result = Q;
                    else if (Q) expr_free(Q);
                    expr_free(logpart);
                } else if (!declined) {
                    /* Substitute t -> kernel(x) and diff-back verify. */
                    Expr* Q = rt_eval_own(expr_new_function(
                        expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ logpart,
                          expr_new_function(expr_new_symbol("List"),
                            (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_copy(tsym),
                                          expr_copy(kernel_back) }, 2) }, 1) }, 2));
                    if (Q && rt_verify_antideriv(Q, f, x)) result = Q;
                    else if (Q) expr_free(Q);
                } else if (logpart) {
                    expr_free(logpart);
                }
            }
        }
    }

    if (num) expr_free(num);
    if (den) expr_free(den);
    expr_free(F);
    expr_free(tsym);
    expr_free(Dt);
    expr_free(kernel_back);
    return result;
}

/* Hermite reduction for a REPEATED transcendental pole, for a logarithmic
 * (theta = Log[u], D theta = u'/u) or exponential (theta = E^u, D theta = u'
 * theta) kernel.  A proper rational function A/D in theta whose denominator D
 * has a repeated factor integrates to  Q = H(theta)/Hden(theta) +
 * sum_j c_j Log(g_j)  where Hden = gcd(D, dD/dtheta) = prod D_m^(m-1) is the
 * repeated part, g_j are the squarefree factors of the radical D/Hden, H is a
 * polynomial in theta of degree < deg(Hden) with bounded-degree polynomial-in-x
 * coefficients, and the c_j are constants.  All unknown constants are solved at
 * once by SolveAlways[Q' - f == 0, {t, x}] — correct by construction (e.g.
 * Integrate[1/(x (1+Log[x])^2), x] = -1/(1+Log[x]),
 * Integrate[E^x/(1+E^x)^2, x] = -1/(1+E^x)).  The exponential kernel is handled
 * only when D is coprime to theta (no theta = 0 pole — that Laurent case is
 * left to the hyperexponential path). */
static Expr* rt_hermite_try(Expr* f, Expr* x, bool is_log) {
    Expr* tsym = expr_new_symbol("rmT");
    Expr* F = NULL; Expr* Dt = NULL; Expr* kback = NULL; Expr* uu = NULL;
    if (is_log) {
        Expr* u = rt_find_log_of_x(f, x);
        if (!u || !rt_kernel_simple(u, x)) { expr_free(tsym); return NULL; }
        uu = expr_copy(u);
        Expr* du = rt_eval2("D", expr_copy(uu), expr_copy(x));
        Expr* invu = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(uu), expr_new_integer(-1) }, 2);
        Dt = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ du, invu }, 2));   /* eta = u'/u */
        kback = expr_new_function(expr_new_symbol("Log"),
            (Expr*[]){ expr_copy(uu) }, 1);
        Expr* rl = expr_new_function(expr_new_symbol("List"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ expr_new_function(expr_new_symbol("Log"),
                    (Expr*[]){ expr_copy(uu) }, 1), expr_new_symbol("rmT") }, 2) }, 1);
        Expr* F0 = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
            (Expr*[]){ expr_copy(f), rl }, 2));
        F = F0 ? rt_eval1("Together", F0) : NULL;
    } else {
        Expr* F0 = rt_exp_kernelize(f, x, &uu);
        if (F0) F = rt_eval1("Together", F0);
        if (uu) {
            Expr* up = rt_eval2("D", expr_copy(uu), expr_copy(x));
            Dt = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ up, expr_new_symbol("rmT") }, 2);   /* u' t */
            kback = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_symbol("E"), expr_copy(uu) }, 2);
        }
    }
    if (!F || !Dt || !uu) {
        expr_free(tsym);
        if (F) { expr_free(F); } if (Dt) { expr_free(Dt); }
        if (kback) { expr_free(kback); } if (uu) { expr_free(uu); }
        return NULL;
    }

    Expr* num = rt_eval1("Numerator", expr_copy(F));
    Expr* den = rt_eval1("Denominator", expr_copy(F));
    Expr* result = NULL;
    bool ok = num && den && rt_is_poly(num, tsym) && rt_is_poly(den, tsym)
        && !rt_free_of_x(den, tsym)
        && rt_find_log_of_x(F, x) == NULL && rt_find_exp_of_x(F, x) == NULL;
    long dnum = ok ? rt_degree(num, tsym) : 0, dden = ok ? rt_degree(den, tsym) : 0;
    ok = ok && dnum < dden;                      /* proper fraction */
    if (ok && !is_log) {
        /* exponential kernel: require D coprime to theta (no theta = 0 pole). */
        Expr* c0 = rt_coeff(den, tsym, 0);
        if (rt_is_zero(c0)) ok = false;
        expr_free(c0);
    }

    Expr* Hden = NULL; Expr* rad = NULL;
    if (ok) {
        Expr* dent = rt_eval2("D", expr_copy(den), expr_copy(tsym));
        Hden = rt_eval_call("PolynomialGCD",
            (Expr*[]){ expr_copy(den), dent }, 2);
        long dH = Hden ? rt_degree(Hden, tsym) : 0;
        if (dH < 1) ok = false;                  /* no repeated pole */
        if (ok) rad = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(Hden), expr_new_integer(-1) }, 2) }, 2));
    }

    if (ok && rad) {
        long dH = rt_degree(Hden, tsym);
        /* distinct t-factors of the squarefree radical (log terms). */
        Expr* factored = rt_eval1("Factor", expr_copy(rad));
        Expr* g[16]; size_t ng = 0; bool bad = (factored == NULL);
        if (factored) {
            Expr** fa; size_t nf; Expr* single[1];
            if (factored->type == EXPR_FUNCTION
                && factored->data.function.head->type == EXPR_SYMBOL
                && factored->data.function.head->data.symbol.name == intern_symbol("Times")) {
                fa = factored->data.function.args; nf = factored->data.function.arg_count;
            } else { single[0] = factored; fa = single; nf = 1; }
            for (size_t i = 0; i < nf && !bad; i++) {
                Expr* term = fa[i]; Expr* base = term;
                if (term->type == EXPR_FUNCTION
                    && term->data.function.head->type == EXPR_SYMBOL
                    && term->data.function.head->data.symbol.name == intern_symbol("Power")
                    && term->data.function.arg_count == 2)
                    base = term->data.function.args[0];
                if (rt_free_of_x(base, tsym)) continue;
                if (ng >= 16) { bad = true; break; }
                g[ng++] = expr_copy(base);
            }
        }
        long dnx = rt_degree(num, x), ddx = rt_degree(den, x);
        long Nx = (dnx > ddx ? dnx : ddx) + 2;   /* derived Hermite x-degree, no cap */
        size_t nh = (size_t)(dH * (Nx + 1));
        if (!bad && (long)(nh + ng) > 0) {
            Expr** hterms = malloc((nh ? nh : 1) * sizeof(Expr*));
            size_t nt = 0;
            for (long p = 0; p < dH; p++)
                for (long k = 0; k <= Nx; k++) {
                    char nm[64]; snprintf(nm, sizeof(nm), "rmHh%ldv%ld", p, k);
                    hterms[nt++] = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_symbol(nm),
                          expr_new_function(expr_new_symbol("Power"),
                            (Expr*[]){ expr_copy(x), expr_new_integer(k) }, 2),
                          expr_new_function(expr_new_symbol("Power"),
                            (Expr*[]){ expr_copy(tsym), expr_new_integer(p) }, 2) }, 3);
                }
            Expr* Hpoly = expr_new_function(expr_new_symbol("Plus"), hterms, nt);
            free(hterms);
            Expr* ratpart = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ Hpoly, expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(Hden), expr_new_integer(-1) }, 2) }, 2);
            Expr** qterms = malloc((1 + ng) * sizeof(Expr*));
            size_t ntq = 0;
            qterms[ntq++] = ratpart;
            for (size_t j = 0; j < ng; j++) {
                char nm[32]; snprintf(nm, sizeof(nm), "rmHc%zu", j);
                qterms[ntq++] = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_symbol(nm),
                      expr_new_function(expr_new_symbol("Log"),
                        (Expr*[]){ expr_copy(g[j]) }, 1) }, 2);
            }
            Expr* Q = expr_new_function(expr_new_symbol("Plus"), qterms, ntq);
            free(qterms);

            Expr* dQx = rt_eval2("D", expr_copy(Q), expr_copy(x));
            Expr* dQt = rt_eval2("D", expr_copy(Q), expr_copy(tsym));
            Expr* Qder = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ dQx, expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_copy(Dt), dQt }, 2) }, 2));
            Expr* diff = expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ Qder, expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(F) }, 2) }, 2);
            Expr* tog = rt_eval1("Together", diff);
            Expr* rnum = tog ? rt_eval1("Numerator", tog) : NULL;
            Expr* sol = NULL;
            if (rnum) {
                Expr* varlist = expr_new_function(expr_new_symbol("List"),
                    (Expr*[]){ expr_copy(tsym), expr_copy(x) }, 2);
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
                    Expr** zero = malloc((nh + ng) * sizeof(Expr*));
                    size_t zi = 0;
                    for (long p = 0; p < dH; p++)
                        for (long k = 0; k <= Nx; k++) {
                            char nm[64]; snprintf(nm, sizeof(nm), "rmHh%ldv%ld", p, k);
                            zero[zi++] = expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_new_symbol(nm), expr_new_integer(0) }, 2);
                        }
                    for (size_t j = 0; j < ng; j++) {
                        char nm[32]; snprintf(nm, sizeof(nm), "rmHc%zu", j);
                        zero[zi++] = expr_new_function(expr_new_symbol("Rule"),
                            (Expr*[]){ expr_new_symbol(nm), expr_new_integer(0) }, 2);
                    }
                    Expr* zl = expr_new_function(expr_new_symbol("List"), zero, nh + ng);
                    free(zero);
                    Qs = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ Qs, zl }, 2));
                    if (Qs) {
                        Expr* back = expr_new_function(expr_new_symbol("List"),
                            (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_copy(tsym), expr_copy(kback) }, 2) }, 1);
                        Expr* rr = rt_eval_own(expr_new_function(
                            expr_new_symbol("ReplaceAll"), (Expr*[]){ Qs, back }, 2));
                        /* Cancel folds away the x-content that Hden carried
                         * (e.g. -x/(x + x Log[x]) -> -1/(1 + Log[x])). */
                        if (rr) result = rt_eval1("Cancel", rr);
                    }
                }
            }
            if (sol) expr_free(sol);
            expr_free(Q);
        }
        for (size_t j = 0; j < ng; j++) expr_free(g[j]);
        if (factored) expr_free(factored);
    }
    if (Hden) expr_free(Hden);
    if (rad) expr_free(rad);
    if (num) expr_free(num);
    if (den) expr_free(den);
    expr_free(Dt);
    expr_free(kback);
    expr_free(uu);
    expr_free(F);
    expr_free(tsym);
    return result;
}

/* Hermite reduction dispatch: logarithmic then exponential kernel. */
Expr* rt_hermite_case(Expr* f, Expr* x) {
    Expr* r = rt_hermite_try(f, x, true);
    if (r) return r;
    r = rt_hermite_try(f, x, false);
    return r;
}

/* Try the fractional log-part for a logarithmic then an exponential kernel. */
Expr* rt_frac_case(Expr* f, Expr* x) {
    Expr* ul = rt_find_log_of_x(f, x);
    if (ul && rt_kernel_simple(ul, x)) {
        Expr* r = rt_frac_try(f, x, ul, true); if (r) return r;
    }
    Expr* ue = rt_find_exp_of_x(f, x);
    if (ue && rt_kernel_simple(ue, x)) {
        Expr* r = rt_frac_try(f, x, ue, false); if (r) return r;
    }
    /* Pure resultant LRT for algebraic (non-rational) residues — runs only
     * when the rational-residue SolveAlways path above declines. */
    if (ul && rt_kernel_simple(ul, x)) {
        Expr* r = rt_frac_lrt(f, x, ul, true); if (r) return r;
    }
    if (ue && rt_kernel_simple(ue, x)) {
        Expr* r = rt_frac_lrt(f, x, ue, false); if (r) return r;
    }
    return NULL;
}

/* Coupled hyperexponential case: a general rational function of a single
 * exponential monomial theta = E^u (u polynomial in x) whose integral mixes a
 * Laurent-polynomial part with a logarithmic part — the case where D(g)/g is
 * itself improper in theta, so the polynomial and log parts do NOT separate
 * (e.g. Tan/Tanh after TrigToExp, or 1/(1+E^x)).  A UNIFIED ansatz
 *     Q = sum_i (sum_k a_{ik} x^k) theta^i + sum_j c_j Log(g_j)
 * (g_j the squarefree factors of the theta-coprime denominator) is solved for
 * all its constant coefficients a_{ik}, c_j at once by requiring Q' = f for all
 * theta and x via SolveAlways[..., {t, x}].  Correct by construction: an exact
 * solution certifies the identity; an imperfect degree bound can only decline.
 * Tried after the pure polynomial and pure-log-part cases (which are cheaper). */
Expr* rt_hyperexp_case(Expr* f, Expr* x) {
    Expr* uexp = NULL;
    Expr* F = rt_exp_kernelize(f, x, &uexp);
    if (!F) return NULL;
    Expr* tsym = expr_new_symbol("rmT");
    Expr* result = NULL;

    Expr* G = rt_eval1("Together", expr_copy(F));
    Expr* num = G ? rt_eval1("Numerator", expr_copy(G)) : NULL;
    Expr* den = G ? rt_eval1("Denominator", expr_copy(G)) : NULL;
    bool ok = num && den && rt_is_poly(num, tsym) && rt_is_poly(den, tsym)
        && !rt_free_of_x(den, tsym)
        && rt_find_exp_of_x(F, x) == NULL && rt_find_log_of_x(F, x) == NULL;

    if (ok) {
        /* a = multiplicity of t in den (first nonzero t-coefficient). */
        long a = 0;
        Expr* cl = rt_eval2("CoefficientList", expr_copy(den), expr_copy(tsym));
        if (cl && cl->type == EXPR_FUNCTION
            && cl->data.function.head->type == EXPR_SYMBOL
            && cl->data.function.head->data.symbol.name == intern_symbol("List")) {
            for (size_t i = 0; i < cl->data.function.arg_count; i++)
                if (!rt_is_zero(cl->data.function.args[i])) { a = (long)i; break; }
        }
        if (cl) expr_free(cl);
        Expr* Dtil = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(tsym), expr_new_integer(-a) }, 2) }, 2));
        /* Split the theta-coprime denominator Dtil into its repeated part
         * Hden = gcd(Dtil, d Dtil/dt) (= prod D_m^(m-1)) and squarefree radical
         * rad = Dtil/Hden.  The distinct t-factors of rad carry constant-residue
         * logarithms; the repeated part is absorbed by a Hermite term
         * H(t)/Hden(t) (H proper in t with x-polynomial coefficients).  When Dtil
         * is squarefree Hden = 1 and the Hermite term is empty, recovering the
         * plain Laurent + log hyperexponential ansatz.  This is the Phase A
         * generalization that couples the Laurent theta = 0 pole with a Hermite
         * repeated pole, closing shapes such as 1/(1 + E^x)^2 and
         * 1/(E^x (1 + E^x)^2) that the pure Hermite / squarefree paths decline. */
        Expr* Hden = NULL; Expr* rad = NULL; long dH = 0;
        if (Dtil) {
            Expr* dDt = rt_eval2("D", expr_copy(Dtil), expr_copy(tsym));
            Hden = rt_eval_call("PolynomialGCD",
                (Expr*[]){ expr_copy(Dtil), dDt }, 2);
            if (Hden) {
                dH = rt_degree(Hden, tsym);
                if (dH < 0) dH = 0;
                rad = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_copy(Dtil), expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(Hden), expr_new_integer(-1) }, 2) }, 2));
            }
        }
        Expr* factored = rad ? rt_eval1("Factor", expr_copy(rad)) : NULL;
        Expr* g[16]; size_t ng = 0; bool bad = (factored == NULL);
        if (factored) {
            Expr** fa; size_t nf; Expr* single[1];
            if (factored->type == EXPR_FUNCTION
                && factored->data.function.head->type == EXPR_SYMBOL
                && factored->data.function.head->data.symbol.name == intern_symbol("Times")) {
                fa = factored->data.function.args; nf = factored->data.function.arg_count;
            } else { single[0] = factored; fa = single; nf = 1; }
            for (size_t i = 0; i < nf && !bad; i++) {
                Expr* term = fa[i]; Expr* base = term;
                if (term->type == EXPR_FUNCTION
                    && term->data.function.head->type == EXPR_SYMBOL
                    && term->data.function.head->data.symbol.name == intern_symbol("Power")
                    && term->data.function.arg_count == 2)
                    base = term->data.function.args[0];   /* rad squarefree: e = 1 */
                if (rt_free_of_x(base, tsym)) continue;    /* FreeQ[base, t] */
                if (ng >= 16) { bad = true; break; }
                g[ng++] = expr_copy(base);
            }
        }

        long dnum = ok ? rt_degree(num, tsym) : 0, dden = ok ? rt_degree(den, tsym) : 0;
        long ihi = dnum - dden; if (ihi < 0) ihi = 0;
        long ilo = -a;
        long degu = rt_degree(uexp, x);
        long dnx = rt_degree(num, x), ddx = rt_degree(den, x);
        /* x-degree bound of the coefficient polynomials — derived from the numerator/
         * denominator degrees and the exponent degree, NO cap.  Generous but a
         * function of the input; a too-large bound only slows SolveAlways (every
         * solution is certified + diff-back verified), never mis-solves. */
        long Nx = (dnx > ddx ? dnx : ddx) + (degu > 0 ? degu : 1) + 1;
        long nwi = ihi - ilo + 1;
        size_t nH_syms = (size_t)(dH * (Nx + 1));   /* Hermite numerator coeffs */
        long nunk = nwi * (Nx + 1) + (long)nH_syms + (long)ng;
        if (!bad && nwi > 0 && nunk > 0) {
            size_t nw_syms = (size_t)(nwi * (Nx + 1));
            Expr** qterms = malloc((nw_syms + 1 + ng) * sizeof(Expr*));
            size_t ntq = 0;
            /* Laurent part: sum_i (sum_k rmW.. x^k) t^i. */
            for (long i = ilo; i <= ihi; i++) {
                for (long k = 0; k <= Nx; k++) {
                    char nm[32]; snprintf(nm, sizeof(nm), "rmW%ldv%ld", i - ilo, k);
                    qterms[ntq++] = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_symbol(nm),
                          expr_new_function(expr_new_symbol("Power"),
                            (Expr*[]){ expr_copy(x), expr_new_integer(k) }, 2),
                          expr_new_function(expr_new_symbol("Power"),
                            (Expr*[]){ expr_copy(tsym), expr_new_integer(i) }, 2) }, 3);
                }
            }
            /* Hermite part: (sum_{p<dH} (sum_k rmHh.. x^k) t^p) / Hden. */
            if (dH >= 1) {
                Expr** hterms = malloc(nH_syms * sizeof(Expr*));
                size_t nh = 0;
                for (long p = 0; p < dH; p++)
                    for (long k = 0; k <= Nx; k++) {
                        char nm[64]; snprintf(nm, sizeof(nm), "rmHh%ldv%ld", p, k);
                        hterms[nh++] = expr_new_function(expr_new_symbol("Times"),
                            (Expr*[]){ expr_new_symbol(nm),
                              expr_new_function(expr_new_symbol("Power"),
                                (Expr*[]){ expr_copy(x), expr_new_integer(k) }, 2),
                              expr_new_function(expr_new_symbol("Power"),
                                (Expr*[]){ expr_copy(tsym), expr_new_integer(p) }, 2) }, 3);
                    }
                Expr* Hpoly = expr_new_function(expr_new_symbol("Plus"), hterms, nh);
                free(hterms);
                qterms[ntq++] = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ Hpoly, expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(Hden), expr_new_integer(-1) }, 2) }, 2);
            }
            /* Log part: sum_j rmHc.. Log(g_j). */
            for (size_t j = 0; j < ng; j++) {
                char nm[32]; snprintf(nm, sizeof(nm), "rmHc%zu", j);
                qterms[ntq++] = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_symbol(nm),
                      expr_new_function(expr_new_symbol("Log"),
                        (Expr*[]){ expr_copy(g[j]) }, 1) }, 2);
            }
            Expr* Q = expr_new_function(expr_new_symbol("Plus"), qterms, ntq);
            free(qterms);

            /* Q'(full) = D[Q,x] + u' t D[Q,t]. */
            Expr* dQx = rt_eval2("D", expr_copy(Q), expr_copy(x));
            Expr* dQt = rt_eval2("D", expr_copy(Q), expr_copy(tsym));
            Expr* up = rt_eval2("D", expr_copy(uexp), expr_copy(x));
            Expr* Qder = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ dQx, expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ up, expr_copy(tsym), dQt }, 3) }, 2));
            Expr* diff = expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ Qder, expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(F) }, 2) }, 2);
            Expr* tog = rt_eval1("Together", diff);
            Expr* rnum = tog ? rt_eval1("Numerator", tog) : NULL;
            Expr* sol = NULL;
            if (rnum) {
                Expr* varlist = expr_new_function(expr_new_symbol("List"),
                    (Expr*[]){ expr_copy(tsym), expr_copy(x) }, 2);
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
                    Expr** zero = malloc((nw_syms + nH_syms + ng) * sizeof(Expr*));
                    size_t zi = 0;
                    for (long i = ilo; i <= ihi; i++)
                        for (long k = 0; k <= Nx; k++) {
                            char nm[32]; snprintf(nm, sizeof(nm), "rmW%ldv%ld", i - ilo, k);
                            zero[zi++] = expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_new_symbol(nm), expr_new_integer(0) }, 2);
                        }
                    for (long p = 0; p < dH; p++)
                        for (long k = 0; k <= Nx; k++) {
                            char nm[64]; snprintf(nm, sizeof(nm), "rmHh%ldv%ld", p, k);
                            zero[zi++] = expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_new_symbol(nm), expr_new_integer(0) }, 2);
                        }
                    for (size_t j = 0; j < ng; j++) {
                        char nm[32]; snprintf(nm, sizeof(nm), "rmHc%zu", j);
                        zero[zi++] = expr_new_function(expr_new_symbol("Rule"),
                            (Expr*[]){ expr_new_symbol(nm), expr_new_integer(0) }, 2);
                    }
                    Expr* zl = expr_new_function(expr_new_symbol("List"), zero,
                                                 nw_syms + nH_syms + ng);
                    free(zero);
                    Qs = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ Qs, zl }, 2));
                    if (Qs) {
                        Expr* back = expr_new_function(expr_new_symbol("List"),
                            (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_copy(tsym),
                                    expr_new_function(expr_new_symbol("Power"),
                                        (Expr*[]){ expr_new_symbol("E"), expr_copy(uexp) }, 2)
                                }, 2) }, 1);
                        result = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                            (Expr*[]){ Qs, back }, 2));
                    }
                }
            }
            if (sol) expr_free(sol);
            expr_free(Q);
        }
        for (size_t j = 0; j < ng; j++) expr_free(g[j]);
        if (factored) expr_free(factored);
        if (rad) expr_free(rad);
        if (Hden) expr_free(Hden);
        if (Dtil) expr_free(Dtil);
    }
    if (G) expr_free(G);
    if (num) expr_free(num);
    if (den) expr_free(den);
    expr_free(F);
    expr_free(tsym);
    expr_free(uexp);
    return result;
}

/* ================================================================== */
/* Case: sum of NON-COMMENSURATE exponentials (Phase B, first tower    */
/* increment) — several independent primitive exponential extensions.  */
/* ================================================================== */
/* Split one multiplicative term T = p(x) * prod_k E^(w_k) into the total
 * x-dependent exponential exponent W = sum_k w_k (over every E^w / Exp[w]
 * factor whose exponent depends on x) and the exponential-free cofactor p.
 * Returns 0 on success with *W_out and *p_out owned; -1 if the cofactor still
 * carries an x-dependent exponential (a shape this simple splitter declines). */
/* Recursively accumulate the factors of a term (descending into nested,
 * un-flattened Times nodes — Expand can leave Times[2, Times[x, E^w]]): every
 * x-dependent exponential exponent is summed into *W; every other factor is
 * appended (as an owned copy) to the growable array *pf. */
static void rt_accum_factors(Expr* T, Expr* x, Expr** W,
                             Expr*** pf, size_t* npf, size_t* cap) {
    if (T->type == EXPR_FUNCTION
        && T->data.function.head->type == EXPR_SYMBOL
        && T->data.function.head->data.symbol.name == intern_symbol("Times")) {
        for (size_t i = 0; i < T->data.function.arg_count; i++)
            rt_accum_factors(T->data.function.args[i], x, W, pf, npf, cap);
        return;
    }
    Expr* w = NULL;   /* borrowed exponent if T is an x-dependent exp kernel */
    const char* h = (T->type == EXPR_FUNCTION
        && T->data.function.head->type == EXPR_SYMBOL)
        ? T->data.function.head->data.symbol.name : NULL;
    if (h == intern_symbol("Exp") && T->data.function.arg_count == 1
        && !rt_free_of_x(T->data.function.args[0], x))
        w = T->data.function.args[0];
    else if (h == intern_symbol("Power") && T->data.function.arg_count == 2
        && T->data.function.args[0]->type == EXPR_SYMBOL
        && T->data.function.args[0]->data.symbol.name == intern_symbol("E")
        && !rt_free_of_x(T->data.function.args[1], x))
        w = T->data.function.args[1];
    if (w) {
        *W = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ *W, expr_copy(w) }, 2));
    } else {
        if (*npf == *cap) {
            *cap = *cap ? *cap * 2 : 4;
            *pf = realloc(*pf, *cap * sizeof(Expr*));
        }
        (*pf)[(*npf)++] = expr_copy(T);
    }
}

static int rt_split_exp_term(Expr* T, Expr* x, Expr** W_out, Expr** p_out) {
    Expr* W = expr_new_integer(0);
    Expr** pf = NULL; size_t npf = 0, cap = 0;
    rt_accum_factors(T, x, &W, &pf, &npf, &cap);

    Expr* p;
    if (npf == 0) p = expr_new_integer(1);
    else if (npf == 1) p = pf[0];                          /* adopt sole factor */
    else p = expr_new_function(expr_new_symbol("Times"), pf, npf);  /* adopts */
    free(pf);

    if (rt_find_exp_of_x(p, x) != NULL) {   /* nested/residual exp in cofactor */
        expr_free(W); expr_free(p);
        return -1;
    }
    *W_out = W; *p_out = p;
    return 0;
}

/* Sum-of-exponentials case.  When TrigToExp (or the raw integrand) yields a sum
 *   f = sum_k p_k(x) E^(W_k)
 * of exponentials whose exponents W_k are NOT integer multiples of a single
 * primitive (e.g. (1 +/- I) x from E^x Sin[x]), the distinct exponentials are
 * independent transcendental extensions.  Because d/dx maps each E^(W_k) to a
 * multiple of itself and never mixes distinct exponents, the terms DECOUPLE:
 *   INT p_k E^(W_k) dx = q_k E^(W_k),   q_k' + W_k' q_k = p_k  (a Risch DE)
 * with W_k polynomial in x, plus any exponential-free (W_k = 0) terms handled by
 * the base-field integral.  Each q_k is found exactly (SolveAlways certificate),
 * so d/dx(sum_k q_k E^(W_k)) = sum_k p_k E^(W_k) = f by linearity — correct by
 * construction.  Declines (whole) if any term is not rational-coefficient over a
 * polynomial exponent (e.g. E^(x^2), left to the Erf recognizer) or carries a
 * residual log / non-elementary cofactor.  Runs after the single-primitive
 * exponential cases, so it only fires on genuinely multi-kernel integrands. */
Expr* rt_expsum_case(Expr* f, Expr* x) {
    if (!rt_find_exp_of_x(f, x)) return NULL;   /* need at least one exp kernel */
    Expr* fe = rt_eval1("Expand", expr_copy(f));
    if (!fe) return NULL;

    Expr** terms; size_t nt; Expr* single[1];
    if (fe->type == EXPR_FUNCTION
        && fe->data.function.head->type == EXPR_SYMBOL
        && fe->data.function.head->data.symbol.name == intern_symbol("Plus")) {
        terms = fe->data.function.args; nt = fe->data.function.arg_count;
    } else { single[0] = fe; terms = single; nt = 1; }

    Expr** outs = malloc((nt ? nt : 1) * sizeof(Expr*));
    size_t no = 0;
    bool fail = false, any_exp = false;
    for (size_t i = 0; i < nt && !fail; i++) {
        Expr* W = NULL; Expr* p = NULL;
        if (rt_split_exp_term(terms[i], x, &W, &p) != 0) { fail = true; break; }
        bool wdep = !rt_free_of_x(W, x);
        Expr* qi = NULL;
        if (!wdep) {
            /* E^W is a constant coefficient: INT p E^W dx = E^W INT p dx. */
            qi = rt_integrate_in_K_with_logs(p, x);
        } else if (rt_is_poly(W, x)) {
            qi = rt_solve_rde(p, 1, W, x);   /* q' + W' q = p */
            if (qi) any_exp = true;
        }   /* rational (non-polynomial) W is Phase C: qi stays NULL -> decline */
        expr_free(p);
        if (!qi) { expr_free(W); fail = true; break; }
        Expr* eW = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_new_symbol("E"), W }, 2);   /* adopts W; E^0 -> 1 */
        outs[no++] = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ qi, eW }, 2);                    /* adopts qi, eW */
    }

    Expr* result = NULL;
    if (!fail && any_exp) {
        result = rt_eval_own(expr_new_function(expr_new_symbol("Plus"), outs, no));
    } else {
        for (size_t i = 0; i < no; i++) expr_free(outs[i]);
    }
    free(outs);
    expr_free(fe);
    return result;
}
