/* cherry_li.c — Cherry's logarithmic integral: LogIntegral (li) engine.
 *
 * G. W. Cherry, "Integration in Finite Terms with Special Functions: The
 * Logarithmic Integral" (SIAM J. Comput. 1986), single-logarithm tower case:
 * theta = Log[w] with w a polynomial in x.  The answer
 *   INT gamma dx = v(x, theta) + Sum_k d_k LogIntegral[w^k]
 * is an elementary Laurent part v = Sum_i b_i(x) theta^i plus logarithmic
 * integrals whose arguments are POWERS w^k (the degree-1 Sigma-decomposition over
 * the single generator w).  It is solved by matching IN THE TOWER: theta is an
 * independent variable, Log[w^k] = k theta by construction, and
 *   D[LogIntegral[w^k]] = (w^k)'/Log[w^k] = w^(k-1) w' / theta,
 * so the derivative of the answer is the rational (in x, theta) function
 *   D_tower[v] + Sum_k d_k w^(k-1) w'/theta,   D_tower = d/dx + (w'/w) d/dtheta.
 * Setting this equal to gamma and clearing denominators gives a linear system
 * over { coeffs(b_i), d_k } solved by SolveAlways over {theta, x}.  The exact
 * tower identity IS the certificate — a plain Simplify diff-back cannot reduce
 * Log[w^k] = k Log[w] (branch-cut safety), so the emitted answer is additionally
 * PowerExpand-verified.  Generalises the single-term rt_try_li recognizer to
 * MULTI-li (Cherry d3: INT x^2/Log[x+1] = li((x+1)^3) - 2 li((x+1)^2) + li(x+1)).
 *
 * Also handled: the transcendental-constant RESCALING (Cherry d2/d4) — a constant
 * nonzero root rho of the theta-denominator gives a rescaled li(e^(-rho) w), since
 * Log[w] - rho = Log[e^(-rho) w], so a simple pole 1/(theta - rho) integrates to
 * e^rho li(e^(-rho) w); the basis term is k w'/(theta - rho).
 *
 * Scope of THIS increment: a single log kernel with a polynomial argument w.
 * Multi-log towers with UNEQUAL factor exponents (a reducible w needing a genuine
 * product Sigma-decomposition Prod f_j^alpha_j, e.g. Log[x^3-x]) and the
 * Sigma-decomposition NON-existence decision (Ex 5.2) are later increments.
 * Correct by construction (tower solve + PowerExpand diff-back), so a
 * mis-generation can only decline.
 */

#include "cherry_li.h"
#include "risch_util.h"

#include "expr.h"
#include "eval.h"
#include "print.h"
#include "sym_intern.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static Expr* mk_sym(const char* s)  { return expr_new_symbol(s); }
static Expr* mk_int(long n)         { return expr_new_integer(n); }
static Expr* mk_pow(Expr* b, Expr* e) {
    return expr_new_function(mk_sym("Power"), (Expr*[]){ b, e }, 2);
}
static Expr* mk_neg(Expr* a) {
    return expr_new_function(mk_sym("Times"), (Expr*[]){ mk_int(-1), a }, 2);
}
static Expr* mk_plus2(Expr* a, Expr* b) {
    return expr_new_function(mk_sym("Plus"), (Expr*[]){ a, b }, 2);
}
static Expr* mk_times2(Expr* a, Expr* b) {
    return expr_new_function(mk_sym("Times"), (Expr*[]){ a, b }, 2);
}

Expr* rt_cherry_li(Expr* f, Expr* x) {
    /* 1. The single log kernel theta = Log[w], w a polynomial in x. */
    Expr* w = rt_find_log_of_x(f, x);                    /* borrowed */
    if (!w || !rt_kernel_simple(w, x) || !rt_is_poly(w, x)) return NULL;
    if (rt_degree(w, x) < 1) return NULL;

    Expr* t = mk_sym("chli$t");                          /* theta */
    Expr* logw = expr_new_function(mk_sym("Log"), (Expr*[]){ expr_copy(w) }, 1);
    Expr* rule = expr_new_function(mk_sym("Rule"),
        (Expr*[]){ logw, expr_copy(t) }, 2);
    Expr* F = rt_eval_own(expr_new_function(mk_sym("ReplaceAll"),
        (Expr*[]){ expr_copy(f), rule }, 2));            /* f with Log[w] -> t */
    if (!F || !rt_free_of_head(F, "Log") || rt_free_of_x(F, t)
        || rt_find_exp_of_x(F, x)) {
        if (F) expr_free(F);
        expr_free(t);
        return NULL;
    }

    /* 2. theta-denominator multiplicity m and x-degree of the numerator. */
    Expr* Ft   = rt_eval1("Together", expr_copy(F));
    Expr* numF = Ft ? rt_eval1("Numerator", expr_copy(Ft)) : NULL;
    Expr* denF = Ft ? rt_eval1("Denominator", expr_copy(Ft)) : NULL;
    if (Ft) expr_free(Ft);
    long m = denF ? rt_degree(denF, t) : 0;
    long dxF = numF ? rt_degree(numF, x) : 0;
    long dw  = rt_degree(w, x);
    if (numF) expr_free(numF);

    /* Constant nonzero roots rho of the theta-denominator give transcendental-
     * constant RESCALED li terms: a simple pole 1/(theta - rho) integrates to
     * e^rho li(e^(-rho) w), since Log[w] - rho = Log[e^(-rho) w] (Cherry d2/d4:
     * 1/(Log[x]+3) -> e^-3 li(e^3 x)). */
    Expr* shifts[16]; size_t ns = 0;
    if (denF && m >= 1) {
        Expr* eqd = expr_new_function(mk_sym("Equal"),
            (Expr*[]){ expr_copy(denF), mk_int(0) }, 2);
        Expr* rs = rt_eval2("Solve", eqd, expr_copy(t));
        if (rs && rs->type == EXPR_FUNCTION && rt_head_is(rs, "List")) {
            for (size_t i = 0; i < rs->data.function.arg_count && ns < 16; i++) {
                Expr* rule = rs->data.function.args[i];
                if (rule->type != EXPR_FUNCTION || !rt_head_is(rule, "List")
                    || rule->data.function.arg_count != 1) continue;
                Expr* rr = rule->data.function.args[0];
                if (rr->type != EXPR_FUNCTION || !rt_head_is(rr, "Rule")) continue;
                Expr* val = rr->data.function.args[1];
                if (!rt_free_of_x(val, x) || !rt_free_of_x(val, t) || rt_is_zero(val))
                    continue;                       /* keep constant, nonzero roots */
                bool dup = false;
                for (size_t j = 0; j < ns; j++)
                    if (expr_eq(shifts[j], val)) { dup = true; break; }
                if (!dup) shifts[ns++] = expr_copy(val);
            }
        }
        if (rs) expr_free(rs);
    }
    if (denF) expr_free(denF);
    if (m < 1) {
        for (size_t i = 0; i < ns; i++) expr_free(shifts[i]);
        expr_free(F); expr_free(t); return NULL;                 /* no 1/theta -> no li */
    }

    /* 3. tower data: w', eta = theta' = w'/w. */
    Expr* wp  = rt_eval2("D", expr_copy(w), expr_copy(x));
    Expr* eta = rt_eval1("Cancel", mk_times2(expr_copy(wp),
        mk_pow(expr_copy(w), mk_int(-1))));

    /* 4. Ansatz bounds (derived, generous — the tower identity certifies soundness).
     *    Laurent theta powers i in [-(m-1), 1]; b_i in x of degree <= Dx; li powers
     *    k in [1, K]. */
    long lo = -(m - 1), hi = 1;
    long Dx = dxF + dw + 2;               if (Dx < 0) Dx = 0;
    long K  = (dw > 0 ? dxF / dw : dxF) + 3; if (K < 2) K = 2;

    size_t nlau = (size_t)(hi - lo + 1);
    size_t nsym = nlau * (size_t)(Dx + 1) + (size_t)K + ns;
    Expr** syms = malloc(nsym * sizeof(Expr*));
    size_t si = 0;

    /* v = Sum_{i=lo}^{hi} b_i(x) theta^i,  b_i = Sum_j a_ij x^j.
     * D_tower[v] = Sum_i ( b_i' theta^i + i b_i eta theta^(i-1) ). */
    Expr** dterms = malloc((2 * nlau + (size_t)K + ns) * sizeof(Expr*));
    size_t nd = 0;
    Expr** vterms = malloc(nlau * sizeof(Expr*));         /* for the emitted answer */
    for (long i = lo; i <= hi; i++) {
        Expr** bt = malloc((size_t)(Dx + 1) * sizeof(Expr*));
        Expr** bt2 = malloc((size_t)(Dx + 1) * sizeof(Expr*));
        for (long j = 0; j <= Dx; j++) {
            char nm[64]; snprintf(nm, sizeof(nm), "chli$a%ld_%ld", i - lo, j);
            Expr* a = mk_sym(nm);
            syms[si++] = expr_copy(a);
            bt[j]  = mk_times2(expr_copy(a), mk_pow(expr_copy(x), mk_int(j)));
            bt2[j] = mk_times2(a, mk_pow(expr_copy(x), mk_int(j)));
        }
        Expr* bi  = expr_new_function(mk_sym("Plus"), bt,  (size_t)(Dx + 1));
        Expr* biA = expr_new_function(mk_sym("Plus"), bt2, (size_t)(Dx + 1));
        free(bt); free(bt2);
        /* emitted v term: b_i theta^i */
        vterms[i - lo] = mk_times2(biA, mk_pow(expr_copy(t), mk_int(i)));
        /* derivative terms */
        Expr* bip = rt_eval2("D", expr_copy(bi), expr_copy(x));
        dterms[nd++] = mk_times2(bip, mk_pow(expr_copy(t), mk_int(i)));
        if (i != 0)
            dterms[nd++] = mk_times2(mk_times2(mk_int(i), bi),
                mk_times2(expr_copy(eta), mk_pow(expr_copy(t), mk_int(i - 1))));
        else
            expr_free(bi);
    }

    /* li terms: contribution d_k w^(k-1) w' / theta, and record the answer terms. */
    Expr** dsyms = malloc((size_t)K * sizeof(Expr*));
    Expr** liargs = malloc((size_t)K * sizeof(Expr*));
    for (long k = 1; k <= K; k++) {
        char nm[32]; snprintf(nm, sizeof(nm), "chli$d%ld", k);
        Expr* d = mk_sym(nm);
        dsyms[k - 1] = expr_copy(d);
        syms[si++] = expr_copy(d);
        liargs[k - 1] = mk_pow(expr_copy(w), mk_int(k));       /* w^k */
        Expr* term = mk_times2(d, mk_times2(
            mk_times2(mk_pow(expr_copy(w), mk_int(k - 1)), expr_copy(wp)),
            mk_pow(expr_copy(t), mk_int(-1))));
        dterms[nd++] = term;
    }

    /* rescaled-li (shift) terms: k_i w'/(theta - rho_i). */
    Expr** ssyms = malloc((ns ? ns : 1) * sizeof(Expr*));
    for (size_t i = 0; i < ns; i++) {
        char nm[32]; snprintf(nm, sizeof(nm), "chli$s%zu", i);
        Expr* k = mk_sym(nm);
        ssyms[i] = expr_copy(k);
        syms[si++] = expr_copy(k);
        Expr* den = mk_plus2(expr_copy(t), mk_neg(expr_copy(shifts[i])));   /* theta - rho */
        dterms[nd++] = mk_times2(k, mk_times2(expr_copy(wp), mk_pow(den, mk_int(-1))));
    }

    Expr* rhs = expr_new_function(mk_sym("Plus"), dterms, nd);
    free(dterms);

    /* resid = F - rhs;  its numerator must vanish identically in {theta, x}. */
    Expr* resid = mk_plus2(expr_copy(F), mk_neg(rhs));
    Expr* rnum = rt_eval1("Numerator", rt_eval1("Together", resid));
    Expr* sol = NULL;
    if (rnum) {
        Expr* eqn = expr_new_function(mk_sym("Equal"), (Expr*[]){ rnum, mk_int(0) }, 2);
        Expr* vl = expr_new_function(mk_sym("List"),
            (Expr*[]){ expr_copy(t), expr_copy(x) }, 2);
        sol = rt_eval2("SolveAlways", eqn, vl);
    }

    /* 5. Assemble and PowerExpand-verify. */
    Expr* result = NULL;
    bool solved = sol && sol->type == EXPR_FUNCTION && rt_head_is(sol, "List")
        && sol->data.function.arg_count >= 1
        && sol->data.function.args[0]->type == EXPR_FUNCTION
        && rt_head_is(sol->data.function.args[0], "List");
    if (solved) {
        Expr* rules = sol->data.function.args[0];
        size_t nans = nlau + (size_t)K + ns;
        Expr** ans = malloc(nans * sizeof(Expr*));
        for (size_t i = 0; i < nlau; i++) ans[i] = expr_copy(vterms[i]);
        for (long k = 1; k <= K; k++)
            ans[nlau + (size_t)(k - 1)] = mk_times2(expr_copy(dsyms[k - 1]),
                expr_new_function(mk_sym("LogIntegral"),
                    (Expr*[]){ expr_copy(liargs[k - 1]) }, 1));
        /* rescaled-li answer: k_i e^rho_i LogIntegral[e^(-rho_i) w] */
        for (size_t i = 0; i < ns; i++) {
            Expr* erho  = mk_pow(mk_sym("E"), expr_copy(shifts[i]));
            Expr* liarg = mk_times2(mk_pow(mk_sym("E"), mk_neg(expr_copy(shifts[i]))),
                                    expr_copy(w));
            Expr* li = expr_new_function(mk_sym("LogIntegral"), (Expr*[]){ liarg }, 1);
            ans[nlau + (size_t)K + i] = mk_times2(expr_copy(ssyms[i]),
                                                  mk_times2(erho, li));
        }
        Expr* Q = expr_new_function(mk_sym("Plus"), ans, nans);
        free(ans);
        /* substitute solution, pin free unknowns to 0, then theta -> Log[w] */
        Q = rt_eval_own(expr_new_function(mk_sym("ReplaceAll"),
            (Expr*[]){ Q, expr_copy(rules) }, 2));
        if (Q) {
            Expr** zero = malloc(nsym * sizeof(Expr*));
            for (size_t j = 0; j < nsym; j++)
                zero[j] = expr_new_function(mk_sym("Rule"),
                    (Expr*[]){ expr_copy(syms[j]), mk_int(0) }, 2);
            Expr* zl = expr_new_function(mk_sym("List"), zero, nsym);
            free(zero);
            Q = rt_eval_own(expr_new_function(mk_sym("ReplaceAll"),
                (Expr*[]){ Q, zl }, 2));
        }
        if (Q) {
            Expr* tr = expr_new_function(mk_sym("Rule"),
                (Expr*[]){ expr_copy(t), expr_new_function(mk_sym("Log"),
                    (Expr*[]){ expr_copy(w) }, 1) }, 2);
            Q = rt_eval_own(expr_new_function(mk_sym("ReplaceAll"),
                (Expr*[]){ Q, tr }, 2));
        }
        /* PowerExpand diff-back (Log[w^k] = k Log[w], the li convention). */
        if (Q && rt_free_of_head(Q, "Integrate")) {
            Expr* diff = mk_plus2(rt_eval2("D", expr_copy(Q), expr_copy(x)),
                                  mk_neg(expr_copy(f)));
            Expr* chk = rt_eval1("Simplify", rt_eval1("PowerExpand", diff));
            if (chk && chk->type == EXPR_INTEGER && chk->data.integer == 0)
                result = Q;
            else
                expr_free(Q);
            if (chk) expr_free(chk);
        } else if (Q) {
            expr_free(Q);
        }
    }
    if (sol) expr_free(sol);

    for (size_t i = 0; i < nlau; i++) expr_free(vterms[i]);
    free(vterms);
    for (long k = 0; k < K; k++) { expr_free(dsyms[k]); expr_free(liargs[k]); }
    free(dsyms); free(liargs);
    for (size_t i = 0; i < ns; i++) { expr_free(shifts[i]); expr_free(ssyms[i]); }
    free(ssyms);
    for (size_t j = 0; j < nsym; j++) expr_free(syms[j]);
    free(syms);
    expr_free(F); expr_free(t); expr_free(wp); expr_free(eta);
    return result;
}
