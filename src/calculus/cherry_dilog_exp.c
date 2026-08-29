/* cherry_dilog_exp.c — exponential-tower dilogarithm integration (Cherry).
 *
 * The exponential-tower mirror of rt_cherry_dilog (cherry_dilog.c).  With the
 * single exponential kernel theta = E^(c x) (theta' = c theta), x = Log[theta]/c
 * is the "root at 0" tower-log and Log[theta - rho_k] are the other tower-logs
 * (so Log[1 + E^x] = Log[theta - (-1)]).  The engine handles integrands that are
 * a rational-in-theta combination of the weight-1 logs {x} u {Log[theta-rho_k]},
 * of total weight 1 — a SINGLE input form covering both x/(E^x-1) (weight-1 log
 * = x, coeff rational in theta) and Log[1+E^x] (weight-1 log = the outer log).
 *
 * The answer is a weight-2 combination
 *   Q = Sum c_ij Log[theta-rho_i] Log[theta-rho_j]
 *       + Sum (a_k x + b_k) Log[theta-rho_k]
 *       + Sum d_g PolyLog[2, g]  +  e x^2/2 + ...,
 * with g a Moebius function of theta between the special points {rho_k} u {0,∞}
 * such that both g and 1-g factor over {theta} u {theta-rho_k}.  Matching IN THE
 * TOWER: differentiate each candidate (theta known to the evaluator via
 * theta = E^(c x)), normalise every log with Log[z] -> PowerExpand[Log[Factor[
 * Together[z]]]] (Together rewrites E^-x = 1/theta so Log[(theta-1)/theta] splits
 * into Log[theta-1] - Log[theta]), substitute Log[theta-rho_k] -> u_k and
 * Log[theta] -> c x, and require the residual to vanish as a polynomial identity
 * in {theta, x, u_1,...,u_m}.  The exact tower identity is the certificate,
 * re-checked by a PowerExpand diff-back.  Any mis-generation declines (NULL).
 *
 * Scope: single commensurable exponential kernel, rational roots rho_k,
 * weight-1 integrand.  Declines cleanly: weight >= 2 (-> PolyLog[3]),
 * incommensurate exponentials, algebraic roots, and negative-argument outer logs
 * (Log[1-E^x]) whose reversed-pair branch injects i pi (dropped by has_bad_log).
 */

#include "cherry_dilog_exp.h"
#include "risch_util.h"
#include "risch_singleext.h"     /* rt_exp_kernelize */
#include "integrate.h"           /* g_integrate_depth (top-level-only gating) */

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
static Expr* mk_log(Expr* a) { return expr_new_function(mk_sym("Log"), (Expr*[]){ a }, 1); }
static Expr* mk_fn2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(mk_sym(h), (Expr*[]){ a, b }, 2);
}

/* Collect the distinct theta-dependent Log[p] arguments in e into ws[] (owned
 * copies, deduped), up to cap. */
static void collect_logs_th(Expr* e, Expr* th, Expr** ws, size_t* n, size_t cap) {
    if (!e || e->type != EXPR_FUNCTION) return;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == intern_symbol("Log")
        && e->data.function.arg_count == 1) {
        Expr* w = e->data.function.args[0];
        if (!rt_free_of_x(w, th)) {
            for (size_t i = 0; i < *n; i++) if (expr_eq(ws[i], w)) return;
            if (*n < cap) ws[(*n)++] = expr_copy(w);
            return;
        }
    }
    collect_logs_th(e->data.function.head, th, ws, n, cap);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        collect_logs_th(e->data.function.args[i], th, ws, n, cap);
}

/* Rational nonzero roots of a polynomial p in th (via Solve), appended deduped. */
static void collect_rational_roots_th(Expr* p, Expr* th, Expr** rs, size_t* n, size_t cap) {
    if (!p || rt_degree(p, th) < 1) return;
    Expr* sols = rt_eval2("Solve",
        mk_fn2("Equal", expr_copy(p), mk_int(0)), expr_copy(th));
    if (sols && sols->type == EXPR_FUNCTION && rt_head_is(sols, "List")) {
        for (size_t i = 0; i < sols->data.function.arg_count && *n < cap; i++) {
            Expr* rule = sols->data.function.args[i];
            if (rule->type != EXPR_FUNCTION || !rt_head_is(rule, "List")
                || rule->data.function.arg_count != 1) continue;
            Expr* rr = rule->data.function.args[0];
            if (rr->type != EXPR_FUNCTION || !rt_head_is(rr, "Rule")) continue;
            Expr* val = rr->data.function.args[1];
            bool rat = (val->type == EXPR_INTEGER || val->type == EXPR_BIGINT)
                || (val->type == EXPR_FUNCTION && rt_head_is(val, "Rational"));
            if (!rat || rt_is_zero(val)) continue;
            bool dup = false;
            for (size_t k = 0; k < *n; k++) if (expr_eq(rs[k], val)) { dup = true; break; }
            if (!dup) rs[(*n)++] = expr_copy(val);
        }
    }
    if (sols) expr_free(sols);
}

/* True if the normalised candidate derivative e carries a disqualifying log:
 *  (a) Complex or Pi (a reversed root pair whose 1-g < 0 leaves Log[-1]=i pi), or
 *  (b) a residual Log of a theta- or x-dependent argument (genuinely unmatched).
 * A Log of a positive constant (real transcendental Log-Log term) is kept. */
static bool has_bad_log_exp(Expr* e, Expr* x, Expr* th) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == intern_symbol("Pi");
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == intern_symbol("Complex")) return true;
        if (h == intern_symbol("Log") && e->data.function.arg_count == 1) {
            Expr* arg = e->data.function.args[0];
            if (!rt_free_of_x(arg, x) || !rt_free_of_x(arg, th)) return true;
        }
    }
    if (has_bad_log_exp(e->data.function.head, x, th)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_bad_log_exp(e->data.function.args[i], x, th)) return true;
    return false;
}

/* norm[e]: E^(k x) -> theta^(k/c), then Log[z] -> PowerExpand[Log[Factor[
 * Together[z]]]], then //. urules (Log[theta-rho_k] -> u_k, Log[theta] -> c x).
 * Adopts e; copies the rule exprs; returns an owned expr or NULL. */
static Expr* norm_exp(Expr* e, Expr* ekrule, Expr* logrule, Expr* urules) {
    Expr* a = rt_eval_own(mk_fn2("ReplaceAll", e, expr_copy(ekrule)));
    Expr* b = a ? rt_eval_own(mk_fn2("ReplaceAll", a, expr_copy(logrule))) : NULL;
    Expr* c = b ? rt_eval_own(mk_fn2("ReplaceRepeated", b, expr_copy(urules))) : NULL;
    return c;
}

Expr* rt_cherry_dilog_exp(Expr* f, Expr* x) {
    /* Fire only on the ORIGINAL top-level integrand (g_integrate_depth <= 1).  A
     * DerivativeDivides substitution (u = Log[x], depth >= 2) can turn a genuine
     * LOG-tower dilog such as Log[x]/(1-x) into the exp-tower shape -(u E^u)/(E^u
     * -1); handling it here would pre-empt the direct log-tower engine
     * rt_cherry_dilog, which gives the cleaner PolyLog[2, 1-x] form.  Declining in
     * nested recursion lets that engine win, exactly as the nonelem message is
     * gated to depth <= 1 (integrate_risch_transcendental.c). */
    if (g_integrate_depth > 1) return NULL;
    /* 1. kernelize E^(c x) -> theta = rmT; derive the rate c from the primitive. */
    if (!rt_find_exp_of_x(f, x)) return NULL;
    Expr* u = NULL;
    Expr* F = rt_exp_kernelize(f, x, &u);          /* F in rmT, u the primitive */
    if (!F || !u) { if (F) expr_free(F); if (u) expr_free(u); return NULL; }
    Expr* th = mk_sym("rmT");
    if (!rt_is_poly(u, x) || rt_degree(u, x) != 1) {
        expr_free(F); expr_free(u); expr_free(th); return NULL;
    }
    Expr* u0 = rt_coeff(u, x, 0);
    bool linok = u0 && rt_is_zero(u0);
    if (u0) expr_free(u0);
    Expr* c = rt_coeff(u, x, 1);                    /* rate c (owned) */
    expr_free(u);
    if (!linok || !c || rt_is_zero(c)) {
        if (c) expr_free(c);
        expr_free(F); expr_free(th); return NULL;
    }
    /* F must actually involve theta and be free of exp/log of x other than the
     * kernel we just peeled. */
    if (rt_free_of_x(F, th) || rt_find_exp_of_x(F, x)) {
        expr_free(F); expr_free(th); expr_free(c); return NULL;
    }

    /* 2. outer theta-logs and their roots + denominator roots. */
    Expr* ws[8]; size_t nw = 0;
    collect_logs_th(F, th, ws, &nw, 8);

    Expr* rs[16]; size_t nr = 0;
    /* denominator roots: F with every Log[_] -> 0, then Denominator[Together]. */
    Expr* zero_logs = expr_new_function(mk_sym("RuleDelayed"),
        (Expr*[]){ mk_log(expr_new_function(mk_sym("Blank"), NULL, 0)), mk_int(0) }, 2);
    Expr* Fnl = rt_eval_own(mk_fn2("ReplaceAll", expr_copy(F), zero_logs));
    if (Fnl) {
        Expr* den = rt_eval1("Denominator", rt_eval1("Together", Fnl));
        if (den) { collect_rational_roots_th(den, th, rs, &nr, 16); expr_free(den); }
    }
    for (size_t i = 0; i < nw; i++) collect_rational_roots_th(ws[i], th, rs, &nr, 16);
    if (nr == 0) {
        for (size_t i = 0; i < nw; i++) expr_free(ws[i]);
        expr_free(F); expr_free(th); expr_free(c); return NULL;
    }

    /* 3. kernel logs u_k = Log[theta - rho_k] -> chde$u%zu; and Log[theta] -> c x. */
    Expr** facs = malloc(nr * sizeof(Expr*));      /* theta - rho_k */
    Expr** us   = malloc(nr * sizeof(Expr*));      /* kernel vars */
    for (size_t k = 0; k < nr; k++) {
        facs[k] = rt_eval_own(mk_plus2(expr_copy(th), mk_neg(expr_copy(rs[k]))));
        char nm[24]; snprintf(nm, sizeof(nm), "chde$u%zu", k);
        us[k] = mk_sym(nm);
    }
    /* urules = { Log[theta-rho_k] -> u_k, ..., Log[theta] -> c x } */
    Expr** url = malloc((nr + 1) * sizeof(Expr*));
    for (size_t k = 0; k < nr; k++)
        url[k] = mk_fn2("Rule", mk_log(expr_copy(facs[k])), expr_copy(us[k]));
    url[nr] = mk_fn2("Rule", mk_log(expr_copy(th)), mk_times2(expr_copy(c), expr_copy(x)));
    Expr* urules = expr_new_function(mk_sym("List"), url, nr + 1);
    free(url);

    /* ekrule = Power[E, ee_] :> theta^(Coefficient[ee,x,1]/c) */
    Expr* ee = mk_sym("chde$ee");
    Expr* ekrule = expr_new_function(mk_sym("RuleDelayed"),
        (Expr*[]){ mk_pow(mk_sym("E"),
                     expr_new_function(mk_sym("Pattern"),
                        (Expr*[]){ expr_copy(ee),
                                   expr_new_function(mk_sym("Blank"), NULL, 0) }, 2)),
                   mk_pow(expr_copy(th),
                     mk_times2(expr_new_function(mk_sym("Coefficient"),
                                  (Expr*[]){ expr_copy(ee), expr_copy(x), mk_int(1) }, 3),
                               mk_pow(expr_copy(c), mk_int(-1)))) }, 2);
    expr_free(ee);
    /* logrule = Log[zz_] :> PowerExpand[Log[Factor[Together[zz]]]] */
    Expr* zz = mk_sym("chde$zz");
    Expr* logrule = expr_new_function(mk_sym("RuleDelayed"),
        (Expr*[]){ mk_log(expr_new_function(mk_sym("Pattern"),
                        (Expr*[]){ expr_copy(zz),
                                   expr_new_function(mk_sym("Blank"), NULL, 0) }, 2)),
                   rt_eval1("HoldForm", mk_int(0)) /* placeholder, replaced below */ }, 2);
    /* build the RHS PowerExpand[Log[Factor[Together[zz]]]] unevaluated */
    {
        Expr* rhs = expr_new_function(mk_sym("PowerExpand"),
            (Expr*[]){ mk_log(expr_new_function(mk_sym("Factor"),
                        (Expr*[]){ expr_new_function(mk_sym("Together"),
                                     (Expr*[]){ expr_copy(zz) }, 1) }, 1)) }, 1);
        expr_free(logrule->data.function.args[1]);
        logrule->data.function.args[1] = rhs;
    }
    expr_free(zz);

    /* 4. candidate answer terms + admissibility filter. */
    size_t cap = 2 * nr * nr + 2 * nr + nr * (nr + 1) / 2 + 4;
    Expr** basis  = malloc(cap * sizeof(Expr*));   /* answer term (E^(c x) form) */
    Expr** dbasis = malloc(cap * sizeof(Expr*));   /* norm[D[sym*term]]          */
    Expr** syms   = malloc(cap * sizeof(Expr*));
    size_t nc = 0;

    /* rmT -> E^(c x) substitution for turning a theta-form term into a real
     * differentiable expression. */
    Expr* th2exp = mk_fn2("Rule", expr_copy(th),
        mk_pow(mk_sym("E"), mk_times2(expr_copy(c), expr_copy(x))));

    #define ADD_TERM(theta_term)                                                    \
        do {                                                                        \
            Expr* real = rt_eval_own(mk_fn2("ReplaceAll", (theta_term),             \
                                            expr_copy(th2exp)));                     \
            if (real) {                                                             \
                char sn[24]; snprintf(sn, sizeof(sn), "chde$c%zu", nc);             \
                Expr* s = mk_sym(sn);                                               \
                Expr* term = mk_times2(expr_copy(s), real);                         \
                Expr* d = norm_exp(rt_eval2("D", expr_copy(term), expr_copy(x)),    \
                                   ekrule, logrule, urules);                        \
                if (d && !has_bad_log_exp(d, x, th)) {                              \
                    syms[nc] = s; basis[nc] = term; dbasis[nc] = d; nc++;           \
                } else {                                                            \
                    expr_free(s); expr_free(term); if (d) expr_free(d);             \
                }                                                                   \
            }                                                                       \
        } while (0)

    /* PolyLog[2, g] — Moebius candidates from families A-D. */
    for (size_t i = 0; i < nr; i++) {
        /* (A) rho_i/theta ; (B) theta/rho_i */
        ADD_TERM(expr_new_function(mk_sym("PolyLog"),
            (Expr*[]){ mk_int(2), rt_eval1("Together",
                mk_times2(expr_copy(rs[i]), mk_pow(expr_copy(th), mk_int(-1)))) }, 2));
        ADD_TERM(expr_new_function(mk_sym("PolyLog"),
            (Expr*[]){ mk_int(2), rt_eval1("Together",
                mk_times2(expr_copy(th), mk_pow(expr_copy(rs[i]), mk_int(-1)))) }, 2));
    }
    for (size_t i = 0; i < nr; i++)
        for (size_t j = 0; j < nr; j++) {
            if (i == j) continue;
            /* (C) affine (theta-rho_i)/(rho_j-rho_i) */
            ADD_TERM(expr_new_function(mk_sym("PolyLog"),
                (Expr*[]){ mk_int(2), rt_eval1("Together", mk_times2(
                    mk_plus2(expr_copy(th), mk_neg(expr_copy(rs[i]))),
                    mk_pow(mk_plus2(expr_copy(rs[j]), mk_neg(expr_copy(rs[i]))),
                           mk_int(-1)))) }, 2));
            /* (D) root-to-root (theta-rho_i)/(theta-rho_j) */
            ADD_TERM(expr_new_function(mk_sym("PolyLog"),
                (Expr*[]){ mk_int(2), rt_eval1("Together", mk_times2(
                    mk_plus2(expr_copy(th), mk_neg(expr_copy(rs[i]))),
                    mk_pow(mk_plus2(expr_copy(th), mk_neg(expr_copy(rs[j]))),
                           mk_int(-1)))) }, 2));
        }
    /* cross terms x Log[theta-rho_k] and single Log[theta-rho_k]. */
    for (size_t k = 0; k < nr; k++) {
        ADD_TERM(mk_times2(expr_copy(x), mk_log(expr_copy(facs[k]))));
        ADD_TERM(mk_log(expr_copy(facs[k])));
    }
    /* Log-Log products. */
    for (size_t i = 0; i < nr; i++)
        for (size_t j = i; j < nr; j++)
            ADD_TERM(mk_times2(mk_log(expr_copy(facs[i])), mk_log(expr_copy(facs[j]))));
    /* root-at-0 monomials x^2/2 and x. */
    ADD_TERM(mk_times2(mk_pow(expr_copy(x), mk_int(2)),
                       expr_new_function(mk_sym("Rational"),
                          (Expr*[]){ mk_int(1), mk_int(2) }, 2)));
    ADD_TERM(expr_copy(x));
    #undef ADD_TERM
    expr_free(th2exp);

    /* 5. match: norm[F] == Sum dbasis, as an identity in {theta, x, u_1,...,u_m}. */
    Expr* Fk = norm_exp(expr_copy(F), ekrule, logrule, urules);
    Expr* sol = NULL;
    if (Fk && nc > 0) {
        Expr** rhsarr = malloc(nc * sizeof(Expr*));
        for (size_t i = 0; i < nc; i++) rhsarr[i] = expr_copy(dbasis[i]);
        Expr* rhs_sum = expr_new_function(mk_sym("Plus"), rhsarr, nc);
        free(rhsarr);
        Expr* resid = mk_plus2(expr_copy(Fk), mk_neg(rhs_sum));
        Expr* rnum = rt_eval1("Numerator", rt_eval1("Together", resid));
        if (rnum) {
            Expr** vl = malloc((nr + 2) * sizeof(Expr*));
            vl[0] = expr_copy(th);
            vl[1] = expr_copy(x);
            for (size_t k = 0; k < nr; k++) vl[2 + k] = expr_copy(us[k]);
            Expr* varlist = expr_new_function(mk_sym("List"), vl, nr + 2);
            free(vl);
            Expr* clist = rt_eval2("CoefficientList", rnum, varlist);
            Expr* threaded = rt_eval1("Thread", mk_fn2("Equal",
                rt_eval1("Flatten", clist), mk_int(0)));
            Expr* eqs = rt_eval2("DeleteCases", threaded, mk_sym("True"));
            Expr** ul = malloc(nc * sizeof(Expr*));
            for (size_t i = 0; i < nc; i++) ul[i] = expr_copy(syms[i]);
            Expr* unklist = expr_new_function(mk_sym("List"), ul, nc);
            free(ul);
            sol = rt_eval2("Solve", eqs, unklist);
        }
    }
    if (Fk) expr_free(Fk);

    /* 6. assemble and PowerExpand diff-back verify. */
    Expr* result = NULL;
    bool solved = sol && sol->type == EXPR_FUNCTION && rt_head_is(sol, "List")
        && sol->data.function.arg_count >= 1
        && sol->data.function.args[0]->type == EXPR_FUNCTION
        && rt_head_is(sol->data.function.args[0], "List");
    if (solved) {
        Expr* rulesol = sol->data.function.args[0];
        Expr** ans = malloc(nc * sizeof(Expr*));
        for (size_t i = 0; i < nc; i++) ans[i] = expr_copy(basis[i]);
        Expr* Q = expr_new_function(mk_sym("Plus"), ans, nc);
        free(ans);
        Q = rt_eval_own(mk_fn2("ReplaceAll", Q, expr_copy(rulesol)));
        if (Q) {
            Expr** zero = malloc(nc * sizeof(Expr*));
            for (size_t i = 0; i < nc; i++)
                zero[i] = mk_fn2("Rule", expr_copy(syms[i]), mk_int(0));
            Expr* zl = expr_new_function(mk_sym("List"), zero, nc);
            free(zero);
            Q = rt_eval_own(mk_fn2("ReplaceAll", Q, zl));
        }
        if (Q && rt_free_of_head(Q, "Integrate")) {
            /* diff = (D[Q,x] - f) /. Log[a_] :> Log[Factor[Together[a]]] */
            Expr* diff = mk_plus2(rt_eval2("D", expr_copy(Q), expr_copy(x)),
                                  mk_neg(expr_copy(f)));
            Expr* aa = mk_sym("chde$aa");
            Expr* lfr = expr_new_function(mk_sym("RuleDelayed"),
                (Expr*[]){ mk_log(expr_new_function(mk_sym("Pattern"),
                              (Expr*[]){ expr_copy(aa),
                                         expr_new_function(mk_sym("Blank"), NULL, 0) }, 2)),
                           mk_log(expr_new_function(mk_sym("Factor"),
                              (Expr*[]){ expr_new_function(mk_sym("Together"),
                                           (Expr*[]){ expr_copy(aa) }, 1) }, 1)) }, 2);
            expr_free(aa);
            Expr* diff2 = rt_eval_own(mk_fn2("ReplaceAll", diff, lfr));
            Expr* chk = diff2 ? rt_eval1("Simplify", rt_eval1("PowerExpand", diff2)) : NULL;
            if (chk && chk->type == EXPR_INTEGER && chk->data.integer == 0) result = Q;
            else expr_free(Q);
            if (chk) expr_free(chk);
        } else if (Q) expr_free(Q);
    }
    if (sol) expr_free(sol);

    /* 7. cleanup. */
    for (size_t i = 0; i < nc; i++) { expr_free(basis[i]); expr_free(dbasis[i]); expr_free(syms[i]); }
    free(basis); free(dbasis); free(syms);
    for (size_t k = 0; k < nr; k++) { expr_free(facs[k]); expr_free(us[k]); }
    free(facs); free(us);
    for (size_t i = 0; i < nr; i++) expr_free(rs[i]);
    for (size_t i = 0; i < nw; i++) expr_free(ws[i]);
    expr_free(urules); expr_free(ekrule); expr_free(logrule);
    expr_free(F); expr_free(th); expr_free(c);
    return result;
}

/* Integrate`Cherry`DilogExp[f, x] — direct debuggable surface for the engine,
 * bypassing the cascade (tests and benchmarks call it head-to-head with the .m
 * prototype).  Returns the antiderivative or leaves the call unevaluated.
 * Registered by cherry_builtins_init (cherry_driver.c). */
Expr* builtin_cherry_dilog_exp(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    return rt_cherry_dilog_exp(res->data.function.args[0],
                               res->data.function.args[1]);
}
