/* cherry_dilog.c — dilogarithm integration (Cherry degree-2 Sigma-decomposition).
 *
 * Integrands R(x) Log[w] whose antiderivative is a combination of Log-Log products
 * and dilogarithms PolyLog[2, g], solved by matching IN THE LOG TOWER: each Log
 * kernel is an independent variable, and D[PolyLog[2,g]] = -Log[1-g] g'/g is a
 * rational function of the kernels (Mathilda reduces PolyLog[1,g] = -Log[1-g]
 * automatically).  The dilogarithm arguments are the degree-1 Sigma-decomposition
 * candidates g_ij = (x - r_i)/(r_j - r_i) between the RATIONAL ROOTS r_k of the
 * linear factors of {w, den(R)} (g(r_i)=0, g(r_j)=1), so 1-g is again a linear
 * factor.  Setting D[ansatz] == R Log[w] as a polynomial identity in {kernels, x}
 * gives a linear system over the constants, closed by SolveAlways; the exact tower
 * identity is the certificate, re-checked by a PowerExpand diff-back.
 *
 * Generalises the single-term rt_try_dilog recognizer (K Log[1+p x]/x ->
 * -K PolyLog[2,-p x]) to the Log-Log + PolyLog[2] answer form:
 *   INT Log[x]/(1+x) dx = Log[x] Log[1+x] + PolyLog[2,-x],
 *   INT Log[x]/(1-x) dx = PolyLog[2, 1-x].
 *
 * Scope of THIS increment: a single log kernel Log[w] (w a monic linear factor)
 * over a rational R whose denominator splits into monic linear factors, where the
 * degree-1 Sigma-decomposition closes WITHOUT introducing a transcendental
 * constant Log[r_j - r_i] (root spacing +-1) or a complex Log[-1] = i pi.  The
 * transcendental-constant cases (Log[2+x]/x, root spacing != 1) and the multi-log
 * / higher-degree Sigma-decomposition are later increments; they decline cleanly
 * (the unmatched Log leaves the tower system unsolvable).
 */

#include "cherry_dilog.h"
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
static Expr* mk_log(Expr* a) { return expr_new_function(mk_sym("Log"), (Expr*[]){ a }, 1); }

/* Collect the distinct x-dependent Log[w] arguments in e into ws[] (owned copies,
 * deduped), up to cap.  Returns the count. */
static void collect_logs(Expr* e, Expr* x, Expr** ws, size_t* n, size_t cap) {
    if (!e || e->type != EXPR_FUNCTION) return;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == intern_symbol("Log")
        && e->data.function.arg_count == 1) {
        Expr* w = e->data.function.args[0];
        if (!rt_free_of_x(w, x)) {
            for (size_t i = 0; i < *n; i++) if (expr_eq(ws[i], w)) return;
            if (*n < cap) ws[(*n)++] = expr_copy(w);
            return;
        }
    }
    collect_logs(e->data.function.head, x, ws, n, cap);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        collect_logs(e->data.function.args[i], x, ws, n, cap);
}

/* Rational roots of a polynomial p in x (via Solve), appended to rs[] deduped. */
static void collect_rational_roots(Expr* p, Expr* x, Expr** rs, size_t* n, size_t cap) {
    if (!p || rt_degree(p, x) < 1) return;
    Expr* sols = rt_eval2("Solve",
        expr_new_function(mk_sym("Equal"), (Expr*[]){ expr_copy(p), mk_int(0) }, 2),
        expr_copy(x));
    if (sols && sols->type == EXPR_FUNCTION && rt_head_is(sols, "List")) {
        for (size_t i = 0; i < sols->data.function.arg_count && *n < cap; i++) {
            Expr* rule = sols->data.function.args[i];
            if (rule->type != EXPR_FUNCTION || !rt_head_is(rule, "List")
                || rule->data.function.arg_count != 1) continue;
            Expr* rr = rule->data.function.args[0];
            if (rr->type != EXPR_FUNCTION || !rt_head_is(rr, "Rule")) continue;
            Expr* val = rr->data.function.args[1];
            /* keep only rational roots */
            bool rat = (val->type == EXPR_INTEGER || val->type == EXPR_BIGINT)
                || (val->type == EXPR_FUNCTION && rt_head_is(val, "Rational"));
            if (!rat) continue;
            bool dup = false;
            for (size_t k = 0; k < *n; k++) if (expr_eq(rs[k], val)) { dup = true; break; }
            if (!dup) rs[(*n)++] = expr_copy(val);
        }
    }
    if (sols) expr_free(sols);
}

Expr* rt_cherry_dilog(Expr* f, Expr* x) {
    /* 1. log kernels present, and F = f with each Log[w_i] -> t_i. */
    Expr* ws[8]; size_t nw = 0;
    collect_logs(f, x, ws, &nw, 8);
    if (nw == 0 || nw > 4) { for (size_t i = 0; i < nw; i++) expr_free(ws[i]); return NULL; }

    Expr** ts = malloc(nw * sizeof(Expr*));
    Expr* rules[8];
    for (size_t i = 0; i < nw; i++) {
        char nm[24]; snprintf(nm, sizeof(nm), "chdl$t%zu", i);
        ts[i] = mk_sym(nm);
        rules[i] = expr_new_function(mk_sym("Rule"),
            (Expr*[]){ mk_log(expr_copy(ws[i])), expr_copy(ts[i]) }, 2);
    }
    Expr* rl = expr_new_function(mk_sym("List"), rules, nw);
    Expr* F = rt_eval_own(expr_new_function(mk_sym("ReplaceAll"),
        (Expr*[]){ expr_copy(f), rl }, 2));
    /* require F rational in x and LINEAR (degree <= 1) in every log kernel. */
    bool ok = F && rt_free_of_head(F, "Log");
    for (size_t i = 0; ok && i < nw; i++) ok = rt_is_poly(F, ts[i]) ? (rt_degree(F, ts[i]) <= 1)
        : (rt_degree(rt_eval1("Numerator", rt_eval1("Together", expr_copy(F))), ts[i]) <= 1);
    if (F) expr_free(F);
    if (!ok) { for (size_t i = 0; i < nw; i++) { expr_free(ws[i]); expr_free(ts[i]); }
               free(ts); return NULL; }

    /* 2. rational roots of {w_i, denominator of f}. */
    Expr* rs[16]; size_t nr = 0;
    for (size_t i = 0; i < nw; i++)
        if (rt_is_poly(ws[i], x)) collect_rational_roots(ws[i], x, rs, &nr, 16);
    Expr* den = rt_eval1("Denominator", rt_eval1("Together", expr_copy(f)));
    if (den) { collect_rational_roots(den, x, rs, &nr, 16); expr_free(den); }
    if (nr < 2) { for (size_t i = 0; i < nw; i++) { expr_free(ws[i]); expr_free(ts[i]); }
                  free(ts); for (size_t i = 0; i < nr; i++) expr_free(rs[i]); return NULL; }

    /* 3. kernels: Log[x - r_k] -> u_k  (superset of the ws logs, deduped by value). */
    Expr** facs = malloc(nr * sizeof(Expr*));   /* x - r_k */
    Expr** us   = malloc(nr * sizeof(Expr*));   /* kernel vars */
    for (size_t k = 0; k < nr; k++) {
        facs[k] = rt_eval_own(mk_plus2(expr_copy(x), mk_neg(expr_copy(rs[k]))));  /* x - r_k */
        char nm[24]; snprintf(nm, sizeof(nm), "chdl$u%zu", k);
        us[k] = mk_sym(nm);
    }

    /* 4. candidates + unknown coefficients:
     *    - Log-Log products u_i u_j  (i <= j),
     *    - dilogs PolyLog[2, (x-r_i)/(r_j-r_i)]  (ordered pairs i != j),
     *    - single logs u_k. */
    size_t nll = nr * (nr + 1) / 2, ndl = nr * (nr - 1), nsl = nr;
    size_t nc = nll + ndl + nsl;
    Expr** basis  = malloc(nc * sizeof(Expr*));   /* the answer term (in u_k / PolyLog) */
    Expr** dbasis = malloc(nc * sizeof(Expr*));   /* its derivative, logs -> u_k */
    Expr** syms   = malloc(nc * sizeof(Expr*));
    size_t ci = 0;

    /* the Log[x-r_k] -> u_k substitution used to rationalise derivatives */
    Expr** krules = malloc(nr * sizeof(Expr*));
    for (size_t k = 0; k < nr; k++)
        krules[k] = expr_new_function(mk_sym("Rule"),
            (Expr*[]){ mk_log(expr_copy(facs[k])), expr_copy(us[k]) }, 2);
    Expr* klist = expr_new_function(mk_sym("List"), krules, nr);

    /* A candidate is admissible only if its derivative RATIONALISES in the tower —
     * i.e. after Log[x-r_k] -> u_k no Log remains.  This drops the polluting dilog
     * of a reversed root pair (D[PolyLog[2,1+x]] = -Log[-x]/(1+x) leaves Log[-x],
     * an i pi shift) that would otherwise break the linear system. */
    #define ADD_TERM(answer_term)                                                   \
        do {                                                                        \
            char sn[24]; snprintf(sn, sizeof(sn), "chdl$c%zu", ci);                 \
            Expr* s = mk_sym(sn);                                                   \
            Expr* term = mk_times2(expr_copy(s), (answer_term));                    \
            Expr* d = rt_eval2("D", expr_copy(term), expr_copy(x));                 \
            Expr* dr = rt_eval_own(expr_new_function(mk_sym("ReplaceAll"),          \
                (Expr*[]){ d, expr_copy(klist) }, 2));                              \
            if (dr && rt_free_of_head(dr, "Log")) {                                 \
                syms[ci] = s; basis[ci] = term; dbasis[ci] = dr; ci++;              \
            } else {                                                                \
                expr_free(s); expr_free(term); if (dr) expr_free(dr);               \
            }                                                                       \
        } while (0)

    for (size_t i = 0; i < nr; i++)
        for (size_t j = i; j < nr; j++)
            ADD_TERM(mk_times2(mk_log(expr_copy(facs[i])), mk_log(expr_copy(facs[j]))));
    for (size_t i = 0; i < nr; i++)
        for (size_t j = 0; j < nr; j++) {
            if (i == j) continue;
            Expr* g = rt_eval1("Together", mk_times2(
                mk_plus2(expr_copy(x), mk_neg(expr_copy(rs[i]))),
                mk_pow(mk_plus2(expr_copy(rs[j]), mk_neg(expr_copy(rs[i]))), mk_int(-1))));
            ADD_TERM(expr_new_function(mk_sym("PolyLog"),
                (Expr*[]){ mk_int(2), g }, 2));
        }
    for (size_t k = 0; k < nr; k++)
        ADD_TERM(mk_log(expr_copy(facs[k])));
    #undef ADD_TERM
    expr_free(klist);
    nc = ci;                          /* actual number of admissible candidates */

    /* 5. match: F (logs -> u_k) == sum dbasis, as a polynomial identity in {u_k, x}.
     * Build f with Log[x-r_k] -> u_k for every kernel factor (PowerExpand first so
     * a Log of a product/quotient splits into the linear-factor logs). */
    Expr* frules[16];
    for (size_t k = 0; k < nr; k++)
        frules[k] = expr_new_function(mk_sym("Rule"),
            (Expr*[]){ mk_log(expr_copy(facs[k])), expr_copy(us[k]) }, 2);
    Expr* flist = expr_new_function(mk_sym("List"), frules, nr);
    Expr* Fk = rt_eval_own(expr_new_function(mk_sym("ReplaceAll"),
        (Expr*[]){ rt_eval1("PowerExpand", expr_copy(f)), flist }, 2));

    Expr** rhs = malloc(nc * sizeof(Expr*));
    for (size_t i = 0; i < nc; i++) rhs[i] = expr_copy(dbasis[i]);
    Expr* rhs_sum = expr_new_function(mk_sym("Plus"), rhs, nc);
    free(rhs);
    Expr* resid = Fk ? mk_plus2(Fk, mk_neg(rhs_sum)) : NULL;
    Expr* rnum = resid ? rt_eval1("Numerator", rt_eval1("Together", resid)) : NULL;

    Expr* sol = NULL;
    if (rnum) {
        Expr** vl = malloc((nr + 1) * sizeof(Expr*));
        for (size_t k = 0; k < nr; k++) vl[k] = expr_copy(us[k]);
        vl[nr] = expr_copy(x);
        Expr* varlist = expr_new_function(mk_sym("List"), vl, nr + 1);
        free(vl);
        sol = rt_eval2("SolveAlways",
            expr_new_function(mk_sym("Equal"), (Expr*[]){ rnum, mk_int(0) }, 2), varlist);
    }

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
        Q = rt_eval_own(expr_new_function(mk_sym("ReplaceAll"),
            (Expr*[]){ Q, expr_copy(rulesol) }, 2));
        if (Q) {
            Expr** zero = malloc(nc * sizeof(Expr*));
            for (size_t i = 0; i < nc; i++)
                zero[i] = expr_new_function(mk_sym("Rule"),
                    (Expr*[]){ expr_copy(syms[i]), mk_int(0) }, 2);
            Expr* zl = expr_new_function(mk_sym("List"), zero, nc);
            free(zero);
            Q = rt_eval_own(expr_new_function(mk_sym("ReplaceAll"), (Expr*[]){ Q, zl }, 2));
        }
        if (Q && rt_free_of_head(Q, "Integrate")) {
            Expr* diff = mk_plus2(rt_eval2("D", expr_copy(Q), expr_copy(x)), mk_neg(expr_copy(f)));
            Expr* chk = rt_eval1("Simplify", rt_eval1("PowerExpand", diff));
            if (chk && chk->type == EXPR_INTEGER && chk->data.integer == 0) result = Q;
            else expr_free(Q);
            if (chk) expr_free(chk);
        } else if (Q) expr_free(Q);
    }
    if (sol) expr_free(sol);

    for (size_t i = 0; i < nc; i++) { expr_free(basis[i]); expr_free(dbasis[i]); expr_free(syms[i]); }
    free(basis); free(dbasis); free(syms);
    for (size_t k = 0; k < nr; k++) { expr_free(facs[k]); expr_free(us[k]); expr_free(rs[k]); }
    free(facs); free(us);
    for (size_t i = 0; i < nw; i++) { expr_free(ws[i]); expr_free(ts[i]); }
    free(ts);
    return result;
}
