/* knowles_erf.c — Knowles' error-function integration of transcendental
 * Liouvillian functions.
 *
 * P. H. Knowles, "Integration of a Class of Transcendental Liouvillian Functions
 * with Error-Functions", Part I (J. Symb. Comput. 13, 1992) and Part II (16, 1993);
 * SYMSAC '86.  Cherry (1983/85) integrates transcendental *elementary* integrands
 * in terms of erf; Knowles extends the integrand class to transcendental
 * *Liouvillian* — the integrand may itself already contain erf/Ei/li — over the K0
 * RT_PRIM tower (KNOWLES_DESIGN.md §2.1).
 *
 * The extended-Liouville form (SSC 1985 / Cherry Lemma 5.1) is
 *     INT g dx = v + Sum_i k_i erf(u_i),
 * v the elementary part and each erf(u_i) an error-function term whose derivative
 * is (2/Sqrt[Pi]) e^(-u_i^2) u_i'.  Differentiating and matching gives a linear
 * system over the constants (Mack/SSC Thm A1) — exactly Cherry's engine, here run
 * over the multi-generator K0 tower.
 *
 * THIS increment's erf-argument generator is the perfect-square gate (Part I
 * Lemma 6.2 in its clean polynomial form): for every exponential monomial e^{w_j}
 * of the tower, if -w_j is a perfect square with root u_j (so e^{w_j} = e^{-u_j^2}
 * is exactly the v_i = e^{-u_i^2} an erf term needs), erf(u_j) is a candidate; the
 * dual +w_j = u_j^2 yields an Erfi candidate.  The erf term's derivative
 * contribution in the tower is k_j (2/Sqrt[Pi]) t_j D_tower[u_j].  Solved by
 * SolveAlways over the tower variables against an undetermined-coefficient
 * elementary part v (a bounded-degree polynomial in the tower monomials with
 * constant coefficients), then diff-back verified (rt_verify_antideriv) — so a
 * mis-generation can only DECLINE, never emit a wrong antiderivative.
 *
 * Pins (KNOWLES_DESIGN.md §1.7, Mathilda's classical Erf convention
 * erf(u) = (Sqrt[Pi]/2) Erf[u]):
 *   INT E^(-x^2 - Erf[x]^2) dx           = (Pi/4) Erf[Erf[x]]         (Part II Ex 4.1)
 *   INT Erf[x] E^(-x^2 - Erf[x]^2) dx    = -(Sqrt[Pi]/4) E^(-Erf[x]^2) (Ex 4.3, elementary)
 *   INT (2 E^(-x^2) Erf[x] - 3 E^(-1/x^2)/x^2) dx = Erf[x]^2 + 3 Erf[1/x]  (Ex 4.4)
 *   INT x E^(-x^2 - Erf[x]^2) dx         -> declines (no erf-elementary antideriv, Ex 4.2)
 *
 * Scope of THIS increment: rational (perfect-square) erf arguments with constant
 * elementary-part coefficients.  The quasiquadratic completing-square case (Part I
 * Ex 8.1, radical arguments like 1/Sqrt[log x]) and x-rational v coefficients are
 * later increments (KNOWLES_DESIGN.md §3, K2 sub-milestones).
 */

#include "knowles_erf.h"
#include "risch_tower.h"
#include "risch_util.h"

#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "print.h"
#include "sym_intern.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static Expr* mk_sym(const char* s)   { return expr_new_symbol(s); }
static Expr* mk_int(long n)          { return expr_new_integer(n); }
static Expr* mk_pow(Expr* b, Expr* e){ return expr_new_function(mk_sym("Power"), (Expr*[]){ b, e }, 2); }
static Expr* mk_times2(Expr* a, Expr* b){ return expr_new_function(mk_sym("Times"), (Expr*[]){ a, b }, 2); }
static Expr* mk_neg(Expr* a)         { return mk_times2(mk_int(-1), a); }
static Expr* mk_plus2(Expr* a, Expr* b){ return expr_new_function(mk_sym("Plus"), (Expr*[]){ a, b }, 2); }

/* 2/Sqrt[Pi] — the Mathilda-Erf derivative constant. */
static Expr* two_over_sqrt_pi(void) {
    return mk_times2(mk_int(2), mk_pow(mk_sym("Pi"),
        expr_new_function(mk_sym("Rational"), (Expr*[]){ mk_int(-1), mk_int(2) }, 2)));
}

static bool is_failed(Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == intern_symbol("$Failed");
}

/* A candidate erf/erfi term: head Erf|Erfi, argument u (in tower variables), and
 * the tower variable expvar = e^{w} = e^{-/+ u^2}. */
typedef struct { const char* head; Expr* u; Expr* expvar; } ErfCand;

/* PolynomialSqrt[target]; return the owned root if target is a perfect square
 * (verified Cancel[root^2 - target] == 0), else NULL.  TAKES OWNERSHIP of target. */
static Expr* perfect_sqrt(Expr* target) {
    if (!target) return NULL;
    Expr* r = rt_eval_call("PolynomialSqrt", (Expr*[]){ expr_copy(target) }, 1);
    if (!r || is_failed(r) || rt_head_is(r, "PolynomialSqrt")) {
        if (r) expr_free(r);
        expr_free(target);
        return NULL;
    }
    Expr* chk = rt_eval1("Cancel", mk_plus2(mk_pow(expr_copy(r), mk_int(2)),
                                            mk_neg(target)));   /* consumes target */
    bool ok = chk && rt_is_zero(chk);
    if (chk) expr_free(chk);
    if (!ok) { expr_free(r); return NULL; }
    return r;
}

Expr* knowles_erf_liouvillian(Expr* f, Expr* x) {
    /* 1. Pre-normalise (split prim-bearing exponentials) and build the K0 tower. */
    Expr* fx = rt_expand_exp_sums(f);
    RtTower T;
    if (!rt_tower_build_min(fx, x, &T, 1)) { rt_tower_free(&T); expr_free(fx); return NULL; }

    /* Erf integration over a tangent sub-tower is out of scope for this increment. */
    for (size_t i = 0; i < T.n; i++)
        if (T.kind[i] == RT_TAN) { rt_tower_free(&T); expr_free(fx); return NULL; }

    /* 2. Integrand in tower variables; must be rational there (fully substituted). */
    Expr* F = rt_subst_kernels(fx, &T);
    if (!F || rt_find_exp_of_x(F, x) || rt_find_log_of_x(F, x)
        || !rt_free_of_head(F, "Erf")  || !rt_free_of_head(F, "Erfi")
        || !rt_free_of_head(F, "Erfc") || !rt_free_of_head(F, "ExpIntegralEi")
        || !rt_free_of_head(F, "LogIntegral")) {
        if (F) expr_free(F);
        rt_tower_free(&T); expr_free(fx); return NULL;
    }

    /* 3. Erf/Erfi candidates from the perfect-square gate over the exp kernels.
     *    For e^{w}: Erf needs w = -u^2 (target -w), Erfi needs w = +u^2 (target w).
     *    Exactly one target is a REAL perfect square for real w — prefer that (so
     *    e^{-x^2} -> Erf[x], e^{x^2} -> Erfi[x], and e^{-erf^2 x} -> Erf[erf x]),
     *    keeping the answer I-free.  A lone I-laden root is still sound (diff-back). */
    ErfCand cand[32]; size_t nc = 0;
    Expr* Isym = parse_expression("I");
    for (size_t j = 0; j < T.n && nc < 30; j++) {
        if (T.kind[j] != RT_EXP) continue;
        Expr* wsub = rt_subst_kernels(T.arg[j], &T);          /* exponent w in tower vars */
        if (!wsub) continue;
        Expr* u_erf  = perfect_sqrt(rt_eval1("Cancel", mk_neg(expr_copy(wsub))));
        Expr* u_erfi = perfect_sqrt(rt_eval1("Cancel", expr_copy(wsub)));
        expr_free(wsub);
        /* choose head/u: prefer the real (I-free) root */
        const char* head = NULL; Expr* u = NULL;
        Expr* fq_erf  = u_erf  ? rt_eval2("FreeQ", expr_copy(u_erf),  expr_copy(Isym)) : NULL;
        Expr* fq_erfi = u_erfi ? rt_eval2("FreeQ", expr_copy(u_erfi), expr_copy(Isym)) : NULL;
        if (u_erf && rt_is_true(fq_erf))       { head = "Erf";  u = u_erf;  }
        else if (u_erfi && rt_is_true(fq_erfi)){ head = "Erfi"; u = u_erfi; }
        else if (u_erf)                        { head = "Erf";  u = u_erf;  }
        else if (u_erfi)                       { head = "Erfi"; u = u_erfi; }
        if (fq_erf)  expr_free(fq_erf);
        if (fq_erfi) expr_free(fq_erfi);
        /* reject a CONSTANT u (free of x AND every tower variable). */
        bool u_const = u && rt_free_of_x(u, x);
        for (size_t g = 0; u_const && g < T.n; g++)
            if (!rt_free_of_x(u, T.t[g])) u_const = false;
        if (u && !u_const) {
            bool dup = false;
            for (size_t c = 0; c < nc; c++)
                if (cand[c].head == intern_symbol(head) && expr_eq(cand[c].u, u)) { dup = true; break; }
            if (!dup) { cand[nc].head = intern_symbol(head);
                        cand[nc].u = expr_copy(u); cand[nc].expvar = expr_copy(T.t[j]); nc++; }
        }
        if (u_erf)  expr_free(u_erf);
        if (u_erfi) expr_free(u_erfi);
    }
    if (Isym) expr_free(Isym);
    if (nc == 0) {                                            /* no erf term possible here */
        expr_free(F); rt_tower_free(&T); expr_free(fx);
        for (size_t c = 0; c < nc; c++) { expr_free(cand[c].u); expr_free(cand[c].expvar); }
        return NULL;
    }

    /* 4. Elementary-part basis: tower monomials Prod t_i^{d_i} with per-kind degree
     *    bounds and CONSTANT coefficients (this increment).  Bounded and capped. */
    long dlo[24], dhi[24];
    if (T.n > 24) { /* absurdly deep tower — decline */
        expr_free(F); for (size_t c=0;c<nc;c++){expr_free(cand[c].u);expr_free(cand[c].expvar);}
        rt_tower_free(&T); expr_free(fx); return NULL;
    }
    size_t total = 1;
    for (size_t i = 0; i < T.n; i++) {
        if (T.kind[i] == RT_EXP)       { dlo[i] = -1; dhi[i] = 1; }
        else if (T.kind[i] == RT_PRIM) { dlo[i] =  0; dhi[i] = 2; }
        else                           { dlo[i] =  0; dhi[i] = 1; } /* RT_LOG */
        total *= (size_t)(dhi[i] - dlo[i] + 1);
        if (total > 512) {                                   /* keep SolveAlways tractable */
            expr_free(F); for (size_t c=0;c<nc;c++){expr_free(cand[c].u);expr_free(cand[c].expvar);}
            rt_tower_free(&T); expr_free(fx); return NULL;
        }
    }

    size_t nsym = total + nc;
    Expr** syms = malloc(nsym * sizeof(Expr*));
    size_t si = 0;

    /* v = Sum_b c_b * monomial_b ; also keep the monomials for the emitted answer. */
    Expr** vterms = malloc(total * sizeof(Expr*));
    Expr** monos  = malloc(total * sizeof(Expr*));
    long decode[24];
    for (size_t b = 0; b < total; b++) {
        size_t r = b;
        for (size_t i = 0; i < T.n; i++) {
            size_t sz = (size_t)(dhi[i] - dlo[i] + 1);
            decode[i] = dlo[i] + (long)(r % sz); r /= sz;
        }
        /* build the monomial Prod t_i^{decode_i} (skip exponent 0) */
        Expr* mono = mk_int(1);
        for (size_t i = 0; i < T.n; i++) {
            if (decode[i] == 0) continue;
            Expr* fac = mk_pow(expr_copy(T.t[i]), mk_int(decode[i]));
            mono = mk_times2(mono, fac);
        }
        char nm[32]; snprintf(nm, sizeof(nm), "kerf$c%zu", b);
        Expr* c = mk_sym(nm);
        syms[si++] = expr_copy(c);
        monos[b]  = mono;                                    /* owned */
        vterms[b] = mk_times2(c, expr_copy(mono));
    }
    Expr* v_ansatz = expr_new_function(mk_sym("Plus"), vterms, total);  /* consumes vterms */
    free(vterms);

    /* 5. D_tower[v] + Sum_j k_j (2/Sqrt[Pi]) t_j D_tower[u_j]. */
    Expr** kd = malloc(nc * sizeof(Expr*));                  /* the k_j symbols */
    Expr** dterms = malloc((nc + 1) * sizeof(Expr*)); size_t nd = 0;
    dterms[nd++] = rt_tower_deriv(v_ansatz, &T, x);
    for (size_t c = 0; c < nc; c++) {
        char nm[32]; snprintf(nm, sizeof(nm), "kerf$k%zu", c);
        Expr* k = mk_sym(nm);
        kd[c] = expr_copy(k);
        syms[si++] = expr_copy(k);
        Expr* du = rt_tower_deriv(cand[c].u, &T, x);
        dterms[nd++] = mk_times2(k, mk_times2(two_over_sqrt_pi(),
            mk_times2(expr_copy(cand[c].expvar), du)));
    }
    Expr* Dexpr = expr_new_function(mk_sym("Plus"), dterms, nd);
    free(dterms);

    /* 6. residual numerator == 0 identically in the tower variables + x. */
    Expr* resid = mk_plus2(expr_copy(F), mk_neg(Dexpr));
    Expr* rnum  = rt_eval1("Numerator", rt_eval1("Together", resid));
    Expr* sol = NULL;
    if (rnum) {
        Expr** vl = malloc((T.n + 1) * sizeof(Expr*));
        for (size_t i = 0; i < T.n; i++) vl[i] = expr_copy(T.t[T.n - 1 - i]);
        vl[T.n] = expr_copy(x);
        Expr* varlist = expr_new_function(mk_sym("List"), vl, T.n + 1);
        free(vl);
        Expr* eqn = expr_new_function(mk_sym("Equal"), (Expr*[]){ rnum, mk_int(0) }, 2);
        sol = rt_eval2("SolveAlways", eqn, varlist);
    }

    /* 7. Assemble the answer, back-substitute tower vars -> kernels, diff-back verify. */
    Expr* result = NULL;
    bool solved = sol && sol->type == EXPR_FUNCTION && rt_head_is(sol, "List")
        && sol->data.function.arg_count >= 1
        && sol->data.function.args[0]->type == EXPR_FUNCTION
        && rt_head_is(sol->data.function.args[0], "List");
    if (solved) {
        Expr* rules = sol->data.function.args[0];
        size_t nans = total + nc;
        Expr** ans = malloc(nans * sizeof(Expr*)); size_t na = 0;
        for (size_t b = 0; b < total; b++) {
            char nm[32]; snprintf(nm, sizeof(nm), "kerf$c%zu", b);
            ans[na++] = mk_times2(mk_sym(nm), expr_copy(monos[b]));
        }
        for (size_t c = 0; c < nc; c++)
            ans[na++] = mk_times2(expr_copy(kd[c]),
                expr_new_function(mk_sym(cand[c].head == intern_symbol("Erf") ? "Erf" : "Erfi"),
                    (Expr*[]){ expr_copy(cand[c].u) }, 1));
        Expr* Q = expr_new_function(mk_sym("Plus"), ans, nans);
        free(ans);
        Q = rt_eval_own(expr_new_function(mk_sym("ReplaceAll"),
            (Expr*[]){ Q, expr_copy(rules) }, 2));            /* substitute solution */
        if (Q) {                                             /* pin free unknowns to 0 */
            Expr** zero = malloc(nsym * sizeof(Expr*));
            for (size_t jj = 0; jj < nsym; jj++)
                zero[jj] = expr_new_function(mk_sym("Rule"),
                    (Expr*[]){ expr_copy(syms[jj]), mk_int(0) }, 2);
            Expr* zl = expr_new_function(mk_sym("List"), zero, nsym);
            free(zero);
            Q = rt_eval_own(expr_new_function(mk_sym("ReplaceAll"), (Expr*[]){ Q, zl }, 2));
        }
        if (Q) {                                             /* tower vars -> kernels */
            Expr** br = malloc((T.n ? T.n : 1) * sizeof(Expr*));
            for (size_t i = 0; i < T.n; i++)
                br[i] = expr_new_function(mk_sym("Rule"),
                    (Expr*[]){ expr_copy(T.t[i]), expr_copy(T.kernel[i]) }, 2);
            Expr* bl = expr_new_function(mk_sym("List"), br, T.n);
            free(br);
            Q = rt_eval_own(expr_new_function(mk_sym("ReplaceAll"), (Expr*[]){ Q, bl }, 2));
        }
        if (Q && rt_free_of_head(Q, "Integrate") && rt_verify_antideriv(Q, f, x))
            result = Q;
        else if (Q) expr_free(Q);
    }
    if (sol) expr_free(sol);

    /* 8. cleanup */
    for (size_t b = 0; b < total; b++) expr_free(monos[b]);
    free(monos);
    for (size_t c = 0; c < nc; c++) { expr_free(cand[c].u); expr_free(cand[c].expvar); expr_free(kd[c]); }
    free(kd);
    for (size_t jj = 0; jj < nsym; jj++) expr_free(syms[jj]);
    free(syms);
    expr_free(v_ansatz); expr_free(F);
    rt_tower_free(&T); expr_free(fx);
    return result;
}
