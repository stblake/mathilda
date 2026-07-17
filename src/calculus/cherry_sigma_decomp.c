/*
 * cherry_sigma_decomp.c — Cherry 1986 Theorem 4.4 restricted Sigma-decomposition
 * (degree-1, all-equal restriction: the logarithmic-integral case).
 *
 * With g(r) = (r,...,r) every term is b_i * P^{r_i}, P = Prod_j f_j, so a
 * decomposition Phi = Sum_i b_i P^{r_i} exists iff Phi is a polynomial in P with
 * constant coefficients (T = Z_{>=0}).  The procedure follows Thm 4.4 exactly:
 *
 *   loop, extracting terms in strictly INCREASING first-exponent:
 *     - a1 = multiplicity of f_1 in Phi_cur;  a1 < 0  => STOP, non-existent
 *       (a pole cannot come from a positive power of P);
 *     - CONSISTENCY (step 7, all-equal so p_j(a1) = a1): every f_j must have
 *       multiplicity a1 in Phi_cur; any mismatch => STOP, non-existent;
 *     - b = (p mod f_1)/(q mod f_1) where Phi_cur / P^{a1} = p/q;  b not in K
 *       => STOP, non-existent;
 *     - Phi_cur <- Phi_cur - b P^{a1};  when it reaches 0 the decomposition is
 *       complete (EXISTS);
 *     - increasing-case degree overshoot (Case ii): a valid decomposition has
 *       at most deg_x(Phi)/deg_x(P) + 1 terms with exponents that cannot exceed
 *       that bound, so if the iteration runs past it with a nonzero residual the
 *       decomposition provably does not exist.
 *
 * The non-existence verdict is the decision property: Integrate[x^2/Log[x^2-1]]
 * declines because (x/2) has no Sigma-decomposition over (x-1, x+1).
 */

#include "cherry_sigma_decomp.h"
#include "cherry_li.h"

#include "risch_util.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "print.h"
#include "expr.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ===================================================================== */
/* Small helpers                                                         */
/* ===================================================================== */

static Expr* mk_int(long n) { return expr_new_integer(n); }

static Expr* mk_pow(Expr* b, long e) {
    return expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ b, mk_int(e) }, 2);
}

/* Multiplicity of the irreducible `f` in the polynomial `poly` (in x): the
 * largest k with f^k | poly, by repeated exact division. */
static long poly_multiplicity(Expr* poly, Expr* f, Expr* x) {
    if (!poly || rt_is_zero(poly)) return 0;
    long k = 0;
    Expr* cur = expr_copy(poly);
    while (k < 100000) {
        Expr* rem = rt_eval3("PolynomialRemainder",
                             expr_copy(cur), expr_copy(f), expr_copy(x));
        bool divides = rem && rt_is_zero(rem);
        if (rem) expr_free(rem);
        if (!divides) break;
        Expr* q = rt_eval3("PolynomialQuotient",
                           expr_copy(cur), expr_copy(f), expr_copy(x));
        expr_free(cur);
        cur = q;
        if (!cur) return k;
        k++;
    }
    expr_free(cur);
    return k;
}

/* Multiplicity of f in a rational Phi = num/den:  mult(num) - mult(den). */
static long phi_multiplicity(Expr* Phi, Expr* f, Expr* x) {
    Expr* tg  = rt_eval1("Together", expr_copy(Phi));
    if (!tg) return 0;
    Expr* num = rt_eval1("Numerator", expr_copy(tg));
    Expr* den = rt_eval1("Denominator", tg);      /* consumes tg */
    long mn = poly_multiplicity(num, f, x);
    long md = poly_multiplicity(den, f, x);
    if (num) expr_free(num);
    if (den) expr_free(den);
    return mn - md;
}

/* ===================================================================== */
/* The decomposition                                                     */
/* ===================================================================== */

SigmaDecomp cherry_sigma_decompose(Expr* Phi, Expr** factors, size_t m, Expr* x) {
    SigmaDecomp d = { SIGMA_UNKNOWN, 0, m, NULL, NULL };
    if (!Phi || !factors || m == 0 || !x) return d;

    /* P = Expand[Prod_j f_j]; its x-degree governs the term bound. */
    Expr* P = mk_int(1);
    for (size_t j = 0; j < m; j++)
        P = rt_eval_own(expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ P, expr_copy(factors[j]) }, 2));
    P = rt_eval_own(expr_new_function(expr_new_symbol(SYM_Expand),
                                      (Expr*[]){ P }, 1));
    long degP = rt_degree(P, x);
    if (degP < 1) { expr_free(P); return d; }        /* degenerate factors */

    Expr* cur = rt_eval1("Cancel", expr_copy(Phi));
    if (!cur) { expr_free(P); return d; }

    /* Term bound: a valid decomposition has exponents in [0, deg(num)/degP]. */
    Expr* num0 = rt_eval1("Numerator", rt_eval1("Together", expr_copy(cur)));
    long degNum = num0 ? rt_degree(num0, x) : 0;
    if (num0) expr_free(num0);
    if (degNum < 0) degNum = 0;
    long iter_cap = degNum / degP + 4;

    Expr** coeffs = NULL; long** exps = NULL;
    size_t n = 0, cap = 0;
    SigmaStatus status = SIGMA_UNKNOWN;
    long prev_a = -1;
    bool first = true;

    for (long it = 0; ; it++) {
        if (rt_is_zero(cur)) { status = SIGMA_EXISTS; break; }
        if (it > iter_cap)   { status = SIGMA_NONEXISTENT; break; }  /* overshoot */

        long a1 = phi_multiplicity(cur, factors[0], x);
        if (a1 < 0)                    { status = SIGMA_NONEXISTENT; break; }
        if (!first && a1 <= prev_a)    { status = SIGMA_NONEXISTENT; break; }

        /* Consistency (Thm 4.4 step 7): every f_j has multiplicity a1. */
        bool consistent = true;
        for (size_t j = 1; j < m; j++) {
            if (phi_multiplicity(cur, factors[j], x) != a1) { consistent = false; break; }
        }
        if (!consistent) { status = SIGMA_NONEXISTENT; break; }

        /* reduced = cur / P^{a1};  b = (p mod f_1) / (q mod f_1). */
        Expr* reduced = (a1 == 0)
            ? rt_eval1("Cancel", expr_copy(cur))
            : rt_eval1("Cancel", expr_new_function(expr_new_symbol(SYM_Times),
                  (Expr*[]){ expr_copy(cur), mk_pow(expr_copy(P), -a1) }, 2));
        Expr* rtg = reduced ? rt_eval1("Together", expr_copy(reduced)) : NULL;
        Expr* rp  = rtg ? rt_eval1("Numerator", expr_copy(rtg)) : NULL;
        Expr* rq  = rtg ? rt_eval1("Denominator", rtg) : NULL;   /* consumes rtg */
        if (reduced) expr_free(reduced);
        Expr* bp = rp ? rt_eval3("PolynomialRemainder",
                                 rp, expr_copy(factors[0]), expr_copy(x)) : NULL;
        Expr* bq = rq ? rt_eval3("PolynomialRemainder",
                                 rq, expr_copy(factors[0]), expr_copy(x)) : NULL;
        Expr* b = (bp && bq) ? rt_eval1("Cancel",
                     expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ bp, mk_pow(bq, -1) }, 2)) : NULL;
        if (bp && !b) expr_free(bp);
        if (bq && !b) expr_free(bq);

        if (!b || !rt_free_of_x(b, x) || rt_is_zero(b)) {
            if (b) expr_free(b);
            status = SIGMA_NONEXISTENT;
            break;
        }

        /* Record (b, [a1, ..., a1]). */
        if (n == cap) {
            cap = cap ? cap * 2 : 8;
            coeffs = realloc(coeffs, cap * sizeof(Expr*));
            exps   = realloc(exps,   cap * sizeof(long*));
        }
        coeffs[n] = expr_copy(b);
        exps[n]   = malloc(m * sizeof(long));
        for (size_t j = 0; j < m; j++) exps[n][j] = a1;
        n++;

        /* cur <- Cancel[cur - b P^{a1}]. */
        Expr* term = (a1 == 0)
            ? expr_copy(b)
            : expr_new_function(expr_new_symbol(SYM_Times),
                  (Expr*[]){ expr_copy(b), mk_pow(expr_copy(P), a1) }, 2);
        Expr* nxt = rt_eval1("Cancel", expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ expr_copy(cur),
                       expr_new_function(expr_new_symbol(SYM_Times),
                           (Expr*[]){ mk_int(-1), term }, 2) }, 2));
        expr_free(b);
        expr_free(cur);
        cur = nxt;
        prev_a = a1;
        first = false;
        if (!cur) { status = SIGMA_UNKNOWN; break; }
    }

    expr_free(cur);
    expr_free(P);

    if (status == SIGMA_EXISTS) {
        d.status = SIGMA_EXISTS; d.n = n; d.coeffs = coeffs; d.exps = exps;
    } else {
        for (size_t i = 0; i < n; i++) { expr_free(coeffs[i]); free(exps[i]); }
        free(coeffs); free(exps);
        d.status = status; d.n = 0; d.coeffs = NULL; d.exps = NULL;
    }
    return d;
}

void cherry_sigma_free(SigmaDecomp* d) {
    if (!d) return;
    for (size_t i = 0; i < d->n; i++) {
        if (d->coeffs) expr_free(d->coeffs[i]);
        if (d->exps)   free(d->exps[i]);
    }
    free(d->coeffs); free(d->exps);
    d->coeffs = NULL; d->exps = NULL; d->n = 0;
}

/* ===================================================================== */
/* Integrate`SigmaDecomposition[Phi, {f1,...,fm}, x]                      */
/* ===================================================================== */

Expr* builtin_sigma_decomposition(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    Expr* Phi  = res->data.function.args[0];
    Expr* flst = res->data.function.args[1];
    Expr* x    = res->data.function.args[2];
    if (!rt_head_is(flst, "List") || flst->data.function.arg_count == 0) return NULL;
    if (x->type != EXPR_SYMBOL) return NULL;

    size_t m = flst->data.function.arg_count;
    Expr** factors = flst->data.function.args;

    SigmaDecomp d = cherry_sigma_decompose(Phi, factors, m, x);

    Expr* out;
    if (d.status == SIGMA_EXISTS) {
        Expr** pairs = malloc(d.n * sizeof(Expr*));
        for (size_t i = 0; i < d.n; i++) {
            Expr** vec = malloc(m * sizeof(Expr*));
            for (size_t j = 0; j < m; j++) vec[j] = mk_int(d.exps[i][j]);
            Expr* vecl = expr_new_function(expr_new_symbol(SYM_List), vec, m);
            free(vec);
            pairs[i] = expr_new_function(expr_new_symbol(SYM_List),
                          (Expr*[]){ expr_copy(d.coeffs[i]), vecl }, 2);
        }
        out = expr_new_function(expr_new_symbol(SYM_List), pairs, d.n);
        free(pairs);
    } else if (d.status == SIGMA_NONEXISTENT) {
        out = expr_new_symbol("$Failed");
    } else {
        out = NULL;   /* out of scope: leave unevaluated */
    }
    cherry_sigma_free(&d);
    return out;
}

/* ===================================================================== */
/* Integrate`LiElementaryQ[f, x] — the logarithmic-integral decision       */
/* ===================================================================== */

/* True  iff INT f is li-elementary (elementary + LogIntegral) — exhibited by
 *          the Cherry li engine producing an antiderivative;
 * False iff the essential form A/Log[w] is PROVEN not li-elementary by the
 *          Thm 4.4 Sigma-decomposition non-existence certificate;
 * unevaluated otherwise (outside the certified single-logarithm scope).
 *
 * Distinct from Risch`ElementaryIntegralQ: Ex 5.1 (x^3/Log[x^2-1]) is
 * li-elementary (True here) but NOT elementary (False there); Ex 5.2
 * (x^2/Log[x^2-1]) is neither (False in both). */
static Expr* builtin_lielementaryq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;

    Expr* anti = rt_cherry_li(f, x);          /* borrowed f */
    if (anti) { expr_free(anti); return expr_new_symbol("True"); }

    if (rt_cherry_li_nonelem(f, x)) return expr_new_symbol("False");

    return NULL;                              /* out of scope: leave unevaluated */
}

void cherry_sigma_init(void) {
    symtab_add_builtin("Integrate`SigmaDecomposition", builtin_sigma_decomposition);
    symtab_get_def("Integrate`SigmaDecomposition")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("Integrate`LiElementaryQ", builtin_lielementaryq);
    symtab_get_def("Integrate`LiElementaryQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Integrate`LiElementaryQ",
        "Integrate`LiElementaryQ[f, x] decides whether f has an antiderivative "
        "that is li-elementary (elementary functions together with the "
        "logarithmic integral LogIntegral): True when the Cherry 1986 engine "
        "exhibits one, False when the Thm 4.4 Sigma-decomposition PROVES none "
        "exists (e.g. x^2/Log[x^2-1]), and unevaluated outside the certified "
        "single-logarithm scope.  Unlike Risch`ElementaryIntegralQ this counts "
        "LogIntegral answers as integrable (x^3/Log[x^2-1] -> True).");
    symtab_set_docstring("Integrate`SigmaDecomposition",
        "Integrate`SigmaDecomposition[Phi, {f1, ..., fm}, x] gives the degree-1 "
        "all-equal restricted Sigma-decomposition of Phi over the irreducibles fi "
        "(Cherry 1986 Thm 4.4): a list of {b_i, {a_i1, ..., a_im}} pairs with "
        "Phi = Sum_i b_i Prod_j fj^(a_ij), or $Failed when Thm 4.4's termination "
        "proves no such decomposition exists.  The engine behind the "
        "logarithmic-integral non-existence decision.");
}
