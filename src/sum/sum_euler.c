/*
 * sum_euler.c -- Sum`Euler: closed forms for infinite linear Euler sums.
 *
 * A linear Euler sum is
 *
 *     sigma_h(p, q) = Sum_{k=1}^Infinity  H_k^{(p)} / k^q ,
 *
 * where H_k^{(p)} = Sum_{j=1}^k 1/j^p is the generalized harmonic number
 * (H_k^{(1)} = H_k).  It converges for q >= 2.  This stage recognises a summand
 * that is (a constant times) a single HarmonicNumber factor over k^q and, when a
 * closed form in Riemann zeta values is known, returns it -- otherwise it leaves
 * the sum unevaluated (returns NULL), never a fabricated value.
 *
 * Closed forms implemented (all rigorous, no numeric fitting):
 *
 *   p = 1  (Euler 1775), for q >= 2:
 *       Sum H_k / k^q = (1 + q/2) zeta(q+1)
 *                       - (1/2) Sum_{j=1}^{q-2} zeta(j+1) zeta(q-j).
 *     e.g. q = 2 -> 2 zeta(3);  q = 3 -> Pi^4/72.
 *
 *   p = q  (diagonal), from the reflection sigma_h(p,q) + sigma_h(q,p)
 *          = zeta(p) zeta(q) + zeta(p+q):
 *       Sum H_k^{(p)} / k^p = (zeta(p)^2 + zeta(2p)) / 2.
 *     e.g. p = 2 -> 7 Pi^4/360 (= 7/4 zeta(4)).
 *
 * Even zeta values reduce to powers of Pi automatically via the Zeta builtin, so
 * the results come out elementary.  General non-diagonal sigma_h(p, q) with
 * p, q >= 2 (odd weight) is not yet reduced and stays unevaluated.
 */

#include "sum_internal.h"
#include "sum_euler_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "sym_names.h"
#include "arithmetic.h"   /* is_infinity_sym */
#include <stdlib.h>
#include <stdio.h>

/* Zeta[n] as an unevaluated node. */
Expr* eu_zeta(int n) {
    return expr_new_function(expr_new_symbol(SYM_Zeta), (Expr*[]){ sum_int(n) }, 1);
}

/* 1/2 as Power[2, -1]. */
Expr* eu_half(void) {
    return expr_new_function(expr_new_symbol(SYM_Power),
               (Expr*[]){ sum_int(2), sum_int(-1) }, 2);
}

/* Euler's formula for Sum_{k>=1} H_k / k^q  (q >= 2), returned unevaluated. */
Expr* euler_order1(int q) {
    /* term1 = (q+2)/2 * Zeta[q+1] */
    Expr* coeff = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ sum_int(q + 2), eu_half() }, 2);
    Expr* acc = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ coeff, eu_zeta(q + 1) }, 2);
    /* - (1/2) Sum_{j=1}^{q-2} Zeta[j+1] Zeta[q-j] */
    for (int j = 1; j <= q - 2; j++) {
        Expr* prod = expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){ eu_zeta(j + 1), eu_zeta(q - j) }, 2);
        Expr* neg = expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){ sum_int(-1), eu_half(), prod }, 3);
        acc = expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ acc, neg }, 2);
    }
    return acc;
}

/* Diagonal Sum_{k>=1} H_k^{(p)} / k^p = (Zeta[p]^2 + Zeta[2p]) / 2, unevaluated. */
static Expr* euler_diag(int p) {
    Expr* zp2 = expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ eu_zeta(p), sum_int(2) }, 2);
    Expr* s = expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ zp2, eu_zeta(2 * p) }, 2);
    return expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ s, eu_half() }, 2);
}

/* --- helpers for the non-diagonal / quadratic reductions --- */

static bool is_harmonic_factor(Expr* g, Expr* var, int* p_out);   /* defined below */

Expr* eu_times2(Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ a, b }, 2);
}
Expr* eu_plus2(Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){ a, b }, 2);
}
/* Exact integer binomial C(n, k) (n, k small: weights <= ~15). */
int64_t eu_binom(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    int64_t r = 1;
    for (int i = 1; i <= k; i++) r = r * (n - k + i) / i;
    return r;
}

/* The double zeta value Z(s,t) = Sum_{m>n>=1} m^{-s} n^{-t} for ODD s, in
 * Riemann-zeta values (Borwein-Borwein-Girgensohn), returned unevaluated:
 *   Z(s,t) = 1/2 ((-1)^s C(w,s) - 1) Zeta[w]
 *          + Sum_{k=1}^{(w-3)/2} [C(w-2k-1,s-1)+C(w-2k-1,t-1)] Zeta[2k] Zeta[w-2k].
 */
Expr* euler_Z_odd(int s, int t) {
    int w = s + t;
    int64_t sign = (s % 2 == 0) ? 1 : -1;
    int64_t Anum = sign * eu_binom(w, s) - 1;         /* A = Anum/2 */
    Expr* acc = eu_times2(make_rational(Anum, 2), eu_zeta(w));
    for (int k = 1; k <= (w - 3) / 2; k++) {
        int O = w - 2 * k;                            /* odd >= 3 */
        int64_t coef = eu_binom(w - 2 * k - 1, s - 1) + eu_binom(w - 2 * k - 1, t - 1);
        Expr* term = eu_times2(sum_int(coef), eu_times2(eu_zeta(2 * k), eu_zeta(O)));
        acc = eu_plus2(acc, term);
    }
    return acc;
}

/* Non-diagonal linear Euler sum sigma(p,q) = Sum H_k^{(p)}/k^q, p != q, p,q >= 2.
 * Uses sigma(p,q) = Z(q,p) + Zeta[w].  Reducible to zeta values only for ODD
 * weight w = p+q; returns NULL (stays unevaluated) for even weight. */
static Expr* euler_nondiag(int p, int q) {
    int w = p + q;
    if (w % 2 == 0) return NULL;
    if (q % 2 == 1) {                                 /* outer q odd: direct */
        return eu_plus2(euler_Z_odd(q, p), eu_zeta(w));
    }
    /* q even => p odd.  Reflection sigma(p,q) = Zeta[p]Zeta[q] + Zeta[w] - sigma(q,p),
     * with sigma(q,p) = Z(p,q) + Zeta[w]  (Z first index p is odd). */
    Expr* sigma_qp = eu_plus2(euler_Z_odd(p, q), eu_zeta(w));
    Expr* refl = eu_plus2(eu_times2(eu_zeta(p), eu_zeta(q)), eu_zeta(w));
    return eu_plus2(refl, eu_times2(sum_int(-1), sigma_qp));
}

/* Quadratic Euler sum Sum_{k>=1} H_k^2 / k^q in zeta values.  Rigorous per-weight
 * table for q = 2..5 (weight w = q+2 <= 7); an irreducible MZV appears at w = 8,
 * so q >= 6 returns NULL. */
static Expr* euler_quadratic(int q) {
    switch (q) {
    case 2:  /* 17/4 Zeta[4] */
        return eu_times2(make_rational(17, 4), eu_zeta(4));
    case 3:  /* 7/2 Zeta[5] - Zeta[2] Zeta[3] */
        return eu_plus2(eu_times2(make_rational(7, 2), eu_zeta(5)),
                        eu_times2(sum_int(-1), eu_times2(eu_zeta(2), eu_zeta(3))));
    case 4:  /* 97/24 Zeta[6] - 2 Zeta[3]^2 */
        return eu_plus2(eu_times2(make_rational(97, 24), eu_zeta(6)),
                        eu_times2(sum_int(-2),
                            expr_new_function(expr_new_symbol(SYM_Power),
                                (Expr*[]){ eu_zeta(3), sum_int(2) }, 2)));
    case 5:  /* 6 Zeta[7] - Zeta[2] Zeta[5] - 5/2 Zeta[3] Zeta[4] */
        return eu_plus2(
                   eu_plus2(eu_times2(sum_int(6), eu_zeta(7)),
                            eu_times2(sum_int(-1), eu_times2(eu_zeta(2), eu_zeta(5)))),
                   eu_times2(make_rational(-5, 2), eu_times2(eu_zeta(3), eu_zeta(4))));
    default: return NULL;
    }
}

/* g == Power[HarmonicNumber[var], 2]?  (H_k squared, order-1 harmonic only.) */
static bool is_harmonic_sq(Expr* g, Expr* var) {
    if (g->type != EXPR_FUNCTION || g->data.function.head->type != EXPR_SYMBOL
        || g->data.function.head->data.symbol.name != SYM_Power
        || g->data.function.arg_count != 2) return false;
    Expr* e2 = g->data.function.args[1];
    if (e2->type != EXPR_INTEGER || e2->data.integer != 2) return false;
    int pp;
    return is_harmonic_factor(g->data.function.args[0], var, &pp) && pp == 1;
}

/* g == HarmonicNumber[var] or HarmonicNumber[var, p] with p a positive integer?
 * Sets *p_out to the order (1 for the one-argument form). */
static bool is_harmonic_factor(Expr* g, Expr* var, int* p_out) {
    if (g->type != EXPR_FUNCTION) return false;
    Expr* h = g->data.function.head;
    if (h->type != EXPR_SYMBOL || h->data.symbol.name != SYM_HarmonicNumber) return false;
    size_t ac = g->data.function.arg_count;
    if (ac < 1 || ac > 2) return false;
    if (!expr_eq(g->data.function.args[0], var)) return false;
    if (ac == 1) { *p_out = 1; return true; }
    Expr* pa = g->data.function.args[1];
    if (pa->type != EXPR_INTEGER || pa->data.integer < 1) return false;
    *p_out = (int)pa->data.integer;
    return true;
}

Expr* builtin_sum_euler(Expr* res);

Expr* builtin_sum_euler(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!sum_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (!definite) return NULL;
    if (!is_infinity_sym(imax)) return NULL;
    /* H_k is defined for k >= 1; require the standard lower bound. */
    if (imin->type != EXPR_INTEGER || imin->data.integer != 1) return NULL;

    /* Canonicalise the held summand (1/k^q parses as (k^q)^(-1)). */
    Expr* fc = expr_copy(f);
    Expr* fn = evaluate(fc);
    expr_free(fc);   /* evaluate() borrows; free the copy we made */

    bool is_times = (fn->type == EXPR_FUNCTION
                     && fn->data.function.head->type == EXPR_SYMBOL
                     && fn->data.function.head->data.symbol.name == SYM_Times);
    size_t n = is_times ? fn->data.function.arg_count : 1;

    int p = 0, hcount = 0, hsqcount = 0;
    Expr* rest = sum_int(1);     /* product of the non-harmonic factors */
    for (size_t i = 0; i < n; i++) {
        Expr* g = is_times ? fn->data.function.args[i] : fn;
        int pp;
        if (is_harmonic_factor(g, var, &pp)) {
            hcount++; p = pp;
        } else if (is_harmonic_sq(g, var)) {
            hsqcount++;                              /* Power[HarmonicNumber[k], 2] */
        } else {
            Expr* t = expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ rest, expr_copy(g) }, 2);
            rest = evaluate(t);
            expr_free(t);
        }
    }
    expr_free(fn);
    /* Exactly one linear H factor (linear/non-diagonal sums), OR exactly one
     * squared H_k (quadratic sums).  Anything else is out of scope. */
    if (!((hcount == 1 && hsqcount == 0) || (hcount == 0 && hsqcount == 1))) {
        expr_free(rest); return NULL;
    }

    /* rest must be c * var^(-q) with c free of var and q >= 2 (convergence). */
    Expr* one = sum_int(1);
    Expr* c = sum_subst(rest, var, one);       /* c = rest|_{var->1} */
    expr_free(one);
    if (!sum_free_of(c, var)) { expr_free(rest); expr_free(c); return NULL; }
    Expr* cinv = expr_new_function(expr_new_symbol(SYM_Power),
                     (Expr*[]){ expr_copy(c), sum_int(-1) }, 2);
    Expr* mt = expr_new_function(expr_new_symbol(SYM_Times),
                   (Expr*[]){ rest, cinv }, 2);      /* consumes rest, cinv */
    Expr* mono = evaluate(mt);
    expr_free(mt);

    int q = 0; bool okmono = false;
    if (mono->type == EXPR_FUNCTION && mono->data.function.head->type == EXPR_SYMBOL
        && mono->data.function.head->data.symbol.name == SYM_Power
        && mono->data.function.arg_count == 2
        && expr_eq(mono->data.function.args[0], var)
        && mono->data.function.args[1]->type == EXPR_INTEGER) {
        int64_t e = mono->data.function.args[1]->data.integer;
        if (e < 0) { q = (int)(-e); okmono = true; }
    }
    expr_free(mono);
    if (!okmono || q < 2) { expr_free(c); return NULL; }

    /* Select the closed form.  Even-weight non-diagonal and q >= 6 quadratic
     * have no zeta reduction -> euler_nondiag/euler_quadratic return NULL. */
    Expr* sig;
    if (hsqcount == 1)   sig = euler_quadratic(q);
    else if (p == 1)     sig = euler_order1(q);
    else if (p == q)     sig = euler_diag(p);
    else                 sig = euler_nondiag(p, q);
    if (!sig) { expr_free(c); return NULL; }

    Expr* raw = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ c, sig }, 2);   /* consumes c, sig */
    Expr* out = evaluate(raw);
    expr_free(raw);
    return out;
}

void sum_euler_init(void) {
    symtab_add_builtin("Sum`Euler", builtin_sum_euler);
    symtab_get_def("Sum`Euler")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Sum`Euler",
        "Sum`Euler[f, k, 1, Infinity] gives the closed form of an infinite Euler "
        "sum: a constant times HarmonicNumber[k], HarmonicNumber[k, p], or "
        "HarmonicNumber[k]^2 over k^q (q >= 2), in Riemann zeta values. Handles "
        "order 1 (Euler's formula), the diagonal p = q, non-diagonal linear sums "
        "of odd weight (double-zeta reduction), and quadratic sums H_k^2/k^q for "
        "q <= 5. Returns unevaluated when no zeta reduction exists.");
}
