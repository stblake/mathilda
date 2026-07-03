/*
 * product_rational_infinite.c -- Product`RationalInfinite: convergent infinite
 * rational products with complex-conjugate roots, via the Gamma canonical form.
 *
 * A convergent product with unit step,
 *   prod_{k=a}^Inf  c * prod_i (k-alpha_i)^{m_i} / prod_j (k-beta_j)^{n_j},
 * converges to a nonzero limit iff c==1, the degrees balance, and the sums of
 * roots balance (the 1/k term vanishes).  Then
 *   prod_{k=a}^Inf R(k) = prod_j Gamma(a-beta_j)^{n_j} / prod_i Gamma(a-alpha_i)^{m_i}.
 *
 * The resulting product of Gammas at rational / complex-conjugate arguments is
 * reduced to elementary closed form in C:
 *   - integer shift  Gamma(z+1) = z Gamma(z)  to bring Re(z) into (0,1];
 *   - cancellation of equal arguments (opposite exponent);
 *   - conjugate pairs  Gamma(1+ib)Gamma(1-ib) = Pi b / Sinh(Pi b),
 *                      Gamma(1/2+ib)Gamma(1/2-ib) = Pi / Cosh(Pi b);
 *   - real reflection  Gamma(r)Gamma(1-r) = Pi / Sin(Pi r).
 * A residual complex Gamma that does not collapse -> bail (NULL).  The candidate
 * is numerically checked against the raw Gamma product before it is returned.
 *
 * Examples:  prod_{k>=2} (k^2-1)/(k^2+1) = Pi/Sinh[Pi];
 *            prod_{k>=2} (k^3-1)/(k^3+1) = 2/3.
 *
 * This stage runs only when a genuine COMPLEX root is present; all-real-root
 * products are left to the (already working) Product`Infinite limit route.
 *
 * Memory contract: takes ownership of res but must not free it; returns an owned
 * closed form or NULL to fall through.
 */

#include "product_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "arithmetic.h"   /* arith_warnings_muted */
#include "sym_names.h"
#include "poly.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

static bool is_inf_sym(const Expr* e) {
    return e->type == EXPR_SYMBOL && e->data.symbol == SYM_Infinity;
}

static void warn_div(void) {
    if (!arith_warnings_muted())
        fprintf(stderr, "Product::div: Product does not converge.\n");
}

/* evaluate(Head[a]) for a single copied arg. */
static Expr* fn1(const char* head, Expr* a) {
    Expr* e = expr_new_function(expr_new_symbol(head), (Expr*[]){ expr_copy(a) }, 1);
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

/* PossibleZeroQ[e] == True (consumes e). */
static bool zeroq(Expr* e) {
    Expr* r = expr_new_function(expr_new_symbol("PossibleZeroQ"), (Expr*[]){ e }, 1);
    Expr* v = evaluate(r);
    expr_free(r);
    bool yes = (v && v->type == EXPR_SYMBOL && v->data.symbol == SYM_True);
    if (v) expr_free(v);
    return yes;
}

/* evaluate boolean cmp(a,b) and test for literal True.  a is copied; b is
 * ADOPTED (freed here) so callers may pass a fresh temporary. */
static bool cmp_true(const char* head, Expr* a, Expr* b) {
    Expr* e = expr_new_function(expr_new_symbol(head),
                  (Expr*[]){ expr_copy(a), b }, 2);   /* adopts b */
    Expr* v = evaluate(e);
    expr_free(e);
    bool yes = (v && v->type == EXPR_SYMBOL && v->data.symbol == SYM_True);
    if (v) expr_free(v);
    return yes;
}

/* pref <- evaluate(Times[pref, factor^exp]); factor copied, pref adopted. */
static Expr* mul_pow(Expr* pref, Expr* factor, int exp) {
    Expr* p = expr_new_function(expr_new_symbol(SYM_Power),
                  (Expr*[]){ expr_copy(factor), expr_new_integer(exp) }, 2);
    Expr* t = expr_new_function(expr_new_symbol(SYM_Times),
                  (Expr*[]){ pref, p }, 2);   /* adopts pref, p */
    Expr* r = evaluate(t);
    expr_free(t);
    return r;
}

typedef struct { Expr* arg; int exp; } GT;

static void free_roots(Expr** r, int* m, size_t n) {
    if (r) { for (size_t i = 0; i < n; i++) expr_free(r[i]); free(r); }
    free(m);
}

Expr* builtin_product_rational_infinite(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!product_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (!definite || !is_inf_sym(imax)) return NULL;
    if (prod_has_symbolic_power(f, var)) return NULL;

    /* Rational split. */
    Expr* tog = prod_eval("Together",    (Expr*[]){ expr_copy(f) }, 1);
    Expr* num = prod_eval("Numerator",   (Expr*[]){ expr_copy(tog) }, 1);
    Expr* den = prod_eval("Denominator", (Expr*[]){ tog }, 1);   /* adopts tog */

    Expr *leadN = NULL, *leadD = NULL, **rootsN = NULL, **rootsD = NULL;
    int *multsN = NULL, *multsD = NULL;
    size_t nN = 0, nD = 0;
    bool okN = false, okD = false;
    bool gotN = prod_rational_roots(num, var, &leadN, &rootsN, &multsN, &nN, &okN);
    bool gotD = prod_rational_roots(den, var, &leadD, &rootsD, &multsD, &nD, &okD);
    expr_free(num); expr_free(den);
    if (!gotN || !gotD || !okN || !okD) {
        if (gotN) free_roots(rootsN, multsN, nN);
        if (gotD) free_roots(rootsD, multsD, nD);
        if (leadN) expr_free(leadN); if (leadD) expr_free(leadD);
        return NULL;
    }

    #define RI_CLEAN() do { free_roots(rootsN, multsN, nN); free_roots(rootsD, multsD, nD); \
                            expr_free(leadN); expr_free(leadD); } while (0)

    /* Require at least one genuine complex root; all-real cases are handled by
     * the existing Product`Infinite limit route.  (zeroq consumes its arg.) */
    bool any_complex = false;
    for (size_t i = 0; i < nN; i++) {
        Expr* im = fn1("Im", rootsN[i]);
        bool z = zeroq(im);   /* consumes im */
        if (!z) { any_complex = true; break; }
    }
    if (!any_complex) for (size_t i = 0; i < nD; i++) {
        Expr* im = fn1("Im", rootsD[i]);
        bool z = zeroq(im);
        if (!z) { any_complex = true; break; }
    }
    if (!any_complex) { RI_CLEAN(); return NULL; }

    /* --- Convergence gate. --- */
    /* c = leadN/leadD == 1 */
    Expr* c = prod_div(leadN, leadD);
    Expr* cm1 = expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ c, expr_new_integer(-1) }, 2);   /* adopts c */
    if (!zeroq(cm1)) { RI_CLEAN(); return NULL; }   /* cm1 consumed */

    /* degree balance */
    int degN = 0, degD = 0;
    for (size_t i = 0; i < nN; i++) degN += multsN[i];
    for (size_t i = 0; i < nD; i++) degD += multsD[i];
    if (degN != degD) { warn_div(); RI_CLEAN(); return NULL; }

    /* root-sum balance: sum m_i alpha_i == sum n_j beta_j */
    Expr* sumN = expr_new_integer(0);
    for (size_t i = 0; i < nN; i++) {
        Expr* t = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ expr_new_integer(multsN[i]), expr_copy(rootsN[i]) }, 2);
        Expr* s = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){ sumN, t }, 2);
        sumN = evaluate(s); expr_free(s);
    }
    Expr* sumD = expr_new_integer(0);
    for (size_t j = 0; j < nD; j++) {
        Expr* t = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ expr_new_integer(multsD[j]), expr_copy(rootsD[j]) }, 2);
        Expr* s = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){ sumD, t }, 2);
        sumD = evaluate(s); expr_free(s);
    }
    Expr* diffS = expr_new_function(expr_new_symbol(SYM_Plus),
                      (Expr*[]){ sumN, expr_new_function(expr_new_symbol(SYM_Times),
                                     (Expr*[]){ expr_new_integer(-1), sumD }, 2) }, 2); /* adopts */
    if (!zeroq(diffS)) { warn_div(); RI_CLEAN(); return NULL; }

    /* --- Build the Gamma multiset:  Gamma(imin - beta)^+n / Gamma(imin - alpha)^+... --- */
    size_t ntcap = nN + nD;
    GT* T = malloc(sizeof(GT) * (ntcap ? ntcap : 1));
    size_t nt = 0;
    Expr* rawG = expr_new_integer(1);   /* raw product for the numeric guard */
    bool bail = false;

    for (size_t j = 0; j < nD && !bail; j++) {   /* denominator -> +exp */
        Expr* arg = expr_new_function(expr_new_symbol(SYM_Plus),
                        (Expr*[]){ expr_copy(imin),
                                   expr_new_function(expr_new_symbol(SYM_Times),
                                       (Expr*[]){ expr_new_integer(-1), expr_copy(rootsD[j]) }, 2) }, 2);
        Expr* ae = evaluate(arg); expr_free(arg);
        T[nt].arg = ae; T[nt].exp = multsD[j]; nt++;
        Expr* g = fn1("Gamma", ae);
        rawG = mul_pow(rawG, g, multsD[j]); expr_free(g);
    }
    for (size_t i = 0; i < nN && !bail; i++) {   /* numerator -> -exp */
        Expr* arg = expr_new_function(expr_new_symbol(SYM_Plus),
                        (Expr*[]){ expr_copy(imin),
                                   expr_new_function(expr_new_symbol(SYM_Times),
                                       (Expr*[]){ expr_new_integer(-1), expr_copy(rootsN[i]) }, 2) }, 2);
        Expr* ae = evaluate(arg); expr_free(arg);
        T[nt].arg = ae; T[nt].exp = -multsN[i]; nt++;
        Expr* g = fn1("Gamma", ae);
        rawG = mul_pow(rawG, g, -multsN[i]); expr_free(g);
    }
    RI_CLEAN();
    #undef RI_CLEAN

    Expr* pref = expr_new_integer(1);

    /* (a) integer shift each arg into Re in (0,1]. */
    for (size_t i = 0; i < nt && !bail; i++) {
        int guard = 0;
        while (guard++ < 4096) {
            Expr* re = fn1("Re", T[i].arg);
            if (cmp_true("Greater", re, expr_new_integer(1) /* leaked copy below */)) {
                /* Gamma(z) = (z-1) Gamma(z-1). */
                expr_free(re);
                Expr* z1 = expr_new_function(expr_new_symbol(SYM_Plus),
                               (Expr*[]){ expr_copy(T[i].arg), expr_new_integer(-1) }, 2);
                Expr* z1e = evaluate(z1); expr_free(z1);
                pref = mul_pow(pref, z1e, T[i].exp);
                expr_free(T[i].arg); T[i].arg = z1e;
                continue;
            }
            bool le0 = cmp_true("LessEqual", re, expr_new_integer(0));
            expr_free(re);
            if (le0) {
                /* Gamma(z) = Gamma(z+1)/z -> factor z^(-exp); z must be nonzero. */
                if (zeroq(expr_copy(T[i].arg))) { bail = true; break; }
                pref = mul_pow(pref, T[i].arg, -T[i].exp);
                Expr* z1 = expr_new_function(expr_new_symbol(SYM_Plus),
                               (Expr*[]){ expr_copy(T[i].arg), expr_new_integer(1) }, 2);
                Expr* z1e = evaluate(z1); expr_free(z1);
                expr_free(T[i].arg); T[i].arg = z1e;
                continue;
            }
            break;
        }
    }

    /* Cache Re/Im per term. */
    Expr** RE = malloc(sizeof(Expr*) * (nt ? nt : 1));
    Expr** IM = malloc(sizeof(Expr*) * (nt ? nt : 1));
    for (size_t i = 0; i < nt; i++) { RE[i] = fn1("Re", T[i].arg); IM[i] = fn1("Im", T[i].arg); }

    /* (b) cancel equal args (opposite exponent). */
    for (size_t i = 0; i < nt; i++) {
        if (T[i].exp == 0) continue;
        for (size_t j = i + 1; j < nt; j++) {
            if (T[j].exp == 0) continue;
            Expr* d = expr_new_function(expr_new_symbol(SYM_Plus),
                          (Expr*[]){ expr_copy(T[i].arg),
                                     expr_new_function(expr_new_symbol(SYM_Times),
                                         (Expr*[]){ expr_new_integer(-1), expr_copy(T[j].arg) }, 2) }, 2);
            if (zeroq(d)) { T[i].exp += T[j].exp; T[j].exp = 0; if (T[i].exp == 0) break; }
        }
    }

    /* (c) conjugate pairs -> hyperbolic. */
    for (size_t i = 0; i < nt && !bail; i++) {
        if (T[i].exp == 0) continue;
        if (zeroq(expr_copy(IM[i]))) continue;   /* real term */
        /* find conjugate partner with equal exponent. */
        size_t part = nt;
        for (size_t j = 0; j < nt; j++) {
            if (j == i || T[j].exp != T[i].exp) continue;
            bool re_eq = zeroq(expr_new_function(expr_new_symbol(SYM_Plus),
                             (Expr*[]){ expr_copy(RE[i]),
                                        expr_new_function(expr_new_symbol(SYM_Times),
                                            (Expr*[]){ expr_new_integer(-1), expr_copy(RE[j]) }, 2) }, 2));
            bool im_opp = zeroq(expr_new_function(expr_new_symbol(SYM_Plus),
                             (Expr*[]){ expr_copy(IM[i]), expr_copy(IM[j]) }, 2));
            if (re_eq && im_opp) { part = j; break; }
        }
        if (part == nt) continue;   /* no partner -> residual (bails later) */
        Expr* bpos = fn1("Abs", IM[i]);
        Expr* pib = expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ expr_new_symbol("Pi"), expr_copy(bpos) }, 2);
        Expr* pibe = evaluate(pib); expr_free(pib);
        Expr* factor;
        if (cmp_true("Equal", RE[i], expr_new_integer(1))) {
            /* Gamma(1+ib)Gamma(1-ib) = Pi b / Sinh(Pi b). */
            Expr* sh = fn1("Sinh", pibe);
            Expr* pib2 = expr_new_function(expr_new_symbol(SYM_Times),
                             (Expr*[]){ expr_new_symbol("Pi"), expr_copy(bpos) }, 2);
            factor = prod_div(pib2, sh);
            expr_free(sh); expr_free(pib2);
        } else {
            /* Re == 1/2: Gamma(1/2+ib)Gamma(1/2-ib) = Pi / Cosh(Pi b). */
            Expr* half = expr_new_function(expr_new_symbol(SYM_Rational),
                             (Expr*[]){ expr_new_integer(1), expr_new_integer(2) }, 2);
            bool is_half = zeroq(expr_new_function(expr_new_symbol(SYM_Plus),
                              (Expr*[]){ expr_copy(RE[i]),
                                         expr_new_function(expr_new_symbol(SYM_Times),
                                             (Expr*[]){ expr_new_integer(-1), half }, 2) }, 2));
            if (!is_half) { expr_free(pibe); expr_free(bpos); bail = true; break; }
            Expr* ch = fn1("Cosh", pibe);
            Expr* pi = expr_new_symbol("Pi");
            factor = prod_div(pi, ch);
            expr_free(ch); expr_free(pi);
        }
        expr_free(pibe); expr_free(bpos);
        pref = mul_pow(pref, factor, T[i].exp);
        expr_free(factor);
        T[i].exp = 0; T[part].exp = 0;
    }

    /* (d) real reflection r & 1-r. */
    for (size_t i = 0; i < nt && !bail; i++) {
        if (T[i].exp == 0) continue;
        if (!zeroq(expr_copy(IM[i]))) continue;   /* only real terms */
        for (size_t j = i + 1; j < nt; j++) {
            if (T[j].exp != T[i].exp) continue;
            if (!zeroq(expr_copy(IM[j]))) continue;
            /* arg_i + arg_j == 1 ? */
            bool refl = zeroq(expr_new_function(expr_new_symbol(SYM_Plus),
                           (Expr*[]){ expr_copy(T[i].arg), expr_copy(T[j].arg),
                                      expr_new_integer(-1) }, 3));
            if (refl) {
                Expr* pir = expr_new_function(expr_new_symbol(SYM_Times),
                                (Expr*[]){ expr_new_symbol("Pi"), expr_copy(T[i].arg) }, 2);
                Expr* pire = evaluate(pir); expr_free(pir);
                Expr* sn = fn1("Sin", pire); expr_free(pire);
                Expr* pi = expr_new_symbol("Pi");
                Expr* factor = prod_div(pi, sn);
                expr_free(pi); expr_free(sn);
                pref = mul_pow(pref, factor, T[i].exp);
                expr_free(factor);
                T[i].exp = 0; T[j].exp = 0;
                break;
            }
        }
    }

    for (size_t i = 0; i < nt; i++) { expr_free(RE[i]); expr_free(IM[i]); }
    free(RE); free(IM);

    /* (e) residual: any leftover COMPLEX Gamma -> bail; keep real Gammas. */
    Expr* result = expr_copy(pref);
    for (size_t i = 0; i < nt && !bail; i++) {
        if (T[i].exp == 0) continue;
        Expr* im = fn1("Im", T[i].arg);
        bool real = zeroq(im);
        if (!real) { bail = true; break; }
        Expr* g = fn1("Gamma", T[i].arg);
        Expr* r2 = mul_pow(result, g, T[i].exp);   /* adopts result */
        expr_free(g);
        result = r2;
    }

    for (size_t i = 0; i < nt; i++) expr_free(T[i].arg);
    free(T);
    expr_free(pref);

    if (bail) { expr_free(result); expr_free(rawG); return NULL; }

    /* Collapse conjugate-factor products (e.g. (1/2-ib)(1/2+ib) -> real). */
    Expr* result_e = prod_eval("Simplify", (Expr*[]){ result }, 1);   /* adopts result */

    /* --- Numeric guard: |N[result] - N[rawG]| < 1e-9 (1 + |N[rawG]|). --- */
    Expr* diff = expr_new_function(expr_new_symbol(SYM_Plus),
                     (Expr*[]){ expr_copy(result_e),
                                expr_new_function(expr_new_symbol(SYM_Times),
                                    (Expr*[]){ expr_new_integer(-1), expr_copy(rawG) }, 2) }, 2);
    Expr* absd = expr_new_function(expr_new_symbol("Abs"), (Expr*[]){ diff }, 1); /* adopts */
    Expr* nabs = fn1("N", absd); expr_free(absd);
    bool ok = cmp_true("Less", nabs, expr_new_real(1e-9));   /* adopts tol */
    expr_free(nabs); expr_free(rawG);

    if (!ok) { expr_free(result_e); return NULL; }
    return result_e;
}

void product_rational_infinite_init(void) {
    symtab_add_builtin("Product`RationalInfinite", builtin_product_rational_infinite);
    symtab_get_def("Product`RationalInfinite")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Product`RationalInfinite",
        "Product`RationalInfinite[f, i, imin, Infinity] evaluates convergent "
        "infinite rational products with complex-conjugate roots via the Gamma "
        "canonical form (reflection / conjugate->hyperbolic reduction). Fires "
        "only when a complex root is present; returns unevaluated otherwise.");
}
