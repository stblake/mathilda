#include "facpoly.h"
#include "common.h"           /* alg_thread_call: Factor threads into Rule/relations */
#include "poly.h"
#include "eval.h"
#include "expand.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"
#include "internal.h"
#include "parse.h"
#include "rationalize.h"
#include "zupoly.h"
#include "bpoly.h"
#include "mpoly.h"
#include "mvfactor.h"
#include "mvfactor3.h"
#include "sym_names.h"
#include "qafactor.h"
#include "flint_bridge.h"
#include "flint_trager.h"
#include "print.h"
#include "message.h"
#include "numeric.h"
#include "checked_int.h"      /* ci_mul_i64 / ci_sub_i64: upoly_div_exact_z */
#include "core.h"             /* tc_deadline_passed(): bz_ddf / recombination */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <gmp.h>

/* facpoly.c is split into a set of .inc fragments included into this single
 * translation unit.  All in-fragment statics keep file scope; the forward
 * declarations below let the fragments be #included in any order.
 *
 *   facpoly_squarefree.inc        FactorSquareFree pipeline.
 *   facpoly_simple.inc            int_root, binomial, degree-one heuristics.
 *   facpoly_monomial_content.inc  Variable-monomial GCD extraction.
 *   facpoly_irreducibility.inc    Image-based irreducibility probe + factor_roots.
 *   facpoly_bv_hensel.inc         Bivariate Hensel lift (Stages 1-3).
 *   facpoly_z_split.inc           n-variate specialise-and-divide.
 *   facpoly_tv_hensel.inc         Trivariate Hensel via MPoly.
 *   facpoly_heuristic.inc         heuristic_factor recursive dispatcher.
 *   facpoly_memo.inc              Per-Simplify Factor result cache.
 *   facpoly_factor_builtin.inc    builtin_factor entry point.
 *   facpoly_factorterms.inc       FactorTerms / FactorTermsList.
 *   facpoly_bz_uni.inc            Univariate Berlekamp-Zassenhaus pipeline.
 */
static Expr* heuristic_factor(Expr* P);
static Expr* factor_monomial_content(Expr* P);
static Expr* factor_binomial(Expr* P);
static Expr* factor_degree_one(Expr* P, Expr** vars, size_t v_count);
static Expr* factor_bivariate_via_hensel(Expr* P, Expr** vars);
static Expr* factor_via_z_independent_split(Expr* P, Expr** vars, size_t v_count);
static Expr* factor_trivariate_via_mhensel(Expr* P, Expr** vars, size_t v_count);
static Expr* factor_roots(Expr* P, Expr** vars, size_t v_count);
static bool is_likely_irreducible_multivariate(Expr* P, Expr** vars, size_t v_count);
static Expr* eval_others_at_alpha(Expr* P, Expr** vars, size_t v_count,
                                  size_t main_idx, int64_t alpha);
static bool univariate_squarefree(Expr* u, Expr* var);
static FactorMemo* factor_memo_top(void);
static const Expr* factor_memo_get(FactorMemo* m, Expr* key);
static void factor_memo_put(FactorMemo* m, Expr* key, Expr* value);

Expr* bz_factor_to_expr(Expr* P, Expr* var);

/* True if `e` still carries a genuine denominator -- a Power[base, k] with k a
 * negative integer (or a Rational with negative numerator) over a non-numeric
 * base (a symbol or compound; a numeric base's negative power is just a rational
 * constant, not a denominator to combine).  builtin_factor derives num/den from
 * Together, which is contracted to return a single reduced fraction, so a clean
 * num/den is denominator-free.  When Together declines to combine (e.g. it bails
 * on >=2 dependent trig/hyperbolic kernels over Q(i) to avoid a GCD blow-up,
 * leaving a sum of fractions), the "numerator" still contains such a Power and is
 * NOT a polynomial -- feeding it to the polynomial factorer treats 1/(x-I) as a
 * variable and yields a NON-EQUIVALENT result (e.g.
 * Factor[(Cosh[x]+Sinh[x])/(-I+x) + (Cosh[x]+Sinh[x])/(I+x)] -> 0).  Detecting
 * this lets Factor decline safely (return the input unchanged) instead of
 * emitting a wrong answer. */
static bool factor_has_leftover_denominator(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h && h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Power &&
        e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* ex   = e->data.function.args[1];
        bool neg = (ex->type == EXPR_INTEGER && ex->data.integer < 0);
        if (!neg && ex->type == EXPR_FUNCTION && ex->data.function.head &&
            ex->data.function.head->type == EXPR_SYMBOL &&
            ex->data.function.head->data.symbol.name == SYM_Rational &&
            ex->data.function.arg_count == 2 &&
            ex->data.function.args[0]->type == EXPR_INTEGER &&
            ex->data.function.args[0]->data.integer < 0)
            neg = true;
        if (neg && base && (base->type == EXPR_SYMBOL || base->type == EXPR_FUNCTION))
            return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (factor_has_leftover_denominator(e->data.function.args[i])) return true;
    return false;
}

#include "expsplit.h"   /* FactorSquareFree reads symbolic powers over generators */
#include "facpoly_options.inc"
#include "facpoly_squarefree.inc"
#include "facpoly_simple.inc"
#include "facpoly_monomial_content.inc"
#include "facpoly_irreducibility.inc"
#include "facpoly_bv_hensel.inc"
#include "facpoly_z_split.inc"
#include "facpoly_tv_hensel.inc"
#include "facpoly_heuristic.inc"
#include "facpoly_memo.inc"
#include "expsplit.h"
#include "facpoly_factor_builtin.inc"
#include "facpoly_factorterms.inc"
#include "facpoly_bz_uni.inc"

/* ===================================================================== */
/* PolynomialSqrt[p] / PolynomialSqrt[p, x] — the polynomial square root.  */
/* ===================================================================== */
/* Returns s with s^2 == p exactly when p is a perfect square (every
 * non-constant irreducible factor occurs to an even multiplicity), else the
 * sentinel `$Failed`.  The numeric content is carried through Sqrt (so
 * PolynomialSqrt[2(1+x)^2] = Sqrt[2] (1+x) over the reals / F-bar), which keeps
 * the routine total over an algebraically closed constant field — exactly the
 * setting Cherry's completing-square Erf-argument step needs (r_i = Sqrt[p+beta q]).
 * Correct by construction: the candidate s is accepted only after the exact
 * Expand[s^2 - p] == 0 check, so a mis-factorisation can never return a wrong root. */

static bool ps_is_numeric(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_BIGINT: case EXPR_REAL: return true;
        case EXPR_FUNCTION:
            return e->data.function.head->type == EXPR_SYMBOL &&
                   (e->data.function.head->data.symbol.name == SYM_Rational ||
                    e->data.function.head->data.symbol.name == SYM_Complex);
        default: return false;
    }
}

static Expr* ps_call1(const char* h, Expr* a) {
    return eval_and_free(expr_new_function(expr_new_symbol(h), &a, 1));
}
static Expr* ps_call2(const char* h, Expr* a, Expr* b) {
    return eval_and_free(expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2));
}

/* The square root of a univariate p in x by undetermined coefficients, from
 * the top: s_d = Sqrt[a_2d], then
 *     s_k = (a_(d+k) - Sum_{i+j = d+k, k < i,j < d} s_i s_j) / (2 s_d),
 * each simplified (Together), and the candidate accepted only after the exact
 * Expand[s^2 - p] == 0 certificate.  This is the route for coefficients the
 * factoriser over Q cannot split -- algebraic numbers (x^4 + 2 Sqrt[2] x^2 + 2
 * is (x^2 + Sqrt[2])^2, irreducible over Q as written).  NULL when it fails. */
static Expr* ps_by_coefficients(Expr* p, Expr* x) {
    Expr* cl = ps_call2("CoefficientList", expr_copy(p), expr_copy(x));
    if (!cl || cl->type != EXPR_FUNCTION || cl->data.function.head->type != EXPR_SYMBOL ||
        cl->data.function.head->data.symbol.name != SYM_List) { if (cl) expr_free(cl); return NULL; }
    size_t len = cl->data.function.arg_count;
    if (len < 3 || len % 2 == 0 || len > 2001) { expr_free(cl); return NULL; }
    for (size_t i = 0; i < len; i++)
        if (contains_any_symbol_from(cl->data.function.args[i], x)) { expr_free(cl); return NULL; }
    size_t d = (len - 1) / 2;
    Expr** a = cl->data.function.args;
    Expr** sk = calloc(d + 1, sizeof(Expr*));
    sk[d] = ps_call2("Power", expr_copy(a[2 * d]),
                     ps_call2("Rational", expr_new_integer(1), expr_new_integer(2)));
    Expr* two_sd = ps_call2("Times", expr_new_integer(2), expr_copy(sk[d]));
    for (size_t k = d; k-- > 0; ) {
        Expr* acc = expr_copy(a[d + k]);
        for (size_t i = k + 1; i < d; i++) {
            size_t j = d + k - i;
            if (j <= k || j >= d) continue;
            acc = ps_call2("Plus", acc, ps_call2("Times", expr_new_integer(-1),
                           ps_call2("Times", expr_copy(sk[i]), expr_copy(sk[j]))));
        }
        sk[k] = ps_call1("Together", ps_call2("Times", acc,
                         ps_call2("Power", expr_copy(two_sd), expr_new_integer(-1))));
    }
    expr_free(two_sd);
    expr_free(cl);
    Expr** terms = malloc((d + 1) * sizeof(Expr*));
    for (size_t k = 0; k <= d; k++)
        terms[k] = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ sk[k], expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(x), expr_new_integer((int64_t)k) }, 2) }, 2);
    free(sk);
    Expr* s = eval_and_free(expr_new_function(expr_new_symbol("Plus"), terms, d + 1));
    free(terms);
    Expr* chk = ps_call1("Expand", ps_call2("Plus", ps_call2("Power", expr_copy(s), expr_new_integer(2)),
                         ps_call2("Times", expr_new_integer(-1), expr_copy(p))));
    bool good = chk && is_zero_poly(chk);
    if (chk) expr_free(chk);
    if (!good) { expr_free(s); return NULL; }
    return s;
}

Expr* builtin_polynomialsqrt(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1
        || res->data.function.arg_count > 2)
        return NULL;
    Expr* p = res->data.function.args[0];
    /* In the 2-arg form PolynomialSqrt[p, x], any factor free of the named
     * variable x is CONSTANT content — carry it through Sqrt just like the
     * numeric content, so a symbolic/algebraic leading coefficient
     * (Sqrt[-c] for -c x^2, Cherry's Erf argument over an algebraically
     * closed constant field) is handled rather than rejected as a "bare
     * irreducible".  The 1-arg form has no variable, so its numeric-only
     * content behaviour is unchanged. */
    Expr* psx = (res->data.function.arg_count == 2)
        ? res->data.function.args[1] : NULL;

    /* 0 -> 0. */
    Expr* pe = expr_expand(p);
    if (is_zero_poly(pe)) { expr_free(pe); return expr_new_integer(0); }
    expr_free(pe);

    /* Factor into  constant * prod(base^exp). */
    Expr** fa = malloc(sizeof(Expr*)); fa[0] = expr_copy(p);
    Expr* f = internal_factor(fa, 1); free(fa);
    if (!f) return NULL;

    Expr** items; size_t nit; Expr* single[1];
    if (f->type == EXPR_FUNCTION && f->data.function.head->type == EXPR_SYMBOL
        && f->data.function.head->data.symbol.name == SYM_Times) {
        items = f->data.function.args; nit = f->data.function.arg_count;
    } else { single[0] = f; items = single; nit = 1; }

    Expr* konst = expr_new_integer(1);                       /* numeric content */
    Expr** halves = malloc((nit ? nit : 1) * sizeof(Expr*)); size_t nh = 0;
    bool ok = true;
    for (size_t i = 0; i < nit && ok; i++) {
        Expr* t = items[i];
        bool t_content = ps_is_numeric(t)
            || (psx && !contains_any_symbol_from(t, psx));  /* x-free => content */
        if (t_content) {
            konst = eval_and_free(expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ konst, expr_copy(t) }, 2));
        } else if (t->type == EXPR_FUNCTION
                   && t->data.function.head->type == EXPR_SYMBOL
                   && t->data.function.head->data.symbol.name == SYM_Power
                   && t->data.function.arg_count == 2
                   && t->data.function.args[1]->type == EXPR_INTEGER) {
            int64_t e = t->data.function.args[1]->data.integer;
            Expr* base = t->data.function.args[0];
            if (ps_is_numeric(base)
                || (psx && !contains_any_symbol_from(base, psx))) {
                konst = eval_and_free(expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ konst, expr_copy(t) }, 2));
            } else if (e % 2 == 0) {
                halves[nh++] = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(base), expr_new_integer(e / 2) }, 2);
            } else {
                ok = false;                       /* odd multiplicity, non-constant */
            }
        } else {
            ok = false;                           /* bare irreducible (multiplicity 1) */
        }
    }
    expr_free(f);

    Expr* result = NULL;
    if (ok) {
        Expr* sk = eval_and_free(expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ konst, expr_new_function(expr_new_symbol("Rational"),
                (Expr*[]){ expr_new_integer(1), expr_new_integer(2) }, 2) }, 2));
        Expr** facs = malloc((nh + 1) * sizeof(Expr*)); size_t nf = 0;
        facs[nf++] = sk;
        for (size_t i = 0; i < nh; i++) facs[nf++] = halves[i];
        Expr* s = (nf == 1) ? facs[0]
                : eval_and_free(expr_new_function(expr_new_symbol("Times"), facs, nf));
        free(facs);
        /* Exact certificate: Expand[s^2 - p] == 0. */
        Expr* chk = eval_and_free(expr_new_function(expr_new_symbol("Expand"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(s), expr_new_integer(2) }, 2),
                  expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(p) }, 2) }, 2) }, 1));
        if (chk && is_zero_poly(chk)) result = s; else expr_free(s);
        if (chk) expr_free(chk);
    } else {
        expr_free(konst);
        for (size_t i = 0; i < nh; i++) expr_free(halves[i]);
    }
    free(halves);
    if (!result) {
        /* The factoriser could not show it (algebraic-number coefficients):
         * the univariate coefficient method, in the named variable or the
         * only one p has. */
        Expr* var = psx ? expr_copy(psx) : NULL;
        if (!var) {
            Expr* vs = ps_call1("Variables", expr_copy(p));
            if (vs && vs->type == EXPR_FUNCTION && vs->data.function.arg_count == 1)
                var = expr_copy(vs->data.function.args[0]);
            if (vs) expr_free(vs);
        }
        if (var) {
            result = ps_by_coefficients(p, var);
            expr_free(var);
        }
    }
    return result ? result : expr_new_symbol("$Failed");
}

void facpoly_init(void) {
    /* The per-Simplify memo stack: a timeout unwinding out of a Simplify must
     * not leave its memo active for the rest of the session (facpoly_memo.inc). */
    tc_register_unwind_counter(factor_memo_depth_counter());
    symtab_add_builtin("FactorSquareFree", builtin_factorsquarefree);
    symtab_get_def("FactorSquareFree")->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
    symtab_add_builtin("Factor", builtin_factor);
    symtab_get_def("Factor")->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
    symtab_add_builtin("FactorTerms", builtin_factorterms);
    symtab_get_def("FactorTerms")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("FactorTermsList", builtin_factortermslist);
    symtab_get_def("FactorTermsList")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("PolynomialSqrt", builtin_polynomialsqrt);
    symtab_get_def("PolynomialSqrt")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("PolynomialSqrt",
        "PolynomialSqrt[p] gives a polynomial s with s^2 == p when p is a perfect "
        "square (every non-constant irreducible factor has even multiplicity; the "
        "numeric content is carried through Sqrt), and $Failed otherwise. "
        "PolynomialSqrt[p, x] treats p as a polynomial in x.");
}
